import pathogenprofiler as pp
import argparse
import ete3
import sys
from tqdm import tqdm
from collections import defaultdict
import csv
import os
from Bio import SeqIO
import json
import os.path


def get_conf(prefix):
    conf = {}
    for f,s in [("ref",".fa"),("barcode",".barcode.bed"),("gff",".gff")]:
        conf[f] = "%s%s" % (prefix,s)
    return conf


def get_sample_meta(samples):
    pp.run_cmd("esearch -db nucleotide -query '%s' | efetch -format gb  > temp.gb" % ",".join(samples))
    data = []
    for seq_record in SeqIO.parse(open("temp.gb"), "gb"):
        sample = seq_record.id.split(".")[0]
        source = [feat for feat in seq_record.features if feat.type=="source"][0]
        country = "NA"
        date = "NA"
        if "country" in source.qualifiers:
            country = source.qualifiers["country"][0].split(":")[0]
        if "collection_date" in source.qualifiers:
            date = source.qualifiers["collection_date"][0]
        data.append({"id":sample,"country":country,"date":date})

    return data

def get_variant_data(vcf_file,ref_file,gff_file):
    pp.run_cmd("samtools faidx %s" % ref_file)
    results = {}
    for l in pp.cmd_out("bcftools view %s | bcftools csq -f %s -g %s  | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%BCSQ\\n'" % (vcf_file,ref_file,gff_file)):
        pos,ref,alts_str,csq_str = l.strip().split()
        if csq_str==".":
            results[int(pos)] = {"alts":alts_str, "types":"intergenic","changes":"NA","gene":"NA"}
        else:
            csqs = csq_str.split(",")
            types = []
            changes = []
            genes = []
            for i in range(len(csqs)):
                csq = csqs[i].split("|")
                types.append(csq[0])
                changes.append(csq[5])
                genes.append(csq[1])
            results[int(pos)] = {"alts":alts_str, "types":",".join(types), "changes":",".join(changes),"gene":genes[0]}
    return results

def main_preprocess(args):
    os.chdir(args.dir)
    conf = get_conf(args.db)
    refseq = pp.fasta(conf["ref"]).fa_dict
    refseqname = list(refseq.keys())[0]

    pp.run_cmd("curl 'https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?q=*:*&fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&fq=%7B!tag=Flags_csv%7DFlags_csv:%22complete%22&cmd=download&sort=&dlfmt=fasta&fl=id,Definition_s,Nucleotide_seq' > temp.fa")
    pp.run_cmd("samtools faidx temp.fa")

    seqs = pp.fasta("temp.fa")
    samples = list(seqs.fa_dict.keys())
    for sample in samples:
        fname = pp.get_random_file()
        open(fname,"w").write(">%s\n%s\n" % (sample,seqs.fa_dict[sample]))
        fasta_obj = pp.fasta(fname)
        vcf_obj = pp.vcf(fasta_obj.get_ref_variants(conf["ref"], sample))
        pp.run_cmd("rm %s" % fname)

    vcf_files = ["%s.vcf.gz" % s  for s in samples]
    vcf_csi_files = ["%s.vcf.gz.csi" % s  for s in samples]
    pp.run_cmd("bcftools merge -0  %s | bcftools view -V indels -Oz -o merged.vcf.gz" % (" ".join(vcf_files)))
    pp.run_cmd("rm %s" % (" ".join(vcf_files)))
    pp.run_cmd("rm %s" % (" ".join(vcf_csi_files)))
    pp.run_cmd("vcf2fasta.py --vcf merged.vcf.gz --ref %s" % conf["ref"])
    if os.path.isfile("merged.fa.log"):
        pp.run_cmd("rm merged.fa.*")
    pp.run_cmd("iqtree -s merged.fa -bb 1000 -nt AUTO -czb -redo")
    variant_data = get_variant_data("merged.vcf.gz",conf["ref"],conf["gff"])
    sample_data = get_sample_meta(samples)
    with open("%s.meta.csv" % args.out,"w") as O:
        writer = csv.DictWriter(O,fieldnames=["id","country","date"])
        writer.writeheader()
        for row in sample_data:
            writer.writerow(row)

    seqs = pp.fasta("merged.fa").fa_dict
    tree = ete3.Tree("merged.fa.treefile",format=1)
    # Reroot tree at S/L types
    tree.set_outgroup(tree.get_common_ancestor("MN996527","MT106053"))
    outgroup_leaf_names = [s for s in tree.get_leaf_names() if seqs[s][8782-1]=="T"]
    tree.set_outgroup(tree.get_common_ancestor(outgroup_leaf_names))

    tree.write(format=1, outfile=args.out+".tree")
    if os.path.isfile("merged.fa.asr.log"):
        pp.run_cmd("rm merged.fa.asr.*")
    pp.run_cmd("iqtree -s merged.fa -te %s.tree -nt AUTO -czb -pre merged.fa.asr -asr" % (args.out))



    tree = ete3.Tree("merged.fa.asr.treefile",format=1)
    node_names = set([tree.name] + [n.name.split("/")[0] for n in tree.get_descendants()])
    leaf_names = set(tree.get_leaf_names())
    internal_node_names = node_names - leaf_names

    for n in tree.traverse():
        if n.name.split("/")[0] in node_names:
            if "Node" in n.name:
                tmp = n.name.split("/")
                n.name = tmp[0]
                if len(tmp)>1:
                    n.support = tmp[1]



    states = defaultdict(dict)
    sites = set()
    sys.stderr.write("Loading states\n")
    for l in tqdm(open("merged.fa.asr.state")):
        if l[0]=="#": continue
        row = l.strip().split()
        if row[0]=="Node": continue
        states[int(row[1])][row[0]] = row[2]
        sites.add(int(row[1]))

    sys.stderr.write("Loading alignment sites\n")
    for site in tqdm(sites):
        for sample in seqs:
            states[site][sample] = seqs[sample][site-1]

    barcoding_sites = []
    convergent_sites = []
    mutations = []
    for site in tqdm(sites):
        nucleotides = set([states[site][n] for n in node_names])
        if len(nucleotides)==1: continue

        # Set up storage objjects
        origins = []
        internal_node_change = False

        tree.add_feature("state",states[site][tree.name])
        for n in tree.traverse():
            if n == tree: continue
            node_state = states[site][n.name]
            if node_state!=n.get_ancestors()[0].state:
                origins.append(n.name)
                if n.name in internal_node_names:
                     internal_node_change = True
            n.add_feature("state",node_state)

        type = "unique"
        if internal_node_change and len(origins)==1:
            type = "barcoding"
            barcoding_sites.append(site)
        if len(origins)>1:
            type = "convergent"
            convergent_sites.append(site)

        tmp_data = {
            "position": site,
            "mutation_type":type,
            "origins": len(origins),
            "branches":",".join(origins),
            "gene": variant_data[site]["gene"],
            "alts": variant_data[site]["alts"],
            "functional_types": variant_data[site]["types"],
            "changes": variant_data[site]["changes"]
            }
        for sample in leaf_names:
            tmp_data[sample] = states[site][sample]
        mutations.append(tmp_data)

    print("Barcoding sites: ",barcoding_sites)
    print("Convergent sites: ",convergent_sites)



    with open(args.out+".barcode.bed","w") as O:
        for pos in barcoding_sites:
            for allele in set(list(states[pos].values())):
                tmp_samps = [x for x in leaf_names if states[pos][x]==allele]
                O.write("%s\t%s\t%s\t%s\t%s\n" % (refseqname,pos-1,pos,allele,tree.get_common_ancestor(tmp_samps).name))

    with open(args.out+".mutation_summary.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames=["position","mutation_type","origins","branches","gene","alts","functional_types","changes"] + list(leaf_names))
        writer.writeheader()
        for row in mutations:
            writer.writerow(row)


def main_position_sample(args):
    os.chdir(args.dir)
    conf = get_conf(args.db)
    refseq = pp.fasta(conf["ref"]).fa_dict
    refseqname = list(refseq.keys())[0]

    tree = ete3.Tree(args.tree, format=1)
    barcoding_sites = {}

    for l in open(args.barcode_bed):
        row = l.strip().split()
        barcoding_sites[row[4]] = (row[0],int(row[2]),row[3])

    fasta_obj = pp.fasta(args.fasta)
    sample_name = list(fasta_obj.fa_dict.keys())[0]
    vcf_obj = pp.vcf(fasta_obj.get_ref_variants(conf["ref"], sample_name))
    bed_gt = vcf_obj.get_bed_gt(args.barcode_bed,conf["ref"])

    closest_node = None
    for node in tree.traverse():
        if node.name in barcoding_sites:
            tmp = barcoding_sites[node.name]
            if tmp[2] in bed_gt[tmp[0]][tmp[1]]:
                closest_node = node

    sys.stdout.write("%s\n" % closest_node.name)
    open(args.outfile,"w").write("%s\n" % closest_node.get_ascii(attributes=["name"], show_internal=False))


def main_profile(args):
    os.chdir(args.dir)
    conf = get_conf(args.db)
    refseq = pp.fasta(conf["ref"]).fa_dict
    refseqname = list(refseq.keys())[0]

    fasta_obj = pp.fasta(args.fasta)
    sample_name = list(fasta_obj.fa_dict.keys())[0]
    wg_vcf_obj = pp.vcf(fasta_obj.get_ref_variants(conf["ref"], sample_name))
    mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])

    barcode = pp.barcode(mutations,conf["barcode"])
    sys.stdout.write("%s\t%s\n" % (sample_name,";".join([d["annotation"] for d in barcode])))


parser = argparse.ArgumentParser(description='Covid pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")


parser_sub = subparsers.add_parser('preprocess', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--dir',default="/tmp/",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.add_argument('--out',default="covid_public",help='First read file')
parser_sub.set_defaults(func=main_preprocess)


parser_sub = subparsers.add_parser('position_isolate', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--fasta',help='First read file',required=True)
parser_sub.add_argument('--tree',default="covid_public.tree",help='First read file')
parser_sub.add_argument('--barcode-bed',default="covid_public.barcode.bed",help='First read file')
parser_sub.add_argument('--dir',default="/tmp/",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.add_argument('--outfile',default="temp.tree.txt",help='First read file')
parser_sub.set_defaults(func=main_position_sample)

parser_sub = subparsers.add_parser('profile', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--fasta',help='First read file',required=True)
parser_sub.add_argument('--dir',default="/tmp/",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.set_defaults(func=main_profile)


args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
