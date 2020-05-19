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
import re
import covid_profiler as cp


def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","version":".version.json","proteins":".proteins.csv"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf


def get_codon_num(mutation):
    re_obj = re.search("([\d]+)[A-Z\*]",mutation)
    if re_obj:
        return int(re_obj.group(1))
    else:
        return None

def get_conf(prefix):
    conf = {}
    for f,s in [("ref",".fasta"),("barcode",".barcode.bed"),("gff",".gff"),("proteins",".proteins.csv")]:
        conf[f] = "%s%s" % (prefix,s)
    return conf


def get_sample_meta(samples,debug=False):
    if not args.debug:
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

def change_codon_number(mut,num):
    re_obj = re.search("(\d+)([A-Z\*]+)>(\d+)([A-Z\*]+)",mut)
    if re_obj:
        return "%s%s>%s%s" % (num,re_obj.group(2),num,re_obj.group(4))
    re_obj = re.search("(\d+)([A-Z\*]+)",mut)
    if re_obj:
        return "%s%s" % (num,re_obj.group(2))
    else:
        import pdb; pdb.set_trace()
def get_variant_data(vcf_file,ref_file,gff_file,protein_file):
    nsp_data = {}
    gene_info = {}
    for row in csv.DictReader(open(protein_file)):
        row["Start"] = int(row["Start"])
        row["End"] = int(row["End"])
        gene_info[row["Gene"]] = {"function":row["Putative function"],"DOI":row["DOI"]}
        if row["Region"]!="nsp": continue
        for i in range(row["Start"],row["End"]+1):
            nsp_data[i] = row

    pp.run_cmd("samtools faidx %s" % ref_file)
    results = {}
    for l in pp.cmd_out("bcftools view %s | bcftools csq -f %s -g %s  | correct_covid_csq.py | bcftools +fill-tags | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%AF\\t%%BCSQ\\n'" % (vcf_file,ref_file,gff_file)):
        pos,ref,alts_str,af_str,csq_str = l.strip().split()
        alt_af = sum([float(x) for x in af_str.split(",")])
        csqs = csq_str.split(",")
        types = []
        changes = []
        genes = []
        pos = int(pos)
        for i in range(len(csqs)):
            if csqs[i][0]=="@":
                results[pos] = results[int(csqs[i][1:])]

            elif csqs[i]==".":
                results[pos] = {"pos":pos, "alts":alts_str, "alt_af":alt_af, "types":"intergenic","changes":"NA","gene":"NA","gene_function":"NA","gene_reference":"NA"}

            else:
                csq = csqs[i].split("|")
                types.append(csq[0].replace("*",""))

                if csq[1]=="orf1ab":
                    codon_pos = get_codon_num(csq[5])
                    if codon_pos in nsp_data:
                        genes.append(nsp_data[codon_pos]["Gene"])
                        codon_pos = codon_pos-nsp_data[codon_pos]["Start"]+1
                        changes.append(change_codon_number(csq[5],codon_pos))
                    else:
                        genes.append("orf1ab")
                        changes.append(csq[5])
                else:
                    changes.append(csq[5])
                    genes.append(csq[1])
                if len(set(types))==1:
                    types = list(set(types))
                results[pos] = {"pos":pos, "alts":alts_str, "alt_af":alt_af, "types":",".join(types), "changes":",".join(changes),"gene":genes[0], "gene_function":gene_info[genes[0]]["function"], "gene_reference":gene_info[genes[0]]["DOI"]}
    return results

def fasta2vcf(fasta_file,outfile):
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    refseq = pp.fasta(conf["ref"]).fa_dict
    seqs = pp.fasta(fasta_file)
    samples = list(seqs.fa_dict.keys())

    for sample in samples:
        fname = pp.get_random_file()
        open(fname,"w").write(">%s\n%s\n" % (sample,seqs.fa_dict[sample]))
        fasta_obj = pp.fasta(fname)
        vcf_obj = pp.vcf(fasta_obj.get_ref_variants(conf["ref"], sample))
        pp.run_cmd("rm %s" % fname)

    sample_chunks = [samples[i:i+200] for i in range(0,len(samples),200)]
    tmp_vcfs = []
    for tmp_samps in sample_chunks:
        tmp_list = pp.get_random_file()
        tmp_vcf = pp.get_random_file()
        open(tmp_list,"w").write("\n".join(["%s.vcf.gz" % x for x in tmp_samps]))
        pp.run_cmd("bcftools merge -0 -l %s -Oz -o %s" % (tmp_list,tmp_vcf))
        pp.run_cmd("bcftools index %s" % tmp_vcf)
        tmp_vcfs.append(tmp_vcf)
        pp.rm_files([tmp_list])

    pp.run_cmd("bcftools merge -0  %s | bcftools view -V indels -Oz -o %s" % (" ".join(tmp_vcfs),outfile))

    vcf_files = ["%s.vcf.gz" % s  for s in samples]
    vcf_csi_files = ["%s.vcf.gz.csi" % s  for s in samples]
    pp.rm_files(vcf_files + vcf_csi_files + tmp_vcfs)

def main_fasta2vcf(args):
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    print(conf)
    fasta2vcf(args.seqs,args.outfile)

def main_asr(args):
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    variant_data = get_variant_data(args.vcf,conf["ref"],conf["gff"],conf["proteins"])
    mutations = cp.find_ancestral_mutations(args.msa,args.tree,args.states,variant_sites = set(list(variant_data)))
    for i in range(len(mutations)):
        for key in variant_data[mutations[i]["position"]]:
            mutations[i][key] = variant_data[mutations[i]["position"]][key]

    with open(args.out+".csv","w") as O:
        tmp = list(mutations[0].keys())
        fieldnames = tmp[:4] + tmp[-6:] + tmp[4:-6]
        writer = csv.DictWriter(O,fieldnames = fieldnames)
        writer.writeheader()
        writer.writerows(mutations)



def main_preprocess(args):
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    refseq = pp.fasta(conf["ref"]).fa_dict
    refseqname = list(refseq.keys())[0]
    if args.seqs:
        seqs = pp.fasta(args.seqs)

    samples = list(seqs.fa_dict.keys())


    fasta2vcf(args.seqs,"merged.vcf.gz")
    #
    #
    # pp.run_cmd("vcf2fasta.py --vcf merged.vcf.gz --ref %s" % conf["ref"])

    # pp.run_cmd("iqtree -s merged.fa -m GTR+F+R2 -bb 1000 -nt AUTO -asr -czb -redo")
    variant_data = get_variant_data("merged.vcf.gz",conf["ref"],conf["gff"],conf["proteins"])
    #### variant_data = {29366: {'alts': 'T', 'types': 'missense', 'changes': '365P>365S', 'gene': 'N', 'gene_function': 'Nucleocapsid protein', 'gene_reference': '10.1186/s40779-020-00240-0'}}


    #### mutations  = [{'position': 29868, 'mutation_type': 'convergent', 'origins': 14, 'branches': 'Node927'}]
    for i in range(len(mutations)):
        for key in variant_data[mutations[i]["position"]]:
            mutations[i][key] = variant_data[mutations[i]["position"]][key]

    with open(args.out+".csv","w") as O:
        tmp = list(mutations[0].keys())
        fieldnames = tmp[:4] + tmp[-6:] + tmp[4:-6]
        writer = csv.DictWriter(O,fieldnames = fieldnames)
        writer.writeheader()
        writer.writerows(mutations)
    quit()
    seqs = pp.fasta("merged.fa").fa_dict

    tree = ete3.Tree("merged.fa.treefile",format=1)
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
    for l in tqdm(open("merged.fa.state")):
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
        if site not in variant_data: continue ####### Not sure if this is right yet
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
            "gene_function": variant_data[site]["gene_function"],
            "gene_reference": variant_data[site]["gene_reference"],
            "alts": variant_data[site]["alts"],
            "functional_types": variant_data[site]["types"],
            "changes": variant_data[site]["changes"]
            }
        for sample in leaf_names:
            tmp_data[sample] = states[site][sample]
        mutations.append(tmp_data)

    print("Barcoding sites: ",barcoding_sites)
    print("Convergent sites: ",convergent_sites)

    # Reroot tree at S/L types
    outgroup_leaf_names = [s for s in leaf_names if seqs[s][8782-1]=="T"]
    tree.set_outgroup(tree.get_common_ancestor(outgroup_leaf_names))


    tree.write(format=1, outfile=args.out+".tree")

    with open(args.out+".barcode.bed","w") as O:
        for pos in barcoding_sites:
            for allele in set(list(states[pos].values())):
                tmp_samps = [x for x in leaf_names if states[pos][x]==allele]
                O.write("%s\t%s\t%s\t%s\t%s\n" % (refseqname,pos-1,pos,allele,tree.get_common_ancestor(tmp_samps).name))

    with open(args.out+".mutation_summary.csv","w") as O:
        writer = csv.DictWriter(O,fieldnames=["position","mutation_type","origins","branches","gene","gene_function","gene_reference","alts","functional_types","changes"] + list(leaf_names))
        writer.writeheader()
        for row in mutations:
            writer.writerow(row)


def main_load_library(args):
    lib_prefix = args.prefix.split("/")[-1]
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","version":".version.json"}
    if pp.nofolder(sys.base_prefix+"/share/covidprofiler"):
        pp.run_cmd("mkdir %s " % (sys.base_prefix+"/share/covidprofiler/"))
    for key in files:
        new_file_location = sys.base_prefix+"/share/covidprofiler/"+lib_prefix+files[key]
        pp.run_cmd("cp %s %s" % (args.prefix+files[key],new_file_location))
    pp.run_cmd("samtools faidx %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".fasta")
    pp.run_cmd("bwa index %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".fasta")
    if os.path.isfile("%s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".dict"):
        pp.run_cmd("rm %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".dict")
    pp.run_cmd("gatk CreateSequenceDictionary -R %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".fasta")
    pp.log("Sucessfully imported library")


def main_get_mutation_data(args):
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    variant_data = list(get_variant_data(args.vcf,conf["ref"],conf["gff"],conf["proteins"]).values())
    with open(args.vcf+".annotations.csv","w") as O:
        fieldnames = list(variant_data[0].keys())
        writer = csv.DictWriter(O,fieldnames)
        writer.writeheader()
        writer.writerows(variant_data)







parser = argparse.ArgumentParser(description='Covid pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")


parser_sub = subparsers.add_parser('preprocess', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--seqs',help='First read file',required=True)
parser_sub.add_argument('--dir',default="/tmp/",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.add_argument('--out',default="covid_public",help='First read file')
parser_sub.add_argument('--debug',action="store_true",help='First read file')
parser_sub.set_defaults(func=main_preprocess)

parser_sub = subparsers.add_parser('fasta2vcf', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--seqs',help='First read file',required=True)
parser_sub.add_argument('--outfile',default="gisaid.vcf.gz",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.set_defaults(func=main_fasta2vcf)

parser_sub = subparsers.add_parser('asr', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--msa',help='First read file',required=True)
parser_sub.add_argument('--tree',help='First read file',required=True)
parser_sub.add_argument('--states',help='First read file',required=True)
parser_sub.add_argument('--vcf',help='First read file',required=True)
parser_sub.add_argument('--out',default="gisaid.asr",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.set_defaults(func=main_asr)


parser_sub = subparsers.add_parser('mutation_data', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--vcf',help='First read file',required=True)
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.set_defaults(func=main_get_mutation_data)

args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
