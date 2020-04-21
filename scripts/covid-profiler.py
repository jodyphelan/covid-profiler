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

def parse_blast_result(result):
    results = json.load(open(result))["BlastOutput2"][0]["report"]["results"]["bl2seq"][0]["hits"][0]["hsps"]
    return results

def blast_seq(seq,ref,word_size=4):
    tmpseqfile = pp.get_random_file()
    open(tmpseqfile,"w").write(">query\n%s\n" % seq)
    tmpresult = pp.get_random_file()
    pp.run_cmd("blastn -task blastn -word_size %s -subject %s -query %s -outfmt 15 -evalue 10000000 -max_hsps 1 > %s" % (word_size,ref,tmpseqfile,tmpresult))
    blast_results = parse_blast_result(tmpresult)
    pp.rm_files([tmpresult,tmpseqfile])
    return blast_results


def get_conf_dict(library_prefix):
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","version":".version.json"}
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

def get_variant_data(vcf_file,ref_file,gff_file,protein_file):
    nsp_data = {}
    gene_info = {}
    for row in csv.DictReader(open(protein_file)):
        gene_info[row["Gene"]] = {"function":row["Putative function"],"DOI":row["DOI"]}
        if row["Region"]!="nsp": continue
        for i in range(int(row["Start"]),int(row["End"])+1):
            nsp_data[i] = row

    pp.run_cmd("samtools faidx %s" % ref_file)
    results = {}
    for l in pp.cmd_out("bcftools view %s | bcftools csq -f %s -g %s  | bcftools query -f '%%POS\\t%%REF\\t%%ALT\\t%%BCSQ\\n'" % (vcf_file,ref_file,gff_file)):
        pos,ref,alts_str,csq_str = l.strip().split()
        pos = int(pos)
        if csq_str==".":
            results[pos] = {"alts":alts_str, "types":"intergenic","changes":"NA","gene":"NA","gene_function":"NA","gene_reference":"NA"}
        else:
            csqs = csq_str.split(",")
            types = []
            changes = []
            genes = []
            for i in range(len(csqs)):
                csq = csqs[i].split("|")
                types.append(csq[0])
                changes.append(csq[5])
                if csq[1]=="orf1ab":
                    codon_pos = get_codon_num(csq[5])
                    if codon_pos in nsp_data:
                        genes.append(nsp_data[codon_pos]["Gene"])
                    else:
                        genes.append("orf1ab")
                else:
                    genes.append(csq[1])

            results[pos] = {"alts":alts_str, "types":",".join(types), "changes":",".join(changes),"gene":genes[0], "gene_function":gene_info[genes[0]]["function"], "gene_reference":gene_info[genes[0]]["DOI"]}
    return results


def main_preprocess(args):
    os.chdir(args.dir)
    conf = get_conf(args.db)
    refseq = pp.fasta(conf["ref"]).fa_dict
    refseqname = list(refseq.keys())[0]
    if not args.debug:
        if args.seqs:
            seqs = pp.fasta(args.seqs)
        else:
            pp.run_cmd("curl 'https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?q=*:*&fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&fq=%7B!tag=Flags_csv%7DFlags_csv:%22complete%22&cmd=download&sort=&dlfmt=fasta&fl=id,Definition_s,Nucleotide_seq' > temp.fa")
            pp.run_cmd("samtools faidx temp.fa")
            seqs = pp.fasta("temp.fa")


    samples = list(seqs.fa_dict.keys())
    if not args.debug:
        for sample in samples:
            fname = pp.get_random_file()
            open(fname,"w").write(">%s\n%s\n" % (sample,seqs.fa_dict[sample]))
            fasta_obj = pp.fasta(fname)
            vcf_obj = pp.vcf(fasta_obj.get_ref_variants(conf["ref"], sample))
            pp.run_cmd("rm %s" % fname)

    vcf_files = ["%s.vcf.gz" % s  for s in samples]
    vcf_csi_files = ["%s.vcf.gz.csi" % s  for s in samples]
    if not args.debug:
        pp.run_cmd("bcftools merge -0  %s | bcftools view -V indels -Oz -o merged.vcf.gz" % (" ".join(vcf_files)))
        pp.run_cmd("rm %s" % (" ".join(vcf_files)))
        pp.run_cmd("rm %s" % (" ".join(vcf_csi_files)))
        pp.run_cmd("vcf2fasta.py --vcf merged.vcf.gz --ref %s" % conf["ref"])
        pp.run_cmd("iqtree -s merged.fa -bb 1000 -nt AUTO -asr -redo")
    variant_data = get_variant_data("merged.vcf.gz",conf["ref"],conf["gff"],conf["proteins"])
    sample_data = get_sample_meta(samples,args.debug)
    with open("%s.meta.csv" % args.out,"w") as O:
        writer = csv.DictWriter(O,fieldnames=["id","country","date"])
        writer.writeheader()
        for row in sample_data:
            writer.writerow(row)

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
    if tree.get_common_ancestor(outgroup_leaf_names).name=="Node1":
        tree.set_outgroup(tree.get_common_ancestor("MN996527","MT106053"))
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


def main_position_sample(args):
    os.chdir(args.dir)
    conf = get_conf(args.db)
    refseq = pp.fasta(conf["ref"]).fa_dict
    refseqname = list(refseq.keys())[0]

    tree = ete3.Tree(args.tree, format=1)
    barcoding_sites = {}

    for l in open(args.barcode_bed):
        row = l.strip().split()
        if len(row)<5: continue
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




def main_profile(args):
    if pp.nofolder(args.dir):
        os.mkdir(args.dir)
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)


    ### Setup prefix for files ###
    files_prefix = args.dir+"/"+args.prefix

    if args.fasta:
        if args.read1 or args.read2:
            sys.stderr.write("Please use --fasta or --read1/2 but not both... Exiting!\n")
            quit()
        fasta_obj = pp.fasta(args.fasta)
        wg_vcf_obj = pp.vcf(fasta_obj.get_ref_variants(conf["ref"], prefix=args.prefix, file_prefix=files_prefix))
    else:
        if not args.read1:
            sys.stderr.write("Please provide assembly using --fasta or at least one read file using --read1... Exiting!\n")
            quit()
        ### Create bam file if fastq has been supplied ###
        if args.read1 and args.read2 and args.no_trim:
            # Paired + no trimming
            fastq_obj = pp.fastq(args.read1,args.read2)
        elif args.read1 and args.read2 and not args.no_trim:
            # Paired + trimming
            untrimmed_fastq_obj = pp.fastq(args.read1,args.read2)
            fastq_obj = untrimmed_fastq_obj.trim(files_prefix,threads=args.threads)
        elif args.read1 and not args.read2 and args.no_trim:
            # Unpaired + trimming
            fastq_obj = pp.fastq(args.read1,args.read2)
        elif args.read1 and not args.read2 and not args.no_trim:
            # Unpaired + trimming
            untrimmed_fastq_obj = pp.fastq(args.read1)
            fastq_obj = untrimmed_fastq_obj.trim(files_prefix,threads=args.threads)
        bam_obj = fastq_obj.map_to_ref(
            ref_file=conf["ref"], prefix=files_prefix,sample_name=args.prefix,
            aligner=args.mapper, platform=args.platform, threads=args.threads
        )
        wg_vcf_obj = bam_obj.call_variants(conf["ref"],args.caller)
        if not args.no_trim:
            pp.run_cmd("rm -f %s" % " ".join(fastq_obj.files))
    refseq = pp.fasta(conf["ref"]).fa_dict
    refseqname = list(refseq.keys())[0]

    mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])

    barcode = pp.barcode(mutations,conf["barcode"])
    sys.stdout.write("%s\t%s\n" % (args.prefix,";".join([d["annotation"] for d in barcode])))





def main_phylogeny(args):
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    samples = [x.replace(".vcf.gz","") for x in os.listdir(args.dir) if x[-len(".vcf.gz"):]==".vcf.gz"]
    vcf_files = ["%s/%s.vcf.gz" % (args.dir,x) for x in samples]
    for vcf in vcf_files:
        pp.index_bcf(vcf)
    args.tmp_vcf = pp.get_random_file()
    pp.run_cmd("bcftools merge -0  %s | bcftools view -V indels -Oz -o %s" % (" ".join(vcf_files),args.tmp_vcf))
    pp.index_bcf(args.tmp_vcf)
    args.sample_file = pp.get_random_file()
    args.ref = conf["ref"]
    open(args.sample_file,"w").write("\n".join(samples)+"\n")
    pp.run_cmd('cat %(sample_file)s | parallel  --bar -j %(threads)s "bcftools consensus -f %(ref)s -s {} %(tmp_vcf)s | sed \'s/^>.*/>{}/\' > {}.tmp.fasta"' % vars(args))
    pp.run_cmd('cat %s > %s.fa' % (" ".join(["%s.tmp.fasta" % s for s in samples]), args.prefix))
    pp.run_cmd('rm %s %s' % (" ".join(["%s.tmp.fasta" % s for s in samples]), args.sample_file))
    pp.run_cmd("iqtree -s %(prefix)s.fa -nt AUTO -czb -bb 1000" % vars(args))



def primer_evaluation(args):
    msa = pp.fasta(args.msa).fa_dict

    forward_results = cp.run_fuzznuc(args.msa, args.primerF, pmismatch=args.mismatch)
    reverse_results = cp.run_fuzznuc(args.msa, args.primerR, pmismatch=args.mismatch)
    probe_results = cp.run_fuzznuc(args.msa, args.probe, pmismatch=args.mismatch)

    rows = []
    print(list(forward_results.values())[0])
    for s in tqdm(msa):
        amplicon = cp.find_amplicon(forward_results[s], reverse_results[s],probe_results[s] )
        rows.append(amplicon)
    with open(args.out,"w") as O:
        writer = csv.DictWriter(O,fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

def main_parse_asr(args):
    mutations = cp.find_ancestral_mutations(args.msa,args.tree,args.states)
    with open(args.out,"w") as O:
        writer = csv.DictWriter(O,fieldnames = list(mutations[0].keys()))
        writer.writeheader()
        writer.writerows(mutations)

parser = argparse.ArgumentParser(description='Covid pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")


parser_sub = subparsers.add_parser('preprocess', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--seqs',help='First read file')
parser_sub.add_argument('--dir',default="/tmp/",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.add_argument('--out',default="covid_public",help='First read file')
parser_sub.add_argument('--debug',action="store_true",help='First read file')
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
parser_sub.add_argument('--fasta','-f',help='Fasta file')
parser_sub.add_argument('--read1','-1',help='First read file')
parser_sub.add_argument('--read2','-2',help='Second read file')
parser_sub.add_argument('--prefix','-p',help='Prefix for output files',required=True)
parser_sub.add_argument('--dir',default="covidv_profiler_results",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.add_argument('--no_trim',action="store_true",help="Don't trim files using trimmomatic")
parser_sub.add_argument('--threads','-t',default=1,help='Threads to use',type=int)
parser_sub.add_argument('--mapper',default="bwa", choices=["bwa","minimap2","bowtie2"],help="Mapping tool to use. If you are using nanopore data it will default to minimap2",type=str)
parser_sub.add_argument('--caller',default="bcftools", choices=["bcftools","gatk","freebayes"],help="Variant calling tool to use.",type=str)
parser_sub.add_argument('--platform','-m',choices=["illumina","nanopore"],default="illumina",help='NGS Platform used to generate data')
parser_sub.set_defaults(func=main_profile)

parser_sub = subparsers.add_parser('phylogeny', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--prefix','-p',help='Prefix for output files',required=True)
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.add_argument('--dir',default="covidv_profiler_results",help='First read file')
parser_sub.add_argument('--threads','-t',default=1,help='Threads to use',type=int)
parser_sub.set_defaults(func=main_phylogeny)

parser_sub = subparsers.add_parser('load_library', help='Load new library', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('prefix',type=str,help='Prefix to the library files')
parser_sub.set_defaults(func=main_load_library)


parser_sub = subparsers.add_parser('primer', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--msa',help='First read file',required=True)
parser_sub.add_argument('--primerF',help='First read file',required=True)
parser_sub.add_argument('--primerR',help='First read file',required=True)
parser_sub.add_argument('--probe',help='First read file',required=True)
parser_sub.add_argument('--out',help='First read file',required=True)
parser_sub.add_argument('--mismatch',type=int,default=3,help='First read file')
parser_sub.set_defaults(func=primer_evaluation)

parser_sub = subparsers.add_parser('asr', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--msa',help='First read file',required=True)
parser_sub.add_argument('--tree',help='First read file',required=True)
parser_sub.add_argument('--states',help='First read file',required=True)
parser_sub.add_argument('--out',default="out.csv",help='First read file')
parser_sub.set_defaults(func=main_parse_asr)

args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
