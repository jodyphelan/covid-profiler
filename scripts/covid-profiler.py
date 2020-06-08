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
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","version":".version.json","proteins":".proteins.csv","msa":".msa.fa"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
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
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","version":".version.json","proteins":".proteins.csv"}
    if pp.nofolder(sys.base_prefix+"/share/covidprofiler"):
        pp.run_cmd("mkdir %s " % (sys.base_prefix+"/share/covidprofiler/"))
    for key in files:
        new_file_location = sys.base_prefix+"/share/covidprofiler/"+lib_prefix+files[key]
        pp.run_cmd("cp %s %s" % (args.prefix+files[key],new_file_location))
    pp.run_cmd("samtools faidx %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".fasta")
    pp.run_cmd("bwa index %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".fasta")
    if os.path.isfile("%s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".dict"):
        pp.run_cmd("rm %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".dict")
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

    results = {}
    barcode_mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])
    barcode = pp.barcode(barcode_mutations,conf["barcode"])
    clade = ";".join(sorted([d["annotation"] for d in barcode]))
    sys.stdout.write("%s\t%s\n" % (args.prefix,clade))
    results["clade"] = clade

    variant_data = cp.get_variant_data(wg_vcf_obj.filename,conf["ref"],conf["gff"],conf["proteins"])
    results["variants"] = variant_data
    json.dump(results,open("%s.results.json" % files_prefix,"w"))

def main_collate(args):
    samples = [x.replace(args.suffix,"") for x in os.listdir(args.dir) if x[-len(args.suffix):]==args.suffix]
    lineages = {}
    sys.stderr.write("Loading data\n")
    for s in tqdm(samples):
        data = json.load(open(pp.filecheck("%s/%s%s" % (args.dir,s,args.suffix))))
        lineages[s] = data["clade"]
    with open(args.out+".clades.csv","w") as O:
        O.write("isolate,clade\n")
        for s in samples:
            O.write("%s,%s\n" % (s,lineages[s]))


def primer_evaluation(args):
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    msa_obj = pp.fasta(conf["msa"]).fa_dict
    forward_results = cp.run_fuzznuc(conf["msa"], args.primerF, pmismatch=args.mismatch)
    reverse_results = cp.run_fuzznuc(conf["msa"], args.primerR, pmismatch=args.mismatch)
    probe_results = cp.run_fuzznuc(conf["msa"], args.probe, pmismatch=args.mismatch)

    rows = []
    for s in tqdm(msa_obj):
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
parser_sub.add_argument('--dir',default="covid_profiler_results",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.add_argument('--no_trim',action="store_true",help="Don't trim files using trimmomatic")
parser_sub.add_argument('--threads','-t',default=1,help='Threads to use',type=int)
parser_sub.add_argument('--mapper',default="bwa", choices=["bwa","minimap2","bowtie2"],help="Mapping tool to use. If you are using nanopore data it will default to minimap2",type=str)
parser_sub.add_argument('--caller',default="bcftools", choices=["bcftools","gatk","freebayes"],help="Variant calling tool to use.",type=str)
parser_sub.add_argument('--platform','-m',choices=["illumina","nanopore"],default="illumina",help='NGS Platform used to generate data')
parser_sub.set_defaults(func=main_profile)

parser_sub = subparsers.add_parser('collate', help='Load new library', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--dir',type=str,default="covid_profiler_results",help='Prefix to the library files')
parser_sub.add_argument('--suffix',type=str,default=".results.json", help='Prefix to the library files')
parser_sub.add_argument('--out',type=str,default="covid_profiler", help='Prefix to the library files')
parser_sub.set_defaults(func=main_collate)

parser_sub = subparsers.add_parser('load_library', help='Load new library', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('prefix',type=str,help='Prefix to the library files')
parser_sub.set_defaults(func=main_load_library)


parser_sub = subparsers.add_parser('primer', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--primerF',help='First read file',required=True)
parser_sub.add_argument('--primerR',help='First read file',required=True)
parser_sub.add_argument('--probe',help='First read file',required=True)
parser_sub.add_argument('--out',help='First read file',required=True)
parser_sub.add_argument('--mismatch',type=int,default=3,help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
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
