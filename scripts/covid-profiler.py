from uuid import uuid4
import pathogenprofiler as pp
import argparse
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
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","version":".version.json","proteins":".proteins.csv","msa":".msa.fa","non_coding_bed":".non_coding.bed"}
    conf = {}
    for key in files:
        sys.stderr.write("Using %s file: %s\n" % (key,library_prefix+files[key]))
        conf[key] = pp.filecheck(library_prefix+files[key])
    return conf



def main_load_library(args):
    lib_prefix = args.prefix.split("/")[-1]
    files = {"gff":".gff","ref":".fasta","barcode":".barcode.bed","version":".version.json","proteins":".proteins.csv","non_coding_bed":".non_coding.bed"}
    if pp.nofolder(sys.base_prefix+"/share/covidprofiler"):
        pp.run_cmd("mkdir %s " % (sys.base_prefix+"/share/covidprofiler/"))
    pp.run_cmd("cp %s %s" % (args.msa,"%s/share/covidprofiler/%s.msa.fa" % (sys.base_prefix,lib_prefix)))
    pp.run_cmd("cp %s %s" % (args.meta,"%s/share/covidprofiler/%s.msa.meta.csv" % (sys.base_prefix,lib_prefix)))
    for key in files:
        new_file_location = sys.base_prefix+"/share/covidprofiler/"+lib_prefix+files[key]
        pp.run_cmd("cp %s %s" % (args.prefix+files[key],new_file_location))
    pp.run_cmd("samtools faidx %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".fasta")
    pp.run_cmd("bwa index %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".fasta")
    if os.path.isfile("%s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".dict"):
        pp.run_cmd("rm %s" % sys.base_prefix+"/share/covidprofiler/"+lib_prefix+".dict")
    pp.log("Sucessfully imported library")


def get_pangolin_lineage(fasta_file):
    tmpfile = str(uuid4())
    pp.run_cmd(f"pangolin --outfile {tmpfile} {fasta_file}")
    for row in csv.DictReader(open(tmpfile)):
        return row["lineage"]


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
        wg_vcf_obj = bam_obj.call_variants(conf["ref"],args.caller,remove_missing=True)
        cp.vcf2consensus(bam_obj.bam_file,wg_vcf_obj.filename,conf["ref"],wg_vcf_obj.samples[0],wg_vcf_obj.prefix+".consensus.fasta")
        if not args.no_trim:
            pp.run_cmd("rm -f %s" % " ".join(fastq_obj.files))
    refseq = pp.fasta(conf["ref"]).fa_dict
    refseqname = list(refseq.keys())[0]

    results = {}
    if not args.fasta:
        for l in pp.cmd_out("bedtools genomecov -ibam %s | datamash median 3" % (bam_obj.bam_file)):
            results["mean_depth"] = int(l.strip())
    # barcode_mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])
    # barcode = pp.barcode(barcode_mutations,conf["barcode"])
    # clade = ";".join(sorted([d["annotation"] for d in barcode]))
    results["clade"] = get_pangolin_lineage(args.fasta if args.fasta else wg_vcf_obj.prefix+".consensus.fasta")

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

def main_aln(args):
    """
    mafft --auto --thread -1 --keeplength --addfragments gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.fasta ~/covid/cvdb.fasta > gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.aln
python ~/gisaid_scripts/get_fasta_stats.py  --fasta gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.aln  --bed ~/covid/static/coding.bed --out gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.stats
awk '$3<=2.5 && $4<=3 && $5<=50' gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.aln.stats | cut -f1 > seq_filtered_samples.txt
seqtk subseq gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.aln seq_filtered_samples.txt > gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.aln
python ~/gisaid_scripts/mask_fasta.py  --fasta gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.aln --bed ~/covid/static/non_coding_mask.bed --out gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.bed_masked.aln
python ~/gisaid_scripts/mask_fasta_non_acgt.py --fasta gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.bed_masked.aln --out gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.bed_masked.acgt.aln
snp-sites -v gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.bed_masked.acgt.aln  | python ~/gisaid_scripts/vcf_fix_ref.py --ref ~/covid/cvdb.fasta | python ~/gisaid_scripts/vcf_mask_non_acgt.py  | tqdm | bcftools view -a -Oz -o gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.bed_masked.acgt.vcf.gz
bcftools norm --threads 4 -m - gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.bed_masked.acgt.vcf.gz -Oz -o gisaid_hcov-19_2020_07_22_09.meta_filtered.filtered.site_filtered.bed_masked.acgt.multi_split.vcf.gz

    """
    conf = get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    pp.run_cmd("mafft --auto --thread %s --keeplength --addfragments %s %s  > %s.aln" % (args.threads,args.fasta,conf["ref"],args.prefix))
    pp.run_cmd("covid_profiler_mask_fasta.py  --fasta %s.aln --bed %s --out %s.bed_masked.aln" % (args.prefix,conf["non_coding_bed"],args.prefix))
    pp.run_cmd("covid_profiler_mask_fasta_non_acgt.py --fasta %s.bed_masked.aln --out %s.bed_masked.acgt.aln" % (args.prefix,args.prefix))
    pp.run_cmd("iqtree -m GTR+F+R2 -s %s.bed_masked.acgt.aln -nt %s -czb -pre %s" % (args.prefix,args.threads,args.prefix))


parser = argparse.ArgumentParser(description='Covid pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")


parser_sub = subparsers.add_parser('profile', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--fasta','-f',help='Fasta file')
parser_sub.add_argument('--read1','-1',help='First read file')
parser_sub.add_argument('--read2','-2',help='Second read file')
parser_sub.add_argument('--prefix','-p',help='Prefix for output files',required=True)
parser_sub.add_argument('--dir',default="covid_profiler_results",help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.add_argument('--no-trim',action="store_true",help="Don't trim files using trimmomatic")
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
parser_sub.add_argument('--prefix',type=str,help='Prefix to the library files',required=True)
parser_sub.add_argument('--msa',type=str,help='Prefix to the MSA file',required=True)
parser_sub.add_argument('--meta',type=str,help='Prefix to the meta file',required=True)
parser_sub.set_defaults(func=main_load_library)


parser_sub = subparsers.add_parser('primer', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--primerF',help='First read file',required=True)
parser_sub.add_argument('--primerR',help='First read file',required=True)
parser_sub.add_argument('--probe',help='First read file',required=True)
parser_sub.add_argument('--out',help='First read file',required=True)
parser_sub.add_argument('--mismatch',type=int,default=3,help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.set_defaults(func=primer_evaluation)


parser_sub = subparsers.add_parser('aln', help='Output program version and exit', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--fasta',help='First read file',required=True)
parser_sub.add_argument('--prefix',help='First read file',required=True)
parser_sub.add_argument('--threads',default=1,help='First read file')
parser_sub.add_argument('--db',default="cvdb",help='First read file')
parser_sub.set_defaults(func=main_aln)

args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
