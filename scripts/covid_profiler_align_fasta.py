import sys
import argparse
import pathogenprofiler as pp
import covid_profiler
import pyfastx
from uuid import uuid4
from tqdm import tqdm
import csv

def main(args):

    args.uuid = str(uuid4())
    conf = covid_profiler.get_conf_dict(sys.base_prefix+"/share/covidprofiler/%s" % args.db)
    vars(args).update(conf)

    args.final_aln = "%s/%s.aln" % (args.working_dir,args.out)
    args.final_vcf = "%s/%s.vcf.gz" % (args.working_dir,args.out)
    args.final_csv = "%s/%s.variant_info.csv" % (args.working_dir,args.out)

    for name,seq in pyfastx.Fasta(args.ref,build_index=False):
        ref_seq = seq
    pp.run_cmd("mafft --auto --thread -1 --keeplength --addfragments %(fasta)s %(ref)s > %(working_dir)s/%(uuid)s.aln" % vars(args))

    troublesome_sites = set()
    if args.mask_troublesome_sites:
        from urllib.request import urlopen
        with urlopen('https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf') as response:
            for l in response:
                row = l.decode().strip().split()
                if row[0][0]=="#": continue
                troublesome_sites.add(int(row[1]))
    print(troublesome_sites)
    with open(args.final_aln,"w") as O:
        for entry in tqdm(pyfastx.Fasta("%(working_dir)s/%(uuid)s.aln" % vars(args), full_name=True)):
            masked_seq = list(entry.seq.upper())
            for start,end in [(1,265),(29675,29903)]:
                for i in range(start-1,end):
                    masked_seq[i] = "N"
            for pos in troublesome_sites:
                masked_seq[pos-1] = "N"
            acgt = set(["A","C","G","T"])
            for pos in [i for i,n in enumerate(masked_seq) if n not in acgt]:
                masked_seq[pos] = "N"
            O.write(">%s\n%s\n" % (entry.name,"".join(masked_seq)))


    pp.run_cmd("snp-sites -v %(final_aln)s  | covid_profiler_vcf_fix_ref.py --ref %(ref)s | covid_profiler_vcf_mask_non_acgt.py  | tqdm | bcftools view -a -Oz -o %(uuid)s.bed_masked.vcf.gz" % vars(args))
    pp.run_cmd("bcftools stats -s - %(uuid)s.bed_masked.vcf.gz > %(uuid)s.bed_masked.vcf.gz.stats" % vars(args))
    pp.run_cmd("bcftools norm -m -  %(uuid)s.bed_masked.vcf.gz -Oz -o %(final_vcf)s" % vars(args))

    variant_data = covid_profiler.get_variant_data(args.final_vcf,conf["ref"],conf["gff"],conf["proteins"])
    with open(args.final_csv,"w") as O:
        fieldnames = list(variant_data[0].keys())
        writer = csv.DictWriter(O,fieldnames)
        writer.writeheader()
        writer.writerows(variant_data)

    pp.run_cmd("iqtree -s %(final_aln)s -m GTR+F+G4 -nt 3" % vars(args))
    sys.stderr.write("\n\n----------------\n")
    sys.stderr.write("Program complete\n")
    sys.stderr.write("----------------\n")
    sys.stderr.write("Alignment: %s\n" % args.final_aln)
    sys.stderr.write("VCF: %s\n" % args.final_vcf)
    sys.stderr.write("Variant summary csv: %s\n" % args.final_csv)

parser = argparse.ArgumentParser(description='Covid alignments',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',type=str,help='File with samples',required=True)
parser.add_argument('--out',type=str,help='File with samples',required=True)
parser.add_argument('--working-dir',default=".",type=str,help='File with samples')
parser.add_argument('--db',default="cvdb",type=str,help='File with samples')
parser.add_argument('--mask-troublesome-sites',action='store_true',help="Mask sites in alignment reported here: https://github.com/W-L/ProblematicSites_SARS-CoV2")
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
