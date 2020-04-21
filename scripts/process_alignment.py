import pathogenprofiler as pp
import argparse
from collections import Counter

refseq = "NC_045512.2"
def main(args):
    seqs = pp.fasta(args.msa).fa_dict
    ref_pos = 0
    print("ref_pos\talignment_pos\tA\tC\tG\tT\tN\tgap")
    for i in range(len(list(seqs.values())[0])):
        alignment_pos = i+1
        if seqs[refseq][i]!="-":
            ref_pos+=1
        allele_count = Counter([seqs[s][i].upper() for s in seqs])
        num_N = allele_count["N"] if "N" in allele_count else 0
        num_gap = allele_count["-"] if "-" in allele_count else 0
        num_A = allele_count["A"] if "A" in allele_count else 0
        num_C = allele_count["C"] if "C" in allele_count else 0
        num_G = allele_count["G"] if "G" in allele_count else 0
        num_T = allele_count["T"] if "T" in allele_count else 0


        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (ref_pos,alignment_pos,num_A,num_C,num_G,num_T,num_N,num_gap))





parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--msa',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
