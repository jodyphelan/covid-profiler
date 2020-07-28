import pathogenprofiler as pp
import argparse
from tqdm import tqdm
from collections import Counter
import sys

def main(args):
    acgtn = set(["A","C","G","T","N"])
    sys.stderr.write("Loading sequences\n")
    seqs = pp.fasta(args.fasta).fa_dict
    sys.stderr.write("Masking sequences\n")
    for seq in tqdm(seqs):
        nucs = Counter(list(seqs[seq]))
        for nuc in nucs:
            if nuc not in acgtn:
                seqs[seq] = seqs[seq].replace(nuc,"N")

    sys.stderr.write("Writing sequences\n")
    with open(args.out,"w") as O:
        for seq in seqs:
            O.write(">%s\n%s\n" % (seq,"".join(seqs[seq])))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',help='VCF file',required=True)
parser.add_argument('--out',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
