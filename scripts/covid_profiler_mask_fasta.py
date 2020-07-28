import pathogenprofiler as pp
import argparse
from tqdm import tqdm


def main(args):
    seqs = pp.fasta(args.fasta).fa_dict
    for seq in seqs:
        seqs[seq] = list(seqs[seq])
    for l in open(args.bed):
        row = l.strip().split()
        for i in tqdm(range(int(row[1])-1,int(row[2]))):
            for seq in seqs:
                seqs[seq][i] = "N"

    with open(args.out,"w") as O:
        for seq in seqs:
            O.write(">%s\n%s\n" % (seq,"".join(seqs[seq])))


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',help='VCF file',required=True)
parser.add_argument('--bed',help='VCF file',required=True)
parser.add_argument('--out',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
