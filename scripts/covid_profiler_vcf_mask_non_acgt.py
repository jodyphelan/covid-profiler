import sys
import argparse


def main(args):
    acgt = set(["A","C","G","T"])
    for l in sys.stdin:
        if l[0]=="#":
            sys.stdout.write(l)
            continue
        row = l.strip().split()
        alleles = [row[3]] + row[4].split(",")
        gt_str = "".join(row[9:])
        for i,a in enumerate(alleles):
            if a not in acgt:
                gt_str = gt_str.replace(str(i),".")
        gt_str = "\t".join(list(gt_str))
        sys.stdout.write("%s\t%s\n" % ("\t".join(row[:9]),gt_str))



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
