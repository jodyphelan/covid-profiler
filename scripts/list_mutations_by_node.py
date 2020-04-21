import csv
import sys
import argparse


def main(args):
    fields = ['position', 'mutation_type', 'origins', 'branches', 'alts', 'types', 'changes', 'gene', 'gene_function', 'gene_reference']
    for node in args.nodes.split(","):
        for row in csv.DictReader(open(args.csv)):
            nodes = row["branches"].split(",")
            # print(nodes)
            if node in nodes:
                sys.stdout.write("%s\t%s\n" % (node,"\t".join([row[x] for x in fields])))



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--csv',help='VCF file',required=True)
parser.add_argument('--nodes',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
