import sys
import argparse
import pathogenprofiler as pp
import ete3
import sqlite3

def main(args):
    # pp.run_cmd("covid-profiler.py preprocess")
    tree_text = open("%s/merged.fa.treefile" % args.dir).readline().strip()
    tree = ete3.Tree("%s/merged.fa.treefile" % args.dir,format=1)
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    c.execute("INSERT INTO tree (newick) VALUES (?)",(tree_text,))
    conn.commit()

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--db',help='Number of threads for parallel operations',required=True)
parser.add_argument('--dir',default="/tmp/",help='VCF file')

parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
