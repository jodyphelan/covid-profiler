import sys
import argparse
import pathogenprofiler as pp
import ete3
import sqlite3
import csv
def main(args):
    # pp.run_cmd("covid-profiler.py preprocess")
    tree_text = open("%s/covid_public.tree" % args.dir).readline().strip()
    tree = ete3.Tree("%s/covid_public.tree" % args.dir,format=1)
    samples = tree.get_leaf_names()

    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    c.execute("INSERT INTO tree (newick) VALUES (?)",(tree_text,))

    for row in csv.DictReader(open("%s/covid_public.meta.csv" % args.dir)):
        try:
            c.execute("INSERT INTO tree_data (id, country, collection_date) VALUES (?,?,?)", (row["id"],row["country"],row["date"],))
        except:
            pass

    c.execute("DROP TABLE IF EXISTS mutations")
    c.execute("CREATE TABLE mutations (position INT, mutation_type TEXT, gene TEXT, alts TEXT, functional_types TEXT, changes TEXT, origins INT, branches TEXT, %s )" % (",".join(["%s TEXT" % s for s in samples])))

    for row in csv.DictReader(open("%s/covid_public.mutation_summary.csv" % args.dir)):
        c.execute("INSERT INTO mutations (position, mutation_type, gene, alts, functional_types, changes, origins, branches) VALUES (%(position)s,'%(mutation_type)s','%(gene)s','%(alts)s','%(functional_types)s','%(changes)s',%(origins)s,'%(branches)s')" % row)
        for s in samples:
            c.execute("UPDATE mutations SET %s = '%s' where position = %s" % (s,row[s],row["position"]))
    conn.commit()

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--db',help='Number of threads for parallel operations',required=True)
parser.add_argument('--dir',default="/tmp/",help='VCF file')

parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
