import sys
import argparse
import pandas as pd
import datetime as dt
from uuid import uuid4
import  subprocess as sp

def main(args):
    df = pd.read_csv(args.meta)
    df.date = pd.to_datetime(df.date, infer_datetime_format=True)
    t = dt.datetime.now() - dt.timedelta(days=args.days)
    df = df[df.date > t]
    df = df.groupby("iso_a3",group_keys=False).apply(lambda x: x.sample(min(len(x),args.max_country_seqs)))
    tmp_file = str(uuid4())
    open(tmp_file,"w").write("\n".join(list(df.id)))
    df.to_csv("%s.csv" % args.out)
    sp.call(f"seqtk subseq {args.fasta} {tmp_file} > {args.out}.fasta",shell=True)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df.groupby("iso_a3").apply(lambda x:len(x)))


parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fasta',type=str,help='CSV file')
parser.add_argument('--meta',type=str,help='CSV file')
parser.add_argument('--out',type=str,help='Output file with sequence names')
parser.add_argument('--max_country_seqs',type=int,default=5000,help='Output file with sequence names')
parser.add_argument('--days',type=int,default=90,help='Output file with sequence names')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)