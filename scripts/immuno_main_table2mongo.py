import json
from pymongo import MongoClient
import csv
from collections import Counter
import argparse
from tqdm import tqdm


def main(args):
    client = MongoClient()
    db = client.test_database
    db.drop_collection("immuno")
    immuno_db = db.immuno

    for row in tqdm(csv.DictReader(open(args.csv))):
        row["Protein AA Position"] = int(row["Protein AA Position"])
        id = immuno_db.insert_one(row)



parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--csv',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
