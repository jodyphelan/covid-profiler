import json
from pymongo import MongoClient
import csv
from collections import Counter
import argparse
from tqdm import tqdm


def main(args):
    client = MongoClient()
    db = client.test_database
    db.drop_collection("mutations")
    db.drop_collection("tree")
    db.drop_collection("meta")
    mutation_db = db.mutations

    for row in tqdm(csv.DictReader(open(args.asr))):
        genotypes = {"_".join(x.split("_")[:3]):row[x] for x in row.keys() if x[:3]=="EPI"}
        allele_counts = Counter(list(genotypes.values()))
        major_allele = allele_counts.most_common(1)[0][0]
        minor_alleles_dict = {s:genotypes[s] for s in genotypes if genotypes[s]!=major_allele}
        results = {c.replace(".","_"):row[c] for c in row.keys() if c[:3]!="EPI"}
        results["major_allele"] = major_allele
        results["minor_allele_genotypes"] = minor_alleles_dict
        results["_id"] = row["position"]
        id = mutation_db.insert_one(results).inserted_id

    db.tree.insert_one({"tree": open(args.tree).readline().strip()})
    meta_db = db.meta
    meta = {}
    for row in csv.DictReader(open(args.meta)):
        meta[row["id"]] = {}
        meta[row["id"]]["id"] = row["id"]
        meta[row["id"]]["country"] = row["country"]
        meta[row["id"]]["date"] = row["date"]
    meta_db.insert_one(meta)


parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--asr',help='VCF file',required=True)
parser.add_argument('--tree',help='VCF file',required=True)
parser.add_argument('--meta',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
