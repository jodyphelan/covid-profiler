from Bio.PDB import *
from collections import defaultdict
import argparse
import json


def main(args):
    p = PDBParser()
    chains = json.loads(args.chain_json)
    structure = p.get_structure('X', args.pdb)
    residues = defaultdict(list)
    for model in structure:
        for chain in model:
            for residue in chain:
                # import pdb; pdb.set_trace()
                if chain.id in chains:
                    residues[chains[chain.id]].append(residue._id[1])

    json.dump({"mapping":chains,"residues":residues},open(args.pdb+".available_residues.json","w"))

parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--pdb',help='VCF file',required=True)
parser.add_argument('--chain-json',help='VCF file',required=True)
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
