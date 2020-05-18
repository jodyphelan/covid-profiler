import sys
import pathogenprofiler as pp
import ete3
from collections import defaultdict
from tqdm import tqdm


def find_ancestral_mutations(msa_file,tree_file,states_file):
    seqs = pp.fasta(msa_file).fa_dict

    tree = ete3.Tree(tree_file,format=1)
    node_names = set([tree.name] + [n.name.split("/")[0] for n in tree.get_descendants()])
    leaf_names = set(tree.get_leaf_names())
    internal_node_names = node_names - leaf_names


    for n in tree.traverse():
        if n.name.split("/")[0] in node_names:
            if "Node" in n.name:
                tmp = n.name.split("/")
                n.name = tmp[0]
                if len(tmp)>1:
                    n.support = tmp[1]


    tree.write(format=1, outfile=tree_file+".reformatted.tree")

    states = defaultdict(dict)
    sites = set()
    sys.stderr.write("Loading states\n")
    for l in tqdm(open(states_file)):
        if l[0]=="#": continue
        row = l.strip().split()
        if row[0]=="Node": continue
        if row[0] not in internal_node_names: continue
        states[int(row[1])][row[0]] = row[2]
        sites.add(int(row[1]))

    sys.stderr.write("Loading alignment sites\n")
    for site in tqdm(sites):
        for sample in seqs:
            states[site][sample] = seqs[sample][site-1]

    barcoding_sites = []
    convergent_sites = []
    mutations = []
    for site in tqdm(sites):
        nucleotides = set([states[site][n] for n in node_names])
        if len(nucleotides)==1: continue

        # Set up storage objjects
        origins = []
        internal_node_change = False

        tree.add_feature("state",states[site][tree.name])
        for n in tree.traverse():
            if n == tree: continue
            node_state = states[site][n.name]
            if node_state!=n.get_ancestors()[0].state:
                # if site==14408: import pdb; pdb.set_trace()
                origins.append(n.name)
                if n.name in internal_node_names:
                     internal_node_change = True
            n.add_feature("state",node_state)

        type = "unique"
        if internal_node_change and len(origins)==1:
            type = "barcoding"
            barcoding_sites.append(site)
        if len(origins)>1:
            type = "convergent"
            convergent_sites.append(site)

        tmp_data = {
            "position": site,
            "mutation_type":type,
            "origins": len(origins),
            "branches":",".join(origins),
            }
        for sample in leaf_names:
            tmp_data[sample] = states[site][sample]
        mutations.append(tmp_data)

    return mutations
