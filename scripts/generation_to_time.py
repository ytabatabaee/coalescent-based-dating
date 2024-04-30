import dendropy
import argparse
import re
import math
import numpy as np
from itertools import combinations


def convert_to_time_tree(args):
    tns = dendropy.TaxonNamespace()
    t = dendropy.Tree.get(path=args.tree, schema='newick', taxon_namespace=tns, rooting='force-rooted')

    root = t.find_node_with_taxon_label('0').parent_node
    t.prune_taxa_with_labels(["0"])
    t.prune_taxa([root])

    for node in t.postorder_node_iter():
        if node.edge.length:
            node.edge.length = node.edge.length * args.generationtime / 1000000

    internal_cnt = 0
    for node in t.preorder_node_iter():
        if node.taxon is None:
            internal_cnt += 1
            node.taxon = dendropy.Taxon('i' + str(internal_cnt))

    with open(args.output, 'w') as f:
        f.write(str(t) + ';\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert generation-unit tree to time-unit tree")
    parser.add_argument("-t", "--tree", type=str,  required=True,
                        help="Generation unit tree in newick format")
    parser.add_argument("-gt", "--generationtime", type=int, required=True,
                        help="Generation time (in years)")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Time tree in newick format")
    args = parser.parse_args()
    convert_to_time_tree(args)
