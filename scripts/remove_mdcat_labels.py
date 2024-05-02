import dendropy
import argparse
import re
import math
import numpy as np
from itertools import combinations


def convert_to_time_tree(args):
    tns = dendropy.TaxonNamespace()
    t = dendropy.Tree.get(path=args.tree, schema='newick', taxon_namespace=tns, rooting='force-rooted')

    for node in t.preorder_node_iter():
        if node.label:
            node.label = re.sub("[\(\[].*?[\)\]]", "", node.label)
        if node.taxon:
            node.taxon.label = re.sub("[\(\[].*?[\)\]]", "", node.taxon.label)

    with open(args.output, 'w') as f:
        f.write(str(t) + ';\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove mdcat labels")
    parser.add_argument("-t", "--tree", type=str,  required=True,
                        help="MDCAT tree in newick format")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Time tree in newick format")
    args = parser.parse_args()
    convert_to_time_tree(args)
