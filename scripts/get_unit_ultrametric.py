import dendropy
import argparse
import re
import math
import numpy as np
from itertools import combinations


def convert_to_unit_tree(args):
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(path=args.tree, schema='newick', taxon_namespace=tns, rooting='force-rooted')
    t2 = dendropy.Tree.get(path=args.tree, schema='newick', taxon_namespace=tns, rooting='force-rooted')
    tree_height = 0
    for node in t1.postorder_node_iter():
        if node.taxon is not None:
            continue
        elif node.edge.length:
            max_child_length = 0
            for child in node._child_nodes:
                if child.edge.length > max_child_length: 
                    max_child_length = child.edge.length
            print(node.edge.length)
            node.edge.length = node.edge.length + max_child_length
            tree_height = max(tree_height, node.edge.length)

    print('tree height:', tree_height)
    for node in t2.postorder_node_iter():
        if node.edge.length:
            node.edge.length = node.edge.length / tree_height

    with open(args.output, 'w') as f:
        f.write(str(t2) + ';\n')
        print('\nUnit-height tree written to', args.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert tree to a tree with unit height")
    parser.add_argument("-t", "--tree", type=str,  required=True,
                        help="Ultrametric tree in newick format")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Unit-height tree in newick format")
    args = parser.parse_args()
    convert_to_unit_tree(args)
