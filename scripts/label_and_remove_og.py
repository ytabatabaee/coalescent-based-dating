import dendropy
import argparse
import re
import math
import numpy as np
from itertools import combinations


def label_tree(args):
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(path=args.tree1, schema='newick', taxon_namespace=tns, rooting='force-rooted')
    t2 = dendropy.Tree.get(path=args.tree2, schema='newick', taxon_namespace=tns, rooting='force-rooted')

    #root = t.find_node_with_taxon_label('0').parent_node
    #t.prune_taxa_with_labels(["0"])
    #t.prune_taxa([root])

    for node in t2.postorder_node_iter():
        if node.is_leaf():
            node.value = node.taxon.label
        else:
            left = node._child_nodes[0]
            right = node._child_nodes[1]
            node.value = left.value + ',' + right.value
            tl = node.value.split(',')
            print(tl)
            mrca = t1.mrca(taxon_labels=tl)
            print(mrca.label)
            node.taxon = dendropy.Taxon(mrca.label)


    internal_cnt = 0
    #for node in t.preorder_node_iter():
    #    if node.taxon is None:
    #        internal_cnt += 1
    #        node.taxon = dendropy.Taxon('i' + str(internal_cnt))

    root = t2.find_node_with_taxon_label('i1')
    root.edge.length = 0

    with open(args.output, 'w') as f:
        f.write(str(t2) + ';\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Label internal nodes")
    parser.add_argument("-t1", "--tree1", type=str,  required=True,
                        help="True tree in newick format")
    parser.add_argument("-t2", "--tree2", type=str, required=True,
                        help="Other tree in newick format")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Labeled tree without outgroup in newick format")
    args = parser.parse_args()
    label_tree(args)
