import dendropy
import argparse
import re
import math
import random
import numpy as np
from itertools import combinations


def convert_to_time_tree(args):
    tns = dendropy.TaxonNamespace()
    t = dendropy.Tree.get(path=args.tree, schema='newick', taxon_namespace=tns, rooting='force-rooted')

    #root = t.find_node_with_taxon_label('0').parent_node
    #t.prune_taxa_with_labels(["0"])
    #t.prune_taxa([root])


    internal_cnt = 0
    for node in t.preorder_node_iter():
        if node.parent_node:
            node.value = node.parent_node.value + 1 # num ancestors
            node.age = node.edge.length + node.parent_node.age # time
        else: # root node
            node.age = 0
            node.value = 0
        if node.is_leaf():
            node.label = node.taxon.label
        print(node.value)
       # if node.taxon is None:
        #    internal_cnt += 1
        #    node.taxon = dendropy.Taxon('i' + str(internal_cnt))
    
    tree_height = t.find_node_with_taxon_label('1').age

    # backward time (present leafs 0, root oldest)
    for node in t.preorder_node_iter():
        node.age = tree_height - node.age

    calibration_points = []
    root = t.find_node_with_label('i1')
    root.edge.length = 0
    calibration_points.append((root.label, root.age))

    rem_calibrations = args.numcalibrations - 1
    selected_taxa = []
    # select taxa with probability 1/#ancestors
    while rem_calibrations > 0:
        for node in t.leaf_node_iter():
            if rem_calibrations > 0 and random.random() < 1 / node.value and \
                    node.parent_node.label != root.label: # select leaf with probability 1/#ancestors
                rem_calibrations -= 1
                selected_taxa.append(node.label)

    print('selected taxa', selected_taxa)

    for taxa in selected_taxa:
        node = t.find_node_with_taxon_label(taxa)
        flag = True
        while flag:
            parent = node
            k = random.randint(1, node.value)
            x = 0
            while x < k:
                parent = node.parent_node
                x += 1
            if not parent.label in [x[0] for x in calibration_points]:
                flag = False
        calibration_points.append((parent.label, parent.age))

    print('calibration nodes:', calibration_points)

    with open(args.output + '.txt', 'w') as f:
        for cp in calibration_points:
            f.write(cp[0] + '\t' + str(round(cp[1], 3)) + '\n')
    #with open(args.output + '.trees', 'w') as f:
    #    f.write(str(t) + ';\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate time trees with calibration points")
    parser.add_argument("-t", "--tree", type=str,  required=True,
                        help="Tree in time units")
    parser.add_argument("-n", "--numcalibrations", type=int, required=True,
                        help="Number of calibration points")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output name")
    args = parser.parse_args()
    convert_to_time_tree(args)
