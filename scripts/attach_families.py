import argparse
import dendropy
import os
import sys
import pandas as pd


def main(args):
    tax = dendropy.TaxonNamespace()
    tr = dendropy.Tree.get(path=args.tree,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
    #annotations = pd.read_csv(args.annotation,sep=',')#sep='\t')
    annotations = pd.read_csv(args.annotation,sep='\t')
    #ant_dict = annotations.set_index('tipnamecodes')['howardmoore.order'].to_dict() # change family to order
    ant_dict = annotations.set_index('taxa')['family'].to_dict() # change family to order
    for node in tr.preorder_node_iter():
        if node.is_leaf():
            node.taxon.label = node.taxon.label.replace(" ", "_")
            node.taxon.label = str(ant_dict[node.taxon.label]) + '_' + node.taxon.label

    with open(args.tree+'.family', 'w') as f: # change family to order
        f.write(str(tr) + ';\n')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="annotate")
    parser.add_argument("-ant", "--annotation", type=str,  required=True,
                        help="Annotation file for taxa")
    parser.add_argument("-t", "--tree", type=str, required=True,
                        help="File containing newick string for tree")
    main(parser.parse_args())
