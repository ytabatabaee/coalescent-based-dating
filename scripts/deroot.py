import dendropy
import sys
import numpy as np
import argparse


def main(args):
    tns = dendropy.TaxonNamespace()
    st = dendropy.Tree.get(path=args.tree, schema='newick', taxon_namespace=tns)
    st.deroot()
    with open(args.tree + '.derooted', 'w') as f:
        f.write(str(st)+';\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="deroot a tree")
    parser.add_argument("-t", "--tree", type=str,  required=True,
                        help="tree file in newick format")
    main(parser.parse_args())
