#!/usr/bin/env python
# modified from https://github.com/kwongj/nw_multi2bifurcation/tree/master

# Usage
import argparse
from ete3 import Tree
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputtree", type=str,  required=True,
                        help="Input unresolved tree in newick format")
parser.add_argument("-o", "--outputtree", type=str, required=True,
                        help="Output resolved tree in newick format")
parser.add_argument('--nosupport', action='store_true', help='Remove internal node support values from tree')
args = parser.parse_args()


t = Tree(args.inputtree)
t.resolve_polytomy(recursive=True)
with open(args.outputtree, 'w') as f:
        f.write(str(t.write()))


# Remove support values from tree
if args.nosupport:
	tn = re.sub(r'\)[0-9\.]+:', r'):', tn)
