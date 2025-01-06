"""
This code is modified from the following script written by Erin K. Molloy in NJMerge software
https://github.com/ekmolloy/njmerge/blob/master/python/compare_trees.py

All rights reserved.
License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import argparse
import dendropy


def clade_distance(tr1, tr2):
    from dendropy.calculate.treecompare \
        import symmetric_difference

    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])

    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)

        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)

        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)

    tr1.update_bipartitions()
    tr2.update_bipartitions()

    nl = len(com)
    cd = symmetric_difference(tr1, tr2)# / (2*nl - 4)

    return nl, cd


def main(args):
    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(path=args.tree1,
                            schema='newick',
                            rooting='force-rooted',
                            taxon_namespace=tax)
    tr2 = dendropy.Tree.get(path=args.tree2,
                            schema='newick',
                            rooting='force-rooted',
                            taxon_namespace=tax)

    nl, cd = clade_distance(tr1, tr2)
    print("Clade distance on %d shared leaves: %.3f" % (nl, cd))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clade distance for two rooted trees")
    parser.add_argument("-t1", "--tree1", type=str,  required=True,
                        help="File containing newick string for tree 1")
    parser.add_argument("-t2", "--tree2", type=str, required=True,
                        help="File containing newick string for tree 2")
    main(parser.parse_args())
