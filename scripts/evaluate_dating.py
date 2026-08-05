import argparse
import matplotlib.pyplot as plt
import dendropy
import math
import pandas as pd
import seaborn as sns
import numpy as np


def compare_bl(tree_path1, tree_path2, df_branches, condition, method, replicate, ad, gtee):
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(path=tree_path1, schema='newick', taxon_namespace=tns)
    t2 = dendropy.Tree.get(path=tree_path2, schema='newick', taxon_namespace=tns)

    t1.deroot()
    t2.deroot()

    length_diffs = dendropy.calculate.treecompare._get_length_diffs(t1, t2)

    idx = 0
    for node in t1.postorder_node_iter():
        if not node.parent_node:
            continue
        (l1, l2) = length_diffs[idx]
        node_type = 'terminal' if node.is_leaf() else 'internal'
        node_label = node.taxon.label if node.is_leaf() else ''
        df_branches.loc[len(df_branches.index)] = [condition, method, replicate, ad, gtee, node_type, l1, l2]
        idx += 1

    return df_branches


if __name__ == "__main__":
    sns.set_theme()
    num_taxa = [50, 100, 200, 500, 1000, 2000, 5000, 10000]
    ru_txts = ['']#, 'root_unfixed_']
    ru_states = [True]#, False]
    methods = ['TreePL+CASTLES-Pro', 'TreePL+Concat(RAxML)']
    df_branches = pd.DataFrame(columns=["Condition", "Method", "replicate", "AD", "GTEE", "Branch Type", "l-true", "l-est"])
    gt_type = 'estimatedgenetre'
    dataset_path = '/scratch/users/syt3/dating-data/large/'

    #if calib == '':
    #footer = ['_s_tree.trees.rooted.labeled']*2
    #else:
    footer = ['_s_tree.trees.rooted.labeled.normalized']*2
    #for u in range(len(ru_txts)):
    #if calib == '' and not ru_states[u]:
    #continue
    #header = ['treepl_n'+calib+'castlespro_', 'treepl_'+calib+'RAxML_result.concat_align']
    #if calib == '':
    #    header += ['s_tree_tu_5_'+calib+ru_txts[u]+'mcmctree.date.nwk']
    #else:
    #    header += ['s_tree_tu_5_calib_'+calib+ru_txts[u]+'mcmctree.date.nwk']
    for j in range(len(num_taxa)):
        print(num_taxa[j])
        condition = str(num_taxa[j])
        calib = str(int(math.sqrt(int(condition))/2))
        header = ['treepl_n'+calib+'_castlespro_', 'treepl_n'+calib+'_RAxML_result.concat']
        for i in range(len(methods)):
            method = methods[i]
            for r in range(1, 21):
                #print(condition, method, r)
                replicate = str(r).zfill(2)
                true_tree_path = dataset_path + condition + '/' +replicate+'/s_tree_unit_ultrametric.trees'
                #true_tree_path = dataset_path + condition + '/' +replicate+'/s_tree_tu_5.trees'
                ad_path = dataset_path + condition + '/' +replicate + '/ad.txt'
                gtee_path = dataset_path + condition + '/' +replicate+ '/gtee.txt'
                try:
                    with open(ad_path) as f:
                        ad = float(f.read())
                    with open(gtee_path) as g:
                        gtee = float(g.read())
                    if 'Concat(RAxML)' in method:
                        est_tree_path = dataset_path + condition + '/' + replicate + '/' + header[i] + footer[i]
                    else:
                        est_tree_path = dataset_path + condition + '/' + replicate + '/' + header[i] + gt_type + footer[i]
                    df_branches = compare_bl(true_tree_path, est_tree_path, df_branches, condition, method, replicate, ad, gtee)
                except Exception as e:
                    print(e)

    df_branches.to_csv('large_estgt_dating_normalized.csv')
