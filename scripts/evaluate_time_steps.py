import argparse
import dendropy
import numpy as np
import pandas as pd
from utils import *
import math

if __name__ == "__main__":
    num_taxa = [50, 100, 200, 500, 1000, 2000, 5000, 10000]
    ru_txts = ['']#, 'root_unfixed_']
    ru_states = [True]#, False]
    genes = [1000]
    methods = ['TreePL+CASTLES-Pro', 'TreePL+ConBL']
    df_branches = pd.DataFrame(columns=["Condition", "Method", "replicate", "AD", "GTEE", "Branch Type", "l-true", "l-est"])
    gt_type = 'estimatedgenetre'
    dataset_path = '/scratch/users/syt3/dating-data/large/'

    OG_flag = False
    OG_txt = '.no_OG' if OG_flag else ''
    #header = ['treepl_'+calib+'_castlespro_', 'treepl_'+calib+'_RAxML_result.concat_']
    footer = ['_s_tree.trees.rooted.labeled'] * 2
    df_time = pd.DataFrame(columns=["Condition", "genes", "Method", "replicate", "Step", "time_s", "mem_gb"])
    dataset_path = '/scratch/users/syt3/dating-data/large/'
    for j in range(len(num_taxa)):
        condition = str(num_taxa[j])
        calib = str(int(math.sqrt(int(condition))/2))
        header = ['treepl_n'+calib+'_castlespro_estimatedgenetre', 'treepl_n'+calib+'_RAxML_result.concat']
        for i in range(len(methods)):
            method = methods[i]
            if condition=='truegenetrees':
                continue
            for g in genes:
                g_str = '' if g == 1000 else '.'+str(g)
                for r in range(1, 21):
                    replicate = str(r).zfill(2)
                    true_tree_path = dataset_path+condition+'/'+replicate+'/s_tree_unit_ultrametric.trees'
                    try:
                        if 'ConBL' in method:
                            tree_path = dataset_path +condition+'/'+ replicate + '/' + 'RAxML_result.concat_s_tree.trees'
                            est_tree_path = dataset_path +condition+'/'+ replicate + '/' + header[i] +g_str+ footer[i]
                            tree_time, tree_mem = get_time_memory(tree_path + '.stat')
                            dating_time, dating_mem = get_time_memory(est_tree_path + '.stat')
                        else:
                            tree_path = dataset_path+condition+'/'+replicate+'/'+'castlespro_estimatedgenetre_s_tree.trees'
                            est_tree_path = dataset_path+condition+'/'+replicate+'/'+header[i]+g_str+footer[i]
                            tree_time, tree_mem = get_time_memory(tree_path + '.stat')
                            dating_time, dating_mem = get_time_memory(est_tree_path + '.stat')
                        #df_time.loc[len(df_time.index)] = [condition, g, method, r, tree_time + dating_time, max(tree_mem, dating_mem)]
                        df_time.loc[len(df_time.index)] = [condition, g, method, r, 'Branch length estimation', tree_time, tree_mem]
                        df_time.loc[len(df_time.index)] = [condition, g, method, r, 'Dating', dating_time, dating_mem]
                    except Exception as e:
                        print(e)
    df_time.to_csv('large_time_treepl_step.csv')
