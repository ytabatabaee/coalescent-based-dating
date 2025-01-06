import dendropy
import numpy as np
import pandas as pd
from compare_trees import *
from clade_distance import *


def get_time_memory(f_name):
    '''
    Compute time and memory from stat file
    :param f_name: file name
    :return: float, float: wall time in seconds, peak memory usage in GB
    '''
    with open(f_name) as f:
        stat = f.read()
        t_start = stat.find('Elapsed (wall clock) time (h:mm:ss or m:ss):')
        t_end = stat.find('Average shared text size (kbytes):')
        time_stat = stat[t_start:t_end - 1].split()[-1]
        time_hms = time_stat.split(':')
        h, m, s = 0, 0, 0
        s = float(time_hms[-1])
        if len(time_hms) > 1:
            m = float(time_hms[-2])
        if len(time_hms) > 2:
            h = float(time_hms[-3])
        time_sec = s + m * 60 + h * 3600
        m_start = stat.find('Maximum resident set size (kbytes):')
        m_end = stat.find('Average resident set size (kbytes):')
        mem_stat = stat[m_start:m_end - 1].split()[-1]
        peak_mem_gb = int(mem_stat) / (1000 * 1000)
        return time_sec, peak_mem_gb


def compute_ad(s_tree_path, g_trees_path):
    '''
    Compute average RF distance (AD level) between model species tree and true gene trees
    :param s_tree_path: path to file containing true species tree
    :param g_trees_path: path to file containing true gene trees
    :return: float: AD level
    '''
    tns = dendropy.TaxonNamespace()
    g_trees = dendropy.TreeList.get(path=g_trees_path, schema='newick', rooting="force-unrooted", taxon_namespace=tns)
    ad = 0
    for g in g_trees:
        #try:
        s_tree = dendropy.Tree.get(path=s_tree_path, schema='newick', rooting="force-unrooted", taxon_namespace=tns)
        _, _, _, _, _, rf = compare_trees(g, s_tree)
        ad += rf
        #except:
        #    print('remove gene')
    return ad / len(g_trees)


def compute_gtee(true_g_trees_path, est_g_trees_path):
    '''
    Compute mean GTEE between true gene trees and estimated gene trees
    :param true_g_trees_path: path to file containing true gene trees
    :param est_g_trees_path: path to file containing estimated gene trees
    :return: float: mean gtee rate
    '''
    tns = dendropy.TaxonNamespace()
    true_g_trees = dendropy.TreeList.get(path=true_g_trees_path, schema='newick', rooting="force-unrooted", taxon_namespace=tns)
    est_g_trees = dendropy.TreeList.get(path=est_g_trees_path, schema='newick', rooting="force-unrooted", taxon_namespace=tns)
    mgte = 0
    k = len(true_g_trees)
    for idx in range(k):
        _, _, _, _, _, rf = compare_trees(true_g_trees[idx], est_g_trees[idx])
        mgte += rf
    return mgte / k
