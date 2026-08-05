import multiprocessing as mp
import subprocess
import sys
import os
import math

r = sys.argv[1].zfill(2) # replicate
#g = '' if sys.argv[2] == '1000' else '_'+sys.argv[2] # number of genes
dataset_path = '/scratch/users/syt3/dating-data/large'
species_tree_name = 's_tree.trees'
#l = int(sys.argv[3])


def simulate(condition):
    #r = str(idx).zfill(2)
    OG_flag = False
    OG_txt = '.no_OG' if OG_flag else ''
    #data_path = dataset_path + '/' + r + '/' 
    #input_path = data_path + 's_tree_tu_5.trees' +OG_txt
    #output_path = data_path + 's_tree_unit_ultrametric.trees' +OG_txt
    output_path = dataset_path + '/' + condition + '/' + r + '/' 
    cmd = 'python3 get_unit_ultrametric.py -t ' + output_path + 's_tree_tu_5.trees' + ' -o ' + output_path + 's_tree_unit_ultrametric.trees'
    #calib = 10
    #species_tree_path = 'castlespro_' + g_path + g + '_' + species_tree_name + '.rooted'
    #species_tree_path = 'RAxML_result.concat_for_fasttree' + '_' + str(l) + g + '_s_tree.trees.rooted' 
    #cmd = 'python3 generation_to_time.py -gt 5 -t ' + output_path + '/s_tree_gu.trees' + ' -o ' + output_path +  's_tree_tu_5.trees'
    #cmd = 'python3 label_and_remove_og.py -t1 ' + output_path +  '/s_tree_tu_5.trees -t2 ' + output_path + species_tree_path + ' -o ' + output_path + 'RAxML_result.concat_for_fasttree' + '_' + str(l) + g + '_s_tree_labeled.trees.rooted'
    #cmd = 'python3 label_and_remove_og.py -t1 ' + output_path +  '/s_tree_tu_5.trees -t2 ' + output_path + species_tree_path + ' -o ' + output_path + species_tree_path + '.labeled'
    #cmd = 'python3 simulate_calibs_unfixed_root.py -n ' + str(calib) + ' -t ' + output_path + '/s_tree_tu_5.trees'+ OG_txt + ' -o ' + output_path +  's_tree_tu_5_calibrations_n'+ str(calib) + '_root_unfixed' + OG_txt
    #print(cmd)
    #calib = int(math.sqrt(int(condition))/2)
    #cmd = 'python3 simulate_calibrations.py -n ' + str(calib) + ' -t ' + output_path + '/s_tree_tu_5.trees' + ' -o ' + output_path +  's_tree_tu_5_calibrations_n' + str(calib)
    #full_cmd = '/usr/bin/time -v -o ' + output_path + '.stat' + ' -f "QR*\t%e\t%M" ' + cmd
    #if os.path.exists(output_path +  's_tree_tu_5_calibrations_n'+ str(calib) +OG_txt+'.txt'):
    #    return
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    _, _ = p.communicate()


if __name__ == '__main__':
    #gene_tree_paths = ['fasttree_genetrees_1600_non',
    #                   'fasttree_genetrees_200_non', 'fasttree_genetrees_400_non', 'fasttree_genetrees_800_non']
    #gene_tree_paths = ['']
    conditions = ['50', '100', '200', '500', '1000', '2000', '5000', '10000']
    with mp.Pool(mp.cpu_count()) as p:
        p.map(simulate, conditions)
