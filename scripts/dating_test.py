import multiprocessing as mp
import subprocess
import sys
import os
import math
import dendropy

r = sys.argv[1].zfill(2) # replicate
g = '' #if sys.argv[2] == '1000' else '.'+sys.argv[2] # number of genes
dataset_path = '/u/syt3/scratch/dating-data/large'
species_tree_name = 's_tree.trees'
#l = int(sys.argv[3])


def run_dating(condition):
    calib_num = int(math.sqrt(int(condition))/2)
    OG_flag = False
    OG_txt = '.no_OG' if OG_flag else ''
    #g_path_str = 'fasttree_genetrees_'+str(g_path)+'_non'
    
    data_path = dataset_path + '/' + condition + '/' + r + '/'
    species_tree_path = 'castlespro_estimatedgenetre' + g + '_' + species_tree_name + '.rooted'
    #species_tree_path = 'RAxML_result.concat_' + species_tree_name + '.rooted'
    #species_tree_path = 'castlespro_' + g_path_str + g + '_s_tree.trees.rooted'
    #species_tree_path = 'RAxML_result.concat_for_fasttree' + '_' + str(g_path) + g + '_s_tree.trees.rooted'
    #species_tree_path = 's_tree.trees.rooted'
    input_path = data_path + species_tree_path + '.labeled' + OG_txt
    calib_path = data_path + 's_tree_tu_5_calibrations_n'+str(calib_num)+ OG_txt +'.txt'
    
    '''
    # MD-CAT
    if calib_num == 1:
        output_path = data_path + 'mdcat_' + species_tree_path + '.labeled'  + OG_txt
        cmd = 'python3 ~/scratch/phylo_software/MD-Cat/md_cat.py -i ' + input_path + ' -o ' + output_path + ' -p 10 -v > ' + output_path + '.log'
    else:
        output_path = data_path + 'mdcat_' + 'n' + str(calib_num) + '_' + species_tree_path + '.labeled' + OG_txt
        cmd = 'python3 ~/scratch/phylo_software/MD-Cat/md_cat.py -i ' + input_path + ' -t ' + calib_path + ' -o ' + output_path + ' -b -v -p 10 > ' + output_path + '.log'
    '''
    '''
    # LOG-DATE
    if calib_num == 1:
        output_path = data_path + 'wlogdate_' + species_tree_path + '.labeled'  + OG_txt
        cmd = 'python3 ~/scratch/phylo_software/wLogDate/launch_wLogDate.py -i ' + input_path + ' -o ' + output_path + ' > ' + output_path + '.log'
    else:
        output_path = data_path + 'wlogdate_' + 'n' + str(calib_num) + '_' + species_tree_path + '.labeled' + OG_txt
        cmd = 'python3 ~/scratch/phylo_software/wLogDate/launch_wLogDate.py -i ' + input_path + ' -t ' + calib_path + ' -o ' + output_path  + ' -b ' + ' > ' + output_path + '.log'
    '''
    # LSD 
    #cmd = 'lsd -i ' + input_path + ' -c -a 0 -z 1 -v 1 -o ' + output_path
    #cmd = 'python3 label_and_remove_og.py -t1 ' + data_path +  '/s_tree_tu_5.trees -t2 ' + input_path + ' -o ' + input_path + '.labeled'


    # LSD2
    # unit ultrametric
    '''
    if calib_num == 1:
        output_path = data_path + 'lsd2_' + species_tree_path + '.labeled'  + OG_txt
        cmd = 'lsd2 -i ' + input_path + ' -a 0 -z 1 -s ' + str(g_path*1000) + ' -u 0.001 -o ' + output_path
    else:
        output_path = data_path + 'lsd2_' + 'n' + str(calib_num) + '_' + species_tree_path + '.labeled' + OG_txt
        lsd_calib_path = data_path + 'lsd_s_tree_tu_5_calibrations_n'+str(calib_num)+OG_txt+'.txt'
        true_sp_path = data_path + 's_tree_tu_5.trees'
        taxa_count = 100 if OG_flag else 101
        if not os.path.exists(lsd_calib_path):
            calib_text = str(taxa_count+calib_num)+'\n'
            tns = dendropy.TaxonNamespace()
            t = dendropy.Tree.get(path=true_sp_path, schema='newick', taxon_namespace=tns, rooting='force-rooted')
            tree_height = 0
            for node in t.postorder_node_iter():
                if node.taxon is not None:
                    continue
                elif node.edge.length:
                    left = node._child_nodes[0]
                    node.edge.length = node.edge.length + left.edge.length
                    tree_height = node.edge.length

            t = dendropy.Tree.get(path=input_path, schema='newick', taxon_namespace=tns, rooting='force-rooted')
            with open(calib_path, 'r') as f:
                lines = f.readlines()
            for line in lines:
                label, age = line.split()
                for node in t.postorder_node_iter():
                    if node.is_leaf():
                        node.value = node.taxon.label
                    else:
                        left = node._child_nodes[0]
                        right = node._child_nodes[1]
                        node.value = left.value + ',' + right.value
                        if node.label == label:
                            tl = node.value.split(',')
                            mrca = label
                            calib_text += 'mrca('
                            for taxa_label in tl:
                                calib_text += taxa_label + ','
                            calib_text = calib_text[:-1]
                            calib_text += ')\t' +  str(round(tree_height-float(age),3)) + '\n'    

            x = 1 if OG_flag else 0
            for i in range(x, 101):
                calib_text += str(i) + '\t' + str(round(tree_height, 3)) + '\n'

            with open(lsd_calib_path, 'w') as f:
                f.write(calib_text)
    
        cmd = 'lsd2 -i ' + input_path + ' -d ' + lsd_calib_path  + '  -s ' + str(g_path*1000) + ' -o -u 0.001 ' + output_path
    '''
    # treePL
    
    config_text = 'treefile = ' + input_path + '\n'
    if calib_num == 1:
        output_path = data_path + 'treepl_DC_' + species_tree_path + '.labeled'  + OG_txt
    else:
        output_path = data_path + 'treepl_DC_' + 'n' + str(calib_num) + '_' + species_tree_path + '.labeled' + OG_txt
    config_text += 'outfile = ' + output_path + '\n'
    config_text += 'smooth = 100\n' + 'numsites = 500000\n\n# these are the result of running prime\n'
    config_text += 'opt = 2\nmoredetail\noptad = 2\nmoredetailad\noptcvad = 2\n\n'
    if calib_num == 1:
        config_text += 'mrca = ROOT '
        taxa_count = int(condition)
        x = 1 if OG_flag else 0
        for i in range(x, taxa_count+1):
            config_text += str(i) + ' '
        config_text += '\n'
        config_text += 'min = ROOT 1\nmax = ROOT 1\n'
    else:
        tns = dendropy.TaxonNamespace()
        t = dendropy.Tree.get(path=input_path, schema='newick', taxon_namespace=tns, rooting='force-rooted')
        with open(calib_path, 'r') as f:
           lines = f.readlines()
           for line in lines:
               label, age = line.split()
               #print(label, age)
               for node in t.postorder_node_iter():
                   if node.is_leaf():
                       node.value = node.taxon.label
                   else:
                       left = node._child_nodes[0]
                       right = node._child_nodes[1]
                       node.value = left.value + ',' + right.value
                       #print(node.label, label)
                       if node.label == label:
                           tl = node.value.split(',')
                           mrca = label
                           config_text += 'mrca = ' + mrca + ' '
                           for taxa_label in tl:
                               config_text += taxa_label + ' '
                           config_text += '\n'
                           config_text += 'min = ' + mrca + ' ' +  age + '\nmax = ' + mrca + ' ' +  age + '\n'     
    
    config_text += 'thorough\nnthreads = 16\n#prime\n'
    config_path = output_path + '.config'
    with open(config_path, 'w') as f:
         f.write(config_text)
    cmd = 'treePL ' + config_path + ' > ' + output_path + '.log' 
    

    #if os.path.exists(output_path):
    #      return

    full_cmd = '/usr/bin/time -v -o ' + output_path + '.stat' + ' -f "QR*\t%e\t%M" ' + cmd
    p = subprocess.Popen(full_cmd, stdout=subprocess.PIPE, shell=True)
    _, _ = p.communicate()


def remove_labels(comb):
    g_path, calib_num = comb
    OG_flag = False
    OG_txt = '.no_OG' if OG_flag else ''
    data_path = dataset_path + '/' + r + '/'
    g_path_str = 'fasttree_genetrees_'+str(g_path)+'_non'
    #species_tree_path = 'castlespro_' + g_path_str + g + '_' + species_tree_name + '.rooted'
    species_tree_path = 'RAxML_result.concat_for_fasttree' + '_' + str(g_path) + g + '_s_tree.trees.rooted' 
    #species_tree_path = 's_tree.trees.rooted'
    #input_path = dataset_path + '/' + r + '/' + 'md_'+ species_tree_path# + '.labeled'
    if calib_num == 1:
        output_path = data_path + 'lsd2_' + species_tree_path + '.labeled'  + OG_txt
    else:
        output_path = data_path + 'lsd2_' + 'n' + str(calib_num) + '_' + species_tree_path + '.labeled' + OG_txt
    #output_path = dataset_path + '/' + r + '/' + 'lsd2_n3_' + species_tree_path# + '.labeled'
    #cmd = 'python3 ~/scratch/phylo_software/MD-Cat/md_cat.py -i ' + input_path + ' -o ' + output_path + ' > ' + output_path + '.log'
    cmd = 'python3 remove_lsd_labels.py -t ' + output_path + '.date.nexus' +  ' -o ' +  output_path + '.date.nwk'
    #cmd = 'python3 remove_mdcat_labels.py -t ' + input_path + ' -o ' + output_path 
    #cmd = 'python3 get_unit_ultrametric.py -t ' + output_path + ' -o ' + output_path + '.normalized'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    _, _ = p.communicate()


def normalize(condition):
    calib_num = int(math.sqrt(int(condition))/2)
    OG_flag = False
    OG_txt = '.no_OG' if OG_flag else ''
    #species_tree_path = 'treepl_n' + str(calib_num) + '_castlespro_estimatedgenetre' + g + '_' + species_tree_name + '.rooted'
    species_tree_path = 'treepl_n' + str(calib_num) + '_RAxML_result.concat_' + species_tree_name + '.rooted'
    #species_tree_path = 'castlespro_' + g_path_str + g + '_' + species_tree_name + '.rooted'
    #species_tree_path = 'RAxML_result.concat_for_fasttree' + '_' + str(g_path) + g + '_s_tree.trees.rooted'
    output_path = dataset_path + '/' + condition+ '/' + r + '/' + species_tree_path + '.labeled' + OG_txt
    #if 'lsd' in tree_name:
    #    output_path += '.date.nwk'
    cmd = 'python3 get_unit_ultrametric.py -t ' + output_path + ' -o ' + output_path + '.normalized'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    _, _ = p.communicate()


if __name__ == '__main__':
    #gene_tree_paths = ['truegenetrees', 'fasttree_genetrees_1600_non',
    #                  'fasttree_genetrees_200_non', 'fasttree_genetrees_400_non', 'fasttree_genetrees_800_non']
    #num_taxa = ['50', '100', '200', '500', '1000', '2000', '5000', '10000']
    #tree_names = ['mdcat_']#, 'mdcat_n3_']
    #tree_names = ['lsd2_', 'wlogdate_', 'mdcat_', 'treepl_', 'lsd2_n10_', 'lsd2_n3_','treepl_n3_', 'treepl_n10_', 'wlogdate_n10_', 'wlogdate_n3_', 'mdcat_n10_', 'mdcat_n3_']
    #tree_names = ['lsd2_n10_root_unfixed_', 'lsd2_n3_root_unfixed_','treepl_n3_root_unfixed_', 'treepl_n10_root_unfixed_', 'wlogdate_n10_root_unfixed_', 'wlogdate_n3_root_unfixed_', 'mdcat_n10_root_unfixed_', 'mdcat_n3_root_unfixed_']
    #calib_num = [1]
    #combinations = []
    #for gt in gene_tree_paths:
    #   for tn in tree_names:
    #       combinations.append((gt, tn))
    #for gt in gene_tree_paths:
    #   for c in calib_num:
    #       combinations.append((gt, c))
    run_dating(sys.argv[2])
    #with mp.Pool(mp.cpu_count()) as p:
    #    p.map(run_dating, num_taxa)
