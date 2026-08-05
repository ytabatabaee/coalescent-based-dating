import multiprocessing as mp
import subprocess
import sys
import os


condition = sys.argv[1]#.zfill(2) # replicate
#g = '' if sys.argv[2] == '1000' else '_'+sys.argv[2] # number of genes
dataset_path = '/scratch/users/syt3/dating-data/large'
species_tree_name = 's_tree.trees'


def run_indelible(i):
    r = str(i).zfill(2)
    g_path = 'truegenetrees'
    gene_tree_path = dataset_path + '/' + condition + '/' + r + '/' + g_path
    output_path = dataset_path  + '/' + condition + '/' + r + '/'
    cmd = 'python3 run_indelible.py -s 1 -e 1000 -t ' + gene_tree_path + ' -p /scratch/users/syt3/scripts/dating_large/parameters_500.csv -o ' + output_path
    #full_cmd = '/usr/bin/time -v -o ' + output_path + '.stat' + ' -f "QR*\t%e\t%M" ' + cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    _, _ = p.communicate()
    os.system("cat "+output_path+"*phy | sed '/^ *$/d' > " + output_path + "all-genes.phylip")
    os.system("rm "+output_path+"*.phy")


if __name__ == '__main__':
    #model_conditions = ['200', '500', '1000', '2000', '5000', '10000']#['50', '100']
    #for condition in model_conditions:
    #    run_aster(condition)
    with mp.Pool(mp.cpu_count()) as p:
        p.map(run_indelible, range(1, 21))
