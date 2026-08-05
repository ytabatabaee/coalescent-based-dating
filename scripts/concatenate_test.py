import multiprocessing as mp
import subprocess
import sys
import os


condition = sys.argv[1]
#r = sys.argv[1].zfill(2) # replicate
#g = '' if sys.argv[2] == '1000' else '_'+sys.argv[2] # number of genes
dataset_path = '/scratch/users/syt3/dating-data/large'
species_tree_name = 's_tree.trees'


def run_concat(i):
    r = str(i).zfill(2) 
    all_genes_path = dataset_path + '/' + condition + '/' + r + '/all-genes.phylip'
    output_path = dataset_path  + '/' + condition + '/' + r + '/concat.fasta'
    cmd = 'python2 concatenate.py ' + all_genes_path + ' > ' + output_path
    #full_cmd = '/usr/bin/time -v -o ' + output_path + '.stat' + ' -f "QR*\t%e\t%M" ' + cmd
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    _, _ = p.communicate()


if __name__ == '__main__':
    with mp.Pool(mp.cpu_count()) as p:
        p.map(run_concat, range(1,21))
