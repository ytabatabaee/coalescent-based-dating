import multiprocessing as mp
import sys
from utils import *

dataset_path = '/home/syt3/scratch/dating-data/large/'
r = sys.argv[1].zfill(2)
m = sys.argv[2]


def compute_stats():
    #r = str(i).zfill(2)
    s_tree_path = dataset_path  + '/' + m + '/' + r + '/s_tree.trees'
    dir_path = dataset_path  + '/' + m + '/' + r 
    ad = compute_ad(s_tree_path, dir_path + '/truegenetrees')
    gtee = compute_gtee(dir_path + '/truegenetrees', dir_path + '/estimatedgenetre')
    print(ad)
    print(gtee)
    if not os.path.exists(dir_path + '/ad.txt'):
        with open(dir_path + '/ad.txt', 'w') as f:
            f.write(str(ad) + '\n')
    if not os.path.exists(dir_path + '/gtee.txt'):
        with open(dir_path + '/gtee.txt', 'w') as f:
            f.write(str(gtee) + '\n')


if __name__ == '__main__':
    #for i in range(1, 21):
    #    compute_stats(i)
    #condition = ['50','100','200','500','1000','2000','5000','10000']
    compute_stats()
    #with mp.Pool(mp.cpu_count()) as p:
    #    p.map(compute_stats, range(1, 21))
