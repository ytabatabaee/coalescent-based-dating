import multiprocessing as mp
import sys
from utils import *


if __name__ == '__main__':
    s_tree_path = sys.argv[1]    
    g_tree_path = sys.argv[2]    
    ad = compute_ad(s_tree_path, g_tree_path)
    print(ad)