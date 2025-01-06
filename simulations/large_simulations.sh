#!/bin/bash

s=20
g=1000
min_pp=20000
max_pp=2000000
gen=20000000
sps=0.000001

for sp in 50 100 200 500 1000 2000 5000 10000; do
  ./simphy_mac64 -rs "$s" -rl f:"$g" -rg 1 -sb lu:0.0000001,0.000001 -sd lu:0.0000001,sb -st ln:16,1 -sl f:"$sp" -so f:1 -si f:1 -sp u:$min_pp,$max_pp -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 -cs 14907 -v 3 -o $sp -ot 0 -op 1 -od 1 > log_$sp.txt
  for r in `ls -d $sp/*/`; do
    cat $r/g_trees*.trees|sed -e "s/_0_0//g" >> $r/truegenetrees
    rm  $r/g_trees*.trees;
    #cat $r/g_trees*.ralpha > $r/truegenetrees.ralpha;
    #rm  $r/g_trees*.ralpha;
    #python3 ../scripts/compute_discord.py $r/s_tree.trees $r/truegenetrees
  done
done
