#!/bin/bash

s=20
g=1000
t=16.11685
pp=2500000
gen=20000000
sps=0.000001

for sp in 100; do
  ./simphy -rs "$s" -rl f:"$g" -rg 1 -sb f:0.0000001 -sd f:0 -st f:$gen -sl f:"$sp" -so f:1 -si f:1 -sp f:$pp -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 -cs 14907 -v 3 -o $sp -ot 0 -op 1 -od 1 > log_$sp.txt
  for r in `ls -d $sp/*/`; do
    cat $r/g_trees*.trees|sed -e "s/_0_0//g" >> $r/truegenetrees
    rm  $r/g_trees*.trees;
    cat $r/g_trees*.ralpha > $r/truegenetrees.ralpha;
    rm  $r/g_trees*.ralpha;
  done
done
