#!/bin/bash

p=~/Downloads/ASTRALIII/
for y in `seq -w 1 50`; do 
	while read line; do  
		nw_distance -ml -n $p/$y/s_tree.trees  $line; 
	done < <(grep "^$y " ./lcas.txt|cut -f 2,3 -d' '|uniq|sed -e "s/_0_0//g") | \
		 awk 'NR%2{a=$1;next}{print "'$y'", $1, a, $3, $2}';  
done|tee ./lcas-st.txt

for y in `seq -w 1 50`; do 
  for x in {0..100}; do  
    cat ASTRALIII/$y/ultrametric-genetrees.tre| nw_distance -n -ml - $(( ( RANDOM % 100 )  + 1 ))_0_0 $(( ( RANDOM % 100 )  + 1 ))_0_0  |awk 'NR%2{a=$1;next}{print "'$y'", $1, a, $3, $2}'; 
  done; 
done |tee ./lcas.tx

 cat lcas.txt |sed -e "s/_0//g" -e "s/  / /g"|awk '{print($1,($2<$3 ? $2 : $3)","($2<$3 ? $3 : $2),$4)}' > lcas-pairs.txt
cat lcas-st.txt |sed -e "s/_0//g" -e "s/  / /g"|awk '{print($1,($2<$3 ? $2 : $3)","($2<$3 ? $3 : $2),$4)}' >lcas-st-pairs.txt
