#!/bin/bash

touch mem_res.txt
touch time_res.txt

declare -a array=("10000/" "20000/" "50000/" "100000/" "200000/" "400000/" "600000/" "800000/" "1000000/" "./")
arraylength=${#array[@]}

for (( i=0; i<${arraylength}; i++ ));
do
  tail -n6 "${array[$i]}"method.Rout | head -n1 >> mem_res.txt
  tail -n1 "${array[$i]}"method.Rout >> time_res.txt
done
