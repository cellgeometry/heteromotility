#!/bin/bash
rm all_tracksX.csv
rm all_tracksY.csv

for exp in 2016*
do
  for well in ${exp}/Well*
  do
    for i in $(seq -f '%02g' 1 25)
    do
      cat ${well}/tracksX_xy${i}.csv >> all_tracksX.csv
      cat ${well}/tracksY_xy${i}.csv >> all_tracksY.csv
    done
  done
done

echo 'X tracks : ' `wc -l all_tracksX.csv`
echo 'Y tracks : ' `wc -l all_tracksY.csv`
