#!/bin/bash
# Combine individual XY point motility statistics into exp level statistics

cat header.csv > exp_motility_statistics.csv
for well in Well*
do
  for xy in $(seq -f '%02g' 1 25)
  do
    for f in ${well}/motility_statistics_xy${xy}.csv
    do
      tail -n +2 ${f} >> exp_motility_statistics.csv
    done
  done
done

echo 'Total cells: ' `wc -l exp_motility_statistics.csv`
