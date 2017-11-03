#!/bin/bash
# Combine individual XY point motility statistics into exp level statistics

rm tmp*.csv

for s in $(seq 15 15)
do
  for well in Well*
  do
    rm ${well}/well_motility_statistics_split*sliding*.csv
    for tau in $(seq -f '%02g' 20 5 30)
    do
      for xy in $(seq -f '%02g' 1 25)
      do
        tail -n +2 ${well}/motility_statistics_split_${tau}_xy${xy}_sliding_s${s}.csv >> ${well}/well_motility_statistics_split_${tau}_sliding_s${s}.csv
      done
    done
  done

  for tau in $(seq -f '%02g' 20 5 30)
  do
    for well in Well*
    do
      cat ${well}/well_motility_statistics_split_${tau}_sliding_s${s}.csv >> tmp_${tau}_sliding_s${s}.csv
    done
    head -n1 header.csv > exp_motility_statistics_split_${tau}_sliding_s${s}.csv
    cat tmp_${tau}_sliding_s${s}.csv >> exp_motility_statistics_split_${tau}_sliding_s${s}.csv
  done

done
rm tmp*.csv

echo '------------------------------------------------------------'
echo 'Total cells: ' `wc -l exp_motility_statistics_sliding_s15.csv`
echo '------------------------------------------------------------'
