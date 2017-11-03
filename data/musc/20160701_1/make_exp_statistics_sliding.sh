#!/bin/bash
# Combine individual XY point motility statistics into exp level statistics

rm tmp*.csv

for s in $(seq 1 14 15)
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
    for well in WellF*
    do
      cat ${well}/well_motility_statistics_split_${tau}_sliding_s${s}.csv >> tmp_fgf${tau}_sliding_s${s}.csv
    done
    head -n1 header.csv > fgf2_exp_motility_statistics_split_${tau}_sliding_s${s}.csv
    cat tmp_fgf${tau}_sliding_s${s}.csv >> fgf2_exp_motility_statistics_split_${tau}_sliding_s${s}.csv
  done

  for tau in $(seq -f '%02g' 20 5 30)
  do
    for well in WellG*
    do
      cat ${well}/well_motility_statistics_split_${tau}_sliding_s${s}.csv >> tmp_nofgf${tau}_sliding_s${s}.csv
    done
    head -n1 header.csv > nofgf2_exp_motility_statistics_split_${tau}_sliding_s${s}.csv
    cat tmp_nofgf${tau}_sliding_s${s}.csv >> nofgf2_exp_motility_statistics_split_${tau}_sliding_s${s}.csv
  done
done
rm tmp*.csv

echo '------------------------------------------------------------'
echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics.csv`
echo '------------------------------------------------------------'
echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics_split_*.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics_split_*.csv`
echo '------------------------------------------------------------'
