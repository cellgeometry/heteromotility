#!/bin/bash
# Combine individual XY point motility statistics into exp level statistics

MAXTAU=20

rm tmp*.csv

for well in Well*
do
  rm ${well}/well_motility_statistics_split*.csv
  for tau in $(seq -f '%02g' 20 $MAXTAU)
  do
    for xy in $(seq -f '%02g' 1 25)
    do
      tail -n +2 ${well}/motility_statistics_split_${tau}_xy${xy}.csv >> ${well}/well_motility_statistics_split_${tau}.csv
    done
  done
done

for tau in $(seq -f '%02g' 20 $MAXTAU)
do
  for well in Well*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_${tau}.csv
  done
  head -n1 header.csv > exp_motility_statistics_split_${tau}.csv
  cat tmp_${tau}.csv >> exp_motility_statistics_split_${tau}.csv
done

rm tmp*.csv

MAXTAU=25

rm tmp*.csv

for well in Well*
do
  rm ${well}/well_motility_statistics_split*.csv
  for tau in $(seq -f '%02g' 25 $MAXTAU)
  do
    for xy in $(seq -f '%02g' 1 25)
    do
      tail -n +2 ${well}/motility_statistics_split_${tau}_xy${xy}.csv >> ${well}/well_motility_statistics_split_${tau}.csv
    done
  done
done

for tau in $(seq -f '%02g' 25 $MAXTAU)
do
  for well in Well*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_${tau}.csv
  done
  head -n1 header.csv > exp_motility_statistics_split_${tau}.csv
  cat tmp_${tau}.csv >> exp_motility_statistics_split_${tau}.csv
done

rm tmp*.csv

MAXTAU=30

rm tmp*.csv

for well in Well*
do
  rm ${well}/well_motility_statistics_split*.csv
  for tau in $(seq -f '%02g' 30 $MAXTAU)
  do
    for xy in $(seq -f '%02g' 1 25)
    do
      tail -n +2 ${well}/motility_statistics_split_${tau}_xy${xy}.csv >> ${well}/well_motility_statistics_split_${tau}.csv
    done
  done
done

for tau in $(seq -f '%02g' 30 $MAXTAU)
do
  for well in Well*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_${tau}.csv
  done
  head -n1 header.csv > exp_motility_statistics_split_${tau}.csv
  cat tmp_${tau}.csv >> exp_motility_statistics_split_${tau}.csv
done

rm tmp*.csv

echo '------------------------------------------------------------'
echo 'Total cells: ' `wc -l exp_motility_statistics.csv`
echo '------------------------------------------------------------'
echo 'Total cells: ' `wc -l exp_motility_statistics_split_*.csv`
echo '------------------------------------------------------------'
