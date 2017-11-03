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
  for well in WellB*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_fgf${tau}.csv
  done
  for well in WellD*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_fgf${tau}.csv
  done
  head -n1 header.csv > fgf2_exp_motility_statistics_split_${tau}.csv
  cat tmp_fgf${tau}.csv >> fgf2_exp_motility_statistics_split_${tau}.csv
done

for tau in $(seq -f '%02g' 20 $MAXTAU)
do
  for well in WellC*
  do
    cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_nofgf${tau}.csv
  done
  for well in WellE*
  do
    cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_nofgf${tau}.csv
  done
  head -n1 header.csv > nofgf2_exp_motility_statistics_split_${tau}.csv
  cat tmp_nofgf${tau}.csv >> nofgf2_exp_motility_statistics_split_${tau}.csv
done

rm tmp*.csv

echo '------------------------------------------------------------'
echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics.csv`
echo '------------------------------------------------------------'
echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics_split_*.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics_split_*.csv`
echo '------------------------------------------------------------'

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
  for well in WellB*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_fgf${tau}.csv
  done
  for well in WellD*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_fgf${tau}.csv
  done
  head -n1 header.csv > fgf2_exp_motility_statistics_split_${tau}.csv
  cat tmp_fgf${tau}.csv >> fgf2_exp_motility_statistics_split_${tau}.csv
done

for tau in $(seq -f '%02g' 25 $MAXTAU)
do
  for well in WellC*
  do
    cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_nofgf${tau}.csv
  done
  for well in WellE*
  do
    cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_nofgf${tau}.csv
  done
  head -n1 header.csv > nofgf2_exp_motility_statistics_split_${tau}.csv
  cat tmp_nofgf${tau}.csv >> nofgf2_exp_motility_statistics_split_${tau}.csv
done

rm tmp*.csv

echo '------------------------------------------------------------'
echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics.csv`
echo '------------------------------------------------------------'
echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics_split_*.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics_split_*.csv`
echo '------------------------------------------------------------'

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
  for well in WellB*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_fgf${tau}.csv
  done
  for well in WellD*
  do
      cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_fgf${tau}.csv
  done
  head -n1 header.csv > fgf2_exp_motility_statistics_split_${tau}.csv
  cat tmp_fgf${tau}.csv >> fgf2_exp_motility_statistics_split_${tau}.csv
done

for tau in $(seq -f '%02g' 30 $MAXTAU)
do
  for well in WellC*
  do
    cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_nofgf${tau}.csv
  done
  for well in WellE*
  do
    cat ${well}/well_motility_statistics_split_${tau}.csv >> tmp_nofgf${tau}.csv
  done
  head -n1 header.csv > nofgf2_exp_motility_statistics_split_${tau}.csv
  cat tmp_nofgf${tau}.csv >> nofgf2_exp_motility_statistics_split_${tau}.csv
done

rm tmp*.csv

echo '------------------------------------------------------------'
echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics.csv`
echo '------------------------------------------------------------'
echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics_split_*.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics_split_*.csv`
echo '------------------------------------------------------------'
