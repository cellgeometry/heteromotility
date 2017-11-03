#!/bin/bash
# Combine individual XY point motility statistics into exp level statistics

for well in Well*
do
  rm ${well}/well_motility_statistics.csv
  for xy in $(seq -f '%02g' 1 25)
  do
    tail -n +2 ${well}/motility_statistics_xy${xy}.csv >> ${well}/well_motility_statistics.csv
  done
done

for well in WellD*
do
  cat ${well}/well_motility_statistics.csv >> tmp0.csv
done
cat header.csv > fgf2_exp_motility_statistics.csv
cat tmp0.csv >> fgf2_exp_motility_statistics.csv

for well in WellE*
do
  cat ${well}/well_motility_statistics.csv >> tmp1.csv
done
cat header.csv > nofgf2_exp_motility_statistics.csv
cat tmp1.csv >> nofgf2_exp_motility_statistics.csv

rm tmp*.csv

echo 'FGF2+ cells: ' `wc -l fgf2_exp_motility_statistics.csv`
echo 'FGF2- cells: ' `wc -l nofgf2_exp_motility_statistics.csv`
