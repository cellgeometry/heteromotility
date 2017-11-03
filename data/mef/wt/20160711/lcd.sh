rm Well*/*lcdensity*.csv

for well in Well*
do
  python3 ~/src/hm_analysis/python/local_cell_density.py ${well}
done

for well in Well*
do
  for xy in $(seq -f '%02g' 1 25)
  do
    cat ${well}/lcdensity_xy${xy}.csv >> ${well}/well_lcdensity.csv
  done
done

for well in Well*
do
  cat ${well}/well_lcdensity.csv >> exp_lcdensity.csv
done

wc -l exp_lcdensity.csv
