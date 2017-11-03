#!/bin/bash
# Run Heteromotility in parallel on
# all XY tracks files in a Well* folder
HM=/home/jkimmel/src/heteromotility/heteromotility/heteromotility.py
cpus=$( ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w )

rm Well*/motility_statistics_split_*.csv
for well in Well*
do
  for xy in $(seq -f '%02g' 1 25)
  do
    $HM ./ ${well}/ --tracksX ${well}/tracksX_xy${xy}.csv --tracksY ${well}/tracksY_xy${xy}.csv --output_suffix xy${xy} --detailedbalance 20 --dbmax 20 &
  done
  wait
done

for well in Well*
do
  for xy in $(seq -f '%02g' 1 25)
  do
    $HM ./ ${well}/ --tracksX ${well}/tracksX_xy${xy}.csv --tracksY ${well}/tracksY_xy${xy}.csv --output_suffix xy${xy} --detailedbalance 25 --dbmax 25 &
  done
  wait
done

for well in Well*
do
  for xy in $(seq -f '%02g' 1 25)
  do
    $HM ./ ${well}/ --tracksX ${well}/tracksX_xy${xy}.csv --tracksY ${well}/tracksY_xy${xy}.csv --output_suffix xy${xy} --detailedbalance 30 --dbmax 30 &
  done
  wait
done
