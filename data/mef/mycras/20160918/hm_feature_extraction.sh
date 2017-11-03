#!/bin/bash
# Run Heteromotility in parallel on
# all XY tracks files in a Well* folder
HM=/home/jkimmel/src/heteromotility/heteromotility/heteromotility.py
cpus=$( ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w )

for well in Well*
do
  for xy in $(seq -f '%02g' 1 25)
  do
    $HM ./ ${well}/ --tracksX ${well}/tracksX_xy${xy}.csv --tracksY ${well}/tracksY_xy${xy}.csv --output_suffix xy${xy} &
  done
  wait
done
