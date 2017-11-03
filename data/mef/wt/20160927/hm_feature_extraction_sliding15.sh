#!/bin/bash
# Run Heteromotility in parallel on
# all XY tracks files in a Well* folder
HM=/home/jkimmel/src/heteromotility/heteromotility/heteromotility.py
cpus=$( ls -d /sys/devices/system/cpu/cpu[[:digit:]]* | wc -w )

rm Well*/motility_statistics_split_*sliding15.csv
# for well in Well*
# do
#   for xy in $(seq -f '%02g' 1 25)
#   do
#     $HM ./ ${well}/ --tracksX ${well}/tracksX_xy${xy}.csv --tracksY ${well}/tracksY_xy${xy}.csv --output_suffix xy${xy}_sliding_s15 --detailedbalance 20 --sliding_window 15 --max_windows 3 &
#   done
#   wait
# done

for well in Well*
do
  for xy in $(seq -f '%02g' 1 25)
  do
    $HM ./ ${well}/ --tracksX ${well}/tracksX_xy${xy}.csv --tracksY ${well}/tracksY_xy${xy}.csv --output_suffix xy${xy}_sliding_s15 --detailedbalance 25 --sliding_window 15 --max_windows 3 &
  done
  wait
done

# for well in Well*
# do
#   for xy in $(seq -f '%02g' 1 25)
#   do
#     $HM ./ ${well}/ --tracksX ${well}/tracksX_xy${xy}.csv --tracksY ${well}/tracksY_xy${xy}.csv --output_suffix xy${xy}_sliding_s15 --detailedbalance 30 --sliding_window 15 --max_windows 3 &
#   done
#   wait
# done
