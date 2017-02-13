#!/bin/sh

dds="mc_qcdem1 mc_qcdmu mc_wjetht"

for nn in `seq 1 1 100`;
do 
  for dd in $dds;
  do
    ./check.sh $dd &> ck_${dd}.log &
  done
  #sleep 1800
  wait
done

