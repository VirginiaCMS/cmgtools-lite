#!/bin/sh

dds="dt_elb dt_eld dt_ele dt_elf dt_elg dt_elh dt_mub dt_mud dt_muh  "
#dds="dt_elb dt_elc dt_eld dt_ele dt_elf dt_elg dt_elh "
#dds="dt_mub dt_muc dt_mud dt_mue dt_muf dt_mug dt_muh "
#dds="dt_phb dt_phc dt_phd dt_phe dt_phf dt_phg dt_phh"
#dds="$dds dt_elc dt_mub dt_muc dt_mud dt_mue dt_muf "

for nn in `seq 1 1 1000`;
do 
  for dd in $dds;
  do
    ./check.sh $dd &> ck_${dd}.log &
  done
  #sleep 1800
  wait
done

