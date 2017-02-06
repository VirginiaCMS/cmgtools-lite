#!/bin/sh

tag=$1
#tag="SinglePhoton_Run2016C_23Sep2016"

njob="0"
for dir in $tag*; 
do
  echo $dir
  cd $dir
  rm -rf Loop*
  #log=`grep processed log.txt`; 
  #if [ -z "$log" ]; then 
  #  echo "run the job now" ; 
  python $CMSSW_BASE/src/PhysicsTools/HeppyCore/python/framework/looper.py pycfg.py config.pck --options=options.json &> local.log && cp -rp Loop/* . && rm -rf Loop &
  #else 
  #  echo "processed before"; 
  #fi
  cd ../
  njob=$(( njob + 1 ))
  if [ "$njob" -eq "20" ]; then
    wait
    njob="0"
  fi
done

