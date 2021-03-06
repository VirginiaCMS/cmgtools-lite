#!/bin/sh

# ask Hengne.Li@cern.ch for questions.

# Warning!! 
# always test before really run this script, 
# it has potential to delete your good and finished job folders!!!!

# privide chunks running directory, 
# and directory to copy jobs out

dir=$1
#dir=dt_b2h_mu_d
#dir=dt_b2h_mu_c
#dir=dt_b2h_mu
#dir=gjetsdt_resub
out=/data2/XZZ2/80X_20170202_Chunks
#out=/datac/heli/XZZ2/80X_20170202_Chunks
#out=/datad/heli/XZZ2/80X_20170202_Chunks

mkdir -p $out

# need to check your jobs with the expected n root files and n pck files
# to verify if the job is finished sucessfully.
n_root_files="3"
n_pck_files="14"
#n_pck_files="13"


if [ ! -e "$out" ]; then
  echo "ERROR:: Do not exist output directory $out, exist... "
  exit 0
fi

#echo "#### first do a overall sync"
#rsync -var $dir $out/

echo "#### go to check each directory "
cd $dir

list=`ls -1`

for job in $list;
do
  echo "## check $job "
  if [ -e ${job}/vvTreeProducer/tree.root ]; then 
    n1=`rootls ${job}/*/*.root | grep root | wc -l`;
    n2=`ls -l ${job}/* | grep pck |wc -l`;

    if [ "$n1" -eq "$n_root_files" -a  "$n2" -eq "$n_pck_files" ]; then
      echo "- job is done correctly with ${n1} root files and ${n2} pck files."
      echo "  - copy out and delete .. "
      echo " > rsync -var $job $out/$dir/"
      rsync -var $job $out/$dir/ && rm -rf $job &
      #echo " > rm -rf $job"
      #rm -rf $job
    fi
  else
    echo "- job is not finished or has problem to be resubmitted .. with ${n1} root files and ${n2} pck files. "
  fi
done

cd ../


