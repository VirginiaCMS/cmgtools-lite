#!/bin/sh

outdir=mc_phdynlo
config=run_xzz2l2nu_80x_cfg_photon_mc.py
#config=run_xzz2l2nu_80x_cfg_loose_mc.py
otherfiles=" pogRecipes.py "

heppy_batch.py -o ${outdir} ${config}  -b 'bsub -q sssss < ./batchScript.sh' -n
for dd in ${outdir}/* ; 
do cp $otherfiles ${dd}/;
done

cd $outdir
heppy_check.py * -b "bsub -q 1nd"

cd ../

#heppy_batch.py -o gjets run_xzz2l2nu_80x_cfg_photon_mc.py -b 'bsub -q ssssss < ./batchScript.sh'

#cd mc1 
#heppy_check.py * -b "bsub -q 2nd"
#cd -
