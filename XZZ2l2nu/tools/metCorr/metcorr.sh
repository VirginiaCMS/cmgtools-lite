#!/bin/sh

#inputs
inputdir=/home/heli/XZZ/80X_20161029_light
outputdir=/home/heli/XZZ/80X_20161029_light_Skim
#inputdir=/datab/heli/XZZ/80X_20161029_light
#inputdir=/data2/XZZ2/80X_20161018_light
#inputdir=/home/heli/XZZ/80X_20161018_light
#outputdir=/home/heli/XZZ/80X_20161018_light_Skim
#inputdir=/home/heli/XZZ/80X_20161006_light
#inputdir=/home/heli/XZZ/80X_20160825_light
#outputdir=/home/heli/XZZ/80X_20161006_light_Skim
#outputdir=/home/heli/XZZ/80X_20160825_light_Skim
config=config/parameters_light

mkdir -p ${outputdir}

gmake all

njob="0"

#for infile in $inputdir/*/vvTreeProducer/tree.root ; 
#for infile in $(ls $inputdir/Single*/vvTreeProducer/tree.root ); 
#for infile in $(ls $inputdir/Single*/vvTreeProducer/tree.root ); 
#for infile in $inputdir/DYJetsToLL_M50_MGMLM/vvTreeProducer/tree.root ; 
#for infile in $inputdir/DYJetsToLL_M50_MGMLM_Ext1/vvTreeProducer/tree.root ; 
#for infile in $inputdir/DYJetsToLL_M50/vvTreeProducer/tree.root ; 
#for infile in $inputdir/*/vvTreeProducer/tree.root ;
#for infile in $inputdir/SingleEMU_*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/DYJetsToLL*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/SingleEMU_Run2016B2G_*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/SingleEMU_Run2016B2H29fbinv_*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/DYJetsToLL_M50_reHLT/vvTreeProducer/tree.root ; 
#for infile in $inputdir/Bulk*/vvTreeProducer/tree.root ; 
#for infile in $(ls $inputdir/*/vvTreeProducer/tree.root | grep -v Single | grep -v DYJets | grep -v Bulk); 
#for infile in $inputdir/SingleEMU_Run2016B2H_ReReco_33fbinv/vvTreeProducer/tree.root ; 
#for infile in $inputdir/ZZTo2L2Nu/vvTreeProducer/tree.root ; 
#for infile in $inputdir/SingleEMU_Run2016B2G_ReReco_27fbinv/vvTreeProducer/tree.root ; 
#for infile in $inputdir/DYJetsToLL_M50_BIG_JEC/vvTreeProducer/tree.root ; 
#for infile in $inputdir/SingleEMU_Run2016B2H_ReReco_36p1fbinv/vvTreeProducer/tree.root ; 
#for infile in $(ls $inputdir/*/vvTreeProducer/tree.root | grep -v Single | grep -v DYJets ); 
#for infile in $inputdir/SingleEMU_Run2016B2H_ReReco_36p22fbinv/vvTreeProducer/tree.root ; 
for infile in $inputdir/DYJetsToLL_M50_BIG/vvTreeProducer/tree.root ; 
do
  echo "+++ skimming $infile +++"
  outfile="${outputdir}/${infile/$inputdir\//}"

  # options for outputs
  #outfile="${outfile/\/vvTreeProducer\/tree/_ResBosRefit_NoRecoil}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_NoRecoil}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_ResBosRefit_Rc36p46}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_Rc36p46DtReCalib}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_DtReCalib}"
  #outfile="${outfile/\/vvTreeProducer\/tree/}"
  outfile="${outfile/\/vvTreeProducer\/tree/_test}"

  inSkimFile=${infile/vvTreeProducer\/tree.root/skimAnalyzerCount\/SkimReport.txt}

  #echo $inSkimFile
  AllEvents=`grep "All Events" ${inSkimFile} | awk {'print $3'}`
  SumWeights=`grep "Sum Weights" ${inSkimFile} | awk {'print $3'}`

  if [ -z $AllEvents ]; then
    echo Fail to get All Events from file ${inSkimFile}
    continue
  fi
  if [ -z $SumWeights ]; then
    SumWeights=$AllEvents
  fi

  echo -- Input file: $infile
  echo -- Output file: $outfile
  echo -- AllEvents: $AllEvents , SumWeights: $SumWeights
  echo -- Selection: $selection
  echo -- Command: ./bin/metcorr.exe $config $infile $outfile $AllEvents $SumWeights

  ./bin/metcorr.exe $config $infile $outfile $AllEvents $SumWeights &> ${outfile}.skim.log &

  njob=$(( njob + 1 ))
  if [ "$njob" -eq "100" ]; then
   # wait
    njob="0"
  fi

done

