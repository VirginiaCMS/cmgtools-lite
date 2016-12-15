#!/bin/sh

#inputs
inputdir=/home/heli/XZZ/80X_20161029_GJets_light
outputdir=/home/heli/XZZ/80X_20161029_GJets_light_Skim
#outputdir=/home/heli/XZZ/80X_20161029_light_Skim
#inputdir=/home/heli/XZZ/80X_20161029_GJets_light
#outputdir=/home/heli/XZZ/80X_20161029_light_Skim
#inputdir=/home/heli/XZZ/80X_20161018_light
#outputdir=/home/heli/XZZ/80X_20161018_light_Skim
#inputdir=/home/heli/XZZ/80X_20161006_light
#outputdir=/home/heli/XZZ/80X_20161006_light_Skim
#inputdir=/home/heli/XZZ/80X_20160927_light
#outputdir=/home/heli/XZZ/80X_20160927_light_Skim
config=config/parameters_light_gjets

mkdir -p ${outputdir}

gmake all

njob="0"

#for infile in $inputdir/SinglePhoton_Run2016BCD_PromptReco/vvTreeProducer/tree.root ; 
#for infile in $inputdir/SinglePhoton_Run2016B2G_PromptReco/vvTreeProducer/tree.root ; 
#for infile in $inputdir/SinglePhoton_Run2016B2H29fbinv_PromptReco/vvTreeProducer/tree.root ; 
#for infile in $inputdir/SinglePhoton_Run2016B2H_ReReco_33fbinv/vvTreeProducer/tree.root ; 
#for infile in $inputdir/GJet_Pt_20toInf_DoubleEMEnriched/vvTreeProducer/tree.root ; 
#for infile in $(ls $inputdir/*/vvTreeProducer/tree.root | grep -v SinglePhoton | grep -v GJet_Pt_ ); 
#for infile in $inputdir/GJet_Pt_*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/ZNuNuGJets*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/GJets_HT40to*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/WJetsToLNu_HT*_BIG/vvTreeProducer/tree.root ; 
#for infile in $inputdir/GJets_HT*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/QCD_*/vvTreeProducer/tree.root ; 
#for infile in $inputdir/SinglePhoton_Run2016B2H_ReReco_36p22fbinv/vvTreeProducer/tree.root ; 
#for infile in $(ls $inputdir/*/vvTreeProducer/tree.root | grep -v SinglePhoton ); 
#for infile in $inputdir/SinglePhoton_Run2016B2H_ReReco_36p46/vvTreeProducer/tree.root ; 
#for infile in $inputdir/WGToLNuG/vvTreeProducer/tree.root ; 
for infile in $(ls $inputdir/*/vvTreeProducer/tree.root | grep -v "/GJet" | grep -v "/QCD" | grep -v "/SinglePhoton" ); 
do
  echo "+++ skimming $infile +++"
  outfile="${outputdir}/${infile/$inputdir\//}"

  # options for outputs
  #outfile="${outfile/\/vvTreeProducer\/tree/_ResBos_Rc36p46ReCalibwHLT}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_ResBos_Rc36p46wHLT}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_Rc36p46ReCalibErr}"
  outfile="${outfile/\/vvTreeProducer\/tree/}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_NoRecoil}"
  #outfile="${outfile/\/vvTreeProducer\/tree/_test}"

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

