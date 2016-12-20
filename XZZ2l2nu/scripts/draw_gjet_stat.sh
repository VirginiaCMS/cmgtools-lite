#!/bin/sh


tag="GJets_"

channels="mu el"
cutChains="tightzpt100met50"
logdir="log_gjets_36p46"

mkdir -p $logdir

for cutChain in $cutChains;
do
   for channel in $channels;
   do

#      python stack_dtmc_skim_gjets_stat.py -l -b -q  --tag="$tag" --cutChain="$cutChain" --channel="$channel" --dyGJets --LogY --muoneg --Blind &> ${logdir}/${tag}${cutChain}_log_${channel}_blind_plot.log &

      python stack_dtmc_skim_gjets_stat.py -l -b -q  --tag="$tag" --cutChain="$cutChain" --channel="$channel" --dyGJets --LogY --muoneg &> ${logdir}/${tag}${cutChain}_log_${channel}_unblind_plot.log &

      python stack_dtmc_skim_gjets_stat.py -l -b -q  --tag="$tag" --cutChain="$cutChain" --channel="$channel" --dyGJets --LogY --muoneg --doSys &> ${logdir}/${tag}${cutChain}_log_${channel}_unblind_sys_plot.log &

   done
done
