#!/usr/bin/env python

import ROOT
import os,sys, string, math, pickle
from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter
ROOT.gROOT.SetBatch()
if len(sys.argv)!=2:
    sys.exit()
tag="test_{0}_".format(sys.argv[1])
nonresscale='0.672679894023'
channel='mu'
LogY=("_log_" in tag)
Blind=("_blind_" in tag)
test=True
DrawLeptons=True
doRhoScale=True

if test: DrawLeptons = False

lepsf="trgsf*isosf*idsf*trksf*1.02942"

if doRhoScale: 
    tag+="RhoWt_"
    lepsf=lepsf+"*(0.32+0.42*TMath::Erf((rho-4.16)/4.58)+0.31*TMath::Erf((rho+115.00)/29.58))" # b2h rereco 36.1 fb-1

outdir='plots'

indir='/home/heli/XZZ/80X_20170124_light_Skim/'

lumi=36.814
sepSig=True
doRatio=True
FakeData=False
UseMETFilter=True

puWeight='puWeightsummer16'
ZJetsZPtWeight=True

k=1 # signal scale

elChannel='((abs(llnunu_l1_l1_pdgId)==11||abs(llnunu_l1_l2_pdgId)==11)&&llnunu_l1_l1_pt>115&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)'
muChannel='((abs(llnunu_l1_l1_pdgId)==13||abs(llnunu_l1_l2_pdgId)==13)&&llnunu_l1_l1_pt>50&&abs(llnunu_l1_l1_eta)<2.4&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4&&(llnunu_l1_l1_highPtID>0.99||llnunu_l1_l2_highPtID>0.99))'

if not os.path.exists(outdir): os.system('mkdir '+outdir)

tag = tag+puWeight+'_'


if UseMETFilter: tag = tag+'metfilter_'

if not Blind: tag = tag+'unblind_'

tag = tag+channel+'_'

paveText="#sqrt{s} = 13 TeV 2016 L = "+"{:.3}".format(float(lumi))+" fb^{-1}"

metfilter='(Flag_EcalDeadCellTriggerPrimitiveFilter&&Flag_HBHENoiseIsoFilter&&Flag_goodVertices&&Flag_HBHENoiseFilter&&Flag_globalTightHalo2016Filter&&Flag_eeBadScFilter)'

cuts_zmass="(llnunu_l1_mass>50&&llnunu_l1_mass<180)&&!(llnunu_l1_mass>70&&llnunu_l1_mass<110)"
cuts_zmassin="(llnunu_l1_mass>50&&llnunu_l1_mass<180)"

cuts=''
if "metzpt30" in tag:cuts=cuts_zmass+'&&(llnunu_l1_pt>30)&&(llnunu_l2_pt>30)'
elif "metzptCR" in tag:cuts=cuts_zmass+'&&(llnunu_l1_pt>60)'
elif "metzpt100" in tag:cuts=cuts_zmass+'&&(llnunu_l1_pt>100)&&(llnunu_l2_pt>100)'
elif "zveto" in tag:cuts=cuts_zmass
elif "full" in tag:cuts=cuts_zmassin
elif "sigall" in tag: cuts="llnunu_l1_mass>70&&llnunu_l1_mass<110"
elif "sigzpt100" in tag: cuts="llnunu_l1_mass>70&&llnunu_l1_mass<110&&(llnunu_l1_pt>100)"
elif "sigpt100" in tag: cuts="llnunu_l1_mass>70&&llnunu_l1_mass<110&&(llnunu_l1_pt>100)&&(llnunu_l2_pt>100)"
else:
    raise RuntimeError, "tag not defined"
if channel=='el': cuts = cuts+'&&'+elChannel
elif channel=='mu': cuts = cuts+'&&'+muChannel
else: cuts = cuts+'&&({0}||{1})'.format(elChannel,muChannel)
if UseMETFilter:
    cuts = '('+cuts+')'

ROOT.gROOT.ProcessLine('.x tdrstyle.C') 

vvPlotters=[]
vvSamples = ['WZTo2L2Q','WZTo3LNu',
'ZZTo2L2Nu',
'ZZTo2L2Q','ZZTo4L',
'ggZZTo2e2nu','ggZZTo2mu2nu',
'TTZToLLNuNu']
for sample in vvSamples:
    vvPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    vvPlotters[-1].addCorrectionFactor('1/SumWeights','norm')
    vvPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    vvPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    vvPlotters[-1].addCorrectionFactor(lepsf, 'lepsf')
    vvPlotters[-1].addCorrectionFactor('xsec','xsec')
    if sample == 'ZZTo2L2Nu' : vvPlotters[-1].addCorrectionFactor("(ZZEwkCorrWeight*ZZQcdCorrWeight)", 'nnlo')
    vvPlotters[-1].setAlias('passMuHLT', '((llnunu_l1_l1_trigerob_HLTbit>>3&1)||(llnunu_l1_l1_trigerob_HLTbit>>4&1)||(llnunu_l1_l2_trigerob_HLTbit>>3&1)||(llnunu_l1_l2_trigerob_HLTbit>>4&1))');
    vvPlotters[-1].setAlias('passElHLT', '((llnunu_l1_l1_trigerob_HLTbit>>1&1)||(llnunu_l1_l2_trigerob_HLTbit>>1&1))');
    vvPlotters[-1].addCorrectionFactor('(passMuHLT||passElHLT)','HLT')
VV = MergedPlotter(vvPlotters)
VV.setFillProperties(1001,ROOT.kMagenta)

nonresPlotters=[]
nonresSamples = ['muonegtree_light_skim_38']
for sample in nonresSamples:
    nonresPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    nonresPlotters[-1].addCorrectionFactor('mtrgsf', 'mtrgsf')
    nonresPlotters[-1].addCorrectionFactor(nonresscale, 'norm')
    nonresPlotters[-1].addCorrectionFactor('mscale', 'mscale')
    nonresPlotters[-1].addCorrectionFactor('1./({0}*1000)'.format(lumi), 'lumi')
NONRES = MergedPlotter(nonresPlotters)
NONRES.setFillProperties(1001,ROOT.kOrange)

zjetsPlotters=[]
zjetsSamples = ['DYJetsToLL_M50_MGMLM_BIG',]
for sample in zjetsSamples:
    zjetsPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    zjetsPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    if ZJetsZPtWeight: zjetsPlotters[-1].addCorrectionFactor('ZPtWeight','ZPtWeight')
    zjetsPlotters[-1].addCorrectionFactor('xsec','xsec')
    zjetsPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    zjetsPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    zjetsPlotters[-1].addCorrectionFactor(lepsf,'lepsf')
    zjetsPlotters[-1].setAlias('passMuHLT', '((llnunu_l1_l1_trigerob_HLTbit>>3&1)||(llnunu_l1_l1_trigerob_HLTbit>>4&1)||(llnunu_l1_l2_trigerob_HLTbit>>3&1)||(llnunu_l1_l2_trigerob_HLTbit>>4&1))');
    zjetsPlotters[-1].setAlias('passElHLT', '((llnunu_l1_l1_trigerob_HLTbit>>1&1)||(llnunu_l1_l2_trigerob_HLTbit>>1&1))');
    zjetsPlotters[-1].addCorrectionFactor('(passMuHLT||passElHLT)','HLT') 
ZJets = MergedPlotter(zjetsPlotters)
ZJets.setFillProperties(1001,ROOT.kGreen+2)

dataPlotters=[]
dataSamples = [
'SingleEMU_Run2016Full_ReReco_v1_DtReCalib'
]
for sample in dataSamples:
    dataPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    dataPlotters[-1].setAlias('passMuHLT', '((llnunu_l1_l1_trigerob_HLTbit>>3&1)||(llnunu_l1_l1_trigerob_HLTbit>>4&1)||(llnunu_l1_l2_trigerob_HLTbit>>3&1)||(llnunu_l1_l2_trigerob_HLTbit>>4&1))');
    dataPlotters[-1].setAlias('passElHLT', '((llnunu_l1_l1_trigerob_HLTbit>>1&1)||(llnunu_l1_l2_trigerob_HLTbit>>1&1))');
    dataPlotters[-1].addCorrectionFactor('(passMuHLT||passElHLT)','HLT')
Data = MergedPlotter(dataPlotters)



Stack = StackPlotter(outTag=tag, outDir=outdir)
Stack.setPaveText(paveText)
Stack.addPlotter(Data, "data_obs", "Data", "data")
#Stack.addPlotter(WJets, "WJets","W+Jets", "background")
Stack.addPlotter(NONRES, "NONReso","non reson.", "background")
Stack.addPlotter(VV, "VVZReso","ZZ WZ reson.", "background")
#Stack.addPlotter(ggZZ, "ggZZ","ggZZ", "background")
Stack.addPlotter(ZJets, "ZJets","ZJets", "background")

Stack.setLog(LogY)
Stack.doRatio(doRatio)



tag+='_'
print cuts
if test: 
    Stack.drawStack('llnunu_l1_mass', cuts, str(lumi*1000), 65, 50, 180, titlex = "M(Z)", units = "GeV",output=tag+'zmass',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 30, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt_low',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_mt', cuts, str(lumi*1000), 40, 0.0, 2000.0, titlex = "M_{T}", units = "GeV",output=tag+'mt_high3',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=300)
    Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 30, 0, 1500, titlex = "MET", units = "GeV",output=tag+'met_low2',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)

Stack.closePSFile()
