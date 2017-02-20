#!/usr/bin/env python

import ROOT
import os,sys, string, math, pickle
from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter
ROOT.gROOT.SetBatch()
if len(sys.argv)!=2:
    sys.exit()
tag="subt_{0}_".format(sys.argv[1])
#nonresscale='0.700104002035'
nonresscale='0.688794241579'
channel='mu'
LogY=("_log_" in tag)
Blind=("_blind_" in tag)
dtnonres=("_dtnonres_" in tag)
test=True
DrawLeptons=True

if test: DrawLeptons = False

lepsf="trgsf*isosf*idsf*trksf"

outdir='plots'

indir='/home/heli/XZZ/80X_20170202_light_Skim/'

lumi=36.814
sepSig=True
doRatio=True
FakeData=False
UseMETFilter=True

puWeight='puWeightsummer16'


k=1 # signal scale

elChannel='((abs(llnunu_l1_l1_pdgId)==11||abs(llnunu_l1_l2_pdgId)==11)&&llnunu_l1_l1_pt>120&&abs(llnunu_l1_l1_eta)<2.5&&llnunu_l1_l2_pt>35&&abs(llnunu_l1_l2_eta)<2.5)'
muChannel='((abs(llnunu_l1_l1_pdgId)==13||abs(llnunu_l1_l2_pdgId)==13)&&llnunu_l1_l1_pt>55&&abs(llnunu_l1_l1_eta)<2.4&&llnunu_l1_l2_pt>20&&abs(llnunu_l1_l2_eta)<2.4&&(llnunu_l1_l1_highPtID>0.99||llnunu_l1_l2_highPtID>0.99))'

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
#elif "metzptCR" in tag:cuts=cuts_zmass+'&&(llnunu_l1_pt>60)'
elif "metzptCR" in tag:cuts=cuts_zmass
elif "metzpt100" in tag:cuts=cuts_zmass+'&&(llnunu_l1_pt>100)&&(llnunu_l2_pt>100)'
elif "zveto" in tag:cuts=cuts_zmass
elif "full" in tag:cuts=cuts_zmassin
elif "sigall" in tag: cuts="llnunu_l1_mass>70&&llnunu_l1_mass<110"
elif "sigzpt100" in tag: cuts="llnunu_l1_mass>70&&llnunu_l1_mass<110&&(llnunu_l1_pt>100)"
elif "sigzpt50" in tag: cuts="llnunu_l1_mass>70&&llnunu_l1_mass<110&&(llnunu_l1_pt>50)"
elif "sigpt100" in tag: cuts="llnunu_l1_mass>70&&llnunu_l1_mass<110&&(llnunu_l1_pt>100)&&(llnunu_l2_pt>50)"
else:
    raise RuntimeError, "tag not defined"
if channel=='el': cuts = cuts+'&&'+elChannel
elif channel=='mu': cuts = cuts+'&&'+muChannel
else: cuts = cuts+'&&({0}||{1})'.format(elChannel,muChannel)
if UseMETFilter:
    cuts = '('+cuts+')'

ROOT.gROOT.ProcessLine('.x tdrstyle.C') 

wwSamples = ['WWTo2L2Nu','WWToLNuQQ_BIG','WZTo1L1Nu2Q']
wwPlotters=[]
for sample in wwSamples:
    wwPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    wwPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    wwPlotters[-1].addCorrectionFactor('xsec','xsec')
    wwPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    wwPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    wwPlotters[-1].addCorrectionFactor(lepsf,'lepsf')
    wwPlotters[-1].setAlias('passMuHLT', '((llnunu_l1_l1_trigerob_HLTbit>>3&1)||(llnunu_l1_l1_trigerob_HLTbit>>4&1)||(llnunu_l1_l2_trigerob_HLTbit>>3&1)||(llnunu_l1_l2_trigerob_HLTbit>>4&1))');
    wwPlotters[-1].setAlias('passElHLT', '((llnunu_l1_l1_trigerob_HLTbit>>1&1)||(llnunu_l1_l2_trigerob_HLTbit>>1&1))');
    wwPlotters[-1].addCorrectionFactor('(passMuHLT||passElHLT)','HLT')
WW = MergedPlotter(wwPlotters)
WW.setFillProperties(1001,ROOT.kOrange)

ttSamples = ['TTTo2L2Nu_noSC','TTWJetsToLNu_BIG', 'T_tWch', 'T_tch_powheg', 'TBar_tWch', 'TBar_tch_powheg']
ttPlotters=[]
for sample in ttSamples:
    ttPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    ttPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    ttPlotters[-1].addCorrectionFactor('xsec','xsec')
    ttPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    ttPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    ttPlotters[-1].addCorrectionFactor(lepsf,'lepsf')
    ttPlotters[-1].setAlias('passMuHLT', '((llnunu_l1_l1_trigerob_HLTbit>>3&1)||(llnunu_l1_l1_trigerob_HLTbit>>4&1)||(llnunu_l1_l2_trigerob_HLTbit>>3&1)||(llnunu_l1_l2_trigerob_HLTbit>>4&1))');
    ttPlotters[-1].setAlias('passElHLT', '((llnunu_l1_l1_trigerob_HLTbit>>1&1)||(llnunu_l1_l2_trigerob_HLTbit>>1&1))');
    ttPlotters[-1].addCorrectionFactor('(passMuHLT||passElHLT)','HLT')
TT = MergedPlotter(ttPlotters)
TT.setFillProperties(1001,ROOT.kAzure-9)

nonresPlotters=[]
nonresSamples = ['muonegtree_light_skim_38_skim']
for sample in nonresSamples:
    nonresPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    nonresPlotters[-1].addCorrectionFactor('mtrgsf', 'mtrgsf')
    nonresPlotters[-1].addCorrectionFactor(nonresscale, 'norm')
    nonresPlotters[-1].addCorrectionFactor('mscale', 'mscale')
    nonresPlotters[-1].addCorrectionFactor('1./({0}*1000)'.format(lumi), 'lumi')
NONRES = MergedPlotter(nonresPlotters)
NONRES.setFillProperties(1001,ROOT.kOrange)

dataPlotters=[]
dataSamples = [
'SingleEMU_Run2016Full_ReReco_v2_DtReCalib',
'WZTo2L2Q','WZTo3LNu',
'ZZTo2L2Nu',
'ZZTo2L2Q','ZZTo4L',
'ggZZTo2e2nu','ggZZTo2mu2nu',
'TTZToLLNuNu',
'DYJetsToLL_M50_MGMLM_BIG'
]
for sample in dataSamples:
    dataPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    dataPlotters[-1].setAlias('passMuHLT', '((llnunu_l1_l1_trigerob_HLTbit>>3&1)||(llnunu_l1_l1_trigerob_HLTbit>>4&1)||(llnunu_l1_l2_trigerob_HLTbit>>3&1)||(llnunu_l1_l2_trigerob_HLTbit>>4&1))');
    dataPlotters[-1].setAlias('passElHLT', '((llnunu_l1_l1_trigerob_HLTbit>>1&1)||(llnunu_l1_l2_trigerob_HLTbit>>1&1))');
    dataPlotters[-1].addCorrectionFactor('(passMuHLT||passElHLT)','HLT')
    if 'SingleEMU' not in sample:
        dataPlotters[-1].addCorrectionFactor('(-{0}*1000)'.format(lumi), 'lumi')
        dataPlotters[-1].addCorrectionFactor('1/SumWeights','norm')
        dataPlotters[-1].addCorrectionFactor('genWeight','genWeight')
        dataPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
        dataPlotters[-1].addCorrectionFactor(lepsf, 'lepsf')
        if sample == 'WZTo3LNu': dataPlotters[-1].addCorrectionFactor('4.4297','xsec') 
        else: dataPlotters[-1].addCorrectionFactor('xsec','xsec')
        if sample == 'ZZTo2L2Nu' : dataPlotters[-1].addCorrectionFactor("(ZZEwkCorrWeight*ZZQcdCorrWeight)", 'nnlo')
#        if sample=='DYJetsToLL_M50_MGMLM_BIG': dataPlotters[-1].addCorrectionFactor('ZPtWeight*1.02509','ZPtWeight')
        if sample=='DYJetsToLL_M50_MGMLM_BIG': dataPlotters[-1].addCorrectionFactor('ZPtWeight*1.24','ZPtWeight')
Data = MergedPlotter(dataPlotters)



Stack = StackPlotter(outTag=tag, outDir=outdir)
Stack.setPaveText(paveText)
Stack.addPlotter(Data, "data_obs", "Data", "data")
if dtnonres:Stack.addPlotter(NONRES, "NONReso","non reson.", "background")
else:
    Stack.addPlotter(WW, "WW","WW/WZ/WJets non-reson.", "background")
    Stack.addPlotter(TT, "TT","TT", "background")


Stack.setLog(LogY)
Stack.doRatio(doRatio)



tag+='_'
print cuts
if test: 
    Stack.drawStack('llnunu_l1_mass', cuts, str(lumi*1000), 65, 50, 180, titlex = "M(Z)", units = "GeV",output=tag+'zmass',outDir=outdir,separateSignal=sepSig)
#    Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 30, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt_low',outDir=outdir,separateSignal=sepSig)
#    Stack.drawStack('llnunu_mt', cuts, str(lumi*1000), 40, 0.0, 2000.0, titlex = "M_{T}", units = "GeV",output=tag+'mt_high3',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=300)
#    Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 30, 0, 1500, titlex = "MET", units = "GeV",output=tag+'met_low2',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
#    Stack.drawStack('llnunu_l1_l1_pt', cuts, str(lumi*1000), 20, 50, 250, titlex = "Pt_{l1}", units = "GeV",output=tag+'pt_l1',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
#    Stack.drawStack('llnunu_l1_l2_pt', cuts, str(lumi*1000), 20, 20, 220, titlex = "Pt_{l2}", units = "GeV",output=tag+'pt_l2',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)

Stack.closePSFile()
