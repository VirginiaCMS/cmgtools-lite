#!/usr/bin/env python

import ROOT
import os,sys, string, math, pickle
from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter
ROOT.gROOT.SetBatch()
if len(sys.argv)!=2:
    sys.exit()
tag="nonres_{0}_".format(sys.argv[1])
nonresscale='0.371917507815'
channel='el'
LogY=("_log_" in tag)
Blind=("_blind_" in tag)
test=True
DrawLeptons=True
doRhoScale=True

if test: DrawLeptons = False

lepsf="trgsf*isosf*idsf*trksf"

if doRhoScale: 
    tag+="RhoWt_"
    lepsf=lepsf+"*(0.366*TMath::Gaus(rho,8.280,5.427)+0.939*TMath::Gaus(rho,18.641,10.001)+0.644*TMath::Gaus(rho,40.041,10.050))" # 2016 rereco/summer16 81.81 fb-1

outdir='plots'

indir='/home/heli/XZZ/80X_20170202_light_Skim/'

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
elif "metzptCR" in tag:cuts=cuts_zmass+'&&(llnunu_l1_pt>50)'
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

dataPlotters=[]
dataSamples = [
    'muonegtree_light_skim_38_skim'
]
for sample in dataSamples:
    dataPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    dataPlotters[-1].addCorrectionFactor('etrgsf', 'etrgsf')
    dataPlotters[-1].addCorrectionFactor('escale', 'escale')
    dataPlotters[-1].addCorrectionFactor(nonresscale, 'norm')
Data = MergedPlotter(dataPlotters)

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

Stack = StackPlotter(outTag=tag, outDir=outdir)
Stack.setPaveText(paveText)
Stack.addPlotter(Data, "data_obs", "data-driven nonreson", "data")
Stack.addPlotter(WW, "NonReso","WW/WZ/WJets non-reson.", "background")
Stack.addPlotter(TT, "TT","TT", "background")

Stack.setLog(LogY)
Stack.doRatio(doRatio)

tag+='_'
print cuts
if test: 
    Stack.drawStack('llnunu_l1_mass', cuts, str(lumi*1000), 65, 50, 180, titlex = "M(Z)", units = "GeV",output=tag+'zmass',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 30, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt_low',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_mt', cuts, str(lumi*1000), 40, 0.0, 2000.0, titlex = "M_{T}", units = "GeV",output=tag+'mt_high3',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=300)
    Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 30, 0, 1500, titlex = "MET", units = "GeV",output=tag+'met_low2',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
    Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 30, 0.0, 200.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt_low',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_pt*cos(llnunu_l2_phi-llnunu_l1_phi)', cuts, str(lumi*1000), 50, -500, 500.0, titlex = "MET_{#parallel}", units = "GeV",output=tag+'met_para',outDir=outdir,separateSignal=sepSig)
Stack.closePSFile()
