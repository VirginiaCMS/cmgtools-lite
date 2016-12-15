#!/usr/bin/env python

import optparse
import ROOT
import os,sys, string, math, pickle
from CMGTools.XZZ2l2nu.plotting.TreePlotter import TreePlotter
from CMGTools.XZZ2l2nu.plotting.MergedPlotter import MergedPlotter
from CMGTools.XZZ2l2nu.plotting.StackPlotter import StackPlotter

grootargs = []
def callback_rootargs(option, opt, value, parser):
    grootargs.append(opt)

parser = optparse.OptionParser()
parser.add_option("-t","--tag",dest="tag",default='DataB2G_ICHEPcfg_',help="")
parser.add_option("--cutChain",dest="cutChain",default='tight',help="")
parser.add_option("--channel",dest="channel",default='mu',help="only works with --ZWeight")
parser.add_option("--LogY",action="store_true", dest="LogY", default=False, help="")
parser.add_option("--Blind",action="store_true", dest="Blind", default=False,help="")
parser.add_option("--test",action="store_true", dest="test", default=False,help="")
parser.add_option("--ZWeight",action="store_true", dest="ZWeight", default=False,help="")
parser.add_option("--SeparateProcess",action="store_true", dest="SeparateProcess", default=False,help="")

parser.add_option("-l",action="callback",callback=callback_rootargs)
parser.add_option("-q",action="callback",callback=callback_rootargs)
parser.add_option("-b",action="callback",callback=callback_rootargs)





(options,args) = parser.parse_args()


tag=options.tag
cutChain=options.cutChain

channel=options.channel
ZWeight=options.ZWeight
LogY=options.LogY
test=options.test
SeparateProcess=options.SeparateProcess
doRhoScale=True
doVtxScale=False
doEtaScale=True
doPhPtScale=True



if SeparateProcess: tag +='SepProc_'

scale='(1)'

if doRhoScale: 
    tag+="RhoWt_"
    scale=scale+"*(0.32+0.42*TMath::Erf((rho-4.16)/4.58)+0.31*TMath::Erf((rho+115.00)/29.58))" # b2h rereco 36.1 fb-1

if doVtxScale:
    tag+="VtxWt_"
    scale=scale+"*(0.807+0.007*nVert+-3.689e-05*nVert*nVert+6.730e-04*exp(2.500e-01*nVert))" # b2h rereco 33.59fb-1

if doEtaScale:
    tag+="EtaWt_"
    scale=scale+"*(0.87*TMath::Gaus(llnunu_l1_eta,0.65,0.56)+0.87*TMath::Gaus(llnunu_l1_eta,-0.65,0.56)+0.65*TMath::Gaus(llnunu_l1_eta,1.90,0.25)+0.65*TMath::Gaus(llnunu_l1_eta,-1.90,0.25))"


if doPhPtScale:
    tag+="PhPtWt_"
    scale=scale+"*((-1.06624+0.0580113*pow(llnunu_l1_pt,1)-5.09328e-4*pow(llnunu_l1_pt,2)+2.28513e-6*pow(llnunu_l1_pt,3)-6.03131e-9*pow(llnunu_l1_pt,4)+9.84946e-12*pow(llnunu_l1_pt,5)-1.00558e-14*pow(llnunu_l1_pt,6)+6.244e-18*pow(llnunu_l1_pt,7)-2.15543e-21*pow(llnunu_l1_pt,8)+3.17021e-25*pow(llnunu_l1_pt,9))*(llnunu_l1_pt<=1000)+(0.688060)*(llnunu_l1_pt>1000))"


outdir='plots_ph_36p46'

indir='/home/heli/XZZ/80X_20161029_GJets_light_Skim'
lumi=36.46
sepSig=True
doRatio=True
Blind=options.Blind
puWeight='puWeightmoriondMC'

if not os.path.exists(outdir): os.system('mkdir -p '+outdir)

tag = tag+cutChain+'_'
tag = tag+puWeight+'_'


if not Blind: tag = tag+'unblind_'

if LogY: tag = tag+'log_'

if ZWeight:
    tag +='ZWeight_'
    if   channel=='el': tag+='el'
    elif channel=='mu': tag+='mu'
    else:  tag+='all'

paveText="#sqrt{s} = 13 TeV 2016 L = "+"{:.4}".format(float(lumi))+" fb^{-1}"


cuts_loose='(1)'
cuts_tight='llnunu_l1_pt>50'

if cutChain=='loosecut': cuts=cuts_loose
elif cutChain=='tight': cuts=cuts_tight
else : cuts=cuts_loose

cuts = '('+cuts+')'


ROOT.gROOT.ProcessLine('.x tdrstyle.C') 


# parameters for GJets
gdataLumi=36.46*1000
gdataYield=3451449849.011390686
gdataFidXsec=gdataYield/gdataLumi

zjetsFidXsecAll = 151.06068438939382759
zjetsFidXsecEl =  1.8318217140038339785
zjetsFidXsecMu =  149.22886267539001892
zjetsFidXsecAll_up = 151.85715853322426483
zjetsFidXsecAll_dn = 150.26421019455997907
zjetsFidXsecEl_up = 1.8728979304188486665
zjetsFidXsecEl_dn = 1.7907454975888201787
zjetsFidXsecMu_up = 149.98426060280544903
zjetsFidXsecMu_dn = 148.47346469697114912
zjetsFidXsecLowLptAll = 807.42655018368884612
zjetsFidXsecLowLptEl = 229.77648821257676559
zjetsFidXsecLowLptMu = 577.65006197098625762

############################################
#
#  Physical MET
#
############################################

################
# DYJets->LL
################
zllSamples = ['DYJetsToLL_M50_reHLT']

zllPlotters=[]
for sample in zllSamples:
    zllPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    zllPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    zllPlotters[-1].addCorrectionFactor('xsec','xsec')
    zllPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    zllPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    zllPlotters[-1].addCorrectionFactor(scale,'scale')



##################
# DYJets -> NuNu
##################
znnSamples = [
'ZJetsToNuNu_HT100to200_BIG',
'ZJetsToNuNu_HT200to400_BIG',
'ZJetsToNuNu_HT400to600_BIG',
'ZJetsToNuNu_HT600to800_BIG',
'ZJetsToNuNu_HT800t1200_BIG',
'ZJetsToNuNu_HT1200to2500_BIG',
'ZJetsToNuNu_HT2500toInf_BIG',
]

znnPlotters=[]
for sample in znnSamples:
    znnPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    znnPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    znnPlotters[-1].addCorrectionFactor('xsec','xsec')
    znnPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    znnPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    znnPlotters[-1].addCorrectionFactor(scale,'scale')



#######################
# Zgamma->nu nu gamma
#######################

znngSamples = [
'ZNuNuGJetsGt40Lt130',
'ZNuNuGJetsGt130',
]

znngPlotters=[]
for sample in znngSamples:
    znngPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    znngPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    znngPlotters[-1].addCorrectionFactor('xsec','xsec')
    znngPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    znngPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    znngPlotters[-1].addCorrectionFactor(scale,'scale')



###########################
# Wgamma-> l nu gamma
###########################
wlngSamples = [
'WGToLNuG',
]

wlngPlotters=[]
for sample in wlngSamples:
    wlngPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    wlngPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    wlngPlotters[-1].addCorrectionFactor('xsec','xsec')
    wlngPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    wlngPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    wlngPlotters[-1].addCorrectionFactor(scale,'scale')



###############################
# WJets->lnu
###############################

wlnSamples = [
#'WJetsToLNu',
'WJetsToLNu_HT100to200_BIG',
'WJetsToLNu_HT1200to2500_BIG',
'WJetsToLNu_HT200to400_BIG',
'WJetsToLNu_HT2500toInf_BIG',
'WJetsToLNu_HT400to600_BIG',
'WJetsToLNu_HT600to800_BIG',
'WJetsToLNu_HT800to1200_BIG',
]

wlnPlotters=[]
for sample in wlnSamples:
    wlnPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    wlnPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    wlnPlotters[-1].addCorrectionFactor('xsec','xsec')
    wlnPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    wlnPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    wlnPlotters[-1].addCorrectionFactor(scale,'scale')


################################
# single top, ttbar
################################

tSamples = [
'TToLeptons_tch_powheg',
'TBarToLeptons_tch_powheg',
'T_tWch',
'TBar_tWch',
'TGJets_BIG',
'TTGJets',
]

tPlotters=[]
for sample in tSamples:
    tPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    tPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    tPlotters[-1].addCorrectionFactor('xsec','xsec')
    tPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    tPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    tPlotters[-1].addCorrectionFactor(scale,'scale')


#########################################
# all above  physical met adding together
##########################################

phymetPlotters = zllPlotters+znnPlotters+znngPlotters+wlngPlotters+wlnPlotters+tPlotters




############################################
#
#  Instrumental MET
#
############################################

##########################
# GJets
##########################

gjetsSamples = [
#'GJets_HT40to100',
'GJets_HT100to200',
'GJets_HT200to400',
'GJets_HT400to600',
'GJets_HT600toInf',
]

gjetsPlotters=[]
for sample in gjetsSamples:
    gjetsPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    gjetsPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    gjetsPlotters[-1].addCorrectionFactor('xsec','xsec')
    gjetsPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    gjetsPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    gjetsPlotters[-1].addCorrectionFactor(scale,'scale')



##########################
# QCD
##########################

qcdSamples = [
'QCD_Pt20to30_EMEnriched',
'QCD_Pt30to50_EMEnriched',
'QCD_Pt50to80_EMEnriched',
'QCD_Pt80to120_EMEnriched',
'QCD_Pt120to170_EMEnriched',
'QCD_Pt170to300_EMEnriched',
'QCD_Pt300toInf_EMEnriched',
]

qcdPlotters=[]
for sample in qcdSamples:
    qcdPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    qcdPlotters[-1].addCorrectionFactor('1./SumWeights','norm')
    qcdPlotters[-1].addCorrectionFactor('xsec','xsec')
    qcdPlotters[-1].addCorrectionFactor('genWeight','genWeight')
    qcdPlotters[-1].addCorrectionFactor(puWeight,'puWeight')
    qcdPlotters[-1].addCorrectionFactor(scale,'scale')




##############################################
# all above instrumental met adding together
##############################################

insmetPlotters = gjetsPlotters+qcdPlotters



############################################
#
# All above Instrumental + Physical MET
#
############################################

allmcPlotters = phymetPlotters+insmetPlotters

for i in range(len(allmcPlotters)) :
    if ZWeight:
        allmcPlotters[i].addCorrectionFactor(str(1/gdataFidXsec),'frac') # divided by g data fid-xsec
        allmcPlotters[i].addCorrectionFactor('GJetsRhoWeight','GJetsRhoWeight')
        if channel=='el' :
            allmcPlotters[i].addCorrectionFactor('GJetsZPtWeightEl','GJetsZPtWeight')
            allmcPlotters[i].addCorrectionFactor(str(zjetsFidXsecEl),'zjetsFidXsecEl')
            allmcPlotters[i].addCorrectionFactor('(1.16725)','GJetsNorm')
        elif channel=='mu' :
            allmcPlotters[i].addCorrectionFactor('GJetsZPtWeightMu','GJetsZPtWeight')
            allmcPlotters[i].addCorrectionFactor(str(zjetsFidXsecMu),'zjetsFidXsecMu')
            allmcPlotters[i].addCorrectionFactor('(1.15559)','GJetsNorm')
        else :
            allmcPlotters[i].addCorrectionFactor('GJetsZPtWeight','GJetsZPtWeight')
            allmcPlotters[i].addCorrectionFactor(str(zjetsFidXsecAll),'zjetsFidXsecAll')
            allmcPlotters[i].addCorrectionFactor('(1.15559)','GJetsNorm')






############################################
#
# gamma  Data 
#
############################################

gdataSamples = [
#'SinglePhoton_Run2016B2H_ReReco_36p46_ResBos_Rc36p46ReCalib',
#'SinglePhoton_Run2016B2H_ReReco_36p46_ResBosRefit_Rc36p46ReCalib',
'SinglePhoton_Run2016B2H_ReReco_36p46_Rc36p46ReCalib',
]

gdataPlotters=[]
for sample in gdataSamples:
    gdataPlotters.append(TreePlotter(sample, indir+'/'+sample+'.root','tree'))
    gdataPlotters[-1].addCorrectionFactor('GJetsPreScaleWeight','prescale')
    if ZWeight:
        gdataPlotters[-1].addCorrectionFactor(str(lumi*1000),'GJetsLumi')
        gdataPlotters[-1].addCorrectionFactor(str(1/gdataYield),'GJetsNorm0')
        gdataPlotters[-1].addCorrectionFactor('GJetsRhoWeight','GJetsRhoWeight')
        if channel=='el' :
            gdataPlotters[-1].addCorrectionFactor('GJetsZPtWeightEl','GJetsZPtWeight')
            gdataPlotters[-1].addCorrectionFactor(str(zjetsFidXsecEl),'zjetsFidXsecEl')
            gdataPlotters[-1].addCorrectionFactor('(1.16725)','GJetsNorm')
        elif channel=='mu' :
            gdataPlotters[-1].addCorrectionFactor('GJetsZPtWeightMu','GJetsZPtWeight')
            gdataPlotters[-1].addCorrectionFactor(str(zjetsFidXsecMu),'zjetsFidXsecMu')
            gdataPlotters[-1].addCorrectionFactor('(1.15559)','GJetsNorm')
        else :
            gdataPlotters[-1].addCorrectionFactor('GJetsZPtWeight','GJetsZPtWeight')
            gdataPlotters[-1].addCorrectionFactor(str(zjetsFidXsecAll),'zjetsFidXsecAll')
            gdataPlotters[-1].addCorrectionFactor('(1.15559)','GJetsNorm')








############################################
#
#  Create Plotters
#
############################################


############################################
# Create MergedPlotters
############################################

#########################
# Physical MET
#########################

ZLL = MergedPlotter(zllPlotters)
ZLL.setFillProperties(1001,ROOT.kOrange)
ZNN = MergedPlotter(znnPlotters)
ZNN.setFillProperties(1001,ROOT.kMagenta)
ZNNG = MergedPlotter(znngPlotters)
ZNNG.setFillProperties(1001,ROOT.kMagenta+2)
WLNG = MergedPlotter(wlngPlotters)
WLNG.setFillProperties(1001,ROOT.kRed)
WLN = MergedPlotter(wlnPlotters)
WLN.setFillProperties(1001,ROOT.kYellow)
T = MergedPlotter(tPlotters)
T.setFillProperties(1001,ROOT.kAzure-9)

#########################
# All Physical MET
#########################

PhyMET = MergedPlotter(phymetPlotters)
PhyMET.setFillProperties(1001, ROOT.kOrange)


#########################
# Instrumental MET
#########################

GJETS = MergedPlotter(gjetsPlotters)
GJETS.setFillProperties(1001,ROOT.kBlue-6)
QCD = MergedPlotter(qcdPlotters)
QCD.setFillProperties(1001,ROOT.kGreen+2)

#########################
# All Instrumental MET
#########################

InsMET = MergedPlotter(insmetPlotters)
InsMET.setFillProperties(1001, ROOT.kGreen)



########################
# gamma Data
########################
GData = MergedPlotter(gdataPlotters)



############################################
# Create StackPlotter
############################################

Stack = StackPlotter(outTag=tag, outDir=outdir)
Stack.setPaveText(paveText)
Stack.addPlotter(GData, "data_obs", "#gamma Data", "data")
if not SeparateProcess:
    Stack.addPlotter(PhyMET, "PhyMET","Physical MET", "background")
    Stack.addPlotter(InsMET, "InsMET","Instrumental MET", "background")
else:
    Stack.addPlotter(ZNN, "ZNN","Z->#nu#nu", "background")
    Stack.addPlotter(ZNNG, "ZNNG","Z#gamma->#nu#nu#gamma", "background")
    Stack.addPlotter(T, "T","Top", "background")
    Stack.addPlotter(ZLL, "ZLL","Z->ll", "background")
    Stack.addPlotter(WLNG, "WGToLNuG","W#gamma->l#nu#gamma", "background")
    Stack.addPlotter(WLN, "WLN","W->l#nu", "background")
    Stack.addPlotter(GJETS, "GJETS","#gamma+jets", "background")
    Stack.addPlotter(QCD, "QCD","QCD", "background")

 
Stack.setLog(LogY)
Stack.doRatio(doRatio)



tag+='_'

if test: 
    Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 75, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt_low',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 25, 0, 1000, titlex = "MET", units = "GeV",output=tag+'met_low',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
else: 
    Stack.drawStack('nVert', cuts, str(lumi*1000), 100, 0.0, 100.0, titlex = "N vertices", units = "",output=tag+'nVert',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('rho', cuts, str(lumi*1000), 55, 0.0, 55.0, titlex = "#rho", units = "",output=tag+'rho',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_pt', cuts, str(lumi*1000), 75, 0.0, 1500.0, titlex = "P_{T}(Z)", units = "GeV",output=tag+'zpt_low',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_eta', cuts, str(lumi*1000), 100, -2.5, 2.5, titlex = "#eta(Z) ", units = "",output=tag+'zeta',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l1_phi', cuts, str(lumi*1000), 64, -3.2, 3.2, titlex = "#phi(Z)", units = "",output=tag+'zphi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_pt', cuts, str(lumi*1000), 25, 0, 1000, titlex = "MET", units = "GeV",output=tag+'met_low',outDir=outdir,separateSignal=sepSig,blinding=Blind,blindingCut=200)
    Stack.drawStack('llnunu_l2_pt*cos(llnunu_l2_phi-llnunu_l1_phi)', cuts, str(lumi*1000), 100, -500, 500.0, titlex = "MET_{#parallel}", units = "GeV",output=tag+'met_para',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_pt*sin(llnunu_l2_phi-llnunu_l1_phi)', cuts, str(lumi*1000), 100, -500, 500.0, titlex = "MET_{#perp}", units = "GeV",output=tag+'met_perp',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_phi', cuts, str(lumi*1000), 100, -3.2, 3.2, titlex = "#phi(MET)", units = "",output=tag+'metPhi',outDir=outdir,separateSignal=sepSig)
    Stack.drawStack('llnunu_l2_sumEt', cuts, str(lumi*1000), 80, 0.0, 3000.0, titlex = "sumE_{T}", units = "GeV",output=tag+'metSumEt',outDir=outdir,separateSignal=sepSig)


Stack.closePSFile()
Stack.closeROOTFile()

