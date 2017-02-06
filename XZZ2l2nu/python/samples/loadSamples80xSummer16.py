#####################
# load samples 
#####################

import PhysicsTools.HeppyCore.framework.config as cfg
import os

# Load backgrounds
from CMGTools.XZZ2l2nu.samples.samples_13TeV_RunIISummer16MiniAODv2 import *
# Load signals
from CMGTools.XZZ2l2nu.samples.samples_13TeV_signal80XSummer16 import *
# Load Data 
from CMGTools.XZZ2l2nu.samples.samples_13TeV_DATA2016 import *
# Load triggers
from CMGTools.XZZ2l2nu.samples.triggers_13TeV_Spring16 import *

# backgrounds
backgroundSamples = [
#DYJetsToLL_M50,
DYJetsToLL_M50_Ext,
#DYJetsToLL_M50_MGMLM,
DYJetsToLL_M50_MGMLM_Ext1,
WJetsToLNu,
#TTTo2L2Nu,
TTTo2L2Nu_noSC,
TTTo2L2Nu_forTTH,
ZZTo2L2Nu,
ZZTo4L,
ZZTo2L2Q,
WWTo2L2Nu,
WWToLNuQQ,
WWToLNuQQ_Ext1,
WZTo1L1Nu2Q,
WZTo2L2Q,
WZTo3LNu,
TTZToLLNuNu,
#TTWJetsToLNu,
TTWJetsToLNu_ext1,
TTWJetsToLNu_ext2,
ggZZTo2e2nu,
ggZZTo2mu2nu,
]

extraBackgroundMC = [
DY1JetsToLL_M50_MGMLM,
DY2JetsToLL_M50_MGMLM,
DY3JetsToLL_M50_MGMLM,
DY4JetsToLL_M50_MGMLM,
#DYBJetsToLL_M50_MGMLM,
#DYJetsToLL_M50_Ext,
]



GJetsMC= WJetsToLNuHT + ZJetsToNuNuHT + SingleTop + ZNuNuGJets + ttgjets + wgjets
GJetsMC +=GJet_Pt_EMEnriched + GJetsHT
GJetsMC += QCDPtEMEnriched + QCDPt + QCDHT



# signals
signalSamples = [
BulkGravToZZToZlepZinv_narrow_600,
BulkGravToZZToZlepZinv_narrow_800,
BulkGravToZZToZlepZinv_narrow_1000,
BulkGravToZZToZlepZinv_narrow_1200,
BulkGravToZZToZlepZinv_narrow_1400,
BulkGravToZZToZlepZinv_narrow_1600,
BulkGravToZZToZlepZinv_narrow_1800,
BulkGravToZZToZlepZinv_narrow_2000,
BulkGravToZZToZlepZinv_narrow_2500,
BulkGravToZZToZlepZinv_narrow_3000,
BulkGravToZZToZlepZinv_narrow_3500,
BulkGravToZZToZlepZinv_narrow_4000,
BulkGravToZZToZlepZinv_narrow_4500,
]

BulkGravToZZToZlepZhad = [
BulkGravToZZToZlepZhad_narrow_600,
BulkGravToZZToZlepZhad_narrow_800,
BulkGravToZZToZlepZhad_narrow_1000,
BulkGravToZZToZlepZhad_narrow_1200,
BulkGravToZZToZlepZhad_narrow_1400,
BulkGravToZZToZlepZhad_narrow_1600,
BulkGravToZZToZlepZhad_narrow_1800,
BulkGravToZZToZlepZhad_narrow_2000,
BulkGravToZZToZlepZhad_narrow_2500,
BulkGravToZZToZlepZhad_narrow_3000,
BulkGravToZZToZlepZhad_narrow_3500,
BulkGravToZZToZlepZhad_narrow_4000,
BulkGravToZZToZlepZhad_narrow_4500,
]

# MC samples
otherMcSamples = BulkGravToZZToZlepZhad + ggHZZ2L2Nu 
mcSamples = signalSamples + backgroundSamples + extraBackgroundMC + GJetsMC
#mcSamples = backgroundSamples

# data
SingleMuon=[SingleMuon_Run2016B_PromptReco,
            SingleMuon_Run2016B_PromptReco_v2,
            SingleMuon_Run2016C_PromptReco_v2,
            SingleMuon_Run2016D_PromptReco_v2,
            SingleMuon_Run2016E_PromptReco_v2,
            SingleMuon_Run2016F_PromptReco_v1,
            SingleMuon_Run2016G_PromptReco_v1,
            SingleMuon_Run2016H_PromptReco_v1,
            SingleMuon_Run2016H_PromptReco_v2,
            SingleMuon_Run2016H_PromptReco_v3,
            ]
SingleElectron=[SingleElectron_Run2016B_PromptReco,
                SingleElectron_Run2016B_PromptReco_v2,
                SingleElectron_Run2016C_PromptReco_v2,
                SingleElectron_Run2016D_PromptReco_v2,
                SingleElectron_Run2016E_PromptReco_v2,
                SingleElectron_Run2016F_PromptReco_v1,
                SingleElectron_Run2016G_PromptReco_v1,
                SingleElectron_Run2016H_PromptReco_v1,
                SingleElectron_Run2016H_PromptReco_v2,
                SingleElectron_Run2016H_PromptReco_v3,
               ]
SinglePhoton=[SinglePhoton_Run2016B_PromptReco,
              SinglePhoton_Run2016B_PromptReco_v2,
              SinglePhoton_Run2016C_PromptReco_v2,
              SinglePhoton_Run2016D_PromptReco_v2,
              SinglePhoton_Run2016E_PromptReco_v2,
              SinglePhoton_Run2016F_PromptReco_v1,
              SinglePhoton_Run2016G_PromptReco_v1,
              SinglePhoton_Run2016H_PromptReco_v1,
              SinglePhoton_Run2016H_PromptReco_v2,
              SinglePhoton_Run2016H_PromptReco_v3,
             ]
MuonEG=[MuonEG_Run2016B_PromptReco,
      MuonEG_Run2016B_PromptReco_v2,
      MuonEG_Run2016C_PromptReco_v2,
      MuonEG_Run2016D_PromptReco_v2,
      MuonEG_Run2016E_PromptReco_v2,
      MuonEG_Run2016F_PromptReco_v1,
      MuonEG_Run2016G_PromptReco_v1,
      MuonEG_Run2016H_PromptReco_v1,
      MuonEG_Run2016H_PromptReco_v2,
      MuonEG_Run2016H_PromptReco_v3,
      ]

MET= [MET_Run2016B_PromptReco,
      MET_Run2016B_PromptReco_v2,
      MET_Run2016C_PromptReco_v2,
      MET_Run2016D_PromptReco_v2,
      MET_Run2016E_PromptReco_v2,
      MET_Run2016F_PromptReco_v1,
      MET_Run2016G_PromptReco_v1,
      MET_Run2016H_PromptReco_v1,
      MET_Run2016H_PromptReco_v2,
      MET_Run2016H_PromptReco_v3,
      ]

SingleMuon23Sep2016=[
            SingleMuon_Run2016B_23Sep2016,
            SingleMuon_Run2016B_23Sep2016_v2,
            SingleMuon_Run2016C_23Sep2016,
            SingleMuon_Run2016D_23Sep2016,
            SingleMuon_Run2016E_23Sep2016,
            SingleMuon_Run2016F_23Sep2016,
            SingleMuon_Run2016G_23Sep2016,
            SingleMuon_Run2016H_PromptReco_v1,
            SingleMuon_Run2016H_PromptReco_v2,
            SingleMuon_Run2016H_PromptReco_v3,
            ]
SingleElectron23Sep2016=[
                SingleElectron_Run2016B_23Sep2016,
                SingleElectron_Run2016B_23Sep2016_v2,
                SingleElectron_Run2016C_23Sep2016,
                SingleElectron_Run2016D_23Sep2016,
                SingleElectron_Run2016E_23Sep2016,
                SingleElectron_Run2016F_23Sep2016,
                SingleElectron_Run2016G_23Sep2016,
                SingleElectron_Run2016H_PromptReco_v1,
                SingleElectron_Run2016H_PromptReco_v2,
                SingleElectron_Run2016H_PromptReco_v3,
               ]
SinglePhoton23Sep2016=[
              SinglePhoton_Run2016B_23Sep2016,
              SinglePhoton_Run2016B_23Sep2016_v2,
              SinglePhoton_Run2016C_23Sep2016,
              SinglePhoton_Run2016D_23Sep2016,
              SinglePhoton_Run2016E_23Sep2016,
              SinglePhoton_Run2016F_23Sep2016,
              SinglePhoton_Run2016G_23Sep2016,
              SinglePhoton_Run2016H_PromptReco_v1,
              SinglePhoton_Run2016H_PromptReco_v2,
              SinglePhoton_Run2016H_PromptReco_v3,
             ]
MuonEG23Sep2016=[
      MuonEG_Run2016B_23Sep2016,
      MuonEG_Run2016B_23Sep2016_v2,
      MuonEG_Run2016C_23Sep2016,
      MuonEG_Run2016D_23Sep2016,
      MuonEG_Run2016E_23Sep2016,
      MuonEG_Run2016F_23Sep2016,
      MuonEG_Run2016G_23Sep2016,
      MuonEG_Run2016H_PromptReco_v1,
      MuonEG_Run2016H_PromptReco_v2,
      MuonEG_Run2016H_PromptReco_v3,
      ]

MET23Sep2016=[
      MET_Run2016B_23Sep2016,
      MET_Run2016B_23Sep2016_v2,
      MET_Run2016C_23Sep2016,
      MET_Run2016D_23Sep2016,
      MET_Run2016E_23Sep2016,
      MET_Run2016F_23Sep2016,
      MET_Run2016G_23Sep2016,
      MET_Run2016H_PromptReco_v1,
      MET_Run2016H_PromptReco_v2,
      MET_Run2016H_PromptReco_v3,
      ]

Data23Sep2016 = SingleMuon23Sep2016+SingleElectron23Sep2016+SinglePhoton23Sep2016+MuonEG23Sep2016+MET23Sep2016



for s in SingleMuon:
    #s.triggers = triggers_1mu_noniso
    s.triggers = [] 
    s.vetoTriggers = []
for s in SingleElectron:
    #s.triggers = triggers_1e_noniso
    s.triggers = [] 
    #s.vetoTriggers = triggers_1mu_noniso
    s.vetoTriggers = []
for s in SinglePhoton:
    #s.trigers = triggers_all_photons
    s.trigers = []
    s.vetoTriggers = []

for s in MuonEG:
    s.trigers = []
    s.vetoTriggers = []

for s in MET:
    s.trigers = []
    s.vetoTriggers = []

for s in Data23Sep2016:
    s.triggers = []
    s.vetoTriggers = []

dataSamples=SingleMuon+SingleElectron+Data23Sep2016

otherDataSamples=MuonEG+MET 

allDataSamples=dataSamples+otherDataSamples

# JSON
silverJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276384_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-279588_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-279931_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T_v2.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-282037_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-283685_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
#goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
goldenJson = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
run_range = (271036,284044)

jsonFile = goldenJson


from CMGTools.XZZ2l2nu.setup.Efficiencies import *
dataDir = "$CMSSW_BASE/src/CMGTools/XZZ2l2nu/data"

for comp in mcSamples+otherMcSamples:
    comp.isMC = True
    comp.isData = False
    comp.splitFactor = 250
    comp.puFileMC=dataDir+"/pileup_MC_80x_271036-276811_68075.root"
    comp.puFileData=dataDir+"/pileup_DATA_80x_271036-276811_68075.root"
    comp.eSFinput=dataDir+"/CutBasedID_LooseWP_76X_18Feb.txt_SF2D.root"
    comp.efficiency = eff2012
    #comp.triggers=triggers_1mu_noniso+triggers_1e_noniso
    comp.triggers= []
    comp.globalTag = "Summer15_25nsV6_MC"

for comp in allDataSamples:
    comp.splitFactor = 250
    comp.isMC = False
    comp.isData = True
    comp.json = jsonFile
    comp.run_range = run_range
    comp.globalTag = "Summer15_25nsV6_DATA"

