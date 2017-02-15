import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator
kreator = ComponentCreator()

# Photon+Jets
GJet_Pt_15To6000 = kreator.makeMCComponent("GJet_Pt_15To6000",
"/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_20M/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",154500, useAAA=False) # 8M evt

#GJet_Pt_20to40_DoubleEMEnriched = kreator.makeMCComponent("GJet_Pt_20to40_DoubleEMEnriched", 
#"/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
#"CMS", ".*root", 137751*0.001587, useAAA=False)   # 24M 

GJet_Pt_40toInf_DoubleEMEnriched = kreator.makeMCComponent("GJet_Pt_40toInf_DoubleEMEnriched", 
"/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 16792*0.0514, useAAA=False)  # 70M

GJet_Pt_20toInf_DoubleEMEnriched = kreator.makeMCComponent("GJet_Pt_20toInf_DoubleEMEnriched",
"/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCUETP8M1_13TeV_Pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",154500*0.02155, useAAA=False)  # 38M

GJet_Pt_EMEnriched=[
GJet_Pt_20toInf_DoubleEMEnriched,
GJet_Pt_40toInf_DoubleEMEnriched,
#GJet_Pt_20to40_DoubleEMEnriched,
]





# DY HT bins:
#https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns#DY_Z
DYJetsToLL_M50_HT100to200 = kreator.makeMCComponent("DYJetsToLL_M50_HT100to200", "/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM", "CMS", ".*root",139.4*1.23)
DYJetsToLL_M50_HT200to400 = kreator.makeMCComponent("DYJetsToLL_M50_HT200to400", "/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM", "CMS", ".*root",42.75*1.23)
DYJetsToLL_M50_HT400to600 = kreator.makeMCComponent("DYJetsToLL_M50_HT400to600", "/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM", "CMS", ".*root",5.497*1.23)
DYJetsToLL_M50_HT600toInf = kreator.makeMCComponent("DYJetsToLL_M50_HT600toInf", "/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM", "CMS", ".*root",2.21*1.23)

# cross-section:
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToGenXSecAnalyzer
# DY inclusive, NLO RunIISpring16MiniAODv2 
# 28M, x-sec recalculated using FEWZ using z_m50_nnlo_inclusive_NNPDF30_nlo_as_0118 QCD NNLO, QED NLO, including ISR, no FSR (because xsec reduction due to FSR is coming from the M50 mass cut)
DYJetsToLL_M50 = kreator.makeMCComponent("DYJetsToLL_M50", 
"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM", 
"CMS", ".*root", 1921.8*3, useAAA=False) 

DYJetsToLL_M50_Ext = kreator.makeMCComponent("DYJetsToLL_M50_Ext", 
"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM", 
"CMS", ".*root", 1921.8*3, useAAA=True)  #the "ext4" set with 129M evts 


 

DYJetsToLL_M50_MGMLM = kreator.makeMCComponent("DYJetsToLL_M50_MGMLM",
"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv1-PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/MINIAODSIM", 
"CMS", ".*root", 1921.8*3) # 9M


DYJetsToLL_M50_MGMLM_Ext1 = kreator.makeMCComponent("DYJetsToLL_M50_MGMLM_Ext1",
"/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM", 
"CMS", ".*root", 1921.8*3) # 50M

# LO* NLO kfactor 1.16261343013
# LO xsec calculated from miniAOD
# NLO/LO = njets NLO calculation / LO from miniAOD
#        = 1921.8*3/4959.0 = 1.16261343013
DY1JetsToLL_M50_MGMLM = kreator.makeMCComponent("DY1JetsToLL_M50_MGMLM",
"/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 1013.0*1.16261343013) 
DY2JetsToLL_M50_MGMLM = kreator.makeMCComponent("DY2JetsToLL_M50_MGMLM",
"/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 334.7*1.16261343013) 
DY3JetsToLL_M50_MGMLM = kreator.makeMCComponent("DY3JetsToLL_M50_MGMLM",
"/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 102.4*1.16261343013) 
DY4JetsToLL_M50_MGMLM = kreator.makeMCComponent("DY4JetsToLL_M50_MGMLM",
"/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 54.45*1.16261343013) 
DYBJetsToLL_M50_MGMLM = kreator.makeMCComponent("DYBJetsToLL_M50_MGMLM",
"/DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 88.2771*1.16261343013) 


# W+Jets
WJetsToLNu = kreator.makeMCComponent("WJetsToLNu",
"/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 3* 20508.9)

### DiBosons

# cross section from https://twiki.cern.ch/twiki/bin/view/CMS/SummaryTable1G25ns#Diboson

ZZTo2L2Nu = kreator.makeMCComponent("ZZTo2L2Nu", 
"/ZZTo2L2Nu_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 0.564) # 226 files, 8.8M evts, 208GB
ZZTo4L = kreator.makeMCComponent("ZZTo4L", 
"/ZZTo4L_13TeV_powheg_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 1.212,useAAA=False) #  6.7M evts
ZZTo2L2Q = kreator.makeMCComponent("ZZTo2L2Q", 
"/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 3.22) # 8.6M.
WWTo2L2Nu = kreator.makeMCComponent("WWTo2L2Nu", 
"/WWTo2L2Nu_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 12.178) #1.9M evts
WWToLNuQQ = kreator.makeMCComponent("WWToLNuQQ", 
"/WWToLNuQQ_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 49.997) # 1.9M evts
WWToLNuQQ_Ext1 = kreator.makeMCComponent("WWToLNuQQ_Ext1", 
"/WWToLNuQQ_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root", 49.997) # 1.9M evts
WZTo1L1Nu2Q = kreator.makeMCComponent("WZTo1L1Nu2Q", 
"/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM", 
"CMS", ".*root", 10.71) # 1.7M
WZTo2L2Q = kreator.makeMCComponent("WZTo2L2Q", 
"/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 5.595) # 25M
WZTo3LNu = kreator.makeMCComponent("WZTo3LNu", 
"/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS",  ".*root", 5.26,useAAA=False ) # 12.5M NLO up to 1 jet in ME #4.4297 from HZZ2l2nu note
### top

# previous: 104M: gen seg: https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/HIG-RunIIWinter15wmLHE-00518/0 
#TTTo2L2Nu = kreator.makeMCComponent("TTTo2L2Nu", 
#"/TTTo2L2Nu_13TeV-powheg/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/MINIAODSIM", 
#"CMS", ".*root", 831.76*((3*0.108)**2) ) # 104M
# noSC: 9.8M gen seg: https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/TOP-RunIISummer15wmLHEGS-00049/0 
TTTo2L2Nu_noSC = kreator.makeMCComponent("TTTo2L2Nu_noSC", 
"/TTTo2L2Nu_noSC_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 831.76*((3*0.108)**2) ) # 9.8M
# forTTH: 79M Gen frag https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_fragment/HIG-RunIISummer15wmLHEGS-00481/0
TTTo2L2Nu_forTTH = kreator.makeMCComponent("TTTo2L2Nu_forTTH", 
"/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 831.76*((3*0.108)**2) ) 

TTZToLLNuNu = kreator.makeMCComponent("TTZToLLNuNu", 
"/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root", 0.2529) # 2M
# old:
#TTWJetsToLNu = kreator.makeMCComponent("TTWJetsToLNu", 
#"/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM", 
#"CMS", ".*root", 0.2043) # 250k evt
TTWJetsToLNu_ext1 = kreator.makeMCComponent("TTWJetsToLNu_ext1", 
"/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/MINIAODSIM", 
"CMS", ".*root", 0.2043) # 2.2M evt
TTWJetsToLNu_ext2 = kreator.makeMCComponent("TTWJetsToLNu_ext2", 
"/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM", 
"CMS", ".*root", 0.2043) # 3.1M evt

### gamma+jets
### GJets Xsec: https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Gamma_jets
GJets_HT40to100 = kreator.makeMCComponent("GJets_HT40to100", 
"/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",20790,useAAA=False) #4M
GJets_HT100to200 = kreator.makeMCComponent("GJets_HT100to200", 
"/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",9238,useAAA=False) # 5 M
GJets_HT200to400 = kreator.makeMCComponent("GJets_HT200to400", 
"/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",2305,useAAA=False) # 10M 
GJets_HT400to600 = kreator.makeMCComponent("GJets_HT400to600", 
"/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",274.4,useAAA=False) # 2.4M
GJets_HT600toInf = kreator.makeMCComponent("GJets_HT600toInf", 
"/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",93.46,useAAA=False) # 2.45M

GJets_HT40to100_ext = kreator.makeMCComponent("GJets_HT40to100_ext", 
"/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",20790,useAAA=False) #4M
GJets_HT100to200_ext = kreator.makeMCComponent("GJets_HT100to200_ext", 
"/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",9238,useAAA=False) # 5 M
GJets_HT200to400_ext = kreator.makeMCComponent("GJets_HT200to400_ext", 
"/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",2305,useAAA=False) # 10M 
GJets_HT400to600_ext = kreator.makeMCComponent("GJets_HT400to600_ext", 
"/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM",
"CMS", ".*root",274.4,useAAA=False) # 2.4M
GJets_HT600toInf_ext = kreator.makeMCComponent("GJets_HT600toInf_ext", 
"/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",93.46,useAAA=False) # 2.45M

GJetsHT = [
GJets_HT40to100,
GJets_HT100to200,
GJets_HT200to400,
GJets_HT400to600,
GJets_HT600toInf,
GJets_HT40to100_ext,
GJets_HT100to200_ext,
GJets_HT200to400_ext,
GJets_HT400to600_ext,
GJets_HT600toInf_ext,
]

### ggZZ
ggZZTo2e2nu = kreator.makeMCComponent("ggZZTo2e2nu", 
"/GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 0.01898 ) #0.01898 
# xsec from McM : https://cms-pdmv.cern.ch/mcm/requests?prepid=HIG-RunIIWinter15pLHE-00083&page=0&shown=4063359
ggZZTo2mu2nu = kreator.makeMCComponent("ggZZTo2mu2nu", 
"/GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 0.01898 ) #0.01898


TTGJets = kreator.makeMCComponent("TTGJets",
"/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 3.697)
TTGJets_ext = kreator.makeMCComponent("TTGJets_ext",
"/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root", 3.697)

ttgjets = [TTGJets, TTGJets_ext]

ZNuNuGJetsGt130 = kreator.makeMCComponent("ZNuNuGJetsGt130", 
"/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 0.1903)
ZNuNuGJetsGt40Lt130 = kreator.makeMCComponent("ZNuNuGJetsGt40Lt130", 
"/ZNuNuGJets_MonoPhoton_PtG-40to130_TuneCUETP8M1_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",2.816)

ZNuNuGJets = [ZNuNuGJetsGt130,ZNuNuGJetsGt40Lt130]

WGToLNuG = kreator.makeMCComponent("WGToLNuG", 
"/WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 585.8)

WGJetsPt130 = kreator.makeMCComponent("WGJetsPt130", 
"/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 585.8)

wgjets = [WGToLNuG, WGJetsPt130]

# QCD HT bins (cross sections from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#QCD
QCD_HT50to100 = kreator.makeMCComponent("QCD_HT50to100", 
"/QCD_HT50to100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 2.465e8) # extracted from miniAOD files.
QCD_HT100to200 = kreator.makeMCComponent("QCD_HT100to200", 
"/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 2.799e7)
QCD_HT200to300 = kreator.makeMCComponent("QCD_HT200to300", 
"/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root",1.712e6)
QCD_HT200to300_ext = kreator.makeMCComponent("QCD_HT200to300_ext", 
"/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root",1.712e6)
QCD_HT300to500 = kreator.makeMCComponent("QCD_HT300to500", 
"/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",3.477e5)
QCD_HT300to500_ext = kreator.makeMCComponent("QCD_HT300to500_ext", 
"/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",3.477e5)
QCD_HT500to700 = kreator.makeMCComponent("QCD_HT500to700", 
"/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",3.21e4)
QCD_HT500to700_ext = kreator.makeMCComponent("QCD_HT500to700_ext", 
"/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM",
"CMS", ".*root",3.21e4)
QCD_HT700to1000 = kreator.makeMCComponent("QCD_HT700to1000", 
"/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",6831)
QCD_HT700to1000_ext = kreator.makeMCComponent("QCD_HT700to1000_ext", 
"/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",6831)
QCD_HT1000to1500 = kreator.makeMCComponent("QCD_HT1000to1500", 
"/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",1207)
QCD_HT1000to1500_ext = kreator.makeMCComponent("QCD_HT1000to1500_ext", 
"/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",1207)
QCD_HT1500to2000 = kreator.makeMCComponent("QCD_HT1500to2000", 
"/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",119.9)
QCD_HT1500to2000_ext = kreator.makeMCComponent("QCD_HT1500to2000_ext", 
"/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",119.9)
QCD_HT2000toInf = kreator.makeMCComponent("QCD_HT2000toInf", 
"/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",25.24)
QCD_HT2000toInf_ext = kreator.makeMCComponent("QCD_HT2000toInf_ext", 
"/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",25.24)


QCDHT = [
QCD_HT50to100,
QCD_HT100to200,
QCD_HT200to300,
QCD_HT200to300_ext,
QCD_HT300to500,
QCD_HT300to500_ext,
QCD_HT500to700,
QCD_HT500to700_ext,
QCD_HT700to1000,
QCD_HT700to1000_ext,
QCD_HT1000to1500,
QCD_HT1000to1500_ext,
QCD_HT1500to2000,
QCD_HT1500to2000_ext,
QCD_HT2000toInf,
QCD_HT2000toInf_ext
]


#QCD_Pt10to15     = kreator.makeMCComponent("QCD_Pt10to15"     , 
#"/QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"     , 
#"CMS" , ".*root" , 5887580000)
QCD_Pt15to30     = kreator.makeMCComponent("QCD_Pt15to30", 
"/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
 "CMS" , ".*root" , 1831718000)
QCD_Pt30to50     = kreator.makeMCComponent("QCD_Pt30to50"     , 
"/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS" , ".*root" , 139803000)
QCD_Pt50to80     = kreator.makeMCComponent("QCD_Pt50to80",    
"/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 19222500,useAAA=False)
QCD_Pt80to120    = kreator.makeMCComponent("QCD_Pt80to120", 
"/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",  2758420, useAAA=False)
QCD_Pt80to120_ext    = kreator.makeMCComponent("QCD_Pt80to120_ext", 
"/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM",
"CMS", ".*root",  2758420, useAAA=False)
QCD_Pt120to170   = kreator.makeMCComponent("QCD_Pt120to170"   , 
"/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",  
"CMS" , ".*root" , 469797,useAAA=False)
QCD_Pt120to170_ext   = kreator.makeMCComponent("QCD_Pt120to170_ext"   ,
"/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS" , ".*root" , 469797,useAAA=False)
QCD_Pt170to300   = kreator.makeMCComponent("QCD_Pt170to300"   , 
"/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS" , ".*root" , 117989, useAAA=False)
QCD_Pt170to300_ext   = kreator.makeMCComponent("QCD_Pt170to300_ext"   ,
"/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS" , ".*root" , 117989, useAAA=False)
QCD_Pt300to470   = kreator.makeMCComponent("QCD_Pt300to470",  
"/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 7820.25,useAAA=False)
QCD_Pt300to470_ext   = kreator.makeMCComponent("QCD_Pt300to470_ext",  
"/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root", 7820.25,useAAA=False)
QCD_Pt470to600   = kreator.makeMCComponent("QCD_Pt470to600",  
"/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 645.528 ,useAAA=False)
QCD_Pt470to600_bkup   = kreator.makeMCComponent("QCD_Pt470to600_bkup",  
"/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 645.528 ,useAAA=False)
QCD_Pt600to800   = kreator.makeMCComponent("QCD_Pt600to800"   , 
"/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS" , ".*root" , 187.109,useAAA=False)
QCD_Pt600to800_ext   = kreator.makeMCComponent("QCD_Pt600to800_ext"   , 
"/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS" , ".*root" , 187.109,useAAA=False)
QCD_Pt600to800_bkup   = kreator.makeMCComponent("QCD_Pt600to800_bkup"   ,
"/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS" , ".*root" , 187.109,useAAA=False)
QCD_Pt800to1000  = kreator.makeMCComponent("QCD_Pt800to1000"  , 
"/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS" , ".*root" , 32.3486 ,useAAA=False)
QCD_Pt800to1000_ext  = kreator.makeMCComponent("QCD_Pt800to1000_ext"  , 
"/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS" , ".*root" , 32.3486 ,useAAA=False)
QCD_Pt1000to1400 = kreator.makeMCComponent("QCD_Pt1000to1400" , 
"/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS" , ".*root" ,  10.4305,useAAA=False)
QCD_Pt1000to1400_ext = kreator.makeMCComponent("QCD_Pt1000to1400_ext" , 
"/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS" , ".*root" ,  10.4305,useAAA=False)
QCD_Pt1400to1800 = kreator.makeMCComponent("QCD_Pt1400to1800" ,
"/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS" , ".*root" , 0.84265,useAAA=False)
QCD_Pt1400to1800_ext = kreator.makeMCComponent("QCD_Pt1400to1800_ext" ,
"/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS" , ".*root" , 0.84265,useAAA=False)
QCD_Pt1800to2400 = kreator.makeMCComponent("QCD_Pt1800to2400" ,
"/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS" , ".*root" , 0.114943,useAAA=False)
QCD_Pt1800to2400_ext = kreator.makeMCComponent("QCD_Pt1800to2400_ext" ,
"/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS" , ".*root" , 0.114943,useAAA=False)
QCD_Pt2400to3200 = kreator.makeMCComponent("QCD_Pt2400to3200" ,
"/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS" , ".*root" , 0.00682981,useAAA=False)
QCD_Pt2400to3200_ext = kreator.makeMCComponent("QCD_Pt2400to3200_ext" ,
"/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS" , ".*root" , 0.00682981,useAAA=False)
QCD_Pt3200toInf  = kreator.makeMCComponent("QCD_Pt3200toInf", 
"/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM",
"CMS" , ".*root" , 0.000165445)

QCDPt = [
#QCD_Pt10to15,
QCD_Pt15to30,
QCD_Pt30to50,
QCD_Pt50to80,
QCD_Pt80to120,
QCD_Pt80to120_ext,
QCD_Pt120to170,
QCD_Pt120to170_ext,
QCD_Pt170to300,
QCD_Pt170to300_ext,
QCD_Pt300to470,
QCD_Pt300to470_ext,
QCD_Pt470to600,
QCD_Pt470to600_bkup,
QCD_Pt600to800,
QCD_Pt600to800_ext,
QCD_Pt600to800_bkup,
QCD_Pt800to1000,
QCD_Pt800to1000_ext,
QCD_Pt1000to1400,
QCD_Pt1000to1400_ext,
QCD_Pt1400to1800,
QCD_Pt1400to1800_ext,
QCD_Pt1800to2400,
QCD_Pt1800to2400_ext,
QCD_Pt2400to3200,
QCD_Pt2400to3200_ext,
QCD_Pt3200toInf
]

# qcd emenr
#QCD_Pt15to20_EMEnriched   = kreator.makeMCComponent("QCD_Pt15to20_EMEnriched"  ,"/QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM"  , "CMS", ".*root", 1273000000*0.0002)
QCD_Pt20to30_EMEnriched   = kreator.makeMCComponent("QCD_Pt20to30_EMEnriched"  ,
"/QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 557600000*0.0096)
QCD_Pt30to50_EMEnriched   = kreator.makeMCComponent("QCD_Pt30to50_EMEnriched"  ,
"/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 136000000*0.073)
QCD_Pt30to50_EMEnriched_ext   = kreator.makeMCComponent("QCD_Pt30to50_EMEnriched_ext"  ,
"/QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root", 136000000*0.073)
QCD_Pt50to80_EMEnriched   = kreator.makeMCComponent("QCD_Pt50to80_EMEnriched", 
"/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", # in production
"CMS", ".*root", 19800000*0.146)
QCD_Pt50to80_EMEnriched_ext   = kreator.makeMCComponent("QCD_Pt50to80_EMEnriched_ext", 
"/QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root", 19800000*0.146)
QCD_Pt80to120_EMEnriched  = kreator.makeMCComponent("QCD_Pt80to120_EMEnriched" ,
"/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 2800000*0.125)
QCD_Pt80to120_EMEnriched_ext  = kreator.makeMCComponent("QCD_Pt80to120_EMEnriched_ext" ,
"/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root", 2800000*0.125)
QCD_Pt120to170_EMEnriched = kreator.makeMCComponent("QCD_Pt120to170_EMEnriched",
"/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 477000*0.132)
QCD_Pt170to300_EMEnriched = kreator.makeMCComponent("QCD_Pt170to300_EMEnriched",
"/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 114000*0.165)
QCD_Pt300toInf_EMEnriched = kreator.makeMCComponent("QCD_Pt300toInf_EMEnriched",
"/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 9000*0.15)

QCDPtEMEnriched = [
#QCD_Pt15to20_EMEnriched,
QCD_Pt20to30_EMEnriched,
QCD_Pt30to50_EMEnriched,
QCD_Pt30to50_EMEnriched_ext,
QCD_Pt50to80_EMEnriched,
QCD_Pt50to80_EMEnriched_ext,
QCD_Pt80to120_EMEnriched,
QCD_Pt80to120_EMEnriched_ext,
QCD_Pt120to170_EMEnriched,
QCD_Pt170to300_EMEnriched,
QCD_Pt300toInf_EMEnriched
]


QCD_Pt15to20_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt15to20_MuEnrichedPt5",
"/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",1273190000*0.003, useAAA=False)
QCD_Pt20to30_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt20to30_MuEnrichedPt5",
"/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",558528000*0.0053, useAAA=False)
QCD_Pt30to50_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt30to50_MuEnrichedPt5",
"/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 139803000* 0.01182 , useAAA=False)
QCD_Pt50to80_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt50to80_MuEnrichedPt5",
"/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 19222500* 0.02276 , useAAA=False)
QCD_Pt80to120_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt80to120_MuEnrichedPt5",
"/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",2758420*0.03844 , useAAA=False)
QCD_Pt80to120_MuEnrichedPt5_ext1 = kreator.makeMCComponent("QCD_Pt80to120_MuEnrichedPt5_ext1",
"/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/MINIAODSIM",
"CMS", ".*root",2758420*0.03844, useAAA=False)
QCD_Pt120to170_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt120to170_MuEnrichedPt5",
"/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",469797*0.05362 , useAAA=False)
QCD_Pt120to170_MuEnrichedPt5_backup = kreator.makeMCComponent("QCD_Pt120to170_MuEnrichedPt5_backup",
"/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",469797*0.05362, useAAA=False)
QCD_Pt170to300_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt170to300_MuEnrichedPt5",
"/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",117989*0.07335, useAAA=False)
QCD_Pt170to300_MuEnrichedPt5_ext1 = kreator.makeMCComponent("QCD_Pt170to300_MuEnrichedPt5_ext1",
"/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",117989*0.07335, useAAA=False)
QCD_Pt170to300_MuEnrichedPt5_backup = kreator.makeMCComponent("QCD_Pt170to300_MuEnrichedPt5_backup",
"/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",117989*0.07335, useAAA=False)
QCD_Pt300to470_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt300to470_MuEnrichedPt5",
"/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 7820.25*0.10196, useAAA=False)
QCD_Pt300to470_MuEnrichedPt5_ext1 = kreator.makeMCComponent("QCD_Pt300to470_MuEnrichedPt5_ext1",
"/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",7820.25*0.10196, useAAA=False)
QCD_Pt300to470_MuEnrichedPt5_ext2 = kreator.makeMCComponent("QCD_Pt300to470_MuEnrichedPt5_ext2",
"/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM",
"CMS", ".*root",7820.25*0.10196, useAAA=False)
QCD_Pt470to600_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt470to600_MuEnrichedPt5",
"/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 645.528*0.12242, useAAA=False)
QCD_Pt470to600_MuEnrichedPt5_ext1 = kreator.makeMCComponent("QCD_Pt470to600_MuEnrichedPt5_ext1",
"/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",645.528*0.12242, useAAA=False)
QCD_Pt470to600_MuEnrichedPt5_ext2 = kreator.makeMCComponent("QCD_Pt470to600_MuEnrichedPt5_ext2",
"/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM",
"CMS", ".*root",645.528*0.12242, useAAA=False)
QCD_Pt600to800_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt600to800_MuEnrichedPt5",
"/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",187.109*0.13412, useAAA=False)
QCD_Pt600to800_MuEnrichedPt5_ext1 = kreator.makeMCComponent("QCD_Pt600to800_MuEnrichedPt5_ext1",
"/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",187.109*0.13412, useAAA=False)
QCD_Pt600to800_MuEnrichedPt5_backup = kreator.makeMCComponent("QCD_Pt600to800_MuEnrichedPt5_backup",
"/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_backup_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",187.109*0.13412, useAAA=False)
QCD_Pt800to1000_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt800to1000_MuEnrichedPt5",
"/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",32.3486*0.14552, useAAA=False)
QCD_Pt800to1000_MuEnrichedPt5_ext1 = kreator.makeMCComponent("QCD_Pt800to1000_MuEnrichedPt5_ext1",
"/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root", 32.3486*0.14552, useAAA=False)
QCD_Pt800to1000_MuEnrichedPt5_ext2 = kreator.makeMCComponent("QCD_Pt800to1000_MuEnrichedPt5_ext2",
"/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM",
"CMS", ".*root", 32.3486*0.14552, useAAA=False)
QCD_Pt1000toInf_MuEnrichedPt5 = kreator.makeMCComponent("QCD_Pt1000toInf_MuEnrichedPt5",
"/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",10.4305*0.15544, useAAA=False)
QCD_Pt1000toInf_MuEnrichedPt5_ext1 = kreator.makeMCComponent("QCD_Pt1000toInf_MuEnrichedPt5_ext1",
"/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v3/MINIAODSIM",
"CMS", ".*root",10.4305*0.15544, useAAA=False)
QCD_Pt20toInf_MuEnrichedPt15 = kreator.makeMCComponent("QCD_Pt20toInf_MuEnrichedPt15",
"/QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",720648000*0.00042, useAAA=False)

QCDPtMuEnriched=[
QCD_Pt15to20_MuEnrichedPt5,
QCD_Pt20to30_MuEnrichedPt5,
QCD_Pt30to50_MuEnrichedPt5,
QCD_Pt50to80_MuEnrichedPt5,
QCD_Pt80to120_MuEnrichedPt5,
QCD_Pt80to120_MuEnrichedPt5_ext1,
QCD_Pt120to170_MuEnrichedPt5,
QCD_Pt120to170_MuEnrichedPt5_backup,
QCD_Pt170to300_MuEnrichedPt5,
QCD_Pt170to300_MuEnrichedPt5_ext1,
QCD_Pt170to300_MuEnrichedPt5_backup,
QCD_Pt300to470_MuEnrichedPt5,
QCD_Pt300to470_MuEnrichedPt5_ext1,
QCD_Pt300to470_MuEnrichedPt5_ext2,
QCD_Pt470to600_MuEnrichedPt5,
QCD_Pt470to600_MuEnrichedPt5_ext1,
QCD_Pt470to600_MuEnrichedPt5_ext2,
QCD_Pt600to800_MuEnrichedPt5,
QCD_Pt600to800_MuEnrichedPt5_ext1,
QCD_Pt600to800_MuEnrichedPt5_backup,
QCD_Pt800to1000_MuEnrichedPt5,
QCD_Pt800to1000_MuEnrichedPt5_ext1,
QCD_Pt800to1000_MuEnrichedPt5_ext2,
QCD_Pt1000toInf_MuEnrichedPt5,
QCD_Pt1000toInf_MuEnrichedPt5_ext1,
QCD_Pt20toInf_MuEnrichedPt15,
]

#TToLeptons_tch_powheg = kreator.makeMCComponent("TToLeptons_tch_powheg", 
#"/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM", 
#"CMS", ".*root", (136.02)*0.108*3)
#TBarToLeptons_tch_powheg = kreator.makeMCComponent("TBarToLeptons_tch_powheg", 
#"/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/MINIAODSIM", 
#"CMS", ".*root", 80.95*0.108*3)

T_tch_powheg = kreator.makeMCComponent("T_tch_powheg", 
"/ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 136.02)
TBar_tch_powheg = kreator.makeMCComponent("TBar_tch_powheg",
"/ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root", 80.95)

TBar_tWch = kreator.makeMCComponent("TBar_tWch", 
"/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",35.6)
T_tWch= kreator.makeMCComponent("T_tWch", 
"/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",35.6)

TGJets = kreator.makeMCComponent("TGJets",
"/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root", 2.967)
TGJets_ext = kreator.makeMCComponent("TGJets_ext", 
"/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root", 2.967)

SingleTop = [
T_tch_powheg,
TBar_tch_powheg,
T_tWch,
TBar_tWch,
TGJets,
TGJets_ext
]

### Zinv
ZJetsToNuNu_HT100to200 = kreator.makeMCComponent("ZJetsToNuNu_HT100to200", 
"/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",280.47*1.23)
ZJetsToNuNu_HT100to200_ext = kreator.makeMCComponent("ZJetsToNuNu_HT100to200_ext", 
"/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",280.47*1.23)
ZJetsToNuNu_HT200to400 = kreator.makeMCComponent("ZJetsToNuNu_HT200to400", 
"/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",78.36*1.23)
ZJetsToNuNu_HT200to400_ext = kreator.makeMCComponent("ZJetsToNuNu_HT200to400_ext", 
"/ZJetsToNuNu_HT-200To400_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",78.36*1.23)
ZJetsToNuNu_HT400to600 = kreator.makeMCComponent("ZJetsToNuNu_HT400to600", 
"/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",10.94*1.23)
ZJetsToNuNu_HT600to800 = kreator.makeMCComponent("ZJetsToNuNu_HT600to800", 
"/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",3.221*1.23)
ZJetsToNuNu_HT800to1200 = kreator.makeMCComponent("ZJetsToNuNu_HT800t1200", 
"/ZJetsToNuNu_HT-800To1200_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",1.474*1.23)
ZJetsToNuNu_HT1200to2500 = kreator.makeMCComponent("ZJetsToNuNu_HT1200to2500", 
"/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",0.3586*1.23 )
ZJetsToNuNu_HT1200to2500_ext = kreator.makeMCComponent("ZJetsToNuNu_HT1200to2500_ext", 
"/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM",
"CMS", ".*root",0.3586*1.23 )
ZJetsToNuNu_HT2500toInf = kreator.makeMCComponent("ZJetsToNuNu_HT2500toInf", 
"/ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM",
"CMS", ".*root",0.008203*1.23)

#280.47*1.23+78.36*1.23+10.94*1.23+1.474*1.23+0.3586*1.23+0.008203*1.23


ZJetsToNuNuHT = [
ZJetsToNuNu_HT100to200,
ZJetsToNuNu_HT100to200_ext,
ZJetsToNuNu_HT200to400,
ZJetsToNuNu_HT200to400_ext,
ZJetsToNuNu_HT400to600,
#ZJetsToNuNu_HT400to600_ext,
#ZJetsToNuNu_HT600toInf,
ZJetsToNuNu_HT600to800,
#ZJetsToNuNu_HT600to800_ext,
ZJetsToNuNu_HT800to1200,
#ZJetsToNuNu_HT800to1200ext,
ZJetsToNuNu_HT1200to2500,
ZJetsToNuNu_HT1200to2500_ext,
ZJetsToNuNu_HT2500toInf,
#ZJetsToNuNu_HT2500toInf_ext,
]



### W+jets
# xsec from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#W_jets
WJetsToLNu_HT100to200 = kreator.makeMCComponent("WJetsToLNu_HT100to200", 
"/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root",1345*1.21,useAAA=False)
WJetsToLNu_HT100to200_ext1 = kreator.makeMCComponent("WJetsToLNu_HT100to200_ext1", 
"/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root",1345*1.21,useAAA=False)
WJetsToLNu_HT100to200_ext2 = kreator.makeMCComponent("WJetsToLNu_HT100to200_ext2", 
"/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM", 
"CMS", ".*root",1345*1.21,useAAA=False)
WJetsToLNu_HT200to400 = kreator.makeMCComponent("WJetsToLNu_HT200to400", 
"/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root",359.7*1.21,useAAA=False)
WJetsToLNu_HT200to400_ext1 = kreator.makeMCComponent("WJetsToLNu_HT200to400_ext1", 
"/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root",359.7*1.21,useAAA=False)
WJetsToLNu_HT200to400_ext2 = kreator.makeMCComponent("WJetsToLNu_HT200to400_ext2", 
"/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM", 
"CMS", ".*root",359.7*1.21,useAAA=False)
WJetsToLNu_HT400to600 = kreator.makeMCComponent("WJetsToLNu_HT400to600", 
"/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root",48.91*1.21,useAAA=False)
WJetsToLNu_HT400to600_ext1 = kreator.makeMCComponent("WJetsToLNu_HT400to600_ext1", 
"/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root",48.91*1.21,useAAA=False)
WJetsToLNu_HT600to800 = kreator.makeMCComponent("WJetsToLNu_HT600to800", 
"/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root",12.05*1.21,useAAA=False)
WJetsToLNu_HT600to800_ext1 = kreator.makeMCComponent("WJetsToLNu_HT600to800_ext1", 
"/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root",12.05*1.21,useAAA=False)
WJetsToLNu_HT800to1200 = kreator.makeMCComponent("WJetsToLNu_HT800to1200", 
"/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root",5.501*1.21,useAAA=False)
WJetsToLNu_HT800to1200_ext1 = kreator.makeMCComponent("WJetsToLNu_HT800to1200_ext1", 
"/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root",5.501*1.21,useAAA=False)
WJetsToLNu_HT1200to2500 = kreator.makeMCComponent("WJetsToLNu_HT1200to2500", 
"/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root",1.329*1.21,useAAA=False)
WJetsToLNu_HT1200to2500_ext1 = kreator.makeMCComponent("WJetsToLNu_HT1200to2500_ext1", 
"/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root",1.329*1.21,useAAA=False)
WJetsToLNu_HT2500toInf = kreator.makeMCComponent("WJetsToLNu_HT2500toInf", 
"/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM", 
"CMS", ".*root",0.03216*1.21,useAAA=False)
WJetsToLNu_HT2500toInf_ext1 = kreator.makeMCComponent("WJetsToLNu_HT2500toInf_ext1", 
"/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM", 
"CMS", ".*root",0.03216*1.21,useAAA=False)

# 1345*1.21+359.7*1.21+48.91*1.21+12.05*1.21+5.501*1.21+1.329*1.21+0.03216*1.21


WJetsToLNuHT = [
WJetsToLNu_HT100to200,
WJetsToLNu_HT100to200_ext1,
WJetsToLNu_HT100to200_ext2,
WJetsToLNu_HT200to400,
WJetsToLNu_HT200to400_ext1,
WJetsToLNu_HT200to400_ext2,
WJetsToLNu_HT400to600,
WJetsToLNu_HT400to600_ext1,
WJetsToLNu_HT600to800,
WJetsToLNu_HT600to800_ext1,
WJetsToLNu_HT800to1200,
WJetsToLNu_HT800to1200_ext1,
WJetsToLNu_HT1200to2500,
WJetsToLNu_HT1200to2500_ext1,
WJetsToLNu_HT2500toInf,
WJetsToLNu_HT2500toInf_ext1
]

