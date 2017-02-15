##########################################################
##      configuration for XZZ2l2nu 
##########################################################

import CMGTools.XZZ2l2nu.fwlite.Config as cfg
from CMGTools.XZZ2l2nu.fwlite.Config import printComps
from CMGTools.XZZ2l2nu.RootTools import *
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption


#Load all common analyzers
from CMGTools.XZZ2l2nu.analyzers.coreXZZ_cff import *

#-------- SAMPLES AND TRIGGERS -----------
from CMGTools.XZZ2l2nu.samples.loadSamples80xSummer16 import *
selectedComponents = mcSamples+dataSamples

triggerFlagsAna.triggerBits ={
    "ISOMU":triggers_1mu_iso,
    "MU":triggers_1mu_noniso,
    "MUv2":triggers_1mu_noniso_v2,
    "MU50":triggers_1mu_noniso_M50,
    "TkMU50":triggers_1mu_noniso_tkM50,
    "ISOELE":triggers_1e,
    "ELE":triggers_1e_noniso,
    "ELEv2":triggers_1e_noniso_v2,
    "ELE115":triggers_1e_noniso_E115,
    "MUMU": triggers_mumu,
    "MUMUNOISO":triggers_mumu_noniso,
    "ELEL": triggers_ee,
    "HT800":triggers_HT800,
    "HT900":triggers_HT900,
    "JJ":triggers_dijet_fat,
    "MET90":triggers_met90_mht90+triggers_metNoMu90_mhtNoMu90,
    "MET120":triggers_metNoMu120_mhtNoMu120,
    "ALLPHOTON": triggers_all_photons,
    "PHOTONIDISO": triggers_photon_idiso,
    "PHOTONIDISOMETEB": triggers_photon_idisometeb,
    "PHOTONIDISOVBFEB": triggers_photon_idisovbfeb,
    "HALOCLEAN": triggers_halo_clean,
    "PHOTONIDISO22": triggers_photon22_idiso,
    "PHOTONIDISO30": triggers_photon30_idiso,
    "PHOTONIDISO36": triggers_photon36_idiso,
    "PHOTONIDISO50": triggers_photon50_idiso,
    "PHOTONIDISO75": triggers_photon75_idiso,
    "PHOTONIDISO90": triggers_photon90_idiso,
    "PHOTONIDISO120": triggers_photon120_idiso,
    "PHOTONIDISO165": triggers_photon165_idiso,
#    "PHOTONIDISOMETEB22": triggers_photon22_idisometeb,
#    "PHOTONIDISOMETEB30": triggers_photon30_idisometeb,
#    "PHOTONIDISOMETEB36": triggers_photon36_idisometeb,
#    "PHOTONIDISOMETEB50": triggers_photon50_idisometeb,
#    "PHOTONIDISOMETEB75": triggers_photon75_idisometeb,
#    "PHOTONIDISOMETEB90": triggers_photon90_idisometeb,
#    "PHOTONIDISOMETEB120": triggers_photon120_idisometeb,
#    "PHOTONIDISOVBFEB22": triggers_photon22_idisovbfeb,
#    "PHOTONIDISOVBFEB30": triggers_photon30_idisovbfeb,
#    "PHOTONIDISOVBFEB36": triggers_photon36_idisovbfeb,
#    "PHOTONIDISOVBFEB50": triggers_photon50_idisovbfeb,
#    "PHOTONIDISOVBFEB75": triggers_photon75_idisovbfeb,
#    "PHOTONIDISOVBFEB90": triggers_photon90_idisovbfeb,
#    "PHOTONIDISOVBFEB120": triggers_photon120_idisovbfeb,
}

#-------- Analyzer
from CMGTools.XZZ2l2nu.analyzers.treeXZZ_cff import *

multiStateAna.processTypes = ['PhotonJets']
multiStateAna.selectPhotonJets = (lambda x: x.leg1.pt()>20.0 and x.leg2.pt()>-0.0)
vvSkimmer.required = ['PhotonJets']

vvTreeProducer.collections = {
	 "jets"       : NTupleCollection("jet",JetType,100, help="all jets in miniaod"),
         "selectedPhotons"       : NTupleCollection("photon",photonType,100, help="selected photons in miniaod"),
         "selectedLeptons" : NTupleCollection("lep",leptonType,10, help="selected leptons"),
         "PhotonJets"     : NTupleCollection("gjet",PhotonJetType ,100, help="photon and MET"),
         "badMuons" : NTupleCollection("badmuon",leptonType,100, help="bad muons"),
     }




#-------- SEQUENCE
coreSequence = [
    skimAnalyzer,
    genAna,
    jsonAna,
    triggerAna,
    pileUpAna,
    vertexAna,
    lepAna,
    photonAna, 
    jetAna,
    metAna,
    multiStateAna,
    eventFlagsAna,
    triggerFlagsAna,
]
multtrg.photonjet=True
multtrg.HLTlist=[
    "HLT_Photon22_R9Id90_HE10_IsoM_v",
    "HLT_Photon30_R9Id90_HE10_IsoM_v",
    "HLT_Photon36_R9Id90_HE10_IsoM_v",
    "HLT_Photon50_R9Id90_HE10_IsoM_v",
    "HLT_Photon75_R9Id90_HE10_IsoM_v",
    "HLT_Photon90_R9Id90_HE10_IsoM_v",
    "HLT_Photon120_R9Id90_HE10_IsoM_v",
    "HLT_Photon165_R9Id90_HE10_IsoM_v",
    ]    
'''
README

The branch gjet_l1_trigerob_HLTbit in the output ntuple can be used to determin which HLT object could be matched with the photon. For example, "HLT_Photon50_R9Id90_HE10_IsoM" is the 3rd element in the list multtrg.HLTlist, then by doing "(gjet_l1_trigerob_HLTbit>>3)&1", one can tell if HLT_Photon50_R9Id90_HE10_IsoM is matched with the photon (1 for yes 0 for no)

'''
#for ihlt in multtrg.HLTlist:
#    vvTreeProducer.globalVariables.append(
#         NTupleVariable("PreScale"+ihlt[3:-2],  lambda ev: getattr(ev,ihlt+"PS"), int, help="Photon HLT prescale")
#)
vvTreeProducer.globalVariables.append(NTupleVariable("PreScale22",  lambda ev: getattr(ev,"HLT_Photon22_R9Id90_HE10_IsoM_vPS"), int, help="Photon HLT prescale"))
vvTreeProducer.globalVariables.append(NTupleVariable("PreScale30",  lambda ev: getattr(ev,"HLT_Photon30_R9Id90_HE10_IsoM_vPS"), int, help="Photon HLT prescale"))
vvTreeProducer.globalVariables.append(NTupleVariable("PreScale36",  lambda ev: getattr(ev,"HLT_Photon36_R9Id90_HE10_IsoM_vPS"), int, help="Photon HLT prescale"))
vvTreeProducer.globalVariables.append(NTupleVariable("PreScale50",  lambda ev: getattr(ev,"HLT_Photon50_R9Id90_HE10_IsoM_vPS"), int, help="Photon HLT prescale"))
vvTreeProducer.globalVariables.append(NTupleVariable("PreScale75",  lambda ev: getattr(ev,"HLT_Photon75_R9Id90_HE10_IsoM_vPS"), int, help="Photon HLT prescale"))
vvTreeProducer.globalVariables.append(NTupleVariable("PreScale90",  lambda ev: getattr(ev,"HLT_Photon90_R9Id90_HE10_IsoM_vPS"), int, help="Photon HLT prescale"))
vvTreeProducer.globalVariables.append(NTupleVariable("PreScale120",  lambda ev: getattr(ev,"HLT_Photon120_R9Id90_HE10_IsoM_vPS"), int, help="Photon HLT prescale"))
vvTreeProducer.globalVariables.append(NTupleVariable("PreScale165",  lambda ev: getattr(ev,"HLT_Photon165_R9Id90_HE10_IsoM_vPS"), int, help="Photon HLT prescale"))

#sequence = cfg.Sequence(coreSequence)
#sequence = cfg.Sequence(coreSequence+[vvSkimmer,vvTreeProducer])
sequence = cfg.Sequence(coreSequence+[vvSkimmer,multtrg,vvTreeProducer])
#sequence = cfg.Sequence(coreSequence+[vvSkimmer,fullTreeProducer])
 

#-------- HOW TO RUN
test = 1
if test==1:
    # test a single component, using a single thread.
    #selectedComponents = dataSamples
    #selectedComponents = mcSamples
    #selectedComponents = SinglePhoton
    #selectedComponents = SinglePhoton23Sep2016+[SinglePhoton_Run2016H_PromptReco_v1,SinglePhoton_Run2016H_PromptReco_v2]
    #selectedComponents = SinglePhoton23Sep2016+[SinglePhoton_Run2016H_PromptReco_v2]
    #selectedComponents = [SinglePhoton_Run2016C_23Sep2016]
    #selectedComponents = [SinglePhoton_Run2016D_23Sep2016]
    #selectedComponents = [SinglePhoton_Run2016B_PromptReco,
    #                      SinglePhoton_Run2016B_PromptReco_v2,
    #                      SinglePhoton_Run2016C_PromptReco_v2,
    #                      SinglePhoton_Run2016D_PromptReco_v2,
    #                      SinglePhoton_Run2016E_PromptReco_v2,
    #                      SinglePhoton_Run2016F_PromptReco_v1,
    #                      SinglePhoton_Run2016G_PromptReco_v1,]
 
    #selectedComponents = [SinglePhoton_Run2016D_PromptReco_v2]
    #selectedComponents = [SinglePhoton_Run2016B_23Sep2016, SinglePhoton_Run2016B_23Sep2016_v2,]
    #selectedComponents = [SinglePhoton_Run2016C_23Sep2016,]
    #selectedComponents = [SinglePhoton_Run2016D_23Sep2016,]
    #selectedComponents = [SinglePhoton_Run2016E_23Sep2016,]
    #selectedComponents = [SinglePhoton_Run2016F_23Sep2016,]
    #selectedComponents = [SinglePhoton_Run2016G_23Sep2016,]
    #selectedComponents = [SinglePhoton_Run2016H_PromptReco_v1, SinglePhoton_Run2016H_PromptReco_v2, SinglePhoton_Run2016H_PromptReco_v3]
    #selectedComponents = [SinglePhoton_Run2016H_PromptReco_v2]
    #selectedComponents = [SinglePhoton_Run2016H_PromptReco_v2, SinglePhoton_Run2016H_PromptReco_v3]
    #selectedComponents = [GJet_Pt_20toInf_DoubleEMEnriched]
 #   selectedComponents = [GJet_Pt_40toInf_DoubleEMEnriched]

    #selectedComponents = [SinglePhoton_Run2016F_PromptReco_v1]
    #selectedComponents = [GJet_Pt_20toInf_DoubleEMEnriched]
    #selectedComponents = [GJet_Pt_20to40_DoubleEMEnriched, GJet_Pt_40toInf_DoubleEMEnriched]
    #selectedComponents = [SingleMuon_Run2015D_Promptv4,SingleElectron_Run2015D_Promptv4]
    #selectedComponents = [SingleMuon_Run2015C_25ns_16Dec]
    #selectedComponents = [SingleMuon_Run2016B_PromptReco_v2] 
    #selectedComponents = SingleMuon+SingleElectron
    #selectedComponents = [SingleMuon_Run2016B_PromptReco_v2,SingleElectron_Run2016B_PromptReco_v2] 
    #selectedComponents = [SingleMuon_Run2016D_PromptReco_v2,SingleElectron_Run2016D_PromptReco_v2] 
    #selectedComponents = [MuonEG_Run2015D_16Dec] #MuEG
    #selectedComponents = [RSGravToZZToZZinv_narrow_800]
    #selectedComponents = [DYJetsToLL_M50]
    #selectedComponents = [DYJetsToLL_M50_MGMLM_Ext1]
    #selectedComponents = [BulkGravToZZToZlepZinv_narrow_600] 
    #selectedComponents = signalSamples
    #selectedComponents = [TTTo2L2Nu]
    #selectedComponents = [BulkGravToZZ_narrow_800]
    #selectedComponents = [BulkGravToZZToZlepZhad_narrow_800]
    #selectedComponents = WJetsToLNuHT + ZJetsToNuNuHT + SingleTop + ZNuNuGJets + ttgjets + wgjets + GJetsHT + QCDPtEMEnriched
    selectedComponents = [DYJetsToLL_M50_Ext]
    #selectedComponents = GJet_Pt_EMEnriched + GJetsHT
    #selectedComponents = QCDPtEMEnriched + QCDPt + QCDHT
    #selectedComponents = [DYJetsToLL_M50_MGMLM_Ext1]
    for c in selectedComponents:
#        c.files = c.files[:2]
        #c.splitFactor = (len(c.files)/10 if len(.files)>10 else 1)
        c.splitFactor = len(c.files)
#        c.splitFactor = 1
        #c.triggers=triggers_1mu_noniso
        #c.triggers=triggers_1e_noniso

## output histogram
outputService=[]
from PhysicsTools.HeppyCore.framework.services.tfile import TFileService
output_service = cfg.Service(
    TFileService,
    'outputfile',
    name="outputfile",
    fname='vvTreeProducer/tree.root',
    option='recreate'
    )
outputService.append(output_service)

from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
preprocessor = CmsswPreprocessor("pogRecipes.py")

from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
event_class = Events
config = cfg.Config( components = selectedComponents,
                     sequence = sequence,
                     services = [],
                     preprocessor=preprocessor, #this would run cmsRun before running Heppy
                     events_class = event_class)

# and the following runs the process directly if running as with python filename.py  
if __name__ == '__main__':
    from PhysicsTools.HeppyCore.framework.looper import Looper
    looper = Looper( 'Loop', config, nPrint = 5,nEvents=300)
    looper.loop()
    looper.write()





