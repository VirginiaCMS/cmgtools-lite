// Hengne Li, 2016 @ CERN, initial version.

#include "metcorr.h"


int main(int argc, char** argv) {

  if( argc<5 ) {
     std::cout << argv[0] << ":  " << std::endl ;
     std::cout << " Functionality: met kin-fit framework, actually a full analysis framework  ... "  << std::endl;
     std::cout << "                 "  << std::endl;
     std::cout << " usage: " << argv[0] << " config.file inputfile.root outputfile.root Nevts SumWeights " << std::endl ;
     exit(1) ;
  }

  // config file name
  _file_config_name = std::string((const char*)argv[1]);

  // input file name
  _file_in_name = std::string((const char*)argv[2]);

  // output file name
  _file_out_name = std::string((const char*)argv[3]);

  // sum events
  _SumEvents = atof(argv[4]);

  // sum weights
  _SumWeights = atof(argv[5]);


  // input file
  _file_in = TFile::Open(_file_in_name.c_str());

  // output file 
  _file_out = TFile::Open(_file_out_name.c_str(), "recreate");



  // read config file
  readConfigFile();

  // prepare the trees
  if (!prepareTrees()) return 1;

  // prepare inputs for pu weights
  if (_addPUWeights && !_isData) preparePUWeights();

  // prepare inputs for muon re-calib
  if (_doMuonPtRecalib) prepareMuonPtRecalib();


  // prepare inputs for addDyZPtWeight
  if (_addDyZPtWeight && !_isData) prepareDyZPtWeight();

  // prepare inputs for simple met recoil tune.
  if (_doRecoil && !_isData) prepareRecoil();


  // loop
  int n_pass = 0;

  // get n selected entries to loop 
  Long64_t n_selected_entries = _selected_entries->GetN();

  for (Int_t i_slist = _n_start; i_slist < n_selected_entries; i_slist++) {

    // get _tree_in's entry number from selected list
    Int_t i = _selected_entries->GetEntry(i_slist);

    // get the selected entry from tree_in
    _tree_in->GetEntry(i);

    // beak if runs more events than  you want to test
    n_pass++;
    if (_n_test>0 && n_pass>_n_test) break;

    //
    if (_debug || i%_n_interval == 0) {
      std::cout << "#############################################" << std::endl;
      std::cout << "##  Entry " << i << ", Run " << _run << ", Event " << _evt << std::endl;
      std::cout << "#############################################" << std::endl;
    }


    // add pu weights
    if (_addPUWeights && !_isData) addPUWeights();

    // do muon re-calib
    if (_doMuonPtRecalib) doMuonPtRecalib(); 

    //  addDyZPtWeight
    if (_addDyZPtWeight && !_isData) addDyZPtWeight();

    // simple met recoil tune.
    if (_doRecoil && !_isData) doRecoil();



    // fill output tree
    _tree_out->Fill(); 
  }  


  // store output tree
  _file_out->cd();
  _tree_out->Write();
  _file_out->Close();

  return 0;

}





//======================================================
// ╦═╗╔═╗╔═╗╔╦╗  ╔═╗╔═╗╔╗╔╔═╗╦╔═╗  ╔═╗╦╦  ╔═╗
// ╠╦╝║╣ ╠═╣ ║║  ║  ║ ║║║║╠╣ ║║ ╦  ╠╣ ║║  ║╣ 
// ╩╚═╚═╝╩ ╩═╩╝  ╚═╝╚═╝╝╚╝╚  ╩╚═╝  ╚  ╩╩═╝╚═╝
//======================================================

void readConfigFile() 
{

  // init paramter file
  PParameterReader parm(_file_config_name.c_str());

  // print it
  parm.Print();

  //===============================================
  // NOTE: 
  //  - See notes above in parameter definitions
  //     for the meanings and usage of them.
  //  - Please always use the same variable name
  //     in the c++ codes and in the config files,
  //     with the c++ one has "_" in front, 
  //      e.g.:   
  //       bool _AAA = parm.GetBool("AAA", kFALSE);
  //     This rule can minimize mistakes.
  //==============================================

  //======================
  // general parameters
  //======================
  _debug = parm.GetBool("debug", kFALSE);
  _n_start = parm.GetInt("n_start", 0);
  _n_test = parm.GetInt("n_test", -1);
  _n_interval = parm.GetInt("n_interval", 3000);
  _selection = parm.GetString("selection", "(1)");
  _storeOldBranches = parm.GetBool("storeOldBranches", kFALSE);


  //==========================
  // use slimmed tree or not
  //=========================
  _useLightTree = parm.GetBool("useLightTree", kTRUE);


  //=========================
  // add PU weights
  //=========================
  _addPUWeights = parm.GetBool("addPUWeights", kTRUE);
  _PUTags = parm.GetVString("PUTags");
  _PUInputDir = parm.GetString("PUInputDir", "data/pileup");
  _PUInputFileNames = parm.GetVString("PUInputFileNames");
  _PUWeightHistName = parm.GetString("PUWeightHistName", "puweight_dtmc");


  //==========================
  // muon pT recalibration
  //==========================
  _doMuonPtRecalib = parm.GetBool("doMuonPtRecalib", kFALSE);


  //========================
  // Add DYJet gen reweight 
  //========================
  _addDyZPtWeight = parm.GetBool("addDyZPtWeight", kTRUE);
  _addDyZPtWeightUseFunction = parm.GetBool("addDyZPtWeightUseFunction", kTRUE);
  _addDyZPtWeightLOUseFunction = parm.GetBool("addDyZPtWeightLOUseFunction", kTRUE);
  _addDyNewGenWeight = parm.GetBool("addDyNewGenWeight", kTRUE);;


  //==============================================
  // Do simple version of the MET recoil tuning
  //==============================================
  _doRecoil = parm.GetBool("doRecoil", kFALSE);
  _doRecoilUseSmooth = parm.GetBool("doRecoilUseSmooth", kTRUE);
  _doRecoilUseSmoothGraph = parm.GetBool("doRecoilUseSmoothGraph", kTRUE);



}




// prepare the trees
bool  prepareTrees() 
{

  // 1.) prepare in put tree:

  // input tree
  _tree_in = (TTree*)_file_in->Get("tree");

  // selection 
  _tree_in->Draw(">>_selected_entries", _selection.c_str(), "entrylist");
  _selected_entries = (TEntryList*)gDirectory->Get("_selected_entries");

  // isData  
  _tree_in->SetBranchAddress("isData",&_isData);

  // check if tree has events
  if (_tree_in->GetEntries()<=0) {
    std::cout << "prepareTrees():: input tree has no event pass selection. Quit." <<std::endl;
    return false;
  }

  // get isData info
  _tree_in->GetEntry(0);

  // set common branches
  _tree_in->SetBranchAddress("run", &_run);
  _tree_in->SetBranchAddress("lumi", &_lumi);
  _tree_in->SetBranchAddress("evt", &_evt);

  _tree_in->SetBranchAddress("llnunu_mt", &_llnunu_mt);
  _tree_in->SetBranchAddress("llnunu_l1_mass",&_llnunu_l1_mass);
  _tree_in->SetBranchAddress("llnunu_l1_mt", &_llnunu_l1_mt);
  _tree_in->SetBranchAddress("llnunu_l1_pt", &_llnunu_l1_pt);
  _tree_in->SetBranchAddress("llnunu_l1_phi", &_llnunu_l1_phi);
  _tree_in->SetBranchAddress("llnunu_l1_eta", &_llnunu_l1_eta);
  _tree_in->SetBranchAddress("llnunu_l1_deltaPhi", &_llnunu_l1_deltaPhi);
  _tree_in->SetBranchAddress("llnunu_l1_deltaR", &_llnunu_l1_deltaR);
  _tree_in->SetBranchAddress("llnunu_l1_rapidity", &_llnunu_l1_rapidity);

  _tree_in->SetBranchAddress("llnunu_l2_pt", &_llnunu_l2_pt);
  _tree_in->SetBranchAddress("llnunu_l2_phi", &_llnunu_l2_phi);

  _tree_in->SetBranchAddress("llnunu_l1_l1_pt", &_llnunu_l1_l1_pt);
  _tree_in->SetBranchAddress("llnunu_l1_l1_eta", &_llnunu_l1_l1_eta);
  _tree_in->SetBranchAddress("llnunu_l1_l1_phi", &_llnunu_l1_l1_phi);
  _tree_in->SetBranchAddress("llnunu_l1_l1_rapidity", &_llnunu_l1_l1_rapidity);
  _tree_in->SetBranchAddress("llnunu_l1_l1_mass", &_llnunu_l1_l1_mass);
  _tree_in->SetBranchAddress("llnunu_l1_l1_pdgId", &_llnunu_l1_l1_pdgId);
  _tree_in->SetBranchAddress("llnunu_l1_l1_charge", &_llnunu_l1_l1_charge);
  _tree_in->SetBranchAddress("llnunu_l1_l1_ptErr", &_llnunu_l1_l1_ptErr);

  _tree_in->SetBranchAddress("llnunu_l1_l2_pt", &_llnunu_l1_l2_pt);
  _tree_in->SetBranchAddress("llnunu_l1_l2_eta", &_llnunu_l1_l2_eta);
  _tree_in->SetBranchAddress("llnunu_l1_l2_phi", &_llnunu_l1_l2_phi);
  _tree_in->SetBranchAddress("llnunu_l1_l2_rapidity", &_llnunu_l1_l2_rapidity);
  _tree_in->SetBranchAddress("llnunu_l1_l2_mass", &_llnunu_l1_l2_mass);
  _tree_in->SetBranchAddress("llnunu_l1_l2_pdgId", &_llnunu_l1_l2_pdgId);
  _tree_in->SetBranchAddress("llnunu_l1_l2_charge", &_llnunu_l1_l2_charge);
  _tree_in->SetBranchAddress("llnunu_l1_l2_ptErr", &_llnunu_l1_l2_ptErr);


  // other branches for not light weight tree   
  if (!_useLightTree) {
    std::cout << " Warning: for not _useLightTree, to be implemented." << std::endl;
  }

  // 
  // MC only
  if (!_isData) {
    _tree_in->SetBranchAddress("nTrueInt", &_nTrueInt );
    _tree_in->SetBranchAddress("genWeight",&_genWeight);
    _tree_in->SetBranchAddress("ngenZ", &_ngenZ);
    _tree_in->SetBranchAddress("genZ_pt", _genZ_pt);
  }


  // 2.) Output tree

  // output tree
  _file_out->cd();
  _tree_out = _tree_in->CloneTree(0);

  // sum of events and weights
  _tree_out->Branch("SumEvents", &_SumEvents, "SumEvents/D");
  _tree_out->Branch("SumWeights", &_SumWeights, "SumWeights/D");



  // keep a copy of old branches
  if (_storeOldBranches) {
    _tree_out->Branch("llnunu_mt_old", &_llnunu_mt_old, "llnunu_mt_old/F");
    _tree_out->Branch("llnunu_l1_mass_old", &_llnunu_l1_mass_old, "llnunu_l1_mass_old/F");
    _tree_out->Branch("llnunu_l1_pt_old", &_llnunu_l1_pt_old, "llnunu_l1_pt_old/F");
    _tree_out->Branch("llnunu_l1_phi_old", &_llnunu_l1_phi_old, "llnunu_l1_phi_old/F");
    _tree_out->Branch("llnunu_l1_eta_old", &_llnunu_l1_eta_old, "llnunu_l1_eta_old/F");
    _tree_out->Branch("llnunu_l2_pt_old", &_llnunu_l2_pt_old, "llnunu_l2_pt_old/F");
    _tree_out->Branch("llnunu_l2_phi_old", &_llnunu_l2_phi_old, "llnunu_l2_phi_old/F");
    _tree_out->Branch("llnunu_l1_l1_pt_old", &_llnunu_l1_l1_pt_old, "llnunu_l1_l1_pt_old/F");
    _tree_out->Branch("llnunu_l1_l1_eta_old", &_llnunu_l1_l1_eta_old, "llnunu_l1_l1_eta_old/F");
    _tree_out->Branch("llnunu_l1_l1_phi_old", &_llnunu_l1_l1_phi_old, "llnunu_l1_l1_phi_old/F");
    _tree_out->Branch("llnunu_l1_l1_ptErr_old", &_llnunu_l1_l1_ptErr_old, "llnunu_l1_l1_ptErr_old/F");
    _tree_out->Branch("llnunu_l1_l2_pt_old", &_llnunu_l1_l2_pt_old, "llnunu_l1_l2_pt_old/F");
    _tree_out->Branch("llnunu_l1_l2_eta_old", &_llnunu_l1_l2_eta_old, "llnunu_l1_l2_eta_old/F");
    _tree_out->Branch("llnunu_l1_l2_phi_old", &_llnunu_l1_l2_phi_old, "llnunu_l1_l2_phi_old/F");
    _tree_out->Branch("llnunu_l1_l2_ptErr_old", &_llnunu_l1_l2_ptErr_old, "llnunu_l1_l2_ptErr_old/F");
  }


  return true;
}



// prepare inputs for pu weights
void preparePUWeights()
{
  
  // for each puWeight, read input weight hist and create branch
  for (int i=0; i<(int)_PUTags.size(); i++) {
    sprintf(name, "%s/%s", _PUInputDir.c_str(), _PUInputFileNames.at(i).c_str());
    _PUFiles.push_back(new TFile(name));
    _PUHists.push_back((TH1D*)_PUFiles.at(i)->Get(_PUWeightHistName.c_str()));
    _PUWeights.push_back(new Float_t(1.0));
    sprintf(name, "puWeight%s", _PUTags.at(i).c_str());
    sprintf(name1, "puWeight%s/F", _PUTags.at(i).c_str());
    _tree_out->Branch(name, _PUWeights.at(i),name1);
  }

}

// add more pileup weights
void addPUWeights() 
{
  // for each puWeight, get the weight
  for (int i=0; i<(int)_PUTags.size(); i++){
    *(_PUWeights.at(i)) = _PUHists.at(i)->GetBinContent(_PUHists.at(i)->FindBin(_nTrueInt));
  }
}


// prepare inputs for muon re-calib
void prepareMuonPtRecalib()
{
  std::cout << "addPileupWeights:: to be implemented" << std::endl;
}


// do muon re-calib
void doMuonPtRecalib()
{
  std::cout << "addPileupWeights:: to be implemented" << std::endl;
}


// prepare inputs for addDyZPtWeight
void prepareDyZPtWeight()
{

  // needed output branches
  _tree_out->Branch("ZPtWeight", &_ZPtWeight, "ZPtWeight/F");
  _tree_out->Branch("ZPtWeight_up", &_ZPtWeight_up, "ZPtWeight_up/F");
  _tree_out->Branch("ZPtWeight_dn", &_ZPtWeight_dn, "ZPtWeight_dn/F");
  _tree_out->Branch("ZJetsGenWeight", &_ZJetsGenWeight, "ZJetsGenWeight/F");

}

// addDyZPtWeight
void addDyZPtWeight()
{
  std::cout << "addDyZPtWeight:: to be implemented" << std::endl;
}

// prepare inputs for simple met recoil tune.
void prepareRecoil()
{
  std::cout << "prepareRecoil:: to be implemented" << std::endl;
}

// do simple met recoil fine tuning
void doRecoil()
{
  std::cout << "addPileupWeights:: to be implemented" << std::endl;
}



