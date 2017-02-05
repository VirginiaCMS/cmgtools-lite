#! /usr/bin/env python
from ROOT import *
from array import array
gROOT.SetBatch()
gStyle.SetOptStat(000000)
gStyle.SetLegendBorderSize(0)
c=TCanvas('c','c',2000,800)
leg=TLegend(.7,.7,.89,.89)
inputfile=TFile("/afs/cern.ch/work/y/yanchu/public/analysis/Moriond/CMSSW_8_0_25/src/CMGTools/XZZ2l2nu/tools/metCorr/all.root")

tree=inputfile.Get("tree")
xbin=[20,30,40,50,60,70,80,100,130,200]
ehist=TH1F("ehist","ehist",len(xbin)-1,array('d',xbin))
#ehist=TH1F("ehist","ehist",30,20,200)
mhist=ehist.Clone("mhist")
tree.Draw("llnunu_l1_l1_pt>>ehist","(abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>20&&llnunu_l1_l2_pt>20)*mscale")
tree.Draw("llnunu_l1_l2_pt>>+ehist","(abs(llnunu_l1_l2_pdgId)==11&&llnunu_l1_l1_pt>20&&llnunu_l1_l2_pt>20)*mscale")
tree.Draw("llnunu_l1_l1_pt>>mhist","(abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>20&&llnunu_l1_l2_pt>20)*mscale")
tree.Draw("llnunu_l1_l2_pt>>+mhist","(abs(llnunu_l1_l2_pdgId)==13&&llnunu_l1_l1_pt>20&&llnunu_l1_l2_pt>20)*mscale")

mhist.SetLineColor(2)
mhist.SetTitle(";pT/GeV")
mhist.SetMaximum(2e5)
leg.AddEntry(ehist,'electron')
leg.AddEntry(mhist,'muon')
c.Divide(2,1)
c.cd(1).SetLogy()
mhist.Draw()
ehist.Draw('same')
leg.Draw('same')
meratio=ehist.Clone("meratio")
meratio.Sumw2()
meratio.Divide(mhist)
meratio.SetTitle(";pT/GeV;electron/muon")
c.cd(2)
meratio.Draw("e")
c.Print('~/www/meratioweight.pdf')

