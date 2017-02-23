from ROOT import *
#gROOT.SetBatch()
c=TCanvas()
f1=TFile("/datab/yanchu/backup80X/nonresmc/tt/TTTo2L2Nu/vvTreeProducer/tree.root")
t1=f1.Get("tree")
f2=TFile("/datac/heli/XZZ2/80X_20161029_light_Skim/TTTo2L2Nu.root")
t2=f2.Get("tree")
t1.Draw("llnunu_l1_mass>>h1(26,50,180)","llnunu_l1_l1_pt>55&&llnunu_l1_l2_pt>35")
t2.Draw("llnunu_l1_mass>>h2(26,50,180)","abs(llnunu_l1_l1_pdgId)==11&&llnunu_l1_l1_pt>55&&llnunu_l1_l2_pt>35")
t2.Draw("llnunu_l1_mass>>h3(26,50,180)","abs(llnunu_l1_l1_pdgId)==13&&llnunu_l1_l1_pt>55&&llnunu_l1_l2_pt>35")
h1.Scale(1./h1.Integral())
h2.Scale(1./h2.Integral())
h3.Scale(1./h3.Integral())
h2.SetLineColor(2)
h3.SetLineColor(3)
h1.Draw()
h2.Draw('same')
h3.Draw('same')
c.Print("closure.pdf")