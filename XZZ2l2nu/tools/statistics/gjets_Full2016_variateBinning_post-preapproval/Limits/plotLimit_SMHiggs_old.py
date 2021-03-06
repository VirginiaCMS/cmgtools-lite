from ROOT import *
from tdrStyle import *
setTDRStyle()
        
from array import array
import ROOT, os, string
import math
from math import *

import sys, os, pwd, commands
from subprocess import *
import optparse, shlex, re

# load signal modules
sys.path.append('./SigInputs')

grootargs = []
def callback_rootargs(option, opt, value, parser):
    grootargs.append(opt)

def parseOptions():

    usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
    parser = optparse.OptionParser(usage)
    parser.add_option("--cutChain", dest="CutChain",default='zpt100met50',help="kinematic cuts on zpt and met")
    parser.add_option("--dyMC",action="store_true", dest="DyMC", default=False,help="whether use DY MC to predict Z+jets MT spectrum")
    parser.add_option("--perBinStatUnc",action="store_true", dest="PerBinStatUnc", default=False,help="whether do per-bin statistical uncertainty")
    parser.add_option("--obs", dest="Observable",default='mT',help="template observable")

    parser.add_option("-l",action="callback",callback=callback_rootargs)
    parser.add_option("-q",action="callback",callback=callback_rootargs)
    parser.add_option("-b",action="callback",callback=callback_rootargs)

    return parser

def plotLimit(parser):

    mass = array('d',[])
    zeros = array('d',[])
    exp_p2 = array('d',[])
    exp_p1 = array('d',[])
    exp = array('d',[])
    exp_m1 = array('d',[])
    exp_m2 = array('d',[])

    exp_VBF = array('d',[])

    (options,args) = parser.parse_args()

    cut=options.CutChain
    perBinStatUnc=options.PerBinStatUnc
    isDyMC=options.DyMC
    isGJets= (not isDyMC)

    zjetsMethod = 'mczjet'
    if( not isDyMC) :
      zjetsMethod = 'gjet'

    observable=options.Observable

    outdir_pre = zjetsMethod+'_'+cut
    if(perBinStatUnc) :
      outdir = outdir + '_perBinStatUnc'

    mH = [500,600,700,800,900,1000,1500,2500]
    for m in mH:
      #scale = 24.7549/1000 #1/((3.363+3.366+3.370)*0.01*0.2*2)
      scale = 1/((3.363+3.366)*0.01*0.2*2)/1000
      mass.append(float(m))
      zeros.append(0.0)

      sigType = 'ggH'
      outdir = outdir_pre + '_' + sigType + '_' + observable

      f = TFile("Datacards/"+outdir+"/higgsCombineTest.Asymptotic.mH"+str(m)+".root","READ")
      t = f.Get("limit")
    
      t.GetEntry(2)
      thisexp = t.limit*scale
      exp.append(thisexp)
      t.GetEntry(0)
      exp_m2.append(thisexp-t.limit*scale)
      t.GetEntry(1)
      exp_m1.append(thisexp-t.limit*scale)
      t.GetEntry(3)
      exp_p1.append(t.limit*scale-thisexp)
      t.GetEntry(4)
      exp_p2.append(t.limit*scale-thisexp)

      sigType = 'VBF'
      outdir = outdir_pre + '_' + sigType + '_' + observable

      f = TFile("Datacards/"+outdir+"/higgsCombineTest.Asymptotic.mH"+str(m)+".root","READ")
      t = f.Get("limit")

      t.GetEntry(2)
      thisexp = t.limit*scale
      exp_VBF.append(thisexp)

    print 'mass',mass
    print 'exp',exp
 
    v_mass = TVectorD(len(mass),mass)
    v_zeros = TVectorD(len(zeros),zeros)
    v_exp_p2 = TVectorD(len(exp_p2),exp_p2)
    v_exp_p1 = TVectorD(len(exp_p1),exp_p1)
    v_exp = TVectorD(len(exp),exp)
    v_exp_m1 = TVectorD(len(exp_m1),exp_m1)
    v_exp_m2 = TVectorD(len(exp_m2),exp_m2)

    v_exp_VBF = TVectorD(len(exp_VBF),exp_VBF)

    c = TCanvas("c","c",800, 800)
    c.SetLogy()
    #c.SetLogx()
    c.SetGridx()
    c.SetGridy()

    c.SetRightMargin(0.06)
    c.SetLeftMargin(0.2)

    dummy = TH1D("dummy","dummy", 1, 500,2500)
    dummy.SetBinContent(1,0.0)
    dummy.GetXaxis().SetTitle('m_{H} [GeV]')   
    dummy.GetYaxis().SetTitle('#sigma(pp#rightarrowH)#timesBR(H#rightarrowZZ) [pb]')   
    dummy.SetLineColor(0)
    dummy.SetLineWidth(0)
    dummy.SetFillColor(0)
    dummy.SetMinimum(0.001)
    dummy.SetMaximum(10.0)
    dummy.GetXaxis().SetMoreLogLabels(kTRUE)
    dummy.GetXaxis().SetNoExponent(kTRUE)
    dummy.Draw()

    gr_exp2 = TGraphAsymmErrors(v_mass,v_exp,v_zeros,v_zeros,v_exp_m2,v_exp_p2)
    gr_exp2.SetLineColor(kYellow)
    gr_exp2.SetFillColor(kYellow)
    gr_exp2.Draw("e3same")

    gr_exp1 = TGraphAsymmErrors(v_mass,v_exp,v_zeros,v_zeros,v_exp_m1,v_exp_p1)
    gr_exp1.SetLineColor(kGreen)
    gr_exp1.SetFillColor(kGreen)
    gr_exp1.Draw("e3same")

    gr_exp = TGraphAsymmErrors(v_mass,v_exp,v_zeros,v_zeros,v_zeros,v_zeros)
    gr_exp.SetLineColor(1)
    gr_exp.SetLineWidth(2)
    gr_exp.SetLineStyle(2)
    gr_exp.Draw("csame")

    gr_exp_VBF = TGraphAsymmErrors(v_mass,v_exp_VBF,v_zeros,v_zeros,v_zeros,v_zeros)
    gr_exp_VBF.SetLineColor(kRed)
    gr_exp_VBF.SetLineWidth(2)
    gr_exp_VBF.SetLineStyle(2)
    gr_exp_VBF.Draw("csame")

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.5*c.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    latex2.DrawLatex(0.87, 0.95,"35.87 fb^{-1} (13 TeV)")
    latex2.SetTextSize(0.7*c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.25, 0.85, "CMS")
    latex2.SetTextSize(0.5*c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    latex2.DrawLatex(0.25, 0.8, "Preliminary")

    legend = TLegend(.45,.70,.90,.90)
    legend.AddEntry(gr_exp_VBF , "100% VBF production", "l")
    legend.AddEntry(gr_exp , "100% ggH production", "l")
    legend.AddEntry(gr_exp1 , "#pm 1#sigma", "f")
    legend.AddEntry(gr_exp2 , "#pm 2#sigma", "f")
    legend.SetShadowColor(0)
    legend.SetFillColor(0)
    legend.SetLineColor(0)            
    legend.Draw("same")
                                                            
    gPad.RedrawAxis()

    c.SaveAs("xzz2l2nu_limit_13TeV_channel_SR_SMHiggs.pdf")
    c.SaveAs("xzz2l2nu_limit_13TeV_channel_SR_SMHiggs.png")


def Run():

    parser=parseOptions()

    plotLimit(parser)

if __name__ == "__main__":
    Run()
