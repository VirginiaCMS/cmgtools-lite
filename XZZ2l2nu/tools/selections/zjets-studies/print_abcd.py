#!/usr/bin/env python

import ROOT
import os
from math import *
from python.TreePlotter import TreePlotter
from python.MergedPlotter import MergedPlotter
from python.StackPlotter import StackPlotter
from python.mylib import *
import python.SetCuts

Channel=raw_input("Please choose a channel (el or mu): \n")
tag0='ZJstudy'
outdir='test'
indir="../AnalysisRegion"
lumi=2.169126704526
LogY=False
doprint=True

if not os.path.exists(outdir): os.system('mkdir '+outdir)

tag = tag0+'_'+'test'
if LogY: tag = tag+'_log'
outTag=outdir+'/'+tag

#  A | C
# ------- dPhi(Z,MET) = 2.5
#  B | D
tex_dic=python.SetCuts.Tex_dic
whichregion=raw_input("Please choose a benchmarck Region (SR or VR): \n")
cuts=python.SetCuts.Cuts(Channel, whichregion)
print cuts

### ----- Initialize (samples):

zjetsPlotters=[]
zjetsSamples = ['DYJetsToLL_M50'] # M50_BIG = M50 + M50_Ext, 150M evts
for sample in zjetsSamples:
    zjetsPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    zjetsPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    zjetsPlotters[-1].addCorrectionFactor('xsec','tree')
    zjetsPlotters[-1].addCorrectionFactor('genWeight','tree')
    zjetsPlotters[-1].addCorrectionFactor('puWeight','tree')
ZJets = MergedPlotter(zjetsPlotters)
ZJets.setFillProperties(1001,ROOT.kGreen+2)

dataPlotters=[]
dataSamples = ['SingleElectron_Run2015C_25ns_16Dec',
               'SingleElectron_Run2015D_16Dec',
               'SingleMuon_Run2015C_25ns_16Dec',
               'SingleMuon_Run2015D_16Dec']
for sample in dataSamples:
    dataPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
Data = MergedPlotter(dataPlotters)

otherBGPlotters=[]
otherBGSamples=['WZTo2L2Q','WZTo3LNu',
                 'ZZTo2L2Nu','ZZTo2L2Q','ZZTo4L',
                 'TTTo2L2Nu',
                 'WWTo2L2Nu','WWToLNuQQ','WZTo1L1Nu2Q']
for sample in otherBGSamples:
    otherBGPlotters.append(TreePlotter(indir+'/'+sample+'.root','tree'))
    otherBGPlotters[-1].addCorrectionFactor('1./SumWeights','tree')
    otherBGPlotters[-1].addCorrectionFactor('xsec','tree')
    otherBGPlotters[-1].addCorrectionFactor('genWeight','tree')
    otherBGPlotters[-1].addCorrectionFactor('puWeight','tree')
otherBG = MergedPlotter(otherBGPlotters)
otherBG.setFillProperties(1001,ROOT.kMagenta)

### ----- Execute (plotting):

ROOT.gROOT.ProcessLine('.x tdrstyle.C')
ROOT.gStyle.SetPadBottomMargin(0.2)
ROOT.gStyle.SetPadLeftMargin(0.15)

#ROOT.TH1.AddDirectory(ROOT.kFALSE) #in this way you could close the TFile after you registe the histograms
yields={}
err={}
err['data']={}
yields['data']={}
err['Zjets']={}
yields['Zjets']={}
err['non-Zjets']={}
yields['non-Zjets']={}

nan=0
fout = ROOT.TFile(outTag+'.root','recreate')

canvas = ROOT.TCanvas('c1', 'c1', 600,630)
canvas.Print(outTag+'.ps[')

for key in tex_dic:
    # MET (data):
    h_met_dt=Data.drawTH1('met_pt',cuts[key],"1",50,0,500,titlex='E_T^{miss}',units='GeV',drawStyle="HIST")
    err['data'][key]=ROOT.Double(0.0)
    yields['data'][key]=h_met_dt.IntegralAndError(0,1+h_met_dt.GetNbinsX(),err['data'][key])
    
    # MET (non-zjets):
    h_met_nonZ=otherBG.drawTH1('met_pt',cuts[key],str(lumi*1000),50,0,500,titlex='E_T^{miss}',units='GeV',drawStyle="HIST")
    err['non-Zjets'][key]=ROOT.Double(0.0)
    yields['non-Zjets'][key]=h_met_nonZ.IntegralAndError(0,1+h_met_nonZ.GetNbinsX(),err['non-Zjets'][key])
    
    # MET (zjets):
    h_met=ZJets.drawTH1('met_pt',cuts[key],str(lumi*1000),50,0,500,titlex='E_T^{miss}',units='GeV',drawStyle="HIST")
    err['Zjets'][key]=ROOT.Double(0.0)
    yields['Zjets'][key]=h_met.IntegralAndError(0,1+h_met.GetNbinsX(),err['Zjets'][key])
    
    canvas.Clear()
    h_met.Draw()
    canvas.Print(outTag+'.ps')
    h_met.SetName("zjets_MET_"+key);
    h_met.Write()

    # deltaR_ll (zjets):
    hdRZmm=ZJets.drawTH1('llnunu_l1_deltaR',cuts[key],str(lumi*1000),50,0,5,titlex='#Delta R(#mu,#mu)',units='',drawStyle="HIST")
    hdRZmm.SetName("zjets_hdRZmm_"+key)
    hdRZmm.GetYaxis().SetTitle("Events")

    canvas.Clear()
    hdRZmm.Draw()
    canvas.Print(outTag+'.ps')
    hdRZmm.Write()

### ----- Finalizing:
canvas.Print(outTag+'.ps]')
os.system('ps2pdf '+outTag+'.ps '+outTag+'.pdf')    
fout.Close()

print yields
print err

### ----- Execute (print):
if doprint:
    err_ac=GetError(yields['Zjets']['SR'], yields['Zjets']['CRc'], err['Zjets']['SR'], err['Zjets']['CRc'])
    err_bd=GetError(yields['Zjets']['CRb'], yields['Zjets']['CRd'], err['Zjets']['CRb'], err['Zjets']['CRd'])
    
    B_data=yields['data']['CRb']-yields['non-Zjets']['CRb']
    D_data=yields['data']['CRd']-yields['non-Zjets']['CRd']
    err_B_data=math.sqrt(err['data']['CRb']**2+err['non-Zjets']['CRb']**2)
    err_D_data=math.sqrt(err['data']['CRd']**2+err['non-Zjets']['CRd']**2)
    err_bd_dt=GetError(B_data, D_data, err_B_data, err_D_data)

    sf={}
    sf['sf_ac_mc']=yields['Zjets']['SR']/yields['Zjets']['CRc']
    sf['sf_bd_mc']=yields['Zjets']['CRb']/yields['Zjets']['CRd']
    sf['sf_bd_dt']=B_data/D_data
    sys_sf={}
    sys_sf['sys_var_corr']=abs(sf['sf_ac_mc']-sf['sf_bd_mc'])/(sf['sf_ac_mc']+sf['sf_bd_mc'])
    sys_sf['sys_nonZ_sub']=abs(sf['sf_bd_mc']-sf['sf_bd_dt'])/(sf['sf_bd_mc']+sf['sf_bd_dt'])
    sys_bd_dt=math.sqrt(sys_sf['sys_var_corr']**2+sys_sf['sys_nonZ_sub']**2)

    
    print '*'*40,'\n{:^40}'.format(Channel+'channel, '+whichregion)
    print '*'*40,'\n{:^40}\n'.format('Scale Factors'),'-'*40 # MC stands for Z(ll)+jets only
    print ' A/C (MC)   | %5.2f +- %5.2f' % (sf['sf_ac_mc'], err_ac)
    print ' B/D (MC)   | %5.2f +- %5.2f' % (sf['sf_bd_mc'], err_bd)
    print ' B/D (data) | %5.2f +- %5.2f +- %5.2f' % (sf['sf_bd_dt'], err_bd_dt, sys_bd_dt)
    print '*'*40,'\n  separate sys. uncert. of B/D (data) \n','-'*40
    print 'sys_var_corr| %.2f \nsys_nonZ_sub| %.2f' % (sys_sf['sys_var_corr'], sys_sf['sys_nonZ_sub'])
    
    # To print the Zjets MC, closure test and data-driven result in SR:
    zjets_mc=yields['Zjets']['CRc']*yields['Zjets']['CRb']/yields['Zjets']['CRd']
    zjets_mc_err=math.sqrt((yields['Zjets']['CRc']*err_bd)**2 + (err['Zjets']['SR']*yields['Zjets']['CRb']/yields['Zjets']['CRd'])**2 )

    C_data=yields['data']['CRc']-yields['non-Zjets']['CRc']
    err_C_data=math.sqrt(err['data']['CRc']**2+err['non-Zjets']['CRc']**2)
    zjets_data=C_data*B_data/D_data
    zjets_data_err=math.sqrt( (C_data*GetError(B_data, D_data))**2+ (err_C_data*B_data/D_data)**2 )

    zjets_data_sys=sys_bd_dt*zjets_data # FIXME
    
    print '*'*40,'\n{:^40}\n'.format('Z+jets estimation'),'-'*40 #'\n Z+jets estimation: \n'
    print '          MC| %5.2f +- %5.2f' % (yields['Zjets']['SR'], err['Zjets']['SR'])
    print 'closure test| %5.2f +- %5.2f' % (zjets_mc, zjets_mc_err)
    print ' data-driven| %5.2f +- %5.2f +- %5.2f' % (zjets_data, zjets_data_err, zjets_data_sys)
    print '*'*40

    # To print the mc and data yields in each region:
    print '-'*40
    print '{:^40}'.format('A | C')
    print '{:^40}'.format('Zjets %6.2f | %6.2f Zjets' % (yields['Zjets']['SR'], yields['Zjets']['CRc']) )
    print '{:^40}'.format('non-Z %6.2f | %6.2f non-Z' % (yields['non-Zjets']['SR'], yields['non-Zjets']['CRc']))
    print '{:^40}'.format('data %6.2f | %6.2f data' % (yields['data']['SR'], yields['data']['CRc']))
    print '-'*40
    print '{:^40}'.format('B | D')
    print '{:^40}'.format('Zjets %6.2f | %6.2f Zjets' % (yields['Zjets']['CRb'], yields['Zjets']['CRd']))
    print '{:^40}'.format('non-Z %6.2f | %6.2f non-Z' % (yields['non-Zjets']['CRb'], yields['non-Zjets']['CRd']))
    print '{:^40}'.format('data %6.2f | %6.2f data' % (yields['data']['CRb'], yields['data']['CRd']))
    print '-'*40