import ROOT as r
from   ROOT import gROOT, TCanvas, TFile, TGraphErrors, SetOwnership
import math, sys, optparse, array, copy, os
import gc, inspect, __main__
import numpy as np
import time
import shutil

import include.Sample as Sample
import include.Launcher as Launcher
import include.helper as helper
import include.Canvas as Canvas
import include.CutManager as CutManager


################################# GLOBAL VARIABLES DEFINITION ####################################

runningfile = os.path.abspath(__file__)
WORKPATH = ''
for level in runningfile.split('/')[:-1]: 
    WORKPATH += level
    WORKPATH += '/'

if __name__ == "__main__":

    parser = optparse.OptionParser(usage='usage: %prog [opts] FilenameWithSamples', version='%prog 1.0')
    parser.add_option('-i', '--input', action='store', type=str, dest='input', default='', help='Target directory')
    (opts, args) = parser.parse_args()

    ############# Set the TDR plot style
    r.gROOT.LoadMacro(WORKPATH + 'include/tdrstyle.C')
    r.gROOT.SetBatch(1)
    r.setTDRStyle()

    ############# Dat file
    filename = 'dat/Samples_cern_Muon_provisionalRun3.dat'

    ############# Data definition
    Muons = ['Muon0_Run2023C']
    Background = ['DYJetsToLL_M-50_22EE']

    ############# Luminosity definition
    lumi = 35.9 # fb-1

    ############# Load Tree
    treeDATA = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Muons, 'DATA'), name = 'DATA', isdata = 1 )
    treeMC   = Sample.Tree( fileName = helper.selectSamples(WORKPATH + filename, Background, 'DATA'), name = 'MC', isdata = 0 )
    treeDATA.setDefinitions(config = 'definitions/defs_MuonVariables.json')
    treeMC.setDefinitions(config = 'definitions/defs_MuonVariables.json')

    ############# Fill histograms
    #histo = treeDATA.getTH1F(lumi = lumi, name = 'Muon_pt', var = 'Muon_pt', nbin = 100, xmin = 0.0, xmax = 100., cut = '', options = '', xlabel = 'Muon p_{T} (GeV)')
    #histoDATA = treeDATA.getTH1F(lumi = lumi, name = 'DiMuon_mass', var = 'DiMuon_mass', cut = '(nMuon > 1 && Muon_charge[0]*Muon_charge[1] < 0.0)', options = '', xlabel = 'Dimuon mass (GeV)', xaxis = np.logspace(-1, 3, 1000))
    histoDATA = treeDATA.getTH1F(lumi = lumi, name = 'DiMuon_mass', var = 'DiMuon_mass', nbin = 100, xmin = 40.0, xmax = 150.0, cut = '(nMuon > 1 && Muon_charge[0]*Muon_charge[1] < 0.0)', options = '', xlabel = 'Dimuon mass (GeV)')
    histoMC   = treeMC.getTH1F(lumi = lumi, name = 'DiMuon_mass', var = 'DiMuon_mass', nbin = 100, xmin = 40.0, xmax = 150.0, cut = '(nMuon > 1 && Muon_charge[0]*Muon_charge[1] < 0.0)', options = '', xlabel = 'Dimuon mass (GeV)')

    ############# Plot histogram
    canvas = Canvas.Canvas(histoMC.GetName(), 'png,pdf', 0.16, 0.72, 0.56, 0.82, 1)
    canvas.addHisto(histoMC,'HIST', 'MC', 'l', r.kBlue, True, 0)
    canvas.addHisto(histoDATA,'P, SAME', 'Data', 'l', r.kBlack, True, 1)
    canvas.save(1, 0, 0, '', '', outputDir = './', isPrivate = False, xlog = False)


