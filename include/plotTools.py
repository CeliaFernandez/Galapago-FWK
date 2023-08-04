from ROOT import TCanvas, TLegend,TPie,  TPad, TLine, TLatex, TGraphAsymmErrors, TH1F, THStack, TGraphErrors, TLine, TPaveStats, TGraph, TArrow, TEllipse
import ROOT as r
import os, copy, math, array


#####################
#####
###
###   Function to plot histograms with signal simulation and estimated background
###   (could accept either data-driven or Monte Carlo)
###
###     - Plots are drawn without ratio nor data in signal region (by the moment)
###
#####
#####################

def makeBlindYieldPlot(lumi, hname_SI, hname_bkg, ylog, treeDATA, inputdir, treeSI, xsec = 1.0, rebin = False, lines = [], xlabel = '', outtag = '', ymax = 0.0, LLlabel = ''):


    ### Get histograms
    luminosity = lumi

    print(inputdir, hname_bkg)
    hbkg_ = treeDATA.getLoopTH1F(inputdir, hname_bkg)

    ### rebinins:
    if type(rebin) != bool:
        if len(rebin) > 1:
            hbkg = hbkg_.Rebin(len(rebin)-1, hbkg_.GetName() + '_rebined', rebin)
        else:
            hbkg = hbkg_.Rebin(rebin)
    else:
        hbkg = hbkg_.Clone()

    hbkg.SetFillColorAlpha(r.kCyan-6, 0.8)
    hbkg.SetLineColor(r.kCyan-2)
    hbkg.GetXaxis().SetLabelSize(0.0)
    hbkg.GetXaxis().SetNdivisions(0)
    hbkg.GetXaxis().SetTitleSize(0.045)
    hbkg.GetYaxis().SetTitleSize(0.045)

    ### Signal histograms
    s_histos = []

    hSIS = treeSI.getLoopStack(inputdir, hname_SI)

    for _i, _h in enumerate(hSIS.GetHists()):
        s_histos.append(copy.deepcopy(_h))
        s_histos[-1].Scale(xsec)

    ### Get maximum
    maxValbkg = hbkg.GetMaximum()
    maxValSI = max([s_histos[i].GetMaximum() for i in range(0, len(s_histos))])
    minValSI = min([s_histos[i].GetMinimum() for i in range(0, len(s_histos))])
    maxVal = max([maxValSI, maxValbkg])
    minVal = max([minValSI, 0.0])

    ### Set Maximum
    if not ylog:
        hbkg.SetMaximum(1.3*maxVal)
        hbkg.SetMaximum(0.7*minVal)
    else:
        hbkg.SetMaximum(1e4*maxVal)
        hbkg.SetMinimum(0.01*minVal)

    if ymax: hbkg.SetMaximum(ymax)


    ### Canvas object
    plot = Canvas.Canvas('BlindedYields_'+hname_SI, 'png,pdf', 0.35, 0.65, 0.7, 0.89, 1, lsize = 0.028)
    plot.addHisto(hbkg, 'HIST', 'Background (predicted)', 'f', '', 1, 0)

    ### Add signals:
    for i,_h in enumerate(s_histos):

        _h.SetLineWidth(2) # provisional
        masses = eval(_h.GetTitle()[3:])
        legend = 'H#rightarrowSS (%d GeV,%d GeV,%d mm)'%(masses[0], masses[1], masses[2])
        plot.addHisto(_h, 'HIST, SAME', legend, 'l', _h.GetFillColor(), 1, i+1) # Signal

    if LLlabel == 'EE':
        plot.addLatex(0.13, 0.93, 'Electron channel')
    elif LLlabel == 'MM':
        plot.addLatex(0.13, 0.93, 'Muon channel')

    ### X label
    plot.addLatex(0.44, 0.06, hname_SI.split('_')[1] + ' region', size = 0.045)


    ### Save it
    outdir = os.path.dirname(os.path.abspath(__main__.__file__)) + '/BlindYieldsPlots_' + outtag + '/'
    plot.save(1, 1, ylog, luminosity, '', outputDir = outdir, xlog = False)


