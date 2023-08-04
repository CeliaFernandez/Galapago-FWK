from ROOT import TCanvas, TLegend,TPie,  TPad, TLine, TLatex, TGraphAsymmErrors, TH1F, THStack, TGraphErrors, TLine, TPaveStats, TGraph, TArrow, TEllipse
import ROOT as r
import os, copy, math, array
from array import array
import time
import numpy as np

class Canvas:
   'Common base class for all Samples'

   def __init__(self, name, _format, x1, y1, x2, y2, c, ww=0, hh=0, lsize=0):
      self.name = name
      self.format = _format
      self.plotNames    = [name + "." + i for i in _format.split(',')]
      self.plotNamesLog = [name + "_log." + i for i in _format.split(',')]
      self.myCanvas = TCanvas(name, name) if not ww else TCanvas(name, name, ww, hh)
      self.ToDraw = []
      self.orderForLegend = []
      self.histos = []
      self.lines = []
      self.arrows= []
      self.ellipses = []
      self.latexs= []
      self.bands = []
      self.options = []
      self.labels = []      
      self.labelsOption = []
      self.myLegend = TLegend(x1, y1, x2, y2)
      self.myLegend.SetFillStyle(0)
      self.myLegend.SetTextFont(42)
      if lsize:
         self.myLegend.SetTextSize(lsize)
      else:
         self.myLegend.SetTextSize(0.03)
      self.myLegend.SetLineWidth(0)
      self.myLegend.SetBorderSize(0)
      self.myLegend.SetNColumns(c)              

   def changeLabelsToNames(self):
      newlabels = []
      for il,lab in enumerate(self.labels):
         print('changing label %s to %s'%(lab, self.histos[il].GetName()))
         newlabels.append(self.histos[il].GetName())
      self.labels = newlabels

   def banner(self, isData, lumi, scy, inProgress = False, isPrivate = False):
     
      latex = TLatex()
      latex.SetNDC();                         
      latex.SetTextAngle(0);                  
      latex.SetTextColor(r.kBlack);           
      latex.SetTextFont(42);                  
      latex.SetTextAlign(31);                 
      latex.SetTextSize(0.06);                
      if isPrivate:
          latex.SetTextSize(0.04);                
          latex.SetTextAlign(11);                 
          if not scy:
              latex.DrawLatex(0.13, 0.93, "#bf{Private work}") 
          else:               
              latex.DrawLatex(0.22, 0.93, "#bf{Private work}") 
      else:
          if not scy:
              latex.DrawLatex(0.25, 0.93, "#bf{CMS}") 
          else:               
              latex.DrawLatex(0.34, 0.93, "#bf{CMS}") 

         
      latexb = TLatex()                      
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(11);
      latexb.SetTextSize(0.04);            

      if isPrivate:
         latexb.SetTextSize(0.03);            
         if not scy:
             latexb.DrawLatex(0.35, 0.93, "#it{CMS data/simulation}")
         else:
             latexb.DrawLatex(0.44, 0.93, "#it{CMS data/simulation}")
      elif inProgress:
         if not scy:
             latexb.DrawLatex(0.26, 0.93, "#it{Work in progress}")
         else:
             latexb.DrawLatex(0.35, 0.93, "#it{Work in progress}")
      elif(isData):
         if not scy:
             latexb.DrawLatex(0.26, 0.93, "#it{Preliminary}")
         else:
             latexb.DrawLatex(0.35, 0.93, "#it{Preliminary}")
      else:
         if not scy:
             latexb.DrawLatex(0.26, 0.93, "#it{Simulation}")
         else:
             latexb.DrawLatex(0.35, 0.93, "#it{Simulation}")

      """
      if(isData):
         if not scy:
             latexb.DrawLatex(0.43, 0.93, "#it{Preliminary}")
         else:
             latexb.DrawLatex(0.53, 0.93, "#it{Preliminary}")
      else:
         if not inProgress:
             if not scy:
                 latexb.DrawLatex(0.43, 0.93, "#it{Simulation}")
             else:
                 latexb.DrawLatex(0.53, 0.93, "#it{Simulation}")
         else:
             if not scy:
                 latexb.DrawLatex(0.54, 0.93, "#it{Work in progress}")
             else:
                 latexb.DrawLatex(0.63, 0.93, "#it{Work in progress}")
      """

      text_lumi = ''
      if lumi: text_lumi = str(lumi)+" fb^{-1}  (13 TeV)"
     
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.04);
      latexc.DrawLatex(0.90, 0.93, text_lumi)                

   def bannerInFrame(self, isData, lumi, inProgress = False, isPrivate = False):
     
      latex = TLatex()
      latex.SetNDC();                         
      latex.SetTextAngle(0);                  
      latex.SetTextColor(r.kBlack);           
      latex.SetTextFont(42);                  
      latex.SetTextAlign(11);                 
      latex.SetTextSize(0.05);                
      if isPrivate:
          latex.DrawLatex(0.17, 0.84, "#bf{Private work}") 
      else:
          latex.DrawLatex(0.17, 0.84, "#bf{CMS}") 

         
      latexb = TLatex()                      
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(11);
      latexb.SetTextSize(0.033);            

      if inProgress:
          latexb.DrawLatex(0.17, 0.8, "#it{Work in progress}")
      elif isPrivate:
          latexb.DrawLatex(0.17, 0.8, "#it{CMS data/simulation}")
      elif(isData):
          latexb.DrawLatex(0.17, 0.8, "#it{Preliminary}")
      else:
          latexb.DrawLatex(0.17, 0.8, "#it{Simulation}")

      text_lumi = ''
      if lumi: text_lumi = str(lumi)+" fb^{-1}  (13 TeV)"
     
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.04);
      latexc.DrawLatex(0.90, 0.93, text_lumi)                

   def bannerRatio(self, isData, lumi, scy = False, inProgress = False, isPrivate = False):
    
      latex = TLatex()
      latex.SetNDC();
      latex.SetTextAngle(0);
      latex.SetTextColor(r.kBlack);
      latex.SetTextFont(42);
      latex.SetTextAlign(11);

      """
      if not scy:
          latex.DrawLatex(0.23, 0.88, "#bf{CMS}")
      else:
          latex.DrawLatex(0.30, 0.88, "#bf{CMS}")
      """

      if isPrivate:
          latex.SetTextSize(0.055);
          if not scy:
              latex.DrawLatex(0.13, 0.88, "#bf{Private work}") 
          else:               
              latex.DrawLatex(0.19, 0.88, "#bf{Private work}") 
      else:
          latex.SetTextSize(0.068);
          if not scy:
              latex.DrawLatex(0.13, 0.88, "#bf{CMS}") 
          else:               
              latex.DrawLatex(0.19, 0.88, "#bf{CMS}") 

      latexb = TLatex()
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.045);            

      if(isData) and not inProgress and not isPrivate:
         if not scy:
             latexb.DrawLatex(0.39, 0.88, "#it{Preliminary}")
         else:
             latexb.DrawLatex(0.46, 0.88, "#it{Preliminary}")
      else:
         if not inProgress and not isPrivate:
             if not scy:
                 latexb.DrawLatex(0.37, 0.88, "#it{Simulation}")
             else:
                 latexb.DrawLatex(0.44, 0.88, "#it{Simulation}")
         elif inProgress:
             if not scy:
                 latexb.DrawLatex(0.46, 0.88, "#it{Work in progress}")
             else:
                 latexb.DrawLatex(0.53, 0.88, "#it{Work in progress}")
         elif isPrivate:
             if not scy:
                 latexb.DrawLatex(0.61, 0.88, "#it{CMS data/simulation}")
             else:
                 latexb.DrawLatex(0.67, 0.88, "#it{CMS data/simulation}")



      text_lumi =str(lumi)+" fb^{-1} (13 TeV)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);
      if lumi != '': latexc.DrawLatex(0.91, 0.88, text_lumi)


   def banner3(self, isData, lumi):
     
      latex = TLatex()
      latex.SetNDC();                         
      latex.SetTextAngle(0);                  
      latex.SetTextColor(r.kBlack);           
      latex.SetTextFont(42);                  
      latex.SetTextAlign(31);                 
      latex.SetTextSize(0.07);                
      latex.DrawLatex(0.1, 1.22, "#bf{CMS}") 
               
      latexb = TLatex()                      
      latexb.SetNDC();
      latexb.SetTextAngle(0);
      latexb.SetTextColor(r.kBlack);
      latexb.SetTextFont(42);
      latexb.SetTextAlign(31);
      latexb.SetTextSize(0.05);            
                                                             
      #if(isData):
      #  latexb.DrawLatex(0.34, 1.22, "#it{Preliminary}")
      #else:
      #  latexb.DrawLatex(0.34, 1.22, "#it{Simulation}")
                                                             
      text_lumi = str(lumi) + " fb^{-1} (13 TeV)"
      latexc = TLatex()
      latexc.SetNDC();
      latexc.SetTextAngle(0);
      latexc.SetTextColor(r.kBlack);
      latexc.SetTextFont(42);
      latexc.SetTextAlign(31);
      latexc.SetTextSize(0.05);

   def addBand(self, x1, y1, x2, y2, color, opacity):

      grshade = TGraph(4)
      grshade.SetPoint(0,x1,y1)
      grshade.SetPoint(1,x2,y1)
      grshade.SetPoint(2,x2,y2)
      grshade.SetPoint(3,x1,y2)
      #grshade.SetFillStyle(3001)
      grshade.SetFillColorAlpha(color, opacity)
      self.bands.append(grshade)

   def addLine(self, x1, y1, x2, y2, color, thickness = 0.):
      line = TLine(x1,y1,x2,y2)
      line.SetLineColor(color)
      line.SetLineStyle(2)
      if thickness:
          line.SetLineWidth(thickness)
      self.lines.append(line)

   def addArrow(self, x1, y1, x2, y2, color, option, thickness = 0.):
      arrow = TArrow(x1,y1,x2,y2, 0.05, option)
      arrow.SetLineColor(color)
      if thickness:
          arrow.SetLineWidth(thickness)
      self.arrows.append(arrow)

   def addEllipse(self, x0, y0, r1, r2, color, fillcolor = False, thickness = 0.):
      ellipse = TEllipse(x0, y0, r1, r2)
      ellipse.SetLineColor(color)
      if fillcolor:
          ellipse.SetFillColor(fillcolor)
      else:
          ellipse.SetFillStyle(0)
      if thickness:
          ellipse.SetLineWidth(thickness)
      self.ellipses.append(ellipse)


   def addLatex(self, x1, y1, text, font=42, size = 0.04, align = 11, color = r.kBlack):
      lat = [x1, y1, text, font, size, align, color]
      self.latexs.append(lat)

   def makeOFHisto(self, histo):
      nbinsx = histo.GetNbinsX()
      xmin = histo.GetXaxis().GetXmin(); xmax = histo.GetXaxis().GetXmax();
      newhisto = r.TH1F(histo.GetName() +'_withOFBin', 'withOFBin' + histo.GetTitle(), nbinsx+1, xmin, xmax+(xmax-xmin)/nbinsx)
      newhisto.Sumw2()
      newhisto.SetMarkerColor(histo.GetMarkerColor())
      newhisto.SetMarkerStyle(histo.GetMarkerStyle())
      newhisto.SetMarkerSize (histo.GetMarkerSize ())
      newhisto.SetLineColor(histo.GetLineColor())
      newhisto.SetLineStyle(histo.GetLineStyle())
      newhisto.SetLineWidth(histo.GetLineWidth())
      newhisto.SetMaximum(histo.GetMaximum())
      newhisto.SetMinimum(histo.GetMinimum())
      newhisto.GetXaxis().SetTitle(histo.GetXaxis().GetTitle())
      newhisto.GetYaxis().SetTitle(histo.GetYaxis().GetTitle())
      for i in range(1,nbinsx+2):
         newhisto.SetBinContent(i,histo.GetBinContent(i))
         newhisto.SetBinError  (i,histo.GetBinError  (i))
      return newhisto
        
   def makeRate(self, eff, option, is2d = False, ymin = 0.0, ymax = 1.2, ratio = False):

      eff.Draw(option)
      self.myCanvas.Update()

      _g = eff.GetPaintedGraph()

      _g.SetMinimum(ymin)
      _g.SetMaximum(ymax)
      _h = eff.GetTotalHistogram()
      xmax = _h.GetXaxis().GetBinUpEdge(_h.GetNbinsX())
      xmin = _h.GetXaxis().GetBinLowEdge(1)
      _g.GetXaxis().SetLimits(xmin,xmax)
      if ratio:
          _g.GetXaxis().SetLabelSize(0)
          _g.GetYaxis().SetTitleSize(0.045)
          _g.GetYaxis().SetTitleOffset(1.25);
          _g.GetYaxis().SetLabelSize(0.04)

      return eff

 
   def addHisto(self, h_, option, label, labelOption, color, ToDraw, orderForLegend, marker = False, doOF = False, normed = False):

      h = copy.deepcopy(h_)

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
          h.SetFillColorAlpha(r.kWhite, 0)
      if(label == ""):
          label = h.GetTitle()

      if normed:
          h.Scale(1.0/h.Integral())
      if marker:
          h.SetMarkerStyle(marker)
      h.SetLineWidth(2)
      self.histos.append(h if not doOF else self.makeOFHisto(h))
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)

   def addRate(self, eff, option, label, labelOption, color, ToDraw, orderForLegend, marker = False):

      if(label == ""):
          label = eff.GetTitle()

      _eff = copy.deepcopy(eff)
      if marker:
          _eff.SetMarkerStyle(marker)
      else:
          _eff.SetMarkerStyle(21)

      if color:
          _eff.SetMarkerColor(color)

      _eff.SetLineWidth(1)
      _eff.SetLineColor(color)
      _eff.SetMarkerSize(1.0)

      self.histos.append(_eff)
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend) 

   def add2DRate(self, eff, option, zmin, zmax):

      _eff = copy.deepcopy(eff)

      _h = _eff.CreateHistogram()
      if zmax:
          _h.GetZaxis().SetRangeUser(zmin, zmax)
      self.histos.append(_h)
      self.options.append(option)
      self.labels.append('')
      self.labelsOption.append('')
      self.ToDraw.append(True)
      self.orderForLegend.append(0) 


   def addProf(self, prof, option, label, labelOption, color, ToDraw, orderForLegend, marker = False):

      if(label == ""):
          label = prof.GetTitle()

      _prof = copy.deepcopy(prof)
      if marker: 
          _prof.SetMarkerStyle(marker)
      else:
          _prof.SetMarkerStyle(21)
      _prof.SetMarkerColor(color)
      _prof.SetLineWidth(2)
      _prof.SetLineColor(color)
      _prof.SetMarkerSize(1.0)

      self.histos.append(_prof)
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)



   def addGraph(self, h, option, label, labelOption, color, ToDraw, orderForLegend):

      if(color != ""):
          h.SetLineColor(color)
          h.SetMarkerColor(color)
      if(label == ""):
          label = h.GetTitle()

      self.histos.append(h)
      self.options.append(option)
      self.labels.append(label)
      self.labelsOption.append(labelOption)
      self.ToDraw.append(ToDraw)
      self.orderForLegend.append(orderForLegend)


   def addStack(self, h, option, ToDraw, orderForLegend):

      legendCounter = orderForLegend
      if(orderForLegend < len(self.orderForLegend)):
          legendCounter = len(self.orderForLegend)

      self.addHisto(h, option, "", "", "", ToDraw, -1)  
      for h_c_ in h.GetHists():
          h_c = copy.deepcopy(h_c_)
          self.addHisto(h_c, "H", h_c.GetTitle(), "F", "", 0, legendCounter)
          legendCounter = legendCounter + 1                                          
       
   def addPies(self, h, option, ToDraw, orderForLegend):
   
      legendCounter = orderForLegend
      if(orderForLegend < len(self.orderForLegend)):
          legendCounter = len(self.orderForLegend)
   
      self.addHisto(h, option, "", "", "", ToDraw, -1)  
      for h_c in h.GetHists():
          self.addHisto(h_c, "H", h_c.GetTitle(), "F", "", 0, legendCounter)
          legendCounter = legendCounter + 1                                       

   def makeLegend(self):

      for i in range(0, len(self.histos)):
          for j in range(0, len(self.orderForLegend)):
              if(self.orderForLegend[j] != -1 and self.orderForLegend[j] == i):
                  self.myLegend.AddEntry(self.histos[j], self.labels[j], self.labelsOption[j])
   

   def ensurePath(self, _path):

      ## Enter a while loop to avoid race conditions
      while True:
          d = os.path.dirname(_path)
          try:
              if not os.path.exists(d):
                  os.makedirs(d)
              break
          except OSError as e:
              if e.errno != os.errno.EEXIST:
                  raise
              print("Sleeping...")
              time.sleep(1.0)
              pass

   def saveRatio(self, legend, isData, log, lumi, hdata, hMC, r_ymin=0, r_ymax=2, r_xmin=0, r_xmax=0, label ="Data/Prediction", hsys = None, outputDir = 'plots/', xlog = False, maxYnumbers = False, inProgress = False, isPrivate = False):

      self.myCanvas.cd()

      ## Create tha pads
      pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0) # for the plot
      pad1.SetBottomMargin(0.015)
      pad1.SetTopMargin(0.13)
      pad1.Draw()                                     
      pad2 = TPad("pad2", "pad2", 0, 0.01, 1, 0.3) # for the ratio
      pad2.SetTopMargin(0.05);
      pad2.SetBottomMargin(0.4);
      #pad2.SetGridy(1);
      pad2.Draw();                                      

      ## Set log scale
      pad1.cd()
      if(log):
          pad1.SetLogy(1)
      if(xlog):
          pad1.SetLogx(1)
          pad2.SetLogx(1)

      ## Draw histograms
      maxYAxisValue = 0.0
      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):
              if str(type(self.histos[i])) == "<class 'ROOT.TEfficiency'>":
                  self.makeRate(self.histos[i], self.options[i], ratio = True)
              else:
                  self.histos[i].GetYaxis().SetTitleSize(0.045)
                  self.histos[i].GetYaxis().SetTitleOffset(1.25);
                  self.histos[i].GetXaxis().SetLabelSize(0)
                  self.histos[i].GetYaxis().SetLabelSize(0.04)
                  self.histos[i].Draw(self.options[i])
                  if self.histos[i].GetMaximum() > maxYAxisValue:
                      maxYAxisValue = self.histos[i].GetMaximum()

      ## Draw frame 1 again
      pad1.Update()
      pad1.RedrawAxis()
      aux_frame = TLine()
      aux_frame.SetLineWidth(2) 
      #aux_frame.DrawLine(pad1.GetUxmax(), pad1.GetUymin(), pad1.GetUxmax(), maxYAxisValue);
      aux_frame.DrawLine(pad1.GetUxmax(), pad1.GetUymin(), pad1.GetUxmax(), pad1.GetUymax());


      if(legend):
          self.makeLegend()
          self.myLegend.SetTextSize(0.04) # Modify the legend size
          self.myLegend.Draw()

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextColor(latex[-1])
          lat.SetTextAlign(latex[-2])
          lat.SetTextSize(latex[-3])
          lat.SetTextFont(latex[-4])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      
      ## Ratios
      if type(hMC) != list:
          hMClist = [hMC]
      else: hMClist = hMC

      ratios = []

      for tmp_hMC in hMClist:
          ind = hMClist.index(tmp_hMC)

          if str(type(tmp_hMC)) == "<class 'ROOT.TEfficiency'>":
              tmp_den = tmp_hMC.GetTotalHistogram().Clone()
              tmp_num = hdata.GetTotalHistogram().Clone()
              for n in range(0,tmp_num.GetNbinsX()):
                  tmp_den.SetBinContent(n+1, tmp_hMC.GetEfficiency(n+1))
                  tmp_num.SetBinContent(n+1, hdata.GetEfficiency(n+1))
                  tmp_den.SetBinError(n+1, tmp_hMC.GetEfficiencyErrorLow(n+1))
                  tmp_num.SetBinError(n+1, hdata.GetEfficiencyErrorLow(n+1))
              tmp_ratio = tmp_num.Clone(tmp_hMC.GetName()+'_ratio')
              tmp_ratio.Divide(tmp_den)
          else:
              tmp_ratio = hdata.Clone(tmp_hMC.GetName()+'_ratio')
              tmp_ratio.Divide(tmp_hMC)


          ## Ratio tunning
          tmp_ratio.SetTitle("")
          tmp_ratio.GetYaxis().SetRangeUser(r_ymin, r_ymax);
          tmp_ratio.GetYaxis().SetTitle(label);
          tmp_ratio.GetYaxis().CenterTitle();
          tmp_ratio.GetYaxis().SetLabelSize(0.10);
          tmp_ratio.GetYaxis().SetNdivisions(4);
          tmp_ratio.GetYaxis().SetTitleOffset(0.5);
          tmp_ratio.GetXaxis().SetLabelSize(0.10);
          tmp_ratio.GetYaxis().SetTitleSize(0.11);
          tmp_ratio.GetXaxis().SetTitleSize(0.12);
          tmp_ratio.GetXaxis().SetLabelOffset(0.02);
          if 'TEfficiency' not in str(type(self.histos[0])): 
              tmp_ratio.GetXaxis().SetTitle(self.histos[0].GetXaxis().GetTitle());
          else:
              tmp_ratio.GetXaxis().SetTitle(self.histos[0].GetTotalHistogram().GetXaxis().GetTitle());
          tmp_ratio.SetMarkerStyle(20);
          tmp_ratio.SetMarkerColor(r.kBlack);
          tmp_ratio.SetMarkerSize(0.8);
          tmp_ratio.SetMarkerColor(r.kBlack if len(hMClist) == 1 else tmp_hMC.GetMarkerColor());
          tmp_ratio.SetLineColor  (r.kBlack if len(hMClist) == 1 else tmp_hMC.GetLineColor  ());
          tmp_ratio.SetLineColor(r.kBlack);
          tmp_ratio.SetLineWidth(2);
          tmp_ratio.SetLineStyle(tmp_hMC.GetLineStyle())

          ratios.append(tmp_ratio)
          xmin = tmp_ratio.GetBinLowEdge(1)
          xmax = tmp_ratio.GetBinLowEdge(tmp_ratio.GetNbinsX()+1)

      pad2.cd();  

      ## Draw systematics (if included)
      ratios[0].Draw('AXIS')
      if hsys is not None:
          hsys.GetYaxis().SetTitle(label);
          hsys.GetYaxis().CenterTitle();
          hsys.GetYaxis().SetLabelSize(0.10);
          hsys.GetYaxis().SetNdivisions(4);
          hsys.GetYaxis().SetTitleOffset(0.5);
          hsys.GetXaxis().SetLabelSize(0.10);
          hsys.GetYaxis().SetTitleSize(0.11);
          hsys.GetXaxis().SetTitleSize(0.12);
          hsys.GetXaxis().SetLabelOffset(0.02);
          hsys.SetLineColor(r.kBlue-10)
          hsys.SetFillColor(r.kBlue-10)
          hsys.SetMarkerSize(0)
          #hsys.SetFillStyle(3013)
          hsys.GetYaxis().SetRangeUser(r_ymin, r_ymax);
          hsys.Draw('E2,same')

      ## Lines
      if (not r_xmin) and (not r_xmax):
          r_xmin = xmin
          r_xmax = xmax
      line = TLine(r_xmin, 1, r_xmax, 1)
      line.SetLineColor(r.kGray+2);
      line.SetLineWidth(2);
      line.Draw('');

      ## Draw ratio
      for rat in ratios:
          rat.Draw('P,same');

      ## Draw frame 2 again
      pad2.Update()
      pad2.RedrawAxis()
      aux_frame2 = TLine()
      aux_frame2.SetLineWidth(2) 
      #aux_frame2.DrawLine(pad2.GetUxmax(), pad2.GetUymin(), pad2.GetUxmax(), r_ymax);
      aux_frame2.DrawLine(pad2.GetUxmax(), pad2.GetUymin(), pad2.GetUxmax(), pad2.GetUymax());

      ## Banner
      pad1.cd()
      if maxYnumbers:
          r.TGaxis().SetMaxDigits(maxYnumbers)
          self.bannerRatio(isData, lumi, scy = True, inProgress = inProgress, isPrivate = isPrivate)
      else:
          self.bannerRatio(isData, lumi, scy = False, inProgress = inProgress, isPrivate = isPrivate)


      if not outputDir[-1] == '/': dirName = outputDir + '/'
      else: dirName = outputDir

      for i,plotName in enumerate(self.plotNames):
          pad1.cd()
          pad1.SetLogy(0)
          path    = dirName+plotName
          pathlog = dirName+self.plotNamesLog[i]
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)
          if not '.root' in pathlog:
              if 'TEfficiency' not in str(type(self.histos[0])) and self.histos[0].GetMinimum() == 0:
                  self.histos[0].SetMinimum(0.1) ### log y axis consistency
              #self.histos[0].SetMaximum(100.0*self.histos[0].GetMaximum()) ### log y axis consistency
              pad1.cd()
              pad1.SetLogy()
              self.myCanvas.SaveAs(pathlog)

          

      pad1.IsA().Destructor(pad1) 
      pad2.IsA().Destructor(pad2) 
      self.myLegend.IsA().Destructor(self.myLegend)
      self.myCanvas.IsA().Destructor(self.myCanvas)                                                                                                                                            

   def saveRatio2(self, legend, isData, log, lumi, hdata, hMC, r_ymin=0, r_ymax=2, r_xmin=0, r_xmax=0, label ="Data/Prediction", sys = None, outputDir = 'plots/', xlog = False, maxYnumbers = False, inProgress = False, isPrivate = False):

      self.myCanvas.cd()

      ## Create tha pads
      pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0) # for the plot
      pad1.SetBottomMargin(0.015)
      pad1.SetTopMargin(0.13)
      pad1.Draw()                                     
      pad2 = TPad("pad2", "pad2", 0, 0.01, 1, 0.3) # for the ratio
      pad2.SetTopMargin(0.05);
      pad2.SetBottomMargin(0.4);
      #pad2.SetGridy(1);
      pad2.Draw();                                      

      ## Set log scale
      pad1.cd()
      if(log):
          pad1.SetLogy(1)
      if(xlog):
          pad1.SetLogx(1)
          pad2.SetLogx(1)

      ## Draw histograms
      maxYAxisValue = 0.0
      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):
              if str(type(self.histos[i])) == "<class 'ROOT.TEfficiency'>":
                  self.makeRate(self.histos[i], self.options[i], ratio = True)
              else:
                  self.histos[i].GetYaxis().SetTitleSize(0.045)
                  self.histos[i].GetYaxis().SetTitleOffset(1.25);
                  self.histos[i].GetXaxis().SetLabelSize(0)
                  self.histos[i].GetYaxis().SetLabelSize(0.04)
                  self.histos[i].Draw(self.options[i])
                  if self.histos[i].GetMaximum() > maxYAxisValue:
                      maxYAxisValue = self.histos[i].GetMaximum()

      ## Draw frame 1 again
      pad1.Update()
      pad1.RedrawAxis()
      aux_frame = TLine()
      aux_frame.SetLineWidth(2) 
      #aux_frame.DrawLine(pad1.GetUxmax(), pad1.GetUymin(), pad1.GetUxmax(), maxYAxisValue);
      aux_frame.DrawLine(pad1.GetUxmax(), pad1.GetUymin(), pad1.GetUxmax(), pad1.GetUymax());


      if(legend):
          self.makeLegend()
          self.myLegend.SetTextSize(0.04) # Modify the legend size
          self.myLegend.Draw()

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextColor(latex[-1])
          lat.SetTextAlign(latex[-2])
          lat.SetTextSize(latex[-3])
          lat.SetTextFont(latex[-4])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      
      ## Ratios
      if type(hMC) != list:
          hMClist = [hMC]
      else: hMClist = hMC

      ratios = []

      for tmp_hMC in hMClist:
          ind = hMClist.index(tmp_hMC)

          if str(type(tmp_hMC)) == "<class 'ROOT.TEfficiency'>":
              tmp_den = tmp_hMC.GetTotalHistogram().Clone()
              tmp_num = hdata.GetTotalHistogram().Clone()
              for n in range(0,tmp_num.GetNbinsX()):
                  tmp_den.SetBinContent(n+1, tmp_hMC.GetEfficiency(n+1))
                  tmp_num.SetBinContent(n+1, hdata.GetEfficiency(n+1))
                  tmp_den.SetBinError(n+1, tmp_hMC.GetEfficiencyErrorLow(n+1))
                  tmp_num.SetBinError(n+1, hdata.GetEfficiencyErrorLow(n+1))
              tmp_ratio = tmp_num.Clone(tmp_hMC.GetName()+'_ratio')
              tmp_ratio.Divide(tmp_den)
          else:
              tmp_ratio = hdata.Clone(tmp_hMC.GetName()+'_ratio')
              tmp_ratio.Divide(tmp_hMC)
              for n in range(1, tmp_ratio.GetNbinsX()+1):
                  if tmp_hMC.GetBinContent(n) == 0 or hdata.GetBinContent(n) == 0: continue
                  tmp_ratio.SetBinContent(n, hdata.GetBinContent(n)/tmp_hMC.GetBinContent(n))
                  tmp_ratio.SetBinError(n, hdata.GetBinError(n)/hdata.GetBinContent(n) * hdata.GetBinContent(n)/tmp_hMC.GetBinContent(n))
                       

          ## Ratio tunning
          tmp_ratio.SetTitle("")
          tmp_ratio.GetYaxis().SetRangeUser(r_ymin, r_ymax);
          tmp_ratio.GetYaxis().SetTitle(label);
          tmp_ratio.GetYaxis().CenterTitle();
          tmp_ratio.GetYaxis().SetLabelSize(0.10);
          tmp_ratio.GetYaxis().SetNdivisions(4);
          tmp_ratio.GetYaxis().SetTitleOffset(0.5);
          tmp_ratio.GetXaxis().SetLabelSize(0.10);
          tmp_ratio.GetYaxis().SetTitleSize(0.11);
          tmp_ratio.GetXaxis().SetTitleSize(0.12);
          tmp_ratio.GetXaxis().SetLabelOffset(0.02);
          if 'TEfficiency' not in str(type(self.histos[0])): 
              tmp_ratio.GetXaxis().SetTitle(self.histos[0].GetXaxis().GetTitle());
          else:
              tmp_ratio.GetXaxis().SetTitle(self.histos[0].GetTotalHistogram().GetXaxis().GetTitle());
          tmp_ratio.SetMarkerStyle(20);
          tmp_ratio.SetMarkerColor(r.kBlack);
          tmp_ratio.SetMarkerSize(0.8);
          tmp_ratio.SetMarkerColor(r.kBlack if len(hMClist) == 1 else tmp_hMC.GetMarkerColor());
          tmp_ratio.SetLineColor  (r.kBlack if len(hMClist) == 1 else tmp_hMC.GetLineColor  ());
          tmp_ratio.SetLineColor(r.kBlack);
          tmp_ratio.SetLineWidth(2);
          tmp_ratio.SetLineStyle(tmp_hMC.GetLineStyle())

          ratios.append(tmp_ratio)
          xmin = tmp_ratio.GetBinLowEdge(1)
          xmax = tmp_ratio.GetBinLowEdge(tmp_ratio.GetNbinsX()+1)

      pad2.cd();  

      ## Draw systematics (if included)
      ratios[0].Draw('AXIS')
      if sys is not None:
          hsys = tmp_ratio.Clone("sys_up")
          hsys.Reset()
          for n in range(1, tmp_hMC.GetNbinsX() + 1):
              if tmp_hMC.GetBinContent(n) == 0: continue
              hsys.SetBinContent(n, 1.0 )
              hsys.SetBinError(n, (math.sqrt(tmp_hMC.GetBinError(n)*tmp_hMC.GetBinError(n) + sys*tmp_hMC.GetBinContent(n)*sys*tmp_hMC.GetBinContent(n) ))/ tmp_hMC.GetBinContent(n))
          hsys.GetYaxis().SetTitle(label);
          hsys.GetYaxis().CenterTitle();
          hsys.GetYaxis().SetLabelSize(0.10);
          hsys.GetYaxis().SetNdivisions(4);
          hsys.GetYaxis().SetTitleOffset(0.5);
          hsys.GetXaxis().SetLabelSize(0.10);
          hsys.GetYaxis().SetTitleSize(0.11);
          hsys.GetXaxis().SetTitleSize(0.12);
          hsys.GetXaxis().SetLabelOffset(0.02);
          hsys.SetFillStyle(3244)
          hsys.SetFillColor(r.kCyan+2)
          hsys.SetMarkerSize(0)
          #hsys.SetFillStyle(3013)
          hsys.GetYaxis().SetRangeUser(r_ymin, r_ymax);
          hsys.Draw('E2,same')

      ## Lines
      if (not r_xmin) and (not r_xmax):
          r_xmin = xmin
          r_xmax = xmax
      line = TLine(r_xmin, 1, r_xmax, 1)
      line.SetLineColor(r.kGray+2);
      line.SetLineWidth(2);
      line.Draw('');

      ## Draw ratio
      for rat in ratios:
          rat.Draw('P,same');

      ## Draw frame 2 again
      pad2.Update()
      pad2.RedrawAxis()
      aux_frame2 = TLine()
      aux_frame2.SetLineWidth(2) 
      #aux_frame2.DrawLine(pad2.GetUxmax(), pad2.GetUymin(), pad2.GetUxmax(), r_ymax);
      aux_frame2.DrawLine(pad2.GetUxmax(), pad2.GetUymin(), pad2.GetUxmax(), pad2.GetUymax());

      ## Banner
      pad1.cd()
      if maxYnumbers:
          r.TGaxis().SetMaxDigits(maxYnumbers)
          self.bannerRatio(isData, lumi, scy = True, inProgress = inProgress, isPrivate = isPrivate)
      else:
          self.bannerRatio(isData, lumi, scy = False, inProgress = inProgress, isPrivate = isPrivate)


      if not outputDir[-1] == '/': dirName = outputDir + '/'
      else: dirName = outputDir

      for i,plotName in enumerate(self.plotNames):
          pad1.cd()
          pad1.SetLogy(0)
          path    = dirName+plotName
          pathlog = dirName+self.plotNamesLog[i]
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)
          if not '.root' in pathlog:
              if 'TEfficiency' not in str(type(self.histos[0])) and self.histos[0].GetMinimum() == 0:
                  self.histos[0].SetMinimum(0.1) ### log y axis consistency
              #self.histos[0].SetMaximum(100.0*self.histos[0].GetMaximum()) ### log y axis consistency
              pad1.cd()
              pad1.SetLogy()
              self.myCanvas.SaveAs(pathlog)

          

      pad1.IsA().Destructor(pad1) 
      pad2.IsA().Destructor(pad2) 
      self.myLegend.IsA().Destructor(self.myLegend)
      self.myCanvas.IsA().Destructor(self.myCanvas)                                                                                                                                            

   def save(self, legend, isData, log, lumi, labelx, ymin=0, ymax=0, outputDir = 'plots/', xlog = False, zlog = False, maxYnumbers = False, inProgress = False, is2d = False, labelz = False, isPrivate = False):

      self.myCanvas.cd()

      if(log):
          self.myCanvas.GetPad(0).SetLogy(1)
      if(xlog):
          self.myCanvas.GetPad(0).SetLogx(1)
      if(zlog):
          self.myCanvas.GetPad(0).SetLogz(1)
     
      for i in range(0, len(self.histos)):
          if(self.ToDraw[i] != 0):        
              if 'TEfficiency' in str(type(self.histos[i])):
                  if ymax:
                      self.makeRate(self.histos[i], self.options[i], False, ymin, ymax)                   
                  else:
                      self.makeRate(self.histos[i], self.options[i], False)                   
              else:
                  if ymax:
                      self.histos[i].GetYaxis().SetRangeUser(ymin, ymax)
                  self.histos[i].Draw(self.options[i])

      ## Draw axis:
      #self.histos[0].Draw('same axis')

      for band in self.bands:
          band.Draw('f')
  
      for line in self.lines:
          line.Draw()
  
      for ellipse in self.ellipses:
          ellipse.Draw()

      for arrow in self.arrows:
          arrow.Draw()
  
      for latex in self.latexs:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextColor(latex[-1])
          lat.SetTextAlign(latex[-2])
          lat.SetTextSize(latex[-3])
          lat.SetTextFont(latex[-4])
          lat.DrawLatex(latex[0], latex[1], latex[2])
  
      if labelz:
          lat = TLatex()
          lat.SetNDC()
          lat.SetTextColor(r.kBlack)
          lat.SetTextAlign(31)
          lat.SetTextSize(0.05)
          lat.SetTextFont(42)
          lat.SetTextAngle(90)
          lat.DrawLatex(0.99, 0.92, labelz)

      if(legend):
          self.makeLegend()
          self.myLegend.Draw()

      lat = TLatex()
      lat.SetNDC()
      lat.SetTextSize(0.05)
      lat.SetTextFont(42)
      lat.DrawLatex(0.46, 0.04, labelx)
      
      r.gPad.RedrawAxis()

      if not is2d:
          if maxYnumbers:
              r.TGaxis().SetMaxDigits(maxYnumbers) 
              self.bannerInFrame(isData, lumi, inProgress = inProgress, isPrivate = isPrivate)
          else:
              self.bannerInFrame(isData, lumi, inProgress = inProgress, isPrivate = isPrivate)
      else:
          if maxYnumbers:
              r.TGaxis().SetMaxDigits(maxYnumbers) 
              self.banner(isData, lumi, scy = True, inProgress = inProgress, isPrivate = isPrivate)
          else:
              self.banner(isData, lumi, scy = False, inProgress = inProgress, isPrivate = isPrivate)

      if not outputDir[-1] == '/': dirName = outputDir + '/'
      else: dirName = outputDir


      for plotName in self.plotNames:
          path = dirName+plotName
          self.ensurePath(path)
          self.myCanvas.SaveAs(path)

      #for _h in self.histos: del _h
      self.myLegend.IsA().Destructor(self.myLegend)
      self.myCanvas.IsA().Destructor(self.myCanvas)

   def savePie(self, legend, lumi, labelx):

      cpie = TCanvas("cpie","TPie test",700,700)
      
      pad1 = TPad("pad1", "pad1", 0.1, 0.1, 0.75, 0.75)
      pad1.SetBottomMargin(0.12)
      pad1.Draw()   
      
      pad1.cd()
    
      colors = []
      names = []
      vals = []
      
      for i in range(0, len(self.histos)):
          vals.append(self.histos[i].Integral())
          colors.append(self.histos[i].GetLineColor())
          names.append(self.histos[i].GetName())
      v = array('d', vals)
      c = array('i', colors)
      pie4 = TPie("p4","",len(v),v,c);
      
      pie4.SetRadius(.45);
      pie4.SetLabelsOffset(.01);
      pie4.SetLabelsOffset(100);
      pie4.SetEntryLineWidth(1,2);
      pie4.SetEntryLineWidth(2,2);
      pie4.SetEntryLineWidth(3,2);
      pie4.SetEntryLineWidth(4,2);
      pie4.SetEntryLineWidth(5,2);
      pie4.SetEntryLineWidth(6,2);
      pie4.SetEntryLineWidth(7,2);
      pie4.Draw();                                      

      self.makeLegend()
      self.myLegend.Draw()              

      lat = TLatex()
      lat.SetNDC()
      lat.SetTextSize(0.06)
      lat.SetTextFont(42)
      lat.DrawLatex(0.12, -0.1, "Slepton signal region, "+labelx)

      self.banner3(True, lumi)
      for plotName in self.plotNames:
          path = 'plots/'+plotName
          #self.ensurePath(path)
          cpie.SaveAs(path)

      self.myLegend.IsA().Destructor(self.myLegend)
      #cpie.Destructor(self.myCanvas)

 #del self.myCanvas                                                                                      
