import ROOT as r
from array import array
from ROOT import TTree, TFile, TCut, TH1F, TH2F, TH3F, THStack, TCanvas, SetOwnership
import copy
import os
import __main__
import json
import glob

########################################################################################
########################################################################################
################################### Sample Class #######################################
########################################################################################
########################################################################################
class Sample:
   'Common base class for all Samples'
   def __init__(self, name, label, color, location, xsection, isdata):

      self.name = name
      self.label = label
      self.color = eval(color)
      self.location = location if location[-1] == '/' else location + '/'
      self.xSection = xsection
      self.isData = isdata
      self.ftpaths = []
      self.ftfiles = []
      self.ttrees = []
      self.tchain = r.TChain('Events')
      self.count = 0.0

      ## Load samples
      if '*' in location:
        for file in glob.glob(location):
          if '.root' not in file: continue
          self.ftpaths.append(file)
          self.tchain.Add(file)
      else:
        for subdir, dirs, files in os.walk(self.location):
          for file in files:
            if '.root' not in file: continue
            self.ftpaths.append(os.path.join(subdir, file))
            self.tchain.Add(os.path.join(subdir, file))
      self.rdf = r.RDataFrame(self.tchain)

      ## Set weight structure
      if not self.isData:
        self.count = self.rdf.Sum('genWeight').GetValue() 
      else:
        self.count = self.tchain.GetEntries()

      if not self.isData:
        self.rdf = self.rdf.Define('lumWeight', str(self.xSection / self.count))
        self.rdf = self.rdf.Define('mcWeight', 'lumWeight * genWeight')

      print("Loaded sample " + name + ' with a total of ' + str(self.count) + ' entries')

   def printSample(self):
      print("#################################")
      print("Sample Name: ", self.name)
      print("Sample Location: ", self.location)
      print("Sample XSection: ", self.xSection)
      print("Sample IsData: ", self.isData)
      print("Sample LumWeight: ", self.lumWeight)
      print("#################################")


   def closeFiles(self):
       for _file in self.ftfiles:
              _file.Close()


   def getTH1F(self, lumi, name, var, nbin = 0, xmin = 0.0, xmax = 0.0, cut = 'true', options = '', xlabel = '', xaxis = []):

      
      ## revisit:
      #if(self.isData == 0):
      #   cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight/abs(genWeight) " + " )" 

      ## revisit: Need to redefine weight structure
      redRFP = self.rdf.Filter(cut)
      if not self.isData:
          redRFP = redRFP.Define('weight', '(' + str(lumi) + ' * mcWeight )')
      else:
          redRFP = redRFP.Define('weight', '1.0')

      if len(xaxis) > 0:
          h = redRFP.Histo1D((name + 'noOF', ";" + xlabel + "; Counts / bin width", len(xaxis) - 1, xaxis), var, 'weight')
      else:
          h = redRFP.Histo1D((name + 'noOF', ";" + xlabel + "; Counts / bin width", nbin, xmin, xmax), var, 'weight')
      h.Sumw2()

      ## revisit: OF histogram
      ofBin = False
      h_of = h.Clone()
      h_of.Sumw2()
      for _bin in range(1, h.GetNbinsX()+2):
          h_of.SetBinContent(_bin, h.GetBinContent(_bin))
          h_of.SetBinError  (_bin, h.GetBinError  (_bin))

      return (h_of if ofBin else h)

    
   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
   
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight/abs(genWeight) )" 
     self.ttree.Project(name, var, cut, options) 
     return h

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel):
   
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx), len(nbiny)-1, array('d', nbiny), len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     if(self.isData == 0):
        cut = cut + "* ( " + str(self.lumWeight*lumi) + " * genWeight/abs(genWeight) )"
     self.ttree.Project(name, var, cut, options) 
     return h


########################################################################################
########################################################################################
#################################### Block Class #######################################
########################################################################################
########################################################################################
class Block:
   'Common base class for all Sample Blocks'

   def __init__(self, name, label, color, isdata):
      self.name  = name
      self.color = color
      self.isData = isdata
      self.label = label
      self.samples = []

   def printBlock(self):

      print("####################")
      print("Block Name: ", self.name)
      print("Block Color: ", self.color)
      print("Block IsData: ", self.isData)
      print("####################")
      print("This block contains the following Samples")

      for l in self.samples:
        l.printSample()
     

   def addSample(self, s):
      self.samples.append(s)

   def getTH1F(self, lumi, name, var, nbin = 0, xmin = 0.0, xmax = 0.0, cut = 'true', options = '', xlabel = '', xaxis = []):
     for _is,s in enumerate(self.samples):
       
       AuxName = "auxT1_sample" + s.name
       haux = s.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, xaxis)
       if not _is:
          h = haux.Clone(name+'_blockHisto')
       else:
          haux_ = haux.Clone(name+'_toAdd')
          h.Add(haux_)
       del haux

     h.SetLineColor(self.color)
     h.SetMarkerColor(self.color)
     h.GetYaxis().SetTitle('Events')
     h.SetTitle(self.label)

     return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)

     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for s in self.samples:
     
       AuxName = "auxT2_block" + s.name
       haux = s.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel)
       h.Add(haux)
       del haux

     return h   

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel):
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny),len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)

     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     for s in self.samples:
     
       AuxName = "auxT3_block" + s.name
       haux = s.getTH3F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel)
       h.Add(haux)
       del haux

     return h   

########################################################################################
########################################################################################
#################################### Tree Class ########################################
########################################################################################
########################################################################################
class Tree:
   'Common base class for a physics meaningful tree'

   def __init__(self, fileName, name, isdata, close = False):
      #print fileName
      self.name  = name
      self.isData = isdata
      self.blocks = []
      self.parseFileName(fileName, close)

   def parseFileName(self, fileName, close = False):
      f = open(fileName)

      for l in f.readlines():
        if (l[0] == "#" or len(l) < 2):
          continue

        splitedLine = str.split(l)
        block       = splitedLine[0]
        theColor    = splitedLine[1]
        name        = splitedLine[2]
        label       = splitedLine[3]
        flocation   = splitedLine[4]
        xsection    = float(splitedLine[5])
        isdata      = int(splitedLine[6])

        color = 0
        plusposition = theColor.find("+")
        if(plusposition == -1):
          color = eval(theColor)
        else:
          color = eval(theColor[0:plusposition])
          color = color + int(theColor[plusposition+1:len(theColor)])

        sample = Sample(name, label, theColor, flocation, xsection, isdata)
        coincidentBlock = [l for l in self.blocks if l.name == block]
        if(coincidentBlock == []):
          newBlock = Block(block, label, color, isdata)
          newBlock.addSample(sample)
          self.addBlock(newBlock)
        else:
          coincidentBlock[0].addSample(sample)

        if close:
            sample.closeFiles()



   def printTree(self):
      print("######")
      print("Tree Name: ", self.name)
      print("Tree IsData: ", self.isData)
      print("######")
      print("This Tree contains the following Blocks")

      for l in self.blocks:
        l.printBlock()
     

   def addBlock(self, b):
      self.blocks.append(b)


   def loadUtils(self, library):
      r.gROOT.LoadMacro(library)


   def setDefinitions(self, config = False):

      with open(config) as f:
            configuration = json.load(f)
      for d in configuration['definitions']:
        for b in self.blocks:
          for s in b.samples:
            s.rdf = s.rdf.Define(d[0], d[1])


   def setSelection(self, selection):

       for b in self.blocks:
          for s in b.samples:
            s.rdf = s.rdf.Filter(selection)


   def getYields(self, lumi, var, xmin, xmax, cut):
      h = self.getTH1F(lumi, "yields", var, 1, xmin, xmax, cut, "", "")
      nbinmin = h.FindBin(xmin)
      nbinmax = h.FindBin(xmax)
      error = r.Double()
      value = h.IntegralAndError(nbinmin, nbinmax, error)
      y = [value, error]
      
      del h
      return y


   def getStack(self, lumi, name, var, nbin, xmin, xmax, cut, options, xlabel):
     if cut == '':
       cut = '(1)'
     hs = THStack(name, "")
     SetOwnership(hs, 0 )
     for b in self.blocks:
     
       AuxName = "auxStack_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel)
       haux.SetFillColor(b.color)
       haux.SetLineColor(r.kBlack)
       haux.SetTitle(b.label)
       hs.Add(haux)
       del haux


     can_aux = TCanvas("can_%s_%s"%(name, b.name))
     can_aux.cd()
     hs.Draw()

     del can_aux

     ylabel = "Events"
     if xmax != xmin:
       hs.GetXaxis().SetTitle(xlabel)
       b = int((xmax-xmin)/nbin)
       ylabel = "Events / " + str(b) + " GeV"
     else:     
       ylabel = "Events"
   
     hs.GetYaxis().SetTitle(ylabel)
     return hs   


   def getTH1F(self, lumi, name, var, nbin = 0, xmin = 0.0, xmax = 0.0, cut = '', options = '', xlabel = '', xaxis = []):
     if cut == '':
       cut = '(1)'
     for ib,b in enumerate(self.blocks):
       AuxName = "auxh1_block_" + name + "_" + b.name
       haux = b.getTH1F(lumi, AuxName, var, nbin, xmin, xmax, cut, options, xlabel, xaxis)
       if not ib:
          h = haux.Clone(name+'_treeHisto')
       else:
          h.Add(haux)
       del haux
       
     return h

   def getTH2F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel):
     if cut == '':
       cut = '(1)'
     if(xmin == xmax) and (ymax == ymin):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny))
     elif (xmin == xmax):
        h = TH2F(name, "", len(nbinx)-1, array('d', nbinx),nbiny,ymin,ymax)
     elif (ymin == ymax):
        h = TH2F(name, "", nbinx,xmin,xmax,len(nbiny)-1, array('d', nbiny))
     else: 
        h = TH2F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
        
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     
     for b in self.blocks:
     
       AuxName = "aux_block" + name + "_" + b.name
       haux = b.getTH2F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, cut, options, xlabel, ylabel)
       h.Add(haux)
       del haux

     return h   

   def getTH3F(self, lumi, name, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel):
     if cut == '':
       cut = '(1)'
     if(xmin == xmax) and (ymax == ymin) and (zmax == zmin):
        h = TH3F(name, "", len(nbinx)-1, array('d', nbinx),len(nbiny)-1, array('d', nbiny), len(nbinz)-1, array('d', nbinz))
     else: 
        h = TH3F(name, "", nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax)
        
     h.Sumw2()
     h.GetXaxis().SetTitle(xlabel)
     h.GetYaxis().SetTitle(ylabel)
     h.GetZaxis().SetTitle(zlabel)
     
     for b in self.blocks:
     
       AuxName = "aux_block" + name + "_" + b.name
       haux = b.getTH3F(lumi, AuxName, var, nbinx, xmin, xmax, nbiny, ymin, ymax, nbinz, zmin, zmax, cut, options, xlabel, ylabel, zlabel)
       h.Add(haux)
       del haux

     return h   



