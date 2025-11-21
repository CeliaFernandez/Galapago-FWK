import ROOT as r
from array import array
from ROOT import TTree, TFile, TCut, TH1F, TH2F, TH3F, THStack, TCanvas, SetOwnership
import copy
import os
import __main__
import json
import glob
import subprocess
import hashlib

# XRootD redirector for CMS data access
XROOTD_REDIRECTOR = "root://cms-xrd-global.cern.ch/"

def queryDAS(dataset, limit=0):
    """
    Query CMS DAS for files in a dataset.

    Args:
        dataset: DAS dataset name (e.g., '/Muon0/Run2023C-PromptNanoAODv12-v1/NANOAOD')
        limit: Maximum number of files to return (0 = all)

    Returns:
        List of xrootd file paths
    """
    query = f"file dataset={dataset}"
    cmd = ["dasgoclient", "--query", query]
    if limit > 0:
        cmd.extend(["--limit", str(limit)])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        files = [f for f in result.stdout.strip().split('\n') if f]
        # Convert LFN to xrootd URL
        xrootd_files = [XROOTD_REDIRECTOR + f for f in files]
        print(f"DAS query returned {len(xrootd_files)} files for {dataset}")
        return xrootd_files
    except subprocess.CalledProcessError as e:
        print(f"DAS query failed: {e.stderr}")
        return []
    except FileNotFoundError:
        print("ERROR: dasgoclient not found. Please setup CMS environment (cmsenv)")
        return []


def getCachedDASFiles(dataset, cache_dir=".das_cache"):
    """
    Get DAS files with local caching to avoid repeated queries.
    """
    os.makedirs(cache_dir, exist_ok=True)
    cache_key = hashlib.md5(dataset.encode()).hexdigest()
    cache_file = os.path.join(cache_dir, f"{cache_key}.txt")

    if os.path.exists(cache_file):
        with open(cache_file, 'r') as f:
            files = [line.strip() for line in f if line.strip()]
            print(f"Using cached DAS query for {dataset} ({len(files)} files)")
            return files

    files = queryDAS(dataset)
    if files:
        with open(cache_file, 'w') as f:
            f.write('\n'.join(files))
    return files


def isDASDataset(location):
    """Check if location is a DAS dataset path."""
    return location.startswith('/') and location.count('/') >= 3 and '/store/' not in location

########################################################################################
########################################################################################
################################### Sample Class #######################################
########################################################################################
########################################################################################
class Sample:
   'Common base class for all Samples'
   def __init__(self, name, label, color, location, xsection, isdata, file_limit=0):

      self.name = name
      self.label = label
      self.color = eval(color) if isinstance(color, str) else color
      self.location = location
      self.xSection = xsection
      self.isData = isdata
      self.ftpaths = []
      self.tchain = r.TChain('Events')
      self.count = 0.0

      ## Load samples based on location type
      self._loadFiles(location, file_limit)

      if not self.ftpaths:
          raise RuntimeError(f"No ROOT files found for sample {name} at {location}")

      self.rdf = r.RDataFrame(self.tchain)

      ## Set weight structure
      if not self.isData:
        self.count = self.rdf.Sum('genWeight').GetValue()
      else:
        self.count = self.tchain.GetEntries()

      if not self.isData:
        self.rdf = self.rdf.Define('lumWeight', str(self.xSection / self.count))
        self.rdf = self.rdf.Define('mcWeight', 'lumWeight * genWeight')

      print(f"Loaded sample {name} with {len(self.ftpaths)} files and {self.count:.0f} entries")

   def _loadFiles(self, location, file_limit=0):
      """Load ROOT files from DAS dataset, xrootd, EOS, or local path."""

      # Case 1: DAS dataset (e.g., /Muon0/Run2023C-PromptNanoAODv12-v1/NANOAOD)
      if isDASDataset(location):
          files = getCachedDASFiles(location)
          if file_limit > 0:
              files = files[:file_limit]
          for f in files:
              self.ftpaths.append(f)
              self.tchain.Add(f)
          return

      # Case 2: XRootD URL or EOS path
      if location.startswith('root://') or location.startswith('/eos/'):
          if location.startswith('/eos/'):
              # Convert EOS path to xrootd URL for remote access
              location = "root://eoscms.cern.ch/" + location
          if '*' in location:
              # Glob pattern - requires local access
              for f in glob.glob(location.replace('root://eoscms.cern.ch/', '')):
                  if '.root' in f:
                      xrd_path = "root://eoscms.cern.ch/" + f if not f.startswith('root://') else f
                      self.ftpaths.append(xrd_path)
                      self.tchain.Add(xrd_path)
          else:
              self.ftpaths.append(location)
              self.tchain.Add(location)
          return

      # Case 3: Local directory or glob pattern
      location = location if location[-1] == '/' else location + '/'
      if '*' in location:
          for f in glob.glob(location):
              if '.root' in f:
                  self.ftpaths.append(f)
                  self.tchain.Add(f)
      else:
          for subdir, dirs, files in os.walk(location):
              for f in files:
                  if '.root' not in f:
                      continue
                  filepath = os.path.join(subdir, f)
                  self.ftpaths.append(filepath)
                  self.tchain.Add(filepath)
                  if file_limit > 0 and len(self.ftpaths) >= file_limit:
                      return

   def printSample(self):
      print("#################################")
      print(f"Sample Name: {self.name}")
      print(f"Sample Location: {self.location}")
      print(f"Sample XSection: {self.xSection}")
      print(f"Sample IsData: {self.isData}")
      print(f"Sample Files: {len(self.ftpaths)}")
      print(f"Sample Entries: {self.count:.0f}")
      print("#################################")


   def closeFiles(self):
       # With RDataFrame, files are managed automatically
       # This method is kept for backward compatibility
       pass


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

   def __init__(self, fileName, name, isdata, close=False, file_limit=0):
      self.name = name
      self.isData = isdata
      self.blocks = []
      self.file_limit = file_limit
      self.parseFileName(fileName, close, file_limit)

   def parseFileName(self, fileName, close=False, file_limit=0):
      with open(fileName) as f:
          lines = f.readlines()

      for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
          continue

        splitedLine = line.split()
        if len(splitedLine) < 7:
          continue

        block = splitedLine[0]
        theColor = splitedLine[1]
        name = splitedLine[2]
        label = splitedLine[3]
        flocation = splitedLine[4]
        xsection = float(splitedLine[5])
        isdata = int(splitedLine[6])

        color = 0
        plusposition = theColor.find("+")
        if plusposition == -1:
          color = eval(theColor)
        else:
          color = eval(theColor[0:plusposition])
          color = color + int(theColor[plusposition+1:])

        sample = Sample(name, label, theColor, flocation, xsection, isdata, file_limit)
        coincidentBlock = [b for b in self.blocks if b.name == block]
        if not coincidentBlock:
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



