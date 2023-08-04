import ROOT as r
import numpy as np

###
### --- color arrays definition
###
gcolors = {}
gcolors['red']       = r.TColor.GetColor('#ea4335');
gcolors['orange']    = r.TColor.GetColor('#ff6d01');
gcolors['yellow']    = r.TColor.GetColor('#fbbc04');
gcolors['green']     = r.TColor.GetColor('#34a853');
gcolors['blue']      = r.TColor.GetColor('#4285f4');
gcolors['violet']    = r.TColor.GetColor('#8f00ff');
gcolors['magenta']    = r.TColor.GetColor('#ff00ff');

acolors = {}
acolors['1']    = r.TColor.GetColor('#219ebc');
acolors['2']    = r.TColor.GetColor('#023047');
acolors['3']    = r.TColor.GetColor('#ffb703');
acolors['4']    = r.TColor.GetColor('#fb8500');

bcolors = {}
bcolors['1']    = r.TColor.GetColor('#177e89');
bcolors['2']    = r.TColor.GetColor('#084c61');
bcolors['3']    = r.TColor.GetColor('#db3a34');
bcolors['4']    = r.TColor.GetColor('#ffc857');

dcolors = {}
dcolors['1mm'] = r.kAzure+10
dcolors['10mm'] = r.kAzure-4
dcolors['100mm'] = r.kBlue-3
dcolors['1000mm'] = r.kMagenta-3
dcolors['10000mm'] = r.kMagenta+2


###
### --- 2D plot palette definition
###

# Signal palettes

sigc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),
                                            np.array([1.00, 204./255.]),
                                            np.array([1.00, 0.00]),
                                            np.array([1.00, 0.00]), 255);
sigpalette_ = []
for i in range(0, 255):
    sigpalette_.append(sigc + i)
sigpalette = np.array(sigpalette_, dtype=np.int32)


# Background (predicted) palette

bkgc = r.TColor.CreateGradientColorTable(2, np.array([0.00, 1.00]),
                                            np.array([1.00, 0.00]),
                                            np.array([1.00, 153/255.]),
                                            np.array([1.00, 153./255.]), 255);

bkgpalette_ = []
for i in range(0, 255):
    sigpalette_.append(sigc + i)

bkgpalette = np.array(bkgpalette_, dtype=np.int32)


