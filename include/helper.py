import math, sys, os, copy, re
import ROOT as r
from   ROOT import TGraphErrors, gROOT, TCanvas, TFile



def ensurePath( _path):
    d = os.path.dirname(_path)
    if not os.path.exists(d):
        os.makedirs(d)             

def ensureDirectory(_path):
   if not os.path.exists(_path):
      os.makedirs(_path)

def selectSamples(inputfile, selList, sType = 'DATA'):
    f = open(inputfile, 'r')
    tmp_file = open('.tmp_sampleFile%s.txt' %sType, 'w')

    checkedList = []
    typeList    = []
    print("selectSamples for ", sType, ": List Of Samples:" , selList)


    for line in f.readlines():
        if '#'==line[0] or not len(line.rstrip('\r')): continue
        for _sample in selList:
            _sample = _sample.replace('*','.*')

            cond_std = _sample == line.split()[2]
            cond_re = re.search(_sample, line.split()[2]) != 0
            if _sample == line.split()[2]:
                print("---> Found a match for", _sample, ":",  line.split()[2], " ", line.split()[3], line.split()[4], line.split()[5], ", matchesName_=", cond_std, ", matchesRegExp", cond_re) 
                if not sType == 'SYNCH':
                    tmp_file.write(line)
                else:
                    tmp_splitline = line.split()
                    tmp_splitline[0] = 'synching'
                    tmp_file.write('  '.join(tmp_splitline+['\n']))
                checkedList.append(_sample)
                typeList   .append(int(line.split()[-1]))
    for _selSample in selList:
        if not '*' in _selSample:
            if _selSample not in checkedList:
                print('ERROR: some samples weren\'t selected, check all sample names!')
                print('check this sample:', _selSample)
                sys.exit('exiting...')
        else:
            print('you used some wildcards in selecting the samples. be careful with that!')
    if not len(set(typeList)) == 1:
            print('ERROR: you\'re mixing DATA and MC!')
            sys.exit('exiting...')
            
    return tmp_file.name



##################################################################################
### These functions were created by Aachen in order to define the proper style ###
##################################################################################
def createMyColors():
    iIndex = 2000

    containerMyColors = []
    for color in defineMyColors.keys():
       	tempColor = r.TColor(iIndex,
       	float(defineMyColors[color][0]) / 255, float(defineMyColors[color][1]) / 255, float(defineMyColors[color][2]) / 255)
       	containerMyColors.append(tempColor)

       	myColors.update({ color: iIndex })
       	iIndex += 1

    return containerMyColors


defineMyColors = {
        'Black' : (0, 0, 0),
        'White' : (255, 255, 255),
        'Red' : (255, 0, 0),
        'DarkRed' : (128, 0, 0),
        'Green' : (0, 255, 0),
        'Blue' : (0, 0, 255),
        'Yellow' : (255, 255, 0),
        'Orange' : (255, 128, 0),
        'DarkOrange' : (255, 64, 0),
        'Magenta' : (255, 0, 255),
        'KDEBlue' : (64, 137, 210),
        'Grey' : (128, 128, 128),
        'DarkGreen' : (0, 128, 0),
        'DarkSlateBlue' : (72, 61, 139),
        'Brown' : (70, 35, 10),

        'MyBlue' : (36, 72, 206),
        'MyDarkBlue' : (18, 36, 103),
        'MyGreen' : (70, 164, 60),
        'AnnBlueTitle' : (29, 47, 126),
        'AnnBlue' : (55, 100, 255),
#        'W11AnnBlue' : (0, 68, 204),
#        'W11AnnBlue' : (63, 122, 240),
    }


myColors = {
            'W11ttbar':  855,
            'W11singlet':  854,
            'W11ZLightJets':  401,
            'W11ZbJets':  400,
            'W11WJets':  842,
            'W11Diboson':  920,
            'W11AnnBlue': 856,
            'W11Rare':  630,
            }


