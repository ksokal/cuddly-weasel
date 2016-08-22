import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import Moog960
#import Moog960_origish as Moog960
import SpectralTools
import astropy.io.fits as pyfits
import random
import string
import sys

def parseParams(labels):
    keys = labels[0].parameters.keys()
    params = {}
    for key in keys:
        params[key] = []
        
    for label in labels:
        for key in keys:
            params[key].append(label.parameters[key])
            
    for key in keys:
        params[key] = numpy.unique(params[key])
        
    return params


def computeMs(orchestra):
    #compute the Interaction and Control Matrix for each grid point

    #should be a list of the labels of all of the spectra
    convolved = orchestra.convolved_labels
    print 'CONVOLVED', convolved
    
    #in music, mastering is taking all of the different inputs and recording them onto a single take
    #this takes the 4 wavelength windows and turns it into 1

    mastered = orchestra.master(selectedLabels=convolved, keySignature='CONVOLVED')
    #looks at the labels, convolved contains a list of all of the labels of the melodies (all the Teff, surface gravities, W start, etc)
    #this is a dictionary
    
    print 'MASTERED', mastered
    
    params = parseParams(mastered)
    print 'PARAMS', params


    wlStart = numpy.min(params["WLSTART"])
    wlStop = numpy.max(params["WLSTOP"])

    for blah in mastered:
        blah.parameters["SELECTED"] = True

    #orchestra.selectExcerpt(wlRange =[wlStart, wlStop], exact=True)

    TeffStroke = 50.0
    loggStroke = 0.5
    BfieldStroke = 1.0
    Tmin = numpy.min(params["TEFF"])
    Tmax = numpy.max(params["TEFF"])
    Gmin = numpy.min(params["LOGG"])
    Gmax = numpy.max(params["LOGG"])
    Bmin = numpy.min(params["BFIELD"])
    Bmax = numpy.max(params["BFIELD"])

    print 'So we are actually ready to     make the matrices!'
    for T in params["TEFF"]:
        for G in params["LOGG"]:
            for B in params["BFIELD"]:
                IM = []
                factor = []
                plusT = numpy.min([T+TeffStroke, Tmax])
                plusMelody, plusLabels = orchestra.blend(desiredParameters = {"TEFF":plusT, 
                    "LOGG":G, "BFIELD":B})
                plus, plusLabel = plusMelody.perform(label=plusLabels[0])

                minusT = numpy.max([T-TeffStroke, Tmin])
                minusMelody, minusLabels = orchestra.blend(desiredParameters = {"TEFF":minusT, 
                    "LOGG":G, "BFIELD":B})
                minus, minusLabel = minusMelody.perform(label=minusLabels[0])

                IM.append(plus-minus)
                factor.append(plusLabel.parameters["TEFF"] - minusLabel.parameters["TEFF"])

                print "lets start with changing the Temp"           
                print plusT, minusT
                ax1.clear()
                plusLabel.Spectrum.plot(ax=ax1)
                minusLabel.Spectrum.plot(ax=ax1)
                ax1.figure.show()
                #raw_input()

                plusG = numpy.min([G+loggStroke, Gmax])
                plusMelody, plusLabels = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":plusG, "BFIELD":B})
                plus, plusLabel = plusMelody.perform(label=plusLabels[0])

                minusG = numpy.max([G-loggStroke, Gmin])
                minusMelody, minusLabels = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":minusG, "BFIELD":B})
                minus, minusLabel = minusMelody.perform(label=minusLabels[0])

                IM.append(plus-minus)
                factor.append(plusLabel.parameters["LOGG"] - minusLabel.parameters["LOGG"])

                print "Now lets try log g"
                print plusG, minusG
                ax1.clear()
                plusLabel.Spectrum.plot(ax=ax1)
                minusLabel.Spectrum.plot(ax=ax1)
                ax1.figure.show()
                raw_input() #this is what you remove so you don't have to enter every time

                plusB = numpy.min([B+BfieldStroke, Bmax])
                plusMelody, plusLabels = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":G, "BFIELD":plusB})
                plus, plusLabel = plusMelody.perform(label=plusLabels[0])

                minusB = numpy.max([B-BfieldStroke, Bmin])
                minusMelody, minusLabels = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":G, "BFIELD":minusB})
                minus, minusLabel = minusMelody.perform(label=minusLabels[0])

                IM.append(plus-minus)
                factor.append(plusLabel.parameters["BFIELD"] - minusLabel.parameters["BFIELD"])

                print "Finally, for the Bfield"
                print plusB, minusB
                ax1.clear()
                plusLabel.Spectrum.plot(ax=ax1)
                minusLabel.Spectrum.plot(ax=ax1)
                ax1.figure.show()
                #raw_input()

    print blah
    return


class ControlMatrix( object ):
    def __init__(self, IM, factor, observed, Gain, ax=None):
        self.IM = IM
        self.factor = factor
        self.observed = observed
        self.ax = ax
        self.calcCM()
        self.Gain = Gain

    def calcCM(self):
        lineIMs = self.IM[:-4:,]
        dims = lineIMs.shape
        U,S,V = scipy.linalg.svd(lineIMs)
        if self.ax != None:
            self.ax.clear()
            self.ax.plot(numpy.log10(S))
            self.ax.figure.show()
            blah = raw_input("Enter number of modes to keep: ")
            try:
                nFiltModes = dims[0] - int(blah)
                self.nFiltModes = nFiltModes
            except:
                pass
            self.ax.clear()
        if self.nFiltModes == 0:
            D = 1.0/S
        else:
            D = 1.0/(S[0:-self.nFiltModes])
            S[-self.nFiltModes:] = 0.0
        newS = numpy.zeros((dims[0], dims[1]))
        I = [i for i in range(dims[1])]
        for i in range(len(D)):
            newS[i][i] = D[i]

        S = newS.copy()
        lineCM = numpy.array(scipy.matrix(V.T.dot(S.T.dot(U.T)),dtype=numpy.float32)).T

        globalIMs = self.IM[-4:,]
        dims = globalIMs.shape
        U,S,V = scipy.linalg.svd(globalIMs)
        D = 1.0/(S[0:-1])
        S[-1:] = 0.0
        newS = numpy.zeros((dims[0], dims[1]))
        I = [i for i in range(dims[1])]
        for i in range(len(D)):
            newS[i][i] = D[i]

        S = newS.copy()
        globalCM = numpy.array(scipy.matrix(V.T.dot(S.T.dot(U.T)),dtype=numpy.float32)).T

        self.CM = numpy.zeros((lineCM.shape[0]+globalCM.shape[0], globalCM.shape[1]))
        for i in range(lineCM.shape[0]):
            self.CM[i] = lineCM[i]
        for j in range(globalCM.shape[0]):
            self.CM[i+j+1] = globalCM[j]
        

    def dot(self, difference):
        overlap_start = numpy.max([numpy.min(self.observed.wl), numpy.min(difference.wl)])
        overlap_stop = numpy.min([numpy.max(self.observed.wl), numpy.max(difference.wl)])
        overlap = scipy.where((self.observed.wl >= overlap_start) & (self.observed.wl <= overlap_stop))

        diff = numpy.zeros(len(self.observed.wl))

        diff[overlap] = difference.flux_I

        return self.CM.dot(diff)

    def getCommand(self, synthesized, ax=None):
        difference = self.observed - synthesized
        #tot = []
        #for i in range(len(difference.wl)):
        #    diff = difference.flux_I[i]
        #    cm = self.CM[3][i]
        #    tot.append(diff*cm)
        #    print i, diff, cm, diff*cm, tot

        #diff = self.dot(difference)
        if ax!=None:
            ax.clear()
            difference.plot(ax=ax)
            ax.plot(self.observed.wl, self.CM[-2])
            ax.plot(self.observed.wl, self.IM[-2])
            ax.figure.show()
            raw_input()
        command =  self.Gain*(self.dot(difference))/self.factor
        command[factor < 1e-3] = 0.0
        return command

    
'''
This program is to take your grid of synthetic spectra then try to find a best fit model. 

i am editing this to see what works.

>run theremin_kim.py alpha.cfg 'Alpha'

     '''
     ###
     #1. read in the stuff you want to define and parse it
     ###
     
#read in alpha.cfg
configFile = sys.argv[1]
#'Alpha' or 'Bravo' or whatever instance it is, only important point is to run only one instance at a time (could run one on Alpha, the other on the other)
moogInstance = sys.argv[2]
try:
    #optional parameter, allow the continuum to float. so you aren't exactly certain of it, and this is mostly a good idea
    contFloat = bool(sys.argv[3]== 'True')
except:
    contFloat = False

    #now put that info into some use
config = AstroUtils.parse_config(configFile)


#now read in your grid of synthetic models
#NOTE THAT THIS IS THE START OF YOUR FILE NAMES AS WELL
directory_grid='/Users/sokal/Desktop/TWHya Kband Grid/TWHydra_Kband/MoogStokesSynthSpec_T5*B1*'
print "Grid from ", directory_grid


observed=config["inFile"]
print 'observed?', observed
orchestra = Moog960.Score(directory=directory_grid, observed=observed, suffix='')



#orchestra = Moog960.Score(directory=directory_grid, observed=None, suffix='')

print "Loaded into a score"


###
#2. Now lets figure out what happens when make changes to parameters, for each point in the grid
#     i - make interaction matrix 
#     ii - then invert it to a control matrix
###

#let make a plot axis so we can check that the right things are happening

f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])




#this is where it is interesting, compute the interaction matrix = what happens when you change T and other parameters?
Ms = computeMs(orchestra) # and the CMs



### so for now lets end here

print 'I AM DonE NOW'
"""
twhya = observed[0][0]
twhya.wl += wlOffset
twhya.flux_I += continuum

window = ((twhya.wl > wlRange[0]) & (twhya.wl < wlRange[1]))

twhya.wl = twhya.wl[window]
twhya.flux_I = twhya.flux_I[window]

Synth.run()
nominal = Synth.Spectra[0]
smoothed = nominal.resample(R=resolution)
resampled = nominal.resample(R=resolution, observedWl=twhya.wl, pad=1.0)
Spectra = [resampled]

if calcIM:
    calculateIM(Synth, resampled, resolution)
"""

Gain = 0.15
IM = pyfits.getdata(Synth.config["outputDir"]+Synth.config["baseName"]+"_IM.fits")
factor = pyfits.getdata(Synth.config["outputDir"]+Synth.config["baseName"]+"_factors.fits")
gfIndices = pyfits.getdata(Synth.config["outputDir"]+Synth.config["baseName"]+'_gf.fits')
vdWIndices = pyfits.getdata(Synth.config["outputDir"]+Synth.config["baseName"]+'_vdW.fits')
CM = ControlMatrix(IM, factor, solar, Gain, ax=ax1)
#CM = ControlMatrix(IM, factor, solar, nFiltModes, Gain)

hdu = pyfits.PrimaryHDU(CM.CM)
hdu.writeto(Synth.config["outputDir"]+Synth.config["baseName"]+"_CM.fits", clobber=True)

cb = ax1.matshow(CM.CM, aspect='auto')
blah = pyplot.colorbar(cb)
f1.show()

f2 = pyplot.figure(1)
f2.clear()
ax2 = f2.add_axes([0.1, 0.1, 0.8, 0.8])

f3 = pyplot.figure(2)
f3.clear()
ax3 = f3.add_axes([0.1, 0.1, 0.8, 0.8])

nLines = Synth.lineList.numLines

continuum = 0.00
wlOffset = 0.0
rotation = 0.0

for j in range(50):
    #command = CM.getCommand(Spectra[-1], ax=ax3)
    command = CM.getCommand(Spectra[-1])
    
    for i in range(len(gfIndices)):
        if gfIndices[i] != -1:
            Synth.lineList.perturbGf(i, command[gfIndices[i]], push=True)
    for i in range(len(vdWIndices)):
        if vdWIndices[i] != -1:
            Synth.lineList.perturbVdW(i, command[vdWIndices[i]], push=True)

    if contFloat:
        continuum = continuum+command[-4]
    else:
        continuum = 0.0
    wlOffset = wlOffset+command[-3]
    #resolution = resolution*(1.0+command[-2])
    if rotFloat:
        rotation = rotation+command[-1]
    else:
        rotation = 0.0
    Synth.run()
    output = Synth.Spectra[-1].resample(R=resolution).rotate(angle=rotation)
    Synth.Spectra = []
    output.wl += wlOffset
    output.flux_I += continuum
    output.bin(solar.wl)
    Spectra.append(output)


    print j, continuum, wlOffset, resolution, rotation
    power = numpy.mean((command**2.0)**(0.5))
    print Synth.config["baseName"], numpy.std(command), numpy.mean((command**2.0)**(0.5))
    #raw_input()

ax2.clear()
ax2.plot(solar.wl, solar.flux_I, lw=2.0)
ax2.plot(Spectra[0].wl, Spectra[0].flux_I)
ax2.plot(Spectra[-1].wl, Spectra[-1].flux_I)
f2.show()
f2.savefig(Synth.config["outputDir"] + Synth.config["baseName"]+'_results.png')
Synth.lineList.saveLineList(filename=outFile, changed=True)
