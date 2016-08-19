import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import Moog960
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

def computeIMs(orchestra):
    #this should maybe be renamed to compute Control Matrices
    #raw, interpolated, integrated, convolved, observed = orchestra.getMelodyParams()

    convolved = orchestra.convolved_labels
    #in music, mastering is taking all of the different inputs and recording them onto a single take
    #this takes the 4 wavelength windows and turns it into 1
    orchestra.master(selectedLabels=convolved, keySignature='CONVOLVED')

    #looks at the labels, convolved contains a list of all of the labels of the melodies (all the Teff, surface gravities, W start, etc)
    #this is a dictionary
    params = parseParams(convolved)

    wlStart = numpy.min(params["WLSTART"])
    wlStop = numpy.max(params["WLSTOP"])

    #this chooses to just take the models that have been convolved/peiced together to the entire wavelength range (the total 4 windows)
    orchestra.selectExcerpt(wlRange =[wlStart, wlStop], exact=True)
    ##so possibly instead of this, have .master output, and that will just then save the ones you want instead of searching for them

    #so we are just defining the range of all of our parameters
    TeffStroke = 50.0
    loggStroke = 0.5
    BfieldStroke = 1.0
    Tmin = numpy.min(params["TEFF"])
    Tmax = numpy.max(params["TEFF"])
    Gmin = numpy.min(params["LOGG"])
    Gmax = numpy.max(params["LOGG"])
    Bmin = numpy.min(params["BFIELD"])
    Bmax = numpy.max(params["BFIELD"])

    #now we want to compute an interaction matrix. Need to figure out how much of a signal to put in to produce a change
    #do this for each point in our grid, define a +T and -T etc. but this delta is inbetween our gridpoints, and assume that any response to this
    #change is roughly linear (as our models have infinite resolution)
    for T in params["TEFF"]:
        for G in params["LOGG"]:
            for B in params["BFIELD"]:
                #orchestra.selectMelodies(wlRange =[wlStart, wlStop], exact=True)
                IM = []
                factor = []
                plusT = numpy.min([T+TeffStroke, Tmax])
                #what blend is doing, it knows what gridpoints are available, it takes values of Veff etc to interpolate and it finds the 8 models or however many it takes to interpolate to the desired value, and does a weighted average. this creates a new merged spectrum to correspond to the input parameters (that were slightly off the grid), the final parameter it returns is the differential sensitivity to each parameter. 
                plusMelody, plusLabels = orchestra.blend(desiredParameters = {"TEFF":plusT, 
                    "LOGG":G, "BFIELD":B})
                #this extracts the spectrum refered to by this label
                plus, plusLabel = plusMelody.perform(label=plusLabels[0])

                minusT = numpy.max([T-TeffStroke, Tmin])
                minusMelody, minusLabels = orchestra.blend(desiredParameters = {"TEFF":minusT, 
                    "LOGG":G, "BFIELD":B})
                minus, minusLabel = minusMelody.perform(label=minusLabels[0])

                #so subtract the + from the - to get a total change, and add this to your matrix, it will add a 1 column of delta T responses
                IM.append(plus-minus)
                factor.append(plusLabel.parameters["TEFF"] - minusLabel.parameters["TEFF"]) #why is this there, this is just deltaT=2*Teffstroke

                #print asdf

                #now we do this to find responses for delta g and delta B

                plusG = numpy.min([G+loggStroke, Gmax])
                plusMelody, plusLabels = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":plusG, "BFIELD":B})
                plus, plusLabel = plusMelody.perform(label=plusLabels[0])

                minusG = numpy.max([G-loggStroke, Gmin])
                minusMelody, minusLabels = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":minusG, "BFIELD":B})
                minus, minusLabel = minusMelody.perform(label=minusLabels[0])

                IM.append(plus-minus)
                factor.append(plusLabel.parameters["TEFF"] - minusLabel.parameters["TEFF"])

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

                #so for each point in the grid we now have an interaction matrix with three columns, TBG
                #so we need to also do this for rk=veiling (excess / continuum flux at 2.2um), delta lambda for wavelength shift in the lines (go rvs!), and continuum shift
                #don't need to do the blend for those

                #once you have an IM that is 6 by the number of wavelength points, we do the singularly valued decomposition  of the interaction matrix
                #essententially you invert it, which make a control matrix
                print asdf

                #so now compute the control matrix for each point in the grid

                #then save this for each, so the total output will be a CM for each grid point - so if the grid is X x Y x Z, then each point in XYZ is also 6 x N_lambda
    print blah
    return

'''
This program is to take your grid of synthetic spectra then try to find a best fit model.     
     '''

#read in alpha.cfg
configFile = sys.argv[1]
#'False' and I changed the default to that
calcIM = bool(sys.argv[2]=='False')
#'Alpha' or 'Bravo' or whatever instance it is, only important point is to run only one instance at a time (could run one on Alpha, the other on the other)
moogInstance = sys.argv[3]
try:
    #optional parameter, allow the continuum to float. so you aren't exactly certain of it, and this is mostly a good idea
    contFloat = bool(sys.argv[4]== 'True')
except:
    contFloat = False
try:
    #optional parameter, this was used to measure oscillator strengths and allows you to rotate your continuum - specifically for one paper to get a saw tooth continuum. Don't use.
    rotFloat = bool(sys.argv[5] == 'True')
except:
    rotFloat = False

    
f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])

config = AstroUtils.parse_config(configFile)

orchestra = Moog960.Score(directory='../MusicMaker/TWHydra_T3*', observed=config["inFile"], suffix='')

#this is where it is interesting, compute the interaction matrix = what happens when you change T and other parameters?
IMs = computeIMs(orchestra) # and the CMs, which are going to be added here

#why? this is for oscillator strength
"""
Synth = MoogTools.MoogStokes(configFile, moogInstance=moogInstance,
            fileBase=''.join(random.choice(string.ascii_letters)
            for _ in range(3)), progressBar=True)

outFile = Synth.config["outputDir"] + Synth.config["baseName"]+'_results.dat'
inFile = Synth.config["inFile"]

continuum = 0.00
wlOffset = Synth.config["wlShift"]
resolution = Synth.config["resolvingPower"]

twhya = observed[0][0]
twhya.wl += wlOffset #that means add the value and store the result in the value  (x=x+value)
twhya.flux_I += continuum

window = ((twhya.wl > wlRange[0]) & (twhya.wl < wlRange[1]))

twhya.wl = twhya.wl[window]
twhya.flux_I = twhya.flux_I[window]

Synth.run()
nominal = Synth.Spectra[0]
smoothed = nominal.resample(R=resolution)
resampled = nominal.resample(R=resolution, observedWl=twhya.wl, pad=1.0)
Spectra = [resampled]

#if calcIM:
#    calculateIM(Synth, resampled, resolution)
"""

#so start with a basic guess of a synthetic spectra, is on the grid
#then we take the model associated with that grid point

#subtract the observed spectrum from that model, which gives us a spectrum of residuals, which is the same dir
#take the residual spectra and take the dot product with the CM, which will then result in the 6 parameters, which is the correcitons
#that you need to make to your parameters

#now you adjust your guess synthetic spectra, running it through blend with the corrections applied
#now you do this again

#this is then repeated for a few iterations, defined by Gain. gain of .1 is 10 iterations, .01 then 1000 iterations
#could input some sort of thing to say convergence to quit on if you like


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

    """
    ax2.clear()
    if diffplot:
        difference = solar.diff_spectra(Spectra[-1], pad=True)
        difference.flux_I[Spectra[-1].flux_I == 0] = 0.0
        ax2.plot(difference.wl, difference.flux_I)
    else:
        ax2.plot(solar.wl, solar.flux_I, lw=2.0)
        for spec in Spectra:
            ax2.plot(spec.wl, spec.flux_I)
    f2.show()
    #"""
    """
    ax3.clear()
    for line in range(Synth.lineList.nStrong):
        ax3.plot(Synth.lineList.strongLines[line].loggfHistory)
        ax3.plot(Synth.lineList.strongLines[line].VdWHistory)
    for line in range(Synth.lineList.numLines - Synth.lineList.nStrong):
        ax3.plot(Synth.lineList.weakLines[line].loggfHistory)
        ax3.plot(Synth.lineList.weakLines[line].VdWHistory)
    f3.show()
    #"""

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
