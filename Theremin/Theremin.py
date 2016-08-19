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
    #raw, interpolated, integrated, convolved, observed = orchestra.getMelodyParams()
    convolved = orchestra.convolved_labels

    mastered = orchestra.master(selectedLabels=convolved, keySignature='CONVOLVED')

    params = parseParams(mastered)

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

           
                print plusT, minusT
                ax1.clear()
                plusLabel.Spectrum.plot(ax=ax1)
                minusLabel.Spectrum.plot(ax=ax1)
                ax1.figure.show()
                raw_input()

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

                print plusG, minusG
                ax1.clear()
                plusLabel.Spectrum.plot(ax=ax1)
                minusLabel.Spectrum.plot(ax=ax1)
                ax1.figure.show()
                raw_input()

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

                print plusB, minusB
                ax1.clear()
                plusLabel.Spectrum.plot(ax=ax1)
                minusLabel.Spectrum.plot(ax=ax1)
                ax1.figure.show()
                raw_input()

                


    print blah
    return


configFile = sys.argv[1]
calcIM = bool(sys.argv[2]=='True')
moogInstance = sys.argv[3]
try:
    contFloat = bool(sys.argv[4]== 'True')
except:
    contFloat = False
#try:
#    rotFloat = bool(sys.argv[5] == 'True')
#except:
#    rotFloat = False

f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])

config = AstroUtils.parse_config(configFile)

orchestra = Moog960.Score(directory='../MusicMaker/TWHydra_T3*', observed=config["inFile"], suffix='')

IMs = computeIMs(orchestra)

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
