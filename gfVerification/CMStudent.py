import MoogStokesPy_Alpha
import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import Moog960
import SpectralTools
import pyfits
import random
import string
import sys

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

def calculateIM(SynthObj, nominalSpectrum, resolution):
    solarWl = nominalSpectrum.wl
    nLines = SynthObj.lineList.numLines
    nModes = 2*nLines+3
    nStrong = SynthObj.lineList.nStrong
    nFilt = nLines
    nPts = len(solarWl)

    IM = numpy.zeros((0, nPts))

    factor = []
    continuumStroke = 0.01                  # Continuum Shift
    wlShiftStroke = 0.5                    # WL Shift
    smoothingStroke = 0.01                # Smoothing
    rotationStroke = 1e-6                 # Spectrum rotation

    counter = 0
    gfIndices = numpy.zeros(nLines, dtype=int)-1
    vdWIndices = numpy.zeros(nLines, dtype=int)-1
    wlStart = SynthObj.parameterFile.config["wlStart"]
    wlStop = SynthObj.parameterFile.config["wlStop"]
    for i in range(nLines):
        SynthObj.Spectra = []
        #  log gf
        stroke = SynthObj.lineList.getGf(i)/3.0
        SynthObj.lineList.perturbGf(i, stroke)
        SynthObj.run(lineIndex = i, partial=True)
        plus = SynthObj.Spectra[-1].resample(R=resolution, observedWl=solarWl)
        SynthObj.lineList.perturbGf(i, -2.0*stroke)
        SynthObj.run(lineIndex = i, partial=True)
        minus = SynthObj.Spectra[-1].resample(R=resolution, observedWl=solarWl)
        SynthObj.lineList.perturbGf(i, stroke)
        if i >= nStrong:
            plus.flux_I[plus.flux_I == 0.0] = 1.0
            minus.flux_I[minus.flux_I == 0.0] = 1.0
        IF = (plus - minus).flux_I/(2.0*stroke)
        gff = numpy.max(numpy.abs(IF))
        strength = numpy.min(plus.flux_I)
        if ((gff > 1e-2) & (strength < 0.995 )):
            counter += 1
            IM = numpy.resize(IM, (counter, nPts))

            overlap_start = numpy.max([numpy.min(solarWl), numpy.min(plus.wl)])
            overlap_stop = numpy.min([numpy.max(solarWl), numpy.max(minus.wl)])
            overlap = scipy.where((solarWl >= overlap_start) & (solarWl <= overlap_stop))
            padded = numpy.zeros(len(solarWl))
            padded[overlap] = IF

            IM[counter-1,:] = padded/gff
            factor.append(gff)
            gfIndices[i] = counter-1

        # damping
        stroke = 1.0
        SynthObj.lineList.perturbVdW(i, stroke)
        SynthObj.run(lineIndex=i, partial=True)
        plus = SynthObj.Spectra[-1].resample(R=resolution, observedWl=solarWl)
        SynthObj.lineList.perturbVdW(i, -2.0*stroke)
        SynthObj.run(lineIndex=i, partial=True)
        minus = SynthObj.Spectra[-1].resample(R=resolution, observedWl=solarWl)
        SynthObj.lineList.perturbVdW(i, stroke)
        if i >= nStrong:
            plus.flux_I[plus.flux_I == 0.0] = 1.0
            minus.flux_I[minus.flux_I == 0.0] = 1.0
        IF = (plus - minus).flux_I/(2.0*stroke)
        vdWf = numpy.max(numpy.abs(IF))
        print("%s Line %.4f  - %d out of %d, strength %.4f" % 
                (SynthObj.config["baseName"], SynthObj.lineList.getWl(i),
                i, nLines, strength))
        if (vdWf > 1e-3):
            counter += 1
            overlap_start = numpy.max([numpy.min(solarWl), numpy.min(plus.wl)])
            overlap_stop = numpy.min([numpy.max(solarWl), numpy.max(minus.wl)])
            overlap = scipy.where((solarWl >= overlap_start) & (solarWl <= overlap_stop))
            padded = numpy.zeros(len(solarWl))
            padded[overlap] = IF

            IM = numpy.resize(IM, (counter, nPts))
            IM[counter-1,:] = padded/vdWf
            factor.append(vdWf)
            vdWIndices[i] = counter-1

    #Continuum Level
    plus = nominalSpectrum.flux_I + continuumStroke
    minus = nominalSpectrum.flux_I - continuumStroke

    counter += 1
    IM = numpy.resize(IM, (counter, nPts))
    IM[-1, :] = (plus-minus)/(2.0*continuumStroke)
    factor.append(numpy.max(numpy.abs(IM[-1,:])))
    IM[-1,:] /= factor[-1]

    #Wavelength Shift
    plus = nominalSpectrum.resample(R=resolution)
    plus.wl += wlShiftStroke
    plus.bin(solarWl)

    minus = nominalSpectrum.resample(R=resolution)
    minus.wl -= wlShiftStroke
    minus.bin(solarWl)
    
    IF = (plus - minus).flux_I/(2.0*wlShiftStroke)
    overlap_start = numpy.max([numpy.min(solarWl), numpy.min(plus.wl)])
    overlap_stop = numpy.min([numpy.max(solarWl), numpy.max(minus.wl)])
    overlap = scipy.where((solarWl >= overlap_start) & (solarWl <= overlap_stop))
    padded = numpy.zeros(len(solarWl))
    padded[overlap] = IF
    
    counter += 1
    IM = numpy.resize(IM, (counter, nPts))
    IM[-1,:] = padded
    factor.append(numpy.max(numpy.abs(IM[-1,:])))
    IM[-1,:] /= factor[-1]
    
    #"""
    #Instrumental Smoothing
    plus = nominalSpectrum.resample(R=resolution*(1.0+smoothingStroke),
            observedWl=solarWl)
    minus = nominalSpectrum.resample(R=resolution*(1.0-smoothingStroke),
            observedWl=solarWl)
    IF = (plus - minus).flux_I/(2.0*smoothingStroke)
    overlap_start = numpy.max([numpy.min(solarWl), numpy.min(plus.wl)])
    overlap_stop = numpy.min([numpy.max(solarWl), numpy.max(minus.wl)])
    overlap = scipy.where((solarWl >= overlap_start) & (solarWl <= overlap_stop))
    padded = numpy.zeros(len(solarWl))
    padded[overlap] = IF
    
    counter += 1
    IM = numpy.resize(IM, (counter, nPts))
    IM[-1,:] = padded
    factor.append(numpy.max(numpy.abs(IM[-1,:])))
    IM[-1,:] /= factor[-1]

    #"""
    # Spectrum Rotation
    plus = nominalSpectrum.rotate(angle=rotationStroke)
    minus = nominalSpectrum.rotate(angle=-rotationStroke)
    IF = (plus - minus).flux_I/(2.0*rotationStroke)

    counter += 1
    IM = numpy.resize(IM, (counter, nPts))
    IM[-1, :] = IF
    factor.append(numpy.max(numpy.abs(IM[-1,:])))
    IM[-1,:] /= factor[-1]
    

    factor = numpy.array(factor)

    hdu = pyfits.PrimaryHDU(IM)
    hdu.writeto(SynthObj.config["outputDir"]+SynthObj.config["baseName"]+"_IM.fits", clobber=True)
    hdu = pyfits.PrimaryHDU(factor)
    hdu.writeto(SynthObj.config["outputDir"]+SynthObj.config["baseName"]+'_factors.fits', clobber=True)
    hdu = pyfits.PrimaryHDU(vdWIndices)
    hdu.writeto(SynthObj.config["outputDir"]+SynthObj.config["baseName"]+'_vdW.fits', clobber=True)
    hdu = pyfits.PrimaryHDU(gfIndices)
    hdu.writeto(SynthObj.config["outputDir"]+SynthObj.config["baseName"]+'_gf.fits', clobber=True)


configFile = sys.argv[1]
calcIM = bool(sys.argv[2]=='True')
moogInstance = sys.argv[3]
try:
    contFloat = bool(sys.argv[4]== 'True')
except:
    contFloat = False
try:
    rotFloat = bool(sys.argv[5] == 'True')
except:
    rotFloat = False

f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])

diffplot = True
diffplot = False

Synth = MoogTools.MoogStokes(configFile, moogInstance=moogInstance, fileBase=''.join(random.choice(string.ascii_letters)
                for _ in range(3)))
#Synth.lineList.writeLineLists()
#Synth.parameterFile.writeParFile()

outFile = Synth.config["outputDir"] + Synth.config["baseName"]+'_results.dat'
inFile = Synth.config["inFile"]

solarSpectrum = Moog960.ObservedMelody.fromFile(filename=inFile)
solarSpectrum.loadData()

wlRange = [Synth.config["wlStart"], Synth.config["wlStop"]]
solarSpectrum.selectPhrases(wlRange=wlRange)
observed, labels = solarSpectrum.perform()

continuum = 0.00
wlOffset = Synth.config["wlShift"]
resolution = Synth.config["resolvingPower"]

solar = observed[0][0]
solar.wl += wlOffset
solar.flux_I += continuum

window = ((solar.wl > wlRange[0]) & (solar.wl < wlRange[1]))

solar.wl = solar.wl[window]
solar.flux_I = solar.flux_I[window]

Synth.run()
nominal = Synth.Spectra[0]
smoothed = nominal.resample(R=resolution)
resampled = nominal.resample(R=resolution, observedWl=solar.wl, pad=1.0)
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
