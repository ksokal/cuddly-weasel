import MoogPy
import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import SpectralTools
import pyfits

def calculateCM(SynthObj, solarWl, rebinnedSpectrum, nominalWave, nominalSpectrum):
    IM = numpy.zeros((SynthObj.lineList.numLines+2, len(solarWl)))
    for i in range(SynthObj.lineList.numLines):
        SynthObj.lineList.perturbLine(i, 0.3)
        wave, flux = SynthObj.run()
        plus = SpectralTools.binSpectrum(flux, wave, solarWl)
        SynthObj.lineList.perturbLine(i, -0.3)
        wave, flux = SynthObj.run()
        minus = SpectralTools.binSpectrum(flux, wave, solarWl)
        IM[i,:] = (plus - minus)/0.6

    #Continuum Level
    plus = rebinnedSpectrum.copy()+ 0.005
    minus = rebinnedSpectrum.copy()- 0.005
    IM[-1, :] = (plus-minus)/0.01

    plus = SpectralTools.binSpectrum(nominalSpectrum, nominalWave+0.1, solarWl)
    minus = SpectralTools.binSpectrum(nominalSpectrum, nominalWave-0.1, solarWl)
    edges = (plus !=0) & (minus != 0)
    IM[-2, edges] = (plus[edges]-minus[edges])/(0.2)

    nFilt = Synth.lineList.numLines-10
    dims = IM.shape
    U,S,V = scipy.linalg.svd(IM)
    D = 1.0/(S[0:-nFilt])
    S[-nFilt:] = 0.0
    newS = numpy.zeros((dims[0], dims[1]))
    I = [i for i in range(dims[1])]
    for i in range(len(D)):
        newS[i][i] = D[i]

    S = newS.copy()
    CM = numpy.array(scipy.matrix(V.T.dot(S.T.dot(U.T)),dtype=numpy.float32)).T

    hdu = pyfits.PrimaryHDU(CM)
    hdu.writeto('CommandMatrix_new.fits', clobber=True)



configFile = 'gfVerification.cfg'

Synth = MoogTools.Moog(configFile)
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()
Synth.solarSpectrum.fixDiscontinuities()
Synth.solarSpectrum.flipWavelength()

solarWave = Synth.solarSpectrum.wave + 0.2


wavelengths, nominalSpectrum = Synth.run()

rebinnedNominalSpectrum = SpectralTools.binSpectrum(nominalSpectrum, wavelengths,
        solarWave)

if False:
    calculateCM(Synth, solarWave, rebinnedNominalSpectrum, wavelengths, nominalSpectrum)

CM = pyfits.getdata("CommandMatrix_new.fits")

f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])

ax1.matshow(CM, aspect='auto')
f1.show()

Spectra = [rebinnedNominalSpectrum.copy()]

Gain = 1.0

f2 = pyplot.figure(1)
f2.clear()
ax2 = f2.add_axes([0.1, 0.1, 0.8, 0.8])

continuum = 0.0
wlOffset = 0.0

while True:
    difference = Spectra[-1] - Synth.solarSpectrum.flux
    difference[Spectra[-1] == 0] = 0.0
    command = -Gain*(CM.dot(difference))
    
    Synth.lineList.applyCorrection(command[:-2])
    continuum = continuum+command[-1]
    wlOffset = wlOffset-command[-2]
    wave, flux = Synth.run()
    Spectra.append(SpectralTools.binSpectrum(flux+continuum,
        wave+wlOffset, solarWave))

    ax2.clear()
    ax2.plot(solarWave, Synth.solarSpectrum.flux)
    for spec in Spectra:
        ax2.plot(solarWave, spec)
    f2.show()
    print wlOffset, command[-2]
    raw_input()


#"""


"""

Line Parameter Fitting -
generate line list
synthesize Solar
synthesize Arcuturs
for line in line list
    push gf
    synthesize solar/arcturus
    pull gf
    synthesize solar/arcturus
    compute influence function for that line
    append to matrix




Star Paramter Fitting
initial guess

"""
