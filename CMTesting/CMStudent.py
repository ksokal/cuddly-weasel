import MoogPy
import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import SpectralTools
import pyfits

def calculateCM(SynthObj, nominalWave, nominalSpectrum):
    solarWl = SynthObj.solarSpectrum.wave
    solarFl = SynthObj.solarSpectrum.flux
    IM = numpy.zeros((SynthObj.lineList.numLines+2, len(solarWl)))
    stroke = numpy.ones(SynthObj.lineList.numLines+2)
    factor = numpy.ones(SynthObj.lineList.numLines+2) * 0.3
    factor[-1] = 0.005
    factor[-2] = 0.1
    while ((factor < 0.9) | (factor > 1.1)).any():
        stroke *= factor
        for i in range(SynthObj.lineList.numLines):
            SynthObj.lineList.perturbLine(i, stroke[i])
            wave, plus = SynthObj.run()
            SynthObj.lineList.perturbLine(i, -stroke[i])
            wave, minus = SynthObj.run()
            factor[i] = numpy.abs(0.1/numpy.min(plus-minus))
            IM[i,:] = SpectralTools.interpolate_spectrum(wave, solarWl,
                    (plus - minus)/(2.0*stroke[i]))

        #Continuum Level
        plus = nominalSpectrum + stroke[-1]
        minus = nominalSpectrum - stroke[-1]
        factor[-1] = 1.0
        IM[-1, :] = SpectralTools.interpolate_spectrum(wave, solarWl, 
                (plus-minus)/0.01)

        plus = SpectralTools.interpolate_spectrum(wave, 
                wave - stroke[-2], nominalSpectrum, pad=True)
        minus = SpectralTools.interpolate_spectrum(wave, 
                wave + stroke[-2], nominalSpectrum, pad=True)
        factor[-2] = 1.0
        IM[-2, :] = SpectralTools.interpolate_spectrum(wave, solarWl, 
                (plus-minus)/(0.2))

        print stroke
        print factor

    hdu = pyfits.PrimaryHDU(IM)
    hdu.writeto("InteractionMatrix.fits", clobber=True)
    #nFilt = Synth.lineList.numLines-2
    dims = IM.shape
    U,S,V = scipy.linalg.svd(IM)
    D = 1.0/S
    #D = 1.0/(S[0:-nFilt])
    #S[-nFilt:] = 0.0
    newS = numpy.zeros((dims[0], dims[1]))
    I = [i for i in range(dims[1])]
    for i in range(len(D)):
        newS[i][i] = D[i]

    S = newS.copy()
    CM = numpy.array(scipy.matrix(V.T.dot(S.T.dot(U.T)),dtype=numpy.float32)).T

    hdu = pyfits.PrimaryHDU(CM)
    hdu.writeto('CommandMatrix.fits', clobber=True)



configFile = 'CMStudent.cfg'

Synth = MoogTools.Moog(configFile)
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()

solarWave = Synth.solarSpectrum.wave
solarFlux = Synth.solarSpectrum.flux


wavelengths, nominalSpectrum = Synth.run()


if True:
    calculateCM(Synth, wavelengths, nominalSpectrum)

CM = pyfits.getdata("CommandMatrix.fits")

f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])

cb = ax1.matshow(CM, aspect='auto')
blah = pyplot.colorbar(cb)
f1.show()

Spectra = [SpectralTools.interpolate_spectrum(wavelengths, solarWave, nominalSpectrum)]

Gain = 0.35

f2 = pyplot.figure(1)
f2.clear()
ax2 = f2.add_axes([0.1, 0.1, 0.8, 0.8])

continuum = 0.0
wlOffset = 0.0

while True:
    x, difference = SpectralTools.diff_spectra(solarWave, solarFlux, 
            solarWave, Spectra[-1], pad=True)
    difference[Spectra[-1] == 0] = 0.0
    command = Gain*(CM.dot(difference))
    
    Synth.lineList.applyCorrection(command[:-2])
    continuum = continuum+command[-1]
    wlOffset = wlOffset+command[-2]
    wave, flux = Synth.run()
    #Spectra.append(wave)
    Spectra.append(SpectralTools.interpolate_spectrum(wave+wlOffset, solarWave, 
        flux+continuum, pad=True))

    ax2.clear()
    ax2.plot(solarWave, Synth.solarSpectrum.flux)
    ax2.plot(wave, flux)
    for spec in Spectra:
        ax2.plot(solarWave, spec)
    f2.show()
    print Synth.lineList.weakLines[0].loggf
    print Synth.lineList.weakLines[1].loggf
    #print wlOffset, command[-2]
    raw_input()


