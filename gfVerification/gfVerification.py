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

    resolution = 80000.0
    nLines = SynthObj.lineList.numLines
    nModes = nLines+3
    IM = numpy.zeros((nModes, len(solarWl)))
    stroke = numpy.ones(nModes)
    target = numpy.ones(nModes) * 0.1
    factor = numpy.ones(nModes) * 0.3
    factor[-1] = 0.005          # Continuum Shift
    factor[-2] = 0.01           # WL Shift
    factor[-3] = 10000.0            # Resolution
    f3 = pyplot.figure(3)
    ax3 = f3.add_axes([0.1, 0.1, 0.8, 0.8])
    #"""
    for i in range(nLines):
        SynthObj.lineList.writeLineLists(i)
        wave, flux = SynthObj.run()
        target[i] = (1.0 - numpy.min(flux))/10.0
        #if i < SynthObj.lineList.nStrong:
        #    target[i] = (1.0 - flux[numpy.abs(wave-SynthObj.lineList.strongLines[i].wl).argsort()[0]])/10.0
        #else:
        #    target[i] = (1.0 - flux[numpy.abs(wave-SynthObj.lineList.weakLines[i-SynthObj.lineList.nStrong].wl).argsort()[0]])/10.0
        while ((factor[i] < 0.9) | (factor[i] > 1.1)):
            stroke[i] *= factor[i]
            SynthObj.lineList.perturbGf(i, stroke[i])
            SynthObj.lineList.writeLineLists(i, partial=True)
            wave, plus = SynthObj.run()
            wave, plus = SpectralTools.resample(wave, plus, resolution)
            SynthObj.lineList.perturbGf(i, -2.0*stroke[i])
            SynthObj.lineList.writeLineLists(i, partial=True)
            wave, minus = SynthObj.run()
            wave, minus = SpectralTools.resample(wave, minus, resolution)
            SynthObj.lineList.perturbGf(i, stroke[i])
            factor[i] = numpy.abs(target[i]/(numpy.min(plus)-numpy.min(minus)))
            if factor[i] > 1e3:  # Probably doesn't contribute much
                factor[i] = 1.0
                stroke[i] = 0.01
            """
            print i, target[i], factor[i], stroke[i], SynthObj.lineList.weakLines[i-SynthObj.lineList.nStrong].loggf
            print numpy.min(plus), numpy.min(minus)
            raw_input()
            #"""
        stroke[i] *= factor[i]
        """
        print stroke
        print factor
        print i
        raw_input()
        #"""

    for i in range(nLines):
        SynthObj.lineList.perturbGf(i, stroke[i])
        SynthObj.lineList.writeLineLists(i)
        wave, plus = SynthObj.run()
        newWave, plus = SpectralTools.resample(wave, plus, resolution)
        SynthObj.lineList.perturbGf(i, -2.0*stroke[i])
        SynthObj.lineList.writeLineLists(i)
        wave, minus = SynthObj.run()
        newWave, minus = SpectralTools.resample(wave, minus, resolution)
        IM[i, :] = SpectralTools.interpolate_spectrum(newWave, solarWl,
                (plus-minus)/(2.0*stroke[i]), pad=True)

    #"""
    #Continuum Level
    stroke[-1] *= factor[-1]
    plus = nominalSpectrum + stroke[-1]
    newWave, plus = SpectralTools.resample(wave, plus, resolution)
    minus = nominalSpectrum - stroke[-1]
    newWave, minus = SpectralTools.resample(wave, minus, resolution)
    factor[-1] = 1.0

    IM[-1, :] = SpectralTools.interpolate_spectrum(newWave, solarWl, (plus-minus)/(2.0*stroke[-1]), pad=True)
    
    #Wavelength Shift
    stroke[-2] *= factor[-2]
    plus = SpectralTools.interpolate_spectrum(wave, wave-stroke[-2],
            nominalSpectrum, pad=True)
    newWave, plus = SpectralTools.resample(wave, plus, resolution)
    minus = SpectralTools.interpolate_spectrum(wave, wave+stroke[-2],
            nominalSpectrum, pad=True)
    newWave, minus = SpectralTools.resample(wave, minus, resolution)
    factor[-2] = 1.0
    IM[-2, :] = SpectralTools.interpolate_spectrum(newWave, solarWl, (plus-minus)/(2.0*stroke[-2]), pad=True)

    #Instrumental smoothing
    stroke[-3] *= factor[-3]
    wavePlus, plus = SpectralTools.resample(wave, nominalSpectrum, resolution + stroke[-3])
    waveMinus, minus = SpectralTools.resample(wave, nominalSpectrum, resolution - stroke[-3])
    diffx, diffy = SpectralTools.diff_spectra(wavePlus, plus, waveMinus, minus, pad=True)
    IM[-3,:] = SpectralTools.interpolate_spectrum(wavePlus, solarWl, diffy/(2.0*stroke[-3]), pad=True)

    hdu = pyfits.PrimaryHDU(IM)
    hdu.writeto("./Output/InteractionMatrix.fits", clobber=True)

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
    hdu.writeto('./Output/CommandMatrix.fits', clobber=True)



configFile = 'gfVerification.cfg'

Synth = MoogTools.Moog(configFile)
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()
Synth.solarSpectrum.fixDiscontinuities()
Synth.solarSpectrum.flipWavelength()

solarWave = Synth.solarSpectrum.wave + 0.2
solarFlux = Synth.solarSpectrum.flux

wave, nominalSpectrum = Synth.run()

#if True:
if False:
    calculateCM(Synth, wave, nominalSpectrum)

CM = pyfits.getdata("./Output/CommandMatrix.fits")

f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])

ax1.matshow(CM, aspect='auto')
f1.show()

Spectra = [SpectralTools.interpolate_spectrum(wave, solarWave, nominalSpectrum, pad=True)]

Gain = 0.15

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
    
    Synth.lineList.applyCorrection(command[:-3])
    continuum = continuum+command[-1]
    wlOffset = wlOffset+command[-2]
    resolution = resolution+command[-3]
    wave, flux = Synth.run()
    Spectra.append(SpectralTools.interpolate_spectrum(wave+wlOffset, solarWave,
        flux+continuum, pad=True))

    ax2.clear()
    ax2.plot(solarWave, Synth.solarSpectrum.flux)
    for spec in Spectra:
        ax2.plot(solarWave, spec)
    f2.show()
    print wlOffset, command[-2]
    print continuum, command[-1]
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
