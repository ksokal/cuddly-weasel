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
    nModes = 2*nLines+3
    nFilt = nLines

    IM = numpy.zeros((nModes, len(solarWl)))
    stroke = numpy.ones(nModes)
    stroke[-1] = 0.005          # Continuum Shift
    stroke[-2] = 0.1           # WL Shift
    stroke[-3] = 0.01            # Resolution
    for i in range(nLines):
        SynthObj.lineList.writeLineLists(i)
        wave, flux = SynthObj.run()
        stroke[i*2] = SynthObj.lineList.getGf(i)/5.0
        stroke[i*2+1] = 0.5

    for i in range(nLines):
        #  log gf
        SynthObj.lineList.perturbGf(i, stroke[i*2])
        SynthObj.lineList.writeLineLists(i)
        wave, plus = SynthObj.run()
        newWave, plus = SpectralTools.resample(wave, plus, resolution)
        SynthObj.lineList.perturbGf(i, -2.0*stroke[i*2])
        SynthObj.lineList.writeLineLists(i)
        wave, minus = SynthObj.run()
        newWave, minus = SpectralTools.resample(wave, minus, resolution)
        SynthObj.lineList.perturbGf(i, stroke[i*2])
        IM[i*2, :] = SpectralTools.interpolate_spectrum(newWave, solarWl,
                (plus-minus)/(2.0*stroke[i*2]), pad=0.0)

        #  Damping
        SynthObj.lineList.perturbVdW(i, stroke[i*2+1])
        SynthObj.lineList.writeLineLists(i)
        wave, plus = SynthObj.run()
        newWave, plus = SpectralTools.resample(wave, plus, resolution)
        SynthObj.lineList.perturbVdW(i, -2.0*stroke[i*2+1])
        SynthObj.lineList.writeLineLists(i)
        wave, minus = SynthObj.run()
        newWave, minus = SpectralTools.resample(wave, minus, resolution)
        SynthObj.lineList.perturbGf(i, stroke[i*2+1])
        IM[i*2+1, :] = SpectralTools.interpolate_spectrum(newWave, solarWl,
                (plus-minus)/(2.0*stroke[i*2]), pad=0.0)

    #Continuum Level
    plus = nominalSpectrum + stroke[-1]
    newWave, plus = SpectralTools.resample(nominalWave, plus, resolution)
    minus = nominalSpectrum - stroke[-1]
    newWave, minus = SpectralTools.resample(nominalWave, minus, resolution)

    IM[-1, :] = SpectralTools.interpolate_spectrum(newWave, solarWl, (plus-minus)/(2.0*stroke[-1]), pad=0.0)
    
    #Wavelength Shift
    plus = SpectralTools.interpolate_spectrum(wave, wave-stroke[-2],
            nominalSpectrum, pad=True)
    newWave, plus = SpectralTools.resample(wave, plus, resolution)
    minus = SpectralTools.interpolate_spectrum(wave, wave+stroke[-2],
            nominalSpectrum, pad=True)
    newWave, minus = SpectralTools.resample(wave, minus, resolution)
    IM[-2, :] = SpectralTools.interpolate_spectrum(newWave, solarWl, (plus-minus)/(2.0*stroke[-2]), pad=True)

    #Instrumental smoothing
    wavePlus, plus = SpectralTools.resample(wave, nominalSpectrum, resolution*(1.0 + stroke[-3]))
    waveMinus, minus = SpectralTools.resample(wave, nominalSpectrum, resolution*(1.0 - stroke[-3]))
    diffx, diffy = SpectralTools.diff_spectra(wavePlus, plus, waveMinus, minus, pad=True)
    IM[-3,:] = SpectralTools.interpolate_spectrum(diffx, solarWl, diffy/(2.0*stroke[-3]), pad=0.0)

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
#configFile = 'small.cfg'

Synth = MoogTools.Moog(configFile)
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()
Synth.solarSpectrum.fixDiscontinuities()
Synth.solarSpectrum.flipWavelength()

continuum = 0.0
wlOffset = 0.0
resolution = 580000.0

diffplot = False

solarWave = Synth.solarSpectrum.wave  + 0.2
solarFlux = Synth.solarSpectrum.flux

nominalWavelength, nominalSpectrum = Synth.run()
wave, spectrum = SpectralTools.resample(nominalWavelength, nominalSpectrum, resolution)

#if False:
if True:
    calculateCM(Synth, nominalWavelength, nominalSpectrum)

CM = pyfits.getdata("./Output/CommandMatrix.fits")

f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])

ax1.matshow(CM, aspect='auto')
f1.show()

Spectra = [spectrum]
Wavelengths = [wave]

Gain = 0.35

f2 = pyplot.figure(1)
f2.clear()
ax2 = f2.add_axes([0.1, 0.1, 0.8, 0.8])

raw_input()

while True:
    x, difference = SpectralTools.diff_spectra(solarWave, solarFlux,
            Wavelengths[-1], Spectra[-1], pad=True)
    difference[Spectra[-1] == 0] = 0.0
    command = Gain*(CM.dot(difference))
    
    Synth.lineList.applyCorrection(command[:-2])
    continuum = continuum+command[-1]
    wlOffset = wlOffset+command[-2]
    resolution = resolution*(1.0+command[-3])
    wave, flux = Synth.run()

    wavelength, spectrum = SpectralTools.resample(wave, flux, resolution)
    Spectra.append(spectrum+continuum)
    Wavelengths.append(wavelength+wlOffset)

    ax1.clear()
    ax1.plot(command)
    f1.show()

    ax2.clear()
    if diffplot:
        x, difference = SpectralTools.diff_spectra(solarWave, solarFlux, 
                Wavelengths[-1], Spectra[-1], pad = True)
        difference[Spectra[-1] == 0.0] = 0.0
        ax2.plot(x, difference)
    else:
        ax2.plot(solarWave, solarFlux)
        for w, spec in zip(Wavelengths, Spectra):
            ax2.plot(w, spec)
    f2.show()
    print resolution, command[-3]
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
