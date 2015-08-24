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



configFile = 'CMTeacher.cfg'

Synth = MoogTools.Moog(configFile)
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()

wavelengths, nominalSpectrum = Synth.run()

SpectralTools.write_2col_spectrum("observed.dat", wavelengths, nominalSpectrum)
