import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import Moog960
import SpectralTools
import pyfits

configFile = 'CMTeacher.cfg'

Synth = MoogTools.MoogStokes(configFile, moogInstance='ALPHA')
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()

Synth.run()
nominal = Synth.Spectra[0]
smoothed = nominal.resample(R=80000)
resampled = nominal.resample(R=80000, nyquist=True)

phrases = []

parameters = {}
parameters["LABEL"] = "CM Teacher"
parameters["TEFF"] = 5770
parameters["LOGG"] = 4.43
parameters["BFIELD"] = 0.0
parameters["SELECTED"] = False
parameters["WLSTART"] = resampled.wl[0]
parameters["WLSTOP"] = resampled.wl[-1]
header = resampled.header
header.set('WLSTART', parameters["WLSTART"])
header.set('WLSTOP', parameters["WLSTOP"])
header.set('TEFF', parameters["TEFF"])
header.set('LOGG', parameters["LOGG"])
header.set('BFIELD', parameters["BFIELD"])
resampled.addHistory(spectrum_type="OBSERVED")
#spectrum = SpectralTools.Spectrum(wl=spectrograph_wl, I=new_y,
#                spectrum_type='OBSERVED', header=header)
phrases.append(Moog960.ObservedPhrase(observedData=
            [resampled], observedLabels = [Moog960.Label(parameters)]))

melody = Moog960.ObservedMelody(phrases=phrases)
melody.selectPhrases(selectAll=True)

spectra, labels = melody.perform()

for label in labels:
    melody.record(labels=label, basename="Teacher")



phrase = Moog960.ObservedPhrase(observedData = [resampled], 
        observedLabels=[Moog960.Label(parameters)])
