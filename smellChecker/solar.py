import MoogTools
import pyfits
import Moog960
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

solar = Moog960.ObservedMelody.fromFile(filename='SolarSpectrum.fits')
solar.loadData()

wlStart = 15000
wlStop = 15100

solar.selectPhrases(wlRange=[wlStart, wlStop])
observed, labels = solar.perform()

observed[0][0].plot(ax=ax)


Synth = MoogTools.MoogStokes('Solar.cfg', wlStart=wlStart, wlStop=wlStop)

Synth.run()
Synth.Spectra[0].plot(ax=ax)

ax.set_xbound(wlStart, wlStop)
fig.show()
