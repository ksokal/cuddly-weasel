import matplotlib.pyplot as pyplot
import numpy
import scipy
import sys
import MoogTools
import Moog960
import AstroUtils

fig1 = pyplot.figure(0)
fig1.clear()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])

fig2 = pyplot.figure(1)
fig2.clear()
ax2 = fig2.add_axes([0.1, 0.1, 0.8, 0.8])

baseName = sys.argv[1]

arcConfig = baseName+".arcturus.cfg"
solConfig = baseName+".solar.cfg"

arcturus = MoogTools.MoogStokes(arcConfig, moogInstance = "ALPHA")
solar = MoogTools.MoogStokes(solConfig, moogInstance = "BRAVO")

solarSpectrum = Moog960.ObservedMelody.fromFile(filename=solar.config["inFile"])
arcturusSpectrum = Moog960.ObservedMelody.fromFile(filename=arcturus.config["inFile"])
solarSpectrum.loadData()
arcturusSpectrum.loadData()

wlRange = [solar.config["wlStart"], solar.config["wlStop"]]
solarSpectrum.selectPhrases(wlRange=wlRange)
obsSol, labels = solarSpectrum.perform()
arcturusSpectrum.selectPhrases(wlRange=wlRange)
obsArc, labels = arcturusSpectrum.perform()

solWlShift =  solar.config["wlShift"]
obsSol = obsSol[0][0]
obsSol.wl += solWlShift
window = ((obsSol.wl > wlRange[0]) & (obsSol.wl < wlRange[1]))
obsSol.wl = obsSol.wl[window]
obsSol.flux_I = obsSol.flux_I[window]

arcWlShift = arcturus.config["wlShift"]
obsArc = obsArc[0][0]
obsArc.wl += arcWlShift
window = ((obsArc.wl > wlRange[0]) & (obsArc.wl < wlRange[1]))
obsArc.wl = obsArc.wl[window]
obsArc.flux_I = obsArc.flux_I[window]


arcturus.run()
solar.run()

obsSol.plot(ax = ax1, color = 'b', lw=2.0, label='Observed Solar')
solar.Spectra[-1].plot(ax=ax1, color = 'r', lw=2.0, label='Arcturus gf\'s')
ax1.legend(loc=3)
ax1.set_ybound(0.5, 1.05)
fig1.show()
fig1.savefig("ObservedSolar.png")

obsArc.plot(ax = ax2, color = 'b', lw=2.0, label='Observed Arcturus')
arcturus.Spectra[-1].plot(ax=ax2, color = 'r', lw=2.0, label='Solar gs\'s')
ax2.legend(loc=3)
ax2.set_ybound(0.45, 1.05)
fig2.show()
fig2.savefig("ObservedArcturus.png")

