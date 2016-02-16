import matplotlib.pyplot as pyplot
import Moog960
import MoogTools
import numpy
import glob

correctionFiles = glob.glob('/home/deen/MoogPyData/AbsorptionLines/corrections/*Solar_results.dat')

solarSpectrum = Moog960.ObservedMelody.fromFile(filename='SolarSpectrum.fits')
solarSpectrum.loadData()

arcturusSpectrum = Moog960.ObservedMelody.fromFile(filename='ArcturusSpectrum.fits')
arcturusSpectrum.loadData()

solarSpectrum.selectPhrases(wlRange=[20000, 24000])
arcturusSpectrum.selectPhrases(wlRange=[20000, 24000])
observed, labels = solarSpectrum.perform()
solar = observed[0][0]
observed, labels = arcturusSpectrum.perform()
arcturus = observed[0][0]

wls = []
ycoord = []
left = []
width = []
height = []

counter = 0
for corrections in correctionFiles:
    with open(corrections, 'r') as f:
        lines = f.readlines()
        for line in lines:
            wls.append(float(line.split()[0]))
            ycoord.append(1.0)

    #print counter
    #print wls[counter], wls[-1], wls[-2]
    left.append(wls[counter])
    width.append(wls[-1] - wls[counter])
    height.append(1.0)
    counter = len(wls)



fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

ax.plot(solar.wl, solar.flux_I, color = 'k')
ax.plot(arcturus.wl-5.5, arcturus.flux_I, color ='g')

#ax.scatter(wls, ycoord, marker='o', s=30, color = 'r')
ax.bar(left, height, width=width, alpha=0.5, color = 'r')
ax.set_xbound(20000, 23000)

ax.set_ybound(0.0, 1.1)

fig.show()
