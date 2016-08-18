import scipy
import numpy
import matplotlib.pyplot as pyplot
import Moog960
import SpectralTools
import glob

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

vsini = 5.8
R = 45000.0

datadir = '/home/cdeen/Data/MoogPyGrid/CorrectedRawData/'

files = glob.glob(datadir+'*_raw.fits')

for filename in files:
    print filename

    syntheticMelody = Moog960.SyntheticMelody(filename=filename)
    syntheticMelody.selectPhrases(selectAll=True)
    convolved = syntheticMelody.rehearse(vsini=vsini, R=R, returnLabels=True)
    syntheticMelody.record(labels=convolved, basename='TWHydra')
    del(convolved)
    del(syntheticMelody)

