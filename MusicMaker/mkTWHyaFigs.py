import Moog960
import SpectralTools
import matplotlib.pyplot as pyplot
import scipy
import scipy.optimize
import numpy

def fit(synthetic, observed):
    veiling = 0.1
    slope = 0.0001
    continuum = 1.0

    params = [veiling, slope, continuum]


    def fitfunc(pars, synth):
        xpts = numpy.arange(len(synth))

        retval = (synth + pars[0])/(1.0+pars[0])*pars[2] + pars[1]*xpts

        return retval
        

    def errfunc(pars, synth, obs):
        return numpy.abs(fitfunc(pars, synth) - obs)

    bestFit, success = scipy.optimize.leastsq(errfunc, params, args=(synthetic, observed))

    #return (synthetic+bestFit[0])/(1.0+bestFit[0])*bestFit[2] + bestFit[1]*xpts
    return bestFit

fig = pyplot.figure(0)
fig.clear()
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.4])

#wlStart = 22880.0
#wlStop = 24000.0
#wlStart = 22700.0
#wlStop = 22890.0
#wlStart = 22100.0
#wlStop = 22550.0
#wlStart = 21700.0
#wlStop = 21980.0
#wlStart = 22010.0
#wlStop = 22550.0
#wlStart = 20000.0
#wlStop = 21000.0
#wlStart = 21000.0
#wlStop = 22000.0
wlStart = 22000.0
wlStop = 22800.0

filename = '/home/deen/Investigations/TWHydra/TWHydra.fits'
TWHya = Moog960.ObservedMelody.fromFile(filename=filename, label='IGRINS TWHydra')

TWHya.selectPhrases(wlRange = [wlStart, wlStop])
TWHya.loadData()
observed_spectra, observed_label = TWHya.perform()

for obs in observed_spectra[0]:
    obs.wl = obs.wl - 5.5

    region = (obs.wl > wlStart) & (obs.wl < wlStop)
    obs.wl = obs.wl[region]
    obs.flux_I = obs.flux_I[region]


orchestra = Moog960.Score(directory='./blended', suffix='')
raw, interpolated, integrated, convolved, observed = orchestra.getMelodyParams()

orchestra.selectMelodies(wlRange=[wlStart, wlStop])
orchestra.selectEnsemble(selectedLabels=convolved)

spectra, params, labels = orchestra.perform(selectedLabels = convolved)

colors = ['g', 'm', 'c', 'r']

bestFit= []
veiled = []
#for T, G, B, color in zip(Ts, Gs, Bs, colors):
for spectrum, label, color in zip(spectra, params, colors):
    for phrase, l in zip(spectrum, label):
        print phrase.wl[0], phrase.wl[-1]
        if ((phrase.wl[0] < wlStart) & (phrase.wl[-1] > wlStop)): 
            phrase.bin(obs.wl)
            bestFit.append(fit(phrase.flux_I, obs.flux_I))
            #bestFit[-1][0] = 0.3
            print bestFit
            phrase.flux_I = (phrase.flux_I+bestFit[-1][0])/(1.0+bestFit[-1][0])
            ax1.plot(phrase.wl, phrase.flux_I, color = color, lw=2.0)
            v = {}
            v["spectrum"] = phrase
            v["Teff"] = l.parameters["TEFF"]
            v["log g"] = l.parameters["LOGG"]
            v["B"] = l.parameters["BFIELD"]
            v["color"] = color
            v["veiling"] = bestFit[-1][0]
            veiled.append(v)


bestFit = numpy.average(numpy.array(bestFit), axis=0)
#bestContinuum = 
#bestFit[2] + bestFit[1]*xpts
observed_spectra[0][0].flux_I = observed_spectra[0][0].flux_I/bestFit[2] - bestFit[1]*numpy.arange(len(observed_spectra[0][0].wl))
ax1.plot(observed_spectra[0][0].wl, observed_spectra[0][0].flux_I, color = 'k', lw=2.0)

ax1.set_xbound(wlStart, wlStop)
ax1.set_ybound(0.5, 1.1)

for v in veiled:
    sp = v["spectrum"]
    difference = sp - observed_spectra[0][0]
    ax2.plot(difference.wl, difference.flux_I, color = v["color"], lw = 2.0)
    print v["Teff"], v["B"], v["veiling"]
fig.show()

#del(orchestra)
#del(spectra)
#del(params)
#del(labels)
