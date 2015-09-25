from __future__ import division, print_function

import MoogPy
import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import SpectralTools
import pyfits
import emcee
import george
from george import kernels
#import triangle

startR = 40000.0
minR = 30000.0
maxR = 50000.0
maxWlShift = 1.0


def model(params):
    Synth.setSpectrumParameters(params)
    #print params
    return Synth.compute()

def lnprior_base(p):
    numLines = (len(p)-3)/2
    # Are any gf's waaaaay to small/large?
    if numpy.any(p[:numLines] > 2.0) | numpy.any(p[:numLines] < -7.0):
        return -numpy.inf
    # Are any damping parameters crazy?
    if numpy.any(p[numLines:-3] > -5.0) | numpy.any(p[:numLines] < -8.5):
        return -numpy.inf
    # Is the smoothing/resolution at a reasonable value?
    if not(minR < startR*p[-3] < maxR):
        return -numpy.inf
    # Is the wl shift within reasonable values?
    if (abs(p[-2]) > maxWlShift):   # 1 angstrom is probably a good value
        return -numpy.inf
    # Is the continuum factor within reasonable values?
    if not(0.99 < p[-1] < 1.01):
        return -numpy.inf
    return 0.0

def lnlike_gp(p, x, y, yerr):
    a, tau = numpy.exp(p[:2])
    gp = george.GP(a*kernels.Matern32Kernel(tau))
    gp.compute(x, yerr)
    return gp.lnlikelihood(y-model(p[2:]))

def lnprior_gp(p):
    lna, lntau = p[:2]
    if not -5 < lna < 5:
        return -np.inf
    if not -5 < lntau < 5:
        return -np.inf
    return lnprior_base(p[2:])

def lnprob_gp(p, x, y, yerr):
    lp = lnprior_gp(p)
    if not numpy.isfinite(lp):
        return -numpy.inf
    return lp + lnlike_gp(p, x, y, yerr)

def fit_gp(Synth, nwalkers=32):
    initialGuess = [0.0, 0.0] + Synth.getSpectrumParameters()
    initialRanges = [1.0e-8, 1.0e-8] + Synth.getInitialParameterSpreads()
    data = Synth.getObserved()
    ndim = len(initialGuess)
    #p0 = [numpy.array(initialGuess) + 
    #        numpy.array(initialRanges)*numpy.random.randn(ndim)
    #        for i in xrange(nwalkers)]
    p0 = [numpy.array(initialGuess) + 
            1.0e-8*numpy.random.randn(ndim)
            for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, args=data, live_dangerously=True)

    print("Running burn-in")
    firstChain = open('first.dat', 'w')
    firstChain.close()

    counter = 0
    for result in sampler.sample(p0, iterations=100, storechain=False):
        counter += 1
        print('Step : %4d' % counter)
        position = result[0]
        firstChain = open('first.dat', 'a')
        for k in range(position.shape[0]):
            firstChain.write("{0:4d} {1:s}\n".format(k, 
                " ".join(str(v) for v in position[k])))
        firstChain.close()
    lnp = result[1]
    sampler.reset()

    counter = 0
    print("Running Second burn-in")
    secondChain = open('second.dat', 'w')
    secondChain.close()
    p = p0[numpy.argmax(lnp)]
    p0 = [p + 1e-8*numpy.random.randn(ndim) for i in xrange(nwalkers)]
    for result in sampler.sample(p0, iterations=500, storechain=False):
        counter += 1
        print('Second Burn-in, Step : %4d' % counter)
        position = result[0]
        secondChain = open('second.dat', 'a')
        for k in range(position.shape[0]):
            secondChain.write("{0:4d} {1:s}\n".format(k, 
                " ".join(str(v) for v in position[k])))
        secondChain.close()
    sampler.reset()

    print("Running production")
    finalChain = open('final.dat', 'w')
    finalChain.close()
    counter = 0
    for result in sampler.sample(p0, iterations=1000, storechain=False):
        counter += 1
        print('Final, Step : %4d' % counter)
        position = result[0]
        finalChain = open('final.dat', 'a')
        for k in range(position.shape[0]):
            finalChain.write("{0:4d} {1:s}\n".format(k, 
                " ".join(str(v) for v in position[k])))
        finalChain.close()
    #"""
    return sampler


configFile = 'MCStudent.cfg'

#diffplot = True
diffplot = False

Synth = MoogTools.Moog(configFile)
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()

solarWave = Synth.solarSpectrum.wave+0.1
solarFlux = Synth.solarSpectrum.flux+0.001+ numpy.random.randn(len(solarWave))*0.001

continuum = 0.0
wlOffset = 0.0
resolution = 45000


sampler = fit_gp(Synth)

samples = sampler.flatchain

fig = triangle.corner(samples[:,2:])
fig.show()
fig.savefig("Final.png")
