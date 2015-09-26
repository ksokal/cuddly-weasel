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
import triangle

startR = 40000.0
minR = 30000.0
maxR = 50000.0
maxWlShift = 1.0


def model(params):
    Synth.setSpectrumParameters(params)
    #print params
    return Synth.compute()

def lnprior_base(p):
    numLines = int((len(p)-3)/2)
    # Are any gf's waaaaay to small/large?
    if numpy.any(p[:numLines] > 2.0) | numpy.any(p[:numLines] < -7.0):
        return -numpy.inf
    # Are any damping parameters crazy?
    if numpy.any(p[numLines:-3] > -5.0) | numpy.any(p[numLines:-3] < -8.5):
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
    diff = y-model(p[2:])
    print(" ".join([str(v) for v in p[2:]]))
    return gp.lnlikelihood(diff)

def lnprior_gp(p):
    lna, lntau = p[:2]
    if not -6 < lna < 5:
        return -numpy.inf
    if not -5 < lntau < 5:
        return -numpy.inf
    return lnprior_base(p[2:])

def lnprob_gp(p, x, y, yerr):
    lp = lnprior_gp(p)
    if not numpy.isfinite(lp):
        return -numpy.inf
    return lp + lnlike_gp(p, x, y, yerr)

def fit_gp(Synth, nwalkers=128):
    initialGuess = [-5.0, 0.0] + Synth.getSpectrumParameters()
    initialRanges = [1.0e-1, 1.0e-1] + Synth.getInitialParameterSpreads()
    data = Synth.getObserved()
    ndim = len(initialGuess)
    #p0 = [numpy.array(initialGuess) + 
    #        numpy.array(initialRanges)*numpy.random.randn(ndim)
    #        for i in xrange(nwalkers)]
    p0 = [numpy.array(initialGuess) + 
            1.0e-3*numpy.random.randn(ndim)
            for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, args=data, live_dangerously=True)

    print("Running burn-in")
    #p0, lnp, _ = sampler.run_mcmc(p0, 50)
    #"""
    firstChain = open('first.dat', 'w')
    firstChain.close()

    counter = 0
    for result in sampler.sample(p0, iterations=100):
        counter += 1
        print('Step : %4d' % counter)
        position = result[0]
        firstChain = open('first.dat', 'a')
        for k in range(position.shape[0]):
            firstChain.write("%4d %.4f %s\n" % (k, numpy.mean(sampler.acceptance_fraction), " ".join(str(v) for v in position[k])))
        firstChain.close()
        print('%.4f' % numpy.mean(sampler.acceptance_fraction))
    lnp = result[1]
    p = position[numpy.argmax(lnp)]
    #"""

    samples = sampler.flatchain

    sampler.reset()
    f1 = triangle.corner(samples[:,2:])
    f1.savefig("first.png")

    #f2 = pyplot.figure(0)
    #ax = f2.add_axes([0.1, 0.1, 0.8, 0.8])
    #ax.plot(model(p[2:]))


    print("Running Second burn-in")
    p0 = [p + 1e-6*numpy.random.randn(ndim) for i in xrange(nwalkers)]
    
    p0, _, _ = sampler.run_mcmc(p0, 500)

    """
    secondChain = open('second.dat', 'w')
    secondChain.close()
    counter = 0
    for result in sampler.sample(p0, iterations=500, storechain=False):
        counter += 1
        print('Second Burn-in, Step : %4d' % counter)
        position = result[0]
        secondChain = open('second.dat', 'a')
        for k in range(position.shape[0]):
            secondChain.write("%4d %.4f %s\n" % (k, numpy.mean(sampler.acceptance_fraction), " ".join(str(v) for v in position[k])))
        secondChain.close()
    #"""
    sampler.reset()

    print("Running production")

    p0, _, _ = sampler.run_mcmc(p0, 1000)
    return sampler


configFile = 'MCStudent.cfg'

#diffplot = True
diffplot = False

Synth = MoogTools.Moog(configFile)
Synth.lineList.writeLineLists()
Synth.parameterFile.writeParFile()

gp = george.GP(0.00001*kernels.ExpSquaredKernel(4.3))

solarWave = Synth.solarSpectrum.wave+0.1
solarFlux = Synth.solarSpectrum.flux+0.001+ numpy.random.randn(len(solarWave))*0.001+gp.sample(solarWave)

#fig = pyplot.figure(0)
#fig.clear()
#ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#ax.plot(solarWave, solarFlux)
#fig.show()
#raw_input()
#print( asdf)
continuum = 0.0
wlOffset = 0.0
resolution = 45000


sampler = fit_gp(Synth)

samples = sampler.flatchain

fig = triangle.corner(samples[:,2:])
fig.savefig("Final.png")
