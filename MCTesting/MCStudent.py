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
    print params
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
    p0 = [numpy.array(initialGuess) + 
            numpy.array(initialRanges)*numpy.random.randn(ndim)
            for i in xrange(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob_gp, args=data, live_dangerously=True)

    print("Running burn-in")
    p0, lnp, _ = sampler.run_mcmc(p0, 100)
    """
    sampler.reset()

    print("Running Second burn-in")
    p = p0[numpy.argmax(lnp)]
    p0 = [p + 1e-2*numpy.random.randn(ndim) for i in xrange(nwalkers)]
    p0, _, _ = sampler.run_mcmc(p0, 500)
    raw_input()
    sampler.reset()

    print("Running production")
    p0, _, _ = sampler.run_mcmc(p0, 1000)
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
fig.savefig("First_Burn_In.png")
