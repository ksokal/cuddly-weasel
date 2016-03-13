import scipy
import numpy
import SpectralTools
import Moog960
import matplotlib.pyplot as pyplot
import pyfits

#def argsort(seq):
#    return sorted(range(len(seq)), key=seq.__getitem__)

def parseParams(labels):
    keys = labels[0].parameters.keys()
    params = {}
    for key in keys:
        params[key] = []

    for label in labels:
        for key in keys:
            params[key].append(label.parameters[key])

    for key in keys:
        params[key] = numpy.unique(params[key])
    return params

def getGridPoints(params, convolved, wlStart, wlStop):
    desired = {}
    gridPoints = {}
    for key in ["TEFF", "LOGG", "BFIELD"]:
        points = params[key]
        low = numpy.min(points)
        high = numpy.max(points)
        while True:
            val = raw_input("Enter desired value for %.2f < %s > %.2f: " %
                    (low, key, high))
            val = float(val)
            if (low <= val) & (val <= high):
                break
        desired[key] = val
        gridPoints[key] = [numpy.sort(points[points<=val])[-1], 
                numpy.sort(points[points >=val])[0]]
        print gridPoints[key], key

    selected = []
    for label in convolved:
        if label.parameters["TEFF"] in gridPoints["TEFF"]:
            if label.parameters["LOGG"] in gridPoints["LOGG"]:
                if label.parameters["BFIELD"] in gridPoints["BFIELD"]:
                    selected.append(label)

    return selected, desired

def interpSpectra(sp1, sp2, param1, param2, desired):
    retval = []
    #newp = p1[0].copy()
    newParameters = []
    fraction = -1.0
    order1 = numpy.argsort(numpy.array(param1))
    order2 = numpy.argsort(numpy.array(param2))
    print order1
    print order2
    param1 = numpy.array(param1)[order1].tolist()
    param2 = numpy.array(param2)[order2].tolist()
    sp1 = numpy.array(sp1)[order1].tolist()
    sp2 = numpy.array(sp2)[order2].tolist()
    for p1, p2 in zip(param1, param2):
        newp = p1.copy()
        for key in desired.keys():
            if p1.parameters[key] != p2.parameters[key]:
                distance = p1.parameters[key] - p2.parameters[key]
                fraction = (p1.parameters[key] - desired[key])/distance
                newp.parameters[key] = desired[key]
        newParameters.append(newp)
    if fraction != -1.0:
        for s1, s2 in zip(sp1, sp2):
            retval.append(s1.blend(s2, fraction))
    else:
        for s1 in sp1:
            retval.append(s1)

    return retval, newParameters

def interpolateSpectra(spectra, params, labels, desired):
    # Keep Teff the same, interpolate B-field
    Bspectra = []
    Bparams = []
    if len(params) == 1:
        Bspectra.append(spectra[0])
        Bparams.append(params[0])
    else:
        while len(params) > 0:
            sp1 = spectra.pop(0)
            p1 = params.pop(0)
            p2 = p1
            for i in range(len(params)):
                if ( (params[i][0].parameters["TEFF"] == p1[0].parameters["TEFF"]) & 
                        (params[i][0].parameters["LOGG"] == p1[0].parameters["LOGG"])):
                    sp2 = spectra.pop(i)
                    p2 = []
                    p2 = params.pop(i)
                    break

            if p2 != p1:
                newsp, newp = interpSpectra(sp1, sp2, p1, p2, desired)
                Bspectra.append(newsp)
                Bparams.append(newp)
            else:
                Bspectra.append(sp1)
                Bparams.append(p1)

    Gspectra = []
    Gparams = []
    if len(Bparams) == 1:
        Gspectra.append(Bspectra[0])
        Gparams.append(Bparams[0])
    else:
        while len(Bparams) > 0:
            sp1 = Bspectra.pop(0)
            p1 = Bparams.pop(0)
            p2 = p1
            for i in range(len(Bparams)):
                if ( (Bparams[i][0].parameters["TEFF"] == p1[0].parameters["TEFF"]) ):
                    sp2 = Bspectra.pop(i)
                    p2 = []
                    p2 = Bparams.pop(i)
                    break
            if p2 != p1:
                newsp, newp = interpSpectra(sp1, sp2, p1, p2, desired)
                Gspectra.append(newsp)
                Gparams.append(newp)
            else:
                Gspectra.append(sp1)
                Gparams.append(p1)

    if len(Gparams) == 1:
        interpolatedSpectrum = Gspectra[0]
        interpolatedParams = Gparams[0]
    else:
        interpolatedSpectrum, interpolatedParams = interpSpectra(Gspectra[0],
            Gspectra[1], Gparams[0], Gparams[1], desired)

    phrases = []
    for sp, p in zip(interpolatedSpectrum, interpolatedParams):
        phrases.append(Moog960.SyntheticPhrase(convolvedData=[sp], convolvedLabels=[p], diskInt='BEACHBALL'))
    header = pyfits.Header()
    for key in desired.keys():
        header.set(key, desired[key])
    header.set("SPECTRUM_CONTENTS", "CONVOLVED")
    blended = Moog960.SyntheticMelody(phrases=phrases, header=header)
    return blended, interpolatedParams
    
#blended = interpolateSpectra(spectra, params, labels, desired)
fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

orchestra = Moog960.Score(directory = './TWHydra', suffix = '')
raw, interpolated, integrated, convolved, observed = orchestra.getMelodyParams()

parameters = parseParams(convolved)

wlStart = 19000
wlStop = 23800

while True:
    selected, desired = getGridPoints(parameters, convolved, wlStart, wlStop)

    orchestra.selectMelodies(wlRange=[wlStart, wlStop])
    orchestra.selectEnsemble(selectedLabels= selected)
    """
TODO: Make the Score.blend(desired) function.

This function should do the blending, and save the output in the orchestra object.
For now, this should be considered a hack.
    """
    spectra, params, labels= orchestra.perform(selectedLabels=selected)

    blended, labels = interpolateSpectra(spectra[:], params[:], labels, desired)

    ax.clear()
    for sp in spectra:
        for s in sp:
            ax.plot(s.wl, s.flux_I)

    labels = []
    for phrase in blended.phrases:
        for convLabel in phrase.convolvedLabels:
            labels.append(convLabel)

    blended.selectPhrases(selectAll=True)

    for blend in labels:
        b, l = blended.perform(label=blend)
        ax.plot(b.wl, b.flux_I, color = 'k', lw=2.0)

    blended.record(labels=labels, basename="blended")
    fig.show()

