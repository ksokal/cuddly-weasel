import scipy
import numpy
import matplotlib.pyplot as pyplot
import pyfits
import MoogTools

baseNames = ['Alpha', 'Bravo', 'Charlie', 'Delta', 'Echo', 'Foxtrot', 'Golf',
        'Hotel', 'India', 'Juliett', 'Kilo', 'Lima', 'Mike', 'November', 'Oscar']

outName = '/home/deen/MoogPyData/AbsorptionLines/corrections/compositeCorrected.dat'

inDir = './Output/'

solar = []
arcturus = []

for baseName in baseNames:
    inFile = inDir+baseName+'_Solar_results.dat'
    for line in open(inFile, 'r'):
        solar.append(MoogTools.MOOG_Line(line))

    inFile = inDir+baseName+'_Arcturus_results.dat'
    for line in open(inFile, 'r'):
        arcturus.append(MoogTools.MOOG_Line(line))

for arc in arcturus:
    if arc in solar:
        if arc.expot_lo < 3.0:
            index = solar.index(arc)
            solar[index].loggf = arc.loggf
            solar[index].VdW = arc.VdW
    else:
        solar.append(arc)


out = open(outName, 'w')
for line in solar:
    line.dump(out=out, mode='MOOGSCALAR')

out.close()

data = open(outName, 'r').readlines()
wl = []
for line in data:
    wl.append(float(line[0:10]))
order = numpy.argsort(wl)
out = open(outName, 'w')
for i in order:
    out.write(data[i])
out.close()
