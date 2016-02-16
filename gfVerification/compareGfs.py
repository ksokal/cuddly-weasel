import matplotlib.pyplot as pyplot
import numpy
import scipy
import sys
import MoogTools
import AstroUtils

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

baseName = sys.argv[1]

arcturusConfig = AstroUtils.parse_config(baseName+".solar.cfg")
solarConfig = AstroUtils.parse_config(baseName+".arcturus.cfg")
originalConfig = AstroUtils.parse_config(baseName+".cfg")



arcturus = MoogTools.LineList(None, arcturusConfig)
solar = MoogTools.LineList(None, solarConfig)
original = MoogTools.LineList(None, originalConfig)

orig = []
sol = []
arc = []
diff = []
dsol = []
darc = []
species = []
expot = []

for i in range(original.nStrong):
    orig.append(original.strongLines[i].loggf)
    sol.append(solar.strongLines[i].loggf)
    arc.append(arcturus.strongLines[i].loggf)
    diff.append(sol[-1] - arc[-1])
    dsol.append(sol[-1] - orig[-1])
    darc.append(arc[-1] - orig[-1])
    species.append(arcturus.strongLines[i].species)
    expot.append(arcturus.strongLines[i].expot_lo)

for i in range(original.numLines - original.nStrong):
    orig.append(original.weakLines[i].loggf)
    sol.append(solar.weakLines[i].loggf)
    arc.append(arcturus.weakLines[i].loggf)
    diff.append(sol[-1] - arc[-1])
    dsol.append(sol[-1] - orig[-1])
    darc.append(arc[-1] - orig[-1])
    species.append(arcturus.weakLines[i].species)
    expot.append(arcturus.weakLines[i].expot_lo)

diff = numpy.array(diff)
orig = numpy.array(orig)
sol = numpy.array(sol)
arc = numpy.array(arc)
dsol = numpy.array(dsol)
darc = numpy.array(darc)
species = numpy.array(species)
expot = numpy.array(expot)

changed = diff != 0.0

#ax.plot([-6, 0], [-6, 0])
ax.scatter(orig[changed], sol[changed], color = 'b', label='Solar', s=species[changed])
ax.scatter(orig[changed], arc[changed], color = 'r', label='Arcturus', s=species[changed])
ax.legend(loc=2)
ax.set_xbound(-4.0, 1.0)
ax.set_ybound(-4.0, 1.0)
ax.set_xlabel("Original log gf")
ax.set_ylabel("New log gf")
#ax.scatter(expot[changed], diff[changed])

#ax.scatter(orig[changed], dsol[changed], color = 'r')
#ax.scatter(orig[changed], darc[changed], color = 'b')
#ax.plot([0.0, 8.0], [0.0, 0.0])
#ax.scatter(species[changed]+1, orig[changed], color = 'k')
#ax.scatter(species[changed], arc[changed], color = 'b')
#ax.scatter(species[changed], sol[changed], color = 'r')
fig.show()

fig.savefig("loggf_changes.png")
