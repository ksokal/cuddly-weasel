import matplotlib.pyplot as pyplot
import numpy
import scipy
import sys
import MoogTools
import AstroUtils

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#baseNames = ['Alpha', 'Bravo', 'Charlie', 'Delta', 'Echo', 'Foxtrot', 'Golf', 'Hotel',
#        'India', 'Juliett', 'Kilo', 'Lima', 'Mike', 'November', 'Oscar']
#baseNames = ['Alpha', 'Bravo', 'Charlie', 'Echo', 'Golf', 'Hotel',
#        'India', 'Juliett', 'Kilo', 'Lima', 'Mike', 'November', 'Oscar']
baseNames = ['Charlie']

arcGfDiff = []
solGfDiff = []
arcEWs = []
arcEWDiffSelf = []
arcEWDiffOther = []
solEWs = []
solEWDiffSelf = []
solEWDiffOther = []

for baseName in baseNames:
    #Original
    arcturus = MoogTools.MoogStokes(baseName+"_Arcturus.cfg", moogInstance = 'ALPHA', fileBase="arcold")
    solar = MoogTools.MoogStokes(baseName+"_Solar.cfg", moogInstance = 'BRAVO', fileBase="solold")
    #New gfs determined by self
    arcturusNew = MoogTools.MoogStokes(baseName+"_Arcturus.cfg", moogInstance = 'CHARLIE', fileBase="arcnew")
    solarNew = MoogTools.MoogStokes(baseName+"_Solar.cfg", moogInstance = 'DELTA', fileBase="solnew")
    #New gfs determined by other
    #arcSol = MoogTools.MoogStokes(baseName+"_Arcturus.cfg", moogInstance = 'CHARLIE', fileBase="arcsol")
    #solArc = MoogTools.MoogStokes(baseName+"_Solar.cfg", moogInstance = 'DELTA', fileBase="solarc")

    originalGfs = solar.lineList

    #arcturusNew.lineList.gf_corrections = solar.lineList.gf_corrections
    arcturusNew.lineList.applyCorrections = True
    arcturusNew.lineList.readInLineLists()
    solarGfs = arcturusNew.lineList
    #solarNew.lineList.gf_corrections = arcturus.lineList.gf_corrections
    solarNew.lineList.applyCorrections = True
    solarNew.lineList.readInLineLists()
    arcturusGfs = solarNew.lineList

    #arcSol.lineList.applyCorrections = True
    #arcSol.lineList.readInLineLists()
    #solArc.lineList.applyCorrections = True
    #solArc.lineList.readInLineLists()


    numLines = arcturusNew.lineList.numLines

    for i in range(numLines):
        arcDiff = arcturusGfs.getGf(i) - originalGfs.getGf(i)
        solarDiff = solarGfs.getGf(i) - originalGfs.getGf(i)
        if ((arcDiff != 0.0) | (solarDiff != 0.0)):
            print("%s Wavelength = %.4f" % (baseName, arcturus.lineList.getWl(i)))
            arcGfDiff.append(arcDiff)
            arcturus.run(lineIndex = i, partial=True)
            oldSp = arcturus.Spectra[-1].resample(R=65000)
            old = oldSp.calc_EW(oldSp.wl[0], oldSp.wl[-1])
            arcturusNew.lineList = arcturusGfs
            arcturusNew.run(lineIndex = i, partial=True)
            newSp = arcturusNew.Spectra[-1].resample(R=65000)
            arcGfEW = newSp.calc_EW(newSp.wl[0], oldSp.wl[-1])
            arcturusNew.lineList = solarGfs
            arcturusNew.run(lineIndex = i, partial=True)
            newSp = arcturusNew.Spectra[-1].resample(R=65000)
            solGfEW = newSp.calc_EW(newSp.wl[0], oldSp.wl[-1])
            arcEWDiffSelf.append(old-arcGfEW)
            arcEWDiffOther.append(old-solGfEW)
            arcEWs.append(old)

            solGfDiff.append(solarDiff)
            solar.run(lineIndex = i, partial=True)
            oldSp = solar.Spectra[-1].resample(R=100000)
            old = oldSp.calc_EW(oldSp.wl[0], oldSp.wl[-1])
            solarNew.lineList = solarGfs
            solarNew.run(lineIndex = i, partial=True)
            newSp = solarNew.Spectra[-1].resample(R=100000)
            solGfEW = newSp.calc_EW(newSp.wl[0], newSp.wl[-1])
            solarNew.lineList = arcturusGfs
            solarNew.run(lineIndex = i, partial=True)
            newSp = solarNew.Spectra[-1].resample(R=100000)
            arcGfEW = newSp.calc_EW(newSp.wl[0], newSp.wl[-1])
            solEWDiffSelf.append(old-solGfEW)
            solEWDiffOther.append(old-arcGfEW)
            solEWs.append(old)

    del(solar)
    del(arcturus)
    del(solarNew)
    del(arcturusNew)

solGfDiff = numpy.array(solGfDiff)
arcGfDiff = numpy.array(arcGfDiff)
arcEWs = numpy.array(arcEWs)
solEWs = numpy.array(solEWs)
arcEWDiffSelf = numpy.array(arcEWDiffSelf)
arcEWDiffOther = numpy.array(arcEWDiffOther)
solEWDiffSelf = numpy.array(solEWDiffSelf)
solEWDiffOther = numpy.array(solEWDiffOther)

ax.scatter(solEWs, solEWs-solEWDiffSelf, color = 'b', label = 'Sun, Solar gf\'s', marker='o', s=40.0)
ax.scatter(solEWs, solEWs-solEWDiffOther, color = 'b', label = 'Sun, Arcturus gf\'s', marker='s', s=40.0)
ax.scatter(arcEWs, arcEWs-arcEWDiffSelf, color = 'r', label = 'Arcturus, Arcturus gf\'s', marker='o', s=40.0)
ax.scatter(arcEWs, arcEWs-arcEWDiffOther, color = 'r', label = 'Arcturus, Solar gf\'s', marker='s', s=40.0)

for sgf, agf, solSelf, solOther, arcSelf, arcOther in zip(solEWs, arcEWs,
        solEWDiffSelf, solEWDiffOther, arcEWDiffSelf, arcEWDiffOther):
    ax.plot([sgf, sgf, agf, agf], [sgf-solSelf, sgf-solOther, 
        agf-arcSelf, agf-arcOther], color = 'k')
#for sg, ag, sew, aew in zip(solGfDiff, arcGfDiff, solEWDiff, arcEWDiff):
#    ax.plot([sg, ag], [sew, aew], color = 'k')

ax.legend()
fig.show()
