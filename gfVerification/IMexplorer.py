import pyfits
import matplotlib.pyplot as pyplot
import scipy
import numpy

PlotIM = False
nFilt = 60

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

IM = pyfits.getdata('./Output/InteractionMatrix.fits')

if PlotIM:
    for row in IM:
        ax.plot(row)
else:
    U,S,V = scipy.linalg.svd(IM)

    ax.plot(numpy.log10(S))

    dims = IM.shape
    D = 1.0/S
    D = 1.0/(S[0:-nFilt])
    S[-nFilt:] = 0.0
    newS = numpy.zeros((dims[0], dims[1]))
    I = [i for i in range(dims[1])]
    for i in range(len(D)):
        newS[i][i] = D[i]

    S = newS.copy()
    CM = numpy.array(scipy.matrix(V.T.dot(S.T.dot(U.T)), dtype=numpy.float32)).T

    hdu = pyfits.PrimaryHDU(CM)
    hdu.writeto('./Output/CommandMatrix.fits', clobber=True)


fig.show()

