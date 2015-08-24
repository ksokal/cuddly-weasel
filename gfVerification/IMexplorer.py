import pyfits
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

IM = pyfits.getdata('./Output/InteractionMatrix.fits')

for row in IM[:-2,:]:
    ax.plot(row)

fig.show()
