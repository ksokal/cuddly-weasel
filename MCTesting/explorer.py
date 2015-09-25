import scipy
import numpy
import matplotlib.pyplot as pyplot

class walker( object ):
    def __init__(self, index, data):
        self.index = index
        self.data = [data]

    def newStep(self, newData):
        self.data.append(newData)

    def __eq__(self, index):
        return index == self.index

df = open('first.dat', 'r')

walkers = []

for line in df.readlines():
    dat = line.split()
    index = int(dat[0])
    if not(index in walkers):
        print index
        walkers.append(walker(index, numpy.array(dat[1:], dtype=numpy.float)))
    else:
        walkers[walkers.index(index)].newStep(numpy.array(dat[1:], dtype=numpy.float))
        
fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

for walker in walkers:
    ax.plot(numpy.array(walker.data)[:,3])
 
fig.show()
