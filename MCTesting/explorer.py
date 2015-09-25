import scipy
import numpy
import matplotlib.pyplot as pyplot
import sys

class walker( object ):
    def __init__(self, index, data):
        self.index = index
        self.data = [data]

    def newStep(self, newData):
        self.data.append(newData)

    def __eq__(self, index):
        return index == self.index

num = numpy.int(sys.argv[1])
df = open('first.dat', 'r')
#df = open('second.dat', 'r')
#df = open('final.dat', 'r')

walkers = []

for line in df.readlines():
    dat = line.split()
    index = int(dat[0])
    if not(index in walkers):
        walkers.append(walker(index, numpy.array(dat[1:], dtype=numpy.float)))
    else:
        walkers[walkers.index(index)].newStep(numpy.array(dat[1:], dtype=numpy.float))
        
fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

for walker in walkers:
    """
    ax.plot(numpy.array(walker.data)[:,3])
    ax.plot(numpy.array(walker.data)[:,4])
    ax.plot(numpy.array(walker.data)[:,5])
    ax.plot(numpy.array(walker.data)[:,6])
    ax.plot(numpy.array(walker.data)[:,7])
    ax.plot(numpy.array(walker.data)[:,8])
    ax.plot(numpy.array(walker.data)[:,9])
    ax.plot(numpy.array(walker.data)[:,10])
    ax.plot(numpy.array(walker.data)[:,11])
    ax.plot(numpy.array(walker.data)[:,12])
    ax.plot(numpy.array(walker.data)[:,13])
    #"""
    ax.plot(numpy.array(walker.data)[:,num])
 
fig.show()
