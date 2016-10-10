import pylab as plt
import numpy as np
import sys

ncols  = 7

class dataprof:
    def __init__(self):
        self.z    = []
        self.p    = []
        self.t    = []
        self.tap  = []
        self.tacb = []
        self.tac  = []
        self.taud = []

indata = dataprof()
        
f = open('profile.out','r')
for ii in np.arange(2): header = f.readline()
while 1:
    line = f.readline()
    if not line: break
    cols = line.rsplit()
    for ii in np.arange(len(cols)): cols[ii] = np.float(cols[ii])

    indata.z.append(cols[0])
    indata.p.append(cols[1])
    indata.t.append(cols[2])
    indata.tap.append(cols[3])
    indata.tacb.append(cols[4])
    indata.tac.append(cols[5])
    indata.taud.append(cols[6])
    
f.close()


plt.plot(indata.t,indata.z,'o-')
plt.show()
