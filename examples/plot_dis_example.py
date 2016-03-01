#/usr/bin/env python

"""
Example of Most Descriptive compound selection.


Code Source: Giuseppe Marco Randazzo
License: BSD 3 clausole

"""
import os
import sys

path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not path in sys.path:
    sys.path.insert(1, path)
del path

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

from optobj.disc import DISC
import time

N = 2000
np.random.seed(N)
mx = np.random.rand(N, 2)

dmx = squareform(pdist(mx, 'euclidean'))

t = time.time()
csel = DISC(dmx, "min", int(0.20*N))
idsel = csel.select()
print "Time: %.3f" % (time.time()-t)

#print(idsel)
print("Selected %d objects in %d" % (len(idsel), float(N)))

#print(idsel)

colors = ["black" for i in range(N)]
for i in idsel:
  colors[i] = "red"

area = [50 for i in range(N)]

x = []
y = []
for i in range(len(mx)):
  x.append(mx[i][0])
  y.append(mx[i][1])

plt.scatter(x, y, s=area, c=colors, alpha=0.8)
plt.show()
