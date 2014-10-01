#/usr/bin/env python

"""
Example of Most Descriptive compound selection.


Code Source: Giuseppe Marco Randazzo
License: BSD 3 clausole

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

from optobj.mdc import MDC


N = 20
np.random.seed(N)
mx = np.random.rand(N, 2)

dmx = squareform(pdist(mx, 'euclidean'))

csel = MDC(dmx)
idsel = csel.select()

print("Selected %d objects in %d" % (len(idsel), float(N)))

print(idsel)

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
