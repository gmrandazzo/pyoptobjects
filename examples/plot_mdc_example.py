#/usr/bin/env python
"""
Example of Most Descriptive compound selection.

"""
import numpy as np
import matplotlib.pyplot as plt

from optobj.mdc import MDC


N = 20
np.random.seed(N)
mx = np.random.rand(N, 2)

csel = MDC()
idsel = csel.select(mx)

print("Selected %d objects in %d" % (len(idsel), float(N)))

colors = [1.0 for i in range(N)]
for i in idsel:
  colors[i] = 0.25

area = [50 for i in range(N)]

x = []
y = []
for i in range(len(mx)):
  x.append(mx[i][0])
  y.append(mx[i][1])

plt.scatter(x, y, s=area, c=colors, alpha=0.8)
plt.show()
