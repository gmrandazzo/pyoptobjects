#/usr/bin/env python
"""
Compare Most Descriptive Compounds Selections
with different metrics using Iris Dataset
"""
# Author: Giuseppe Marco Randazzo gmrandazzo@gmail.com
# License: BSD 3 clause

import os
import sys

path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if not path in sys.path:
    sys.path.insert(1, path)
del path


from sklearn import datasets
from optobj.mdc import MDC
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import math
from collections import Counter

# Calculate the Cosine Similarity between two lists
def get_cosim(list1, list2):
    count1 = Counter(list1)
    count2 = Counter(list2)
    allitems = set(count1.keys()).union(set(count2.keys()))
    vect1 = [count1[k] for k in allitems]
    vect2 = [count2[k] for k in allitems]
    dotprod = np.dot(vect1, vect2)
    magn1 = math.sqrt(np.dot(vect1, vect1))
    magn2 = math.sqrt(np.dot(vect2, vect2))
    return dotprod / (magn1 * magn2)

#import iris dataset
xdata = np.array(datasets.load_iris().data)

# Alternative you can import a random dataset
#N = 500
#np.random.seed(N)
#xdata = np.random.rand(N, 2)

allmetrics = ['braycurtis',
              'canberra',
              'chebyshev',
              'cityblock',
              'correlation',
              'cosine',
              'euclidean',
              'hamming',
              'mahalanobis',
              'matching',
              'minkowski',
              'seuclidean',
              'sokalmichener',
              'sokalsneath',
              'sqeuclidean']

selections = []
for m in allmetrics:
    print "Compute %s" % (m)
    if m == 'minkowski':
        dmx = squareform(pdist(xdata, m, 1))
    elif m == 'seuclidean':
        dmx = squareform(pdist(xdata, m, None))
    elif m == 'mahalanobis':
        dmx = squareform(pdist(xdata, m, None))
    else:
        dmx = squareform(pdist(xdata, m))

    csel = MDC(dmx, 20)
    selections.append(csel.select())


# Plot the results as heat map
cosdmx = np.identity(len(allmetrics))
d = []
for i in range(len(selections)):
    for j in range(i+1, len(selections)):
        cosim = get_cosim(selections[i], selections[j])
        cosdmx[i][j] = cosdmx[j][i] = cosim
        d.append(cosim)

fig1 = plt.figure(0)
ax = fig1.add_subplot(111)
ax.pcolor(cosdmx, cmap=plt.cm.Greens, alpha=1.0)
fig1 = plt.gcf()
fig1.set_size_inches(8, 11)
ax.set_frame_on(False)
ax.set_yticks(np.arange(cosdmx.shape[0])+0.5, minor=False)
ax.set_xticks(np.arange(cosdmx.shape[1])+0.5, minor=False)
ax.invert_yaxis()
ax.xaxis.tick_top()
ax.set_xticklabels(allmetrics, minor=False)
ax.set_yticklabels(allmetrics, minor=False)
label_font_size = 8
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]
             + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(int(label_font_size))
plt.xticks(rotation=90)

ax.grid(False)
ax = plt.gca()
for t in ax.xaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False
for t in ax.yaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False


# Build Dendogram
fig2 = plt.figure(1)
fig2.add_subplot(111)
linkage_matrix = linkage(d, 'complete')
dendrogram(linkage_matrix,
           color_threshold=1,
           labels=allmetrics,
           leaf_rotation=90,
           show_leaf_counts=True,
           distance_sort='ascending')
plt.show()

