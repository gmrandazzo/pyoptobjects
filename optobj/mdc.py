"""
Most Descriptor Compounds: Compound Selection algorithm
for databases.
"""
# Author: Giuseppe Marco Randazzo gmrandazzo@gmail.com

# License: BSD 3 clause

from scipy.spatial.distance import pdist, squareform
from numpy import zeros, delete

from time import sleep

class MDC:
    """Perform Most Descriptor Compound Selection of data

    Parameters
    ----------
    nobjects : int, optional, default: 0
        Number of object to select. 0 means an autostop
        criterion.

    metrics : string, optional, default: euclidean
        Set the metrics function to use in order to calculate the
        distance matrix. See scipy to view how many distace function
        are aivalable at
        http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    Attributes
    ----------
    dmx_ : array, shape (row_, row_)
        Distance Matrix of all x objects

    info_ : array, shape (row_,)
        Information Vector to select the mdc


    Returns
    ------
    mdcids: list
        Return the list of id selected from the algorithm.


    Notes
    -----
    See examples/plot_mdc_example.py for an example.

    References
    ----------
    Brian D. Hudson, Richard M. Hyde, Elizabeth Rahr and John Wood,
    Parameter Based Methods for Compound Selection from Chemical Databases,
    Quant. Struct. Act. Relat. j. 185-289 1996

    """

    def __init__(self, nobjects=0, metric='euclidean'):
        self.nobjects = nobjects
        self.metric = metric
        self.dmx_ = None
        self.info_ = None
        self.mdcids = []

    def select(self, x):
        """ Run the Most Descriptive Compound Selection"""
        self.mdcids = []
        self.dmx_ = squareform(pdist(x, self.metric))
        row = len(x)
        self.info_ = zeros(row)

        self._build_infovector()
        # Start MDC Selection from the info vector
        nmdc = 0
        stopcondition = True
        while stopcondition:
            dist = self.info_[0]
            mdc = 0

            # Select the MDC with the major information
            for i in range(len(self.info_)):
                if self.info_[i] > dist:
                    dist = self.info_[i]
                    mdc = i
                else:
                    continue

            nmdc += 1
            self.mdcids.append(mdc)
            self._rm_mdc_contrib()

            # Check Stop Condition
            if self.nobjects > 0:
                if nmdc < self.nobjects:
                    stopcondition = True
                else:
                    stopcondition = False
            else:
                cc = 0
                for item in self.info_:
                    if item < 1:
                        cc += 1
                    else:
                        continue

                if cc > len(self.mdcids):
                  stopcondition = False

        return self.mdcids

    def _build_infovector(self):
        """ build the information vector """
        row = len(self.dmx_)
        tmp = zeros((row, 2))
        for i in range(row):
            for j in range(row):
                tmp[j][0] = self.dmx_[i][j]
                tmp[j][1] = j
            tmp.sort(axis=0)

            # Reciprocal of the rank
            div = 2.0
            for j in range(row):
                if j == i:
                    self.info_[j] += 1
                else:
                    k = int(tmp[j][1])
                    self.info_[k] += 1/div
                    div += 1.0

    def _rm_mdc_contrib(self):
        """ remove the most descriptive compound contribution """
        mdc = self.mdcids[-1]
        row = len(self.dmx_)
        tmp = zeros((row, 2))
        rank = zeros(row)
        for j in range(row):
            tmp[j][0] = self.dmx_[mdc][j]
            tmp[j][1] = j
        tmp.sort(axis=0)

        div = 2.0
        for i in range(row):
            j = int(tmp[i][1])
            if j == mdc:
                rank[j] = 0.0
            else:
                rank[j] = 1.0 - (1.0/div)
                div += 1.0

        for i in range(row):
            self.info_[i] *= rank[i]

