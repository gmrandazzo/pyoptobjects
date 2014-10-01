"""
Most Descriptor Compounds: Compound Selection algorithm
for databases.
"""
# Author: Giuseppe Marco Randazzo gmrandazzo@gmail.com
# License: BSD 3 clause

from numpy import zeros, array

class MDC(object):
    """Perform Most Descriptor Compound Selection of data

    Parameters
    ----------
    dmx : array, shape(row,row)
        A square distance matrix.
        To build a distance matrix see scipy at:
        http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    nobjects : int, optional, default: 0
        Number of object to select. 0 means an autostop
        criterion.

    Attributes
    ----------
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

    def __init__(self, dmx, nobjects=0):
        self.dmx_ = dmx
        self.nobjects = nobjects
        self.info_ = None
        self._build_infovector()
        self.mdcids = []


    def mdclist(self):
        """ Return the list of most descriptor compounds """
        return self.mdcids


    def getnext(self):
        """ Get the next most descriptor compound """
        self._appendnext()
        return self.mdcids[-1]


    def select(self):
        """ Run the Most Descriptive Compound Selection"""
        stopcondition = True
        while stopcondition:
            self._appendnext()
            self._rm_mdc_contrib()
            # Check Stop Condition
            if self.nobjects > 0:
                if len(self.mdcids) == len(self.dmx_):
                    stopcondition = False
                else:
                    if len(self.mdcids) < self.nobjects:
                        continue
                    else:
                        stopcondition = False
            else:
                ncheck = 0
                for item in self.info_:
                    if item < 1:
                        ncheck += 1
                    else:
                        continue

                if ncheck > len(self.mdcids):
                    stopcondition = False

        return self.mdcids


    def _build_infovector(self):
        """ build the information vector """
        row = len(self.dmx_)
        self.info_ = zeros(row)
        tmp = zeros((row, 2))
        for i in range(row):
            for j in range(row):
                tmp[j][0] = self.dmx_[i][j]
                tmp[j][1] = j
            tmp = array(sorted(tmp, key=lambda item: item[0]))
            # Reciprocal of the rank
            div = 2.0
            for j in range(row):
                if j == i:
                    self.info_[j] += 1
                else:
                    k = int(tmp[j][1])
                    self.info_[k] += 1/div
                    div += 1.0


    def _appendnext(self):
        """ Append the next most descriptive compound to list """
        dist = self.info_[0]
        mdc = 0

        # Select the MDC with the major information
        for i in range(1, len(self.info_)):
            if self.info_[i] > dist:
                dist = self.info_[i]
                mdc = i
            else:
                continue
        self.mdcids.append(mdc)


    def _rm_mdc_contrib(self):
        """ remove the most descriptive compound contribution """
        mdc = self.mdcids[-1]
        row = len(self.dmx_)
        tmp = zeros((row, 2))
        rank = zeros(row)
        for j in range(row):
            tmp[j][0] = self.dmx_[mdc][j]
            tmp[j][1] = j
        tmp = array(sorted(tmp, key=lambda item: item[0]))
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

