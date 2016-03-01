"""
Kennard Stones object selection
"""
# Author: Giuseppe Marco Randazzo gmrandazzo@gmail.com
# License: BSD 3 clause

from random import randrange
import time

from numpy import zeros, array

class KS(object):
    """Perform Kenard-Stoness compound object selection

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
        try:
            self.dmx_ = dmx.tolist() #convert to list to be faster
        except AttributeError:
            self.dmx_ = dmx
        self.nobjects = nobjects
        self.info_ = None
        self._build_infovector()
        self.ksids = []


    def kslist(self):
        """ Return the list of Kennard-Stones selected objects """
        return self.ksids


    def getnext(self):
        """ Get the next compound """
        self._appendnext()
        return self.ksids[-1]


    def select(self):
        """ Run the Kennard-Stones compound Selection """
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



function [model,test]=kenstone(X,k)

[m,n]=size(X);
if k>=m | k<=0
    h=errordlg('Wrongly specified number of objects to be selected to model set.','Error');
    model=[];
    if nargout==2
        test=[];
    end
    waitfor(h)
    return
end

x=[[1:size(X,1)]' X];
n=size(x,2);
[i1,ind1]=min(fastdist(mean(x(:,2:n)),x(:,2:n)));
model(1)=x(ind1,1);
x(ind1,:)=[];

[i2,ind2]=max(fastdist(X(model(1),:),x(:,2:n)));
model(2)=x(ind2,1);
x(ind2,:)=[];

for d=3:k
    [ii,ww]=max(min(fastdist(x(:,2:n),X(model,:))));
	model(d)=x(ww,1);
	x(ww,:)=[];
end

if nargout==2
    test=1:size(X,1);
    test(model)=[];
end


% --->

function D=fastdist(x,y)

% Calculates squared Euclidean distances between two sets of objetcs

D=((sum(y'.^2))'*ones(1,size(x,1)))+(ones(size(y,1),1)*(sum(x'.^2)))-2*(y*x');
