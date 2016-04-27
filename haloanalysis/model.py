import numpy as np

from astropy.table import Table, Column

class Axis(object):

    def __init__(self,name,edges):
        self._edges = edges
        self._name = name

    @property
    def name(self):
        return self._name

    @property
    def edges(self):
        return self._edges
        
    @property
    def nbin(self):
        return len(self._edges)-1


class SpectralModel(object):

    def __init__(self,norm,index):

        self._norm = norm
        self._index = index

    @property
    def norm(self):
        return self._norm

    @property
    def index(self):
        return self._index

    def eflux(self):
        pass
    
class HaloModelCalc(object):

    def __init__(self,mmap,smodel):

        self._mmap = mmap
        self._smodel = smodel

    def eflux(self,eobs,p0,p1):
        """Evaluate the model energy flux at observed energy x for parameters
        p.

        Parameters
        ----------

        eobs : array
           Min/max energy of bin edges in 2 x N.

        p0 : array
           Model parameters of primary spectrum.

        p1 : array
           Model parameters of pair cascade.

        """

        x = np.linspace(3.0,7.0,21)
        etrue = np.vstack((x[:-1],x[1:]))
        
        w = smodel.eflux(etrue,p0)/1E-6

        # Model Map
        # Axis 0 : True Energy
        # Axis 1 : Observed Energy
        # Axis 2 : B_IGMF
        # Axis N ...

        # Interpolate model flux onto grid of true, observed energy
        mmap = interpolate(self._mmap,etrue,eobs,p1)

        detrue = x
        deobs = x
        
        # Summation over true energy
        eflux = np.sum(w*mmap*deobs*detrue,axis=0)

        return eflux

    def th68(self,x,p):
        """Evaluate the 68% containment angle at observed energy x for
        parameters p.

        Parameters
        ----------

        x : array
           Observed energy.

        p : array
           Model parameters.

        """
        pass
        
class HaloModelMap(object):

    def __init__(self,axes):

        data_shape = []
        for axis in axes:
            data_shape += [axis.nbin]        
        self._eflux = np.zeros(data_shape)
        self._th68 = np.zeros(data_shape)
        self._axes = axes

    @property
    def eflux(self):
        return self._eflux
        
    def write_fits(self,filename):

        cols = [Column(name='eflux',dtype='f8',shape=self.eflux.shape,data=self.eflux)]
        for axis in axes:
            cols += [Column(name=axis.name,dtype='f8',shape=axis.nbins,data=axis.edges)]
        tab = Table(cols)        
        tab.write(filename,format='fits',overwrite=True)

    @staticmethod
    def create_from_fits(filename):
        pass
