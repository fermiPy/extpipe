import numpy as np

from astropy.table import Table, Column, join
from astropy.modeling.powerlaws import LogParabola1D, PowerLaw1D

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

    def interpolate(self):
        #scipy.interpolate.LinearNDInterpolator, or RegularGridInterpolator
        pass

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
        
        # primary spectrum start at ~100 GeV approximate as powerlaw
        w = smodel.eflux(etrue,p0)/1E-6

        # Model Map
        # Axis 0 : True Energy
        # Axis 1 : Observed Energy
        # Axis 2 : B_IGMF
        # Axis N ...

        # Interpolate model flux onto grid of true, observed energy
        mmap = interpolate(self._mmap,etrue,eobs,p1) 

        #bin width in energy 
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
        #weight theta by energy flux, not primary spectrum as in w, it's energy flux. 
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

    #def get_eflux(self, energy, p0):
    #here is where we make the map
    #write values into eflux array
    # contains log parabola information 
    # model: LogParabola1D(10**8,1000,2,1)
        
    def write_fits(self,filename):

        cols = [Column(name='eflux',dtype='f8',shape=self.eflux.shape,data=self.eflux)]
        for axis in axes:
            cols += [Column(name=axis.name,dtype='f8',shape=axis.nbins,data=axis.edges)]
        tab = Table(cols)        
        tab.write(filename,format='fits',overwrite=True)

    @staticmethod
    def create_from_fits(filename):
        pass


if __name__ == '__main__':
    
    from astropy.io import fits 
    import stack
    
    def column(matrix, i):
        return[row[i] for row in matrix]

    # Reading in a fits file using the Table class

    hdulist=Table.read("/u/gl/mdwood/ki20/mdwood/fermi/ext_analysis/v9/table_std_all_3fgl_glat050.fits")
    hdu_lnl=Table.read("/u/gl/mdwood/ki20/mdwood/fermi/ext_analysis/v9/table_std_all_3fgl_glat050_lnl.fits")

    joint=join(hdulist,hdu_lnl)
    test_mask=stack.create_mask(joint, {'assoc':['1ES 0229+200']})

    print joint['name','assoc'][test_mask]

    test_source=joint[test_mask]

    # Doing this using astropy.fits just for fun
    #hdulist=fits.open("/u/gl/mdwood/ki20/mdwood/fermi/ext_analysis/v9/table_std_all_3fgl_glat050.fits")
    #hdu_lnl=fits.open("/u/gl/mdwood/ki20/mdwood/fermi/ext_analysis/v9/table_std_all_3fgl_glat050_lnl.fits")
    #data=hdulist[1].data
    #lnl=hdu_lnl[1].data

    #for i in range(len(data)):
    #    if data.field('assoc')[i]=='1ES 0229+200':
    #        print data.field('name')[i], lnl.field('name')[i]
    #        eflux_ul95=data.field('fit_halo_scan_eflux_ul95')[i]
    #        src_eflux=lnl.field('fit_src_sed_scan_eflux')[i]
    #        halo_eflux=lnl.field('fit_halo_sed_scan_eflux')[i]
    #        print eflux_ul95, src_eflux


