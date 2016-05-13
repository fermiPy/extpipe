import numpy as np

from astropy.table import Table, Column, join
from astropy.modeling.models import Gaussian1D
from astropy.convolution import convolve

from fermipy.spectrum import *
import itertools

from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d

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

    @property
    def centers(self):
        #naively averages the bin edges to determine the center
        centers=np.zeros(len(self._edges)-1)
        for i in range(len(self._edges)-1):
            centers[i]=(self._edges[i+1]+self._edges[i])/2.0
        return centers


    
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
        
        # primary spectrum start at ~100 GeV approximate as powerlaw
        #think about this implementation... WHAT is in mind for smodel     
        w = smodel.eflux(etrue,p0)/1E-6


        # Model Map
        # Axis 0 : True Energy
        # Axis 1 : Observed Energy
        # Axis 2 : B_IGMF
        # Axis N ...

        # Interpolate model flux onto grid of true, observed energy
        #data = np.meshgrid(etrue, w, indexing='ij',sparse=True)
        #define what exactly we want to do here. 
        #I interpret this is we have the etrue and eflux grid points. We want
        #the eflux if given an eobs grid point. Output eflux, input point eobs
        #mmap = interpolate(self._mmap,etrue,w) 
        #mmap=RegularGridInterpolator((etrue, w), data)
        interp=interp1d(x,w)
        mmap=interp(eobs)

        print interp, mmap
        
        #bin width in energy 
        detrue = (x[1]-x[0])
        deobs =  (x[1]-x[0])

        print detrue, deobs
        
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

    def get_eflux(self, axes): 
        #axes right now are E, and B
        #some random constants that we can change
        B0=1E-16
        gamma=0.5

        #right now the definition dependent on the E and B fields
        #eflux = PowerLaw.eval_dfde(np.array(x),[1E-10,-2.0],1E3)#*(axes[1].edges/B0)**gamma
        allaxes=[]
        for axis in axes:
            allaxes.append(axis.centers)

        params = np.meshgrid(*allaxes,indexing='ij')
        #the application does make it dependend on the axes...
        ef=PowerLaw.eval_dfde(10**params[0],[1E-10,-2.0],1E3)*(10**params[1]/B0)**gamma
        #print ef
        #print self.eflux
        self._eflux=ef
        #print self.eflux

        pass
 
        
    def write_fits(self,filename,axes):

        cols = [Column(name='eflux',dtype='f8',shape=self.eflux.shape,data=self.eflux)]
        allaxes=[]
        for axis in axes:
            allaxes.append(axis.centers)

        params = np.meshgrid(*allaxes,indexing='ij')
        for i in range(len(axes)):
            cols += [Column(name=axes[i].name,dtype='f8',data=params[i])]
            #print cols
        tab = Table(cols)        
        tab.write(filename,format='fits',overwrite=True)

        return tab

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

    #print joint['name','assoc'][test_mask]

    test_source=joint[test_mask]

    x = np.linspace(3.0,7.0,21)
    eobs = np.vstack((x[:-1],x[1:]))

    E=Axis('E',x)
    Ep=Axis('Ep',convolve(x, [0.2, 0.6, 0.3]))
    B=Axis('B',np.linspace(-13,-16,9))
    axes=[E,B]

    #print E, Ep, B
    test=HaloModelMap(axes)
    test.get_eflux(axes)
    cols=test.write_fits('test.fits',axes)
    
    #test_source_eflux=HaloModelCalc(mmap, test_source['spectrum_type']) 
    #test_source_eflux.eflux(eobs, p0, p1)

    #print test_source_eflux


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


