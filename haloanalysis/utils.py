import numpy as np
from numpy.core import defchararray
from scipy.interpolate import RegularGridInterpolator


def load_source_rows(tab, names, key='assoc'):
    """Load the rows from a table that match a source name.

    Parameters
    ----------
    tab : `astropy.table.Table`
       Table that will be searched.

    names : list
       List of source identifiers.

    key : str
       Name of the table column that will be searched for a source
       matching key.

    Returns
    -------
    outtab : `astropy.table.Table`
       Table containing the subset of rows with matching source identifiers.

    """
    names = [name.lower().replace(' ', '') for name in names]
    tab[key] = defchararray.replace(defchararray.lower(tab[key]),
                                    ' ', '')
    mask = create_mask(tab, {key: names})
    return tab[mask]


def create_mask(tab, target_def):
    """Create a table mask from a target definition."""

    m = np.empty(len(tab), dtype=bool); m.fill(True)

    for k,v in target_def.items():
        
        if isinstance(v,list):

            m0 = np.zeros(len(tab),dtype=bool)

            for t in v:
                m0 |= (tab[k] == t)

            m &= m0
            
        elif isinstance(v,dict):

            m0 = np.empty(len(tab),dtype=bool)
            m0.fill(True)
            if 'min' in v:
                m0 &= (tab[k] >= v['min'])

            if 'max' in v:
                m0 &= (tab[k] <= v['max'])
            m &= m0
            
    return m

def interp_map(z, axis0, axis1,dim=0):

    s0 = z.ndim*[None]
    s1 = z.ndim*[slice(None)]
    s0[idim] = slice(None)
    s1[idim] = slice(0,1)
    z /= axis0.width[s0]

    shape = list(z.shape)
    shape[idim] = len(axis1.centers)
    
    zinterp = np.zeros(shape)
    
    for x, idx in np.ndenumerate(z[s1]):
        zinterp = np.interp(axis1.centers,
                            axis0.centers,
                            z[:,i])
        
class MapND(object):
    """Container class representing an n-dimensional map."""
    
    def __init__(self, axes, data, log_interp=False):
        """
        Parameters
        ----------
        axes : list
           List of `Axis` objects defining the n-dimensional grid.

        data : `~numpy.ndarray`
           
        log_interp : bool
           Perform interpolation in log-space.

        """
        
        self._axes = axes
        self._data = data
        self._log_interp = log_interp
        points = [ax.centers for ax in axes]

        if log_interp:
            self._fn = RegularGridInterpolator(points, np.log(data),
                                               bounds_error=False,
                                               fill_value=None)
        else:
            self._fn = RegularGridInterpolator(points, data,
                                               bounds_error=False,
                                               fill_value=None)
            
    @property
    def fn(self):
        return self._fn

    @property
    def axes(self):
        return self._axes

    @property
    def data(self):
        return self._data

    def marginalize(self, dims):

        data = np.squeeze(np.apply_over_axes(np.sum,self.data,axes=dims))

        axes = []
        for i, axis in enumerate(self.axes):
            if i not in dims:
                axes += [axis]
                
        return MapND(axes, data, self._log_interp)
    
    def slice(self, dims, vals):
        axis_xvals = []
        axes = []
        for i, axis in enumerate(self.axes):
            axis_xvals += [axis.centers]
            if not i in dims:
                axes += [axis]

        for i, v in zip(dims,vals):
            axis_xvals[i] = np.array(v,ndmin=1)
            
        interp_xvals = np.meshgrid(*axis_xvals,indexing='ij',sparse=True)
        data = np.squeeze(self.interp(tuple(interp_xvals)))
        return MapND(axes, data, self._log_interp)
        
    def interp(self, *args):

        if self._log_interp:
            return np.exp(self._fn(*args))
        else:
            return self._fn(*args)

class Axis(object):

    def __init__(self, name, edges, centers=None):
        self._edges = edges
        self._centers = (0.5*(edges[1:] + edges[:-1])
                         if centers is None else centers)
        self._name = name

    @staticmethod
    def create_from_centers(name, centers, logscale=False):
        """Create a new axis object from a sequence of bin centers."""
        
        if logscale:
            delta = np.log(centers[1:])-np.log(centers[:-1])
        else:
            delta = centers[1:]-centers[:-1]

        if len(delta) == 0:
            delta = np.array([1.0])
        else:
            delta = np.insert(delta,0,delta[0])
            
        if logscale:
            edges_lo = np.log(centers) - 0.5*delta
            edges_hi = np.log(centers) + 0.5*delta
            edges = np.exp(np.insert(edges_hi,0,edges_lo[0]))
        else:
            edges_lo = centers - 0.5*delta
            edges_hi = centers + 0.5*delta
            edges = np.insert(edges_hi,0,edges_lo[0])
            
        return Axis(name, edges, centers)
        
    @property
    def name(self):
        return self._name

    @property
    def edges(self):
        return self._edges

    @property
    def lo(self):
        """Return the lower bin edges."""
        return self._edges[:-1]

    @property
    def hi(self):
        """Return the upper bin edges."""
        return self._edges[1:]

    @property
    def nbin(self):
        return len(self._edges)-1

    @property
    def centers(self):
        return self._centers

    @property
    def width(self):
        return self._edges[1:] - self._edges[:-1]
