import numpy as np
import numpy.random as rnd
import sys
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from pixel import Pixel
import cProfile
import time
import pickle

class Histogram(object):
    """ This is a wrapper class around the data from a fits image. It is
    used to get precise containment(radius) and radius(containment) values"""
    def __init__(self, x_values=None, y_values=None):
        """ One of two ways to build the Histogram, this takes a list of x and y values

        The bulk of the work is done by setup

        Args:
            x_values (numpy array, optional): If this is not present, an empty Histogram is constructed.
                Otherwise, it is used to build the internal 2d histogram.
            y_values (numpy array, optional): If this is not present, an empty Histogram is constructed.
                Otherwise, it is used to build the internal 2d histogram.
        """

        if x_values is None or y_values is None:
            return
            
        self._x_bin_max = 400
        self._y_bin_max = 400
        hist, xedges, yedges = np.histogram2d(x_values,  y_values, 
                                bins=[self._x_bin_max, self._y_bin_max], range=[[-2,2],[-2,2]])

        self.setup(hist, xedges, yedges)

    def setup(self, hist, xedges, yedges):
        """ If an empty histogram was constructed from __init__, this must be called
        to fill it. Takes equivalent output as from the numpy function histogram2d, but
        those values need not be constructed using that function.

        Args:
            hist (numpy array): The 2d histogram of data
            xedges (numpy array): The xedges, constituting a nbins+1 array with left and right
                edges for each x-axis bin
            yedges (numpy array): As xedges, but for y.
        """    
            
        self._hist = hist
        self._xedges = xedges
        self._yedges = yedges

        self._total = np.sum(self._hist)
        self._x_del = self._xedges[1] - self._xedges[0]
        self._y_del = self._yedges[1] - self._yedges[0]
        self._search_radius = np.sqrt(self._x_del**2 + self._y_del**2) * 2
        
        self._pixels = []
        self._middle_rad = []
        for i,x in enumerate(self._xedges[:-1]):
            for j,y in enumerate(self._yedges[:-1]):
                v = self._hist[i][j]
                self._pixels.append((Pixel(x,y,self._x_del, self._y_del, v), np.sqrt((x+self._x_del/2.)**2 + (y+self._y_del/2.)**2), v))
                
        self._containment_cache = {}
        self._min_containment = None
        self._max_containment = None
        self._build_containment_cache()


    def _value(self,pixel_set, radius):
        """ A shorthand function, so we don't actually call pixels that are either fully below the radius
            or fully above the radius. If above, 0 is returned. If below, the number of photons in this
            pixel are returned. If in between, we call the pixel to get the value.

            Args:
                pixel_set (list): The first entry is the pixel, the second is the min radius for the pixel,
                    the third is the max radius for the pixel.
                radius (float): The radius at which to search

            Return:
                float: The value corresponding to the weghted number of photons contained in the radius
                    for this pixel.
                    """

        if pixel_set[1] > radius + self._search_radius:
            return 0.
        
        if pixel_set[1] < radius - self._search_radius:
            return pixel_set[2]
        
        return pixel_set[0].value(radius)    
        

    #@do_profile(['Histogram._get_integral'])        
    def _get_integral(self,radius):
        """ Integrates over all the pixels from 0 to the specified radius

        Args:
            radius (float): The max radius to integrate up to

        Return:
            float: The integral of all pixels up to radius
            """
        integral = 0
        #print 'Starting integral'
#        start = time.time()
        for i, pixel in enumerate(self._pixels):
#            if i % 2000 == 0:
#                sys.stdout.write('.')
#                if i > 0 and i % 100000 == 0:
#                    sys.stdout.write(' = %d of %d (%f)\n'%(i,len(self._pixels), time.time()-start))
#                sys.stdout.flush()
            integral += self._value(pixel, radius)
#        if i % 100000 != 0:
#            sys.stdout.write(' = %d of %d (%f)\n'%(i,len(self._pixels), time.time()-start))
#        print
        return integral
    

    def _containment_single(self, radius):
        """ Run a single radius through the containment machinery. Caches values we've seen before to
            (greatly) speed up computation time.

            Args:
                radius (float): Radius for which to determine containment

            Return:
                float: The containment fraction for this radius
                """
        if radius in self._containment_cache:
            return self._containment_cache[radius]
        else:
            int = self._get_integral(radius)
            containment = int/self._total
            self._containment_cache[radius] = containment
            
            return containment
        
    def containment(self, radius):
        """ A vectorized version of the containment search.

        Args:
            radius (numpy array): A 1d array of floats to get the containment for

        Return:
            numpy array: A 1d array of containment fractions corresponding to the radius array
            """
        f_vec = np.vectorize(self._containment_single)
        return f_vec(radius)  
        
    def _step_containment(self, start, stop, check_value, step, done):
        """Identify and cache the containments for a series of radii to build a
            cdf curve. 

        Args:
            start (float): Radius at which to start the stepping
            stop (float): Radius to stop at
            check_value (float): Part of the short-circut machinery. Stops execution of we've gone
                higher than the requested containment value
            step (float): The step value for the loop
            done (function): Part of the short-circuit machinery. Takes two floats, the current 
                containment and the desired, and should return true if we can stop executing.
        
        Return:
            (float, float): The last radius checked, and the last containment value reached

            """

        for r in np.arange(start,stop,step):
            contain = self.containment(r)
            if done(contain, check_value):
                return (r,contain)
        return (stop-step, self.containment(stop-step))
        
    def _setup_containment_fcn(self):
        """ Sets up the containment interpolation curve from the containment_cache."""
        containment_curve = []
        for key,value in self._containment_cache.iteritems():
            containment_curve.append((key,value))
            
        containment_curve.sort(key=lambda x:x[0])
        y,x = zip(*containment_curve)
        
        self._contain_fcn = self._containment_fcn = interp.interp1d(x,y)
            
        
    def _build_containment_cache(self):
        """ Creates the initial containment cache values from containment = (0.34, 0.95) 
            or radius=(0, 2.1), whichever is smaller."""
        self._max_containment = self._step_containment(.6, 2.1,0.95,.1, done=lambda x,y: x>y)
        self._min_containment = self._step_containment(.6, 0,0.34,-.1, done=lambda x,y: x<y)
        self._setup_containment_fcn()
                
    def _radius_single(self,containment):
        """ Gets the radius for a single containment value. If the cache extends beyond the requested
            value, just return the interpolated containment function for the requested value. Otherwise
            build the cache out until that is true. In the case that the requested containment is 
            greater than the size of the image (true if significant photons are in the corners--because
            we stop searching a r=R_0, not x,y=R_0), return the size of the image.

        Args:
            containment (float): The containment for which to find a radius

        Return:
            float: The radius corresponding to the desired containment

            """

        Uses the containment cache, if it contains enough values to 
        if containment < self._min_containment[1]:
            self._min_containment = self._step_containment(self._min_containment[0],0,containment,-.1, done=lambda x,y: x<y)
            self._setup_containment_fcn()      
        elif containment > self._max_containment[1]:
            self._max_containment = self._step_containment(self._max_containment[0],2.1, containment,.1, done=lambda x,y: x>y)
            self._setup_containment_fcn()            
        
        try:
            radius = self._contain_fcn(containment)
        except ValueError:
            return self._max_containment[0]
        return radius
        
    def radius(self, containment):
        """ A vectorized version of the radius search.

        Args:
            containment (numpy array): A 1d array of floats to get the radius for

        Return:
            numpy array: A 1d array of radius fractions corresponding to the containment array
            """
        f_vec = np.vectorize(self._radius_single)
        return f_vec(containment)  
            
            

#def test_radius():
#    dist = Distribution()
#    h = Histogram(*dist.values(10))
    
#    assert ( h._radius(199,199) == h._radius(200,200))
#    assert ( h._radius(0,0) == h._radius(399,399))
    
            
def test_fraction_of_area():
    """ Basic pixel and short_cut funciton testing"""
    dist = Distribution()
    h = Histogram(*dist.values(10))
    
    assert(h._fraction_of_area(.1, 201,201) == 1.)
    assert(h._fraction_of_area(.01, 201,201) == 0.)
  
    r = np.arange(0.01,0.03,0.001)
    fractions = h._fraction_of_area(r,201*np.ones(len(r)),201*np.ones(len(r)))
    assert np.all(fractions - [ 0.,0.,0.,0.,0.,0.00833252,0.03529066, 0.07989413,
                    0.14155475, 0.22046858, 0.31604745, 0.42789923, 0.55641604,
                    0.69199098, 0.80109793, 0.88412901, 0.94382904, 0.98137438,
                    0.99852956, 1.] < 1e-6 )

def test_containment():
    """ Make sure the containment is right"""
    dist = Distribution()
    print dist
    
    dist._rotation = 0
    dist._std_one = .5
    dist._std_two = .5
    dist._set_containment()
    
    h = Histogram(*dist.values(10000))
    x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
    y = h.containment(x)
    
    test_y = np.array([ 0.019119  ,  0.07709066,  0.15880573,  0.26645272,  0.38526369,
        0.50733592,  0.61797834,  0.71550907,  0.79519714,  0.85742165,
        0.90597305,  0.93793006,  0.96381745,  0.97971121,  0.98867297,
        0.99379649,  0.99773116,  0.9987    ,  0.99924042])
        
    assert np.all(test_y == y)


def test_contain_radius():
    """ Make sure the radius value are right"""
    dist = Distribution()
    print dist
    
    dist._rotation = 0
    dist._std_one = .5
    dist._std_two = .5
    dist._set_containment()
    
    h = Histogram(*dist.values(1000000))
    
    radius01 = h.radius(.1)
    radius02 = h.radius(.2)
    radius03 = h.radius(.3)
    radius04 = h.radius(.4)
    radius05 = h.radius(.5)
    radius06 = h.radius(.6)
    radius07 = h.radius(.7)
    radius08 = h.radius(.8)
    radius09 = h.radius(.9)
    
    good = np.array([ 0.22684681,  0.33301551,  0.42216387,  0.50615889,
                    0.59069347,  0.68062232,  0.77999141,  0.90031203,
                    1.07865906])

    trial = np.array([radius01, radius02, radius03, radius04, radius05, 
            radius06, radius07, radius08, radius09,])
                          
    assert np.all( np.abs(good - trial) < 3e-3 ), np.abs(good - trial)
