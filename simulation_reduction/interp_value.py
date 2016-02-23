
#Imports
from argparse import ArgumentParser
import numpy as np
import sys
from scipy.interpolate import RegularGridInterpolator as Interp


class DataStore(object):
    def __init__(self, filename):
        """ Open a numpy array in the appropriate format and read it into this storage class

        Args:
            filename (string): Filename of the data in numpy structured array format

        """
        data = np.load(filename)
        assert data.dtype == np.dtype([('redshift', '<f8'), ('strength', '<f8'), ('power_law_index', '<f8'), ('boost', '<f8'), ('correlation_length', '<f8'), 
                ('angle', '<f8'), ('low_energy', '<f8'), ('high_energy', '<f8'), ('total_photons', '<f8'), ('energy', '<f8'), ('10% Containment', '<f8'), 
                ('20% Containment', '<f8'), ('30% Containment', '<f8'), ('40% Containment', '<f8'), 
                ('50% Containment', '<f8'), ('60% Containment', '<f8'), ('70% Containment', '<f8'), 
                ('80% Containment', '<f8'), ('90% Containment', '<f8'), ('95% Containment', '<f8')])
        num_containments = 10

        # Since we collapse dimensions down below, we can't trust the interpolator to do bounds checking for us
        self.redshift_bounds = (np.min(data['redshift']), np.max(data['redshift']))    
        self.strength_bounds = (np.min(data['strength']), np.max(data['strength']))
        self.power_law_index_bounds = (np.min(data['power_law_index']), np.max(data['power_law_index']))
        self.boost_bounds = (np.min(data['boost']), np.max(data['boost']))
        self.correlation_length_bounds = (np.min(data['correlation_length']), np.max(data['correlation_length']))
        self.angle_bounds = (np.min(data['angle']), np.max(data['angle']))
        
        # Build the uniques array
        uniques = np.array([len(np.unique(data['redshift'])),
                           len(np.unique(data['strength'])),
                           len(np.unique(data['power_law_index'])),
                           len(np.unique(data['boost'])),
                           len(np.unique(data['correlation_length'])),
                           len(np.unique(data['angle']))])
   
        # Convert the data to a usable format
        ndata = data.view((np.double, len(data.dtype.names)))
        ndata[:,1] = np.log(ndata[:,1])/np.log(10) # Make the strength a log interpolation
        energies = np.unique(ndata[:,9])
        
        # Very important! Active dimensions are those with more than one grid point
        # This is needed because RegularGridInterpolator has a bug when dimensions with
        # one point are put in
        self.active_dimensions = np.array([len(np.unique(x))>1 for x in [ndata[:,i] for i in xrange(0,6)]])
        
        # Build the interpolation function table, indexed by energy
        self.interp_table = dict()
        for energy in energies:
            edata = ndata[np.abs(ndata[:,9]-energy)<.01,:]
            e_data_set = {}
            for row in edata:
                point = tuple(row[self.active_dimensions])
                containments = row[num_containments:]
                e_data_set[point] = containments
            
            points = np.array([np.unique(x) for x in [ndata[:,i] for i in xrange(0,6)]])
            points = points[self.active_dimensions]
            points = points.T
            containment_fcns = []
            for i in xrange(len(e_data_set[e_data_set.keys()[0]])):
                value_array = self._build_value_array(points, e_data_set, i)
                containment_fcns.append(Interp(points, value_array, fill_value=None))
            self.interp_table[energy] = containment_fcns                

    
    def _build_value_array(self, points, data_set, data_set_index, value_array = None, point_index = 0, cur_point=None, head=True):
        """ The RegularGridInterpolator requires an odd input format which does not work
        transparently for an unknown number of dimensions. Here we take the known data
        and the axis points (which will be later passed into the interpolator) and 
        recursively produce the correct value array.
        
        Args:
            points (numpy array): The same points array that will go into the interpolator
            data_set (dict of tuple:numpy array): We need a way to store all the 
                data information before putting it into the value array. This is
                that way
            data_set_index (int): A convenience index, since we store many different
                containment radii in the same data_set. This is the one we're using
                to make this value array
            value_array (list): Starts out blank, and is recursively filled
            point_index (int): Which point axis are we iterating over
            cur_point (numpy array): The array contains the current point, built up
                by iterating over the axes in the points array
            head (bool): Just letting us know whether to reshape and return the 
                value_array
                
        Return:
            numpy array: The value array should be exactly what the RegularGridInterpolator
                needs
                """
        if cur_point is None:
            cur_point = np.zeros(len(points))
            
        if value_array is None:
            value_array = []
        
        if point_index >= len(points):
            # Finish
            value_array.append(data_set[tuple(cur_point)][data_set_index])
            return
            
        for coord in points[point_index]:
            cur_point[point_index] = coord
            
            self._build_value_array(points, data_set, data_set_index, value_array, point_index+1, cur_point, False)
            
        if head:
            shape = np.array([len(x) for x in points])
            return np.reshape(np.array(value_array), shape)
    
    def _check_bound(self, point, bound):
        """ Check a single bound
        
        Args:
            point (Point): The Point class to check
            bound ((float, float)): The lower and upper bound to check
            
            """
        if point < bound[0] or point > bound[1]:
            raise OutOfDefinedRegionException(point, bound[0], bound[1])

    def _check_bounds(self, point):
        """ Since the interpolation function isn't going to do it for us after
        we collapsed the dimensions, we have to check our own bounds.
        
        Args:
            point (Point): The Point class containing the desired interpolation point
            """
        self._check_bound(point.redshift, self.redshift_bounds)
        self._check_bound(point.strength, self.strength_bounds)
        self._check_bound(point.power_law_index, self.power_law_index_bounds)
        self._check_bound(point.boost, self.boost_bounds)
        self._check_bound(point.correlation_length, self.correlation_length_bounds)
        self._check_bound(point.angle, self.angle_bounds)
        
    def _collapse(self, point):
        """Changes a point class into numpy array representation used for the internal
        interpolation functions, and also slices down to the active dimensions
        
        Args:
            point (Point): A Point class filled with the desired interpolation point
            
        Return:
            numpy array: An array of the points that will later be used in the interpolation
            """
        array = np.array([point.redshift, point.strength, point.power_law_index,
                    point.boost, point.correlation_length, point.angle])
        array[1] = np.log(array[1])/np.log(10)

        return array[self.active_dimensions]
        
    def _make_row(self, energy, interps):
        """Takes the results of an interpolation, along with an energy, and makes
        a standard structured array row out of them.
        
        The interpolation format should be:
            Radius at 0.1 degree containment, @0.2, @0.3, @0.4, @0.5, @0.6,
            @0.7, @0.8, @0.9, @0.95
            
        Args: 
            energy (float): The energy of this row
            interps (numpy array): A numpy array containing the results of an interpolation 
                row run
                
        Return:
            numpy array: Structured array row containing the input information
            """
        # Order is redshift, strength, power_law_index, boost,
    # correlation_length, angle, low_energy, high_energy, total_photons,
    # energy, rad at 0.1 containment, @0.2, @0.3, @0.4, @0.5, @0.6,
    # @0.7, @0.8, @0.9, @0.95
        row = np.array([(energy, interps[0], interps[1], interps[2], interps[3], interps[4],
                            interps[5], interps[6], interps[7], interps[8], interps[9])],
                        dtype=[('energy', 'd'),('10% Containment','d'),
                                ('20% Containment','d'),('30% Containment','d'),('40% Containment','d'),
                                ('50% Containment','d'),('60% Containment','d'),('70% Containment','d'),
                                ('80% Containment','d'),('90% Containment','d'),('95% Containment','d'),
                                ])
        return row
        
    def __call__(self, point):
        """ Uses the given point to build a table of interpolation values from the
        internal data.
        
        Args:
            point (Point): The point for which to do an interpolation
            
        Return:
            numpy array: This is a series of structured array rows, one for each
                energy bin in the relevant input data, containing the radius data
                for the given point. This data is interpolated using the internal
                data storage.
            """        
        self._check_bounds(point)
            
        point_row = self._collapse(point)
        outrows = []
        for energy, interp_set in self.interp_table.iteritems():
            interp_values = []
            for i, interp in enumerate(interp_set):
                try:
                    value = interp(point_row)
                    interp_values.append(value)
                except ValueError:
                    print "We're working with energy (",energy,") and interp (",i,") when a failure occurred"
                    print "The relevant point values were", point_row
                    print 'The interp grid is:', interp.grid
                    print 'The interp values are:', interp.values
                    raise
            outrows.append(self._make_row(energy, interp_values))
            
        return np.concatenate(outrows)
        
        
class OutOfDefinedRegionException(Exception):
    """A response to trying to interpolate outside the bounds."""
    def __init__ (self, value, low, high):
        """Stores the value, low bound, and high bound for the string's use
        
        Args:
            value (float): The value which was out of bounds
            low (float): The low bound
            high (float): The high bound
            
        """
        self.value = value
        self.high = high
        self.low = low
        
    def __str__ (self):
        """ Returns a string for printing the exception"""
        return "The value %s is outside the bounds (%s, %s)"%(self.value, self.low, self.high)

class Point(object):
    """A holder object for requested interpolation data"""
    def __init__ (self, redshift, strength, power_law_index, boost, correlation_length, angle):
        self.redshift = redshift
        self.strength = strength
        self.power_law_index = power_law_index
        self.boost = boost
        self.correlation_length = correlation_length
        self.angle = angle


# Arguments and parsing
def parse_options(input=None):
    """Parse the options for this script.
    
    Args:
        input (string): If you want to replace the command line input with your own version.
        
    Return:
        Options: Options object with the relevant options present
        """
    parser = ArgumentParser(description='A file to run many instances of the simulation for the given parameters')
    # General Arguments



    # Simulation parameters
    parser.add_argument('--redshift', type=float, default=0.1, help='Desired redshift')
    parser.add_argument('--strength', type=float, default=1e-16, help='Desired magnetic field strength in gauss')
    parser.add_argument('--power-law-index', type=float, default=1.8, help='Desired power law index for the injection spectrum')
    parser.add_argument('--boost', type=float, default=10, help='Desired doppler boost for the source')
    parser.add_argument('--correlation-length', type=float, default=1, help='Desired magnetic field correlation length in Mpc')
    parser.add_argument('--angle', type=float, default=2, help='Desired observation angle in degrees')

    options = parser.parse_args(input)
    return options



if __name__ == '__main__':
    options = parse_options()
    p = Point(options.redshift, options.strength, options.power_law_index, options.boost, options.correlation_length, options.angle)

    ds = DataStore('image_reduction.npy')
    print ds(p)
