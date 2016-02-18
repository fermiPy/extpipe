from astropy.io import fits
from histogram import Histogram
from astropy.wcs import WCS
import numpy as np

class FitsParser:
    """ A fits file parser, customized for the simulation image model files.

    I tend to write my own fits file parsers for each usage, because I want 
    different things from each one. The purpose of this parser is to 
    extract data from the fits file and build a Histogram for each 
    energy bin. These can in turn be used to extract containment information.
    """
    def __init__(self, filename):
        """ Opens the filename and extracts relevant data

        Args:
            filename (string): The fits filename
            """
        self._filename = filename
        self._hdulist = fits.open(filename)

        self._data_header = self._hdulist[0].header
        self._data = self._hdulist[0].data
        self._energies_header = self._hdulist[1].header
        self._energies = self._hdulist[1].data

        self._parse_data_header()
        self._zero_data()

        self._totals = np.full(len(self._energies), -1)

    def _parse_data_header(self):
        """ Extract needed information from the fits header"""
        assert self._data_header['NAXIS'] == 3
        self._x_max = self._data_header['NAXIS1']
        self._y_max = self._data_header['NAXIS2']
        self._z_max = self._data_header['NAXIS3']
        self._rad_max = np.sqrt(self._x_max*self._x_max + self._y_max*self._y_max)

        self._x_ref = 200.5 #self._data_header['CRPIX1'] # The real center point is not what is in the header
        assert self._data_header['CRVAL1'] == 0, self._data_header['CRVAL1']
        self._x_del = self._data_header['CDELT1']

        self._y_ref = 200.5 #self._data_header['CRPIX2']
        assert self._data_header['CRVAL2'] == 0, self._data_header['CRVAL2']
        self._y_del = self._data_header['CDELT2']

    def _zero_data(self):
        """ All data has a small floor, set to several orders of magnitude below
        the lowest actual data value. This is because gtobssim cannot deal with
        0 value in its model file parser. So, here we remove that false floor"""
        min_data = np.min(self._data)
        np.place(self._data,self._data == min_data, [0])
        
    def _get_edges(self, ref, delta):
        """ For building a histogram 
            
        Args:
            ref (float): reference pixel
            delta (float): spacing between fits points
            
        Return:
            numpy array: Length = number of pixels + 1, this constitutes the left and right
                edges of a Histogram along an axis
        """
        return (np.arange(0.5,401.5) - ref)*delta

    def hist_set(self):
        """ A method to generate and yield (energy, histogram) pairs.

        Usage: for energy, hist in fits_parser.hist_set():
            ... do something ...

        yield:
            (float, Histogram): The energy and histogram pairs. This will yield one set
                for each energy in the fits file.
                """
        print 'Generating hist set'
        xedges = self._get_edges(self._x_ref, self._x_del)
        yedges = self._get_edges(self._y_ref, self._y_del)

        for i, energy in enumerate(self._energies):
            print '\tHistogram for', energy
            h = Histogram()
            h.setup(self._data[i], xedges, yedges)
            
            energy = energy[0]

            yield (energy, h)

        print 'Done'


def test_edges():
    """ Test function to make sure the edges are correct"""
    fp = FitsParser('images/image-z0.1_b1e-15_p1.8-0/image_combined.halo.fits')

    gold = np.array([-1.,    -0.995, -0.99,  -0.985, -0.98,  -0.975, -0.97,  -0.965, -0.96,  -0.955,
 -0.95,  -0.945, -0.94,  -0.935, -0.93,  -0.925, -0.92,  -0.915, -0.91,  -0.905, -0.9,
 -0.895, -0.89,  -0.885, -0.88,  -0.875, -0.87,  -0.865, -0.86,  -0.855, -0.85,
 -0.845, -0.84,  -0.835, -0.83,  -0.825, -0.82,  -0.815, -0.81,  -0.805, -0.8,   -0.795,
 -0.79,  -0.785, -0.78,  -0.775, -0.77,  -0.765, -0.76,  -0.755, -0.75,  -0.745,
 -0.74,  -0.735, -0.73,  -0.725, -0.72,  -0.715, -0.71,  -0.705, -0.7,   -0.695,
 -0.69,  -0.685, -0.68,  -0.675, -0.67,  -0.665, -0.66,  -0.655, -0.65,  -0.645,
 -0.64,  -0.635, -0.63,  -0.625, -0.62,  -0.615, -0.61,  -0.605, -0.6,   -0.595,
 -0.59,  -0.585, -0.58,  -0.575, -0.57,  -0.565, -0.56,  -0.555, -0.55,  -0.545,
 -0.54,  -0.535, -0.53,  -0.525, -0.52,  -0.515, -0.51,  -0.505, -0.5,   -0.495,
 -0.49,  -0.485, -0.48,  -0.475, -0.47,  -0.465, -0.46,  -0.455, -0.45,  -0.445,
 -0.44,  -0.435, -0.43,  -0.425, -0.42,  -0.415, -0.41,  -0.405, -0.4,   -0.395,
 -0.39,  -0.385, -0.38,  -0.375, -0.37,  -0.365, -0.36,  -0.355, -0.35,  -0.345,
 -0.34,  -0.335, -0.33,  -0.325, -0.32,  -0.315, -0.31,  -0.305, -0.3,   -0.295,
 -0.29,  -0.285, -0.28,  -0.275, -0.27,  -0.265, -0.26,  -0.255, -0.25,  -0.245,
 -0.24,  -0.235, -0.23,  -0.225, -0.22,  -0.215, -0.21,  -0.205, -0.2,   -0.195,
 -0.19,  -0.185, -0.18,  -0.175, -0.17,  -0.165, -0.16,  -0.155, -0.15,  -0.145,
 -0.14,  -0.135, -0.13,  -0.125, -0.12,  -0.115, -0.11,  -0.105, -0.1,   -0.095,
 -0.09,  -0.085, -0.08,  -0.075, -0.07,  -0.065, -0.06,  -0.055, -0.05,  -0.045,
 -0.04,  -0.035, -0.03,  -0.025, -0.02,  -0.015, -0.01,  -0.005,  0,     0.005,
  0.01,   0.015,  0.02,   0.025,  0.03,   0.035,  0.04,   0.045,  0.05,   0.055,
  0.06,   0.065,  0.07,   0.075,  0.08,   0.085,  0.09,   0.095,  0.1,    0.105,
  0.11,   0.115,  0.12,   0.125,  0.13,   0.135,  0.14,   0.145,  0.15,   0.155,
  0.16,   0.165,  0.17,   0.175,  0.18,   0.185,  0.19,   0.195,  0.2,    0.205,
  0.21,   0.215,  0.22,   0.225,  0.23,   0.235,  0.24,   0.245,  0.25,   0.255,
  0.26,   0.265,  0.27,   0.275,  0.28,   0.285,  0.29,   0.295,  0.3,    0.305,
  0.31,   0.315,  0.32,   0.325,  0.33,   0.335,  0.34,   0.345,  0.35,   0.355,
  0.36,   0.365,  0.37,   0.375,  0.38,   0.385,  0.39,   0.395,  0.4,    0.405,
  0.41,   0.415,  0.42,   0.425,  0.43,   0.435,  0.44,   0.445,  0.45,   0.455,
  0.46,   0.465,  0.47,   0.475,  0.48,   0.485,  0.49,   0.495,  0.5,    0.505,
  0.51,   0.515,  0.52,   0.525,  0.53,   0.535,  0.54,   0.545,  0.55,   0.555,
  0.56,   0.565,  0.57,   0.575,  0.58,   0.585,  0.59,   0.595,  0.6,    0.605,
  0.61,   0.615,  0.62,   0.625,  0.63,   0.635,  0.64,   0.645,  0.65,   0.655,
  0.66,   0.665,  0.67,   0.675,  0.68,   0.685,  0.69,   0.695,  0.7,    0.705,
  0.71,   0.715,  0.72,   0.725,  0.73,   0.735,  0.74,   0.745,  0.75,   0.755,
  0.76,   0.765,  0.77,   0.775,  0.78,   0.785,  0.79,   0.795,  0.8,    0.805,
  0.81,   0.815,  0.82,   0.825,  0.83,   0.835,  0.84,   0.845,  0.85,   0.855,
  0.86,   0.865,  0.87,   0.875,  0.88,   0.885,  0.89,   0.895,  0.9,    0.905,
  0.91,   0.915,  0.92,   0.925,  0.93,   0.935,  0.94,   0.945,  0.95,   0.955,
  0.96,   0.965,  0.97,   0.975,  0.98,   0.985,  0.99,   0.995,  1.   ] )
    trial = fp._get_edges(200.5, 0.005)
    assert np.all(trial - gold < 0.001)

if __name__ == '__main__':
    test_hist_set()

