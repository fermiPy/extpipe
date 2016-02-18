#/bin/python

class PhysicalParameters (object):
    """This class stores the physical properties of the simulations in 
        a particular image directory. It has the following properties

    Properties:
        total_photons: The number of (scaled and added) photons generated to make this image
        low_energy: The lowest generation energy for the image (MeV)
        high_energy: The highest generation energy for the image (MeV)
        redshift: The redshift of the source
        strength: The strength of the B field (Gauss)
        power_law_index: The power law index of the source 
        boost: The lorentz boost of the source
        correlation_length: The correlation length of the B field (Mpc)
        image_index: The index of the image, in the simulation set
        angle: The observation angle of this image
        """
    def __init__(self):
        """ Sets everything to blank values"""
        self.total_photons = 0
        self.low_energy = 0
        self.high_energy = 0
        self.redshift = 0
        self.strength = 0
        self.power_law_index = 0
        self.boost = 0
        self.correlation_length = 0
        self._image_index = -1
        self._angle = -1

    @property
    def image_index(self):
        """The index must be set before reading, and setting it
            also sets the angle"""
        assert self._image_index >= 0
        return self._image_index

    @image_index.setter
    def image_index(self, index):
        self._image_index = index
        # Hardcoded for the moment
        if index == 0:
            self._angle = 0
        elif index == 1:
            self._angle = 2
        elif index == 2:
            self._angle = 5

    @property
    def angle(self):
        """The angle is set by setting the image index"""
        assert self._image_index >= 0
        return self._angle
    
    def image_directory_name(self, format='zbLi'):
        """ Create the image directory name
        
        The run format codes are:
            z = redshift
            b = magnetic field strength
            l = magnetic field correlation length
            G = AGN doppler boost
            L = power law index
            i = image index

        Args
            format (string): A combination of format codes to determine the
                string which makes up the image directory

        Return:
            string: Image directory name
            """

        name = 'image'
        code = format

        while len(code) > 0:
            if code[0] == 'z':
                name += ('-z%g'%self.redshift).strip("0")
            elif code[0] == 'b':
                name += '-b%s'%self.strength
            elif code[0] == 'l':
                name += '-c%g'%self.correlation_length
            elif code[0] == 'G':
                name += '-g%g'%self.boost
            elif code[0] == 'L':
                name += '-p%g'%self.power_law_index
            elif code[0] == 'i':
                name += '-%g'%self.image_index
            
            code = code[1:]
                
        return name



