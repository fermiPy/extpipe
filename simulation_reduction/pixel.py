import numpy as np
import numpy.random as rnd
import sys


class Pixel(object):
    """ A single pixel, used to determine the weighted value enclosed by 
        a circle """
    def __init__(self, x_0, y_0, width, height, value):
        """ Set up the pixel

        Args:
            x_0 (float): Lower left corner's x value
            y_0 (float): Lower left corner's y value
            width (float): Width of the pixel
            height (float): Heigh of the pixel
            value (float): The value of the pixel, usually the number of photons
                it contains
        """

        self._x_0 = x_0
        self._y_0 = y_0        
        self._width = width
        self._height = height
        self._value = value
        
        self._total_area = width*height
    
    def _triangle(self, w, h):
        """ Convenience function returning the area of a triangle

        Args:
            w (float): width
            h (float): height

        Returns:
            float: Area of triangle"""
        return (w*h/2.)/self._total_area
        
    def _removed_triangle(self, w,h):
        """ Convenience function returning the area of a rectangle minus a triangle.
            Rectangle is assumed to be self.width*self.height size
            
        Args:
            w (float): width
            h (float): height

        Returns:
            float: Area of rectangle minus triangle"""
        return (self._total_area - w*h/2.)/self._total_area
    
    def _trapizoid(self, b1, b2, h):
        """ Convenience function returning the area of a trapizoid

        Args:
            b1 (float): base 1
            b2 (float): base 2
            h (float): height

        Returns:
            float: Area of trapizoid"""
        return ((b1+b2)*h / 2.)/self._total_area
                                
    def _rad_triple(self, x,y):
        """ Convenience function to make a radius triple

        Args:
            x (float): x-value
            y (float): y-value

        Return:
            (float, float, float): (x,y,radius)"""
        return (x,y,np.sqrt(x*x+y*y))

    def _get_orientation(self):
        """Returns the corners in (nearest, second nearest, second farthest, farthest)
        In the case of ties, the two in the middle will be the same.
        
        Return:
            (list): As above"""
        
        return sorted([self._rad_triple(self._x_0,self._y_0),
                        self._rad_triple(self._x_0+self._width,self._y_0),
                        self._rad_triple(self._x_0,self._y_0+self._height),
                        self._rad_triple(self._x_0+self._width,self._y_0+self._height),],
                        key=lambda x: x[2])
                        
    def _check_low(self, radius):
        """Is the radius greater than all corners of the pixel?

        Args:
            radius (float): The radius to check

        Return:
            bool: True if the radius is greater than all corners of the pixel,
                False if that is not true
                """
        r2 = radius**2
        if (self._x_0)**2 + (self._y_0)**2  < r2 and \
            (self._x_0+self._width)**2 + (self._y_0)**2  < r2 and \
            (self._x_0)**2 + (self._y_0+self._height)**2  < r2 and \
            (self._x_0+self._width)**2 + (self._y_0+self._height)**2  < r2:
                return True
                
        return False
    
    def _check_high(self, radius):
        """Is the radius less than all corners of the pixel?

        Args:
            radius (float): The radius to check

        Return:
            bool: True if the radius is less than all corners of the pixel,
                False if that is not true
                """
        r2 = radius**2
        if (self._x_0)**2 + (self._y_0)**2  > r2 and \
            (self._x_0+self._width)**2 + (self._y_0)**2  > r2 and \
            (self._x_0)**2 + (self._y_0+self._height)**2  > r2 and \
            (self._x_0+self._width)**2 + (self._y_0+self._height)**2  > r2:
                return True
                
        return False    
                        
    def _assign_to_case(self, radius):
        """ Given the radius, chooses what formula to use to calculate the area,
            and does so.

            I tried to optimize this function, but I'm not sure if I did a good job.

        Args:
            radius (float): Radius to check

        Return:
            float: 
        """
        
        
        if self._check_high(radius):
            return 0.
            
        elif self._check_low(radius):
            return 1.
            
        near, near_next, far_next, far = self._get_orientation()
            
        if near_next[2] > radius and far_next[2] > radius: # Nearby triangle
            x_diff = np.abs(np.abs(near[0]) - np.sqrt(radius**2-near[1]**2))
            y_diff = np.abs(np.abs(near[1]) - np.sqrt(radius**2-near[0]**2))
                        
            return self._triangle(x_diff, y_diff)
            
        elif near_next[2] <= radius and far_next[2] <= radius: # Far triangle
            x_diff = np.abs(np.abs(far[0]) - np.sqrt(radius**2-far[1]**2))
            y_diff = np.abs(np.abs(far[1]) - np.sqrt(radius**2-far[0]**2))
                        
            return self._removed_triangle(x_diff, y_diff)  
            
        elif near_next[2] <= radius and far_next[2] > radius: # Trapizoid
            if near_next[1] == near[1]: # X trapizoid
                h = self._width
                b1 = np.abs(np.abs(self._y_0) - np.sqrt(np.abs(radius**2 - self._x_0**2)))
                b2 = np.abs(np.abs(self._y_0) - np.sqrt(np.abs(radius**2 - (self._x_0+self._width)**2)))
                
                return self._trapizoid(b1, b2, h)

            else: # Y Trapizoid
                h = self._height
                b1 = np.abs(np.abs(self._x_0) - np.sqrt(np.abs(radius**2 - self._y_0**2)))
                b2 = np.abs(np.abs(self._x_0) - np.sqrt(np.abs(radius**2 - (self._y_0+self._height)**2)))
                
                return self._trapizoid(b1, b2, h)
                            
                            
    def _mc_area(self, radius, n = 1e3):
        """ Calculate the area of the portion of the pixel contained within
            the radius using a psuedo monte-carlo approach (really a grid, 
            since it's a lot simplier).

        Args:
            radius (float): The radius to check
            n (int, optional): The number of points to use along one axis

        Return:
            float: The area contained within the radius, to some precision
            """
        
        if self._check_low(radius):
            return 1.
            
        near, near_next, far_next, far = self._get_orientation()
            
        xs = np.linspace(self._x_0,self._width+self._x_0,n+1)**2
        ys = np.linspace(self._y_0,self._height+self._y_0,n+1)**2
        
        total = 0
        r2 = radius ** 2
        for x in xs:
            total += np.sum(ys+x < r2)
                
        return total / float(n**2)
        
                            
                            
    def fraction(self, radius):
        """Given the radius, use the appropriate method to find the fractional
            area of the pixel contained within that value.

        This method uses a grid/monte-carlo approach if the radius is small, and
            a linear approximation if the radius is large.

        Args:
            radius (float): The radius to check

        Return:
            float: Fractional area of the pixel corresponding to the radius.
            """
        if radius > np.max(self._width, self._height)*25.:
            return self._assign_to_case(radius)
        else:        
            return self._mc_area(radius)
            
            
    def value(self, radius=None):
        """Given the radius, return the fractional area of the pixel within that number,
        weighted by the value of the pixel.

        Args:
            radius (float, optional): The radius to check. If none, then the value of
                this pixel will be returned.

        Return:
            float: Weighted fractional area of the pixel corresponding to the radius.
            """
        if radius is None:
            return self._value
        else: 
            return self.fraction(radius) * self._value
        
        
def assert_fl(val, limit):
    """Convenience test function"""
    assert(np.abs(val) < limit), val
        
def test_cases():
    """Basic pixel tests"""
    # Basic cases
    f = Pixel(1.41421 - 0.5, -.5, 1, 1, 1)
    assert_fl(f._assign_to_case(.5), 1e-6)
    assert_fl(f._assign_to_case(3) - 1., 1e-6)
    
    # Trapizoids
    f = Pixel(1.41421 - 0.5, -.5, 1, 1, 1)
    assert_fl (f._assign_to_case(1.5) - .5, 1e-5)
    f = Pixel(-.5, 1.41421 - 0.5, 1, 1, 1)
    assert_fl (f._assign_to_case(1.5) - .5, 1e-5)
    f = Pixel(-1.41421 - 0.5, -.5, 1, 1, 1)
    assert_fl (f._assign_to_case(1.5) - .5, 1e-5)
    f = Pixel(-.5, -1.41421 - 0.5, 1, 1, 1)
    assert_fl (f._assign_to_case(1.5) - .5, 1e-5)
    
    # Near Triangles
    f = Pixel(1, 1, 1, 1, 1)
    assert_fl (f._assign_to_case(1.8027) - .125, 5e-5)
    f = Pixel(-2, 1, 1, 1, 1)
    assert_fl (f._assign_to_case(1.8027) - .125, 5e-5)
    f = Pixel(1, -2, 1, 1, 1)
    assert_fl (f._assign_to_case(1.8027) - .125, 5e-5)
    f = Pixel(-2, -2, 1, 1, 1)
    assert_fl (f._assign_to_case(1.8027) - .125, 5e-5)
    
    # Far triangles
    f = Pixel(1, 1, 1, 1, 1)
    assert_fl (f._assign_to_case(2.5) - .875, 5e-5)
    f = Pixel(-2, 1, 1, 1, 1)
    assert_fl (f._assign_to_case(2.5) - .875, 5e-5)
    f = Pixel(1, -2, 1, 1, 1)
    assert_fl (f._assign_to_case(2.5) - .875, 5e-5)
    f = Pixel(-2, -2, 1, 1, 1)
    assert_fl (f._assign_to_case(2.5) - .875, 5e-5)
    
def test_mc():
    """ Test the monte-carlo generator"""
    def build_f(angle):
        return Pixel(1.5*np.cos(angle/180. * np.pi) - .4, 1.5*np.sin(angle/180. * np.pi)-.4, 1, 1, 1)

    angles = np.empty(13)
    for i,a in enumerate(np.arange(0,360+30, 30)):
        f = build_f(a)
        angles[i] = f._mc_area(1.5, 3e3)
        
    test_angles = np.array([0.36821809999999999, 0.29881760000000002, 0.29883700000000002, 
                            0.36828100000000003, 0.41309459999999998, 0.48250900000000002, 
                            0.56814140000000002, 0.61303399999999997, 0.61308810000000002, 
                            0.5680347, 0.48233340000000002, 0.41288350000000001, 0.36825560000000002])
                                                    
    assert ( np.all(np.abs(angles - test_angles) < 5e-4) ), angles-test_angles
    
    
    
def determine_mc_number():
    """How many monte carlo points to use?"""
    f = Pixel(1.5*np.cos(10/180. * np.pi) - .4, 1.5*np.sin(10/180. * np.pi)-.4, 1, 1, 1)
    
    real = f._mc_area(1.5,1e5)
    
    ns = np.logspace(1,5,40)
    values = [None,]*40
    for i,n in enumerate(ns):
        sys.stdout.write('.')
        sys.stdout.flush()
        ntrials=1
        data = np.empty(ntrials)
        for j in np.arange(ntrials):
            data[j] = f._mc_area(1.5, n)
        values[i] = (n, np.mean(data), np.std(data))
    print
    
    x,y,st = zip(*values)
    print st
        
    import matplotlib.pyplot as plt
    
    plt.figure()
    ax = plt.subplot(211)
    ax.set_xscale('log')
    plt.errorbar(x,y,yerr=st)
    
    ax = plt.subplot(212)
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.plot(x,y-real)
    
    plt.show()
    
def determine_switch_radius():
    """When do we switch from the monte carlo method to the approximation?"""
    size = 1
    max_val = size*100
    angle = 0
    
    import matplotlib.pyplot as plt
    
    n = 201
    radius = np.linspace(0,max_val,n)
    mc_data = np.empty(n)
    approx_data = np.empty(n)
    for i,r in enumerate(radius):      
        sys.stdout.write('.')
        sys.stdout.flush()  
        f = Pixel(r*np.cos(0/180. * np.pi) - size/2., r*np.sin(angle/180. * np.pi)-size/2., size, size, 1)
        
        mc_data[i] = f._mc_area(r)
        approx_data[i] = f._assign_to_case(r)
    print
    
    plt.figure()
    ax = plt.subplot(211)
    
    plt.plot(radius/size, mc_data, 'b', radius/size, approx_data, 'r')
    
    
    ax = plt.subplot(212)
    ax.set_yscale('log')
    plt.plot(radius/size, mc_data - approx_data)
    
    plt.show()
    
    
        
    
#@do_profile(follow=[Pixel._mc_area])                
def test():
    """Debugging test"""
    p = Pixel(1, 1, 1, 1, 1)
    for i in xrange(100):
        print p._mc_area(1.8027, 1e3)
        
    
import cProfile
    
if __name__ == "__main__":
    #test()
    
    #cProfile.run("determine_switch_radius()", sort='cumtime')
    #determine_mc_number()
    test_cases()
    test_mc()
        
