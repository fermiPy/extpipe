import numpy as np
from parameter_set import PhysicalParameters as PhysPar
from fits_parser import FitsParser
import os.path
import pickle
import numpy.lib.recfunctions


def build_row(par, energy, profile):
    """A convenience function which takes data extracted from the image
    directory (physical parameter set and image file), and puts it into
    a numpy structured array.
    
    Args:
        par (PhysicalParameter): For this row. Used to populate the first 9 columns
        energy (float): The energy in MeV
        profile (numpy array, 9 digits long): The array, assumed to be containment values
            at (10%, 20%, 30%, 40%, 50%, 60%, 70%, 80%, 90%, 95%)

    Returns:
        numpy array: A structured array with the dtype:
            dtype=[('redshift','d'), ('strength','d'), ('power_law_index','d'), ('boost','d'),
                    ('correlation_length','d'), ('angle','d'), ('low_energy','d'), ('high_energy','d'),
                    ('total_photons', 'd'), ('energy', 'd'),('10% Containment','d'),
                    ('20% Containment','d'),('30% Containment','d'),('40% Containment','d'),
                    ('50% Containment','d'),('60% Containment','d'),('70% Containment','d'),
                    ('80% Containment','d'),('90% Containment','d'),('95% Containment','d')]
            Columns are as indicated
    """
    # Order is redshift, strength, power_law_index, boost,
    # correlation_length, angle, low_energy, high_energy, total_photons,
    # energy, rad at 0.1 containment, @0.2, @0.3, @0.4, @0.5, @0.6,
    # @0.7, @0.8, @0.9, @0.95
    row = np.array([(par.redshift, par.strength, par.power_law_index, par.boost, par.correlation_length,
            par.angle, par.low_energy, par.high_energy, par.total_photons, energy, profile[0],
            profile[1], profile[2], profile[3], profile[4], profile[5], profile[6], profile[7],
            profile [8], profile[9])],
            dtype=[('redshift','d'), ('strength','d'), ('power_law_index','d'), ('boost','d'),
                    ('correlation_length','d'), ('angle','d'), ('low_energy','d'), ('high_energy','d'),
                    ('total_photons', 'd'), ('energy', 'd'),('10% Containment','d'),
                    ('20% Containment','d'),('30% Containment','d'),('40% Containment','d'),
                    ('50% Containment','d'),('60% Containment','d'),('70% Containment','d'),
                    ('80% Containment','d'),('90% Containment','d'),('95% Containment','d'),
                    ])
    return row

def rows_to_array(rows):
    """ Turns a set of rows produced by build_rows into a full numpy array.

    Args:
        rows (list): A list of numpy rows, created by build_row

    Return:
        numpy array: The above concatenated together
        """
    return np.concatenate(rows)


def handle_directory(directory):
    """This function extracts all relevant information from a given directory, assuming
    standard filenames

    Args:
        directory (string): Absolute or relative path to the directory in question

    Return:
        list: A list of build_row rows, one for each energy in the contained "halo" file
        """
    print '\n', directory
    profiles = {}
    phys_param = None
    for filename in os.listdir(directory):
        if 'halo' in filename:
            fp = FitsParser(os.path.join(directory, filename))
            for energy, hist in fp.hist_set():
                print '\t\t',energy
                profiles[energy] = hist.radius(np.concatenate((np.arange(.1, 1., .1),[.95]),axis=1))

        elif filename == 'phys_param.pkl':
            with open(os.path.join(directory, filename), 'rb') as param_file:
                phys_param = pickle.load(param_file)
                


    rows = []
    for key, value in profiles.iteritems():
        rows.append(build_row(phys_param, key, value))

    return rows


def run_all_directories(base = 'images'):
    """Runs handle_directory for each directory in the base directory

    At the moment, this function caches the output after each handle_directory call,
    in case of crashes. To take advantage of the cache, the function would
    have to be extended--that functionality is not present yet.

    Args:
        base (string): Directory containing the image directories

    Return:
        list: A combined list of the rows returned from handle_directory 
            applied to each directory.
            """
    full_rows = []
    for i,dirname in enumerate(os.listdir(base)):
        rows = handle_directory(os.path.join(base,dirname))
        full_rows += rows
        with open('cache.pkl', 'wb') as cachefile:
            pickle.dump(full_rows, cachefile)

    os.remove('cache.pkl')
    
    return full_rows

def make_array(rows):
    return np.concatenate(rows)


if __name__ == '__main__':
    rows = run_all_directories('images')
    array = rows_to_array(rows)
    np.save('image_reduction',array)

    print array
    


