import glob

import yaml
from fermipy.utils import *

from fermipy.roi_model import ROIModel

import sys
import time, os, stat

def file_age_in_seconds(pathname):
    return time.time() - os.stat(pathname)[stat.ST_MTIME]

def check_log(logfile, string='Successfully', exists=True):
    """ Often logfile doesn't exist because the job hasn't begun
    to run. It is unclear what you want to do in that case...
    logfile : String with path to logfile
    exists  : Is the logfile required to exist
    string  : Value to check for in existing logfile
    """
    if not os.path.exists(logfile):
	return not exists
    return string in open(logfile).read()


overwrite = True


dirs = glob.glob(sys.argv[1] + '/*')

script = os.path.abspath(sys.argv[2])

for dirname in dirs:

    cfgfile = os.path.join(dirname,'config.yaml')
    logfile = os.path.join(dirname,'lsf.log')

    config = yaml.load(open(cfgfile))

    srcname = config['selection']['target'].replace(' ','').lower()
    
    cmd = 'bsub -W 800 -o %s python %s --config=%s --source=%s'%(logfile,script,cfgfile,
                                                                 srcname)

    print cmd
    
    age = 1000
    if os.path.isfile(logfile):
        age = file_age_in_seconds(logfile)/60.
    

#    print age/60., check_log(logfile)
    
    if overwrite or (not check_log(logfile) and age > 60):
        print cmd
        os.system(cmd)
    else:
        print 'Skipping ', dirname, age



    

    
