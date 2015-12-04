import glob

import yaml
from fermipy.utils import *

from fermipy.roi_model import ROIModel

import sys
import time, os, stat

def file_age_in_seconds(pathname):
    return time.time() - os.stat(pathname)[stat.ST_MTIME]

def check_log(logfile, exited='Exited with exit code',
              successful='Successfully completed', exists=True):
    """ Often logfile doesn't exist because the job hasn't begun
    to run. It is unclear what you want to do in that case...
    logfile : String with path to logfile
    exists  : Is the logfile required to exist
    string  : Value to check for in existing logfile
    """
    if not os.path.exists(logfile):
	return not exists

    if exited in open(logfile).read():
        return 'Exited'
    elif successful in open(logfile).read():
        return 'Successful'
    else:
        return 'None' 
    #return string in open(logfile).read()


overwrite = False


dirs = glob.glob(sys.argv[1] + '/*')

for dirname in dirs:

    cfgfile = os.path.join(dirname,'config.yaml')
    logfile = os.path.join(dirname,'lsf.log')
    runscript = os.path.join(dirname,'run.sh')
    
    if not os.path.isfile(cfgfile): continue

    config = yaml.load(open(cfgfile))

    srcname = config['selection']['target'].replace(' ','').lower()
    
    cmd = 'bsub -W 800 -oo %s bash %s'%(logfile,runscript)

    age = 1000
    if os.path.isfile(logfile):
        age = file_age_in_seconds(logfile)/60.

    print check_log(logfile), age
    
    if overwrite or (check_log(logfile)=='Exited'): # and age > 60
        print "Job Exited, resending command:"
        print cmd
        os.system(cmd)
    elif (check_log(logfile)=='None') and age > 60:
        print "Job did not exit, but no activity on log file for > 60 min. Resending command:"
        print cmd
        os.system(cmd)
    elif not check_log(logfile):
        print cmd
        os.system(cmd)
    else:
        print 'Skipping ', dirname, age



    

    
