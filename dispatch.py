import glob

import yaml
from fermipy.utils import *

from fermipy.roi_model import ROIModel
from batch import *

import sys
import time, os, stat

def check_num_jobs():
    pass

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
max_jobs = 500
jobs_per_cycle = 25
time_per_cycle = 10

dirs = glob.glob(sys.argv[1] + '/*')

def collect_jobs(dirs,overwrite=False): 

    jobs = []
    
    for dirname in dirs:
        
        cfgfile = os.path.join(dirname,'config.yaml')
        logfile = os.path.join(dirname,'lsf.log')
        runscript = os.path.join(dirname,'run.sh')
        o = dict(cfgfile = os.path.join(dirname,'config.yaml'),
                 logfile = os.path.join(dirname,'lsf.log'),
                 runscript = os.path.join(dirname,'run.sh'))
        
        if not os.path.isfile(cfgfile):
            continue

        if not os.path.isfile(logfile):
            jobs.append(o)
            continue
            
        age = file_age_in_seconds(logfile)/60.

        print dirname, check_log(logfile), age
        
        if overwrite or (check_log(logfile)=='Exited'): 
            print "Job Exited, resending command:"
            jobs.append(o)
        elif (check_log(logfile)=='None') and age > 60:
            print "Job did not exit, but no activity on log file for > 60 min. Resending command:"
            jobs.append(o)
        elif not check_log(logfile):
            jobs.append(o)

    return jobs

jobs = collect_jobs(dirs,overwrite)

while(1):

    if len(jobs) == 0:
        break
    
    status = get_lsf_status()

    njob_to_submit = min(max_jobs - status['NJOB'],jobs_per_cycle)
    if njob_to_submit > 0:
        
        print 'Submitting ', njob_to_submit, 'jobs'
    
        for job in jobs[:njob_to_submit]:
            cmd = 'bsub -W 800 -oo %s bash %s'%(job['logfile'],
                                                job['runscript'])
            print cmd
            os.system(cmd)

        del jobs[:njob_to_submit]
    
    time.sleep(time_per_cycle)

    
