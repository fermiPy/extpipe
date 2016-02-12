import glob

import yaml
from fermipy.utils import *

from fermipy.roi_model import ROIModel
from batch import *

import sys
import time, os, stat
import datetime
import argparse

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

usage = "usage: %(prog)s [config file]"
description = "Dispatch analysis jobs."
parser = argparse.ArgumentParser(usage=usage,description=description)

parser.add_argument('--config', default = 'sample_config.yaml')
parser.add_argument('--max_jobs', default = 500, type=int)
parser.add_argument('--jobs_per_cycle', default = 20, type=int)
parser.add_argument('--time_per_cycle', default = 15, type=float,
                    help='Time per submission cycle in seconds.')
parser.add_argument('--max_job_age', default = 90, type=float,
                    help='Max job age in minutes.')
parser.add_argument('--dry_run', default = False, action='store_true')
parser.add_argument('--overwrite', default = False, action='store_true')
parser.add_argument('--runscript', default = None, required=True)

parser.add_argument('dirs', nargs='+', default = None,
                    help='Run analyses in all subdirectories of this '
                    'directory.')

args = parser.parse_args()
    
dirs = []
for d in args.dirs:

    if not os.path.isdir(d):
        continue

    #subdirs = glob.glob(d + '/*')
    
    for subdir in os.listdir(d):
        
        subdir = os.path.join(d,subdir)
        
        if not os.path.isdir(subdir):
            continue
         
        dirs += [subdir]
         
    dirs += [d]


    
def collect_jobs(dirs,runscript,overwrite=False): 

    jobs = []
    
    for dirname in dirs:        
        
        o = dict(cfgfile = os.path.join(dirname,'config.yaml'),
                 logfile = os.path.join(dirname,os.path.splitext(runscript)[0] + '.log'),
                 runscript = os.path.join(dirname,runscript))
        
        if not os.path.isfile(o['cfgfile']):
            continue

        if not os.path.isfile(o['runscript']):
            continue

        if not os.path.isfile(o['logfile']):
            jobs.append(o)
            continue
            
        age = file_age_in_seconds(o['logfile'])/60.

        print dirname, check_log(o['logfile']), age
        
        if overwrite or (check_log(o['logfile'])=='Exited'): 
            print "Job Exited, resending command:"
            jobs.append(o)
        elif (check_log(o['logfile'])=='None') and age > 90:
            print "Job did not exit, but no activity on log file for > 90 min. Resending command:"
            jobs.append(o)
        elif not check_log(o['logfile']):
            jobs.append(o)

    return jobs

jobs = collect_jobs(dirs,args.runscript,args.overwrite)

while(1):

    print '-'*80
    print datetime.datetime.now()
    print len(jobs), 'jobs in queue'

    print args.dry_run
    
    if len(jobs) == 0:
        break
    
    status = get_lsf_status()
    
    njob_to_submit = min(args.max_jobs - status['NJOB'],
                         args.jobs_per_cycle)

    import pprint
    pprint.pprint(status)    
    print 'njob_to_submit ', njob_to_submit
    
    if njob_to_submit > 0:
        
        print 'Submitting ', njob_to_submit, 'jobs'
    
        for job in jobs[:njob_to_submit]:
            cmd = 'bsub -W 1000 -oo %s bash %s'%(job['logfile'],
                                                job['runscript'])
            print cmd
            if not args.dry_run:
                print 'submitting'
                os.system(cmd)

        del jobs[:njob_to_submit]

    print 'Sleeping %f seconds'%args.time_per_cycle
    sys.stdout.flush()
    time.sleep(args.time_per_cycle)


    
