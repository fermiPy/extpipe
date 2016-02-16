import os
import subprocess

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

def get_lsf_status():
    
    p = subprocess.Popen(['bjobs'],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    p.stderr.close()

    output = p.stdout.readlines()

    status_count = {'RUN' : 0,
                    'PEND' : 0,
                    'SUSP' : 0,
                    'USUSP': 0,
                    'NJOB' : 0,
                    'UNKNWN' : 0}
    
    for line in output[1:]:
        line = line.strip().split()

        status_count['NJOB'] += 1
        
        for k in status_count.keys():

            if line[2] == k:
                status_count[k] += 1
                
    return status_count
