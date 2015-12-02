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


config = yaml.load(open(sys.argv[1]))
config.update(yaml.load(open(sys.argv[2])))

src_list = yaml.load(open(sys.argv[3]))
basedir = sys.argv[4]

roi = ROIModel(config['model'])
roi.load()

for k in src_list:

    s = roi.get_source_by_name(k,True)
    dirname = os.path.join(basedir,s.name.lower().replace(' ','_'))
    
    mkdir(dirname)

    config['selection']['target'] = s.name    
    cfgfile = os.path.join(dirname,'config.yaml')
    yaml.dump(config,open(cfgfile,'w'),default_flow_style=False)
    

