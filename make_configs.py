import yaml
from fermipy.utils import *

from fermipy.roi_model import ROIModel

import sys
import time, os, stat

import argparse

usage = "usage: %(prog)s [config files]"
description = "Run fermipy analysis chain."
parser = argparse.ArgumentParser(usage=usage,description=description)

parser.add_argument('--basedir', default = None, required=True)
parser.add_argument('--source_list', default = None, required=True,
                    help='YAML file containing a list of sources to be '
                    'analyzed.')
parser.add_argument('--script', default = None, required=True,
                    help='The python script.')
parser.add_argument('configs', nargs='+', default = None,
                    help='One or more configuration files that will be merged '
                    'to construct the analysis configuration.')

args = parser.parse_args()

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

config = {}
for c in args.configs:
    config.update(yaml.load(open(c)))

src_list = yaml.load(open(args.source_list))
basedir = args.basedir

roi = ROIModel(config['model'])
roi.load()

bash_script = """
cat $0
python {script} --config={config} --source="{source}"
"""

for k in src_list:

    s = roi.get_source_by_name(k,True)
    dirname = os.path.join(basedir,s.name.lower().replace(' ','_'))
    
    mkdir(dirname)

    config['selection']['target'] = s.name    
    cfgfile = os.path.abspath(os.path.join(dirname,'config.yaml'))
    yaml.dump(config,open(cfgfile,'w'),default_flow_style=False)

    script = os.path.basename(args.script)
    scriptpath = os.path.abspath(os.path.join(dirname,script))
    
    os.system('cp %s %s'%(args.script,scriptpath))
    
    with open(os.path.join(dirname,'run.sh'),'wt') as f:
        f.write(bash_script.format(source=s.name,config=cfgfile,
                                   script=scriptpath))
    
