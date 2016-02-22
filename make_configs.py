import os
import copy
import yaml
import numpy as np
import healpy as hp
import fermipy.utils as utils
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

config = {}
for c in args.configs:
    config.update(yaml.load(open(c)))

src_list = yaml.load(open(args.source_list))
basedir = args.basedir

bash_script = """
cat $0
python {script} --config={config} --source="{source}"
"""

scriptdir = os.path.join(basedir,'scripts')
utils.mkdir(scriptdir)
os.system('cp %s %s'%(args.script,scriptdir))

for name, v in src_list.items():

#    if isinstance(target,dict):
#        name = 'hp_region_%03i_%04i'%(target['nside'],target['pix'])
#        theta, phi = hp.pix2ang(target['nside'],target['pix'])        
#        config['selection']['glat'] = np.degrees(np.pi/2.-theta)
#        config['selection']['glon'] = np.degrees(phi)
#    else:
#        name = target
#        config['selection']['target'] = name 

    print name
    
#    dirname = os.path.join(basedir,name.lower().replace(' ','_'))
    dirname = os.path.join(basedir,name)    
    utils.mkdir(dirname)

    c = copy.deepcopy(config)
    c = utils.merge_dict(c,v,add_new_keys=True)
    
    cfgfile = os.path.abspath(os.path.join(dirname,'config.yaml'))
    yaml.dump(utils.tolist(c),open(cfgfile,'w'),default_flow_style=False)

    script = os.path.basename(args.script)
    scriptpath = os.path.abspath(os.path.join(dirname,script))
    
    os.system('ln -sf %s %s'%(os.path.abspath(os.path.join(scriptdir,script)),
                             scriptpath))

    runscript = os.path.abspath(os.path.join(dirname,os.path.splitext(script)[0] + '.sh'))
    
    with open(os.path.join(runscript),'wt') as f:
        f.write(bash_script.format(source=name,config=cfgfile,
                                   script=scriptpath))
    
