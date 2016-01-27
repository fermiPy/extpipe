import os
import yaml
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

for name in src_list:

    dirname = os.path.join(basedir,name.lower().replace(' ','_'))    
    utils.mkdir(dirname)

    config['selection']['target'] = name    
    cfgfile = os.path.abspath(os.path.join(dirname,'config.yaml'))
    yaml.dump(utils.tolist(config),open(cfgfile,'w'),default_flow_style=False)

    script = os.path.basename(args.script)
    scriptpath = os.path.abspath(os.path.join(dirname,script))
    
    os.system('cp %s %s'%(args.script,scriptpath))
    
    with open(os.path.join(dirname,'run.sh'),'wt') as f:
        f.write(bash_script.format(source=name,config=cfgfile,
                                   script=scriptpath))
    
