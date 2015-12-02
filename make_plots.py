import sys
from fermipy.gtanalysis import *
from fermipy.config import *
from fermipy.utils import *
import yaml
import pprint
from fermipy.logger import *

        
config = ConfigManager.create(sys.argv[1])
gta = GTAnalysis(config,logging={'verbosity' : 3})

gta.setup()

gta.load_roi('fit1')
gta.write_roi('fit1',make_residuals=True)

