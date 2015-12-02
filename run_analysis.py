import sys
from fermipy.gtanalysis import GTAnalysis

import argparse

usage = "usage: %(prog)s [config file]"
description = "Run fermipy analysis chain."
parser = argparse.ArgumentParser(usage=usage,description=description)

parser.add_argument('--config', default = 'sample_config.yaml')
parser.add_argument('--source', default = None)

args = parser.parse_args()

gta = GTAnalysis(args.config,logging={'verbosity' : 3})

gta.setup()
gta.optimize()
gta.extension(args.source)
gta.sed(args.source)

# Baseline Fit
gta.write_roi('fit0',make_residuals=True)

newname = args.source.replace(' ','').lower() + '_reloc'
gta.localize(args.source,update=True,newname=newname,nstep=6)

gta.optimize()
gta.extension(newname)
gta.sed(newname)

# Pot-Relocalization Fit
gta.write_roi('fit1',make_residuals=True)

skydir = gta.roi.get_source_by_name(newname)[0].skydir

halo_source_dict = {
    'ra' : skydir.ra.deg,
    'dec' : skydir.dec.deg,
    'SpectrumType' : 'PowerLaw', 
    'Index' : 2.0, 
    'Scale' : 1000,
    'Prefactor' : { 'value' : 0.0, 'scale' : 1e-13 },
    'SpatialModel' : 'GaussianSource',
    'SpatialWidth' : 1.0
    }


halo_width = [0.1,0.316,1.0]


for i,w in enumerate(halo_width):
    halo_source_dict['SpatialWidth'] = w

    halo_source_name = 'halo_%06.3f'%w
    
    gta.add_source(halo_source_name,halo_source_dict)
    gta.free_norm(halo_source_name)
    gta.free_norm(newname)
    gta.free_norm('galdiff')
    gta.free_norm('isodiff')
    gta.fit()
    gta.sed(halo_source_name)    
    gta.write_roi('halo_fit_%02i'%i,make_residuals=True)
    gta.delete_source(halo_source_name)    
    gta.load_roi('fit1')
