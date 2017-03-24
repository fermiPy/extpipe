import os
import sys

configdir='haloanalysis/config'

model_lookup = {
    'stdmodel' : None,
    'altmodel' : 'altmodel',
    'model0' : 'Lorimer_z4_Ts100000_v5',
    'model1' : 'Lorimer_z4_Ts150_v5',
    'model2' : 'Lorimer_z10_Ts100000_v5',
    'model3' : 'Lorimer_z10_Ts150_v5',
    'model4' : 'SNR_z4_Ts100000_v5',
    'model5' : 'SNR_z4_Ts150_v5',
    'model6' : 'SNR_z10_Ts100000_v5',
    'model7' : 'SNR_z10_Ts150_v5',
}
    

#models = model_lookup.keys()
models = ['stdmodel']

script=sys.argv[1]
#sourcelist=sys.argv[2]
basedir=sys.argv[2]

#sourcelist = 'haloanalysis/diffuse_syst_list.yaml'
#sourcelist = 'galactic_list.yaml'
#sourcelist = 'haloanalysis/ext_candidates_list.yaml'
#sourcelist = 'haloanalysis/fhes_list.yaml'

sourcelist = ['haloanalysis/fhes_list.yaml',
              'haloanalysis/ext_candidates_list.yaml',
              'haloanalysis/cta1_list.yaml']

#configs = ['config_std.yaml','config_90mo.yaml','config_bigroi.yaml']
configs = ['config_std.yaml','config_psf0123_joint2a.yaml']
#configs = ['config_std.yaml']


for m in models:

    cmd = 'fermipy-clone-configs --script=%s'%(script)
    cmd += ' --basedir=%s_%s '%(basedir,m)

    for s in sourcelist:
        cmd += ' --source_list=%s '%s
    
    #if script == 'run-region-analysis' and ('fhes' in sourcelist or 'cta1' in sourcelist):
    #    cmd += ' --args="--radius=1.5" '    
    #cmd += '  --basedir=bigroi_all_galactic_scan_%s'%m
    for c in configs:
        cmd += ' %s '%(os.path.join(configdir,c))
    if model_lookup[m] is not None:
        cmd += ' %s '%(os.path.join(configdir,'config_%s.yaml'%(model_lookup[m])))
    print cmd
    os.system(cmd)
    
