import os
import sys
import argparse

configdir='haloanalysis/config'

model_lookup = {
    'stdmodel' : None,
    'altmodel' : 'altmodel',
    'p8iem' : 'p8iem',
    'model0' : 'Lorimer_z4_Ts100000_v5',
    'model1' : 'Lorimer_z4_Ts150_v5',
    'model2' : 'Lorimer_z10_Ts100000_v5',
    'model3' : 'Lorimer_z10_Ts150_v5',
    'model4' : 'SNR_z4_Ts100000_v5',
    'model5' : 'SNR_z4_Ts150_v5',
    'model6' : 'SNR_z10_Ts100000_v5',
    'model7' : 'SNR_z10_Ts150_v5',
}
    

def main():

    usage = "usage: %(prog)s"
    description = "Aggregate analysis output."
    parser = argparse.ArgumentParser(usage=usage,description=description)


    parser.add_argument('--config', default = None, required=True)
    parser.add_argument('--script', default = None, required=True)
    parser.add_argument('--args', default = None)
    parser.add_argument('--model', default = 'stdmodel')
    parser.add_argument('sourcelists', nargs='+', default = None,
                        help='')    
    args = parser.parse_args()

    if args.model == 'all':    
        models = model_lookup.keys()
    else:
        models = [args.model]

    #sourcelist = 'haloanalysis/diffuse_syst_list.yaml'
    #sourcelist = 'galactic_list.yaml'
    #sourcelist = 'haloanalysis/ext_candidates_list.yaml'
    #sourcelist = 'haloanalysis/fhes_list.yaml'

    #sourcelists = [
        #'haloanalysis/sourcelists/fhes_infill_list.yaml',]
    #    'haloanalysis/sourcelists/fhes_list.yaml',
    #    'haloanalysis/sourcelists/ext_candidates_list.yaml',
    #    'haloanalysis/sourcelists/cta1_list.yaml',]
        #'haloanalysis/sourcelists/3fgl_srcs_list_glat050.yaml',
        #'haloanalysis/sourcelists/3fhl_srcs_list_glat050.yaml']

    if args.config == 'std_all':
        configs = ['config_std.yaml']
    elif args.config == 'std_psf0123_joint2a':
        configs = ['config_std.yaml','config_psf0123_joint2a.yaml']
    else:
        raise ValueError
        
    for m in models:

        cmd = 'fermipy-clone-configs --script=%s'%(args.script)
        cmd += ' --basedir=%s_%s '%(args.config,m)

        for s in args.sourcelists:
            cmd += ' --source_list=%s '%s

        if args.args:
            cmd += ' --args=%s '%args.args
            
        #cmd += '  --basedir=bigroi_all_galactic_scan_%s'%m
        for c in configs:
            cmd += ' %s '%(os.path.join(configdir,c))
        if model_lookup.get(m) is not None:
            cmd += ' %s '%(os.path.join(configdir,'config_%s.yaml'%(model_lookup[m])))
        print cmd
        os.system(cmd)
    
if __name__ == '__main__':

    main()
