import numpy as np
import glob
import sys
import os
from astropy.table import Table, Column
from fermipy.utils import mkdir
from jinja2 import Template, Environment, FileSystemLoader
from fermipy.utils import collect_dirs
#from haloanalysis.batch import *
import argparse

usage = "usage: %(prog)s"
description = "Aggregate analysis output."
parser = argparse.ArgumentParser(usage=usage,description=description)


parser.add_argument('--output', default = 'test.fits')
parser.add_argument('--table', default = None, required=True)
parser.add_argument('--path', default = None, required=True)
parser.add_argument('dirs', nargs='+', default = None,
                    help='Run analyses in all subdirectories of this '
                    'directory.')

args = parser.parse_args()
#dirs = sorted(collect_dirs(args.dirs))
#dirs = sorted([d for s in _ for _ in collect_dirs(args.dirs)])
dirs = sorted([d for argdir in args.dirs for d in collect_dirs(argdir)])



#dirs = glob.glob(sys.argv[1] + '/3fgl_j0*')
tab = Table.read(args.table)


#rootpath = os.path.expandvars('$HOME/public_html/extanalysis')
rootpath = os.path.abspath(args.path)
cwd = os.getcwd()

#env = Environment(loader=FileSystemLoader(os.path.join(PATH, 'templates')))
env = Environment(loader=FileSystemLoader('/u/gl/mdwood/ki20/mdwood/fermi/ext_analysis'))

files = ['fit[1-5]_pointsource_powerlaw_2.00_tsmap_sqrt_ts.png',
#         'fit*_pointsource_powerlaw_2.00_tsmap_npred.png',
         'fit[1-5]_nosource_pointsource_powerlaw_2.00_tsmap_sqrt_ts.png',
#         'fit*_nosource_pointsource_powerlaw_2.00_tsmap_npred.png',
         'fit[1-5]_*_residmap_data.png',
         'fit[1-5]_*_residmap_model.png']
#         'fit*_pointsource_powerlaw_2.00_residmap_sigma.png']
         
for d in sorted(dirs):

    codename = os.path.basename(d)
    rows = tab[tab['codename']==codename]

    if len(rows) == 0:
        continue
    
    newdir = os.path.join(rootpath,os.path.basename(d))

    print newdir
    
    mkdir(newdir)

    for f in files:

        input_files = sorted(glob.glob(os.path.join(os.path.abspath(d),f)))
        if len(input_files) > 1:
            input_files = [input_files[0],input_files[-1]]
        
#        if not os.path.isfile(os.path.join(d,f)):
#            continue

        for t in input_files:
                                
            dest = os.path.join(newdir,os.path.basename(t))
        
            if os.path.exists(dest):
                os.system('rm %s'%dest)

            print 'ln -s %s %s'%(t,dest)            
            os.system('ln -s %s %s'%(t,dest))
    
    row = rows[0]
    idx = row['fit_nsrc']

    src_dict = dict(source_name=d,
                    name=row['name'],
                    assoc=row['assoc'],
                    src_class=row['class'],
                    spectrum_type=row['spectrum_type'],
                    fit1_tsmap_sqrt_ts='fit1_pointsource_powerlaw_2.00_tsmap_sqrt_ts.png',
                    fit1_nosource_tsmap_sqrt_ts='fit1_nosource_pointsource_powerlaw_2.00_tsmap_sqrt_ts.png',
                    fit1_tsmap_npred='fit1_pointsource_powerlaw_2.00_tsmap_npred.png',
                    fit1_residmap_data='fit1_gaussian_s0.20_powerlaw_2.00_residmap_data.png',
                    fit1_residmap_model='fit1_gaussian_s0.20_powerlaw_2.00_residmap_model.png',
                    fit1_residmap_sigma='fit1_pointsource_powerlaw_2.00_residmap_sigma.png',
                    fit2_tsmap_sqrt_ts='fit%i_pointsource_powerlaw_2.00_tsmap_sqrt_ts.png'%idx,
                    fit2_nosource_tsmap_sqrt_ts='fit%i_nosource_pointsource_powerlaw_2.00_tsmap_sqrt_ts.png'%idx,                    
                    fit2_nosource_tsmap_npred='fit%i_nosource_pointsource_powerlaw_2.00_tsmap_npred.png'%idx,                    
                    fit2_residmap_data='fit%i_gaussian_s0.20_powerlaw_2.00_residmap_data.png'%idx,
                    fit2_residmap_model='fit%i_gaussian_s0.20_powerlaw_2.00_residmap_model.png'%idx,
                    fit2_residmap_sigma='fit%i_pointsource_powerlaw_2.00_residmap_sigma.png'%idx,
                    )
    
#              'ext1_ts','ext1_mle','ext1_ul95','ext1_err',

    fcols2 = ['ts','npred','dfde1000_index']
#              'fit0_ts','fit1_ts','fit2_ts','fit3_ts','fit4_ts',
#              'fit1_dlike','fit2_dlike','fit3_dlike','fit4_dlike']
             
#    fcols3 = ['fit0_offset','fit1_offset','fit2_offset','fit3_offset','fit4_offset']

    fcols3 = ['ra','dec','glon','glat']
    
    gcols = ['flux10000','flux10000_err',
             'eflux10000','eflux10000_err',
             'flux1000','flux1000_err',
             'eflux1000','eflux1000_err',
             'flux100','flux100_err',
             'eflux100','eflux100_err']

#    for c in ['fitn_ext_ts']:
#        for i in range(len(row[c])):
#            if not np.isfinite(row[c][i]):

        
    for i in range(6):
        src_dict['fit%i_offset'%i] = '%10.3f'%row['fitn_offset'][i]
        src_dict['fit_ext_ts'] = '%10.2f'%row['fit_ext_ts']
        src_dict['fit_ext_err'] = '%10.3f'%row['fit_ext_err']
        src_dict['fit_ext_mle'] = '%10.3f'%row['fit_ext_mle']
        src_dict['fit_ext_ul95'] = '%10.3f'%row['fit_ext_ul95']
        src_dict['fit_halo_ts'] = '%10.3f'%row['fit_halo_ts']
        src_dict['fit_halo_width'] = '%10.3f'%row['fit_halo_width']
        src_dict['fit_halo_index'] = '%10.3f'%row['fit_halo_index']
        
        src_dict['ext%i_ts'%i] = '%10.2f'%row['fitn_ext_ts'][i]
        src_dict['ext%i_err'%i] = '%10.3f'%row['fitn_ext_err'][i]
        src_dict['ext%i_mle'%i] = '%10.3f'%row['fitn_ext_mle'][i]
        src_dict['ext%i_ul95'%i] = '%10.3f'%row['fitn_ext_ul95'][i]
        src_dict['fit%i_halo_ts'%i] = '%10.2f'%row['fitn_halo_ts'][i]
        src_dict['fit%i_halo_width'%i] = '%10.3f'%row['fitn_halo_width'][i]
        src_dict['fit%i_halo_index'%i] = '%10.3f'%row['fitn_halo_index'][i]
        src_dict['fit%i_dlike'%i] = '%10.2f'%row['fitn_dlike'][i]
        src_dict['fit%i_dlike_ext'%i] = '%10.2f'%row['fitn_dlike_ps_ext'][i]
        src_dict['fit%i_dlike_halo'%i] = '%10.2f'%row['fitn_dlike_ps_halo'][i]
        src_dict['fit%i_ts'%i] = '%10.2f'%row['fitn_ts'][i]
        src_dict['fit%i_ra'%i] = '%10.3f'%row['fitn_ra'][i]
        src_dict['fit%i_dec'%i] = '%10.3f'%row['fitn_dec'][i]
    
    for c in fcols2:        
        src_dict[c] = '%10.2f'%row[c]

    for c in fcols3:
        src_dict[c] = '%10.3f'%row[c]
        
    for c in gcols:
        src_dict[c] = '%10.3g'%row[c]

    for k, v in src_dict.items():
        if v.strip() == 'nan':
            src_dict[k] = '---'
            
#    for i, index in enumerate(['1.50','1.75','2.00','2.25','2.50','2.75','3.00']):        
#        for j, ext in enumerate(['0.100','0.178','0.316','0.562','0.750','1.000']):
#            src_dict['fit1_halo_idx%02i_ext%02i_ts'%(i,j)] = '%10.2f'%row['fit1_halo_%s_%s_ts'%(index,ext)]
#            src_dict['fit_halo_idx%02i_ext%02i_ts'%(i,j)] = '%10.2f'%row['fit_halo_%s_%s_ts'%(index,ext)]
    
    template = env.get_template('source_template.html')
    template = template.render(**src_dict)

    src_page = os.path.join(newdir,'index.html')
    
    with open(src_page, 'w') as f:
        f.write(template)
