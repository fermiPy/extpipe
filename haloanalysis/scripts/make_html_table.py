import sys
import numpy as np
from lxml import html
from astropy.table import Table, Column
import argparse
import os

#  <script src="http://code.jquery.com/jquery-1.12.0.min.js"></script>
#  <script src="http://cdn.datatables.net/1.10.10/js/jquery.dataTables.min.js"></script>
#  <script type="text/javascript" charset="utf-8">

# <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/r/bs-3.3.5/jq-2.1.4,dt-1.10.8/datatables.min.css"/>

def main():

    usage = "usage: %(prog)s [config files]"
    description = "Merge tables."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--output', default = None)   
    parser.add_argument('file', default=None, nargs=1,
                        help='Catalog file.')

    args = parser.parse_args()
    
    tab = Table.read(args.file[0])

    if args.output is None:
        args.output = os.path.splitext(args.file[0])[0] + '.html'
    
    html_table_cols = ['linkname','assoc','class','ra','dec','glon','glat',
                       'ts','npred',
                       'fit1_ext_ts','fit1_halo_ts',
                       'fit_ext_ts','fit_halo_ts',
    #                   'fit1_dlike','fit_dlike',
                       'fit_nsrc',
                       'fit_dlike_ext','fit_dlike_halo','fit_mean_sep']

    subs = {}

    for row in tab:
        subs[row['codename'].replace('+','p').replace('.','_')] = \
            r'<a href="%s/index.html">%s</a>'%(row['codename'],row['name'])

    tab['fit1_ext_ts'] = tab['fitn_ext_ts'][:,1]
    tab['fit1_halo_ts'] = tab['fitn_halo_ts'][:,1]
    tab['fit1_dlike'] = tab['fitn_dlike'][:,1]

    tab2 = tab[html_table_cols]

    #max_ts = np.max(tab['fit1_halo_ts'][:,:,2],axis=1)
    #max_ts[max_ts<0] = 0    
    #tab2['fit1_halo_2.00_ts'] = max_ts
    #tab2['fit1_halo_2.00_ts'].format='%.2f'

    #max_ts = np.max(tab['fit_halo_ts'][:,:,2],axis=1)
    #max_ts[max_ts<0] = 0   
    #tab2['fit_halo_2.00_ts'] = max_ts
    #tab2['fit_halo_2.00_ts'].format='%.2f'

    tab2.columns['ts'].format='%.2f'
    tab2.columns['npred'].format='%.2f'
    tab2.columns['ra'].format='%.2f'
    tab2.columns['dec'].format='%.2f'
    tab2.columns['glon'].format='%.2f'
    tab2.columns['glat'].format='%.2f'
    tab2.columns['fit1_ext_ts'].format='%.2f'
    tab2.columns['fit1_halo_ts'].format='%.2f'
    #tab2.columns['fit1_dlike'].format='%.2f'
    tab2.columns['fit_ext_ts'].format='%.2f'
    tab2.columns['fit_halo_ts'].format='%.2f'
    #tab2.columns['fit_dlike'].format='%.2f'
    tab2.columns['fit_dlike_ext'].format='%.2f'
    tab2.columns['fit_dlike_halo'].format='%.2f'
    tab2.columns['fit_nsrc'].format='%.2f'
    tab2.columns['fit_mean_sep'].format='%.2f'

    tab2.columns['linkname'].name = 'name'
    tab2.write(args.output,format='html')
    with open(args.output, 'r') as myfile:
        data=myfile.read()
    with open(args.output, 'w') as myfile:
        myfile.write(data.format(**subs))


    scripts = [
        #{'name' : 'script', 'attrib' : {'src' : 'http://code.jquery.com/jquery-1.12.0.min.js'}},
        #{'name' : 'script', 'attrib' : {'src' : 'http://cdn.datatables.net/1.10.10/js/jquery.dataTables.min.js'}},
        {'name' : 'script', 'attrib' : { 'type' : 'text/javascript', 'src' : "https://cdn.datatables.net/r/bs-3.3.5/jqc-1.11.3,dt-1.10.8/datatables.min.js"} },
        {'name' : 'script', 'attrib' : {'src' : 'table.js'}},
        {'name' : 'link', 'attrib' : {'rel' : 'stylesheet', 'type' : 'text/css',
                                      'href' : 'https://cdn.datatables.net/r/bs-3.3.5/jq-2.1.4,dt-1.10.8/datatables.min.css'} },
        {'name' : 'link', 'attrib' : {'rel' : 'stylesheet', 'type' : 'text/css',
                                      'href' : 'style.css'} }
         ]

    root = html.parse(args.output)
    path = root.xpath('head')[0]

    elements = []

    e0 = html.Element('div')
    e0.attrib['class'] = 'container-fluid'


    for s in scripts:

        e = html.Element(s['name'])
        for attrib_name, attrib_val in s['attrib'].items():
            e.attrib[attrib_name] = attrib_val
        e.text = ''
        path.insert(-1,e)
    #    e0.insert(-1,e)

    path0 = root.xpath('body')[0]
    path0.insert(-1,e0)

    path1 = root.xpath('body/table')[0]

    path1.attrib['id']='mytable'
    path1.attrib['class']='table table-striped table-bordered table-hover'

    #path1.insert(0,e0)
    e0.append(path1)

    root.write(args.output)


if __name__ == "__main__":
    main()
