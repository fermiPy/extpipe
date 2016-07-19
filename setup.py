#!/usr/bin/env python

from setuptools import setup, find_packages
import os
import sys

from haloanalysis.version import get_git_version

setup(name='haloanalysis',
      version=get_git_version(),
      license='BSD',
      packages=find_packages(exclude='tests'),
      include_package_data = True, 
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: Implementation :: CPython',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Development Status :: 4 - Beta',
      ],
      entry_points= {'console_scripts': [
            'run-region-analysis = haloanalysis.scripts.region_analysis:main',
            'run-halo-analysis = haloanalysis.scripts.halo_analysis:main',
            'haloanalysis-aggregate = haloanalysis.scripts.aggregate:main',
            'haloanalysis-stack = haloanalysis.scripts.stack:main',
            'haloanalysis-make-html-table = haloanalysis.scripts.make_html_table:main',
            ]},
      install_requires=['numpy >= 1.6.1',
                        'matplotlib >= 1.4.0',
                        'scipy >= 0.14',
                        'astropy >= 1.0',
                        'pyyaml',
                        'healpy',
                        'wcsaxes'])

