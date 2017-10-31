Extension Pipeline
==================

Creating Analysis Directories
-----------------------------

Create directories manually:

.. code-block:: bash

   $ fermipy-clone-configs --source_list=source_list.yaml --script=run-region-analysis --basedir=v0 config.yaml

If it's not working, make sure you ran 

.. code-block:: bash

   $ python setup.py develop 

in the haloanalysis directory. 

To send things to the batch farm, do 

.. code-block:: bash

   $ fermipy-dispatch --runscript=run-region-analysis.sh 

Create directories for standard analysis:

.. code-block:: bash

   $ python scripts/setup_config.py --config=std_psf0123_joint2a --script=run-region-analysis \ 
   haloanalysis/sourcelists/3fgl_srcs_list_glat050.yaml haloanalysis/sourcelists/3fhl_srcs_list_glat050.yaml
   $ python scripts/setup_config.py --config=std_psf0123_joint2a --script=run-halo-analysis \ 
   haloanalysis/sourcelists/3fgl_srcs_list_glat050.yaml haloanalysis/sourcelists/3fhl_srcs_list_glat050.yaml


Create directories for alternative IEMs:

.. code-block:: bash

   $ python scripts/setup_config.py --config=std_psf0123_joint2a --script=run-region-analysis \
   --model=all --script=run-region-analysis haloanalysis/sourcelists/fhes_list.yaml
   $ python scripts/setup_config.py --config=std_psf0123_joint2a --script=run-halo-analysis \
   --model=all --script=run-region-analysis haloanalysis/sourcelists/fhes_list.yaml
