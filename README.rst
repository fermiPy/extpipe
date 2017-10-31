Extension Pipeline
==================

Creating Analysis Directories
-----------------------------

Create directories for standard analysis:

.. code-block:: bash

   $ python scripts/setup_config.yaml --config=std_psf0123_joint2a --script=run-region-analysis \ 
   haloanalysis/sourcelists/3fgl_srcs_list_glat050.yaml haloanalysis/sourcelists/3fhl_srcs_list_glat050.yaml


Create directories for alternative IEMs:

.. code-block:: bash

   $ python scripts/setup_config.yaml --config=std_psf0123_joint2a --script=run-region-analysis \
   --model=all --script=run-region-analysis haloanalysis/sourcelists/fhes_list.yaml &> /dev/null &
