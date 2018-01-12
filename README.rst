Extension Pipeline
==================

The extension pipeline consists of two analysis scripts:

 * run-region-analysis : Tests a given source for extension while
   iteratively adding up to 5 point sources to the inner ROI.  The
   outputs of this step are the extension analysis results and a
   sequence of up to 5 PS models.
 * run-halo-analysis : Using the PS models derived in
   run-region-analysis test the target source for an extended
   component that is superimposed on the source.

After running these you will need to aggregate the analysis results using the extpipe-aggregate script.



Creating Analysis Directories
-----------------------------

Before running the analysis you will need to create the directory
structure.  The ``fermipy-clone-configs`` script will generate a set of
analysis directories for each entry in a source list file.  To create
directories manually for the two analysis scripts:

.. code-block:: bash

   $ fermipy-clone-configs --source_list=source_list.yaml --script=run-region-analysis --basedir=v0 config.yaml
   $ fermipy-clone-configs --source_list=source_list.yaml --script=run-halo-analysis --basedir=v0 config.yaml

where ``config.yaml`` is the baseline configuration and
``source_list.yaml`` is a hash of analysis configurations keyed to the
analysis subdirectory name.

After generating the analysis directories to send things to the batch
farm run fermipy-dispatch with the root analysis directory as its
argument:

.. code-block:: bash

   $ fermipy-dispatch --runscript=run-region-analysis.sh v0

The setup_config.py script can be used to automatically the analysis
directories for a matrix of different analysis configurations. To
create directories for standard analysis:

.. code-block:: bash

   $ python scripts/setup_config.py --config=std_psf0123_joint2a --script=run-region-analysis \ 
   haloanalysis/sourcelists/3fgl_srcs_list_glat050.yaml haloanalysis/sourcelists/3fhl_srcs_list_glat050.yaml
   $ python scripts/setup_config.py --config=std_psf0123_joint2a --script=run-halo-analysis \ 
   haloanalysis/sourcelists/3fgl_srcs_list_glat050.yaml haloanalysis/sourcelists/3fhl_srcs_list_glat050.yaml


To Create directories for alternative IEMs:

.. code-block:: bash

   $ python scripts/setup_config.py --config=std_psf0123_joint2a --script=run-region-analysis \
   --model=all --script=run-region-analysis haloanalysis/sourcelists/fhes_list.yaml
   $ python scripts/setup_config.py --config=std_psf0123_joint2a --script=run-halo-analysis \
   --model=all --script=run-region-analysis haloanalysis/sourcelists/fhes_list.yaml


Running the Analysis
--------------------

First create a conda environment with the releases of fermipy and extpipe that you plan to use:

.. code-block:: bash

   $ conda create --name extanalysis-vXX python=2.7 numpy scipy matplotlib==2.0.0 astropy pytest pyyaml ipython healpy
   $ source activate extanalysis-vXX
   $ cd <path to fermipy>
   $ python setup.py install
   $ cd <path to extpipe>
   $ python setup.py install

For testing purposes you can also run the following instead of install: 

.. code-block:: bash

   $ python setup.py develop 



Aggregating Analysis Results
----------------------------

After the analysis has completed you need to run the aggregation
script to merge the results for each target into a single catalog
file.  This step of the analysis also selects the best-fit model for
each source.

.. code-block:: bash

   $ extpipe-aggregate <basedir> --batch


where ``<basedir>`` is the root analysis directory.  After this step
completes you can build a single catalog file by running
``extpipe-build-catalog``:

.. code-block:: bash

   $ extpipe-build-catalog std_psf0123_joint2a_stdmodel/*_cat.fits --output=catalog.fits --ts_threshold=16
   $ extpipe-build-catalog std_psf0123_joint2a_stdmodel/*_lnl.fits --output=catalog_lnl.fits --ts_threshold=16
