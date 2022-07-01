Installation
================

The instructions below are intended for installation and first-time use of the ``fmri_preproc`` preprocessing pipeline.


Install package
~~~~~~~~~~~~~~~~~~

The ``fmri_preproc`` package can be installed via git by typing:

.. code-block:: bash

    # CLI installation instructions
    git clone https://github.com/AdebayoBraimah/fmri_preproc.git
    cd fmri_preproc
    pip install -r requirements.txt [--user]


Auxillary files
~~~~~~~~~~~~~~~~~~~~~~~~~~

Auxillary files must also be downloaded (e.g. atlases, templates, etc.) for the pipeline to run correctly.
To download the necessary files, enter the below commands into the command line:

.. code-block:: bash

    # Auxillary file download instructions
    #
    # While still in the fmri_preproc git repository
    ./fmri_preproc/utils/fetch_dhcp_resources.py [-h | --help]


The above command should begin downloading all the needed files - which should take some time.


External Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~

External dependencies and relevant information are described in the table shown below.
Relevant information includes the software, environmental variables, and so on (via external links).

The preprocessing pipeline assumes that these dependencies are met. If not, then a :meth:`DependncyError` exception will be
raised - halting execution of the pipeline.

.. csv-table:: Dependencies
    :file: tables/deps.csv
    :widths: 30, 70
    :header-rows: 1


