Troubleshooting
=================

``FileNotFoundError``
^^^^^^^^^^^^^^^^^^^^^^^

Should the :meth:`FileNotFoundError` exception be raised in any part of the pipeline AFTER ``import-data`` - 
then re-running all of the previous pipeline stages should *likely* help.


Pipeline failures
^^^^^^^^^^^^^^^^^^^

Sometimes certain stages of the pipeline fail for one reason or another (e.g. the ``mcdc`` stage is one that comes to mind).

Import-data failures
~~~~~~~~~~~~~~~~~~~~~~~

This failure likely occurs when any one of the required input files are not specifed.

Additionally, should all the necessary files have been specifed and the ``import-data`` stage still fails - then
this could likely be attributed to:

    - Image dimension mis-match between the two input fieldmaps/single band references.
    - Image file truncation of one or both of the input fieldmap(s).

At the moment, there is currently no fix for this issue.


Pre-MCDC failures
~~~~~~~~~~~~~~~~~~~

This could likely be caused by one or more the below issues:

    - Image file truncation of one or both of the input fieldmap(s).
    - Topup error.

At the moment, there is currently no fix for this issue.

MCDC failures
~~~~~~~~~~~~~~~~

This failure could occur for any reason (as this is the most error-prone stage of the pipeline).

Should any issues arise - the solution usually is to decrement the ``--mporder`` option.

.. note:: 

    The ``--mporder`` option should not be 0 - as this disables the slice-to-volume motion and distortion 
    correction steps that are performed by `FSL's eddy <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#A--mporder>`_.


Post-MCDC
~~~~~~~~~~~~~~~

Should failures occur in this stage of the pipeline there are several possible solutions:

    - If there is a registration failure, then ensure that `ANTs <http://stnava.github.io/ANTs/>`_ and `C3D <https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md>`_ are installed, and that ``ANTSPATH`` has been set.
    - If there is an issue with `FSL's FIX <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide#Running_FIX>`_ (which would mostly like be an :meth:`IndexError`) - then set the ``--ica-dim`` option to 150 (or whatever is temporally feasible/allowable) and decrement this value from there. Should issues still persist as the ``ica-dim`` option nears 0 - then the independent components will need to be manually classified (as either signal or noise), followed by re-running this stage of the pipeline.
    - If there is a ``FileNotFoundError`` - see above.
