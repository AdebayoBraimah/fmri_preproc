Import Data
~~~~~~~~~~~~~~~~~

The ``import-data`` subcommand menu contains all the options needed for importing resting-state fMRI data into the
pipeline.

.. note:: 

  ``--sp-??`` options pairs are mutually exclusive, **e.g.** specifying ``--sp-PA`` and ``--sp-LR`` would 
  raise an :mod:`AttributeError` exception.

.. code-block:: text

    usage: fmri_preproc.py import-data [-h] [-s <str>] [-r <str>] [-o <str>] [-e <str>] [-version]
                                      [--write-config <file>] [-v] [--log-level <str>]
                                      [--config <file>] [-a <age>] [--birth-age <age>] [--T2w <file>]
                                      [--T2w-brainmask <file>] [--T2w-tissue-seg <file>]
                                      [--tissue-seg-type <str>] [--T1w <file>] [--wmmask <file>]
                                      [-f <file>] [--func-echospacing <float>] [--func-pedir <PE-dir>]
                                      [--func-slorder <file>] [--func-brainmask <file>]
                                      [--func-inplane-accel <float>] [--func-mbfactor <int>]
                                      [--sbref <file>] [--sbref-brainmask <file>]
                                      [--sbref-echospacing <float>] [--sbref-pedir <PE-dir>]
                                      [--spinecho <file>] [--sp-echospacing <float>]
                                      [--sp-pedir <PE-dir,PE-dir,..,PE-dir>]
                                      [--sp-inplane-accel <float>] [--sp-AP <file>] [--sp-PA <file>]
                                      [--sp-LR <file>] [--sp-RL <file>] [--sp-IS <file>]
                                      [--sp-SI <file>]
                                      {import-data,pre-mcdc,mcdc,post-mcdc,run-all} ...

    positional arguments:
      {import-data,pre-mcdc,mcdc,post-mcdc,run-all}  Type the 'subcommand' name followed by '-h' or '--
                                                    help' for more information.
        import-data                                  Imports the subject's imaging data into a working
                                                    directory for image preprocessing.
        pre-mcdc                                     Pre-Motion Correction and Distortion Correction
                                                    (MCDC).
        mcdc                                         Motion Correction and Distortion Correction (MCDC).
        post-mcdc                                    Post-MCDC preprocessing, which includes: func ->
                                                    template registration, ICA-based denoising, and QC
                                                    (web page).
        run-all                                      Perform/Run all preprocessing steps of the
                                                    pipeline.

    optional arguments:
      -h, --help                                     show this help message and exit
      -s <str>, --sub <str>                          REQUIRED: Subject ID.
      -r <str>, --run <str>                          REQUIRED: Run number.
      -o <str>, --outdir <str>                       REQUIRED: Output parent directory.
      -e <str>, --ses <str>                          OPTIONAL: Session ID.
      -version, --version                            Prints version, then exits.
      --write-config <file>                          Writes a template configuration (JSON) to file.
      -v, --verbose                                  Enables verbose output to the command and log
                                                    files.
      --log-level <str>                              Log level to be written to output log files.
                                                    Acceptable inputs include: 'info', 'debug',
                                                    'critical', 'warning', and 'error'
                                                    [default='info'].
      --config <file>, --settings <file>             Configuration (JSON) file that contains additional
                                                    input arguments not specified at the command line.
      -a <age>, --age <age>                          REQUIRED: Age at scan. Acceptable inputs include:
                                                    'neo', 28 - 44 (PMA, in wks), and 1 or 2 (yrs).
      --birth-age <age>                              Gestational age at birth.
      --T2w <file>                                   REQUIRED: Input (bias corrected, e.g. '_restore')
                                                    T2w image data.
      --T2w-brainmask <file>                         Input T2w brainmask image data (_brainmask_bet).
      --T2w-tissue-seg <file>                        REQUIRED: Input Draw-EM or other neonatal brain
                                                    tissue segmentation (e.g. '_drawem_tissue_labels')
                                                    image data. NOTE: Other brain tissue files can be
                                                    used (e.g. from FreeSurfer or FSL's FAST).
      --tissue-seg-type <str>                        REQUIRED: Input brain tissue segmentation type.
                                                    Acceptable inputs include: 'drawem',
                                                    'freesurfer_aseg', or 'fsl_fast'.
      --T1w <file>                                   Input (bias corrected, e.g. '_restore') T1w image
                                                    data.
      --wmmask <file>                                WM mask image file.
      -f <file>, --func <file>                       REQUIRED: Input rs-fMRI data. NOTE: IF input data
                                                    conforms to the BIDS naming convention (e.g.
                                                    'sub-001[_ses-001]_task-
                                                    rest[_run-01]_bold.nii.gz'), then '--sub', '--ses',
                                                    and '--run' need not be specified.
      --func-echospacing <float>                     REQUIRED: Functional echospacing. NOTE: If unsure
                                                    use 0.0005.
      --func-pedir <PE-dir>                          REQUIRED: Phase encoding direction. Acceptable
                                                    inputs include: 'PA', 'AP', 'RL', 'LR', 'IS', or
                                                    'SI'.
      --func-slorder <file>                          Functional slice order specification file. This
                                                    file specifies the slice acquisition order of the
                                                    input rs-fMRI. If not provided one can be created
                                                    internally IF '--func-mbfactor' is specified.
      --func-brainmask <file>                        Functional brain mask.
      --func-inplane-accel <float>                   Functional inplane acceleration (e.g. SENSE factor
                                                    on Philips MR scanners or GRAPPA, the general term)
                                                    [default=1.00].
      --func-mbfactor <int>                          Functional multi-band acceleration factor.
      --sbref <file>                                 Single band reference MR image data.
      --sbref-brainmask <file>                       Single band reference MR brain mask image data.
      --sbref-echospacing <float>                    Single band reference MR image echospacing. If not
                                                    specified, then it is assumed to match '--func-
                                                    echospacing'.
      --sbref-pedir <PE-dir>                         Single band reference MR image data phase encoding
                                                    direction. NOTE: This should match '--func-pedir'.
                                                    If not specified, then it is assumed to match '--
                                                    func-pedir'. If this is not case, this could result
                                                    in mis-registration and distortion correction
                                                    issues.
      --spinecho <file>                              Spinecho MR image data.
      --sp-echospacing <float>                       Spinecho MR image echospacing. If not specified,
                                                    then it is assumed to match '--func-echospacing'.
      --sp-pedir <PE-dir,PE-dir,..,PE-dir>           Spinecho MR image phase encoding direction(s),
                                                    specified as a comma separated list. Accepatable
                                                    inputs are the same as those specifed for both '--
                                                    func-pedir', and '--sbref-pedir'.
      --sp-inplane-accel <float>                     Spinecho MR image inplane acceleration . If not
                                                    specified, then it is assumed to match '--func-
                                                    inplane-accel'.
      --sp-AP <file>                                 A -> P (anterior to posterior) Spinecho MR image
                                                    data.
      --sp-PA <file>                                 P -> A (posterior to anterior) Spinecho MR image
                                                    data.
      --sp-LR <file>                                 L -> R (left to right) Spinecho MR image data.
      --sp-RL <file>                                 R -> L (right to left) Spinecho MR image data.
      --sp-IS <file>                                 I -> S (inferior to superior) Spinecho MR image
                                                    data.
      --sp-SI <file>                                 S -> I (superior to inferior) Spinecho MR image
                                                    data.
