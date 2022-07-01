MCDC
~~~~~~~~~~~~~~~~~

The ``mcdc`` subcommand menu contains all the options needed for the motion correction, 
distortion correction steps of the preprocessing pipeline.

.. warning:: 

  This stage of the pipeline **REQUIRES** a GPU.
  Running this stage of the pipeline without a GPU will likely result in either:
  
    - The pipeline failing as a :meth:`DependencyError` exception is raised (this happens with versions of FSL <6.0.4).
    - The pipeline will have an *excessively* long runtime (so long in fact, I have no idea how long it would take - and this usually happens with versions of FSL 6.0.5+).

.. code-block:: text
    
    usage: fmri_preproc.py mcdc [-h] [-s <str>] [-r <str>] [-o <str>] [-e <str>] [-version]
                                [--write-config <file>] [-v] [--log-level <str>] [--config <file>]
                                [--s2v] [--dc] [--mbs] [--mcflirt] [--mporder <int>]
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
      --s2v, --slice-to-volume                       Enables slice-to-volume motion correction. NOTE:
                                                    This will LIKELY require EDDY_CUDA and GPU
                                                    processing. This is possible without GPUs in FSL
                                                    v6.0.5+ but is VERY SLOW.
      --dc, --distortion-correction                  Enables distortion correction.
      --mbs                                          Enables movement by susceptibility distortion
                                                    correction.
      --mcflirt                                      Enables FSL's MCFLIRT-based motion correction.
                                                    NOTE: Enabling this option DISABLES: '--s2v', '--
                                                    dc', and '--mbs' options.
      --mporder <int>                                The number of discrete cosine (DCT) basis sets used
                                                    to model the intra-slice movement within a volume
                                                    (used by FSL's eddy). This value should be
                                                    maximally defined as the number of rows in the '--
                                                    func-slorder' option. For reference, this value
                                                    should not exceed the number of excitations per
                                                    volume (e.g. given a MB factor of 3, with 45
                                                    acquired slices, the mporder should not exceed 15)
                                                    [NOTE: this value is automatically computed when
                                                    the '--s2v' option is enabled].
