Pre-MCDC
~~~~~~~~~~~~~~~~~

The ``pre-mcdc`` subcommand menu contains all the options needed for the pre- motion correction, distortion correction
step of the preprocessing pipeline.

.. code-block:: text
    
    usage: fmri_preproc.py pre-mcdc [-h] [-s <str>] [-r <str>] [-o <str>] [-e <str>] [-version]
                                    [--write-config <file>] [-v] [--log-level <str>] [--config <file>]
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
