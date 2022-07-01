Post-MCDC
~~~~~~~~~~~~~~~~~

The ``post-mcdc`` subcommand menu contains all the options needed for the remainder of the pipeline which includes:

  - Registration to the age-matched template(s) (, which also includes templates from the `UNC AAL atlas <https://www.med.unc.edu/bric/ideagroup/software/unc-infant-0-1-2-atlases/>`_) via `ANTs <http://stnava.github.io/ANTs/>`_ and `C3D <https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md>`_
  - ICA-based cleanup via `FSL's FIX <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide#Running_FIX>`_
  - Spatial smoothing and intensity normalization (of a 3D volume, not the entire timeseries).

.. warning:: 

  Be sure that ``ANTSPATH`` has been set in your command line environment.

.. code-block:: text
    
    usage: fmri_preproc.py post-mcdc [-h] [-s <str>] [-r <str>] [-o <str>] [-e <str>] [-version]
                                    [--write-config <file>] [-v] [--log-level <str>] [--config <file>]
                                    [--standard-age <int>] [--template-ages <str,str,..str>]
                                    [--atlas-dir <path>] [--temporal-fwhm <float>] [--ica-dim <int>]
                                    [--rdata <str>] [--fix-threshold <int>] [--smooth <float>]
                                    [--int-norm]
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
      --standard-age <int>                           Standard age to use for final dHCP age-matched
                                                    template (PMA, in weeks) [default=40 (wks)].
      --template-ages <str,str,..str>                Comma separated list of ages. Acceptable inputs
                                                    include: 'neo', 1, 2 (yrs), and/or 28 - 44 (PMA,
                                                    wks).
      --atlas-dir <path>                             Atlas directory path.
      --temporal-fwhm <float>                        Temporal FWHM filtering coefficient [default=150.0
                                                    (s)].
      --ica-dim <int>                                Maximal ICA dimensionality.
      --rdata <str>                                  Trained FIX classifier.
      --fix-threshold <int>                          FIX noise classification threshold. [default=10].
      --smooth <float>                               Smoothing kernel (FWHM, in mm) [default: disabled].
      --int-norm                                     Performs intensity normalization OR grand-mean
                                                    scaling IF '--smooth' is specified and greater than
                                                    0 [default: disabled].
