Adjacency Matrices
~~~~~~~~~~~~~~~~~~~~~

This shell script is located in ``fmri_preproc/resources/misc/gen_ts_mat.sh``.

This script applies pre-computed linear and/or non-linear transforms to template space 
and extracts the mean timeseries from each ROI (supplied in the label file).

.. code-block:: text

    Usage: 
        
        gen_ts_mat.sh -f <nifti> -o <nifti> --label <nifti> --lin-xfm <mat> [--no-cleanup]
        gen_ts_mat.sh -f <nifti> -o <nifti> --label <nifti> --nl-xfm <warp> [--no-cleanup]

    Required arguments
      -f, --func                      Input 4D functional image file
      -o, --out-prefix                Output prefix name
      --label                         Input atlas label file

    Mutually exclusive arguments
      --lin-xfm                       Input LINEAR transformation matrix file for the standard -> native (func) transform
      --nl-xfm                        Input NON-LINEAR transformation warp file for the standard -> native (func) transform
    
    Optional arguments
      --no-cleanup                    DO NOT perform clean-up of temporary working directories
      -h, -help, --help               Prints the help menu, then exits.

