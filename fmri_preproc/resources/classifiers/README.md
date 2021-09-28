# Classifiers
---------------

## Description of parameters and units for FIX classifiers.

`x_size` (mm): X dimension length.        
`y_size` (mm): Y dimension length.        
`z_size` (mm): Z dimension length.        
`t_size`: Number of dynamics (time points) in the functional MR image timeseries.        
`x_vox` (mm): X dimension voxel length.        
`y_vox` (mm): Y dimension voxel length.        
`z_vox` (mm): Z dimension voxel length.        
`TR` (sec): Repetition time.        
`InplaneAccelerationFactor` (SENSE / GRAPPA): Acceleration factor used by the MR scanner for slice acquisition speed.        
`MultibandAccelerationFactor`: Multi-band acceleration factor used by the MR scanner for multi-slice acquisitions.        
`topup`: Topup distortion estimation of input EPI fieldmaps.        
`dseg_type`: Segmentation image type (e.g. ``drawem``, ``fsl_fast``, ``freesurfer_ase``).        
`fieldmap_correction`: Fieldmap correction and registration to structural image.        
`eddy`: Eddy current distortion estimation and correction.        
`dc`: Distortion correction.        
`mbs`: Movement-by-susceptibility estimation.        
`s2v`: Slice-to-volume motion correction.        
`mc`: Motion correction (Performed either by FSL'S ``EDDY`` or ``MCFLIRT``.).        
`icadim`: ICA dimensionality (usually set to automatic dimensionality estimation).        
`temporal_fwhm` (s): Temporal highpass filter.        
`smooth` (mm): Smoothing kernel (FWHM).        
`intnorm`: Intensity normalization.        
`NumberOfSubjects`: Number of subjects used to train the classifier.        
`MeanAge`: Mean age (in the case of neonate, PMA in wks).        