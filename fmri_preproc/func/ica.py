# -*- coding: utf-8 -*-
"""Performs single-subject ICA for the ``fmri_preproc`` fMRI pipeline.
"""
import os
import nibabel as nib

from typing import (
    Dict, 
    Optional,
    Tuple
)

from fmri_preproc.utils.outputs.ica import ICAFiles
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.fileio import NiiFile

from fmri_preproc.utils.fslpy import (
    fslmaths,
    melodic
)

def ica(outdir: str,
        func: str,
        func_brainmask: str,
        func_tr: Optional[float] = None,
        temporal_fwhm: Optional[float] = 150,
        icadim: Optional[int] = None,
        log: Optional[LogFile] = None
       ) -> Tuple[str]:
    """Performs single-subject ICA.
    """
    # Validate input NIFTI files
    with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
        with NiiFile(src=func_brainmask, assert_exists=True, validate_nifti=True) as fb:
            func: str = fn.abspath()
            func_brainmask: str = fb.abspath()

    # Define outputs
    with WorkDir(src=outdir) as od:
        ica_out: ICAFiles = ICAFiles(outdir=od.abspath())
        outputs: Dict[str,str] = ica_out.outputs()

        icadir: str = outputs.get('icadir')
        meldir: str = outputs.get('meldir')
        func_filt: str = outputs.get('func_filt')

    if func_tr:
        func_tr: float = float(func_tr)
    else:
        func_tr: float = nib.load(func).header.get('pixdim','')[4]

    if temporal_fwhm:
        temporal_fwhm: float = float(temporal_fwhm)
        # sigma: float = np.round(temporal_fwhm / (2 * func_tr))  # 150 second filter

        with TmpDir(src=icadir) as tmp:
            func_mean: str = tmp.join('func_mean.nii.gz')
            func_mean: str = fslmaths(img=func).Tmean().run(out=func_mean, log=log)
            func_filt: str = fslmaths(img=func).bptf(high_pass=temporal_fwhm, 
                                                     low_pass=-1, 
                                                     tr=func_tr, 
                                                     input_is_sec=True).add(input=func_mean).run(out=func_filt, log=log)
    else:
        with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
            # func_filt: str = fn.copy(dst=func_filt)
            func_filt: str = fn.sym_link(dst=func_filt, relative=True)
    
    meldir: str = melodic(input=func_filt,
                          mask=func_brainmask,
                          outdir=meldir,
                          dim=icadim,
                          tr=func_tr,
                          log=log)
    return func_filt, meldir
