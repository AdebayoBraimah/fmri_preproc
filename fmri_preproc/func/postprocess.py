# -*- coding: utf-8 -*-
"""Perform post-processing spatial smoothing and intensity norm for the ``fmri_preproc`` fMRI pre-processing pipeline.
"""
import os
import numpy as np
import nibabel as nib

from typing import (
    Dict, 
    Optional,
    Tuple
)

from fmri_preproc.utils.command import Command
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.fileio import NiiFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.util import timeops


from fmri_preproc.utils.fslpy import fslmaths


# Globlally define (temporary) log file object
with TmpDir(src=os.getcwd()) as tmpd:
    with TmpDir.TmpFile(tmp_dir=tmpd.src, ext='.log') as tmpf:
        log: LogFile = LogFile(log_file=tmpf.src)


@timeops(log)
def postprocess(func: str,
                func_mean: str,
                func_brainmask: str,
                outdir: str,
                spatial_fwhm: Optional[float] = 5,
                intnorm: bool = False,
                basename: Optional[str] = None,
                log: Optional[LogFile] = None
               ) -> Tuple[str,str]:
    """Perform post-processing spatial smoothing and intensity normalization. 
    """
    if isinstance(spatial_fwhm, int) or isinstance(spatial_fwhm, float):
        pass
    else:
        if log: log.error(f"RuntimeError: Input spatial smoothing kernel is neither an integer nor a float: {spatial_fwhm}")
        raise RuntimeError(f"Input spatial smoothing kernel is neither an integer nor a float: {spatial_fwhm}")

    with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
        with NiiFile(src=func_mean, assert_exists=True, validate_nifti=True) as fm:
            with NiiFile(src=func_brainmask, assert_exists=True, validate_nifti=True) as fnb:
                func: str = fn.abspath()
                func_mean: str = fm.abspath()
                func_brainmask: str = fnb.abspath()

    if log: log.log("Post-processing: spatial smoothing and/or intensity normalization")

    with WorkDir(src=outdir) as od:
        outdir: str = od.abspath()

    if basename:
        pass
    else:
        with NiiFile(src=func) as fn:
            _, basename, _ = fn.file_parts()
    
    # Define outputs
    outputs: Dict[str,str] = {
                                "func_smooth": os.path.join(outdir,f'{basename}_{spatial_fwhm}mm_smooth.nii.gz'),
                                "func_intnorm": os.path.join(outdir,f'{basename}_intnorm.nii.gz')
                             }
    
    # FEAT-style spatial smoothing
    spatial_fwhm: float = float(spatial_fwhm)

    func_smooth: str = outputs.get('func_smooth')
    func_smooth: str = fslmaths(img=func).mas(func_brainmask).run(out=func_smooth, log=log)

    d: nib.Nifti1Image = nib.load(func_smooth).get_data().astype(float)
    func_q2, func_q98 = np.percentile(d, [2, 98])

    d: nib.Nifti1Image = nib.load(func_smooth).get_data().astype(float)
    m: nib.Nifti1Image = nib.load(func_brainmask).get_data().astype(float)
    median_intensity: float = np.percentile(d[m >= 1], 50)

    if spatial_fwhm > 0.01:
        #  Spatial smoothing using a Gaussian kernel of FWHM filter_width_mm
        if log: log.log(f"Performing spatial smoothing with spatial smoothing kernel (mm): {spatial_fwhm}")

        sigma: float = spatial_fwhm / 2.355
        susan_thr: float = (median_intensity - func_q2) * 0.75

        cmd: Command = Command("susan")

        cmd.opt(f"{func_smooth}")
        cmd.opt(f"{susan_thr}")
        cmd.opt(f"{sigma}")
        cmd.opt(f"3")
        cmd.opt(f"1")
        cmd.opt(f"1")
        cmd.opt(f"{func_mean}")
        cmd.opt(f"{func_smooth}")
        cmd.opt(f"{func_smooth}")

        cmd.run(log=log)

        func_smooth: str = fslmaths(img=func_smooth).mas(func_brainmask).run(out=func_smooth)
        func: str = func_smooth
    
    # FEAT-style intensity normalization
    func_intnorm: str = outputs.get('func_intnorm')
    normmean: int = 10_000

    if intnorm:
        # Multiplicative mean intensity normalization of the volume at each
        #   timepoint
        if log: log.log(f"Performing intensity normalization with a multiplicative mean intensity of: {normmean}")
        func_intnorm: str = fslmaths(img=func).inm(normmean).run(out=func_intnorm, log=log)
    else:
        # Grand-mean intensity normalization of the entire 4D dataset by a
        #   single multiplicative factor
        if log: log.log(f"Performing grand-mean intensity normalization")
        scaling: float = normmean / median_intensity
        func_intnorm: str = fslmaths(img=func).mul(scaling).run(out=func_intnorm, log=log)
    
    return func_smooth, func_intnorm
