# -*- coding: utf-8 -*-
"""Performs ICA-based denoising of fMRI data for the ``fmri_preproc`` preprocessing pipeline.

NOTE: 
    External dependency: FIX - https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide.
"""
import os
import pandas as pd
import numpy as np

from typing import (
    Optional,
    Tuple
)

from fmri_preproc.utils.util import Command
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.func.ica import ica
from fmri_preproc.func.mcdc import brain_extract

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

from fmri_preproc.utils.fslpy import (
    applyxfm,
    invxfm,
    fslmaths,
)

def fix_extract(func_filt: str,
                func_ref: str,
                struct: str,
                struct_brainmask: str,
                struct_dseg: str,
                dseg_type: str,
                func2struct_mat: str,
                mot_param: str,
                outdir: str,
                icadir: Optional[str] = None,
                func_brainmask: Optional[str] = None,
                func_tr: Optional[float] = None,
                temporal_fwhm: Optional[float] = 150,
                log: Optional[LogFile] = None
               ) -> None:
    """FIX feature extraction."""
    with NiiFile(src=func_filt, assert_exists=True, validate_nifti=True) as ff:
        with NiiFile(src=func_ref, assert_exists=True, validate_nifti=True) as fr:
            with NiiFile(src=struct, assert_exists=True, validate_nifti=True) as st:
                with NiiFile(src=struct_brainmask ,assert_exists=True, validate_nifti=True) as sb:
                    with NiiFile(src=struct_dseg, assert_exists=True, validate_nifti=True) as sd:
                        func_filt: str = ff.abspath()
                        func_ref: str = fr.abspath()
                        struct: str = st.abspath()
                        struct_brainmask: str = sb.abspath()
                        struct_dseg: str = sd.abspath()
    
    with File(src=func2struct_mat, assert_exists=True) as fm:
        with File(src=mot_param, assert_exists=True) as mt:
            func2struct_mat: str = fm.abspath()
            mot_param: str = mt.abspaths()
    
    with WorkDir(src=outdir) as od:
        denoisedir: str = os.path.join(od.src,'denoise')
        fixdir: str = os.path.join(denoisedir,'fix')
        with WorkDir(src=fixdir) as fx:
            fx.mkdir()
            denoisedir: str = os.path.abspath(denoisedir)
            fixdir: str = fx.abspath()
    
    # Setup fake FIX directory
    with NiiFile(src=func_filt) as ff:
        tmp_func: str = os.path.join(fixdir,'filtered_func_data.nii.gz')
        tmp_func: str = ff.sym_link(dst=tmp_func, relative=True)
    
    if icadir:
        with WorkDir(src=icadir) as icd:
            tmp_ica: str = os.path.join(fixdir,'filtered_func_data.ica')
            tmp_ica: str = icd.sym_link(dst=tmp_ica, relative=True)
    else:
        icadir: str = os.path.join(fixdir,'filtered_func_data.ica')
        with WorkDir(src=icadir) as icd:
            icd.mkdir()
        
        _: Tuple[str] = ica(outdir=icadir,
                            func=tmp_func,
                            func_brainmask=func_brainmask,
                            func_tr=func_tr,
                            temporal_fwhm=temporal_fwhm,
                            log=log)

    mcdir: str = os.path.join(fixdir,'mc')
    parfile: str = os.path.join(mcdir,'prefiltered_func_data_mcf.par')
    with WorkDir(src=mcdir) as mcd:
        mcd.mkdir()
    
    mp: pd.DataFrame = pd.read_csv(mot_param, delimiter='\t', index_col=None)
    np.savetxt(parfile, 
        mp[['RotX', 'RotY', 'RotZ', 'X', 'Y', 'Z']].values, 
        fmt='%0.6f',
    )

    mean_func: str = os.path.join(fixdir,'mean_func.nii.gz')
    mean_func: str = fslmaths(img=func_filt).Tmean().run(out=mean_func, log=log)

    fixregdir: str = os.path.join(fixdir, 'reg')
    with WorkDir(src=fixregdir) as fr:
        fr.mkdir()
    
    example_func: str = os.path(fixdir, 'example_func.nii.gz')
    highres: str = os.path(fixdir, 'highres.nii.gz')
    highres2examp: str = os.path(fixdir, 'highres2example_func.nii.gz')
    funcmask: str = os.path(fixdir, 'mask.nii.gz')

    with File(src=func_ref) as fr:
        example_func: str = fr.sym_link(dst=example_func, relative=True)
    
    highres: str = fslmaths(img=struct).mul(struct_brainmask).run(out=highres, log=log)
    highres2examp: str = invxfm(inmat=func2struct_mat,outmat=highres2examp, log=log)
    funcmask: str = applyxfm(src=struct_brainmask, ref=example_func, mat=highres2examp, out=funcmask)
    funcmask: str = fslmaths(img=funcmask).thr(0.5).bin().run(out=funcmask, log=log)

    # TODO:
    #   mask convert function

    # Extract FIX features
    cmd: Command = Command("fix")
    cmd.opt("-f")
    cmd.opt(f"{fixdir}")
    cmd.run(log=log)

    return fixdir

def _classify(fixdir: str,
              rdata: str,
              thr: int,
              fix_src: Optional[str],
              log: Optional[LogFile] = None
             ) -> str:
    """Helper function that performs unlabeled ICA classification.
    """
    thr: int = int(thr)

    with File(src=rdata, assert_exists=True) as rd:
        rdata: str = rd.abspath()
    
    with WorkDir(src=fixdir) as fd:
        if not fd.exists():
            raise FileNotFoundError(f"Input FIX directory does not exist.")
        
        fixdir: str = fd.abspath()
        fix_log: str = os.path.join(fd.abspath(),'.fix_2b_predict.log')
    
    if fix_src:
        cmd: Command = Command(f"{fix_src}")
    else:
        cmd: Command = Command("fix")
    
    cmd.opt("-c")
    cmd.opt(f"{fixdir}")
    cmd.opt(f"{rdata}")
    cmd.opt(f"{thr}")

    try:
        cmd.run(log=log)
    except Exception as _:
        with open(fix_log, 'r') as f:
            s: str = f.read()
        raise RuntimeError(s)

    return fixdir



def fix_classify():
    pass

def fix_apply():
    pass
