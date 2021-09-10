# -*- coding: utf-8 -*-
"""Performs ICA-based denoising of fMRI data for the ``fmri_preproc`` preprocessing pipeline.

NOTE: 
    External dependency: FIX - https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FIX/UserGuide.
"""
import os
import pandas as pd
import numpy as np

from typing import (
    Dict,
    Optional,
    Tuple
)

from fmri_preproc.utils.util import timeops
from fmri_preproc.utils.command import Command
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.mask import convert_to_fsl_fast
from fmri_preproc.utils.fixlabels import loadLabelFile
from fmri_preproc.func.ica import ica

from fmri_preproc.utils.outputs.denoise import (
    FIXApply,
    FIXClassify,
    FIXExtract
)

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

from fmri_preproc.utils.fslpy import (
    applyxfm,
    bet,
    invxfm,
    fslmaths,
)


# Globlally define (temporary) log file object
with TmpDir(src=os.getcwd()) as tmpd:
    with TmpDir.TmpFile(tmp_dir=tmpd.src, ext='.log') as tmpf:
        log: LogFile = LogFile(log_file=tmpf.src)
        tmpd.rmdir()


@timeops(log)
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
            mot_param: str = mt.abspath()
    
    # Define outputs
    with WorkDir(src=outdir) as od:
        outfix: FIXExtract = FIXExtract(outdir=od.abspath())
        outputs: Dict[str,str] = outfix.outputs()
        with WorkDir(src=outputs.get('fixdir')) as fd:
            denoisedir: str = outputs.get('denoisedir')
            fixdir: str = fd.abspath()
    
    # Setup fake FIX directory
    if log: log.log("Setting up FIX directory")

    with NiiFile(src=func_filt) as ff:
        tmp_func: str = os.path.join(fixdir,'filtered_func_data.nii.gz')
        tmp_func: str = ff.sym_link(dst=tmp_func, relative=True)
    
    if icadir:
        with WorkDir(src=icadir) as icd:
            tmp_ica: str = os.path.join(fixdir,'filtered_func_data.ica')
            tmp_ica: str = icd.sym_link(dst=tmp_ica, relative=True)
    else:
        icadir: str = os.path.join(fixdir,'filtered_func_data.ica')
        with WorkDir(src=icadir) as _:
            pass
        
        if func_brainmask:
            with NiiFile(src=func_brainmask, assert_exists=True, validate_nifti=True) as fnb:
                func_brainmask: str = fnb.abspath()
        else:
            tmp_func_brain: str = os.path.join(fixdir, 'func_brain.nii.gz')
            tmp_func_brainmask: str = os.path.join(fixdir, 'func_brain_mask.nii.gz')

            _, func_brainmask = brain_extract(img=tmp_func,
                                              out=tmp_func_brain,
                                              mask=tmp_func_brainmask,
                                              robust=True,
                                              frac_int=0.4,
                                              log=log)
        _: Tuple[str] = ica(outdir=icadir,
                            func=tmp_func,
                            func_brainmask=func_brainmask,
                            func_tr=func_tr,
                            temporal_fwhm=temporal_fwhm,
                            log=log)

    mcdir: str = os.path.join(fixdir,'mc')
    parfile: str = os.path.join(mcdir,'prefiltered_func_data_mcf.par')
    with WorkDir(src=mcdir) as _:
        pass
    
    mp: pd.DataFrame = pd.read_csv(mot_param, delimiter='\t', index_col=None)
    np.savetxt(parfile, 
        mp[['RotX', 'RotY', 'RotZ', 'X', 'Y', 'Z']].values, 
        fmt='%0.6f',
    )

    mean_func: str = os.path.join(fixdir,'mean_func.nii.gz')
    mean_func: str = fslmaths(img=func_filt).Tmean().run(out=mean_func, log=log)

    fixregdir: str = os.path.join(fixdir, 'reg')
    with WorkDir(src=fixregdir) as _:
        pass
    
    example_func: str = os.path.join(fixregdir, 'example_func.nii.gz')
    highres: str = os.path.join(fixregdir, 'highres.nii.gz')
    highres2examp: str = os.path.join(fixregdir, 'highres2example_func.mat')
    funcmask: str = os.path.join(fixdir, 'mask.nii.gz')

    with File(src=func_ref) as fr:
        example_func: str = fr.sym_link(dst=example_func, relative=True)
    
    highres: str = fslmaths(img=struct).mul(struct_brainmask).run(out=highres, log=log)
    highres2examp: str = invxfm(inmat=func2struct_mat, outmat=highres2examp, log=log)
    funcmask: str = applyxfm(src=struct_brainmask, ref=example_func, mat=highres2examp, out=funcmask)
    funcmask: str = fslmaths(img=funcmask).thr(0.5).bin().run(out=funcmask, log=log)

    fsl_labels: str = os.path.join(fixregdir,'highres_pveseg.nii.gz')
    fsl_labels: str = convert_to_fsl_fast(dseg=struct_dseg,
                                          seg_type=dseg_type,
                                          out=fsl_labels)

    # Extract FIX features
    if log: log.log("Performing FIX feature extraction")

    cmd: Command = Command("fix")
    cmd.opt("-f")
    cmd.opt(f"{fixdir}")
    cmd.run(log=log, raise_exc=False)

    return fixdir


def _classify(fixdir: str,
              rdata: str,
              thr: int,
              fix_src: Optional[str] = None,
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
            if log: log.error(s)
        raise RuntimeError(s)

    return fixdir


@timeops(log)
def fix_classify(rdata: str,
                 thr: int,
                 outdir: str,
                 log: Optional[LogFile] = None
                ) -> Tuple[str,str]:
    """FIX feature classification.
    """
    # Check inputs
    with File(src=rdata, assert_exists=True) as rd:
        rdata: str = rd.abspath()
    
    with WorkDir(src=outdir) as od:
        denoisedir: str = od.join('denoise')
        fixdir: str = os.path.join(denoisedir,'fix')
        with WorkDir(src=fixdir) as fd:
            if not fd.exists():
                log.error("RuntimeError: FIX feature extraction must be performed before feature classification.")
                raise RuntimeError(f"FIX feature extraction must be performed before feature classification.")
            fixdir: fd.abspath()
    
    # Define outputs
    outfix: FIXClassify = FIXClassify(outdir=outdir)
    outputs: Dict[str,str] = outfix.outputs()
    
    # FIX feature classification
    if log: log.log("Performing FIX feature classification.")
    _: str = _classify(fixdir=fixdir, rdata=rdata, thr=thr)

    # Rename FIX labels
    with File(src=rdata) as rd:
        _, basename, _ = rd.file_parts()
        labels: str = f'fix4melview_{basename}_thr{thr}.txt'
        labels: str = os.path.join(fixdir, labels)
    
    with File(src=labels, assert_exists=True) as f:
        fix_labels: str = f.copy(dst=outputs.get('fix_labels'))
    
    # Create FIX regressors
    mel_mix_dir: str = os.path.join(fixdir, 'filtered_func_data.ica/melodic_mix')
    fix_reg: str = outputs.get('fix_regressors')

    noise_idx: str = loadLabelFile(filename=fix_labels, returnIndices=True)[2]
    noise_idx: np.array = np.array(noise_idx) - 1
    mix: np.array = np.loadtxt(mel_mix_dir)[:, noise_idx]
    df: pd.DataFrame = pd.DataFrame(data=mix, columns=[f'noise_{i}' for i in noise_idx])
    df.to_csv(fix_reg, sep='\t', index=None)
    
    return fix_labels, fix_reg


@timeops(log)
def fix_apply(outdir: str,
              temporal_fwhm: Optional[float] = 150.0,
              log: Optional[LogFile] = None
             ) -> None:
    """FIX regress noise ICs from data.
    """
    with WorkDir(src=outdir) as od:
        denoisedir: str = os.path.join(od.src,'denoise')
        fixdir: str = os.path.join(denoisedir,'fix')
        with WorkDir(src=fixdir) as fx:
            denoisedir: str = os.path.abspath(denoisedir)
            fixdir: str = fx.abspath()

    # Define outputs
    outfix: FIXApply = FIXApply(outdir=outdir)
    outputs: Dict[str,str] = outfix.outputs()
    
    # Output variables
    labels: str = outputs.get('fix_labels')
    fix_clean: str = outputs.get('fix_clean')
    func_clean: str = outputs.get('func_clean')
    func_mean: str = outputs.get('func_clean_mean')
    func_stdev: str = outputs.get('func_clean_std')
    func_tsnr: str = outputs.get('func_clean_tsnr')

    # FIX apply
    if log: log.log("Performing FIX noise/nuissance regression")
    cmd: Command = Command("fix")
    cmd.opt("-a")
    cmd.opt(f"{labels}")
    cmd.opt("-m")
    cmd.opt("-h")
    cmd.opt(f"{temporal_fwhm}")
    cmd.run(log=log)

    # Organize output files
    with NiiFile(src=fix_clean, assert_exists=True, validate_nifti=True) as fxc:
        func_clean: str = fxc.sym_link(dst=func_clean, relative=True)
    
    func_mean: str = fslmaths(img=func_clean).Tmean().run(out=func_mean, log=log)
    func_stdev: str = fslmaths(img=func_clean).Tstd().run(out=func_stdev, log=log)
    func_tsnr: str = fslmaths(img=func_mean).div(func_stdev).run(out=func_tsnr, log=log)
    
    return (func_clean,
            func_mean,
            func_stdev,
            func_tsnr)


@timeops(log)
def brain_extract(img: str,
                  out: str,
                  mask: Optional[str] = None,
                  robust: bool = False,
                  seg: bool = True,
                  frac_int: Optional[float] = None,
                  log: Optional[LogFile] = None
                 ) -> Tuple[str,str]:
    """Performs brain extraction.
    """
    if mask:
        mask_bool: bool = True
    else:
        mask_bool: bool = False

    brain, _ = bet(img=img,
                   out=out,
                   mask=mask_bool,
                   robust=robust,
                   seg=seg,
                   frac_int=frac_int,
                   log=log)
    mask: str = fslmaths(img=brain).bin().run(out=mask, log=log)
    return brain, mask
