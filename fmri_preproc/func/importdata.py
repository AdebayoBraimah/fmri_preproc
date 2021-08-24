# -*- coding: utf-8 -*-
"""Import data for ``fmri_preproc`` resting-state fMRI preprocessing pipeline.
"""
import os
import nibabel as nib

from typing import (
    Any,
    Dict,
    Optional,
    List,
    Tuple,
    Union
)

from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.mask import create_mask
from fmri_preproc.utils.enums import SegType
from fmri_preproc.func.fieldmap import _merge_rpe as merge_rpe

from fmri_preproc.utils.util import (
    dict2json,
    json2dict,
    update_sidecar,
    timeops
)

from fmri_preproc.utils.enums import (
    SegType,
    PhaseEncodeDirection,
    PhaseUnits
)

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

from fmri_preproc.utils.fslpy import (
    bet,
    fslmaths,
    fslreorient2std,
    fslroi
)

from fmri_preproc.func.mcdc import (
    motion_outlier,
    write_slice_order
)

# Globlally define log file object
log: LogFile = None

@timeops(log)
def testfunc(s, log):
    log.log(s)
    return s

@timeops(log)
def import_info(outdir: str,
                subid: str,
                scan_pma: float,
                sesid: Optional[str] = None,
                birth_ga: Optional[float] = None,
                log: Optional[LogFile] = None,
                **kwargs
               ) -> Tuple[str, Dict[Any,Any]]:
    """Import subject related information and data.
    """
    with WorkDir(src=outdir) as od:
        if not od.exists(): 
            if log: log.log("Creating subject information files.")
            od.mkdir()
        outdir: str = od.abspath()
    
    # Define outputs
    if sesid:
        info_data: str = f'sub-{subid}_ses-{sesid}'
    else:
        info_data: str = f'sub-{subid}'
    
    outputs: Dict[str,Union[int,str]] = { "subject_info": os.path.join(outdir, info_data + "_info.json") }
    
    # Define info dictionary
    info_name: str = outputs.get('subject_info')
    info: Dict[Any,Any] = json2dict(jsonfile=info_name) if os.path.exists(info_name) else {}
    info: Dict[Any,Any] = { 
                            **info,
                            "subid": subid,
                            "sesid": sesid,
                            "scan_pma": scan_pma,
                            "birth_ga": birth_ga,
                          }
    info_name: str = dict2json(dict=info, jsonfile=info_name)
    return info_name, info

@timeops(log)
def import_func(outdir: str,
                func: str,
                func_echospacing: float,
                func_pedir: str,
                func_brainmask: Optional[str] = None,
                func_slorder: Optional[str] = None,
                mb_factor: Optional[int] = None,
                sbref: Optional[str] = None,
                sbref_brainmask: Optional[str] = None,
                sbref_echospacing: Optional[float] = None,
                sbref_pedir: Optional[str] = None,
                mask_func: bool = False,
                log: Optional[LogFile] = None
               ) -> Tuple[str]:
    """Import function EP image data into the preprocessing pipeline, in addition to:
        * Slice order file (optional)
        * Single-band reference (optional)
    """
    with WorkDir(src=outdir) as od:
        importdir: str = os.path.join(od.src, 'import')
        with WorkDir(src=importdir) as impd:
            if not impd.exists():
                if log: log.log("Creating import data directory")
                impd.mkdir()
            importdir: str = impd.abspath()
    
    # Check required inputs
    func_echospacing: float = float(func_echospacing)
    func_pedir: str = PhaseEncodeDirection(func_pedir.upper()).name
    
    # Define outputs
    outputs: Dict[str,str] = {
                                "func": os.path.join(importdir,'func.nii.gz'),
                                "func_mean": os.path.join(importdir,'func_mean.nii.gz'),
                                "func_brainmask": os.path.join(importdir,'func_brainmask.nii.gz'),
                                "func_slorder": os.path.join(importdir,'func.slorder'),
                                "sbref": os.path.join(importdir,'sbref.nii.gz'),
                                "sbref_brainmask": os.path.join(importdir,'sbref_brainmask.nii.gz'),
                                "func_std": os.path.join(importdir,'func_std.nii.gz'),
                                "func_tsnr": os.path.join(importdir,'func_tsnr.nii.gz'),
                                "func0": os.path.join(importdir,'func0.nii.gz'),
                                "func_dil_brainmask": os.path.join(importdir,'func_dilated_brainmask.nii.gz'),
                                "func_metrics": os.path.join(importdir,'func_regressors.tsv'),
                                "func_out_plot": os.path.join(importdir,'func_outliers.png'),
                             }

    # Import functional EP image data
    with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
        func: str = fn.abspath()
        func: str = fslreorient2std(img=func, out=outputs.get('func'), log=log)
        _: str = update_sidecar(file=func, 
                                echo_spacing=func_echospacing, 
                                phase_encode_dir=func_pedir)

    # Create func_mean
    func_mean: str = fslmaths(img=func).Tmean().run(out=outputs.get('func_mean'), log=log)
    _: str = fslroi(img=func, out=outputs.get('func0'), tmin=0, tsize=1, log=log)

    # Create or import functional brain mask
    if func_brainmask:
        with NiiFile(src=func_brainmask, assert_exists=True, validate_nifti=True) as fnb:
            func_brainmask: str = fnb.abspath()
            func_brainmask: str = fslreorient2std(img=func_brainmask, out=outputs.get('func_brainmask'), log=log)
    else:
        with NiiFile(src=func) as fn:
            dirname, _, _ = fn.file_parts()

        with TmpDir(src=dirname) as tmp:
            tmp.mkdir()
            brain: str = os.path.join(tmp.src, 'brain.nii.gz')
            brain: str = bet(img=func_mean, out=brain, frac_int=0.4, mask=False, log=log)
            func_brainmask: str = fslmaths(img=brain).bin().run(out=outputs.get('func_brainmask'), log=log)
            tmp.rmdir()
    
    # Create dilated brainmask
    func_dilated_brainmask: str = fslmaths(img=func_brainmask).dilM().dilM().bin().run(out=outputs.get('func_dilated_brainmask'), log=log)

    # Mask func (with dilated brainmask) to reduce file size
    if mask_func:
        func : str = fslmaths(img=func).mas(img=func_dilated_brainmask).run(out=func)
    
    # Create std and tsnr
    func_std: str = fslmaths(img=func).Tstd().run(out=outputs.get('func_std'), log=log)
    _: str = fslmaths(img=func_mean).div(func_std).run(out=outputs.get('func_tsnr'), log=log)

    # Compute motion metrics and outliers
    _: Tuple[str,int,None] = motion_outlier(func=func, 
                                            metric_name=outputs.get('func_metrics'),
                                            plot_name=outputs.get('func_out_plot'))

    # Copy functional slice order, OR recompute if provided multi-band factor
    if func_slorder:
        with File(src=func_slorder, assert_exists=True) as fs:
            func_slorder: str = fs.copy(dst=outputs.get('func_slorder'))
    elif mb_factor:
        slices: int = nib.load(filename=func).header.get('dim','')[3]
        func_slorder: str = write_slice_order(slices=slices,
                                              mb_factor=mb_factor,
                                              mode="single-shot",
                                              out_file=outputs.get('func_slorder'),
                                              return_mat=False)

    if sbref:

        if not sbref_pedir:
            sbref_pedir: str = PhaseEncodeDirection(sbref_pedir.upper()).name
        else:
            if log: log.error("RuntimeError: sbref' phase-encoding direction required.")
            raise RuntimeError("'sbref' phase-encoding direction required.")
        
        if sbref_echospacing:
            sbref_echospacing: float = float(sbref_echospacing)
        else:
            if log: log.error("RuntimeError: sbref' echo-spacing required.")
            raise RuntimeError("'sbref' echo-spacing required.")
        
        with NiiFile(src=sbref, assert_exists=True, validate_nifti=True) as sb:
            sbref: str = fslreorient2std(img=sb.abspath(), out=outputs.get('sbref'), log=log)
            _: str = update_sidecar(file=sbref,
                                    echo_spacing=sbref_echospacing,
                                    phase_encode_dir=sbref_pedir)
        
        if sbref_brainmask:
            with NiiFile(src=sbref_brainmask, assert_exists=True, validate_nifti=True) as sbm:
                sbref_brainmask: str = fslreorient2std(img=sbm.abspath(), out=outputs.get('sbref_brainmask'), log=log)
        else:
            with NiiFile(src=func) as fn:
                dirname, _, _ = fn.file_parts()
            
            with TmpDir(src=dirname) as tmp:
                tmp.mkdir()
                brain: str = os.path.join(tmp.src, 'brain.nii.gz')
                brain: str = bet(img=sbref, out=brain, frac_int=0.4, mask=False, log=log)
                sbref_brainmask: str = fslmaths(img=brain).bin().run(out=outputs.get('sbref_brainmask'), log=log)
    return (func,
            func_mean, 
            func_brainmask, 
            func_slorder, 
            sbref, 
            sbref_brainmask)

@timeops(log)
def import_struct(outdir: str,
                  T2w: str,
                  brainmask: str,
                  dseg: str,
                  dseg_type: Optional[str] = 'drawem',
                  probseg: Optional[str] = None,
                  probseg_type: Optional[str] = 'drawem',
                  wmmask: Optional[str] = None,
                  T1w: Optional[str] = None,
                  log: Optional[LogFile] = None
                 ) -> Tuple[str,str,str]:
    """Import anatomical/structural neuroimages into the preprecossing pipeline.
    """
    with WorkDir(src=outdir) as od:
        importdir: str = os.path.join(od.src, 'import')
        with WorkDir(src=importdir) as impd:
            if not impd.exists():
                if log: log.log("Creating import data directory")
                impd.mkdir()
            importdir: str = impd.abspath()
    
    # Check required inputs
    dseg_type: str = SegType(dseg_type.lower()).name
    probseg_type: str = SegType(probseg_type.lower()).name

    # Define outputs
    outputs: Dict[str,str] = {
                                "T1w": os.path.join(importdir,'T1w.nii.gz'),
                                "T2w": os.path.join(importdir,'T2w.nii.gz'),
                                "T2w_brainmask": os.path.join(importdir,'T2w_brainmask.nii.gz'),
                                "T2w_dseg": os.path.join(importdir,'T2w_dseg.nii.gz'),
                                "T2w_wmmask": os.path.join(importdir,'T2w_wmmask.nii.gz'),
                                "T2w_probseg": os.path.join(importdir,'T2w_probseg.nii.gz'),
                             }

    # Import anatomical data

    with NiiFile(src=T2w, assert_exists=True, validate_nifti=True) as t2:
        T2w: str = fslreorient2std(img=t2.abspath(), out=outputs.get('T2w'), log=log)
    
    with NiiFile(src=brainmask, assert_exists=True, validate_nifti=True) as bm:
        brainmask: str = fslreorient2std(img=bm.abspath(), out=outputs.get('T2w_brainmask'), log=log)
    
    with NiiFile(src=dseg, assert_exists=True, validate_nifti=True) as ds:
        dseg: str = fslreorient2std(img=ds.abspath(), out=outputs.get('T2w_dseg'), log=log)
        _: str = update_sidecar(file=dseg,
                                seg_type=dseg_type)
    
    if probseg:
        with NiiFile(src=probseg, assert_exists=True, validate_nifti=True) as pb:
            probseg: str = fslreorient2std(img=pb.abspath(), out=outputs.get('T2w_probseg'), log=log)
            _: str = update_sidecar(file=probseg,
                                    seg_type=probseg_type)
    
    if wmmask:
        with NiiFile(src=wmmask, assert_exists=True, validate_nifti=True) as wm:
            wmmask: str = fslreorient2std(img=wm.abspath(), out=outputs.get('T2w_wmmask'), log=log)
    else:
        wmmask: str = create_mask(dseg=dseg,
                                  dseg_type=dseg_type,
                                  labels=['wm'],
                                  out=outputs.get('T2w_wmmask'))
    
    if T1w:
        with NiiFile(src=T1w, assert_exists=True, validate_nifti=True) as t1:
            T1w: str = fslreorient2std(img=t1.abspath(), out=outputs.get('T1w'), log=log)   
    
    return T2w, wmmask, dseg

@timeops(log)
def import_spinecho(outdir: str,
                    spinecho: Optional[str] = None,
                    spinecho_echospacing: Optional[float] = 0.1,
                    spinecho_pedir: Union[str,List[str]] = None,
                    ap_dir: Optional[str] = None,
                    pa_dir: Optional[str] = None,
                    lr_dir: Optional[str] = None,
                    rl_dir: Optional[str] = None,
                    is_dir: Optional[str] = None,
                    si_dir: Optional[str] = None,
                    spinecho_epifactor: Optional[int] = None,
                    spinecho_inplaneacc: Optional[float] = 1,
                    log: Optional[LogFile] = None
                   ) -> None:
    """Import reverse phase-encoded (blip-up, blip-down) spin-echo images.
    """
    # TODO: Pick-up from here
    pass

