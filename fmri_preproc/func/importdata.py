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

from fmri_preproc.utils.util import timeops
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.mask import create_mask
from fmri_preproc.utils.enums import SegType
from fmri_preproc.func.fieldmap import _merge_rpe as merge_rpe

from fmri_preproc.utils.outputs.importdata import (
    ImportFunc,
    ImportSpinEcho,
    ImportStruct
)

from fmri_preproc.utils.util import (
    dict2json,
    json2dict,
    update_sidecar,
    timeops
)

from fmri_preproc.utils.enums import (
    SegType,
    LogLevel,
    PhaseEncodeDirection
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


# Globlally define (temporary) log file object
with TmpDir(src=os.getcwd()) as tmpd:
    with TmpDir.TmpFile(tmp_dir=tmpd.src, ext='.log') as tmpf:
        log: LogFile = LogFile(log_file=tmpf.src)
        tmpd.rmdir()


@timeops(log)
def import_info(outdir: str,
                func: str,
                scan_pma: Union[int,float,str],
                birth_ga: Optional[float] = None,
                log: Optional[LogFile] = None,
                verbose: bool = False,
                log_datetime: bool = False,
                log_level: str = 'info',
                **kwargs
               ) -> Tuple[str,str,str,Dict[Any,Any]]:
    """Import subject related information and data.
    """
    with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
        _, fname, _ = fn.file_parts()
        
        sesid: str = None

        for item in fname.split('_'):
            if 'sub' in item: subid: str = item[4:]
            if 'ses' in item: sesid: str = item[4:]
            if 'run' in item: runid: str = item[4:]
        
        if not runid:
            runid: str = '01'
        
        if not subid:
            raise RuntimeError("No subject ID can be inferred from the data. \
                The input data should be in BIDS format, and prefixed with \
                'sub-' in the filename.")
    
    with WorkDir(src=outdir) as od:

        if sesid:
            info_dir: str = od.join(f'sub-{subid}',f'ses-{sesid}',f'run-{runid}')
            info_data: str = f'sub-{subid}_ses-{sesid}'
            sub_info: str = f"\n\nSubject information:\n \nsub: {subid} \nses: {sesid} \nrun: {runid}"
        else:
            info_dir: str = od.join(f'sub-{subid}',f'run-{runid}')
            info_data: str = f'sub-{subid}'
            sub_info: str = f"\n\nSubject information:\n \nsub: {subid} \nrun: {runid}"

        with WorkDir(src=info_dir) as ifd:

            logdir: str = ifd.join('logs')

            with WorkDir(src=logdir) as logdr:
                
                info_dir: str = ifd.abspath()
                logdir: str = logdr.abspath()

                if not log:
                    log_file: str = os.path.join(logdr.src, info_data + ".log")
                    log: LogFile = LogFile(log_file=log_file,
                                           print_to_screen=verbose,
                                           format_log_str=log_datetime,
                                           level=LogLevel(log_level.lower()).name)
                
        log.log(sub_info, use_header=True)
        log.log("Creating subject information files.")
            
    # Define outputs
    outputs: Dict[str,Union[int,str]] = { "subject_info": os.path.join(info_dir, info_data + "_info.json") }
    
    # Define info dictionary
    info_name: str = outputs.get('subject_info')
    info: Dict[Any,Any] = json2dict(jsonfile=info_name) if os.path.exists(info_name) else {}
    info: Dict[Any,Any] = { 
                            **info,
                            "subid": subid,
                            "sesid": sesid,
                            "scan_pma": scan_pma,
                            "birth_ga": birth_ga,
                            **kwargs,
                          }
    info_name: str = dict2json(dict=info, jsonfile=info_name)

    return info_dir, logdir, info_name, info


@timeops(log)
def import_func(outdir: str,
                func: str,
                func_echospacing: float,
                func_pedir: str,
                func_brainmask: Optional[str] = None,
                func_slorder: Optional[str] = None,
                func_inplane_accel: Optional[float] = 1,
                mb_factor: Optional[int] = None,
                sbref: Optional[str] = None,
                sbref_brainmask: Optional[str] = None,
                sbref_echospacing: Optional[float] = None,
                sbref_pedir: Optional[str] = None,
                mask_func: bool = False,
                log: Optional[LogFile] = None
               ) -> Tuple[str,Dict[str,str]]:
    """Import function EP image data into the preprocessing pipeline, in addition to:
        * Slice order file (optional)
        * Single-band reference (optional)
    """
    with WorkDir(src=outdir) as od:
        importdir: str = od.join('import')
        with WorkDir(src=importdir) as impd:
            if log: log.log("Creating import data directory")
            importdir: str = impd.abspath()
    
    # Check required inputs
    func_echospacing: float = float(func_echospacing)
    func_pedir: str = PhaseEncodeDirection(func_pedir.upper()).name
    
    # Define outputs
    out: ImportFunc = ImportFunc(outdir=outdir)
    outputs: Dict[str,str] = out.outputs()

    # Import functional EP image data
    with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
        func: str = fn.abspath()
        _, func, _ = fslreorient2std(img=func, out=outputs.get('func'), log=log)
        _: str = update_sidecar(file=func, 
                                echo_spacing=func_echospacing,
                                mb_factor=mb_factor,
                                phase_encode_dir=func_pedir,
                                inplane_accel=float(func_inplane_accel))

    # Create func_mean
    func_mean: str = fslmaths(img=func).Tmean().run(out=outputs.get('func_mean'), log=log)
    _: str = fslroi(img=func, out=outputs.get('func0'), tmin=0, tsize=1, log=log)

    # Create or import functional brain mask
    if func_brainmask:
        with NiiFile(src=func_brainmask, assert_exists=True, validate_nifti=True) as fnb:
            func_brainmask = fnb.abspath()
            _, func_brainmask, _ = fslreorient2std(img=func_brainmask, out=outputs.get('func_brainmask'), log=log)
    else:
        with NiiFile(src=func) as fn:
            dirname, _, _ = fn.file_parts()

        with TmpDir(src=dirname) as tmp:
            brain: str = os.path.join(tmp.src, 'brain.nii.gz')
            brain, _ = bet(img=func_mean, out=brain, frac_int=0.4, mask=False, log=log)
            func_brainmask: str = fslmaths(img=brain).bin().run(out=outputs.get('func_brainmask'), log=log)
    
    # Create dilated brainmask
    func_dilated_brainmask: str = fslmaths(img=func_brainmask).dilM().dilM().bin().run(out=outputs.get('func_dil_brainmask'), log=log)

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

        if sbref_pedir:
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
            _, sbref, _ = fslreorient2std(img=sb.abspath(), out=outputs.get('sbref'), log=log)
            _: str = update_sidecar(file=sbref,
                                    echo_spacing=sbref_echospacing,
                                    phase_encode_dir=sbref_pedir)
        
        if sbref_brainmask:
            with NiiFile(src=sbref_brainmask, assert_exists=True, validate_nifti=True) as sbm:
                _, sbref_brainmask, _ = fslreorient2std(img=sbm.abspath(), out=outputs.get('sbref_brainmask'), log=log)
        else:
            with NiiFile(src=func) as fn:
                dirname, _, _ = fn.file_parts()
            
            with TmpDir(src=dirname) as tmp:
                brain: str = os.path.join(tmp.src, 'brain.nii.gz')
                brain, _ = bet(img=sbref, out=brain, frac_int=0.4, mask=False, log=log)
                sbref_brainmask: str = fslmaths(img=brain).bin().run(out=outputs.get('sbref_brainmask'), log=log)
    return (func,
            func_mean, 
            func_brainmask, 
            func_slorder, 
            sbref, 
            sbref_brainmask,
            outputs)


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
                 ) -> Tuple[str,str,str,Dict[str,str]]:
    """Import anatomical/structural neuroimages into the preprecossing pipeline.
    """
    with WorkDir(src=outdir) as od:
        importdir: str = od.join('import')
        with WorkDir(src=importdir) as impd:
            if log: log.log("Creating import data directory")
            importdir: str = impd.abspath()
    
    # Check required inputs
    dseg_type: str = SegType(dseg_type.lower()).name
    probseg_type: str = SegType(probseg_type.lower()).name

    # Define outputs
    out: ImportStruct = ImportStruct(outdir=outdir)
    outputs: Dict[str,str] = out.outputs()

    # Import anatomical data

    with NiiFile(src=T2w, assert_exists=True, validate_nifti=True) as t2:
        _, T2w, _ = fslreorient2std(img=t2.abspath(), out=outputs.get('T2w'), log=log)
    
    with NiiFile(src=brainmask, assert_exists=True, validate_nifti=True) as bm:
        _, brainmask, _ = fslreorient2std(img=bm.abspath(), out=outputs.get('T2w_brainmask'), log=log)
    
    with NiiFile(src=dseg, assert_exists=True, validate_nifti=True) as ds:
        _, dseg, _ = fslreorient2std(img=ds.abspath(), out=outputs.get('T2w_dseg'), log=log)
        _: str = update_sidecar(file=dseg,
                                dseg_type=dseg_type)
    
    if probseg:
        with NiiFile(src=probseg, assert_exists=True, validate_nifti=True) as pb:
            _, probseg, _ = fslreorient2std(img=pb.abspath(), out=outputs.get('T2w_probseg'), log=log)
            _: str = update_sidecar(file=probseg,
                                    probseg_type=probseg_type)
    
    if wmmask:
        with NiiFile(src=wmmask, assert_exists=True, validate_nifti=True) as wm:
            _, wmmask, _ = fslreorient2std(img=wm.abspath(), out=outputs.get('T2w_wmmask'), log=log)
    else:
        wmmask: str = create_mask(dseg=dseg,
                                  dseg_type=dseg_type,
                                  labels=['wm'],
                                  out=outputs.get('T2w_wmmask'))
    
    if T1w:
        with NiiFile(src=T1w, assert_exists=True, validate_nifti=True) as t1:
            _, T1w, _ = fslreorient2std(img=t1.abspath(), out=outputs.get('T1w'), log=log)   
    
    return (T2w, 
            wmmask, 
            dseg, 
            outputs)


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
                   ) -> str:
    """Import reverse phase-encoded (blip-up, blip-down) spin-echo images.
    """
    if log: log.log("Importing field map.")

    # Define outputs
    out: ImportSpinEcho = ImportSpinEcho(outdir=outdir)
    outputs: Dict[str,str] = out.outputs()

    outdir: str = os.path.join(outdir,"fmap")
    with WorkDir(src=outdir) as od:
        topup_dir: str = od.join("topup")
        with WorkDir(src=topup_dir) as td:
            if log: log.log(f"Making fieldmap directory: {od.src}.")
            outdir: str = od.abspath()
            topup_dir: str = td.abspath()

    # Check inputs to ensure that they are reversed phase encoded
    with TmpDir(src=topup_dir) as tmp:
        
        # input_spinecho: str = os.path.join(topup_dir,'spinecho.nii.gz')
        input_spinecho: str = outputs.get('spinecho')

        if spinecho:
            if log: log.log(f"Input spinecho fieldmap(_epi): {spinecho}.")

            if not spinecho_pedir:
                if log: log.error(f"TypeError: No phase encoding direction list was passed as input.")
                raise TypeError("No phase encoding direction list was passed as input.")
            elif not isinstance(spinecho_pedir, list):
                spinecho_pedir: List[str] = spinecho_pedir.split(",")

            # Verify input phase encoding directions
            spinecho_pedir: List[str] = [ PhaseEncodeDirection(x.upper()).name for x in spinecho_pedir ]

            if log: log.log(f"Phase encodind directions: {spinecho_pedir}")

            with NiiFile(src=spinecho, assert_exists=True, validate_nifti=True) as sp:
                _, input_spinecho, _ = fslreorient2std(img=sp.abspath(), out=input_spinecho)

        elif ap_dir and pa_dir:
            if log:
                log.log(f"Input AP fieldmap(_epi): {ap_dir}.")
                log.log(f"Input PA fieldmap(_epi): {pa_dir}.")
            
            input_ap: str = os.path.join(tmp.src,'AP_epi.nii.gz')
            input_pa: str = os.path.join(tmp.src,'PA_epi.nii.gz')

            with NiiFile(src=ap_dir, assert_exists=True, validate_nifti=True) as ap:
                with NiiFile(src=pa_dir, assert_exists=True, validate_nifti=True) as pa:
                    _, input_ap, _ = fslreorient2std(img=ap.abspath(), out=input_ap)
                    _, input_pa, _ = fslreorient2std(img=pa.abspath(), out=input_pa)
            
            spinecho_pedir: List[str] = [ "PA", "AP" ]
            input_spinecho: str = merge_rpe(out=input_spinecho, log=log, ap_dir=input_ap, pa_dir=input_pa)

        elif lr_dir and rl_dir:
            if log:
                log.log(f"Input LR fieldmap(_epi): {lr_dir}.")
                log.log(f"Input RL fieldmap(_epi): {rl_dir}.")
            
            input_lr: str = os.path.join(tmp.src,'LR_epi.nii.gz')
            input_rl: str = os.path.join(tmp.src,'RL_epi.nii.gz')

            with NiiFile(src=lr_dir, assert_exists=True, validate_nifti=True) as lrd:
                with NiiFile(src=rl_dir, assert_exists=True, validate_nifti=True) as rld:
                    _, input_lr, _ = fslreorient2std(img=lrd.abspath(), out=input_lr)
                    _, input_rl, _ = fslreorient2std(img=rld.abspath(), out=input_rl)
            
            spinecho_pedir: List[str] = [ "LR", "RL" ]
            input_spinecho: str = merge_rpe(out=input_spinecho, log=log, lr_dir=input_lr, rl_dir=input_rl)

        elif is_dir and si_dir:
            if log:
                log.log(f"Input IS fieldmap(_epi): {is_dir}.")
                log.log(f"Input SI fieldmap(_epi): {si_dir}.")
            
            input_is: str = os.path.join(tmp.src,'IS_epi.nii.gz')
            input_si: str = os.path.join(tmp.src,'SI_epi.nii.gz')

            with NiiFile(src=is_dir, assert_exists=True, validate_nifti=True) as isd:
                with NiiFile(src=si_dir, assert_exists=True, validate_nifti=True) as sid:
                    _, input_is, _ = fslreorient2std(img=isd.abspath(), out=input_is)
                    _, input_si, _ = fslreorient2std(img=sid.abspath(), out=input_si)
            
            spinecho_pedir: List[str] = [ "IS", "SI" ]
            input_spinecho: str = merge_rpe(out=input_spinecho, log=log, is_dir=input_is, si_dir=input_si)

        else:
            if log: log.error(f"AttributeError: Input images are not reversed phase encoded fieldmap EPIs: Spinecho: {spinecho} \nAP: {ap_dir} \nPA: {pa_dir} \nLR: {lr_dir} \nRL: {rl_dir} \nIS: {is_dir} \nSI: {si_dir}")
            raise AttributeError(f"Input images are not reversed phase encoded fieldmap EPIs: Spinecho: {spinecho} \nAP: {ap_dir} \nPA: {pa_dir} \nLR: {lr_dir} \nRL: {rl_dir} \nIS: {is_dir} \nSI: {si_dir}")
    
    # Update JSON sidecar
    _: str = update_sidecar(file=input_spinecho,
                            echo_spacing=spinecho_echospacing,
                            phase_encode_dir=spinecho_pedir,
                            epi_factor=spinecho_epifactor,
                            inplane_accel=spinecho_inplaneacc)
    
    return input_spinecho, spinecho_pedir
