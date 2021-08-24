# -*- coding: utf-8 -*-
"""Fieldmap preprocessing.

NOTE:
    * Check this resource pertaining to the odd number of slices: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;28f6c983.1806
"""
import os
import nibabel as nib

from typing import (
    Dict,
    List,
    Optional,
    Tuple,
    Union
)

from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.acqparam import write_func_params
from fmri_preproc.utils.enums import PhaseEncodeDirection

from fmri_preproc.utils.fslpy import (
    bet,
    FSLDIR,
    fslmerge,
    fslreorient2std,
    topup,
    fslmaths
)

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

def fieldmap(outdir: str,
             echo_spacing: Optional[float] = 0.1,
             spinecho: Optional[str] = None,
             pedir: Union[str,List[str]] = None,
             ap_dir: Optional[str] = None,
             pa_dir: Optional[str] = None,
             lr_dir: Optional[str] = None,
             rl_dir: Optional[str] = None,
             is_dir: Optional[str] = None,
             si_dir: Optional[str] = None,
             epifactor: Optional[int] = None,
             inplane_acc: Optional[float] = 1,
             config: Optional[str] = None,
             verbose: bool = False,
             log: Optional[LogFile] = None
            ) -> Tuple[str, str, str]:
    """Prepare fieldmaps.

    NOTE: 
        Input directory should be parent FEAT/processing directory path.
    """
    # TODO: 
    # Re-configure function to ensure compatibility with import_spinecho and import fieldmap functions
    if log: log.log("Preparing fieldmaps")

    outdir: str = os.path.join(outdir,"fmap")
    with WorkDir(src=outdir) as od:
        topup_dir: str = os.path.join(od.src, "topup")
        with WorkDir(src=topup_dir) as td:
            if not td.exists(): 
                if log: log.log(f"Making fieldmap directory: {od.src}.")
                td.mkdir()
            outdir: str = od.abspath()
            topup_dir: str = td.abspath()
    
    # Define output files
    outputs: Dict[str,str] = {
                                "fmap": os.path.join(outdir,"fieldmap.nii.gz"),
                                "fmap_mag": os.path.join(outdir,"fieldmap_magnitude.nii.gz"),
                                "fmap_mask": os.path.join(outdir,"fieldmap_magnitude.nii.gz"),
                                "topup_out": os.path.join(topup_dir, "topup_dist_corr"),
                             }

    with TmpDir(src=topup_dir) as tmp:
        tmp.mkdir()
        
        input_spinecho: str = os.path.join(topup_dir,'spinecho.nii.gz')

        if spinecho:
            if log: log.log(f"Input spinecho fieldmap(_epi): {spinecho}.")

            if not pedir:
                if log: log.error(f"TypeError: No phase encoding direction list was passed as input.")
                raise TypeError("No phase encoding direction list was passed as input.")
            elif not isinstance(pedir, list):
                pedir: List[str] = pedir.split(",")

            # Verify input phase encoding directions
            pedir: List[str] = [ PhaseEncodeDirection(x.upper()).name for x in pedir ]

            if log: log.log(f"Phase encodind directions: {pedir}")

            with NiiFile(src=spinecho, assert_exists=True, validate_nifti=True) as sp:
                input_spinecho: str = fslreorient2std(img=sp.abspath(), out=input_spinecho)
                tmp.rmdir()

        elif ap_dir and pa_dir:
            if log:
                log.log(f"Input AP fieldmap(_epi): {ap_dir}.")
                log.log(f"Input PA fieldmap(_epi): {pa_dir}.")
            
            input_ap: str = os.path.join(tmp.src,'AP_epi.nii.gz')
            input_pa: str = os.path.join(tmp.src,'PA_epi.nii.gz')

            with NiiFile(src=ap_dir, assert_exists=True, validate_nifti=True) as ap:
                with NiiFile(src=pa_dir, assert_exists=True, validate_nifti=True) as pa:
                    input_ap: str = fslreorient2std(img=ap.abspath(), out=input_ap)
                    input_pa: str = fslreorient2std(img=pa.abspath(), out=input_pa)
            
            pedir: List[str] = [ "PA", "AP" ]
            input_spinecho: str = _merge_rpe(out=input_spinecho, log=log, ap_dir=input_ap, pa_dir=input_pa)
            tmp.rmdir()
        elif lr_dir and rl_dir:
            if log:
                log.log(f"Input LR fieldmap(_epi): {lr_dir}.")
                log.log(f"Input RL fieldmap(_epi): {rl_dir}.")
            
            input_lr: str = os.path.join(tmp.src,'LR_epi.nii.gz')
            input_rl: str = os.path.join(tmp.src,'RL_epi.nii.gz')

            with NiiFile(src=lr_dir, assert_exists=True, validate_nifti=True) as lrd:
                with NiiFile(src=rl_dir, assert_exists=True, validate_nifti=True) as rld:
                    input_lr: str = fslreorient2std(img=lrd.abspath(), out=input_lr)
                    input_rl: str = fslreorient2std(img=rld.abspath(), out=input_rl)
            
            pedir: List[str] = [ "LR", "RL" ]
            input_spinecho: str = _merge_rpe(out=input_spinecho, log=log, lr_dir=input_lr, rl_dir=input_rl)
            tmp.rmdir()
        elif is_dir and si_dir:
            if log:
                log.log(f"Input IS fieldmap(_epi): {is_dir}.")
                log.log(f"Input SI fieldmap(_epi): {si_dir}.")
            
            input_is: str = os.path.join(tmp.src,'IS_epi.nii.gz')
            input_si: str = os.path.join(tmp.src,'SI_epi.nii.gz')

            with NiiFile(src=is_dir, assert_exists=True, validate_nifti=True) as isd:
                with NiiFile(src=si_dir, assert_exists=True, validate_nifti=True) as sid:
                    input_is: str = fslreorient2std(img=isd.abspath(), out=input_is)
                    input_si: str = fslreorient2std(img=sid.abspath(), out=input_si)
            
            pedir: List[str] = [ "IS", "SI" ]
            input_spinecho: str = _merge_rpe(out=input_spinecho, log=log, is_dir=input_is, si_dir=input_si)
            tmp.rmdir()
        else:
            if log: log.error(f"AttributeError: Input images are not reversed phase encoded fieldmap EPIs: Spinecho: {spinecho} \nAP: {ap_dir} \nPA: {pa_dir} \nLR: {lr_dir} \nRL: {rl_dir} \nIS: {is_dir} \nSI: {si_dir}")
            raise AttributeError(f"Input images are not reversed phase encoded fieldmap EPIs: Spinecho: {spinecho} \nAP: {ap_dir} \nPA: {pa_dir} \nLR: {lr_dir} \nRL: {rl_dir} \nIS: {is_dir} \nSI: {si_dir}")

    slices: int = nib.load(filename=input_spinecho).header.get('dim','')[3]

    if config:
        pass
    else:
        config: str = _get_b0_conf(slices=slices, default=True)
    
    dist_acqp: str = os.path.join(topup_dir,'spinecho.acqp')

    dist_acqp: str = write_func_params(epi=input_spinecho,
                                       echospacing=echo_spacing,
                                       pedir=pedir,
                                       out=dist_acqp,
                                       epifactor=epifactor,
                                       inplane_acc=inplane_acc)

    # Run topup
    if log: log.log("Performing topup distortion field estimation")

    _,field_img, mag_img = topup(img=spinecho.src, 
                                 param=dist_acqp, 
                                 out=outputs.get('topup_out'),
                                 config=config,
                                 fout=True,
                                 iout=True,
                                 subsamp="1,1,1,1,1,1,1,1,1",
                                 verbose=verbose,
                                 log=log)

    # Compute field map
    if log: log.log("Computing fieldmaps.")

    fmap: str = fslmaths(img=field_img).mul(6.2832).run(out=outputs.get('fmap'), log=log)
    mag: str = fslmaths(img=mag_img).Tmean().run(out=outputs.get('fmap_mag'), log=log)

    # Create brain mask
    if log: log.log("Creating fieldmap brain mask.")

    with TmpDir(src=os.path.join(topup_dir)) as tmp:
        tmp.mkdir()
        brain: str = os.path.join(tmp.src,"brain")
        fmap_brain, _ = bet(img=fmap,
                            out=brain,
                            mask=False,
                            robust=True,
                            log=log)
        fmap_mask: str = fslmaths(img=fmap_brain).bin().run(out=outputs.get('fmap_mask'),log=log)
        tmp.rmdir()
    return fmap, mag, fmap_mask

def _get_b0_conf(slices: int,
                 default: bool = True,
                ) -> str:
    """Returns the optimal b0 configuration file for use with FSL's ``topup`` and ``eddy``.
    """
    config_parent_dir: str = os.path.join(FSLDIR,'etc','flirtsch')

    if ((slices % 2) == 0) and default:
        config: str = os.path.join(config_parent_dir,'b02b0.cnf')
        config: File = File(src=config, assert_exists=True)
        return config.abspath()

    if (slices % 4) == 0:
        config: str = os.path.join(config_parent_dir,'b02b0_4.cnf')
        config: File = File(src=config, assert_exists=True)
        return config.abspath()
    elif (slices % 2) == 0:
        config: str = os.path.join(config_parent_dir,'b02b0_2.cnf')
        config: File = File(src=config, assert_exists=True)
        return config.abspath()
    else:
        config: str = os.path.join(config_parent_dir,'b02b0_1.cnf')
        config: File = File(src=config, assert_exists=True)
        return config.abspath()

def _merge_rpe(out: str,
               ap_dir: Optional[str] = None,
               pa_dir: Optional[str] = None,
               lr_dir: Optional[str] = None,
               rl_dir: Optional[str] = None,
               is_dir: Optional[str] = None,
               si_dir: Optional[str] = None,
               log: Optional[LogFile] = None
              ) -> str:
    """Helper function that merges reverse phase-encoded (rPE) fieldmaps (EPI fieldmaps).
    """
    # NOTE: Input order for `fslmerge`:
    #   * PA, AP
    #   * LR, RL
    #   * IS, SI

    # Verify and validate input files
    if ap_dir and pa_dir:
        with NiiFile(src=ap_dir, assert_exists=True, validate_nifti=True) as ap:
            with NiiFile(src=pa_dir, assert_exists=True, validate_nifti=True) as pa:
                rpedir1: str = pa.abspath()
                rpedir2: str = ap.abspath()
    elif lr_dir and rl_dir:
        with NiiFile(src=lr_dir, assert_exists=True, validate_nifti=True) as lrd:
            with NiiFile(src=rl_dir, assert_exists=True, validate_nifti=True) as rld:
                rpedir1: str = lrd.abspath()
                rpedir2: str = rld.abspath()
    elif is_dir and si_dir:
        with NiiFile(src=is_dir, assert_exists=True, validate_nifti=True) as isd:
            with NiiFile(src=si_dir, assert_exists=True, validate_nifti=True) as sid:
                rpedir1: str = isd.abspath()
                rpedir2: str = sid.abspath()
    else:
        if log: log.error(f"AttributeError: Input images are not reversed phase encoded: AP: {ap_dir} \nPA: {pa_dir} \nLR: {lr_dir} \nRL: {rl_dir} \nIS: {is_dir} \nSI: {si_dir}")
        raise AttributeError(f"Input images are not reversed phase encoded: AP: {ap_dir} \nPA: {pa_dir} \nLR: {lr_dir} \nRL: {rl_dir} \nIS: {is_dir} \nSI: {si_dir}")

    out: str = fslmerge(out,
                        "t",
                        None,
                        log,
                        rpedir1,
                        rpedir2)
    return out
