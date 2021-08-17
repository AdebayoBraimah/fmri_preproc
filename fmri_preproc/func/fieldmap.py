# -*- coding: utf-8 -*-
"""Fieldmap preprocessing.
"""
import os
import nibabel as nib

from math import pi as PI

from typing import (
    Optional,
    Tuple
)

from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir

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
             spinecho: Optional[str] = None,
             ap_dir: Optional[str] = None,
             pa_dir: Optional[str] = None,
             echo_spacing: Optional[float] = 0.05,
             config: Optional[str] = None,
             verbose: bool = False,
             log: Optional[LogFile] = None
            ) -> Tuple[str, str, str]:
    """Prepare fieldmaps.
    """
    if log:
        log.log("Preparing fieldmaps")

    with WorkDir(src=outdir) as od:
        topup_dir: str = os.path.join(od.src, "topup")
        with WorkDir(src=topup_dir) as td:
            if log:
                log.log(f"Making fieldmap directory: {od.src}.")
            td.mkdir()
            outdir: str = od.abspath()
            topup_dir: str = td.abspath()

    if spinecho:
        if log:
            log.log(f"Input spinecho fieldmap/sbref: {spinecho}.")

        spinecho: NiiFile = NiiFile(src=spinecho, assert_exists=True, validate_nifti=True)
        _, spinecho, _ = fslreorient2std(img=spinecho.src, out=os.path.join(outdir,f"spinecho{spinecho.ext}"))
        spinecho: NiiFile = NiiFile(src=spinecho)

    elif ap_dir and pa_dir:
        if log:
            log.log(f"Input AP fieldmap/sbref: {ap_dir}.")
            log.log(f"Input PA fieldmap/sbref: {pa_dir}.")

        ap_dir: NiiFile = NiiFile(src=ap_dir, assert_exists=True, validate_nifti=True)
        pa_dir: NiiFile = NiiFile(src=pa_dir, assert_exists=True, validate_nifti=True)
        _, ap_dir, _ = fslreorient2std(img=ap_dir.src, out=os.path.join(outdir,f"sbref_ap{ap_dir.ext}"))
        _, pa_dir, _ = fslreorient2std(img=pa_dir.src, out=os.path.join(outdir,f"sbref_pa{pa_dir.ext}"))

        spinecho: str = _merge_rpe(ap_dir=ap_dir,
                                   pa_dir=pa_dir,
                                   out=os.path.join(outdir,"spinecho"),
                                   log=log)
        
        spinecho: NiiFile = NiiFile(src=spinecho, assert_exists=True, validate_nifti=True)

    else:
        raise AttributeError("Neither of the 'spinecho' OR 'ap_dir' and 'pa_dir' options were specified.")

    slices: int = nib.load(filename=spinecho.abspath()).header.get('dim','')[3]

    if config:
        pass
    else:
        config: str = _get_b0_conf(slices=slices, default=True)

    dist_acqp: str = write_fieldmap_acq_params(effective_echo_spacing=echo_spacing, 
                                               out_prefix=os.path.join(topup_dir, "fieldmap"))

    # Run topup
    if log:
        log.log("Running topup")

    (_,
     field_img,
     mag_img) = topup(img=spinecho.src, 
                      param=dist_acqp, 
                      out=os.path.join(topup_dir, "topup_dist_corr"),
                      config=config,
                      fout=True,
                      iout=True,
                      subsamp="1,1,1,1,1,1,1,1,1",
                      verbose=verbose,
                      log=log)

    # Compute field map
    if log:
        log.log("Computing fieldmaps.")

    fmap: str = fslmaths(img=field_img).mul(2 * PI).run(out=os.path.join(outdir,"fmap"), log=log)
    mag: str = fslmaths(img=mag_img).Tmean().run(out=os.path.join(outdir,"fmap_magnitude"), log=log)

    # Create brain mask
    if log:
        log.log("Creating fieldmap brain mask.")

    with TmpDir(tmp_dir=os.path.join(topup_dir)) as tmp:
        tmp.mkdir()
        brain: str = os.path.join(tmp.src,"brain")
        fmap_brain, _ = bet(img=fmap,
                            out=brain,
                            mask=False,
                            robust=True,
                            log=log)
        fmap_mask: str = fslmaths(img=fmap_brain).bin().run(out=os.path.join(outdir,"fmap_brainmask"),log=log)
        tmp.rmdir()
    return (fmap,
            mag,
            fmap_mask)

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

def write_fieldmap_acq_params(effective_echo_spacing: Optional[float] = 0.05, 
                              out_prefix: Optional[str] = 'file'
                             ) -> str:
    """Creates acquisition parameters files for fMRI acquisition provided the number of frames and the 
    effective echo spacing.

    NOTE:
        The corresponding 4D single band reference (sbref) files are assumed to be concatenated into a single file that:
            * The first image was acquired in the posterior -> anterior (P -> A, PA) direction.
            * The second image was acquired in the anterior -> posterior (A -> P, AP) direction.
    
    Usage example:
        >>> dist_acqp = write_fieldmap_acq_params(num_frames=1200, 
        ...                                       effective_echo_spacing=0.05,
        ...                                       out_prefix="dist")
        ...
        
    Arguments:
        num_frames: Number of temporal frames/dynamics.
        effective_echo_spacing: Effective echo spacing.
        out_prefix: Output file prefix.
        
    Returns:
        String that corresponds to the acquisition parameter file for distortion correction.
    
    Raises:
        TypeError: Exception that is raised when ``effective_echo_spacing`` is not an ``int`` or a ``float``.
    """
    if (not isinstance(effective_echo_spacing, int) and 
        not isinstance(effective_echo_spacing, float)):
        raise TypeError(f"Input for effective_echo_spacing: {effective_echo_spacing} is not a float.")
    
    # Generate distortion correction
    out_dist: str = out_prefix + "_distortion_correction.acqp"
    
    # Write distortion correction acqp file
    with open(out_dist,'w') as f:
        f.write(f"0 1 0 {effective_echo_spacing}\n")
        f.write(f"0 -1 0 {effective_echo_spacing}\n")
        f.close()
    
    # Obtain absolute file paths
    out_dist: str = os.path.abspath(out_dist)
    return out_dist

def _merge_rpe(ap_dir: str,
               pa_dir: str,
               out: str,
               log: Optional[LogFile]
              ) -> str:
    """Merges reverse phase-encoded (rPE) fieldmaps or single band references (sbrefs).
    """
    ap_dir: NiiFile = NiiFile(src=ap_dir, assert_exists=True, validate_nifti=True)
    pa_dir: NiiFile = NiiFile(src=pa_dir, assert_exists=True, validate_nifti=True)

    fslmerge(out,
             "t",
             None,
             log,
             pa_dir.abspath(),
             ap_dir.abspath())
    return out
