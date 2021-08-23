# -*- coding: utf-8 -*-
"""Image registration for the fmri_preproce pipeline.

NOTE: 
    Non-FSL dependencies:
        * ``ANTs``: http://stnava.github.io/ANTs/ (**MUST** be compiled natively on target OS)
        * Convert 3D ``c3d``: https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md 
            * Can download pre-compiled binaries: https://sourceforge.net/projects/c3d/files/c3d/
"""
import os

from typing import (
    Dict,
    List,
    Optional, 
    Tuple,
    Union
)

from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.util import Command
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.enums import BBRType

from fmri_preproc.utils.acqparam import (
    _get_axis as get_axis, 
    phase_encode_dict
)

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

from fmri_preproc.utils.fslpy import (
    FSLDIR,
    applywarp,
    applyxfm,
    convertwarp,
    flirt,
    fslmaths,
    fugue,
    invxfm
)

# TODO: ALL registration functions will need:
#   * src_space
#   * ref_space
#   * outdir
# 
# These variables will be used as inputs for all 
# registration functions.

def fmap_to_struct(outdir: str,
                   fmap: str,
                   fmap_magnitude: str,
                   fmap_brainmask: str,
                   struct: str,
                   struct_brainmask: str,
                   struct_boundarymask: Optional[str] = None,
                   src_space: Optional[str] = None,
                   ref_space: Optional[str] = None,
                   bbr: bool = False,
                   bbr_slope: float = 0.5,
                   bbr_type: Optional[str] = "signed",
                   log: Optional[LogFile] = None
                  ) -> None:
    """Perform (rigid-body) fieldmap to structural image registration, with optional BBR.
    """
    # Check and validate input fmap associated NIFTI files
    with NiiFile(src=fmap, assert_exists=True, validate_nifti=True) as fmp:
        with NiiFile(src=fmap_magnitude, assert_exists=True, validate_nifti=True) as fmm:
            with NiiFile(src=fmap_brainmask, assert_exists=True, validate_nifti=True) as fmb:
                fmap: str = fmp.abspath()
                fmap_magnitude: str = fmm.abspath()
                fmap_brainmask: str = fmb.abspath()
                if not src_space:
                    _, src_space, _ = fmp.file_parts()
    
    # Check and validate input structural associated NIFTI images
    with NiiFile(src=struct, assert_exists=True, validate_nifti=True) as stc:
        with NiiFile(src=struct_brainmask, assert_exists=True, validate_nifti=True) as sbm:
            struct: str = stc.abspath()
            struct_brainmask: str = sbm.abspath()
            if not ref_space:
                _, ref_space, _ = stc.file_parts()

    with WorkDir(src=outdir) as od:
        regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
        with WorkDir(src=regdir) as rd:
            rd.mkdir()
            outdir: str = od.abspath()
            regdir: str = rd.abspath()
            regname: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')

    if struct_boundarymask:
        with NiiFile(src=struct_boundarymask, assert_exists=True, validate_nifti=True) as sbm:
            struct_boundarymask: str = sbm.abspath()
    
    if bbr:
        if (struct_boundarymask is None) or (struct_boundarymask == ""):
            if log: log.error("RuntimeError: 'struct_boundarymask' is required for BBR")
            raise RuntimeError("'struct_boundarymask' is required for BBR")

        if isinstance(bbr_slope, float) or isinstance(bbr_slope, int):
            pass
        else:
            if log: log.error(f"TypeError: Input 'bbr_slope': {bbr_slope} is not an 'int' nor 'float'.")
            raise TypeError(f"Input 'bbr_slope': {bbr_slope} is not an 'int' nor 'float'.")
    
    # Define output files
    outputs: Dict[str,str] = {
                                "resampled_image": f"{regname}_img.nii.gz",
                                "affine": f"{regname}_affine.mat",
                                "inv_affine": f"{regname}_invaffine.mat",
                                "init_affine": f"{regname}_init_affine.mat",
                                "resampled_image_init": f"{regname}_init_img.nii.gz"
                             }

    affine: str = outputs.get('affine')
    inv_affine: str = outputs.get('inv_affine')
    src2ref: str = outputs.get('resampled_image')
    init_affine: str = None

    with TmpDir(src=regdir) as tmp:
        tmp.mkdir()

        # Brain extraction
        struct_brain: str = os.path.join(tmp.src, 'struct_brain.nii.gz')
        struct_brain: str = fslmaths(img=struct).mas(img=struct_brainmask).run(out=struct_brain, log=log)

        fmap_brain: str = os.path.join(tmp.src, 'fmap_brain.nii.gz')
        fmap_brain: str = fslmaths(img=fmap_magnitude).mas(img=fmap_brainmask).run(out=fmap_brain, log=log)

        # Perform registration
        kwargs: Dict[str,str] = { 
                                    "outdir": regdir,
                                    "log": log,
                                    "src": fmap_brain,
                                    "ref": struct_brain,
                                    "affine": affine,
                                    "inv_affine": inv_affine,
                                    "src2ref": src2ref,
                                }
        if bbr:
            init_affine: str = outputs.get('init_affine')
            kwargs: Dict[str,str] = {
                                        **kwargs,
                                        "ref_boundarymask": struct_boundarymask,
                                        "bbrslope": bbr_slope,
                                        "bbrtype": BBRType(bbr_type.lower()).name,
                                        "init_affine": init_affine
                                    }

        (affine, 
         inv_affine, 
         src2ref, 
         _) = epireg(**kwargs)

        tmp.rmdir()
    
    return (affine,
            inv_affine,
            src2ref,
            init_affine)

def func_to_sbref():
    """doc-string"""
    pass

def sbref_to_struct():
    """doc-string"""
    pass

def func_to_struct_composite():
    """doc-string"""
    pass

def fmap_to_func_composite():
    """doc-string"""
    pass

def template_to_struct():
    pass

def struct_to_template_composite():
    pass

def func_to_template_composite():
    pass

def epireg(outdir: str,
           src: str,
           ref: str,
           src_space: Optional[str] = None,
           ref_space: Optional[str] = None,
           src_brainmask: Optional[str] = None,
           ref_brainmask: Optional[str] = None,
           src_pedir: Optional[str] = None,
           src_echospacing: Optional[float] = 0.1,
           ref_boundarymask: Optional[str] = None,
           fmap: Optional[str] = None,
           fmap_brainmask: Optional[str] = None,
           fmap2ref_xfm: Optional[str] = None,
           init_xfm: Optional[str] = None,
           bbrslope: float = -0.5,
           bbrtype: str = "signed",
           basename: Optional[str] = None,
           log: Optional[LogFile] = None,
           **kwargs
          ) -> Tuple[str,str,str,str,str]:
    """Linear (FLIRT-based) registration of EPI to ref, with optional BBR and distortion correction.
    """
    with NiiFile(src=src, assert_exists=True, validate_nifti=True) as sr:
        with NiiFile(src=ref, assert_exists=True, validate_nifti=True) as rf:
            src: str = sr.abspath()
            ref: str = rf.abspath()
            if not src_space: _, src_space, _ = sr.file_parts()
            if not ref_space: _, ref_space, _ = rf.file_parts()
    
    if log: log.log(f"Performing epireg: Register {src} to {ref} ")

    if basename:
        pass
    else: 
        with WorkDir(src=outdir) as od:
            regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
            with WorkDir(src=regdir) as rd:
                if not rd.exists(): 
                    rd.mkdir()
                outdir: str = od.abspath()
                regdir: str = rd.abspath()
                basename: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')
                dc_warp: str = os.path.join(regdir, f'{src_space}_dc_warp.nii.gz')
    
    # Define outputs
    outputs: Dict[str,str] = {
                                "affine": f'{basename}_affine.mat',
                                "init_affine": f'{basename}_init_affine.mat',
                                "src2ref_init": f'{basename}_init_img.nii.gz',
                                "warp": f'{basename}_warp.nii.gz',
                                "inv_affine": f'{basename}_invaffine.mat',
                                "inv_warp": f'{basename}_invwarp.nii.gz',
                                "src2ref": f'{basename}_img.nii.gz',
                                "dc_warp": dc_warp
                             }

    # Update outputs/templates from kwargs
    for k,v in kwargs.items():
        if k not in outputs:
            raise ValueError(f'Unknown output template {k}')
        outputs[k] = v
    
    # Logic tests
    _has_fmap: bool = (fmap is not None) and (fmap_brainmask is not None) and (fmap2ref_xfm is not None)
    _do_bbr: bool = (ref_boundarymask is not None)
    _do_dc: bool = _has_fmap and (src_pedir is not None) and (src_echospacing is not None)

    if log: log.debug(f'Has fieldmap: {_has_fmap}')
    if log: log.debug(f'Perform BBR: {_do_bbr}')
    if log: log.debug(f'Perform distortion correction: {_do_dc}')

    # Check fieldmap inputs
    if _has_fmap:
        with NiiFile(src=fmap, assert_exists=True, validate_nifti=True) as fmp:
            with NiiFile(src=fmap_brainmask, assert_exists=True, validate_nifti=True) as fmb:
                fmap: str = fmp.abspath()
                fmap_brainmask: str = fmb.abspath()

    # Phase encode direction
    # TODO: 
    #   * Add phase encoding directions dict
    #   * Add acqp params module
    if src_pedir in ['PA', 'AP', 'IS', 'SI', 'RL', 'LR']:
        src_pedir: str = get_axis(src, src_pedir)[0]
    
    pedict: Dict[str,Tuple[int,str]] = phase_encode_dict()
    if src_pedir:
        src_pedir: Tuple[int,str] = pedict[src_pedir]
    else:
        src_pedir: Tuple[None] = (None, None)
    
    # TODO: Incorporate sigloss
    sigloss2ref: str = None

    with TmpDir(src=regdir) as tmp:
        tmp.mkdir()

        affine: str = None
        inv_affine: str = None
        src2ref: str = None
        warp: str = None
        inv_warp: str = None

        if src_brainmask:
            src_brain: str = os.path.join(tmp.src,'src_brain.nii.gz')
            src_brain: str = fslmaths(img=src).mas(src_brainmask).run(out=src_brain, log=log)
        
        if ref_brainmask:
            ref_brain: str = os.path.join(tmp.src, 'ref_brain.nii.gz')
            ref_brain: str = fslmaths(img=ref).mas(ref_brain).run(out=ref_brain, log=log)
        
        # Initial registration src -> ref (for BBR)
        init_affine: str = outputs.get('init_affine')
        init_img: str = outputs.get('src2ref_init')

        if _do_bbr:
            bbrtype: str = BBRType(bbrtype.lower()).name
            bbrslope: float = float(bbrslope)

            _: Tuple[str] = flirt(
                                src=src if src_brainmask is None else src_brain,
                                ref=ref if ref_brainmask is None else ref_brain,
                                omat=init_affine,
                                out=init_img,
                                dof=6,
                                interp='spline',
                                usesqform=True,
                                nosearch=True,
                                log=log)

        # Main registration: src -> ref
        affine: str = outputs.get('affine')

        flirt_kwargs: Dict[str,str] = {
            "src": src if src_brainmask is None else src_brain,
            "ref": ref if ref_brainmask is None else ref_brain,
            "dof": 6,
            "cost": 'corratio',
            "omat": affine,
            "refweight": sigloss2ref,
            "usesqform": True,
            "nosearch": True,
            "interp": 'spline',
            "log": log,
        }

        if init_xfm:
            with File(src=init_xfm, assert_exists=True) as xfm:
                flirt_kwargs['init'] = xfm.abspath()
        
        if _do_bbr:
            schedule: str = os.path.join(FSLDIR, 'etc', 'flirtsch', 'bbr.sch')

            flirt_kwargs: Dict[str,str] = {
                **flirt_kwargs,
                "cost": 'bbr',
                "wmseg": ref_boundarymask,
                "bbrslope": bbrslope,
                "bbrtype": bbrtype,
                "schedule": schedule,
                "init": init_affine
            }
        
        # Resample fmap and fmap_brainmask -> ref space
        if _do_dc:
            fmap_to_ref: str = os.path.join(tmp.src, 'fmap_to_ref.nii.gz')
            fmap_brainmask_to_ref: str = os.path.join(tmp.src, 'brainmask_to_ref.nii.gz')

            # fmap_to_ref
            fmap_to_ref: str = applyxfm(src=fmap,
                                        ref=ref,
                                        mat=fmap2ref_xfm,
                                        out=fmap_to_ref,
                                        interp='spline',
                                        log=log)
            # fmap_brainmask_to_ref
            fmap_brainmask_to_ref: str = applyxfm(src=fmap_brainmask,
                                                  ref=ref,
                                                  mat=fmap2ref_xfm,
                                                  out=fmap_brainmask_to_ref,
                                                  interp='nearestneighbour',
                                                  log=log)
            
            flirt_kwargs: Dict[str,str] = {
                **flirt_kwargs,
                "echospacing": src_echospacing,
                "pedir": src_pedir[0],
                "fieldmap": fmap_to_ref,
                "fieldmapmask": fmap_brainmask_to_ref,
            }
        
        # Run FLIRT: Perform linear registration
        _: str = flirt(**flirt_kwargs)

        # Invert matrix
        inv_affine: str = outputs.get('inv_affine')
        inv_affine: str = invxfm(inmat=affine, outmat=inv_affine, log=log)

        # If not performing distortion correction return current transforms
        if not _do_dc:
            src2ref: str = outputs.get('src2ref')

            src2ref: str = applyxfm(src=src,
                                    ref=ref,
                                    mat=affine,
                                    out=src2ref,
                                    interp='spline',
                                    log=log)
            tmp.rmdir()
            return affine, inv_affine, src2ref, warp, inv_warp
        
        # Perform distortion correction
        warp: str = outputs.get('warp')

        ## Transform fieldmap_brainmask from ref to src space
        fmap2src_mask: str = os.path.join(tmp.src, "fmap_to_src_brainmask.nii.gz")
        fmap2src_mask: str = applyxfm(src=fmap_brainmask_to_ref,
                                      ref=src,
                                      mat=inv_affine,
                                      out=fmap2src_mask,
                                      interp='nearestneighbour',
                                      log=log)

        ## Transform fieldmap from ref to src space
        fmap2src: str = os.path.join(tmp.src, "fmap_to_src_phase.nii.gz")
        fmap2src: str = applyxfm(src=fmap_to_ref,
                                 ref=src,
                                 mat=inv_affine,
                                 out=fmap2src,
                                 interp='spline',
                                 log=log)
        fmap2src: str = fslmaths(img=fmap2src).mul(fmap2src_mask).run(out=fmap2src, log=log) # Re-mask

        ## Calculate shiftmap and combine transform (src2ref)
        fmap2src_shift: str = os.path.join(tmp.src, "fmap_to_src_shift.nii.gz")

        _: Tuple[str] = fugue(log=log, 
                              loadfmap=fmap2src,
                              mask=fmap2src_mask,
                              saveshift=fmap2src_shift,
                              unmaskshift=True,
                              dwell=src_echospacing,
                              unwarpdir=src_pedir[1])

        _: str = convertwarp(out=warp,
                             ref=ref,
                             postmat=affine,
                             shiftmap=fmap2src_shift,
                             shiftdir=src_pedir[1],
                             abs=True,
                             absout=True,
                             log=log)
        ## Inverse warp
        inv_warp: str = outputs.get('inv_warp')
        inv_warp: str = inv_warp(inwarp=warp,
                                 ref=src,
                                 outwarp=inv_warp,
                                 log=log)

        src2ref: str = outputs.get('src2ref')
        src2ref: str = applywarp(src=src,
                                 ref=ref,
                                 warp=warp,
                                 out=src2ref,
                                 interp='spline',
                                 log=log)
        tmp.rmdir()
    return affine, inv_affine, src2ref, warp, inv_warp

def nonlinear_reg(outdir: str,
                  src: str,
                  ref: str,
                  dim: int = 3,
                  quick: bool = False,
                  ref_brainmask: Optional[str] = None,
                  src_space: Optional[str] = None,
                  ref_space: Optional[str] = None,
                  basename: Optional[str] = None,
                  log: Optional[LogFile] = None,
                  **kwargs: Dict[str]
                 ):
    """Perform ``ANTs`` deformable registration."""
    if isinstance(src, list):
        pass
    elif isinstance(src, str):
        src: List[str] = [ src ]
    else:
        raise TypeError(f"Input for 'src' is not a list or string: {src}")
    
    if isinstance(ref, list):
        pass
    elif isinstance(ref, str):
        ref: List[str] = [ ref ]
    else:
        raise TypeError(f"Input for 'ref' is not a list or string: {ref}")
    
    if basename:
        pass
    else: 
        with NiiFile(src=src[0]) as sr:
            with NiiFile(src=ref[0]) as rf:
                with WorkDir(src=outdir) as od:

                    if not src_space: _, src_space, _ = sr.file_parts()
                    if not ref_space: _, ref_space, _ = rf.file_parts()
                    regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')

                    with WorkDir(src=regdir) as rd:
                        if not rd.exists(): 
                            rd.mkdir()
                        
                        outdir: str = od.abspath()
                        regdir: str = rd.abspath()
                        basename: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')
    
    # Define outputs
    outputs: Dict[str,str] = {
                                "antsout": f'{basename}.ants/output_',
                                "affine": f'{basename}_affine.mat',
                                "warp": f'{basename}_warp.nii.gz',
                                "inv_warp": f'{basename}_invaffine.mat',
                                "src2ref": f'{basename}_img.nii.gz'
                             }

    # Update outputs/templates from kwargs
    for k,v in kwargs.items():
        if k not in outputs:
            raise ValueError(f'Unknown output template {k}')
        outputs[k] = v
    
    # Pick-up from here

def antsRegistrationSyN(fixed: Union[List[str],str],
                        moving: str,
                        dim: int,
                        out_prefix: str,
                        fixed_mask: Optional[str] = None,
                        quick: bool = False,
                        log: Optional[LogFile] = None
                       ) -> Tuple[str]:
    """Python wrapper for ``ANTs`` ``antsRegistrationSyN``."""
    if isinstance(fixed, list):
        pass
    elif isinstance(fixed, str):
        fixed: List[str] = [ fixed ]
    else:
        raise TypeError(f"Input for 'fixed' is not a list or string: {fixed}")
    
    if isinstance(moving, list):
        pass
    elif isinstance(moving, str):
        moving: List[str] = [ moving ]
    else:
        raise TypeError(f"Input for 'moving' is not a list or string: {moving}")
    
    assert len(fixed) == len(moving), 'The same number of fixed and moving images must be provided'

    if quick:
        cmd: Command = Command("antsRegistrationSyNQuick.sh")
    else:
        cmd: Command = Command("antsRegistrationSyN.sh")
    
    cmd.opt("-d"); cmd.opt(f"{dim}")

    for fixed_img, moving_img in zip(fixed, moving):
        cmd.opt("-f"); cmd.opt(f"{fixed_img}")
        cmd.opt("-m"); cmd.opt(f"{moving_img}")
    
    cmd.opt("-o"), cmd.opt(f"{out_prefix}")
    cmd.opt("-t"); cmd.opt("s")
    cmd.opt("-j"); cmd.opt("1")

    if fixed_mask:
        cmd.opt("-x"); cmd.opt(f"{fixed_mask}")
    
    return None
