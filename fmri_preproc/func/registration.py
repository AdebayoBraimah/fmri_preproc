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
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.enums import BBRType

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

from fmri_preproc.utils.fslpy import (
    flirt,
    fslmaths
)

def fmap_to_struct(outdir: str,
                   fmap: str,
                   fmap_magnitude: str,
                   fmap_brainmask: str,
                   ref_space: str,
                   struct: str,
                   struct_brainmask: str,
                   struct_boundarymask: Optional[str] = None,
                   bbr: bool = False,
                   bbr_slope: float = 0.5,
                   bbr_type: Optional[str] = "signed",
                   src_space: str = 'native',
                   log: Optional[LogFile] = None
                  ) -> None:
    """Perform (rigid-body) fieldmap to structural image registration, with optional BBR.
    """
    # TODO: 
    #   Need to set 'src_space' and 'ref_space' variables

    with WorkDir(src=outdir) as od:
        regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
        with WorkDir(src=regdir) as rd:
            rd.mkdir()
            outdir: str = od.abspath()
            regdir: str = rd.abspath()
            regname: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')

    # Check and validate input fmap associated NIFTI files
    with NiiFile(src=fmap, assert_exists=True, validate_nifti=True) as fmp:
        with NiiFile(src=fmap_magnitude, assert_exists=True, validate_nifti=True) as fmm:
            with NiiFile(src=fmap_brainmask, assert_exists=True, validate_nifti=True) as fmb:
                fmap: str = fmp.abspath()
                fmap_magnitude: str = fmm.abspath()
                fmap_brainmask: str = fmb.abspath()

    # Check and validate input structural associated NIFTI images
    with NiiFile(src=struct, assert_exists=True, validate_nifti=True) as stc:
        with NiiFile(src=struct_brainmask, assert_exists=True, validate_nifti=True) as sbm:
            struct: str = stc.abspath()
            struct_brainmask: str = sbm.abspath()

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
                                "resampled_image": regname + "_img.nii.gz",
                                "affine": regname + "_affine.mat",
                                "inv_affine": regname + "_invaffine.mat",
                                "init_affine": regname + "_init_affine.mat",
                                "resampled_image_init": regname + "_init_img.nii.gz"
                             }
    with TmpDir(src=regdir) as tmp:
        tmp.mkdir()

        # Brain extraction
        struct_brain: str = os.path.join(tmp.src, 'struct_brain.nii.gz')
        struct_brain: str = fslmaths(img=struct).mas(img=struct_brainmask).run(out=struct_brain, log=log)

        fmap_brain: str = os.path.join(tmp.src, 'fmap_brain.nii.gz')
        fmap_brain: str = fslmaths(img=fmap_magnitude).mas(img=fmap_brainmask).run(out=fmap_brain, log=log)

        # Perform registration
        kwargs: Dict[str,str] = { 
                                    "src": fmap_brain,
                                    "ref": struct_brain,
                                    "affine": outputs.get('affine'),
                                    "inv_affine": outputs.get('inv_affine'),
                                    "src2ref": outputs.get('resampled_image')
                                }
        if bbr:
            kwargs: Dict[str,str] = {
                                        **kwargs,
                                        'ref_boundarymask': struct_boundarymask,
                                        'bbrslope': bbr_slope,
                                        'bbrtype': BBRType(bbr_type).name,
                                        'init_affine': outputs.get('init_affine')
                                    }

        # TODO: Get outputs of epireg
        out: str = epireg(**kwargs)

        tmp.rmdir()
    return out

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

def epireg(src: str,
           ref: str,
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
          ) -> None:
    """Linear (FLIRT-based) registration of EPI to ref, with optional BBR and distortion correction.
    """
    with NiiFile(src=src, assert_exists=True, validate_nifti=True) as sr:
        with NiiFile(src=ref, assert_exists=True, validate_nifti=True) as rf:
            src: str = sr.abspath()
            ref: str = rf.abspath()
    
    if log: log.log(f"Performing epireg: Register {src} to {ref} ")

    if basename:
        pass
    else: 
        with NiiFile(src=src) as sr:
            with NiiFile(src=ref) as rf:
                dir1, base1, _ = sr.file_parts()
                _, base2, _ = rf.file_parts()
                basename: str = os.path.join(dir1,f'{base1}_to_{base2}')
                dc_warp: str = sr.rm_ext() + '_dc_warp.nii.gz'
    
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
    _has_fmap = (fmap is not None) and (fmap_brainmask is not None) and (fmap2ref_xfm is not None)
    _do_bbr = (ref_boundarymask is not None)
    _do_dc = _has_fmap and (src_pedir is not None) and (src_echospacing is not None)

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
    pass

def nonlinear_reg(src: str,
                  ref: str,
                  dim: int = 3,
                  quick: bool = False,
                  ref_brainmask: Optional[str] = None,
                  basename: Optional[str] = None,
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
    
    # TODO: Pick-up from here

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
