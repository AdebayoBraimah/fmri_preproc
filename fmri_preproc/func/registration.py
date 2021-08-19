# -*- coding: utf-8 -*-
"""Image registration for the fmri_preproce pipeline.

NOTE: 
    Non-FSL dependencies:
        * ``ANTs``: http://stnava.github.io/ANTs/ (**MUST** be compiled natively on target OS)
        * Convert 3D ``c3d``: https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md 
            * Can download pre-compiled binaries: https://sourceforge.net/projects/c3d/files/c3d/
"""
import os
from fmri_preproc.utils.workdir import WorkDir
from typing import (
    Dict,
    List,
    Optional, 
    Tuple,
    Union
)

from fmri_preproc.utils.util import Command
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)


def fmap_to_struct(outdir: str,
                   fmap: str,
                   fmap_magnitude: str,
                   fmap_brainmask: str,
                   struct: str,
                   struct_brainmask: str,
                   struct_boundarymask: Optional[str] = None,
                   bbr: bool = True,
                   bbr_slope: float = 0.5,
                   bbr_type: Optional[str] = None, # This SHOULD be an enum
                   log: Optional[LogFile] = None
                  ) -> None:
    """Perform (rigid-body) fieldmap to structural image registration, with optional BBR.
    """
    with WorkDir(src=outdir) as od:
        regdir: str = os.path.join(od.abspath(),'reg') # This might change as src_to_ref space directories need to be specififed.
        outdir: str = od.abspath()

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
            if log:
                log.error("RuntimeError: struct_boundarymask is required for BBR")
            raise RuntimeError('struct_boundarymask is required for BBR')

    pass

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

def epireg():
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
