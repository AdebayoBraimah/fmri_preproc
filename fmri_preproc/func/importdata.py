# -*- coding: utf-8 -*-
"""Import data for ``fmri_preproc`` resting-state fMRI preprocessing pipeline.
"""
import os
import nibabel as nib

from typing import (
    Optional,
    List
)

from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.mask import create_mask
from fmri_preproc.utils.enums import SegType
from fmri_preproc.utils.util import Command

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

from fmri_preproc.utils.fslpy import (
    bet,
    fslmaths
)

from fmri_preproc.func.mcdc import motion_outlier

def import_info(outdir: str,
                subid: str,
                scan_pma: float,
                sesid: Optional[str],
                birth_ga: Optional[float] = None,
                log: Optional[LogFile] = None,
                **kwargs
               ) -> None:
    """doc
    """
    pass

def import_func(outdir: str,
                func: str,
                func_echospacing: float,
                func_pedir: str,
                func_brainmask: Optional[str] = None,
                func_slorder: Optional[str] = None,
                sbref: Optional[str] = None,
                sbref_brainmask: Optional[str] = None,
                sbref_echospacing: Optional[float] = None,
                sbref_pedir: Optional[str] = None,
                mask_func: bool = False,
                log: Optional[LogFile] = None
               ) -> None:
    """doc
    """
    pass

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
                 ) -> None:
    """doc
    """
    pass

def import_spinecho():
    """doc
    """
    pass
