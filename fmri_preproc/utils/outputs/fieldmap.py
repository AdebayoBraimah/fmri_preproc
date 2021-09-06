# -*- coding: utf-8 -*-
"""Output data filepaths dictionaries for the ``fieldmap`` module.
"""
import os

from typing import (
    Dict,
    Union
)

from fmri_preproc.utils.outputs.outdict import OutDict
from fmri_preproc.utils.workdir import WorkDir

class FmapDict(OutDict):
    """class doc-string.
    """
    def __init__(self,
                 outdir: Union[WorkDir,str]
                ) -> None:
        """Class constructor method.
        NOTE: ``outdir`` should be subject working directory.
        """
        self.outdir: Union[WorkDir,str] = outdir
    
    def outputs(self) -> Dict[str,str]:
        """doc-string.
        """
        if isinstance(self.outdir,WorkDir):
            outdir: str = self.outdir.src
        else:
            outdir: str = self.outdir
        
        with WorkDir(src=outdir) as od:
            self.output: Dict[str,str] = {
                "fmap": od.join("fieldmap.nii.gz"),
                "fmap_mag": od.join("fieldmap_magnitude.nii.gz"),
                "fmap_mask": od.join("fieldmap_brainmask.nii.gz"),
            }
        return self.output