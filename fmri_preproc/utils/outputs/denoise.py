# -*- coding: utf-8 -*-
"""Output data filepaths dictionaries for the ``denoise`` module.
"""
import os

from typing import (
    Dict,
    Union
)

from fmri_preproc.utils.outputs.outdict import OutDict
from fmri_preproc.utils.workdir import WorkDir


class FIXClassify(OutDict):
    """class doc-string.
    """
    def __init__(self,
                 outdir: Union[WorkDir,str]
                ) -> None:
        """Class constructor method.
        NOTE: ``outdir`` should be subject working directory.
        """
        self.outdir: Union[WorkDir,str] = outdir
        super(FIXClassify, self).__init__()
    
    def outputs(self) -> Dict[str,str]:
        """doc-string.
        """
        if isinstance(self.outdir,WorkDir):
            outdir: str = self.outdir.src
        else:
            outdir: str = self.outdir
        
        denoisedir: str = os.path.join(outdir,'denoise')

        with WorkDir(src=denoisedir) as od:
            self.output: Dict[str,str] = {
                "fix_labels": od.join("fix_labels.txt"),
                "fix_regressors": od.join("fix_regressors.tsv")
            }
        return self.output


class FIXApply(OutDict):
    """class doc-string.
    """
    def __init__(self,
                 outdir: Union[WorkDir,str]
                ) -> None:
        """Class constructor method.
        NOTE: ``outdir`` should be subject working directory.
        """
        self.outdir: Union[WorkDir,str] = outdir
        super(FIXClassify, self).__init__()
    
    def outputs(self) -> Dict[str,str]:
        """doc-string.
        """
        if isinstance(self.outdir,WorkDir):
            outdir: str = self.outdir.src
        else:
            outdir: str = self.outdir
        
        denoisedir: str = os.path.join(outdir,'denoise')

        with WorkDir(src=denoisedir) as od:
            fixdir: str = od.join('fix')
            self.output: Dict[str,str] = {
                "fix_labels": od.join("fix_labels.txt"),
                "fix_labels": od.join('fix_labels.txt'),
                "func_clean": od.join('func_clean.nii.gz'),
                "func_clean_mean": od.join('func_clean_mean.nii.gz'),
                "func_clean_std": od.join('func_clean_std.nii.gz'),
                "func_clean_tsnr": od.join('func_clean_tsnr.nii.gz'),
                "fix_clean": os.path.join(fixdir,'filtered_func_data_clean.nii.gz')
            }
        return self.output
