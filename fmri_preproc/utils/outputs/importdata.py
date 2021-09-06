# -*- coding: utf-8 -*-
"""Output data filepaths dictionaries for the ``importdata`` module.
"""
import os

from typing import (
    Dict,
    Union
)

from fmri_preproc.utils.outputs.outdict import OutDict
from fmri_preproc.utils.workdir import WorkDir


class ImportFunc(OutDict):
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
                                "func": od.join('func.nii.gz'),
                                "func_mean": od.join('func_mean.nii.gz'),
                                "func_brainmask": od.join('func_brainmask.nii.gz'),
                                "func_slorder": od.join('func.slorder'),
                                "sbref": od.join('sbref.nii.gz'),
                                "sbref_brainmask": od.join('sbref_brainmask.nii.gz'),
                                "func_std": od.join('func_std.nii.gz'),
                                "func_tsnr": od.join('func_tsnr.nii.gz'),
                                "func0": od.join('func0.nii.gz'),
                                "func_dil_brainmask": od.join('func_dilated_brainmask.nii.gz'),
                                "func_metrics": od.join('func_regressors.tsv'),
                                "func_out_plot": od.join('func_outliers.png'),
                             }
        return self.output


class ImportStruct(OutDict):
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
                "T1w": od.join('T1w.nii.gz'),
                "T2w": od.join('T2w.nii.gz'),
                "T2w_brainmask": od.join('T2w_brainmask.nii.gz'),
                "T2w_dseg": od.join('T2w_dseg.nii.gz'),
                "T2w_wmmask": od.join('T2w_wmmask.nii.gz'),
                "T2w_probseg": od.join('T2w_probseg.nii.gz'),
                }
        return self.output


class ImportSpinEcho(OutDict):
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
            topup_dir: str = od.join('fmap','topup')

            if os.path.exists(topup_dir): 
                remove_dir: bool = False
            else:
                remove_dir: bool = True
            
            with WorkDir(src=topup_dir) as td:
                self.output: Dict[str,str] = {"spinecho": td.join('spinecho.nii.gz')}
                if remove_dir:
                    td.rmdir(rm_parent=True)
        return self.output
