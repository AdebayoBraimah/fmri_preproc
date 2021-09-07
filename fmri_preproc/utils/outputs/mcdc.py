# -*- coding: utf-8 -*-
"""Output data filepaths dictionaries for the ``mcdc`` module.
"""
import os

from typing import (
    Dict,
    Union
)

from fmri_preproc.utils.outputs.outdict import OutDict
from fmri_preproc.utils.workdir import WorkDir

class MCDCFiles(OutDict):
    """class doc-string
    """
    def __init__(self,
                 outdir: Union[WorkDir,str]
                ) -> None:
        """Class constructor method.
        NOTE: ``outdir`` should be subject working directory.
        """
        self.outdir: Union[WorkDir,str] = outdir
        super(MCDCFiles, self).__init__()
    
    def outputs(self,
                dc: bool = False) -> Dict[str,str]:
        """doc-string.
        """
        if isinstance(self.outdir,WorkDir):
            outdir: str = self.outdir.src
        else:
            outdir: str = self.outdir
        
        if dc:
            mc_name: str = 'mcdc'
        else:
            mc_name: str = 'mc'
        
        mcdir: str = os.path.join(outdir,f"{mc_name}")
        
        with WorkDir(src=mcdir) as mc:
            self.output: Dict[str,str] = {
                "mcdir": mc.abspath(),
                "func_mcdc": mc.join(f"func_{mc_name}.nii.gz"),
                "func_mot": mc.join(f"func_{mc_name}_motion.tsv"),
                "func_metrics": mc.join(f"func_{mc_name}_regressors.tsv"),
                "func_out_plot": mc.join(f"func_{mc_name}_outliers.png"),
                "mcdc_mean": mc.join(f"func_{mc_name}_mean.nii.gz"),
                "mcdc_std": mc.join(f"func_{mc_name}_std.nii.gz"),
                "mcdc_tsnr": mc.join(f"func_{mc_name}_tsnr.nii.gz"),
                "mcdc_brainmask": mc.join(f"func_{mc_name}_brainmask.nii.gz"),
                "func_mcdc_fovmask": mc.join(f"func_{mc_name}_fovmask.nii.gz"),
                "func_mcdc_fovpercent": mc.join(f"func_{mc_name}_fovpercent.nii.gz")
            }

        return self.output
