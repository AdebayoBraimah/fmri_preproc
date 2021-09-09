# -*- coding: utf-8 -*-
"""Output data filepaths dictionaries for the ``registration`` module.
"""
import os

from typing import (
    Dict,
    Union
)

from fmri_preproc.utils.outputs.outdict import OutDict
from fmri_preproc.utils.workdir import WorkDir


class MRreg(OutDict):
    """class doc-string.
    """
    def __init__(self,
                 outdir: Union[WorkDir,str]
                ) -> None:
        """Class constructor method.
        NOTE: ``outdir`` should be subject working directory.
        """
        self.outdir: Union[WorkDir,str] = outdir
        super(MRreg, self).__init__()
    
    def outputs(self,
                src_space: str,
                ref_space: str,
               ) -> Dict[str,str]:
        """doc-string.
        """
        if isinstance(self.outdir,WorkDir):
            outdir: str = self.outdir.src
        else:
            outdir: str = self.outdir
        
        with WorkDir(src=outdir) as od:
            regdir: str = od.join('reg',f'{src_space}_to_{ref_space}')

            if os.path.exists(regdir):
                remove_dir: bool = False
            else:
                remove_dir: bool = True

            with WorkDir(src=regdir) as rd:
                regname: str = rd.join(f'{src_space}_to_{ref_space}')

                self.output: Dict[str,str] = {
                    "resampled_image": f"{regname}_img.nii.gz",
                    "warp": f"{regname}_warp.nii.gz",
                    "inv_warp": f"{regname}_invwarp.nii.gz",
                    "affine": f"{regname}_affine.mat",
                    "inv_affine": f"{regname}_invaffine.mat",
                    "init_affine": f"{regname}_init_affine.mat",
                    "resampled_image_init": f"{regname}_init_img.nii.gz",
                    "resampled_brainmask": f"{regname}_brainmask.nii.gz",
                    "dc_warp": rd.join(f'{src_space}_dc_warp.nii.gz'),
                    "dc_image": rd.join(f'{src_space}_dc_img.nii.gz'),
                    "dc_brainmask": rd.join(f'{src_space}_dc_brainmask.nii.gz'),
                }

                if remove_dir:
                    rd.rmdir(rm_parent=False)
        return self.output
