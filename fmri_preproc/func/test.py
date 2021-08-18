#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import os
_mod_path: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc"
sys.path.append(_mod_path)

from fmri_preproc.utils.logutil import LogFile

# # global variables
# outdir: str = '/Users/adebayobraimah/Desktop/projects/fmri_preproc/test.code/test.1'
# log: LogFile = LogFile(os.path.join(outdir),"test.1.log")
# ap: str = '/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-AP_run-01_sbref.nii.gz'
# pa: str = '/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz'
# 
# # Test field map
# from fmri_preproc.func.fieldmap import fieldmap
# 
# fieldmap(outdir=outdir,
#          ap_dir=ap,
#          pa_dir=pa,
#          echo_spacing=0.1,
#          log=log)

# # global variables
# outdir: str = '/Users/adebayobraimah/Desktop/projects/fmri_preproc/test.code/test.2'
# log: LogFile = LogFile("/Users/adebayobraimah/Desktop/projects/fmri_preproc/test.code/test.2.log")
# func: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_bold.nii.gz"
# func_mask: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test.code/func_brain_mask.nii.gz"
# 
# # Test MCDC (mcflirt)
# from fmri_preproc.func.mcdc import mcdc
# 
# mcdc(func=func,
#      outdir=outdir,
#      func_echospacing=0.1,
#      func_brainmask=func_mask,
#      use_mcflirt=True,
#      log=log)

# global variables
outdir: str = '/Users/adebayobraimah/Desktop/projects/fmri_preproc/test.code/test.3'
log: LogFile = LogFile("/Users/adebayobraimah/Desktop/projects/fmri_preproc/test.code/test.3.log")
func: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_bold.nii.gz"
fmap: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test.code/test.1/fmap.nii.gz"
func_mask: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test.code/func_brain_mask.nii.gz"

# Test MCDC (mcflirt)
from fmri_preproc.func.mcdc import mcdc

mcdc(func=func,
     outdir=outdir,
     fmap=fmap,
     func_echospacing=0.1,
     func_brainmask=func_mask,
     mb_factor=5,
     s2v=True,
     mbs=True,
     use_gpu=True,
     use_mcflirt=False,
     log=log)