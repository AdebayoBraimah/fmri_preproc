#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from logging import log
import os
import sys
_mod_path: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc"
sys.path.append(_mod_path)

from fmri_preproc.func.pipeline import Pipeline

# Global variables

# Directory variables
outdir: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test/test.1"

# File variables
func: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_bold.nii.gz"
sbref: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz"

T2w: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_T2w_restore.nii.gz"
T2w_brainmask: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_brainmask_bet.nii.gz"
dseg: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_drawem_tissue_labels.nii.gz"

ap: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-AP_run-01_sbref.nii.gz"
pa: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz"

# Numerical and misc. variables
scan_pma: float = 44
echospacing: float = 0.54
func_pedir: str = "pa"
mb_factor: int = 3
epifactor: int = 63
inplane_acc: float = 1

preproc: Pipeline = Pipeline(outdir=outdir,
                             func=func,
                             scan_pma=44)

preproc.import_data(func_echospacing=echospacing,
                    func_pedir=func_pedir,
                    func_inplane_accel=inplane_acc,
                    T2w=T2w,
                    T2w_brainmask=T2w_brainmask,
                    dseg=dseg,
                    mb_factor=mb_factor,
                    sbref=sbref,
                    sbref_pedir=func_pedir,
                    sbref_echospacing=echospacing,
                    ap_dir=ap,
                    pa_dir=pa,
                    spinecho_inplaneacc=inplane_acc)

preproc.run_all(use_mcflirt=True)