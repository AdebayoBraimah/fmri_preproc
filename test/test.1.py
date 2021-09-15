#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import os
# import sys
# _mod_path: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc"
# sys.path.append(_mod_path)

from fmri_preproc.func.pipeline import Pipeline

# Global variables

# Echospacing (ms) calculation
# 
# echospacing = 1/1844.5 * 1000
# 
# echospacing = 0.54 # (ms)
# echospacing/1000 = 0.00054 (sec)

# Directory variables
# outdir: str = "/users/brac4g/fmri_neonate/FMRI_test_1/test/test.1"
outdir: str = "/users/brac4g/fmri_neonate/FMRI_test_1/test/test.2"

# # (Local) File variables
# outdir: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test/test.1"
# 
# func: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_bold.nii.gz"
# sbref: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz"
# 
# T2w: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_T2w_restore.nii.gz"
# T2w_brainmask: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_brainmask_bet.nii.gz"
# dseg: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_drawem_tissue_labels.nii.gz"
# 
# ap: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-AP_run-01_sbref.nii.gz"
# pa: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz"

# Cluster variables
func: str = "/users/brac4g/fmri_neonate/test_data/func/func/sub-194_task-rest_dir-PA_run-01_bold.nii.gz"
sbref: str = "/users/brac4g/fmri_neonate/test_data/func/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz"

T2w: str = "/users/brac4g/fmri_neonate/test_data/struct/ses-session1/anat/sub-1_ses-1_T2w_restore.nii.gz"
T2w_brainmask: str = "/users/brac4g/fmri_neonate/test_data/struct/ses-session1/anat/sub-1_ses-1_brainmask_bet.nii.gz"
dseg: str = "/users/brac4g/fmri_neonate/test_data/struct/ses-session1/anat/sub-1_ses-1_drawem_tissue_labels.nii.gz"

ap: str = "/users/brac4g/fmri_neonate/test_data/func/func/sub-194_task-rest_dir-AP_run-01_sbref.nii.gz"
pa: str = sbref

# Rdata (for testing purposes)
rdata: str = "/usr/local/fix/1.6.15/training_files/Standard.RData"

# Numerical and misc. variables
scan_pma: float = 44
echospacing: float = 0.54/1000
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

# preproc.import_data()

# preproc.run_all(use_mcflirt=False, s2v=True, rdata=rdata, fix_threshold=20)
# preproc.run_all(use_mcflirt=True, s2v=False, rdata=rdata, fix_threshold=20)

# preproc.run_all(use_mcflirt=False, s2v=True, rdata=rdata, fix_threshold=20, template_ages='neo')
# preproc.run_all(use_mcflirt=True, s2v=False, rdata=rdata, fix_threshold=20, template_ages='neo')

preproc.run_all(use_mcflirt=False, s2v=True, rdata=rdata, fix_threshold=20, template_ages=['neo', 1, 2, 38, 44])
# preproc.run_all(use_mcflirt=True, s2v=False, rdata=rdata, fix_threshold=20, template_ages=['neo', 1, 2, 38, 44])
