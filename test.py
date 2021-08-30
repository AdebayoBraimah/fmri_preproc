#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
_mod_path: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc"
sys.path.append(_mod_path)

from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.func.fieldmap import fieldmap
from fmri_preproc.func.mcdc import mcdc

from fmri_preproc.func.importdata import (
     import_info,
     import_func,
     import_struct,
     import_spinecho
)

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

# Test import info
sub_workdir, logdir, sub_json, sub_dict = import_info(outdir=outdir,
                                                      func=func,
                                                      scan_pma=scan_pma,
                                                      log_datetime=True)

# Test import functional data
_import_log: str = os.path.join(logdir,'import.log')

import_log: LogFile = LogFile(_import_log,
                              format_log_str=True)

(func,
func_mean, 
func_brainmask, 
func_slorder, 
sbref, 
sbref_brainmask) = import_func(outdir=sub_workdir,
                               func=func,
                               func_echospacing=echospacing,
                               func_pedir=func_pedir,
                               mb_factor=mb_factor,
                               sbref=sbref,
                               sbref_pedir=func_pedir,
                               sbref_echospacing=echospacing,
                               log=import_log)

T2w, wmmask, dseg = import_struct(outdir=sub_workdir,
                                  T2w=T2w,
                                  brainmask=T2w_brainmask,
                                  dseg=dseg,
                                  dseg_type='drawem',
                                  log=import_log)

spinecho, spinecho_pedir = import_spinecho(outdir=sub_workdir,
                                           ap_dir=ap,
                                           pa_dir=pa,
                                           spinecho_epifactor=63,
                                           spinecho_inplaneacc=inplane_acc,
                                           log=import_log)

# Test processing fmap data
_fmap_log: str = os.path.join(logdir,'fieldmap.log')

fmap_log: LogFile = LogFile(_fmap_log,
                              format_log_str=True)

fmap, mag, fmap_mask = fieldmap(outdir=sub_workdir,
                                spinecho=spinecho,
                                echo_spacing=echospacing,
                                pedir=spinecho_pedir,
                              #   epifactor=epifactor,
                                inplane_acc=inplane_acc,
                                log=fmap_log)

# Test processing fmap data
_mcdc_log: str = os.path.join(logdir,'mcdc.log')

mcdc_log: LogFile = LogFile(_mcdc_log,
                              format_log_str=True)

(func_mcdc,
mcf,
func_metrics,
mcdc_mean,
mcdc_std,
mcdc_tsnr,
mcdc_brainmask,
fov_mask,
fov_percent) = mcdc(func=func,
                    outdir=sub_workdir,
                    func_echospacing=echospacing,
                    func_pedir="PA",
                    inplane_acc=inplane_acc,
                    func_brainmask=func_brainmask,
                    fmap=fmap,
                    mb_factor=mb_factor,
                    dc=True,
                    s2v=True,
                    mbs=True,
                    use_mcflirt=False,
                    log=mcdc_log)
