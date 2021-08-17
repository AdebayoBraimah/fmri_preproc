#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
_mod_path: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc"
sys.path.append(_mod_path)

from fmri_preproc.utils.logutil import LogFile

log: LogFile = LogFile("test.1.log")

# global variables
outdir = 'test.1'
ap = '/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-AP_run-01_sbref.nii.gz'
pa = '/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz'


# Test field map
from fmri_preproc.func.fieldmap import fieldmap

fieldmap(outdir=outdir,
         ap_dir=ap,
         pa_dir=pa,
         echo_spacing=0.1,
         log=log)