#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

_pkg_path: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc"
sys.path.append(_pkg_path)

from fmri_preproc.func.pipeline import Pipeline

# directories
outdir: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test/test_proc"

# data
func: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test_data/fmri_preproc/sub-148/run-01/import/func.nii.gz"
age: int = 44
config: str = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test/settings.test.json"

func_pedir = "pa"
T2w = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test_data/fmri_preproc/sub-148/run-01/import/T2w.nii.gz"
T2w_brainmask = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test_data/fmri_preproc/sub-148/run-01/import/T2w_brainmask.nii.gz"
dseg = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test_data/fmri_preproc/sub-148/run-01/import/T2w_dseg.nii.gz"

func_echospacing = 0.00054
mb_factor = 3
sbref = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test_data/fmri_preproc/sub-148/run-01/import/sbref.nii.gz"
sbref_pedir = "pa"

spinecho = "/Users/adebayobraimah/Desktop/projects/fmri_preproc/test_data/fmri_preproc/sub-148/run-01/fmap/topup/spinecho.nii.gz"
spinecho_pedir = "pa,ap"


# pipeline
preproc: Pipeline = Pipeline(outdir=outdir,
    subid='148',
    runid='01',
    scan_pma=age,
    settings_json_dict=config)

preproc.import_data(func=func,
    func_echospacing=func_echospacing,
    func_pedir=func_pedir,
    T2w=T2w,
    T2w_brainmask=T2w_brainmask,
    dseg=dseg,
    func_inplane_accel=1.0,
    mb_factor=mb_factor,
    sbref=sbref,
    sbref_pedir=sbref_pedir,
    spinecho=spinecho,
    spinecho_pedir=spinecho_pedir,
    spinecho_inplaneacc=1.0,
    sbref_echospacing=func_echospacing,
    spinecho_echospacing=func_echospacing)

preproc.run_all(use_mcflirt=True,
    quick=True,
    icadim=10,
    group_map=None,
    preproc_only=True,
    spatial_fwhm=1.5,
    intnorm=True)
