#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# Script dir
scripts_dir=$(realpath $(dirname $(realpath ${0}))/..)
fmri_preproc=${scripts_dir}/fmri_preproc.py

# echo ${scripts_dir}

# Add to python path
export PYTHONPATH=${PYTHONPATH}:${scripts_dir}

# Define variables

# # Local
# outdir="/Users/adebayobraimah/Desktop/projects/fmri_preproc/test/TEST001/test.1"
# 
# func="/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_bold.nii.gz"
# sbref="/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz"
# 
# T2w="/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_T2w_restore.nii.gz"
# T2w_brainmask="/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_brainmask_bet.nii.gz"
# dseg="/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/structural/sub-194/ses-session1/anat/sub-1_ses-1_drawem_tissue_labels.nii.gz"
# 
# ap="/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-AP_run-01_sbref.nii.gz"
# pa="/Users/adebayobraimah/Desktop/projects/fmri_preproc/data/sub-194/functional/sub-194/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz"
# 
# config="/Users/adebayobraimah/Desktop/projects/fmri_preproc/test/settigs.test.json"

# HPC
outdir="/users/brac4g/fmri_neonate/FMRI_test_1/test/test.3"

func="/users/brac4g/fmri_neonate/test_data/func/func/sub-194_task-rest_dir-PA_run-01_bold.nii.gz"
sbref="/users/brac4g/fmri_neonate/test_data/func/func/sub-194_task-rest_dir-PA_run-01_sbref.nii.gz"

T2w="/users/brac4g/fmri_neonate/test_data/struct/ses-session1/anat/sub-1_ses-1_T2w_restore.nii.gz"
T2w_brainmask="/users/brac4g/fmri_neonate/test_data/struct/ses-session1/anat/sub-1_ses-1_brainmask_bet.nii.gz"
dseg="/users/brac4g/fmri_neonate/test_data/struct/ses-session1/anat/sub-1_ses-1_drawem_tissue_labels.nii.gz"

ap="/users/brac4g/fmri_neonate/test_data/func/func/sub-194_task-rest_dir-AP_run-01_sbref.nii.gz"
pa=${sbref}

config=$(realpath /users/brac4g/fmri_neonate/fmri_preproc/test/settigs.test.json)

# Rdata (for testing purposes)
rdata="/usr/local/fix/1.6.15/training_files/Standard.RData"

# Numerical and misc. variables
scan_pma=44
echospacing=0.00054
func_pedir="pa"
mb_factor=3
epifactor=63
inplane_acc=1

# Run
${fmri_preproc} \
run-all \
--outdir=${outdir} \
--func=${func} \
--age=${scan_pma} \
--T2w=${T2w} \
--T2w-brainmask=${T2w_brainmask} \
--T2w-tissue-seg=${dseg} \
--tissue-seg-type=drawem \
--func-echospacing=${echospacing} \
--func-inplane-accel=${inplane_acc} \
--func-pedir=PA \
--func-mbfactor=${mb_factor} \
--sbref=${sbref} \
--sp-AP=${ap} \
--sp-PA=${pa} \
--config=$(realpath ${config})
