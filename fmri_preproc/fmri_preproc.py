#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Performs preprocessing of neonatal resting-state functional MR neuroimages.
"""
# TODO: 
#   * Write CLI
#   * Create settings function(s)/module and corresponding settings(.json) file.

import os
import argparse

from typing import (
    Optional,
    Tuple
)

# from fmri_preproc.func.pipeline import Pipeline


def main() -> int:
    """Main function.
    """
    
    args, parser = arg_parser()

    return 0


def preproc_data():
    """Preprocess rs-fMRI data.
    """
    # preproc: Pipeline = Pipeline(outdir=outdir,
    #                              func=func,
    #                              scan_pma=44)
    pass


def arg_parser() -> Tuple[argparse.ArgumentParser.parse_args, argparse.ArgumentParser]:
    """Argument parser.
    """
    # Init parser
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=__doc__)#, formatter_class=argparse.RawTextHelpFormatter)

    # Parse Arguments

    # Required Arguments
    reqoptions = parser.add_argument_group('Required Arguments')
    reqoptions.add_argument('-o','--outdir',
                            type=str,
                            metavar="<path>",
                            dest="outdir",
                            help="REQUIRED: Output location for working directory.")
    reqoptions.add_argument('-f','--func',
                            type=str,
                            metavar="<file>",
                            dest="func",
                            help="REQUIRED: Input rs-fMRI data. NOTE: Input data MUST conform to the BIDS naming convention (e.g. 'sub-001[_ses-001]_task-rest[_run-01]_bold.nii.gz').")
    reqoptions.add_argument('-a','--age',
                            type=str,
                            metavar="<age>",
                            dest="age",
                            help="REQUIRED: Age at scan. Acceptable inputs include: 'neo', 28 - 44 (wks), and 1 or 2 (yrs).")

    importoptions = parser.add_argument_group('Import Data Arguments')

    # Structural import data
    importoptions.add_argument('--T2w',
                               type=str,
                               metavar="<file>",
                               dest="T2w",
                               help="Input (bias corrected, e.g. '_restore') T2w image data.")
    importoptions.add_argument('--T2w-brainmask',
                               type=str,
                               metavar="<file>",
                               dest="T2w_brainmask",
                               help="Input T2w brainmask image data (_brainmask_bet).")
    importoptions.add_argument('--T2w-tissue-seg',
                               type=str,
                               metavar="<file>",
                               dest="dseg",
                               help="Input Draw-EM neonatal brain tissue segmentation (e.g. '_drawem_tissue_labels') image data. NOTE: Other brain tissue files can be used (e.g. from FreeSurfer or FSL's FAST).")
    importoptions.add_argument('--tissue-seg-type',
                               type=str,
                               metavar="<str>",
                               dest="dseg_type",
                               help="Input brain tissue segmentation type. Acceptable inputs include: 'drawem', 'freesurfer_aseg', or 'fsl_fast'.")
    importoptions.add_argument('--T1w',
                               type=str,
                               metavar="<file>",
                               dest="T1w",
                               help="Input (bias corrected, e.g. '_restore') T1w image data.")
    
    # Functional import data
    importoptions.add_argument('--func-echospacing',
                               type=float,
                               metavar="<float>",
                               dest="func_echospacing",
                               help="Functional echospacing.")
    importoptions.add_argument('--func-pedir',
                               type=str,
                               metavar="<PE-dir>",
                               dest="func_pedir",
                               help="Phase encoding direction. Acceptable inputs include: 'PA', 'AP', 'RL', 'LR', 'IS', or 'SI'.")
    importoptions.add_argument('--func-slorder',
                               type=str,
                               metavar="<file>",
                               dest="func_slorder",
                               help="Functional slice order specification file. This file specifies the slice acquisition order of the input rs-fMRI. If not provided one can be created internally IF 'func-mbfactor' is specified.")
    importoptions.add_argument('--func-inplane-accel',
                               type=float,
                               metavar="<float>",
                               dest="func_inplane_accel",
                               default=1.00,
                               help="Functional inplane acceleration (e.g. SENSE factor on Philips MR scanners or GRAPPA, the general term) [default=1.00].")
    importoptions.add_argument('--func-mbfactor',
                               type=int,
                               metavar="<int>",
                               dest="mb_factor",
                               help="Functional multi-band factor.")
    
    # Single-band reference import data
    importoptions.add_argument('--sbref',
                               type=str,
                               metavar="<file>",
                               dest="sbref",
                               help="Single band reference MR image data.")
    importoptions.add_argument('--sbref-brainmask',
                               type=str,
                               metavar="<file>",
                               dest="sbref_brainmask",
                               help="Single band reference MR brain mask image data.")
    importoptions.add_argument('--sbref-echospacing',
                               type=float,
                               metavar="<float>",
                               dest="sbref_echospacing",
                               help="Single band reference MR image echospacing. If not specified, then it is assumed to match '--func-echospacing'.")
    importoptions.add_argument('--sbref-pedir',
                               type=str,
                               metavar="<PE-dir>",
                               dest="sbref_pedir",
                               help="Single band reference MR image data phase encoding direction. NOTE: This should match '--func-pedir'. If not specified, then it is assumed to matched '--func-pedir'. If this is not case, this could result in mis-registration and distortion correction issues.")
    
    # Single-band reference import data
    importoptions.add_argument('--spinecho',
                               type=str,
                               metavar="<file>",
                               dest="spinecho",
                               help="Spinecho MR image data.")
    importoptions.add_argument('--sp-echospacing',
                               type=float,
                               metavar="<float>",
                               dest="spinecho_echospacing",
                               help="Spinecho MR image echospacing. If not specified, then it is assumed to match '--func-echospacing'.")
    importoptions.add_argument('--sp-pedir',
                               type=str,
                               metavar="<PE-dir,PE-dir,..,PE-dir>",
                               dest="spinecho_pedir",
                               help="Spinecho MR image phase encoding direction(s), specified as a comma separated list. Accepatable inputs are the same as those specifed for both '--func-pedir', and '--sp-pedir'.")
    importoptions.add_argument('--sp-inplane-accel',
                               type=float,
                               metavar="<float>",
                               dest="spinecho_inplaneacc",
                               help="Spinecho MR image inplane acceleration . If not specified, then it is assumed to match '--func-inplane-accel'.")
    importoptions.add_argument('--sp-AP',
                               type=str,
                               metavar="<file>",
                               dest="ap_dir",
                               help="A -> P (anterior to posterior) Spinecho MR image data.")
    importoptions.add_argument('--sp-PA',
                               type=str,
                               metavar="<file>",
                               dest="pa_dir",
                               help="P -> A (posterior to anterior) Spinecho MR image data.")
    importoptions.add_argument('--sp-LR',
                               type=str,
                               metavar="<file>",
                               dest="lr_dir",
                               help="L -> R (left to right) Spinecho MR image data.")
    importoptions.add_argument('--sp-RL',
                               type=str,
                               metavar="<file>",
                               dest="rl_dir",
                               help="R -> L (right to left) Spinecho MR image data.")
    importoptions.add_argument('--sp-IS',
                               type=str,
                               metavar="<file>",
                               dest="is_dir",
                               help="I -> S (inferior to superior) Spinecho MR image data.")
    importoptions.add_argument('--sp-SI',
                               type=str,
                               metavar="<file>",
                               dest="si_dir",
                               help="S -> I (superior to inferior) Spinecho MR image data.")

    mcdcoptions = parser.add_argument_group('Motion Correction and Distortion Correction Arguments')

    mcdcoptions.add_argument('--s2v', '--slice-to-volume',
                             dest="s2v",
                             action="store_true",
                             help="Enables slice-to-volume motion correction. NOTE: This will LIKELY require EDDY_CUDA and GPU processing. This is possible without GPUs in FSL v6.0.5+ but is VERY SLOW.")
    mcdcoptions.add_argument('--dc', '--distortion-correction',
                             dest="dc",
                             action="store_true",
                             help="Enables distortion correction.")
    mcdcoptions.add_argument('--mbs',
                             dest="mbs",
                             action="store_true",
                             help="Enables movement by susceptibility distortion correction.")
    mcdcoptions.add_argument('--mcflirt',
                             dest="use_mcflirt",
                             action="store_true",
                             help="Enables FSL's MCFLIRT-based motion correction. NOTE: Enabling this option DISABLES: '--s2v', '--dc', and '--mbs' options.")

    regoptions = parser.add_argument_group('Template Registration Arguments')

    regoptions.add_argument('--standard-age',
                            type=int,
                            metavar="<int>",
                            dest="standard_age",
                            help="Standard age to use for final dHCP age-matched template (PMA, in weeks) [default=40 (wks)].")
    regoptions.add_argument('--template-ages',
                            type=str,
                            metavar="<str,str,..str>",
                            dest="template_ages",
                            help="Comma separated list of ages. Acceptable inputs include: 'neo', 1, 2 (yrs), and/or 28 - 44 (PMA, wks).")
    regoptions.add_argument('--atlas-dir',
                            type=str,
                            metavar="<path>",
                            dest="atlasdir",
                            help="Atlas directory path.")

    icaoptions = parser.add_argument_group('ICA Arguments')

    icaoptions.add_argument('--temporal-fwhm',
                            type=float,
                            metavar="<float>",
                            dest="temporal_fwhm",
                            help="Temporal FWHM filtering coefficient [default=150.0 (s)].")
    icaoptions.add_argument('--ica-dim',
                            type=int,
                            metavar="<int>",
                            dest="icadim",
                            help="Maximal ICA dimensionality.")

    fixoptions = parser.add_argument_group('Denoise/FIX Arguments')

    fixoptions.add_argument('--rdata',
                            type=str,
                            metavar="<str>",
                            dest="rdata",
                            help="Trained FIX classifier (stored as a Rdata file).")
    fixoptions.add_argument('--fix-threshold',
                            type=str,
                            metavar="<str>",
                            dest="rdata",
                            help="FIX noise classification threshold. [default= TODO: SET LATER].")

    pipeoptions = parser.add_argument_group('Pipeline Runtime Behavior Arguments')

    pipeoptions.add_argument('--pre-mcdc',
                             dest="pre_mcdc",
                             action="store_true",
                             help="Perform pipeline preprocessing steps prior to motion correction and distortion correction.")
    pipeoptions.add_argument('--mcdc',
                             dest="mcdc",
                             action="store_true",
                             help="Perform motion correction and distortion correction step exclusively.")
    pipeoptions.add_argument('--post-mcdc',
                             dest="post_mcdc",
                             action="store_true",
                             help="Perform pipeline preprocessing steps after motion correction and distortion correction.")

    optoptions = parser.add_argument_group('Optional Arguments')

    optoptions.add_argument('--config', '--settings',
                            type=str,
                            metavar="<file>",
                            dest="settings",
                            help="Configuration (JSON) file that contains additional input arguments not specified at the command line.")


    args: argparse.ArgumentParser.parse_args = parser.parse_args()

    return args, parser


if __name__ == '__main__':
    main()
