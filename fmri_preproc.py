#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Performs preprocessing of neonatal resting-state functional MR neuroimages.
"""
import sys
import argparse

from typing import (
    Any,
    Dict,
    Literal,
    Tuple
)

from fmri_preproc import (
    DEFAULT_SETTINGS_FILE,
    __version__ as version_id
)

from fmri_preproc.func.pipeline import Pipeline

from fmri_preproc.utils.util import (
    dict2json,
    settings
)


def main() -> Literal[0]:
    """Main function.
    """
    preproc_data()
    return 0


def preproc_data() -> Literal[0]:
    """Preprocess rs-fMRI data.
    """
    args, parser = arg_parser()

    # Print help message in the case of no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    # Version
    if args.version:
        print(f"\nfmri_preproc: v{version_id}\n")
        sys.exit(0)

    # Write settings file
    if args.out_settings:
        settings_dict: Dict[str,Any] = settings(jsonfile=DEFAULT_SETTINGS_FILE)
        _: str = dict2json(dict=settings_dict, jsonfile=args.out_settings, indent=4)
        sys.exit(0)
    
    # Preprocess data

    # Check arguments
    if args.outdir and ((args.sub and args.run) or args.func):

        # Confirm phase-encoding directions
        if not args.sbref_pedir:
            args.sbref_pedir = args.func_pedir
        
        # Confirm echo-spacing paramters
        if not args.sbref_echospacing:
            args.sbref_echospacing = args.func_echospacing
        
        if not args.spinecho_echospacing:
            args.spinecho_echospacing = args.func_echospacing
        
        if not args.spinecho_inplaneacc:
            args.spinecho_inplaneacc = args.func_inplane_accel
        
        # Configure runtime settings
        settings_dict: Dict[str,Any] = settings(jsonfile=args.settings,
                                                func_inplane_accel=args.func_inplane_accel,
                                                func_slorder=args.func_slorder,
                                                mb_factor=args.mb_factor,
                                                sbref_echospacing=args.sbref_echospacing,
                                                dseg_type=args.dseg_type,
                                                probseg_type=None, # Set this to None for now
                                                mask_func=None, # Set this to None for now
                                                spinecho_echospacing=args.spinecho_echospacing,
                                                spinecho_epifactor=None, # Set this to None for now
                                                spinecho_inplaneacc=args.spinecho_inplaneacc,
                                                use_mcflirt=args.use_mcflirt,
                                                s2v=args.s2v,
                                                dc=args.dc,
                                                mbs=args.mbs,
                                                standard_age=args.standard_age,
                                                quick=None,
                                                atlasdir=args.atlasdir,
                                                template_ages=args.template_ages,
                                                temporal_fwhm=args.temporal_fwhm,
                                                icadim=args.icadim,
                                                rdata=args.rdata,
                                                fix_threshold=args.fix_threshold,
                                                group_qc=None,
                                                group_map=None,
                                                standard_res=None,
                                                verbose=args.verbose,
                                                log_level=args.log_level,
                                                smooth=args.smooth,
                                                intnorm=args.intnorm)
        
        # Get subject info
        preproc: Pipeline = Pipeline(outdir=args.outdir, 
                                     subid=args.subid,
                                     sesid=args.sesid,
                                     runid=args.runid,
                                     func=args.func, 
                                     scan_pma=args.age, 
                                     birth_ga=args.birth_age, 
                                     verbose=settings_dict.get('verbose'), 
                                     log_level=settings_dict.get('log_level'), 
                                     settings_json_dict=settings_dict)
        if args.method == 'import':
            preproc.import_data(func=args.func,
                                func_echospacing=args.func_echospacing,
                                func_pedir=args.func_pedir,
                                T2w=args.T2w,
                                T2w_brainmask=args.T2w_brainmask,
                                dseg=args.dseg,
                                func_brainmask=args.func_brainmask,
                                func_slorder=settings_dict.get('func_slorder'),
                                func_inplane_accel=settings_dict.get('func_inplane_accel'),
                                mb_factor=settings_dict.get('mb_factor'),
                                sbref=args.sbref,
                                sbref_brainmask=args.sbref_brainmask,
                                sbref_echospacing=settings_dict.get('sbref_echospacing'),
                                sbref_pedir=args.sbref_pedir,
                                dseg_type=settings_dict.get('dseg_type'),
                                T1w=args.T1w,
                                spinecho=args.spinecho,
                                spinecho_echospacing=settings_dict.get('spinecho_echospacing'),
                                spinecho_pedir=args.spinecho_pedir,
                                ap_dir=args.ap_dir,
                                pa_dir=args.pa_dir,
                                lr_dir=args.lr_dir,
                                rl_dir=args.rl_dir,
                                is_dir=args.is_dir,
                                si_dir=args.si_dir,
                                spinecho_epifactor=settings_dict.get('spinecho_epifactor'),
                                spinecho_inplaneacc=settings_dict.get('spinecho_inplaneacc'))
        elif args.method == 'pre_mcdc':
            preproc.pre_mcdc()
        elif args.method == 'mcdc':
            preproc.mcdc(use_mcflirt=settings_dict.get('use_mcflirt'),
                         s2v=settings_dict.get('s2v'),
                         dc=settings_dict.get('dc'),
                         mbs=settings_dict.get('mbs'))
        elif args.method == 'post_mcdc':
            preproc.post_mcdc(standard_age=settings_dict.get('standard_age'),
                              template_ages=settings_dict.get('template_ages'),
                              temporal_fwhm=settings_dict.get('temporal_fwhm'),
                              quick=settings_dict.get('quick'),
                              icadim=settings_dict.get('icadim'),
                              rdata=settings_dict.get('rdata'),
                              fix_threshold=settings_dict.get('fix_threshold'),
                              group_map=settings_dict.get('group_map'),
                              standard_res=settings_dict.get('standard_res'),
                              group_qc=settings_dict.get('group_qc'),
                              atlasdir=settings_dict.get('atlasdir'),
                              preproc_only=settings_dict.get('preproc_only'),
                              smooth=settings_dict.get('smooth'),
                              intnorm=settings_dict.get('intnorm'))
        elif args.method == 'run_all':
            # Import data
            preproc.import_data(func=args.func,
                                func_echospacing=args.func_echospacing,
                                func_pedir=args.func_pedir,
                                T2w=args.T2w,
                                T2w_brainmask=args.T2w_brainmask,
                                dseg=args.dseg,
                                func_brainmask=args.func_brainmask,
                                func_slorder=settings_dict.get('func_slorder'),
                                func_inplane_accel=settings_dict.get('func_inplane_accel'),
                                mb_factor=settings_dict.get('mb_factor'),
                                sbref=args.sbref,
                                sbref_brainmask=args.sbref_brainmask,
                                sbref_echospacing=settings_dict.get('sbref_echospacing'),
                                sbref_pedir=args.sbref_pedir,
                                dseg_type=settings_dict.get('dseg_type'),
                                T1w=args.T1w,
                                spinecho=args.spinecho,
                                spinecho_echospacing=settings_dict.get('spinecho_echospacing'),
                                spinecho_pedir=args.spinecho_pedir,
                                ap_dir=args.ap_dir,
                                pa_dir=args.pa_dir,
                                lr_dir=args.lr_dir,
                                rl_dir=args.rl_dir,
                                is_dir=args.is_dir,
                                si_dir=args.si_dir,
                                spinecho_epifactor=settings_dict.get('spinecho_epifactor'),
                                spinecho_inplaneacc=settings_dict.get('spinecho_inplaneacc'))

            # RUN-ALL
            preproc.run_all(standard_age=settings_dict.get('standard_age'),
                            template_ages=settings_dict.get('template_ages'),
                            temporal_fwhm=settings_dict.get('temporal_fwhm'),
                            use_mcflirt=settings_dict.get('use_mcflirt'),
                            s2v=settings_dict.get('s2v'),
                            dc=settings_dict.get('dc'),
                            mbs=settings_dict.get('mbs'),
                            quick=settings_dict.get('quick'),
                            icadim=settings_dict.get('icadim'),
                            rdata=settings_dict.get('rdata'),
                            fix_threshold=settings_dict.get('fix_threshold'),
                            group_map=settings_dict.get('group_map'),
                            standard_res=settings_dict.get('standard_res'),
                            group_qc=settings_dict.get('group_qc'),
                            atlasdir=settings_dict.get('atlasdir'),
                            preproc_only=settings_dict.get('preproc_only'),
                            smooth=settings_dict.get('smooth'),
                            intnorm=settings_dict.get('intnorm'))
    else:
        print("\nREQUIRED: '--outdir', '--sub', AND '--run'.\n")
        parser.print_help()

    return 0


def arg_parser() -> Tuple[argparse.ArgumentParser.parse_args, argparse.ArgumentParser]:
    """Argument parser.
    """
    # Init parser
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=__doc__,
                                                              formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=55, width=100)
                                                              #, formatter_class=argparse.RawTextHelpFormatter)
                                                             )

    # Parse Arguments

    mainargparser = parser.add_subparsers(title="subcommands",
                                          description="Data preprocessing pipeline sections for 'fmri_preproc'.",
                                          help="Type the 'subcommand' name followed by '-h' or '--help' for more information.")

    reqoptions = parser.add_argument_group('Required Arguments')

    reqoptions.add_argument('-s','--sub',
                            dest="sub",
                            type=str,
                            metavar="<str>",
                            default=None,
                            help="REQUIRED: Subject ID.")
    reqoptions.add_argument('-r','--run',
                            dest="run",
                            type=str,
                            metavar="<str>",
                            default=None,
                            help="REQUIRED: Run number.")
    reqoptions.add_argument('-o','--outdir',
                            dest="ses",
                            type=str,
                            metavar="<str>",
                            default=None,
                            help="REQUIRED: Output parent directory.")

    optoptions = parser.add_argument_group('Optional Arguments')

    optoptions.add_argument('-e','--ses',
                            dest="ses",
                            type=str,
                            metavar="<str>",
                            default=None,
                            help="OPTIONAL: Session ID.")
    optoptions.add_argument('-version','--version',
                            dest="version",
                            action="store_true",
                            help="Prints version, then exits.")
    optoptions.add_argument('--write-config',
                            type=str,
                            metavar="<file>",
                            dest="out_settings",
                            default=None,
                            help="Writes a template configuration (JSON) to file.")
    
    # Generic input options
    optoptions.add_argument('-v','--verbose',
                            dest="verbose",
                            action="store_true",
                            help="Enables verbose output to the command and log files.")
    optoptions.add_argument('--log-level',
                            type=str,
                            metavar="<str>",
                            dest="log_level",
                            default='info',
                            help="Log level to be written to output log files. Acceptable inputs include: 'info', 'debug', 'critical', 'warning', and 'error' [default='info'].")
    optoptions.add_argument('--config', '--settings',
                            type=str,
                            metavar="<file>",
                            dest="settings",
                            default=DEFAULT_SETTINGS_FILE,
                            help="Configuration (JSON) file that contains additional input arguments not specified at the command line.")

    # IMPORT DATA ARGS
    importoptions = mainargparser.add_parser('import-data', 
                                             parents=[reqoptions],
                                             formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=55, width=100),
                                             add_help=False,
                                             help="Imports the subject's imaging data into a working directory for image preprocessing.")

    importoptions.add_argument('-a','--age',
                               type=str,
                               metavar="<age>",
                               dest="age",
                               default=None,
                               help="REQUIRED: Age at scan. Acceptable inputs include: 'neo', 28 - 44 (PMA, in wks), and 1 or 2 (yrs).")
    importoptions.add_argument('--birth-age',
                               type=float,
                               metavar="<age>",
                               dest="birth_age",
                               default=None,
                               help="Gestational age at birth.")

    # Structural import data
    importoptions.add_argument('--T2w',
                               type=str,
                               metavar="<file>",
                               dest="T2w",
                               default=None,
                               help="REQUIRED: Input (bias corrected, e.g. '_restore') T2w image data.")
    importoptions.add_argument('--T2w-brainmask',
                               type=str,
                               metavar="<file>",
                               dest="REQUIRED: T2w_brainmask",
                               default=None,
                               help="Input T2w brainmask image data (_brainmask_bet).")
    importoptions.add_argument('--T2w-tissue-seg',
                               type=str,
                               metavar="<file>",
                               dest="dseg",
                               default=None,
                               help="REQUIRED: Input Draw-EM or other neonatal brain tissue segmentation (e.g. '_drawem_tissue_labels') image data. NOTE: Other brain tissue files can be used (e.g. from FreeSurfer or FSL's FAST).")
    importoptions.add_argument('--tissue-seg-type',
                               type=str,
                               metavar="<str>",
                               dest="dseg_type",
                               default=None,
                               help="REQUIRED: Input brain tissue segmentation type. Acceptable inputs include: 'drawem', 'freesurfer_aseg', or 'fsl_fast'.")
    importoptions.add_argument('--T1w',
                               type=str,
                               metavar="<file>",
                               dest="T1w",
                               default=None,
                               help="Input (bias corrected, e.g. '_restore') T1w image data.")
    
    # Functional import data
    importoptions.add_argument('-f','--func',
                               type=str,
                               metavar="<file>",
                               dest="func",
                               default=None,
                               help="REQUIRED: Input rs-fMRI data. NOTE: IF input data conforms to the BIDS naming convention (e.g. 'sub-001[_ses-001]_task-rest[_run-01]_bold.nii.gz'), then '--sub', '--ses', and '--run' need not be specified.")
    importoptions.add_argument('--func-echospacing',
                               type=float,
                               metavar="<float>",
                               dest="func_echospacing",
                               default=None,
                               help="REQUIRED: Functional echospacing. NOTE: If unsure use 0.0005.")
    importoptions.add_argument('--func-pedir',
                               type=str,
                               metavar="<PE-dir>",
                               dest="func_pedir",
                               default=None,
                               help="REQUIRED: Phase encoding direction. Acceptable inputs include: 'PA', 'AP', 'RL', 'LR', 'IS', or 'SI'.")
    importoptions.add_argument('--func-slorder',
                               type=str,
                               metavar="<file>",
                               dest="func_slorder",
                               default=None,
                               help="Functional slice order specification file. This file specifies the slice acquisition order of the input rs-fMRI. If not provided one can be created internally IF '--func-mbfactor' is specified.")
    importoptions.add_argument('--func-brainmask',
                               type=str,
                               metavar="<file>",
                               dest="func_brainmask",
                               default=None,
                               help="Functional brain mask.")
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
                               default=None,
                               help="Functional multi-band acceleration factor.")
    
    # Single-band reference import data
    importoptions.add_argument('--sbref',
                               type=str,
                               metavar="<file>",
                               dest="sbref",
                               default=None,
                               help="Single band reference MR image data.")
    importoptions.add_argument('--sbref-brainmask',
                               type=str,
                               metavar="<file>",
                               dest="sbref_brainmask",
                               default=None,
                               help="Single band reference MR brain mask image data.")
    importoptions.add_argument('--sbref-echospacing',
                               type=float,
                               metavar="<float>",
                               dest="sbref_echospacing",
                               default=None,
                               help="Single band reference MR image echospacing. If not specified, then it is assumed to match '--func-echospacing'.")
    importoptions.add_argument('--sbref-pedir',
                               type=str,
                               metavar="<PE-dir>",
                               dest="sbref_pedir",
                               default=None,
                               help="Single band reference MR image data phase encoding direction. NOTE: This should match '--func-pedir'. If not specified, then it is assumed to match '--func-pedir'. If this is not case, this could result in mis-registration and distortion correction issues.")
    
    # Single-band reference import data
    importoptions.add_argument('--spinecho',
                               type=str,
                               metavar="<file>",
                               dest="spinecho",
                               default=None,
                               help="Spinecho MR image data.")
    importoptions.add_argument('--sp-echospacing',
                               type=float,
                               metavar="<float>",
                               dest="spinecho_echospacing",
                               default=None,
                               help="Spinecho MR image echospacing. If not specified, then it is assumed to match '--func-echospacing'.")
    importoptions.add_argument('--sp-pedir',
                               type=str,
                               metavar="<PE-dir,PE-dir,..,PE-dir>",
                               dest="spinecho_pedir",
                               default=None,
                               help="Spinecho MR image phase encoding direction(s), specified as a comma separated list. Accepatable inputs are the same as those specifed for both '--func-pedir', and '--sbref-pedir'.")
    importoptions.add_argument('--sp-inplane-accel',
                               type=float,
                               metavar="<float>",
                               dest="spinecho_inplaneacc",
                               default=None,
                               help="Spinecho MR image inplane acceleration . If not specified, then it is assumed to match '--func-inplane-accel'.")
    importoptions.add_argument('--sp-AP',
                               type=str,
                               metavar="<file>",
                               dest="ap_dir",
                               default=None,
                               help="A -> P (anterior to posterior) Spinecho MR image data.")
    importoptions.add_argument('--sp-PA',
                               type=str,
                               metavar="<file>",
                               dest="pa_dir",
                               default=None,
                               help="P -> A (posterior to anterior) Spinecho MR image data.")
    importoptions.add_argument('--sp-LR',
                               type=str,
                               metavar="<file>",
                               dest="lr_dir",
                               default=None,
                               help="L -> R (left to right) Spinecho MR image data.")
    importoptions.add_argument('--sp-RL',
                               type=str,
                               metavar="<file>",
                               dest="rl_dir",
                               default=None,
                               help="R -> L (right to left) Spinecho MR image data.")
    importoptions.add_argument('--sp-IS',
                               type=str,
                               metavar="<file>",
                               dest="is_dir",
                               default=None,
                               help="I -> S (inferior to superior) Spinecho MR image data.")
    importoptions.add_argument('--sp-SI',
                               type=str,
                               metavar="<file>",
                               dest="si_dir",
                               default=None,
                               help="S -> I (superior to inferior) Spinecho MR image data.")

    importoptions.set_defaults(method='import')

    
    # PRE-MCDC ARGS
    premcdcoptions = mainargparser.add_parser('pre-mcdc', 
                                              parents=[reqoptions],
                                              formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=55, width=100),
                                              add_help=False,
                                              help="Pre-Motion Correction and Distortion Correction (MCDC).")

    premcdcoptions.set_defaults(method='pre_mcdc')

    # MCDC ARGS
    mcdcoptions = mainargparser.add_parser('mcdc', 
                                           parents=[reqoptions],
                                           formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=55, width=100),
                                           add_help=False,
                                           help="Motion Correction and Distortion Correction (MCDC).")
    
    # Slice-to-volume args
    mcdcoptions.add_argument('--s2v', '--slice-to-volume',
                             dest="s2v",
                             action="store_true",
                             default=None,
                             help="Enables slice-to-volume motion correction. NOTE: This will LIKELY require EDDY_CUDA and GPU processing. This is possible without GPUs in FSL v6.0.5+ but is VERY SLOW.")
    mcdcoptions.add_argument('--dc', '--distortion-correction',
                             dest="dc",
                             action="store_true",
                             default=None,
                             help="Enables distortion correction.")
    mcdcoptions.add_argument('--mbs',
                             dest="mbs",
                             action="store_true",
                             default=None,
                             help="Enables movement by susceptibility distortion correction.")
    mcdcoptions.add_argument('--mcflirt',
                             dest="use_mcflirt",
                             action="store_true",
                             default=None,
                             help="Enables FSL's MCFLIRT-based motion correction. NOTE: Enabling this option DISABLES: '--s2v', '--dc', and '--mbs' options.")

    mcdcoptions.set_defaults(method='mcdc')

    # POST-MCDC
    postmcdcoptions = mainargparser.add_parser('post-mcdc', 
                                               parents=[reqoptions],
                                               formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=55, width=100),
                                               add_help=False,
                                               help="Post-MCDC preprocessing, which includes: func -> template registration, ICA-based denoising, and QC (web page).")

    # Registration args
    postmcdcoptions.add_argument('--standard-age',
                                 type=int,
                                 metavar="<int>",
                                 dest="standard_age",
                                 default=40,
                                 help="Standard age to use for final dHCP age-matched template (PMA, in weeks) [default=40 (wks)].")
    postmcdcoptions.add_argument('--template-ages',
                                 type=str,
                                 metavar="<str,str,..str>",
                                 dest="template_ages",
                                 default=None,
                                 help="Comma separated list of ages. Acceptable inputs include: 'neo', 1, 2 (yrs), and/or 28 - 44 (PMA, wks).")
    postmcdcoptions.add_argument('--atlas-dir',
                                 type=str,
                                 metavar="<path>",
                                 dest="atlasdir",
                                 default=None,
                                 help="Atlas directory path.")

    # ICA Args
    postmcdcoptions.add_argument('--temporal-fwhm',
                                 type=float,
                                 metavar="<float>",
                                 dest="temporal_fwhm",
                                 default=150.0,
                                 help="Temporal FWHM filtering coefficient [default=150.0 (s)].")
    postmcdcoptions.add_argument('--ica-dim',
                                 type=int,
                                 metavar="<int>",
                                 dest="icadim",
                                 default=None,
                                 help="Maximal ICA dimensionality.")

    # FIX Args
    postmcdcoptions.add_argument('--rdata',
                                 type=str,
                                 metavar="<str>",
                                 dest="rdata",
                                 default=None,
                                 help="Trained FIX classifier.")
    postmcdcoptions.add_argument('--fix-threshold',
                                 type=int,
                                 metavar="<int>",
                                 dest="fix_threshold",
                                 default=None,
                                 help="FIX noise classification threshold. [default=10].") # Default is in settings function.

    # POST-PROCESS
    postmcdcoptions.add_argument('--smooth',
                                 type=float,
                                 metavar="<float>",
                                 dest="smooth",
                                 default=None,
                                 help="Smoothing kernel (FWHM, in mm) [default: disabled].")
    postmcdcoptions.add_argument('--int-nrom',
                                 dest="intnorm",
                                 action="store_true",
                                 help="Performs intensity normalization OR grand-mean scaling IF '--smooth' is specified and greater than 0 [default: disabled].")

    postmcdcoptions.set_defaults(method='post_mcdc')

    # RUN ALL SECTIONS
    runalloptions = mainargparser.add_parser('run-all', 
                                             parents=[importoptions],
                                             formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=55, width=100),
                                             add_help=False,
                                             help="Perform/Run all preprocessing steps of the pipeline.")

    # Slice-to-volume args
    runalloptions.add_argument('--s2v', '--slice-to-volume',
                             dest="s2v",
                             action="store_true",
                             default=None,
                             help="Enables slice-to-volume motion correction. NOTE: This will LIKELY require EDDY_CUDA and GPU processing. This is possible without GPUs in FSL v6.0.5+ but is VERY SLOW.")
    runalloptions.add_argument('--dc', '--distortion-correction',
                             dest="dc",
                             action="store_true",
                             default=None,
                             help="Enables distortion correction.")
    runalloptions.add_argument('--mbs',
                             dest="mbs",
                             action="store_true",
                             default=None,
                             help="Enables movement by susceptibility distortion correction.")
    runalloptions.add_argument('--mcflirt',
                             dest="use_mcflirt",
                             action="store_true",
                             default=None,
                             help="Enables FSL's MCFLIRT-based motion correction. NOTE: Enabling this option DISABLES: '--s2v', '--dc', and '--mbs' options.")

    # Registration args
    runalloptions.add_argument('--standard-age',
                                 type=int,
                                 metavar="<int>",
                                 dest="standard_age",
                                 default=40,
                                 help="Standard age to use for final dHCP age-matched template (PMA, in weeks) [default=40 (wks)].")
    runalloptions.add_argument('--template-ages',
                                 type=str,
                                 metavar="<str,str,..str>",
                                 dest="template_ages",
                                 default=None,
                                 help="Comma separated list of ages. Acceptable inputs include: 'neo', 1, 2 (yrs), and/or 28 - 44 (PMA, wks).")
    runalloptions.add_argument('--atlas-dir',
                                 type=str,
                                 metavar="<path>",
                                 dest="atlasdir",
                                 default=None,
                                 help="Atlas directory path.")

    # ICA Args
    runalloptions.add_argument('--temporal-fwhm',
                                 type=float,
                                 metavar="<float>",
                                 dest="temporal_fwhm",
                                 default=150.0,
                                 help="Temporal FWHM filtering coefficient [default=150.0 (s)].")
    runalloptions.add_argument('--ica-dim',
                                 type=int,
                                 metavar="<int>",
                                 dest="icadim",
                                 default=None,
                                 help="Maximal ICA dimensionality.")

    # FIX Args
    runalloptions.add_argument('--rdata',
                                 type=str,
                                 metavar="<str>",
                                 dest="rdata",
                                 default=None,
                                 help="Trained FIX classifier (stored as a Rdata file).")
    runalloptions.add_argument('--fix-threshold',
                                 type=int,
                                 metavar="<int>",
                                 dest="fix_threshold",
                                 default=None,
                                 help="FIX noise classification threshold. [default=10].") # Default is in settings function.
    
    # POST-PROCESS
    runalloptions.add_argument('--smooth',
                               type=float,
                               metavar="<float>",
                               dest="smooth",
                               default=None,
                               help="Smoothing kernel (FWHM, in mm) [default: disabled].")
    runalloptions.add_argument('--int-nrom',
                               dest="intnorm",
                               action="store_true",
                               help="Performs intensity normalization OR grand-mean scaling IF '--smooth' is specified and greater than 0 [default: disabled].")

    runalloptions.set_defaults(method='run_all')

    args: argparse.ArgumentParser.parse_args = parser.parse_args()
    return args, parser


if __name__ == '__main__':
    main()
