#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Performs preprocessing of neonatal resting-state functional MR neuroimages.
"""
# TODO: 
#   * Write CLI
#   * Write function for multi-atlas registration
#   * Create settings function(s)/module and corresponding settings(.json) file.

import os
import argparse

from typing import (
    Optional,
    Tuple
)

from fmri_preproc.func.pipeline import Pipeline


def main():
    """Main function.
    """
    pass


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
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

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
                            metavar="<path>",
                            dest="func",
                            help="REQUIRED: Path to input rs-fMRI data. NOTE: Input data MUST conform to the BIDS naming convention (e.g. 'sub-001[_ses-001]_task-rest[_run-01]_bold.nii.gz').")
    reqoptions.add_argument('-a','--age',
                            type=str,
                            metavar="<age>",
                            dest="age",
                            help="REQUIRED: Age at scan. Acceptable inputs include: 'neo', 28 - 44 (wks), and 1 or 2 (yrs).")

    importoptions = parser.add_argument_group('Import Data Arguments')
    importoptions.add_argument('--T2w',
                               type=str,
                               metavar="<path>",
                               dest="T2w",
                               help="Path to input (bias corrected, e.g. '_restore') T2w image data.")
    importoptions.add_argument('--T2w-brainmask',
                               type=str,
                               metavar="<path>",
                               dest="T2w_brainmask",
                               help="Path to input T2w brainmask image data (_brainmask_bet).")
    importoptions.add_argument('--tissue-seg',
                               type=str,
                               metavar="<path>",
                               dest="dseg",
                               help="Path to input Draw-EM neonatal tissue segmentation (e.g. '_drawem_tissue_labels') image data.")
    pass


if __name__ == '__main__':
    main()
