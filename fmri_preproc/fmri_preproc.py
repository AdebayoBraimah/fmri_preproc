#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Performs preprocessing of neonatal resting-state functional MR neuroimages.
"""
# TODO: Recompress UNC.tar.gz on a linux computer
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
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=__doc__)

    # Parse Arguments

    # Required Arguments
    reqoptions = parser.add_argument_group('Required Arguments')
    reqoptions.add_argument('--path',
                            type=str,
                            metavar="<path>",
                            dest="path",
                            help="Location for resources (default = <FMRI_PREPROC_DIR>/resources).")
    pass


if __name__ == '__main__':
    main()
