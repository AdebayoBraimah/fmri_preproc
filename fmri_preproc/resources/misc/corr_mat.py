#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Construct adjacency/correlation matrices for input mean timeseries file of atlas ROIs.
"""
import os
import numpy as np
import sys
import argparse

from typing import (
    Tuple,
    TypeVar
)


file = TypeVar("file", bound=str)


def main() -> None:
    """Main function.
    """
    args, parser = arg_parser()

    # Print help message in the case of no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.mat_file and args.out_file:
        compute_corr_mat(mat_file=args.mat_file,
                         out_file=args.out_file)
    else:
        print("\nREQUIRED: '--matrix-file' and '--output'.\n")
    return None


def compute_corr_mat(mat_file: file,
                     out_file: file
                    ) -> file:
    """Computes correlation matrix for some input mean timeseries of ROIs 
    (represented as a r x c matrix, stored in a text file).
    """
    mat_file: str = os.path.abspath(mat_file)
    mat: np.array = np.loadtxt(mat_file)

    # NOTE: The transpose of the matrix is used
    #   as corrcoef usually returns a matrix of the
    #   dimension r x r, rather than c x c - of which
    #   is expected to convey ROI-to-ROI correlations.
    corr: np.array = np.corrcoef(mat.transpose())
    corr: np.array = np.nan_to_num(corr,nan=0.0)

    np.savetxt(out_file,corr, fmt='%.9f')
    out_file: str = os.path.abspath(out_file)

    return out_file


def arg_parser() -> Tuple[argparse.ArgumentParser.parse_args, argparse.ArgumentParser]:
    """Argument parser.
    """
    # Init parser
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=__doc__,
                                                              formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=55, width=100)
                                                             )
    
    reqoptions = parser.add_argument_group('Required Arguments')

    reqoptions.add_argument('-m','--matrix-file',
                            dest="mat_file",
                            type=str,
                            metavar="<file>",
                            default=None,
                            help="Matrix file of mean timeseries for each ROI")
    reqoptions.add_argument('-o','--output',
                            dest="out_file",
                            type=str,
                            metavar="<file>",
                            default=None,
                            help="Output file name of correlation matrix.")
    args: argparse.ArgumentParser.parse_args = parser.parse_args()
    return args, parser


if __name__ == '__main__':
    main()
