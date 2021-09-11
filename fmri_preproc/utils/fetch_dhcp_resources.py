#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Fetch dHCP volumetric atlas, group maps, and/or group QC files from the dHCP study.

NOTE: This does not include the dHCP FIX training data as it is assumed that the data being preprocessed will not 
    likely match the acquisition parameters used for the dHCP study.
"""
import os
import argparse

from typing import Optional, Tuple

from fmri_preproc.utils.util import (
    fetch_dhcp_group_maps,
    fetch_dhcp_group_qc,
    fetch_dhcp_volumetric_atlas
)


def main() -> int:
    """Main function.
    """
    args, parser = arg_parser()

    # Print help message in the case of no arguments
    try:
        args = parser.parse_args()
    except SystemExit as err:
        if err.code == 2:
            parser.print_help()
    
    download_data(path=args.path,
                  fetch_vol_atlas=args.fetch_vol_atlas,
                  fetch_group_maps=args.fetch_group_maps,
                  fetch_group_qc=args.fetch_group_qc)
    return 0


def download_data(path: Optional[str] = None,
                  fetch_vol_atlas: bool = True,
                  fetch_group_maps: bool = True,
                  fetch_group_qc: bool = True
                 ) -> None:
    """Downloads the specified dHCP related data.
    """
    if path is not None:
        atlas_path: str = os.path.join(path,'atlases')
        group_maps_path: str = os.path.join(path,'group_maps')
        group_qc_path: str = os.path.join(path,'group_qc')
    else:
        atlas_path: str = None
        group_maps_path: str = None
        group_qc_path: str = None

    if fetch_vol_atlas:
        fetch_dhcp_volumetric_atlas(path=atlas_path, extended=True)
    
    if fetch_group_maps:
        fetch_dhcp_group_maps(path=group_maps_path)
    
    if fetch_group_qc:
        fetch_dhcp_group_qc(path=group_qc_path)

    return None


def arg_parser() -> Tuple[argparse.ArgumentParser.parse_args, argparse.ArgumentParser]:
    """Argument parser.
    """
    # Init parser
    parser: argparse.ArgumentParser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

    # Parse Arguments

    # Required Arguments
    reqoptions = parser.add_argument_group('Argument(s)')
    reqoptions.add_argument('--path',
                            type=str,
                            metavar="<path>",
                            dest="path",
                            help="Location for resources (default=<FMRI_PREPROC_DIR>/resources).")

    reqoptions.add_argument('--skip-volumetric-atlas',
                            action='store_false',
                            dest="fetch_vol_atlas")

    reqoptions.add_argument('--skip-group-maps',
                            action='store_false',
                            dest="fetch_group_maps")

    reqoptions.add_argument('--skip-group-qc',
                            action='store_false',
                            dest="fetch_group_qc")
    
    args: argparse.ArgumentParser.parse_args = parser.parse_args()

    return args, parser


if __name__ == '__main__':
    main()