#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Performs preprocessing of neonatal resting-state functional MR neuroimages.
"""
import sys
from typing import Any, Dict, Literal

from fmri_preproc import DEFAULT_SETTINGS_FILE, __version__ as version_id
from fmri_preproc.func.pipeline import Pipeline
from fmri_preproc.utils.parser import arg_parser
from fmri_preproc.utils.util import dict2json, settings


def main() -> Literal[0]:
    """Main function.
    """
    preproc_data()
    return 0


def preproc_data() -> Literal[0]:
    """Preprocess rs-fMRI data.
    """
    parser = arg_parser()
    args = parser.parse_args()

    # Print help message in the case of no arguments
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        args: Dict[str, Any] = vars(args)

    # Version
    if args.get("version"):
        print(f"\nfmri_preproc: v{version_id}\n")
        sys.exit(0)

    # Write settings file
    if args.get("out_settings"):
        settings_dict: Dict[str, Any] = settings(jsonfile=DEFAULT_SETTINGS_FILE)
        _: str = dict2json(
            dict=settings_dict, jsonfile=args.get("out_settings"), indent=4
        )
        sys.exit(0)

    # Preprocess data

    # Check arguments
    if args.get("outdir") and (
        (args.get("sub") and args.get("run")) or args.get("func")
    ):

        # Confirm phase-encoding directions
        if not args.get("sbref_pedir"):
            args["sbref_pedir"] = args.get("func_pedir")

        # Confirm echo-spacing paramters
        if not args.get("sbref_echospacing"):
            args["sbref_echospacing"] = args.get("func_echospacing")

        if not args.get("spinecho_echospacing"):
            args["spinecho_echospacing"] = args.get("func_echospacing")

        if not args.get("spinecho_inplaneacc"):
            args["spinecho_inplaneacc"] = args.get("func_inplane_accel")

        # Configure runtime settings
        settings_dict: Dict[str, Any] = settings(
            jsonfile=args.get("settings"),
            func_inplane_accel=args.get("func_inplane_accel"),
            func_slorder=args.get("func_slorder"),
            mb_factor=args.get("mb_factor"),
            sbref_echospacing=args.get("sbref_echospacing"),
            dseg_type=args.get("dseg_type"),
            probseg_type=None,  # Set this to None for now
            mask_func=None,  # Set this to None for now
            spinecho_echospacing=args.get("spinecho_echospacing"),
            spinecho_epifactor=None,  # Set this to None for now
            spinecho_inplaneacc=args.get("spinecho_inplaneacc"),
            use_mcflirt=args.get("use_mcflirt"),
            s2v=args.get("s2v"),
            dc=args.get("dc"),
            mbs=args.get("mbs"),
            standard_age=args.get("standard_age"),
            quick=None,
            atlasdir=args.get("atlasdir"),
            template_ages=args.get("template_ages"),
            temporal_fwhm=args.get("temporal_fwhm"),
            icadim=args.get("icadim"),
            rdata=args.get("rdata"),
            fix_threshold=args.get("fix_threshold"),
            group_qc=None,
            group_map=None,
            standard_res=None,
            verbose=args.get("verbose"),
            log_level=args.get("log_level"),
            smooth=args.get("smooth"),
            intnorm=args.get("intnorm"),
        )

        # Get subject info
        preproc: Pipeline = Pipeline(
            outdir=args.get("outdir"),
            subid=args.get("sub"),
            sesid=args.get("ses"),
            runid=args.get("run"),
            func=args.get("func"),
            scan_pma=args.get("age"),
            birth_ga=args.get("birth_age"),
            verbose=settings_dict.get("verbose"),
            log_level=settings_dict.get("log_level"),
            settings_json_dict=settings_dict,
        )
        if args.get("method") == "import":
            preproc.import_data(
                func=args.get("func"),
                func_echospacing=args.get("func_echospacing"),
                func_pedir=args.get("func_pedir"),
                T2w=args.get("T2w"),
                T2w_brainmask=args.get("T2w_brainmask"),
                dseg=args.get("dseg"),
                func_brainmask=args.get("func_brainmask"),
                func_slorder=settings_dict.get("func_slorder"),
                func_inplane_accel=settings_dict.get("func_inplane_accel"),
                mb_factor=settings_dict.get("mb_factor"),
                sbref=args.get("sbref"),
                sbref_brainmask=args.get("sbref_brainmask"),
                sbref_echospacing=settings_dict.get("sbref_echospacing"),
                sbref_pedir=args.get("sbref_pedir"),
                dseg_type=settings_dict.get("dseg_type"),
                T1w=args.get("T1w"),
                wmmask=args.get("wmmask"),
                spinecho=args.get("spinecho"),
                spinecho_echospacing=settings_dict.get("spinecho_echospacing"),
                spinecho_pedir=args.get("spinecho_pedir"),
                ap_dir=args.get("ap_dir"),
                pa_dir=args.get("pa_dir"),
                lr_dir=args.get("lr_dir"),
                rl_dir=args.get("rl_dir"),
                is_dir=args.get("is_dir"),
                si_dir=args.get("si_dir"),
                spinecho_epifactor=settings_dict.get("spinecho_epifactor"),
                spinecho_inplaneacc=settings_dict.get("spinecho_inplaneacc"),
            )
        elif args.get("method") == "pre_mcdc":
            preproc.pre_mcdc()
        elif args.get("method") == "mcdc":
            preproc.mcdc(
                use_mcflirt=settings_dict.get("use_mcflirt"),
                s2v=settings_dict.get("s2v"),
                dc=settings_dict.get("dc"),
                mbs=settings_dict.get("mbs"),
                mporder=args.get("mporder"),
            )
        elif args.get("method") == "post_mcdc":
            preproc.post_mcdc(
                standard_age=settings_dict.get("standard_age"),
                template_ages=settings_dict.get("template_ages"),
                temporal_fwhm=settings_dict.get("temporal_fwhm"),
                quick=settings_dict.get("quick"),
                icadim=settings_dict.get("icadim"),
                rdata=settings_dict.get("rdata"),
                fix_threshold=settings_dict.get("fix_threshold"),
                group_map=settings_dict.get("group_map"),
                standard_res=settings_dict.get("standard_res"),
                group_qc=settings_dict.get("group_qc"),
                atlasdir=settings_dict.get("atlasdir"),
                preproc_only=settings_dict.get("preproc_only"),
                spatial_fwhm=settings_dict.get("smooth"),
                intnorm=settings_dict.get("intnorm"),
            )
        elif args.get("method") == "run_all":
            # Import data
            preproc.import_data(
                func=args.get("func"),
                func_echospacing=args.get("func_echospacing"),
                func_pedir=args.get("func_pedir"),
                T2w=args.get("T2w"),
                T2w_brainmask=args.get("T2w_brainmask"),
                dseg=args.get("dseg"),
                func_brainmask=args.get("func_brainmask"),
                func_slorder=settings_dict.get("func_slorder"),
                func_inplane_accel=settings_dict.get("func_inplane_accel"),
                mb_factor=settings_dict.get("mb_factor"),
                sbref=args.get("sbref"),
                sbref_brainmask=args.get("sbref_brainmask"),
                sbref_echospacing=settings_dict.get("sbref_echospacing"),
                sbref_pedir=args.get("sbref_pedir"),
                dseg_type=settings_dict.get("dseg_type"),
                T1w=args.get("T1w"),
                wmmask=args.get("wmmask"),
                spinecho=args.get("spinecho"),
                spinecho_echospacing=settings_dict.get("spinecho_echospacing"),
                spinecho_pedir=args.get("spinecho_pedir"),
                ap_dir=args.get("ap_dir"),
                pa_dir=args.get("pa_dir"),
                lr_dir=args.get("lr_dir"),
                rl_dir=args.get("rl_dir"),
                is_dir=args.get("is_dir"),
                si_dir=args.get("si_dir"),
                spinecho_epifactor=settings_dict.get("spinecho_epifactor"),
                spinecho_inplaneacc=settings_dict.get("spinecho_inplaneacc"),
            )

            # RUN-ALL
            preproc.run_all(
                standard_age=settings_dict.get("standard_age"),
                template_ages=settings_dict.get("template_ages"),
                temporal_fwhm=settings_dict.get("temporal_fwhm"),
                use_mcflirt=settings_dict.get("use_mcflirt"),
                s2v=settings_dict.get("s2v"),
                dc=settings_dict.get("dc"),
                mbs=settings_dict.get("mbs"),
                quick=settings_dict.get("quick"),
                icadim=settings_dict.get("icadim"),
                rdata=settings_dict.get("rdata"),
                fix_threshold=settings_dict.get("fix_threshold"),
                group_map=settings_dict.get("group_map"),
                standard_res=settings_dict.get("standard_res"),
                group_qc=settings_dict.get("group_qc"),
                atlasdir=settings_dict.get("atlasdir"),
                preproc_only=settings_dict.get("preproc_only"),
                spatial_fwhm=settings_dict.get("smooth"),
                intnorm=settings_dict.get("intnorm"),
            )
    else:
        print("\nREQUIRED: '--outdir', '--sub', AND '--run'.\n")
        parser.print_help()

    return 0


if __name__ == "__main__":
    main()
