# -*- coding: utf-8 -*-
"""fMRI preprocessing pipeline - main pipeline.

.. autosummary::
    :nosignatures:

    Pipeline

.. autoclass:: Pipeline
    :members:
"""
import os
import numpy as np
import nibabel as nib
from collections import OrderedDict

from typing import Any, Dict, List, Optional, Set, Tuple, Union

from fmri_preproc import DEFAULT_CLASSIFIER, GROUP_MAP_DIR, GROUP_QC_DIR

from fmri_preproc.utils.fileio import NiiFile
from fmri_preproc.func.qc import Subject
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.func.fieldmap import fieldmap
from fmri_preproc.func.mcdc import mcdc
from fmri_preproc.func.ica import ica

from fmri_preproc.utils.fslpy import FSLDIR, flirt

from fmri_preproc.utils.util import (
    dict2json,
    json2dict,
    load_sidecar,
    settings,
    update_sidecar,
)

from fmri_preproc.utils.outputs.importdata import (
    ImportFunc,
    ImportSpinEcho,
    ImportStruct,
)

from fmri_preproc.utils.outputs.outdict import key_to_none
from fmri_preproc.utils.outputs.fieldmap import FmapFiles
from fmri_preproc.utils.outputs.registration import MRreg
from fmri_preproc.utils.outputs.mcdc import MCDCFiles
from fmri_preproc.utils.outputs.ica import ICAFiles
from fmri_preproc.func.postprocess import postprocess

from fmri_preproc.utils.outputs.denoise import FIXApply, FIXClassify, FIXExtract

from fmri_preproc.func.importdata import (
    import_func,
    import_info,
    import_spinecho,
    import_struct,
)

from fmri_preproc.func.registration import (
    fmap_to_struct,
    fmap_to_func_composite,
    func_to_sbref,
    func_to_struct_composite,
    func_to_template_composite,
    sbref_to_struct,
    _select_atlas as select_atlas,
    struct_to_template_composite,
    template_to_struct,
)

from fmri_preproc.func.denoise import fix_apply, fix_classify, fix_extract

from fmri_preproc.utils.mask import create_mask, get_dseg_labels


class Pipeline:
    """class doc-string
    """

    def __init__(
        self,
        outdir: str,
        func: Optional[str] = None,
        subid: Optional[str] = None,
        sesid: Optional[str] = None,
        runid: Optional[str] = None,
        scan_pma: Optional[float] = None,
        birth_ga: Optional[float] = None,
        verbose: bool = False,
        log_level: str = "info",
        settings_json_dict: Optional[Union[Dict[str, Any], str]] = None,
        **kwargs,
    ) -> None:
        """Constructor for pipeline class.
        """
        # Import information
        self.verbose: bool = verbose

        if isinstance(settings_json_dict, str):
            self.settings: Dict[str, Any] = settings(jsonfile=settings_json_dict)
        elif isinstance(settings_json_dict, dict):
            self.settings: Dict[str, Any] = settings(**settings_json_dict)
        else:
            self.settings: Dict[str, Any] = settings()

        (sub_workdir, logdir, sub_json, sub_dict) = import_info(
            outdir=outdir,
            func=func,
            subid=subid,
            sesid=sesid,
            runid=runid,
            scan_pma=scan_pma,
            birth_ga=birth_ga,
            log=None,
            verbose=verbose,
            log_datetime=True,
            log_level=log_level,
            **kwargs,
        )

        # Define additional class variables
        self.logdir: str = logdir
        self.verbose: bool = verbose
        self.workdir: str = sub_workdir
        self.sub_dict: str = sub_dict
        self.func: str = func
        self.sub_json: str = sub_json

        with WorkDir(src=self.logdir) as lgd:
            _import_log: str = lgd.join("import.log")
            import_log: LogFile = LogFile(
                _import_log, format_log_str=True, print_to_screen=self.verbose
            )
            self.import_log: LogFile = import_log

            subid: str = sub_dict.get("subid")
            sesid: str = sub_dict.get("sesid")

            if sesid:
                self.proc: str = lgd.join(f"sub-{subid}_ses-{sesid}_fmripreproc.json")
            else:
                self.proc: str = lgd.join(f"sub-{subid}_fmripreproc.json")

            if os.path.exists(self.proc):
                self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)
                self.outputs: Dict[str, str] = key_to_none(d=self.outputs)
            else:
                self.outputs: Dict[str, str] = {}
        return None

    def import_data(
        self,
        func: Optional[str] = None,
        func_echospacing: Optional[float] = None,
        func_pedir: Optional[str] = None,
        T2w: Optional[str] = None,
        T2w_brainmask: Optional[str] = None,
        dseg: Optional[str] = None,
        func_brainmask: Optional[str] = None,
        func_slorder: Optional[str] = None,
        func_inplane_accel: Optional[float] = 1,
        mb_factor: Optional[int] = None,
        sbref: Optional[str] = None,
        sbref_brainmask: Optional[str] = None,
        sbref_echospacing: Optional[float] = None,
        sbref_pedir: Optional[str] = None,
        dseg_type: Optional[str] = "drawem",
        probseg: Optional[str] = None,
        probseg_type: Optional[str] = "drawem",
        wmmask: Optional[str] = None,
        T1w: Optional[str] = None,
        mask_func: bool = False,
        spinecho: Optional[str] = None,
        spinecho_echospacing: Optional[float] = None,
        spinecho_pedir: Union[str, List[str]] = None,
        ap_dir: Optional[str] = None,
        pa_dir: Optional[str] = None,
        lr_dir: Optional[str] = None,
        rl_dir: Optional[str] = None,
        is_dir: Optional[str] = None,
        si_dir: Optional[str] = None,
        spinecho_epifactor: Optional[int] = None,
        spinecho_inplaneacc: Optional[float] = None,
    ) -> Dict[Any, str]:
        """Import data class-function.
        """
        # Import information
        sub_workdir: str = self.workdir
        import_log: str = self.import_log
        sub_dict: Dict[Any, Any] = self.sub_dict
        # func: str = self.func
        sub_json: str = self.sub_json
        settings_file: str = os.path.join(sub_workdir, "logs", "settings.json")

        # Check if func has been passed as class arg
        if self.func:
            func: str = self.func

        # Use previuos settings if present
        if os.path.exists(settings_file):
            self.settings: Dict[str, Any] = json2dict(jsonfile=settings_file)

        settings_file: str = dict2json(dict=self.settings, jsonfile=settings_file)

        # Check if func has been imported
        out_func: ImportFunc = ImportFunc(outdir=sub_workdir)
        _: Dict[str, str] = out_func.outputs()
        func_files: Tuple[str] = ("func", "func_mean", "func_brainmask", "func_slorder")

        if not out_func.check_exists(*func_files):
            if func and func_echospacing and func_pedir:
                (
                    func,
                    _,
                    func_brainmask,
                    func_slorder,
                    sbref,
                    sbref_brainmask,
                    func_info_dict,
                ) = import_func(
                    outdir=sub_workdir,
                    func=func,
                    func_echospacing=func_echospacing,
                    func_pedir=func_pedir,
                    func_brainmask=func_brainmask,
                    func_slorder=func_slorder,
                    func_inplane_accel=func_inplane_accel,
                    mb_factor=mb_factor,
                    sbref=sbref,
                    sbref_brainmask=sbref_brainmask,
                    sbref_echospacing=sbref_echospacing,
                    sbref_pedir=sbref_pedir,
                    mask_func=mask_func,
                    log=import_log,
                )
            else:
                raise RuntimeError(
                    "Minimum required information was not provided: 'func_echospacing' and/or 'func_pedir'."
                )
        else:
            func_info_dict: Dict[str, str] = out_func.outputs()

        # Check if struct has been imported
        out_struct: ImportStruct = ImportStruct(outdir=sub_workdir)
        _: Dict[str, str] = out_struct.outputs()
        struct_files: Tuple[str] = ("T2w", "T2w_brainmask", "T2w_dseg")

        if not out_struct.check_exists(*struct_files):
            if T2w and T2w_brainmask and dseg and dseg_type:
                (T2w, wmmask, dseg, struct_info_dict) = import_struct(
                    outdir=sub_workdir,
                    T2w=T2w,
                    brainmask=T2w_brainmask,
                    dseg=dseg,
                    dseg_type=dseg_type,
                    probseg=probseg,
                    probseg_type=probseg_type,
                    wmmask=wmmask,
                    T1w=T1w,
                    log=import_log,
                )
            else:
                raise RuntimeError(
                    "Required inputs were not found and/or provided: 'T2w', 'T2w_brainmask', 'dseg', 'dseg_type'."
                )
        else:
            struct_info_dict: Dict[str, str] = out_struct.outputs()

        # Check if fmap has been imported
        out_fmap: ImportSpinEcho = ImportSpinEcho(outdir=sub_workdir)
        _: Dict[str, str] = out_fmap.outputs()
        fmap_files: Tuple[str] = ("spinecho",)

        if not out_fmap.check_exists(*fmap_files):
            if (
                spinecho
                or (ap_dir and pa_dir)
                or (lr_dir and rl_dir)
                or (is_dir and si_dir)
            ):
                (spinecho, spinecho_pedir) = import_spinecho(
                    outdir=sub_workdir,
                    spinecho=spinecho,
                    spinecho_echospacing=spinecho_echospacing,
                    spinecho_pedir=spinecho_pedir,
                    ap_dir=ap_dir,
                    pa_dir=pa_dir,
                    lr_dir=lr_dir,
                    rl_dir=rl_dir,
                    is_dir=is_dir,
                    si_dir=si_dir,
                    spinecho_epifactor=spinecho_epifactor,
                    spinecho_inplaneacc=spinecho_inplaneacc,
                    log=import_log,
                )
            else:
                spinecho: str = None
                spinecho_pedir: str = None
        else:
            fmap_info_dict: Dict[str, str] = out_fmap.outputs()
            spinecho: str = fmap_info_dict.get("spinecho")

        self.outputs: Dict[str, str] = {
            **sub_dict,
            **func_info_dict,
            **struct_info_dict,
            "settings_file": settings_file,
            "subject_info": sub_json,
            "workdir": sub_workdir,
            "spinecho": spinecho,
            # "subject_info": sub_json.name,
        }

        # Set files that do not exist to None
        self.outputs: Dict[str, str] = key_to_none(d=self.outputs)
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # Remove unneeded class variables
        attrs: Tuple[str] = ("func", "workdir", "import_log", "sub_dict", "sub_json")

        for attr in attrs:
            self.__dict__.pop(attr, None)

        return self.outputs

    def prepare_fieldmap(self) -> Dict[Any, str]:
        """Prepare fieldmaps for additional preprocessing steps.
        """
        with WorkDir(src=self.logdir) as lgd:
            _fmap_log: str = lgd.join("fmap.log")
            fmap_log: LogFile = LogFile(_fmap_log, format_log_str=True)

        self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)

        spinecho: str = self.outputs.get("spinecho")
        metadata: Dict[str, str] = load_sidecar(file=spinecho)

        fmap_out: FmapFiles = FmapFiles(outdir=self.outputs.get("workdir"))
        fmap_info_dict: Dict[str, str] = fmap_out.outputs()
        fmap_files: Tuple[str] = ("fmap", "fmap_mag", "fmap_mask")

        if not fmap_out.check_exists(*fmap_files):
            _: Tuple[str] = fieldmap(
                outdir=self.outputs.get("workdir"),
                spinecho=spinecho,
                echo_spacing=metadata.get("echo_spacing"),
                pedir=metadata.get("phase_encode_dir"),
                epifactor=metadata.get("epi_factor"),
                inplane_acc=metadata.get("inplane_accel", 1.0),
                verbose=self.verbose,
                log=fmap_log,
            )

        self.outputs: Dict[str, str] = {
            **self.outputs,
            **fmap_info_dict,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # fmap -> struct

        src_space: str = "fmap"
        ref_space: str = "struct"

        fm2st_reg: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        fm2st_dict: Dict[str, str] = fm2st_reg.outputs(
            src_space=src_space, ref_space=ref_space
        )
        fm2st_files: Tuple[str] = (
            "affine",
            "inv_affine",
            "resampled_image",
            "init_affine",
        )

        if not fm2st_reg.check_exists(*fm2st_files):
            (
                fmap2struct_affine,
                fmap2struct_inv_affine,
                fmap2struct_src2ref,
                fmap2struct_init_affine,
            ) = fmap_to_struct(
                outdir=self.outputs.get("workdir"),
                fmap=self.outputs.get("fmap"),
                fmap_magnitude=self.outputs.get("fmap_mag"),
                fmap_brainmask=self.outputs.get("fmap_mask"),
                struct=self.outputs.get("T2w"),
                struct_brainmask=self.outputs.get("T2w_brainmask"),
                struct_boundarymask=self.outputs.get("T2w_wmmask"),
                src_space=src_space,
                ref_space=ref_space,
                bbr=True,
                log=fmap_log,
            )
        else:
            fmap2struct_affine: str = fm2st_dict.get("affine")
            fmap2struct_inv_affine: str = fm2st_dict.get("inv_affine")
            fmap2struct_src2ref: str = fm2st_dict.get("resampled_image")
            fmap2struct_init_affine: str = fm2st_dict.get("init_affine")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "fmap2struct_affine": fmap2struct_affine,
            "fmap2struct_inv_affine": fmap2struct_inv_affine,
            "fmap2struct_src2ref": fmap2struct_src2ref,
            "fmap2struct_init_affine": fmap2struct_init_affine,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # func -> sbref (distorted)

        src_space: str = "func"
        ref_space: str = "sbref"

        fn2sb_reg: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        fn2sb_dict: Dict[str, str] = fn2sb_reg.outputs(
            src_space=src_space, ref_space=ref_space
        )
        fn2sb_files: Tuple[str] = ("resampled_image", "affine", "inv_affine")

        if not fn2sb_reg.check_exists(*fn2sb_files):
            (
                func2sbref_src2ref,
                func2sbref_affine,
                func2sbref_inv_affine,
            ) = func_to_sbref(
                outdir=self.outputs.get("workdir"),
                func=self.outputs.get("func0"),
                func_brainmask=self.outputs.get("func_brainmask"),
                sbref=self.outputs.get("sbref"),
                sbref_brainmask=self.outputs.get("sbref_brainmask"),
                src_space=src_space,
                ref_space=ref_space,
                log=fmap_log,
            )
        else:
            func2sbref_src2ref: str = fn2sb_dict.get("resampled_image")
            func2sbref_affine: str = fn2sb_dict.get("affine")
            func2sbref_inv_affine: str = fn2sb_dict.get("inv_affine")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "func2sbref_src2ref": func2sbref_src2ref,
            "func2sbref_affine": func2sbref_affine,
            "func2sbref_inv_affine": func2sbref_inv_affine,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # sbref -> struct (with BBR and DC)

        sbref: str = self.outputs.get("sbref")
        sbref_meta: Dict[str, str] = load_sidecar(file=sbref)

        src_space: str = "sbref"
        ref_space: str = "struct"

        sb2st_reg: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        sb2st_dict: Dict[str, str] = sb2st_reg.outputs(
            src_space=src_space, ref_space=ref_space
        )
        sb2st_files: Tuple[str] = (
            "affine",
            "inv_affine",
            "resampled_image",
            "warp",
            "dc_warp",
            "dc_image",
            "dc_brainmask",
        )

        if not sb2st_reg.check_exists(*sb2st_files):
            (
                sbref2struct_affine,
                sbref2struct_inv_affine,
                sbref2struct_src2ref,
                sbref2struct_warp,
                sbref_dc_warp,
                sbref_dc_img,
                sbref_dc_brainmask,
            ) = sbref_to_struct(
                outdir=self.outputs.get("workdir"),
                sbref=sbref,
                sbref_brainmask=self.outputs.get("sbref_brainmask"),
                struct=self.outputs.get("T2w"),
                struct_brainmask=self.outputs.get("T2w_brainmask"),
                dc=True,
                fmap=self.outputs.get("fmap"),
                fmap_brainmask=self.outputs.get("fmap_mask"),
                fmap2struct_xfm=fmap2struct_affine,
                sbref_pedir=sbref_meta.get("phase_encode_dir"),
                sbref_echospacing=sbref_meta.get("echo_spacing"),
                struct_boundarymask=self.outputs.get("T2w_wmmask"),
                bbr=True,
                src_space=src_space,
                ref_space=ref_space,
                log=fmap_log,
            )
        else:
            sbref2struct_affine: str = sb2st_dict.get("affine")
            sbref2struct_inv_affine: str = sb2st_dict.get("inv_affine")
            sbref2struct_src2ref: str = sb2st_dict.get("resampled_image")
            sbref2struct_warp: str = sb2st_dict.get("warp")
            sbref_dc_warp: str = sb2st_dict.get("dc_warp")
            sbref_dc_img: str = sb2st_dict.get("dc_image")
            sbref_dc_brainmask: str = sb2st_dict.get("dc_brainmask")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "sbref2struct_affine": sbref2struct_affine,
            "sbref2struct_inv_affine": sbref2struct_inv_affine,
            "sbref2struct_src2ref": sbref2struct_src2ref,
            "sbref2struct_warp": sbref2struct_warp,
            "sbref_dc_warp": sbref_dc_warp,
            "sbref_dc_img": sbref_dc_img,
            "sbref_dc_brainmask": sbref_dc_brainmask,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # func (distorted) -> sbref -> struct (composite)

        src_space: str = "func"
        ref_space: str = "struct"

        fn2st_reg: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        fn2st_dict: Dict[str, str] = fn2st_reg.outputs(
            src_space=src_space, ref_space=ref_space
        )
        fn2st_files: Tuple[str] = (
            "affine",
            "inv_affine",
            "warp",
            "inv_warp",
            "dc_warp",
            "dc_image",
            "resampled_image",
        )

        if not fn2st_reg.check_exists(*fn2st_files):
            (
                func2struct_affine,
                func2struct_inv_affine,
                func2struct_warp,
                func2struct_inv_warp,
                func_dc_warp,
                func_dc_img,
                func2struct_resamp_img,
            ) = func_to_struct_composite(
                outdir=self.outputs.get("workdir"),
                func=self.outputs.get("func0"),
                struct=self.outputs.get("T2w"),
                func2sbref_affine=func2sbref_affine,
                sbref2struct_affine=sbref2struct_affine,
                sbref2struct_warp=sbref2struct_warp,
                src_space=src_space,
                ref_space=ref_space,
                log=fmap_log,
            )
        else:
            func2struct_affine: str = fn2st_dict.get("affine")
            func2struct_inv_affine: str = fn2st_dict.get("inv_affine")
            func2struct_warp: str = fn2st_dict.get("warp")
            func2struct_inv_warp: str = fn2st_dict.get("inv_warp")
            func_dc_warp: str = fn2st_dict.get("dc_warp")
            func_dc_img: str = fn2st_dict.get("dc_image")
            func2struct_resamp_img: str = fn2st_dict.get("resampled_image")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "func2struct_affine": func2struct_affine,
            "func2struct_inv_affine": func2struct_inv_affine,
            "func2struct_warp": func2struct_warp,
            "func2struct_inv_warp": func2struct_inv_warp,
            "func_dc_warp": func_dc_warp,
            "func_dc_img": func_dc_img,
            "func2struct_resamp_img": func2struct_resamp_img,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # fmap -> func (composite)

        src_space: str = "fmap"
        ref_space: str = "func"

        fm2fn_reg: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        fm2fn_dict: Dict[str, str] = fm2fn_reg.outputs(
            src_space=src_space, ref_space=ref_space
        )
        fm2fn_files: Tuple[str] = ("affine", "resampled_image")

        if not fm2fn_reg.check_exists(*fm2fn_files):
            (fmap2func_affine, fmap2func_resamp_img) = fmap_to_func_composite(
                outdir=self.outputs.get("workdir"),
                fmap=self.outputs.get("fmap"),
                func=self.outputs.get("func0"),
                fmap2struct_affine=fmap2struct_affine,
                func2struct_invaffine=func2struct_inv_affine,
                src_space=src_space,
                ref_space=ref_space,
                log=fmap_log,
            )
        else:
            fmap2func_affine: str = fm2fn_dict.get("affine")
            fmap2func_resamp_img: str = fm2fn_dict.get("resampled_image")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "fmap2func_affine": fmap2func_affine,
            "fmap2func_resamp_img": fmap2func_resamp_img,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        return self.outputs

    def mcdc(
        self,
        use_mcflirt: bool = False,
        s2v: bool = False,
        dc: bool = False,
        mbs: bool = False,
        mporder: Optional[int] = None,
    ) -> Dict[Any, str]:
        """doc-string
        """
        with WorkDir(src=self.logdir) as lgd:
            _mcdc_log: str = lgd.join("mcdc.log")
            mcdc_log: LogFile = LogFile(_mcdc_log, format_log_str=True)

        self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)

        func: str = self.outputs.get("func")
        metadata: Dict[str, str] = load_sidecar(file=func)

        assert "echo_spacing" in metadata, "func sidecar must contain 'echo_spacing'."
        assert (
            "phase_encode_dir" in metadata
        ), "func sidecar must contain 'phase_encode_dir'."

        kwargs: Dict[str, str] = {
            "func": func,
            "func_brainmask": self.outputs.get("func_brainmask"),
            "fmap": self.outputs.get("fmap2func_resamp_img"),
            # "fmap2func_affine": self.outputs.get('fmap2func_affine'),
            "func_echospacing": metadata.get("echo_spacing"),
            "func_pedir": metadata.get("phase_encode_dir"),
            "func_slorder": self.outputs.get("func_slorder"),
            "inplane_acc": metadata.get("inplane_accel"),
            "mb_factor": metadata.get("mb_factor"),
            "outdir": self.outputs.get("workdir"),
            "s2v": s2v,
            "dc": dc,
            "mbs": mbs,
            "mporder": mporder,
            "use_mcflirt": use_mcflirt,
        }

        # Check if mc(dc) files exist
        outmcdc: MCDCFiles = MCDCFiles(outdir=self.outputs.get("workdir"))
        outmcdc_dict: Dict[str, str] = outmcdc.outputs(dc=(not use_mcflirt))
        outmcdc_files: Tuple[str] = (
            "func_mcdc",
            "mcdc_mean",
            "motparams",
            "mcdc_brainmask",
        )

        if not outmcdc.check_exists(*outmcdc_files):
            (
                func_mcdc,
                motparams,
                func_metrics,
                mcdc_mean,
                mcdc_std,
                mcdc_tsnr,
                mcdc_brainmask,
                fov_mask,
                fov_percent,
            ) = mcdc(**kwargs, log=mcdc_log)
        else:
            func_mcdc: str = outmcdc_dict.get("func_mcdc")
            motparams: str = outmcdc_dict.get("motparams")
            func_metrics: str = outmcdc_dict.get("func_metrics")
            mcdc_mean: str = outmcdc_dict.get("mcdc_mean")
            mcdc_std: str = outmcdc_dict.get("mcdc_std")
            mcdc_tsnr: str = outmcdc_dict.get("mcdc_tsnr")
            mcdc_brainmask: str = outmcdc_dict.get("mcdc_brainmask")
            fov_mask: str = outmcdc_dict.get("fov_mask")
            fov_percent: str = outmcdc_dict.get("fov_percent")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "func_mcdc": func_mcdc,
            "motparams": motparams,
            "func_metrics": func_metrics,
            "mcdc_mean": mcdc_mean,
            "mcdc_std": mcdc_std,
            "mcdc_tsnr": mcdc_tsnr,
            "mcdc_brainmask": mcdc_brainmask,
            "fov_mask": fov_mask,
            "fov_percent": fov_percent,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # func (undistorted) -> sbref (undistorted)

        src_space: str = "func-mcdc"
        ref_space: str = "sbref-dc"

        mc2sb_dc: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        mc2sb_dict: Dict[str, str] = mc2sb_dc.outputs(
            src_space=src_space, ref_space=ref_space
        )
        mc2sb_files: Tuple[str] = ("affine", "resampled_image", "inv_affine")

        if not mc2sb_dc.check_exists(*mc2sb_files):
            (
                func_mcdc2sbref_dc_src2ref,
                func_mcdc2sbref_dc_affine,
                func_mcdc2sbref_dc_inv_affine,
            ) = func_to_sbref(
                outdir=self.outputs.get("workdir"),
                func=mcdc_mean,
                func_brainmask=mcdc_brainmask,
                sbref=self.outputs.get("sbref_dc_img"),
                sbref_brainmask=self.outputs.get("sbref_dc_brainmask"),
                src_space=src_space,
                ref_space=ref_space,
                log=mcdc_log,
            )
        else:
            func_mcdc2sbref_dc_src2ref: str = mc2sb_dict.get("resampled_image")
            func_mcdc2sbref_dc_affine: str = mc2sb_dict.get("affine")
            func_mcdc2sbref_dc_inv_affine: str = mc2sb_dict.get("inv_affine")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "func_mcdc2sbref_dc_src2ref": func_mcdc2sbref_dc_src2ref,
            "func_mcdc2sbref_dc_affine": func_mcdc2sbref_dc_affine,
            "func_mcdc2sbref_dc_inv_affine": func_mcdc2sbref_dc_inv_affine,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # func (undistorted) -> sbref (undistorted) -> struct (composite)

        src_space: str = "func-mcdc"
        ref_space: str = "struct"

        mc2st: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        mc2st_dict: Dict[str, str] = mc2sb_dc.outputs(
            src_space=src_space, ref_space=ref_space
        )
        mc2st_files: Tuple[str] = (
            "affine",
            "inv_affine",
            "warp",
            "inv_warp",
            "dc_warp",
            "dc_image",
            "resampled_image",
        )

        if not mc2st.check_exists(*mc2st_files):
            (
                func_mcdc2struct_affine,
                func_mcdc2struct_inv_affine,
                func_mcdc2struct_warp,
                func_mcdc2struct_inv_warp,
                func_mcdc2struct_dc_warp,
                func_mcdc2struct_dc_img,
                func_mcdc2struct_resamp_img,
            ) = func_to_struct_composite(
                outdir=self.outputs.get("workdir"),
                func=mcdc_mean,
                struct=self.outputs.get("T2w"),
                func2sbref_affine=func_mcdc2sbref_dc_affine,
                sbref2struct_affine=self.outputs.get("sbref2struct_affine"),
                sbref2struct_warp=self.outputs.get("sbref2struct_warp"),
                src_space=src_space,
                ref_space=ref_space,
                log=mcdc_log,
            )
        else:
            func_mcdc2struct_affine: str = mc2st_dict.get("affine")
            func_mcdc2struct_inv_affine: str = mc2st_dict.get("inv_affine")
            func_mcdc2struct_warp: str = mc2st_dict.get("warp")
            func_mcdc2struct_inv_warp: str = mc2st_dict.get("inv_warp")
            func_mcdc2struct_dc_warp: str = mc2st_dict.get("dc_warp")
            func_mcdc2struct_dc_img: str = mc2st_dict.get("dc_image")
            func_mcdc2struct_resamp_img: str = mc2st_dict.get("resampled_image")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "func_mcdc2struct_affine": func_mcdc2struct_affine,
            "func_mcdc2struct_inv_affine": func_mcdc2struct_inv_affine,
            "func_mcdc2struct_warp": func_mcdc2struct_warp,
            "func_mcdc2struct_inv_warp": func_mcdc2struct_inv_warp,
            "func_mcdc2struct_dc_warp": func_mcdc2struct_dc_warp,
            "func_mcdc2struct_dc_img": func_mcdc2struct_dc_img,
            "func_mcdc2struct_resamp_img": func_mcdc2struct_resamp_img,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        return self.outputs

    def _reg_to_standard(
        self,
        template_source: str,
        template_space: str,
        template_age: Union[int, str],
        standard_age: Union[int, str] = 40,
        quick: bool = False,
        atlasdir: Optional[str] = None,
        std_log: Optional[LogFile] = None,
    ) -> Dict[Any, str]:
        """doc-stirng

        NOTE: 
            * Template source should be the aged match template.
            * Template space should be the standard aged template.
        """
        self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)

        try:
            age: int = int(np.round(template_age))
        except (TypeError, ValueError):
            age: Union[int, float, str] = template_age

        #  template -> struct

        # src_space: str = f'template-{age}wks'
        src_space: str = template_source
        ref_space: str = "struct"

        tm2st_reg: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        tm2st_dict: Dict[str, str] = tm2st_reg.outputs(
            src_space=src_space, ref_space=ref_space
        )
        tm2st_files: Tuple[str] = ("warp", "inv_warp", "affine", "resampled_image")

        if not tm2st_reg.check_exists(*tm2st_files):
            (
                template2struct_warp,
                template2struct_inv_warp,
                template2struct_affine,
                template2struct_src2ref,
            ) = template_to_struct(
                outdir=self.outputs.get("workdir"),
                age=age,
                struct_brainmask=self.outputs.get("T2w_brainmask"),
                struct_T2w=self.outputs.get("T2w"),
                struct_T1w=self.outputs.get("T1w"),
                quick=quick,
                src_space=src_space,
                ref_space=ref_space,
                atlasdir=atlasdir,
                log=std_log,
            )
        else:
            template2struct_warp: str = tm2st_dict.get("warp")
            template2struct_inv_warp: str = tm2st_dict.get("inv_warp")
            template2struct_affine: str = tm2st_dict.get("affine")
            template2struct_src2ref: str = tm2st_dict.get("resampled_image")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "template2struct_warp": template2struct_warp,
            "template2struct_inv_warp": template2struct_inv_warp,
            "template2struct_affine": template2struct_affine,
            "template2struct_src2ref": template2struct_src2ref,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # struct -> age-matched template -> standard template (composite)

        src_space: str = "struct"
        # ref_space: str = 'standard'
        ref_space: str = template_space

        st2std_reg: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        st2std_dict: Dict[str, str] = st2std_reg.outputs(
            src_space=src_space, ref_space=ref_space
        )
        st2std_files: Tuple[str] = ("warp", "inv_warp", "resampled_image")

        if not st2std_reg.check_exists(*st2std_files):
            (
                struct2std_warp,
                struct2std_inv_warp,
                struct2std_resamp_img,
            ) = struct_to_template_composite(
                outdir=self.outputs.get("workdir"),
                age=age,
                standard_age=standard_age,
                struct=self.outputs.get("T2w"),
                struct2template_warp=template2struct_inv_warp,
                atlasdir=atlasdir,
                src_space=src_space,
                ref_space=ref_space,
                log=std_log,
            )
        else:
            struct2std_warp: str = st2std_dict.get("warp")
            struct2std_inv_warp: str = st2std_dict.get("inv_warp")
            struct2std_resamp_img: str = st2std_dict.get("resampled_image")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "struct2std_warp": struct2std_warp,
            "struct2std_inv_warp": struct2std_inv_warp,
            "struct2std_resamp_img": struct2std_resamp_img,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # func (undistorted) -> struct -> age-matched template -> standard template (composite)

        src_space: str = "func-mcdc"
        # ref_space: str = 'standard'
        ref_space: str = template_space

        fn2tm_reg: MRreg = MRreg(outdir=self.outputs.get("workdir"))
        fn2tm_dict: Dict[str, str] = fn2tm_reg.outputs(
            src_space=src_space, ref_space=ref_space
        )
        fn2tm_files: Tuple[str] = ("warp", "inv_warp", "resampled_image")

        if not fn2tm_reg.check_exists(*fn2tm_files):
            (
                func2std_warp,
                func2std_inv_warp,
                func2std_resamp_img,
            ) = func_to_template_composite(
                outdir=self.outputs.get("workdir"),
                age=age,
                standard_age=standard_age,
                func=self.outputs.get("mcdc_mean"),
                func2struct_affine=self.outputs.get("func_mcdc2struct_affine"),
                struct2template_warp=struct2std_warp,
                atlasdir=atlasdir,
                standarddir=atlasdir,
                src_space=src_space,
                ref_space=ref_space,
                log=std_log,
            )
        else:
            func2std_warp: str = fn2tm_dict.get("warp")
            func2std_inv_warp: str = fn2tm_dict.get("inv_warp")
            func2std_resamp_img: str = fn2tm_dict.get("resampled_image")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "func2std_warp": func2std_warp,
            "func2std_inv_warp": func2std_inv_warp,
            "func2std_resamp_img": func2std_resamp_img,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        return self.outputs

    def standard(
        self,
        standard_age: Union[int, str] = 40,
        quick: bool = False,
        template_ages: Optional[Union[List[Union[int, str]], int, str]] = None,
        atlasdir: Optional[str] = None,
    ) -> Dict[Any, str]:
        """doc-string
        """
        with WorkDir(src=self.logdir) as lgd:
            _std_log: str = lgd.join("standard.log")
            std_log: LogFile = LogFile(_std_log, format_log_str=True)

        self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)

        if template_ages is not None:
            if not isinstance(template_ages, list):
                template_ages: List[str] = template_ages.split(",")

            # Create set of unique template ages
            age_set: Set = set(template_ages)

            if ("neo" in template_ages) or ("neonate" in template_ages):
                try:
                    age_set.remove("neo")
                except KeyError:
                    age_set.remove("neonate")

                # Sort
                template_ages: List[str] = list(age_set)
                template_ages.sort()
                template_ages.insert(0, "neonate")
                # template_ages.insert(0,'neo')
            else:
                # Sort
                template_ages: List[str] = list(age_set)
                template_ages.sort()
        else:
            template_ages: List[str] = []

        try:
            scan_pma: Union[float, str] = self.outputs.get("scan_pma")
            age: int = int(np.round(scan_pma))
        except (TypeError, ValueError):
            age: Union[int, float, str] = self.outputs.get("scan_pma")

        # Append scan age to age list
        template_ages.append(age)

        # Fill dictionary
        age_dict: OrderedDict = OrderedDict()

        for template_age in template_ages:
            if (template_age == "neo") or (template_age == "neonate"):
                template_source: str = f"UNCAAL-{template_age}"
                template_space: str = f"UNCAAL-{template_age}"
                temp_age: Union[int, str] = template_age
                std_age: Union[int, str] = template_age
            elif (template_age == 1) or (template_age == 2):
                template_source: str = f"UNCAAL-{template_age}yrs"
                template_space: str = f"UNCAAL-{template_age}yrs"
                temp_age: Union[int, str] = template_age
                std_age: Union[int, str] = template_age
            elif template_age == age:
                template_source: str = f"template-{template_age}wks"
                template_space: str = "standard"
                temp_age: Union[int, str] = template_age
                std_age: Union[int, str] = standard_age
            else:
                template_source: str = f"template-{template_age}wks"
                template_space: str = f"template-{template_age}wks"
                temp_age: Union[int, str] = template_age
                std_age: Union[int, str] = template_age

            tmp_dict: Dict[str, str] = {
                template_age: {
                    "template_source": template_source,
                    "template_space": template_space,
                    "template_age": temp_age,
                    "standard_age": std_age,
                }
            }

            age_dict.update(**tmp_dict)

            std_log.log(
                f"Performing registration of structural and functional data to: {age_dict.get(template_age).get('template_space')}",
                use_header=True,
            )

            _: Dict[str, str] = self._reg_to_standard(
                **age_dict.get(template_age),
                quick=quick,
                atlasdir=atlasdir,
                std_log=std_log,
            )
        return self.outputs

    def ica(
        self, temporal_fwhm: float = 150.0, icadim: Optional[int] = None
    ) -> Dict[Any, str]:
        """Perform single-subject ICA.
        """
        with WorkDir(src=self.logdir) as lgd:
            _ica_log: str = lgd.join("ica.log")
            ica_log: LogFile = LogFile(_ica_log, format_log_str=True)

        self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)

        temporal_fwhm: float = float(temporal_fwhm)

        # Check if ICA files exist
        ica_out: ICAFiles = ICAFiles(outdir=self.outputs.get("workdir"))
        ica_dict: Dict[str, str] = ica_out.outputs()
        ica_files: Tuple[str] = ("func_filt", "meldir")

        if not ica_out.check_exists(*ica_files):
            func_filt, icadir = ica(
                outdir=self.outputs.get("workdir"),
                func=self.outputs.get("func_mcdc"),
                func_brainmask=self.outputs.get("mcdc_brainmask"),
                icadim=icadim,
                temporal_fwhm=temporal_fwhm,
                log=ica_log,
            )
        else:
            func_filt: str = ica_dict.get("func_filt")
            icadir: str = ica_dict.get("icadir")

        _: str = update_sidecar(file=func_filt, temporal_fwhm=temporal_fwhm)

        # Update output dictionary
        self.outputs: Dict[str, str] = {
            **self.outputs,
            "func_filt": func_filt,
            "icadir": icadir,
            "meldir": ica_dict.get("meldir"),
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)
        return self.outputs

    def denoise(
        self, rdata: Optional[str] = None, fix_threshold: Optional[int] = None
    ) -> Dict[Any, str]:
        """Perform FIX-based denoising.
        """
        with WorkDir(src=self.logdir) as lgd:
            _fix_log: str = lgd.join("fix.log")
            fix_log: LogFile = LogFile(_fix_log, format_log_str=True)

        self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)

        struct_dseg: str = self.outputs.get("T2w_dseg")
        dseg_type: str = load_sidecar(struct_dseg).get("dseg_type")

        func_filt: str = self.outputs.get("func_filt")
        temporal_fwhm: float = load_sidecar(func_filt).get("temporal_fwhm")

        # Check if FIX feature extraction has been performed

        outfix1: FIXExtract = FIXExtract(outdir=self.outputs.get("workdir"))
        outdict1: Dict[str, str] = outfix1.outputs()
        outfiles1: Tuple[str] = ("fixdir",)

        if not outfix1.check_exists(*outfiles1):
            fixdir: str = fix_extract(
                outdir=self.outputs.get("workdir"),
                func_filt=func_filt,
                func_ref=self.outputs.get("mcdc_mean"),
                struct=self.outputs.get("T2w"),
                struct_brainmask=self.outputs.get("T2w_brainmask"),
                struct_dseg=struct_dseg,
                dseg_type=dseg_type,
                func2struct_mat=self.outputs.get("func_mcdc2struct_affine"),
                mot_param=self.outputs.get("motparams"),
                icadir=self.outputs.get("meldir"),
                temporal_fwhm=temporal_fwhm,
                log=fix_log,
            )
        else:
            fixdir: str = outdict1.get("fixdr")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "fixdir": fixdir,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # Check if FIX feature classification has been performed

        outfix2: FIXClassify = FIXClassify(outdir=self.outputs.get("workdir"))
        outdict2: Dict[str, str] = outfix2.outputs()
        outfiles2: Tuple[str] = ("fix_labels", "fix_regressors")

        if not outfix2.check_exists(*outfiles2):

            if rdata is None:
                rdata: str = DEFAULT_CLASSIFIER

            if fix_threshold is None:
                fix_threshold: int = int(self.settings("fix_threshold"))

            if (rdata is not None) and (fix_threshold is not None):
                fix_labels, fix_regressors = fix_classify(
                    outdir=self.outputs.get("workdir"),
                    rdata=rdata,
                    thr=fix_threshold,
                    log=fix_log,
                )
        else:
            fix_labels: str = outdict2.get("fix_labels")
            fix_regressors: str = outdict2.get("fix_regressors")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "fix_labels": fix_labels,
            "fix_regressors": fix_regressors,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        # Check if FIX ICA data cleaning has been performed

        outfix3: FIXApply = FIXApply(outdir=self.outputs.get("workdir"))
        outdict3: Dict[str, str] = outfix3.outputs()
        outfiles3: Tuple[str] = (
            "func_clean",
            "func_clean_mean",
            "func_clean_std",
            "func_clean_tsnr",
        )

        if not outfix3.check_exists(*outfiles3):
            (func_clean, func_clean_mean, func_clean_std, func_clean_tsnr) = fix_apply(
                outdir=self.outputs.get("workdir"),
                temporal_fwhm=temporal_fwhm,
                log=fix_log,
            )
        else:
            func_clean: str = outdict3.get("func_clean")
            func_clean_mean: str = outdict3.get("func_clean_mean")
            func_clean_std: str = outdict3.get("func_clean_std")
            func_clean_tsnr: str = outdict3.get("func_clean_tsnr")

        self.outputs: Dict[str, str] = {
            **self.outputs,
            "func_clean": func_clean,
            "func_clean_mean": func_clean_mean,
            "func_clean_std": func_clean_std,
            "func_clean_tsnr": func_clean_tsnr,
        }
        _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        return self.outputs

    def report(self, group_qc: Optional[str] = None) -> str:
        """doc-string
        """
        self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)

        if group_qc is None:
            group_qc: str = os.path.join(GROUP_QC_DIR, "grp_qc_512_anon.json")

        template: str = "individual_qc_report_template.html"
        qcdir: str = self.outputs.get("qcdir")

        qc0: Subject = Subject.from_qcdir(workdir=qcdir)
        qc0.parse(template, self.outputs.get("qc_report"), group_json=group_qc)
        return self.outputs.get("qc_report")

    def qc(
        self,
        group_map: Optional[str] = None,
        standard_age: Optional[int] = 40,
        standard_res: Optional[float] = 1.5,
        atlasdir: Optional[str] = None,
        preproc_only: Optional[bool] = None,
    ) -> Dict[Any, Any]:
        """doc-string
        """
        # Set up logging
        with WorkDir(src=self.logdir) as lgd:
            _qc_log: str = lgd.join("qc.log")
            qc_log: LogFile = LogFile(_qc_log, format_log_str=True)

        if group_map is None:
            group_map: str = os.path.join(GROUP_MAP_DIR, "group_maps.nii.gz")

        self.outputs: Dict[str, str] = json2dict(jsonfile=self.proc)
        qcdir: str = os.path.join(self.outputs.get("workdir"), "qc")

        subid: str = json2dict(jsonfile=self.outputs.get("subject_info")).get("subid")
        sesid: str = json2dict(jsonfile=self.outputs.get("subject_info")).get("sesid")

        if sesid:
            fname: str = f"sub-{subid}_ses-{sesid}_qc.html"
        else:
            sesid: str = "000"
            fname: str = f"sub-{subid}_qc.html"

        with WorkDir(src=qcdir) as qdir:
            with WorkDir(src=self.outputs.get("workdir")) as od:
                qcdir: str = qdir.abspath()
                qc_report: str = od.join(fname)
                self.outputs: Dict[str, str] = {
                    **self.outputs,
                    "qcdir": qcdir,
                    "qc_report": qc_report,
                }
                _: str = dict2json(dict=self.outputs, jsonfile=self.proc)

        qc: Subject = Subject(workdir=qcdir)

        # ==========================
        # Add subject information
        # ==========================

        scan_pma: Union[int, float] = json2dict(
            jsonfile=self.outputs.get("subject_info")
        ).get("scan_pma")
        birth_ga: Union[int, float] = json2dict(
            jsonfile=self.outputs.get("subject_info")
        ).get("birth_ga")
        qc.add_subjectinfo(
            subject_id=subid, session_id=sesid, scan_age=scan_pma, birth_age=birth_ga
        )

        # ==========================
        # Add fmap
        # ==========================

        fmap_out: FmapFiles = FmapFiles(outdir=self.outputs.get("workdir"))
        _: Dict[str, str] = fmap_out.outputs()
        fmap_files: Tuple[str] = ("fmap", "fmap_mag")

        if fmap_out.check_exists(*fmap_files):

            out_sp: ImportSpinEcho = ImportSpinEcho(outdir=self.outputs.get("workdir"))
            _: Dict[str, str] = out_sp.outputs()
            sp_files: Tuple[str] = ("spinecho",)

            if out_sp.check_exists(*sp_files):
                spinecho: str = self.outputs.get("spinecho")
            else:
                spinecho: str = None

            qc_log.log("Add fieldmap to QC")

            qc.add_fmap(
                label="fmap",
                fmap=self.outputs.get("fmap"),
                fmap_mag=self.outputs.get("fmap_mag"),
                fmap_brainmask=self.outputs.get("fmap_mask"),
                spinecho=spinecho,
            )
        else:
            qc_log.warning("Fieldmap not found")

        # ==========================
        # Registration: fmap2struct
        # ==========================

        struct_brainmask: str = os.path.join(qcdir, "struct_brainmask.nii.gz")
        sidecar: Dict[str, str] = load_sidecar(file=self.outputs.get("T2w_dseg"))

        struct_brainmask: str = create_mask(
            dseg=self.outputs.get("T2w_dseg"),
            dseg_type=sidecar.get("dseg_type"),
            labels=["wm", "gm", "sc", "cb", "bs"],
            out=struct_brainmask,
        )

        wd0: ImportStruct = ImportStruct(outdir=self.outputs.get("workdir"))
        _: Dict[str, str] = wd0.outputs()
        wd0_files: Tuple[str] = ("T2w", "T2w_brainmask", "T2w_dseg", "T2w_wmmask")

        if wd0.check_exists(*wd0_files) and os.path.exists(
            self.outputs.get("fmap2struct_src2ref")
        ):
            qc_log.log("Add fieldmap2struct to QC")

            qc.add_reg(
                label="fmap2struct",
                source=self.outputs.get("fmap2struct_src2ref"),
                ref=self.outputs.get("T2w"),
                ref_brainmask=struct_brainmask,
                ref_boundarymask=self.outputs.get("T2w_wmmask"),
                ref_dseg=self.outputs.get("T2w_wmmask"),
            )
        else:
            qc_log.warning("fmap2struct not added")

        # ==========================
        # Registration: sbref2struct
        # ==========================

        wd0: ImportStruct = ImportStruct(outdir=self.outputs.get("workdir"))
        _: Dict[str, str] = wd0.outputs()
        wd0_files: Tuple[str] = ("T2w", "T2w_brainmask", "T2w_dseg", "T2w_wmmask")

        if wd0.check_exists(*wd0_files) and os.path.exists(
            self.outputs.get("sbref2struct_src2ref")
        ):
            qc_log.log("Add sbref2struct to QC")

            qc.add_reg(
                label="sbref2struct",
                source=self.outputs.get("sbref2struct_src2ref"),
                ref=self.outputs.get("T2w"),
                ref_brainmask=struct_brainmask,
                ref_boundarymask=self.outputs.get("T2w_wmmask"),
                ref_dseg=self.outputs.get("T2w_wmmask"),
            )
        else:
            qc_log.warning("sbref2struct not added")

        # ==========================
        # Registration: func2sbref
        # ==========================

        wd0: ImportFunc = ImportFunc(outdir=self.outputs.get("workdir"))
        _: Dict[str, str] = wd0.outputs()
        wd0_files: Tuple[str] = ("sbref", "sbref_brainmask")

        if wd0.check_exists(*wd0_files) and os.path.exists(
            self.outputs.get("func2sbref_src2ref")
        ):
            qc_log.log("Add func2sbref to QC")

            qc.add_reg(
                label="func2sbref",
                source=self.outputs.get("func2sbref_src2ref"),
                ref=self.outputs.get("sbref"),
                ref_brainmask=self.outputs.get("sbref_brainmask"),
            )
        else:
            qc_log.warning("func2sbref not added")

        # ==========================
        # Registration: func-mcdc2sbref
        # ==========================

        wd0: ImportFunc = ImportFunc(outdir=self.outputs.get("workdir"))
        _: Dict[str, str] = wd0.outputs()
        wd0_files: Tuple[str] = ("sbref", "sbref_brainmask")

        if wd0.check_exists(*wd0_files) and os.path.exists(
            self.outputs.get("func_mcdc2sbref_dc_src2ref")
        ):
            qc_log.log("Add func-mcdc2sbref to QC")

            qc.add_reg(
                label="func-mcdc2sbref",
                source=self.outputs.get("func_mcdc2sbref_dc_src2ref"),
                ref=self.outputs.get("sbref"),
                ref_brainmask=self.outputs.get("sbref_brainmask"),
            )
        else:
            qc_log.warning("func-mcdc2sbref not added")

        # ==========================
        # Registration: template2struct
        # ==========================

        # scan_pma: int = int(np.round(scan_pma))

        wd0: ImportStruct = ImportStruct(outdir=self.outputs.get("workdir"))
        _: Dict[str, str] = wd0.outputs()
        wd0_files: Tuple[str] = ("T2w", "T2w_brainmask", "T2w_dseg", "T2w_wmmask")

        if wd0.check_exists(*wd0_files) and os.path.exists(
            self.outputs.get("template2struct_src2ref")
        ):
            qc_log.log("Add template2struct to QC")

            qc.add_reg(
                label="template2struct",
                source=self.outputs.get("template2struct_src2ref"),
                ref=self.outputs.get("T2w"),
                ref_brainmask=struct_brainmask,
                ref_boundarymask=self.outputs.get("T2w_wmmask"),
                ref_dseg=self.outputs.get("T2w_wmmask"),
            )
        else:
            qc_log.warning("template2struct not added")

        # ==========================
        # Prepare standard space
        # ==========================

        # Resample standard and standard_dseg to specified resolution (typically down sampled)

        atlas: Dict[str, str] = select_atlas(
            age=standard_age, standard_age=standard_age, atlasdir=atlasdir
        )

        standard: str = atlas.get("template_T2")

        if os.path.exists(standard):
            std_resamp: str = flirt(
                src=standard,
                ref=standard,
                applyisoxfm=standard_res,
                init=os.path.join(FSLDIR, "etc/flirtsch/ident.mat"),
                out=os.path.join(
                    self.outputs.get("qcdir"), f"standard-{standard_res}mm.nii.gz"
                ),
                interp="spline",
                log=qc_log,
            )

            standard: nib.Nifti1Image = nib.load(std_resamp[0])
        else:
            standard: str = None

        standard_dseg: str = atlas.get("seg")

        if os.path.exists(standard_dseg):
            std_resamp: str = flirt(
                src=standard_dseg,
                ref=standard_dseg,
                applyisoxfm=standard_res,
                init=os.path.join(FSLDIR, "etc/flirtsch/ident.mat"),
                out=os.path.join(
                    self.outputs.get("qcdir"), f"standard-desg-{standard_res}mm.nii.gz"
                ),
                interp="nearestneighbour",
                log=qc_log,
            )

            standard_dseg: nib.Nifti1Image = nib.load(std_resamp[0])
        else:
            standard_dseg: str = None

        # Load warps

        func2standard: str = self.outputs.get("func2std_warp")
        standard2func: str = self.outputs.get("func2std_inv_warp")

        if not os.path.exists(func2standard):
            func2standard: str = None

        if not os.path.exists(standard2func):
            standard2func: str = None

        # ==========================
        # Func: Raw
        # ==========================

        wd0: ImportFunc = ImportFunc(outdir=self.outputs.get("workdir"))
        _: Dict[str, str] = wd0.outputs()
        wd0_files: Tuple[str] = ("func", "func_brainmask")

        if wd0.check_exists(*wd0_files):
            qc_log.log("Add func_raw to QC")

            qc.add_func(
                func=self.outputs.get("func"),
                label="raw",
                brainmask=self.outputs.get("func_brainmask"),
                standard=standard,
                func2standard_warp=func2standard,
                template=group_map,
                template2func_warp=standard2func,
                template_dseg=standard_dseg,
                template_dseg_labels=get_dseg_labels(seg_type="drawem"),
            )
        else:
            qc_log.warning("func_raw not added")

        # ==========================
        # Func: MCDC
        # ==========================

        with NiiFile(
            src=self.outputs.get("func_mcdc"), assert_exists=True, validate_nifti=True
        ) as fm:
            _, fname, _ = fm.file_parts()

            if "dc" in fname.lower():
                dc: bool = True
            else:
                dc: bool = False

        wd0: MCDCFiles = MCDCFiles(outdir=self.outputs.get("workdir"))
        _: Dict[str, str] = wd0.outputs(dc=dc)
        wd0_files: Tuple[str] = ("func_mcdc", "mcdc_brainmask")

        if wd0.check_exists(*wd0_files):
            qc_log.log("Add func_mcdc to QC")

            qc.add_func(
                func=self.outputs.get("func_mcdc"),
                label="mcdc",
                brainmask=self.outputs.get("mcdc_brainmask"),
                standard=standard,
                func2standard_warp=func2standard,
                template=group_map,
                template2func_warp=standard2func,
                template_dseg=standard_dseg,
                template_dseg_labels=get_dseg_labels(seg_type="drawem"),
            )
        else:
            qc_log.warning("func_mcdc not added")

        # ==========================
        # Func: Clean
        # ==========================

        # Include FIX cleaned data to QC report
        if (not preproc_only) or (preproc_only is None):
            wd1: FIXApply = FIXApply(outdir=self.outputs.get("workdir"))
            _: Dict[str, str] = wd1.outputs()

            if wd1.check_exists("func_clean"):
                qc_log.log("Add func_clean to QC")

                qc.add_func(
                    func=self.outputs.get("func_clean"),
                    label="clean",
                    brainmask=self.outputs.get("mcdc_brainmask"),
                    standard=standard,
                    func2standard_warp=func2standard,
                    template=group_map,
                    template2func_warp=standard2func,
                    template_dseg=standard_dseg,
                    template_dseg_labels=get_dseg_labels(seg_type="drawem"),
                )
            else:
                qc_log.warning("func_clean not added")
        else:
            # MOCK CODE:
            #
            # This code block is used ONLY on the condition that FIX has not
            #   run. With that being the case - this used as a place holder
            #   for data that will be displayed in the QC report.
            wd1: MCDCFiles = MCDCFiles(outdir=self.outputs.get("workdir"))
            _: Dict[str, str] = wd1.outputs()

            if wd0.check_exists(*wd0_files):
                qc_log.log("Add func_clean [place holder] to QC")

                qc.add_func(
                    func=self.outputs.get("func_filt"),
                    label="clean",
                    brainmask=self.outputs.get("mcdc_brainmask"),
                    standard=standard,
                    func2standard_warp=func2standard,
                    template=group_map,
                    template2func_warp=standard2func,
                    template_dseg=standard_dseg,
                    template_dseg_labels=get_dseg_labels(seg_type="drawem"),
                )
            else:
                qc_log.warning("func_clean [place holder] not added")

        # ==========================
        # Motion parameters
        # ==========================

        if os.path.exists(self.outputs.get("motparams")):
            qc_log.log("Add motion parameters to QC")
            qc.add_motparam(label="motparam", motparams=self.outputs.get("motparams"))
        else:
            qc_log.warning("motion parameters not added")

        # ==========================
        # FIX
        # ==========================

        meldir: str = self.outputs.get("meldir")

        if os.path.exists(meldir):

            with WorkDir(src=meldir) as md:
                mel_ic: str = md.join("melodic_IC.nii.gz")
                mel_mix: str = md.join("melodic_mix")

            wd0: FIXClassify(outdir=self.outputs.get("workdir"))
            _: Dict[str, str] = wd0.outputs()

            if os.path.exists(self.outputs.get("fix_labels", "")):
                labels: str = self.outputs.get("fix_labels")
            else:
                labels: str = None

            qc_log.log("Add FIX to QC")

            qc.add_fix(label="fix", ic=mel_ic, mix=mel_mix, labels=labels)
        else:
            qc_log.warning("FIX/ICA not found")

        return self.outputs

    def pre_mcdc(self) -> Dict[Any, str]:
        """Perform all pre-motion-and-distortion-correction stages of the preprocessing pipeline."""
        self.outputs: Dict[Any, str] = self.prepare_fieldmap()
        return self.outputs

    def post_mcdc(
        self,
        standard_age: Union[int, str] = 40,
        template_ages: Optional[Union[List[Union[int, str]], int, str]] = None,
        temporal_fwhm: float = 150.0,
        quick: bool = False,
        icadim: Optional[int] = None,
        rdata: Optional[str] = None,
        fix_threshold: Optional[int] = None,
        group_map: Optional[str] = None,
        standard_res: Optional[float] = 1.5,
        group_qc: Optional[str] = None,
        atlasdir: Optional[str] = None,
        preproc_only: Optional[bool] = None,
        spatial_fwhm: Optional[float] = None,
        intnorm: bool = False,
    ) -> Dict[Any, str]:
        """Perform all post-motion-and-distortion-correction stages of the preprocessing pipeline."""
        self.outputs: Dict[Any, str] = self.standard(
            standard_age=standard_age,
            quick=quick,
            template_ages=template_ages,
            atlasdir=atlasdir,
        )

        self.outputs: Dict[Any, str] = self.ica(
            temporal_fwhm=temporal_fwhm, icadim=icadim
        )

        if not preproc_only:
            self.outputs: Dict[Any, str] = self.denoise(
                rdata=rdata, fix_threshold=fix_threshold
            )

        self.outputs: Dict[Any, str] = self.qc(
            group_map=group_map,
            standard_res=standard_res,
            standard_age=standard_age,
            preproc_only=preproc_only,
        )

        _: str = self.report(group_qc=group_qc)

        if (spatial_fwhm is not None) and (spatial_fwhm != 0) or intnorm:

            if spatial_fwhm is None:
                spatial_fwhm: int = 0

            with WorkDir(src=self.logdir) as lgd:
                _post_log: str = lgd.join("postprocess.log")
                post_log: LogFile = LogFile(_post_log, format_log_str=True)

            if os.path.exists(self.outputs.get("func_clean", "")):
                func: str = self.outputs.get("func_clean")
            else:
                func: str = self.outputs.get("func_filt")

            _: Tuple[str] = postprocess(
                func=func,
                func_mean=self.outputs.get("mcdc_mean"),
                func_brainmask=self.outputs.get("mcdc_brainmask"),
                outdir=self.outputs.get("workdir"),
                spatial_fwhm=spatial_fwhm,
                intnorm=intnorm,
                log=post_log,
            )

        return self.outputs

    def run_all(
        self,
        standard_age: Union[int, str] = 40,
        template_ages: Optional[Union[List[Union[int, str]], int, str]] = None,
        temporal_fwhm: float = 150.0,
        use_mcflirt: bool = False,
        s2v: bool = False,
        dc: bool = False,
        mbs: bool = False,
        mporder: Optional[int] = None,
        quick: bool = False,
        icadim: Optional[int] = None,
        rdata: Optional[str] = None,
        fix_threshold: Optional[int] = None,
        group_map: Optional[str] = None,
        standard_res: Optional[float] = 1.5,
        group_qc: Optional[str] = None,
        atlasdir: Optional[str] = None,
        preproc_only: Optional[bool] = None,
        spatial_fwhm: Optional[float] = None,
        intnorm: bool = False,
    ) -> None:
        """Perform all stages of the preprocessing pipeline."""
        self.outputs: Dict[Any, str] = self.pre_mcdc()

        self.outputs: Dict[Any, str] = self.mcdc(
            use_mcflirt=use_mcflirt, s2v=s2v, dc=dc, mbs=mbs, mporder=mporder
        )

        self.outputs: Dict[Any, str] = self.post_mcdc(
            standard_age=standard_age,
            template_ages=template_ages,
            temporal_fwhm=temporal_fwhm,
            quick=quick,
            icadim=icadim,
            rdata=rdata,
            fix_threshold=fix_threshold,
            group_map=group_map,
            standard_res=standard_res,
            group_qc=group_qc,
            atlasdir=atlasdir,
            preproc_only=preproc_only,
            spatial_fwhm=spatial_fwhm,
            intnorm=intnorm,
        )

        return self.outputs
