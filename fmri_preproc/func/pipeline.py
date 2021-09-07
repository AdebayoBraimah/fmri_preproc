# -*- coding: utf-8 -*-
"""fMRI preprocessing pipeline - main pipeline.
"""
# NOTE:
# 
# Pipeline overview
# 
# 0. Import/copy data to working directory
# 
# 1. [X] Prepare fieldmap
#   a. [X] Topup [mc]
# 2. [X] Motion correction, distortion correction
#   a. [X] eddy [mc]
#   b. [X] mcflirt [mc]
# 3. Registration (linear, and non-linear)
#   a. flirt ( + mcflirt) [reg]
#       i. sym-link field map files from [mc]
# 4. [X] ICA
#   a. [X] melodic
# 5. [X] Denoise
#   a. [X] FIX
# 6. [X] Post-processing
#   a. [X] Spatial smooting
#   b. [X] Intensity norm
# 7. QC
#   a. PENDING
import os
import numpy as np

from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
    Union
)

from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.util import (
    dict2json,
    json2dict,
    load_sidecar, 
    update_sidecar
)

from fmri_preproc.utils.outputs.importdata import (
    ImportFunc,
    ImportSpinEcho,
    ImportStruct
)

from fmri_preproc.utils.outputs.outdict import key_to_none
from fmri_preproc.utils.outputs.fieldmap import FmapFiles
from fmri_preproc.utils.outputs.registration import MRreg

from fmri_preproc.func.importdata import (
    import_func,
    import_info,
    import_spinecho,
    import_struct
)

from fmri_preproc.func.fieldmap import fieldmap

from fmri_preproc.func.registration import (
    ATLASDIR,
    fmap_to_struct,
    fmap_to_func_composite,
    func_to_sbref,
    func_to_struct_composite,
    func_to_template_composite,
    sbref_to_struct,
    struct_to_template_composite,
    template_to_struct
)

from fmri_preproc.func.mcdc import mcdc
from fmri_preproc.func.ica import ica

from fmri_preproc.func.denoise import (
    fix_apply,
    fix_classify,
    fix_extract
)

# TODO: 
#   * Create outputs module that handles output
#       dictionaries for all of the other modules.
#   * Import the dictionaries from the outputs modules.

class Pipeline:
    def __init__(self) -> None:
        """Constructor for pipeline class.
        """
        self.outputs: Dict[str,str] = {}
        return None
    
    def import_data(self,
                    outdir: str,
                    func: str,
                    func_echospacing: float,
                    func_pedir: str,
                    T2w: str,
                    T2w_brainmask: str,
                    dseg: str,
                    scan_pma: float,
                    birth_ga: Optional[float] = None,
                    verbose: bool = False,
                    log_level: str = 'info',
                    func_brainmask: Optional[str] = None,
                    func_slorder: Optional[str] = None,
                    func_inplane_accel: Optional[float] = 1,
                    mb_factor: Optional[int] = None,
                    sbref: Optional[str] = None,
                    sbref_brainmask: Optional[str] = None,
                    sbref_echospacing: Optional[float] = None,
                    sbref_pedir: Optional[str] = None,
                    dseg_type: Optional[str] = 'drawem',
                    probseg: Optional[str] = None,
                    probseg_type: Optional[str] = 'drawem',
                    wmmask: Optional[str] = None,
                    T1w: Optional[str] = None,
                    mask_func: bool = False,
                    spinecho: Optional[str] = None,
                    spinecho_echospacing: Optional[float] = 0.1,
                    spinecho_pedir: Union[str,List[str]] = None,
                    ap_dir: Optional[str] = None,
                    pa_dir: Optional[str] = None,
                    lr_dir: Optional[str] = None,
                    rl_dir: Optional[str] = None,
                    is_dir: Optional[str] = None,
                    si_dir: Optional[str] = None,
                    spinecho_epifactor: Optional[int] = None,
                    spinecho_inplaneacc: Optional[float] = 1,
                    **kwargs
                   ) -> Dict[Any,str]:
        """Import data class-function.
        """
        # Import information
        self.verbose: bool = verbose

        (sub_workdir, 
        logdir, 
        sub_json, 
        sub_dict) = import_info(outdir=outdir,
                                func=func,
                                scan_pma=scan_pma,
                                birth_ga=birth_ga,
                                log=None,
                                verbose=self.verbose,
                                log_datetime=True,
                                log_level=log_level,
                                **kwargs)
        
        self.logdir: str = logdir
        
        with WorkDir(src=self.logdir) as lgd:
            _import_log: str = lgd.join('import.log')
            import_log: LogFile = LogFile(_import_log, format_log_str=True, print_to_screen=self.verbose)

            subid: str = sub_dict.get('subid')
            sesid: str = sub_dict.get('sesid')

            if sesid:
                self.proc: str = lgd.join(f'sub-{subid}_ses-{sesid}_fmripreproc.json')
            else:
                self.proc: str = lgd.join(f'sub-{subid}_fmripreproc.json')
            
            if os.path.exists(self.proc):
                self.outputs: Dict[str,str] = json2dict(jsonfile=self.proc)
                self.outputs: Dict[str,str] = key_to_none(d=self.outputs)
                return self.outputs
        
        # Check if func has been imported
        out_func: ImportFunc = ImportFunc(outdir=outdir)
        func_files: Tuple[str] = ('func',
                                'func_mean',
                                'func_brainmask',
                                'func_slorder')

        if not out_func.check_exists(*func_files):
            (func,
            func_mean, 
            func_brainmask, 
            func_slorder, 
            sbref, 
            sbref_brainmask,
            func_info_dict) = import_func(outdir=sub_workdir,
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
                                        log=import_log)
        else:
            func_info_dict: Dict[str,str] = out_func.outputs()
        
        # Check if struct has been imported
        out_struct: ImportStruct = ImportStruct(outdir=outdir)
        struct_files: Tuple[str] = ('T2w',
                                    'T2w_brainmask',
                                    'T2w_dseg')

        if not out_struct.check_exists(*struct_files):
            (T2w, 
            wmmask, 
            dseg,
            struct_info_dict) = import_struct(outdir=sub_workdir,
                                            T2w=T2w,
                                            brainmask=T2w_brainmask,
                                            dseg=dseg,
                                            dseg_type=dseg_type,
                                            probseg=probseg,
                                            probseg_type=probseg_type,
                                            wmmask=wmmask,
                                            T1w=T1w,
                                            log=import_log)
        else:
            struct_info_dict: Dict[str,str] = out_struct.outputs()
        
        # Check if fmap has been imported
        out_fmap: ImportSpinEcho = ImportSpinEcho(outdir=outdir)
        fmap_files: Tuple[str] = ('spinecho')

        if not out_fmap.check_exists(*fmap_files):
            if (spinecho or
                (ap_dir and pa_dir) or
                (lr_dir and rl_dir) or
                (is_dir and si_dir)):
                (spinecho, spinecho_pedir) = import_spinecho(outdir=sub_workdir,
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
                                                             log=import_log)
            else:
                spinecho: str = None
                spinecho_pedir: str = None
        else:
            fmap_info_dict: Dict[str,str] = out_fmap.outputs
            spinecho: str = fmap_info_dict.get('spinecho')

        self.outputs: Dict[str,str] = {
            **sub_dict,
            **func_info_dict,
            **struct_info_dict,
            "workdir": sub_workdir,
            "func_json": sub_json,
            "spinecho": spinecho,
            "spinecho_pedir": spinecho_pedir
        }

        # TODO: Add read/writing of outputs dict to json - then read in json
        #   for cases when pipeline has stopped.

        # Set files that do not exist to None
        self.outputs: Dict[str,str] = key_to_none(d=self.outputs)
        self.proc: str = dict2json(dict=self.outputs)
        return self.outputs

    def prepare_fieldmap(self) -> Dict[Any,str]:
        """doc-string
        """
        with WorkDir(src=self.logdir) as lgd:
            _fmap_log: str = lgd.join('fmap.log')
            fmap_log: LogFile = LogFile(_fmap_log, format_log_str=True)
        
        spinecho: str = self.outputs.get('spinecho')
        metadata: Dict[str,str] = load_sidecar(file=spinecho)

        fmap_out: FmapFiles = FmapFiles(outdir=self.outputs.get('workdir'))
        fmap_info_dict: Dict[str,str] = fmap_out.outputs()

        _: Tuple[str] = fieldmap(outdir=self.outputs.get('workdir'),
                                 spinecho=spinecho,
                                 echo_spacing=metadata.get('echo_spacing'),
                                 pedir=metadata.get('phase_encode_dir'),
                                 epifactor=metadata.get('epi_factor'),
                                 inplane_acc=metadata.get('inplane_accel',1.0),
                                 verbose=self.verbose,
                                 log=fmap_log)
        
        self.outputs: Dict[str,str] = {
            **self.outputs,
            **fmap_info_dict,
        }

        # fmap -> struct
        
        (fmap2struct_affine,
        fmap2struct_inv_affine,
        fmap2struct_src2ref,
        fmap2struct_init_affine) = fmap_to_struct(outdir=self.outputs.get('workdir'),
                                                  fmap=self.outputs.get('fmap'),
                                                  fmap_magnitude=self.outputs.get('fmap_mag'),
                                                  fmap_brainmask=self.outputs.get('fmap_mask'),
                                                  struct=self.outputs.get('T2w'),
                                                  struct_brainmask=self.outputs.get('T2w_brainmask'),
                                                  struct_boundarymask=self.outputs.get('T2w_wmmask'),
                                                  src_space='fmap',
                                                  ref_space='struct',
                                                  bbr=True,
                                                  log=fmap_log)

        # func -> sbref (distorted)
        (func2sbref_src2ref, 
        func2sbref_affine, 
        func2sbref_inv_affine) = func_to_sbref(outdir=self.outputs.get('workdir'),
                                               func=self.outputs.get('func0'),
                                               func_brainmask=self.outputs.get('func_brainmask'),
                                               sbref=self.outputs.get('sbref'),
                                               sbref_brainmask=self.outputs.get('sbref_brainmask'),
                                               src_space='func',
                                               ref_space='sbref',
                                               log=fmap_log)

        # sbref -> struct (with BBR and DC)
        sbref: str = self.outputs.get('sbref')
        sbref_meta: Dict[str,str] = load_sidecar(file=sbref)

        (sbref2struct_affine, 
        sbref2struct_inv_affine, 
        sbref2struct_src2ref, 
        sbref2struct_warp,
        sbref_dc_warp, 
        sbref_dc_img, 
        sbref_dc_brainmask) = sbref_to_struct(outdir=self.outputs.get('workdir'),
                                              sbref=sbref,
                                              sbref_brainmask=self.outputs.get('sbref_brainmask'),
                                              struct=self.outputs.get('T2w'),
                                              struct_brainmask=self.outputs.get('T2w_brainmask'),
                                              dc=True,
                                              fmap=self.outputs.get('fmap'),
                                              fmap_brainmask=self.outputs.get('fmap_mask'),
                                              fmap2struct_xfm=fmap2struct_affine,
                                              sbref_pedir=sbref_meta.get('phase_encode_dir'),
                                              sbref_echospacing=sbref_meta.get('echo_spacing'),
                                              struct_boundarymask=self.outputs.get('T2w_wmmask'),
                                              bbr=True,
                                              src_space='sbref',
                                              ref_space='struct',
                                              log=fmap_log)
        
        # func (distorted) -> sbref -> struct (composite)
        (func2struct_affine,
        func2struct_inv_affine,
        func2struct_warp,
        func2struct_inv_warp,
        func_dc_warp,
        func_dc_img,
        func2struct_resamp_img) = func_to_struct_composite(outdir=self.outputs.get('workdir'),
                                                           func=self.outputs.get('func0'),
                                                           struct=self.outputs.get('T2w'),
                                                           func2sbref_affine=func2sbref_affine,
                                                           sbref2struct_affine=sbref2struct_affine,
                                                           sbref2struct_warp=sbref2struct_warp,
                                                           src_space='func',
                                                           ref_space='struct',
                                                           log=fmap_log)

        # fmap -> func (composite)
        (fmap2func_affine, 
        fmap2func_resamp_img) = fmap_to_func_composite(outdir=self.outputs.get('workdir'),
                                                       fmap=self.outputs.get('fmap'),
                                                       func=self.outputs.get('func0'),
                                                       fmap2struct_affine=fmap2struct_affine,
                                                       func2struct_invaffine=func2struct_inv_affine,
                                                       src_space='fmap',
                                                       ref_space='func',
                                                       log=fmap_log)

        # Update output dictionary
        self.outputs: Dict[str,str] = {
            **self.outputs,
            "fmap2struct_affine": fmap2struct_affine,
            "fmap2struct_inv_affine": fmap2struct_inv_affine,
            "fmap2struct_src2ref": fmap2struct_src2ref,
            "fmap2struct_init_affine": fmap2struct_init_affine,
            "func2sbref_src2ref": func2sbref_src2ref,
            "func2sbref_affine": func2sbref_affine,
            "func2sbref_inv_affine": func2sbref_inv_affine,
            "sbref2struct_affine": sbref2struct_affine,
            "sbref2struct_inv_affine": sbref2struct_inv_affine,
            "sbref2struct_src2ref": sbref2struct_src2ref,
            "sbref2struct_warp": sbref2struct_warp,
            "sbref_dc_warp": sbref_dc_warp,
            "sbref_dc_img": sbref_dc_img,
            "sbref_dc_brainmask": sbref_dc_brainmask,
            "func2struct_affine": func2struct_affine,
            "func2struct_inv_affine": func2struct_inv_affine,
            "func2struct_warp": func2struct_warp,
            "func2struct_inv_warp": func2struct_inv_warp,
            "func_dc_warp": func_dc_warp,
            "func_dc_img": func_dc_img,
            "func2struct_resamp_img": func2struct_resamp_img,
            "fmap2func_affine": fmap2func_affine, 
            "fmap2func_resamp_img": fmap2func_resamp_img
        }
        self.proc: str = dict2json(dict=self.outputs)
        return self.outputs

    def mcdc(self,
             use_mcflirt: bool = False
            ) -> Dict[Any,str]:
        """doc-string
        """
        with WorkDir(src=self.logdir) as lgd:
            _mcdc_log: str = lgd.join('mcdc.log')
            mcdc_log: LogFile = LogFile(_mcdc_log, format_log_str=True)

        func: str = self.outputs.get('func')
        metadata: Dict[str,str] = load_sidecar(file=func)

        assert 'echo_spacing' in metadata, "func sidecar must contain 'echo_spacing'."
        assert 'phase_encode_dir' in metadata, "func sidecar must contain 'phase_encode_dir'."

        kwargs: Dict[str,str] = {
            "func": func,
            "func_brainmask": self.outputs.get('func_brainmask'),
            "fmap": self.outputs.get('fmap2func_resamp_img'),
            # "fmap2func_affine": self.outputs.get('fmap2func_affine'),
            "func_echospacing": metadata.get('echo_spacing'),
            "func_pedir": metadata.get('phase_encode_dir'),
            "func_slorder": self.outputs.get('func_slorder'),
            "inplane_acc": self.outputs.get('inplane_accel'),
        }

        (func_mcdc,
        motparams,
        func_metrics,
        mcdc_mean,
        mcdc_std,
        mcdc_tsnr,
        mcdc_brainmask,
        fov_mask,
        fov_percent) = mcdc(**kwargs,
                            use_mcflirt=use_mcflirt,
                            log=mcdc_log)

        # func (undistorted) -> sbref (undistorted)
        (func_mcdc2sbref_dc_src2ref, 
        func_mcdc2sbref_dc_affine, 
        func_mcdc2sbref_dc_inv_affine) = func_to_sbref(outdir=self.outputs.get('workdir'),
                                                       func=mcdc_mean,
                                                       func_brainmask=mcdc_brainmask,
                                                       sbref=self.outputs.get('sbref_dc_img'),
                                                       sbref_brainmask=self.outputs.get('sbref_dc_brainmask'),
                                                       src_space='func-mcdc',
                                                       ref_space='sbref-dc',
                                                       log=mcdc_log)

        # func (undistorted) -> sbref (undistorted) -> struct (composite)
        (func_mcdc2struct_affine,
        func_mcdc2struct_inv_affine,
        func_mcdc2struct_warp,
        func_mcdc2struct_inv_warp,
        func_mcdc2struct_dc_warp,
        func_mcdc2struct_dc_img,
        func_mcdc2struct_resamp_img) = func_to_struct_composite(outdir=self.outputs.get('workdir'),
                                                                func=mcdc_mean,
                                                                struct=self.outputs.get('T2w'),
                                                                func2sbref_affine=func_mcdc2sbref_dc_affine,
                                                                sbref2struct_affine=self.outputs.get('sbref2struct_affine'),
                                                                sbref2struct_warp=self.outputs.get('sbref2struct_warp'),
                                                                src_space='func-mcdc',
                                                                ref_space='struct',
                                                                log=mcdc_log)

        # Update output dictionary
        self.outputs: Dict[str,str] = {
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
            "func_mcdc2sbref_dc_src2ref": func_mcdc2sbref_dc_src2ref,
            "func_mcdc2sbref_dc_affine": func_mcdc2sbref_dc_affine,
            "func_mcdc2sbref_dc_inv_affine": func_mcdc2sbref_dc_inv_affine,
            "func_mcdc2struct_affine": func_mcdc2struct_affine,
            "func_mcdc2struct_inv_affine": func_mcdc2struct_inv_affine,
            "func_mcdc2struct_warp": func_mcdc2struct_warp,
            "func_mcdc2struct_inv_warp": func_mcdc2struct_inv_warp,
            "func_mcdc2struct_dc_warp": func_mcdc2struct_dc_warp,
            "func_mcdc2struct_dc_img": func_mcdc2struct_dc_img,
            "func_mcdc2struct_resamp_img": func_mcdc2struct_resamp_img,
        }
        self.proc: str = dict2json(dict=self.outputs)
        return self.outputs,

    def standard(self,
                 standard_age: Union[int,str] = 40,
                 quick: bool = False,
                #  template_ages: Optional[Union[List[int,str],int,str]] = None,
                 atlasdir: Optional[str] = None
                ) -> Dict[Any,str]:
        """doc-string
        """
        with WorkDir(src=self.logdir) as lgd:
            _std_log: str = lgd.join('standard.log')
            std_log: LogFile = LogFile(_std_log, format_log_str=True)

        # if not isinstance(template_ages, list):
        #     template_ages: List[int,str] = [template_ages]
        
        # for template_age in template_ages:
        #     try:
        #         template_age: int = int(template_age)
        #     except(TypeError,ValueError):
        #         template_age: str = template_age
        
        # try:
        #     template_ages.remove(standard_age)
        # except ValueError:
        #     pass

        # template_ages.append(standard_age)
        
        try:
            scan_pma: Union[float,str] = self.outputs.get('scan_pma')
            age: int = int(np.round(scan_pma))
        except (TypeError,ValueError):
            age: Union[int,float,str] = scan_pma
        
        if atlasdir is None:
            atlasdir: str = ATLASDIR
        
        #  template -> struct
        (template2struct_warp,
        template2struct_inv_warp,
        template2struct_affine,
        template2struct_src2ref) = template_to_struct(outdir=self.outputs.get('workdir'),
                                                      age=age,
                                                      struct_brainmask=self.outputs.get('T2w_brainmask'),
                                                      struct_T2w=self.outputs.get('T2w'),
                                                      struct_T1w=self.outputs.get('T1w'),
                                                      quick=quick,
                                                      src_space=f'template-{age}',
                                                      ref_space='struct',
                                                      atlasdir=atlasdir,
                                                      log=std_log)

        # struct -> age-matched template -> standard template (composite)
        (struct2std_warp, 
        struct2std_inv_warp, 
        struct2std_resamp_img) = struct_to_template_composite(outdir=self.outputs.get('workdir'),
                                                              age=age,
                                                              standard_age=standard_age,
                                                              struct=self.outputs.get('T2w'),
                                                              struct2template_warp=template2struct_inv_warp,
                                                              atlasdir=atlasdir,
                                                              src_space='struct',
                                                              ref_space='standard',
                                                              log=std_log)

        # func (undistorted) -> struct -> age-matched template -> standard template (composite)
        (func2std_warp, 
        func2std_inv_warp, 
        func2std_resamp_img) = func_to_template_composite(outdir=self.outputs.get('workdir'),
                                                          age=age,
                                                          standard_age=standard_age,
                                                          func=self.outputs.get('mcdc_mean'),
                                                          struct2template_warp=struct2std_warp,
                                                          atlasdir=atlasdir,
                                                          standarddir=atlasdir,
                                                          src_space='func-mcdc',
                                                          ref_space='standard',
                                                          log=std_log)
        # Update output dictionary
        self.outputs: Dict[str,str] = {
            **self.outputs,
            "template2struct_warp": template2struct_warp,
            "template2struct_inv_warp": template2struct_inv_warp,
            "template2struct_affine": template2struct_affine,
            "template2struct_src2ref": template2struct_src2ref,
            "struct2std_warp": struct2std_warp,
            "struct2std_inv_warp": struct2std_inv_warp,
            "struct2std_resamp_img": struct2std_resamp_img,
            "func2std_warp": func2std_warp,
            "func2std_inv_warp": func2std_inv_warp,
            "func2std_resamp_img": func2std_resamp_img,
        }
        self.proc: str = dict2json(dict=self.outputs)
        return self.outputs

    def ica(self,
            temporal_fwhm: float = 150.0,
            icadim: Optional[int] = None
           ) -> Dict[Any,str]:
        """Perform single-subject ICA.
        """
        with WorkDir(src=self.logdir) as lgd:
            _ica_log: str = lgd.join('ica.log')
            ica_log: LogFile = LogFile(_ica_log, format_log_str=True)
        
        temporal_fwhm: float = float(temporal_fwhm)

        func_filt, icadir = ica(outdir=self.outputs.get('workdir'),
                                func=self.outputs.get('func_mcdc'),
                                func_brainmask=self.outputs.get('mcdc_brainmask'),
                                icadim=icadim,
                                temporal_fwhm=temporal_fwhm,
                                log=ica_log)
        
        _: str = update_sidecar(file=func_filt,
                                temporal_fwhm=temporal_fwhm)

        # Update output dictionary
        self.outputs: Dict[str,str] = {
            **self.outputs,
            "func_filt": func_filt,
            "icadir": icadir,
        }
        self.proc: str = dict2json(dict=self.outputs)
        return self.outputs

    def denoise(self,
                rdata: Optional[str] = None,
                fix_threshold: Optional[int] = None
               ) -> Dict[Any,str]:
        """Perform FIX-based denoising.
        """
        with WorkDir(src=self.logdir) as lgd:
            _fix_log: str = lgd.join('fix.log')
            fix_log: LogFile = LogFile(_fix_log, format_log_str=True)

        struct_dseg: str = self.outputs.get('T2w_dseg')
        dseg_type: str = load_sidecar(struct_dseg).get('dseg_type')

        func_filt: str = self.outputs.get('func_filt')
        temporal_fwhm: float = load_sidecar(func_filt).get('temporal_fwhm')

        fixdir: str = fix_extract(outdir=self.outputs.get('workdir'),
                                  func_filt=func_filt,
                                  func_ref=self.outputs.get('mcdc_mean'),
                                  struct=self.outputs.get('T2w'),
                                  struct_brainmask=self.outputs.get('T2w_brainmask'),
                                  struct_dseg=struct_dseg,
                                  dseg_type=dseg_type,
                                  func2struct_mat=self.outputs.get('func_mcdc2struct_affine'),
                                  mot_param=self.outputs.get('motparams'),
                                  icadir=self.outputs.get('icadir'),
                                  temporal_fwhm=temporal_fwhm,
                                  log=fix_log)

        if rdata is None:
            # Raise exception for now.
            raise RuntimeError ("Rdata file is required.")
        
        if fix_threshold is None:
            # Raise exception for now - grab this from settings file.
            raise RuntimeError("Must define FIX threshold.")
        
        if (rdata is not None) and (fix_threshold is not None):
            fix_labels, fix_reg = fix_classify(outdir=self.outputs.get('workdir'),
                                               rdata=rdata,
                                               thr=fix_threshold,
                                               log=fix_log)

            (func_clean,
            func_mean,
            func_stdev,
            func_tsnr) = fix_apply(outdir=self.outputs.get('workdir'),
                                   temporal_fwhm=temporal_fwhm,
                                   log=fix_log)
        # Update output dictionary
        self.outputs: Dict[str,str] = {
            **self.outputs,
            "fixdir": fixdir,
            "fix_labels": fix_labels,
            "fix_reg": fix_reg,
            "func_clean": func_clean,
            "func_mean": func_mean,
            "func_stdev": func_stdev,
            "func_tsnr": func_tsnr
        }
        self.proc: str = dict2json(dict=self.outputs)
        return self.outputs
    
    # TODO: Work on these class function later
    # 
    # def report(self) -> None:
    #     pass
    # 
    # def qc(self) -> None:
    #     pass

    def pre_mcdc(self) -> Dict[Any,str]:
        """Perform all pre-motion-and-distortion-correction stages of the preprocessing pipeline."""
        self.outputs: Dict[Any,str] = self.prepare_fieldmap()
        return self.outputs

    def post_mcdc(self,
                  standard_age: Union[int,str] = 40,
                  temporal_fwhm: float = 150.0,
                  icadim: Optional[int] = None,
                  rdata: Optional[str] = None,
                  fix_threshold: Optional[int] = None,
                  atlasdir: Optional[str] = None
                 ) -> Dict[Any,str]:
        """Perform all post-motion-and-distortion-correction stages of the preprocessing pipeline."""
        self.outputs: Dict[Any,str] = self.standard(standard_age=standard_age,
                                                    quick=True, # Remove this later after testing
                                                    atlasdir=atlasdir)

        self.outputs: Dict[Any,str] = self.ica(temporal_fwhm=temporal_fwhm,
                                               icadim=icadim)

        self.outputs: Dict[Any,str] = self.denoise(rdata=rdata,
                                                   fix_threshold=fix_threshold)
        return self.outputs

    def run_all(self,
                standard_age: Union[int,str] = 40,
                temporal_fwhm: float = 150.0,
                use_mcflirt: bool = False,
                icadim: Optional[int] = None,
                rdata: Optional[str] = None,
                fix_threshold: Optional[int] = None,
                atlasdir: Optional[str] = None
               ) -> None:
        """Perform all stages of the preprocessing pipeline."""
        self.outputs: Dict[Any,str] = self.prepare_fieldmap()
        self.outputs: Dict[Any,str] = self.mcdc(use_mcflirt=use_mcflirt)
        self.outputs: Dict[Any,str] = self.standard(standard_age=standard_age,
                                                    atlasdir=atlasdir)
        self.outputs: Dict[Any,str] = self.ica(temporal_fwhm=temporal_fwhm,
                                               icadim=icadim)

        self.outputs: Dict[Any,str] = self.denoise(rdata=rdata,
                                                   fix_threshold=fix_threshold)
        return self.outputs
