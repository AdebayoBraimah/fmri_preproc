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
from fmri_preproc.utils.util import load_sidecar

from fmri_preproc.utils.outputs.importdata import (
    ImportFunc,
    ImportSpinEcho,
    ImportStruct
)

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
    fmap_to_struct,
    fmap_to_func_composite,
    func_to_sbref,
    func_to_struct_composite,
    sbref_to_struct,
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
        return self.outputs

    def prepare_fieldmap(self) -> None:
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
        # fmap -> func (composite)

        # TODO: Add fmap/sbref reg outputs to class dict.
        pass

    def mcdc(self) -> None:
        pass

    def standard(self) -> None:
        pass

    def ica(self) -> None:
        pass

    def denoise(self) -> None:
        pass

    def report(self) -> None:
        pass

    def qc(self) -> None:
        pass

    def pre_mcdc(self) -> None:
        """Perform all pre-motion-and-distortion-correction stages of the preprocessing pipeline."""
        pass

    def post_mcdc(self) -> None:
        """Perform all post-motion-and-distortion-correction stages of the preprocessing pipeline."""
        pass

    def run_all(self) -> None:
        """Perform all stages of the preprocessing pipeline."""
        pass
