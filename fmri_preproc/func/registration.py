# -*- coding: utf-8 -*-
"""Image registration for the fmri_preproce pipeline.

NOTE: 
    Non-FSL dependencies:
        * ``ANTs``: http://stnava.github.io/ANTs/ (**MUST** be compiled natively on target OS)
        * Convert 3D ``c3d``: https://sourceforge.net/p/c3d/git/ci/master/tree/doc/c3d.md 
            * Can download pre-compiled binaries: https://sourceforge.net/projects/c3d/files/c3d/
"""
import os
import glob

from typing import (
    Dict,
    List,
    Optional, 
    Tuple,
    Union
)

from fmri_preproc import _resourcedir
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.command import Command
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.enums import (
    BBRType,
    PhaseEncodeDirection    
)

from fmri_preproc.utils.acqparam import (
    _get_axis as get_axis, 
    phase_encode_dict
)

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

from fmri_preproc.utils.fslpy import (
    FSLDIR,
    applywarp,
    applyxfm,
    concatxfm,
    convertwarp,
    flirt,
    fslmaths,
    fslmerge,
    fslroi,
    fugue,
    invwarp,
    invxfm
)

# TODO: ALL registration functions will need:
#   * src_space
#   * ref_space
#   * outdir
# 
# These variables will be used as inputs for all 
# registration functions.

ATLASDIR: str = os.path.join(_resourcedir,'atlases')

def fmap_to_struct(outdir: str,
                   fmap: str,
                   fmap_magnitude: str,
                   fmap_brainmask: str,
                   struct: str,
                   struct_brainmask: str,
                   struct_boundarymask: Optional[str] = None,
                   src_space: Optional[str] = None,
                   ref_space: Optional[str] = None,
                   bbr: bool = False,
                   bbr_slope: float = 0.5,
                   bbr_type: Optional[str] = "signed",
                   log: Optional[LogFile] = None
                  ) -> Tuple[str,str,str,str]:
    """Perform (rigid-body) fieldmap to structural image registration, with optional BBR.
    """
    # Check and validate input fmap associated NIFTI files
    with NiiFile(src=fmap, assert_exists=True, validate_nifti=True) as fmp:
        with NiiFile(src=fmap_magnitude, assert_exists=True, validate_nifti=True) as fmm:
            with NiiFile(src=fmap_brainmask, assert_exists=True, validate_nifti=True) as fmb:
                fmap: str = fmp.abspath()
                fmap_magnitude: str = fmm.abspath()
                fmap_brainmask: str = fmb.abspath()
                if not src_space:
                    _, src_space, _ = fmp.file_parts()
    
    # Check and validate input structural associated NIFTI images
    with NiiFile(src=struct, assert_exists=True, validate_nifti=True) as stc:
        with NiiFile(src=struct_brainmask, assert_exists=True, validate_nifti=True) as sbm:
            struct: str = stc.abspath()
            struct_brainmask: str = sbm.abspath()
            if not ref_space:
                _, ref_space, _ = stc.file_parts()

    with WorkDir(src=outdir) as od:
        regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
        with WorkDir(src=regdir) as rd:
            rd.mkdir()
            outdir: str = od.abspath()
            regdir: str = rd.abspath()
            regname: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')

    if struct_boundarymask:
        with NiiFile(src=struct_boundarymask, assert_exists=True, validate_nifti=True) as sbm:
            struct_boundarymask: str = sbm.abspath()
    
    if bbr:
        if (struct_boundarymask is None) or (struct_boundarymask == ""):
            if log: log.error("RuntimeError: 'struct_boundarymask' is required for BBR")
            raise RuntimeError("'struct_boundarymask' is required for BBR")

        if isinstance(bbr_slope, float) or isinstance(bbr_slope, int):
            pass
        else:
            if log: log.error(f"TypeError: Input 'bbr_slope': {bbr_slope} is not an 'int' nor 'float'.")
            raise TypeError(f"Input 'bbr_slope': {bbr_slope} is not an 'int' nor 'float'.")
    
    # Define output files
    outputs: Dict[str,str] = {
                                "resampled_image": f"{regname}_img.nii.gz",
                                "affine": f"{regname}_affine.mat",
                                "inv_affine": f"{regname}_invaffine.mat",
                                "init_affine": f"{regname}_init_affine.mat",
                                "resampled_image_init": f"{regname}_init_img.nii.gz"
                             }

    affine: str = outputs.get('affine')
    inv_affine: str = outputs.get('inv_affine')
    src2ref: str = outputs.get('resampled_image')
    init_affine: str = None

    with TmpDir(src=regdir) as tmp:
        tmp.mkdir()

        # Brain extraction
        struct_brain: str = os.path.join(tmp.src, 'struct_brain.nii.gz')
        struct_brain: str = fslmaths(img=struct).mas(img=struct_brainmask).run(out=struct_brain, log=log)

        fmap_brain: str = os.path.join(tmp.src, 'fmap_brain.nii.gz')
        fmap_brain: str = fslmaths(img=fmap_magnitude).mas(img=fmap_brainmask).run(out=fmap_brain, log=log)

        # Perform registration
        kwargs: Dict[str,str] = { 
                                    "outdir": regdir,
                                    "log": log,
                                    "src": fmap_brain,
                                    "ref": struct_brain,
                                    "affine": affine,
                                    "inv_affine": inv_affine,
                                    "src2ref": src2ref,
                                }
        if bbr:
            init_affine: str = outputs.get('init_affine')
            kwargs: Dict[str,str] = {
                                        **kwargs,
                                        "ref_boundarymask": struct_boundarymask,
                                        "bbrslope": bbr_slope,
                                        "bbrtype": BBRType(bbr_type.lower()).name,
                                        "init_affine": init_affine
                                    }

        (affine, 
         inv_affine, 
         src2ref, 
         _) = epireg(**kwargs)

        tmp.rmdir()
    
    return (affine,
            inv_affine,
            src2ref,
            init_affine)


def func_to_sbref(outdir: str,
                  func: str,
                  func_brainmask: str,
                  sbref: str,
                  sbref_brainmask: str,
                  src_space: Optional[str] = None,
                  ref_space: Optional[str] = None,
                  log: Optional[LogFile] = None
                 ) -> Tuple[str,str,str]:
    """Rigid-body register functional 4D EPI to sbref."""
    # Check and validate input structural associated NIFTI images
    with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
        with NiiFile(src=func_brainmask, assert_exists=True, validate_nifti=True) as fb:
            func: str = fn.abspath()
            func_brainmask: str = fb.abspath()
            if not src_space:
                _, src_space, _ = fn.file_parts()
    
    with NiiFile(src=sbref, assert_exists=True, validate_nifti=True) as sb:
        with NiiFile(src=sbref_brainmask, assert_exists=True, validate_nifti=True) as sm:
            sbref: str = sb.abspath()
            sbref_brainmask: str = sm.abspath()
            if not ref_space:
                _, ref_space, _ = sb.file_parts()

    with WorkDir(src=outdir) as od:
        regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
        with WorkDir(src=regdir) as rd:
            if not rd.exists(): rd.mkdir()
            outdir: str = od.abspath()
            regdir: str = rd.abspath()
            regname: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')
    
    # Define outputs
    outputs: Dict[str,str] = {
                                "resampled_image": f"{regname}_img.nii.gz",
                                "affine": f"{regname}_affine.mat",
                                "inv_affine": f"{regname}_invaffine.mat",
                             }

    with TmpDir(src=regdir) as tmp:
        tmp.mkdir()

        func_brain: str = os.path.join(tmp.src, 'func_brain.nii.gz')
        func_brain: str = fslmaths(img=func).mas(func_brainmask).run(out=func_brain, log=log)

        sbref_brain: str = os.path.join(tmp.src, 'sbref_brain.nii.gz')
        sbref_brain: str = fslmaths(img=sbref).mas(sbref_brainmask).run(out=sbref_brain, log=log)

        # Perform registration
        #   1. Brain-extracted
        #   2. Not brain-extracted

        init_affine: str = os.path.join(tmp.src, 'init.mat')

        (init_affine, 
        _, 
        _, 
        _, 
        _) = epireg(outdir=tmp.src,
                    src=func_brain,
                    ref=sbref_brain,
                    affine=init_affine,
                    inv_affine=os.path.join(tmp.src, 'init_inv.mat'),
                    src2ref=os.path.join(tmp.src, 'init.nii.gz'))
        
        (affine, 
        inv_affine, 
        src2ref, 
        _, 
        _) = epireg(outdir=tmp.src,
                    src=func,
                    ref=sbref,
                    init_xfm=init_affine,
                    affine=outputs.get('affine'),
                    inv_affine=outputs.get('inv_affine'),
                    src2ref=outputs.get('resampled_image'))
        tmp.rmdir()
    return src2ref, affine, inv_affine


def sbref_to_struct(outdir: str,
                    sbref: str,
                    sbref_brainmask: str,
                    struct: str,
                    struct_brainmask: str,
                    dc: bool = True,
                    fmap: Optional[str] = None,
                    fmap_brainmask: Optional[str] = None,
                    fmap2struct_xfm: Optional[str] = None,
                    sbref_pedir: Optional[str] = None,
                    sbref_echospacing: Optional[float] = None,
                    struct_boundarymask: Optional[str] = None,
                    bbr: bool = True,
                    bbrslope: float = 0.5,
                    bbrtype: str = "signed",
                    src_space: Optional[str] = None,
                    ref_space: Optional[str] = None,
                    log: Optional[LogFile] = None
                   ) -> Tuple[Union[str,None]]:
    """Register ``sbref`` to structural image with BBR (optional) and distortion correction (optional)."""
    # Check and validate input structural associated NIFTI images
    with NiiFile(src=sbref, assert_exists=True, validate_nifti=True) as sf:
        with NiiFile(src=sbref_brainmask, assert_exists=True, validate_nifti=True) as sb:
            sbref: str = sf.abspath()
            sbref_brainmask: str = sb.abspath()
            if not src_space:
                _, src_space, _ = sf.file_parts()

    with NiiFile(src=struct, assert_exists=True, validate_nifti=True) as st:
        with NiiFile(src=struct_brainmask, assert_exists=True, validate_nifti=True) as sb:
            struct: str = st.abspath()
            struct_brainmask: str = sb.abspath()
            if not ref_space:
                _, ref_space, _ = st.file_parts()
    
    with WorkDir(src=outdir) as od:
        regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
        with WorkDir(src=regdir) as rd:
            outdir: str = od.abspath()
            regdir: str = rd.abspath()
            regname: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')

            dc_warp: str = rd.join(f'{src_space}_dc_warp.nii.gz') 
            dc_img: str = rd.join(f'{src_space}_dc_img.nii.gz')
            dc_brainmask: str = rd.join(f'{src_space}_dc_brainmask.nii.gz')
    
    # Logic tests
    _has_fmap: bool = fmap is not None and fmap_brainmask is not None
    _has_boundary: bool = struct_boundarymask is not None
    _has_sbref_props: bool = sbref_pedir is not None and sbref_echospacing is not None

    if bbr:
        if not _has_boundary:
            raise RuntimeError("'struct_boundarymask' is required for BBR.")
        
        with NiiFile(src=struct_boundarymask, assert_exists=True, validate_nifti=True) as sb:
            struct_boundarymask: str = sb.abspath()
            bbrslope: float = float(bbrslope)
            bbrtype: str = BBRType(bbrtype.lower()).name
    
    if dc:
        if not _has_fmap:
            raise RuntimeError("'fmap' and 'fmap_brainmask' are required for distortion correction.")

        if not _has_sbref_props:
            raise RuntimeError("'sbref_pedir' and 'sbref_echospacing' are required for distortion correction.")

        with NiiFile(src=fmap, assert_exists=True, validate_nifti=True) as fm:
            with NiiFile(src=fmap_brainmask, assert_exists=True, validate_nifti=True) as fb:
                fmap: str = fm.abspath()
                fmap_brainmask: str = fb.abspath()

                sbref_echospacing: float = float(sbref_echospacing)
                sbref_pedir: str = PhaseEncodeDirection(sbref_pedir.lower()).name
    
    # Define outputs
    outputs: Dict[str, str] = {
                                "resampled_image": f"{regname}_img.nii.gz",
                                "affine": f"{regname}_affine.mat",
                                "inv_affine": f"{regname}_invaffine.mat",
                                "init_affine": f'{regname}_init_affine.mat',
                                "resampled_image_init": f"{regname}_init_img.nii.gz",
                                "warp": f'{regname}_warp.nii.gz',
                                "dc_warp": dc_warp,
                                "dc_image": dc_img,
                                "dc_brainmask": dc_brainmask,
                              }

    # Perform registration(s)
    with TmpDir(src=regdir) as tmp:
        sbref_brain: str = tmp.join('sbref_brain.nii.gz')
        struct_brain: str = tmp.join('struct_brain.nii.gz')

        sbref_brain: str = fslmaths(img=sbref).mas(sbref_brainmask).run(out=sbref_brain, log=log)
        struct_brain: str = fslmaths(img=struct).mas(struct_brainmask).run(out=struct_brain, log=log)

        epireg_kwargs: Dict[str,str] = {
            "src": sbref_brain,
            "ref": struct_brain,
            "affine": outputs.get('affine'),
            "inv_affine": outputs.get('inv_affine'),
            "src2ref": outputs.get('resampled_image')
        }

        if bbr:
            epireg_kwargs: Dict[str,str] = {
                **epireg_kwargs,
                "ref_boundarymask": struct_boundarymask,
                "bbrslope": bbrslope,
                "bbrtype": bbrtype,
                "init_affine": outputs.get('init_affine'),
                "src2ref_init": outputs.get('resampled_image_init')
            }
        
        if dc:
            if fmap2struct_xfm is None:
                fmap2struct_xfm = os.path.join(FSLDIR, 'etc/flirtsch/ident.mat')

            with File(src=fmap2struct_xfm, assert_exists=True) as fxm:
                fmap2struct_xfm: str = fxm.abspath()
            
            epireg_kwargs: Dict[str,str] = {
                **epireg_kwargs,
                "fmap": fmap,
                "fmap_brainmask": fmap_brainmask,
                "fmap2ref_xfm": fmap2struct_xfm,
                "src_pedir": sbref_pedir,
                "src_echospacing": sbref_echospacing,
                "warp": outputs.get('warp')
            }
        
        (affine,
        inv_affine,
        src2ref,
        warp,
        _) = epireg(**epireg_kwargs)

        # Create warp to distortion correct func and/or sbref/sbref_brainmask
        if dc:
            dc_warp: str = convertwarp(warp1=warp,
                                       postmat=inv_affine,
                                       ref=sbref,
                                       out=outputs.get('dc_warp'))

            dc_img: str = applywarp(src=sbref,
                                    ref=sbref,
                                    warp=dc_warp,
                                    out=outputs.get('dc_image'),
                                    interp='spline')

            dc_brainmask: str = applywarp(src=sbref_brainmask,
                                          ref=sbref_brainmask,
                                          warp=dc_warp,
                                          out=outputs.get('dc_brainmask'),
                                          interp='trilinear')

            dc_brainmask: str = fslmaths(img=dc_brainmask).thr(0.5).bin().run(out=dc_brainmask, log=log)
        else:
            dc_warp: str = None
            dc_img: str = None
            dc_brainmask: str = None

    return (affine, 
            inv_affine, 
            src2ref, 
            warp,
            dc_warp, 
            dc_img, 
            dc_brainmask)


def func_to_struct_composite(outdir: str,
                             func: str,
                             struct: str,
                             func2sbref_affine: str,
                             sbref2struct_affine: str,
                             sbref2struct_warp: Optional[str] = None,
                             src_space: Optional[str] = None,
                             ref_space: Optional[str] = None,
                             log: Optional[LogFile] = None
                            ) -> Tuple[Union[str,None]]:
    """Create composite transform func -> sbref -> struct.
    Create distortion correction warp for func and sbref.
    """
    with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
        with NiiFile(src=struct, assert_exists=True, validate_nifti=True) as st:
            func: str = fn.abspath()
            struct: str = st.abspath()

            if not src_space: _, src_space, _ = fn.file_parts()
            if not ref_space: _, ref_space, _ = st.file_parts()
    
    with File(src=func2sbref_affine, assert_exists=True) as fsa:
        with File(src=sbref2struct_affine, assert_exists=True) as ssa:
            func2sbref_affine: str = fsa.abspath()
            sbref2struct_affine: str = ssa.abspath()

    with WorkDir(src=outdir) as od:
        regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
        with WorkDir(src=regdir) as rd:
            outdir: str = od.abspath()
            regdir: str = rd.abspath()
            regname: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')

            dc_warp: str = rd.join(f'{src_space}_dc_warp.nii.gz') 
            dc_img: str = rd.join(f'{src_space}_dc_img.nii.gz')

    _has_warp: bool = sbref2struct_warp is not None

    # Define outputs
    outputs: Dict[str,str] = {
                                "resampled_image": f"{regname}_img.nii.gz",
                                "affine": f"{regname}_affine.mat",
                                "inv_affine": f"{regname}_invaffine.mat",
                                "warp": f"{regname}_warp.nii.gz",
                                "inv_warp": f"{regname}_invwarp.nii.gz",
                                "dc_warp": dc_warp,
                                "dc_image": dc_img
                             }

    # Func (undistorted) -> struct :: affine
    affine: str = concatxfm(inmat1=func2sbref_affine, 
                            inmat2=sbref2struct_affine, 
                            outmat=outputs.get('affine'),
                            log=log)

    # Struct -> func (undistorted) :: affine
    inv_affine: str = invxfm(inmat=affine, 
                             outmat=outputs.get('inv_affine'), 
                             log=log)

    if _has_warp:
        with NiiFile(src=sbref2struct_warp, assert_exists=True, validate_nifti=True) as sbs:
            sbref2struct_warp: str = sbs.abspath()
        
        # Func (distorted) -> struct :: warp
        warp: str = convertwarp(premat=func2sbref_affine, 
                                warp1=sbref2struct_warp, 
                                out=outputs.get('warp'), 
                                ref=struct, 
                                log=log)

        # Struct -> func (distorted) :: warp
        inv_warp: str = invwarp(inwarp=warp, 
                                ref=func, 
                                outwarp=outputs.get('inv_warp'),
                                log=log)
        
        # Create warp to distortion correct func
        dc_warp: str = convertwarp(warp1=warp, 
                                   postmat=inv_affine, 
                                   ref=func, 
                                   out=outputs.get('dc_warp'), 
                                   log=log)

        dc_img: str = applywarp(src=func, 
                                ref=func, 
                                warp=dc_warp, 
                                out=outputs.get('dc_image'), 
                                log=log)

        # Resample func -> structural space
        resamp_img: str = applywarp(src=func, 
                                    ref=struct, 
                                    warp=warp, 
                                    out=outputs.get('resampled_image'), 
                                    interp='spline', 
                                    log=log)
    else:
        # Resample func -> structural space
        resamp_img: str = applywarp(src=func, 
                                    ref=struct, 
                                    premat=affine,
                                    out=outputs.get('resampled_image'), 
                                    interp='spline', 
                                    log=log)
        warp: str = None
        inv_warp: str = None
        dc_warp: str = None
        dc_img: str = None

    return (affine,
            inv_affine,
            warp,
            inv_warp,
            dc_warp,
            dc_img,
            resamp_img)


def fmap_to_func_composite(outdir: str,
                           fmap: str,
                           func: str,
                           fmap2struct_affine: str,
                           func2struct_invaffine: str,
                           src_space: Optional[str] = None,
                           ref_space: Optional[str] = None,
                           log: Optional[LogFile] = None
                          ) -> Tuple[str,str]:
    """Create composite transform fieldmap -> struct -> sbref -> func.
    """
    with NiiFile(src=fmap, assert_exists=True, validate_nifti=True) as fm:
        with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
            fmap: str = fm.abspath()
            func: str = fn.abspath()

            if not src_space: _, src_space, _ = fm.file_parts()
            if not ref_space: _, ref_space, _ = fn.file_parts()
    
    with File(src=fmap2struct_affine, assert_exists=True) as fsa:
        with File(src=func2struct_invaffine, assert_exists=True) as fsi:
            fmap2struct_affine: str = fsa.abspath()
            func2struct_invaffine: str = fsi.abspath()
    
    with WorkDir(src=outdir) as od:
        regdir: str = od.join('reg',f'{src_space}_to_{ref_space}')
        with WorkDir(src=regdir) as rd:
            outdir: str = od.abspath()
            regdir: str = rd.abspath()
            regname: str = rd.join(f'{src_space}_to_{ref_space}')
    
    # Define outputs
    outputs: Dict[str,str] = {
                                "resampled_image": f"{regname}_img.nii.gz",
                                "affine": f"{regname}_affine.mat"
                             }

    # Fmap -> func
    affine: str = concatxfm(inmat1=fmap2struct_affine, 
                            inmat2=func2struct_invaffine, 
                            outmat=outputs.get('affine'), 
                            log=log)
    
    resamp_img: str = applyxfm(src=fmap,
                               ref=func,
                               mat=affine,
                               interp='spline',
                               out=outputs.get('resampled_image'),
                               log=log)
    return affine, resamp_img


def _select_atlas(age: Union[int,str],
                  atlasdir: Optional[str] = None
                 ) -> Dict[str,str]:
    """Helper function designed to gather the relevant atlas/template
    files for either the dHCP volumetric atlas or the UNC AAL atlas.
    """
    # Set atlas age
    try:
        age: int = int(age)
        if age == 1:
            age: str = "1yr"
            atlas_name: str = "UNCAAL"
        elif age == 2:
            age: str = "2yr"
            atlas_name: str = "UNCAAL"
        else:
            # Cap or trim age if in PMA (for dHCP)
            age: int = 28 if age < 28 else age
            age: int = 44 if age > 44 else age
            atlas_name: str = "dHCPtemplate"
    except ValueError:
        if (age.lower() == 'neo') or (age.lower() == 'neonate'):
            age: str = 'neo'
            atlas_name: str = "UNCAAL"
        else:
            raise RuntimeError(f"The input for age is insufficient: {age}. Use 'neo' or 'neonate' or the PMA in weeks.")
    
    # Select atlas directory
    if atlasdir:
        pass
    elif ((age == 'neo') or (age == '1yr') or (age == '2yr')):
        atlasdir: str = ' '.join(map(str, glob.glob(os.path.join(ATLASDIR,"UNC*2020*"))))
    elif ((28 <= age) and (age <= 44)):
        atlasdir: str = ' '.join(map(str, glob.glob(os.path.join(ATLASDIR,"dHCPatlas"))))
    else:
        raise FileNotFoundError("'atlasdir' does exists or was not specified.")

    
    with WorkDir(src=atlasdir) as ad:
        atlasdir: str = ad.abspath()
        templatedir: str = ad.join("templates")

        with WorkDir(src=templatedir) as td:
            templatedir: str = td.abspath()
            _cwd: str = os.getcwd()
            os.chdir(templatedir)

            t1: str = ' '.join(map(str, glob.glob(f"*{age}*T1*.nii*")))
            t2: str = ' '.join(map(str, glob.glob(f"*{age}*T2*.nii*")))
            seg: str = ' '.join(map(str, glob.glob(f"*{age}*dseg*.nii*")))
            
            if t2 == "":
                t2: str = ' '.join(map(str, glob.glob(f"*{age}.nii*")))
            
            if seg == "":
                seg: str = ' '.join(map(str, glob.glob(f"*{age}*-seg.nii*")))
            
            gm_prob: str = ' '.join(map(str, glob.glob(f"*{age}*gm*.nii*")))
            wm_prob: str = ' '.join(map(str, glob.glob(f"*{age}*wm*.nii*")))
            probseg: str = ' '.join(map(str, glob.glob(f"*{age}*probseg*.nii*")))

            if probseg == "":
                probseg: str = ' '.join(map(str, glob.glob(f"*{age}*avgseg*.nii*")))

            # Set variables to None if no atlas file is found
            if t1 == "": 
                t1: str = None
            else:
                t1: str = td.join(t1)

            if gm_prob == "": 
                gm_prob: str = None
            else:
                gm_prob: str = td.join(gm_prob)

            if wm_prob == "": 
                wm_prob: str = None
            else:
                wm_prob: str = td.join(wm_prob)

            if probseg == "": 
                probseg: str = None
            else:
                probseg: str = td.join(probseg)

            if seg == "": 
                seg: str = None
            else:
                seg: str = td.join(seg)

            if t2 == "":
                raise FileNotFoundError(f"Unable to find T2w image template for the specified age: {age} and atlas directory: {td.abspath()}")
            
            if 'dHCPatlas' in atlasdir:
                age: str = f"{age}wks"
            
            atlasdict: Dict[str,str] = {
                                            "template_T1": t1,
                                            "template_T2": t2,
                                            "gm_prob": gm_prob,
                                            "wm_prob": wm_prob,
                                            "seg": seg,
                                            "prob": probseg,
                                            "atlas_name": atlas_name,
                                            "atlas_age": age
                                        }
            os.chdir(_cwd)

    return atlasdict


def template_to_struct(outdir: str,
                       age: Union[int,str],
                       struct_brainmask: str,
                       struct_T1w: Optional[str] = None,
                       struct_T2w: Optional[str] = None,
                       struct_gmprob: Optional[str] = None,
                       struct_wmprob: Optional[str] = None,
                       quick: bool = False,
                       src_space: Optional[str] = None,
                       ref_space: Optional[str] = None,
                       atlasdir: Optional[str] = None,
                       log: Optional[LogFile] = None
                      ) -> Tuple[str]:
    """Register (age-matched-) template -> struct.
    """
    atlas: Dict[str,str] = _select_atlas(age=age, atlasdir=atlasdir)

    if (struct_T1w is None) and (struct_T2w is None):
        raise RuntimeError('T2w and/or T1w is required')
    
    with NiiFile(src=struct_brainmask, assert_exists=True, validate_nifti=True) as sbn:
        struct_brainmask: str = sbn.abspath()
    
    if not src_space:
        src_space: str = f"{atlas.get('atlas_name')}-{atlas.get('atlas_age')}"
    
    if not ref_space:
        ref_space: str = 'struct'
    
    with WorkDir(src=outdir) as od:
        regdir: str = od.join('reg',f'{src_space}_to_{ref_space}')
        with WorkDir(src=regdir) as rd:
            outdir: str = od.abspath()
            regdir: str = rd.abspath()
            regname: str = rd.join(f'{src_space}_to_{ref_space}')

    # Define outputs
    outputs: Dict[str,str] = {
                                "resampled_image": f"{regname}_img.nii.gz",
                                "affine": f"{regname}_affine.mat",
                                "warp": f"{regname}_warp.nii.gz",
                                "inv_warp": f"{regname}_invwarp.nii.gz"
                             }

    # Build structural and template arguments
    struct: List[str] = []
    template: List[str] = []

    if struct_T2w and atlas.get('template_T2'):
        with NiiFile(src=struct_T2w, assert_exists=True, validate_nifti=True) as st2:
            with NiiFile(src=atlas.get('template_T2'), assert_exists=True, validate_nifti=True) as tt2:
                struct_T2w: str = st2.abspath()
                template_T2w: str = tt2.abspath()

                struct += [struct_T2w]
                template += [template_T2w]

    if struct_T1w and atlas.get('template_T1'):
        with NiiFile(src=struct_T1w, assert_exists=True, validate_nifti=True) as st1:
            with NiiFile(src=atlas.get('template_T1'), assert_exists=True, validate_nifti=True) as tt1:
                struct_T1w: str = st1.abspath()
                template_T1w: str = tt1.abspath()

                struct += [struct_T1w]
                template += [template_T1w]
    
    if struct_gmprob:
        with NiiFile(src=struct_gmprob, assert_exists=True, validate_nifti=True) as stg:
            struct_gmprob: str = stg.abspath()
            struct += [struct_gmprob]

        if 'dHCP' in atlas.get('prob'):
            with NiiFile(src=atlas.get('prob'), assert_exists=True, validate_nifti=True) as tgmpb:
                tprob: str = tgmpb.abspath()
                gmpb: str = os.path.join(regdir,'template_gmprob.nii.gz')
                gmpb: str = fslroi(img=tprob,
                                   out=gmpb,
                                   tmin=1,
                                   tsize=1,
                                   log=log)
                gmpb: str = fslmaths(img=gmpb).div(100).run(out=gmpb, log=log)
        else:
            with NiiFile(src=atlas.get('gm_prob'), assert_exists=True, validate_nifti=True) as tgmpb:
                gmpb: str = tgmpb.abspath()
        
        template += [gmpb]
    
    if struct_wmprob:
        with NiiFile(src=struct_wmprob, assert_exists=True, validate_nifti=True) as stw:
            struct_wmprob: str = stw.abspath()
            struct += [struct_wmprob]

        if 'dHCP' in atlas.get('prob'):
            with NiiFile(src=atlas.get('prob'), assert_exists=True, validate_nifti=True) as twmpb:
                tprob: str = twmpb.abspath()
                wmpb: str = os.path.join(regdir,'template_wmprob.nii.gz')
                wmpb: str = fslroi(img=tprob,
                                   out=wmpb,
                                   tmin=1,
                                   tsize=1,
                                   log=log)
                wmpb: str = fslmaths(img=wmpb).div(100).run(out=wmpb, log=log)
        else:
            with NiiFile(src=atlas.get('wm_prob'), assert_exists=True, validate_nifti=True) as twmpb:
                wmpb: str = twmpb.abspath()
        
        template += [wmpb]
    
    with NiiFile(src=outputs.get('warp')) as wp:
        antsout_dir: str = wp.rm_ext()

    (warp,
    inv_warp,
    affine,
    src2ref) = nonlinear_reg(src=template,
                            ref=struct,
                            quick=quick,
                            ref_brainmask=struct_brainmask,
                            antsout=antsout_dir + + '.ants/output_',
                            affine=outputs.get('affine'),
                            warp=outputs.get('warp'),
                            inv_warp=outputs.get('inv_warp'),
                            src2ref=outputs.get('resampled_image'),
                            log=log)
    
    inv_warp: str = invwarp(inwarp=warp,
                            ref=template[0],
                            outwarp=inv_warp,
                            log=log)

    return (warp,
            inv_warp,
            affine,
            src2ref)


def struct_to_template_composite():
    """doc-string"""
    pass


def func_to_template_composite():
    """doc-string"""
    pass


def epireg(outdir: str,
           src: str,
           ref: str,
           src_space: Optional[str] = None,
           ref_space: Optional[str] = None,
           src_brainmask: Optional[str] = None,
           ref_brainmask: Optional[str] = None,
           src_pedir: Optional[str] = None,
           src_echospacing: Optional[float] = 0.1,
           ref_boundarymask: Optional[str] = None,
           fmap: Optional[str] = None,
           fmap_brainmask: Optional[str] = None,
           fmap2ref_xfm: Optional[str] = None,
           init_xfm: Optional[str] = None,
           bbrslope: float = -0.5,
           bbrtype: str = "signed",
           basename: Optional[str] = None,
           log: Optional[LogFile] = None,
           **kwargs
          ) -> Tuple[str,str,str,str,str]:
    """Linear (FLIRT-based) registration of EPI to ref, with optional BBR and distortion correction.
    """
    with NiiFile(src=src, assert_exists=True, validate_nifti=True) as sr:
        with NiiFile(src=ref, assert_exists=True, validate_nifti=True) as rf:
            src: str = sr.abspath()
            ref: str = rf.abspath()
            if not src_space: _, src_space, _ = sr.file_parts()
            if not ref_space: _, ref_space, _ = rf.file_parts()
    
    if log: log.log(f"Performing epireg: Register {src} to {ref} ")

    if basename and src_space and ref_space:
        with WorkDir(src=outdir) as od:
            regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
            basename: str = os.path.join(regdir,basename)
            dc_warp: str = os.path.join(regdir, f'{src_space}_dc_warp.nii.gz')
    elif basename:
        with WorkDir(src=outdir) as od:
            if not od.exists(): od.mkdir()
            regdir: str = od.abspath()
            basename: str = os.path.join(regdir,basename)
            dc_warp: str = os.path.join(regdir, f'{src_space}_dc_warp.nii.gz')
    else: 
        with WorkDir(src=outdir) as od:
            regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
            with WorkDir(src=regdir) as rd:
                if not rd.exists(): 
                    rd.mkdir()
                outdir: str = od.abspath()
                regdir: str = rd.abspath()
                basename: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')
                dc_warp: str = os.path.join(regdir, f'{src_space}_dc_warp.nii.gz')
    
    # Define outputs
    outputs: Dict[str,str] = {
                                "affine": f'{basename}_affine.mat',
                                "init_affine": f'{basename}_init_affine.mat',
                                "src2ref_init": f'{basename}_init_img.nii.gz',
                                "warp": f'{basename}_warp.nii.gz',
                                "inv_affine": f'{basename}_invaffine.mat',
                                "inv_warp": f'{basename}_invwarp.nii.gz',
                                "src2ref": f'{basename}_img.nii.gz',
                                "dc_warp": dc_warp
                             }

    # Update outputs/templates from kwargs
    for k,v in kwargs.items():
        if k not in outputs:
            raise ValueError(f'Unknown output template {k}')
        outputs[k] = v
    
    # Logic tests
    _has_fmap: bool = (fmap is not None) and (fmap_brainmask is not None) and (fmap2ref_xfm is not None)
    _do_bbr: bool = (ref_boundarymask is not None)
    _do_dc: bool = _has_fmap and (src_pedir is not None) and (src_echospacing is not None)

    if log: log.debug(f'Has fieldmap: {_has_fmap}')
    if log: log.debug(f'Perform BBR: {_do_bbr}')
    if log: log.debug(f'Perform distortion correction: {_do_dc}')

    # Check fieldmap inputs
    if _has_fmap:
        with NiiFile(src=fmap, assert_exists=True, validate_nifti=True) as fmp:
            with NiiFile(src=fmap_brainmask, assert_exists=True, validate_nifti=True) as fmb:
                fmap: str = fmp.abspath()
                fmap_brainmask: str = fmb.abspath()

    # Get phase-encode direction for FLIRT
    if src_pedir in ['PA', 'AP', 'IS', 'SI', 'RL', 'LR']:
        src_pedir: str = get_axis(src, src_pedir)[0]
    
    pedict: Dict[str,Tuple[int,str]] = phase_encode_dict()
    if src_pedir:
        src_pedir: Tuple[int,str] = pedict[src_pedir]
    else:
        src_pedir: Tuple[None] = (None, None)
    
    # TODO: Incorporate sigloss
    sigloss2ref: str = None

    with TmpDir(src=regdir) as tmp:
        tmp.mkdir()

        affine: str = None
        inv_affine: str = None
        src2ref: str = None
        warp: str = None
        inv_warp: str = None

        if src_brainmask:
            src_brain: str = os.path.join(tmp.src,'src_brain.nii.gz')
            src_brain: str = fslmaths(img=src).mas(src_brainmask).run(out=src_brain, log=log)
        
        if ref_brainmask:
            ref_brain: str = os.path.join(tmp.src, 'ref_brain.nii.gz')
            ref_brain: str = fslmaths(img=ref).mas(ref_brain).run(out=ref_brain, log=log)
        
        # Initial registration src -> ref (for BBR)
        init_affine: str = outputs.get('init_affine')
        init_img: str = outputs.get('src2ref_init')

        if _do_bbr:
            bbrtype: str = BBRType(bbrtype.lower()).name
            bbrslope: float = float(bbrslope)

            _: Tuple[str] = flirt(
                                src=src if src_brainmask is None else src_brain,
                                ref=ref if ref_brainmask is None else ref_brain,
                                omat=init_affine,
                                out=init_img,
                                dof=6,
                                interp='spline',
                                usesqform=True,
                                nosearch=True,
                                log=log)

        # Main registration: src -> ref
        affine: str = outputs.get('affine')

        flirt_kwargs: Dict[str,str] = {
            "src": src if src_brainmask is None else src_brain,
            "ref": ref if ref_brainmask is None else ref_brain,
            "dof": 6,
            "cost": 'corratio',
            "omat": affine,
            "refweight": sigloss2ref,
            "usesqform": True,
            "nosearch": True,
            "interp": 'spline',
            "log": log,
        }

        if init_xfm:
            with File(src=init_xfm, assert_exists=True) as xfm:
                flirt_kwargs['init'] = xfm.abspath()
        
        if _do_bbr:
            schedule: str = os.path.join(FSLDIR, 'etc', 'flirtsch', 'bbr.sch')

            flirt_kwargs: Dict[str,str] = {
                **flirt_kwargs,
                "cost": 'bbr',
                "wmseg": ref_boundarymask,
                "bbrslope": bbrslope,
                "bbrtype": bbrtype,
                "schedule": schedule,
                "init": init_affine
            }
        
        # Resample fmap and fmap_brainmask -> ref space
        if _do_dc:
            fmap_to_ref: str = os.path.join(tmp.src, 'fmap_to_ref.nii.gz')
            fmap_brainmask_to_ref: str = os.path.join(tmp.src, 'brainmask_to_ref.nii.gz')

            # fmap_to_ref
            fmap_to_ref: str = applyxfm(src=fmap,
                                        ref=ref,
                                        mat=fmap2ref_xfm,
                                        out=fmap_to_ref,
                                        interp='spline',
                                        log=log)
            # fmap_brainmask_to_ref
            fmap_brainmask_to_ref: str = applyxfm(src=fmap_brainmask,
                                                  ref=ref,
                                                  mat=fmap2ref_xfm,
                                                  out=fmap_brainmask_to_ref,
                                                  interp='nearestneighbour',
                                                  log=log)
            
            flirt_kwargs: Dict[str,str] = {
                **flirt_kwargs,
                "echospacing": src_echospacing,
                "pedir": src_pedir[0],
                "fieldmap": fmap_to_ref,
                "fieldmapmask": fmap_brainmask_to_ref,
            }
        
        # Run FLIRT: Perform linear registration
        _: str = flirt(**flirt_kwargs)

        # Invert matrix
        inv_affine: str = outputs.get('inv_affine')
        inv_affine: str = invxfm(inmat=affine, outmat=inv_affine, log=log)

        # If not performing distortion correction return current transforms
        if not _do_dc:
            src2ref: str = outputs.get('src2ref')

            src2ref: str = applyxfm(src=src,
                                    ref=ref,
                                    mat=affine,
                                    out=src2ref,
                                    interp='spline',
                                    log=log)
            tmp.rmdir()
            return affine, inv_affine, src2ref, warp, inv_warp
        
        # Perform distortion correction
        warp: str = outputs.get('warp')

        ## Transform fieldmap_brainmask from ref to src space
        fmap2src_mask: str = os.path.join(tmp.src, "fmap_to_src_brainmask.nii.gz")
        fmap2src_mask: str = applyxfm(src=fmap_brainmask_to_ref,
                                      ref=src,
                                      mat=inv_affine,
                                      out=fmap2src_mask,
                                      interp='nearestneighbour',
                                      log=log)

        ## Transform fieldmap from ref to src space
        fmap2src: str = os.path.join(tmp.src, "fmap_to_src_phase.nii.gz")
        fmap2src: str = applyxfm(src=fmap_to_ref,
                                 ref=src,
                                 mat=inv_affine,
                                 out=fmap2src,
                                 interp='spline',
                                 log=log)
        fmap2src: str = fslmaths(img=fmap2src).mul(fmap2src_mask).run(out=fmap2src, log=log) # Re-mask

        ## Calculate shiftmap and combine transform (src2ref)
        fmap2src_shift: str = os.path.join(tmp.src, "fmap_to_src_shift.nii.gz")

        _: Tuple[str] = fugue(log=log, 
                              loadfmap=fmap2src,
                              mask=fmap2src_mask,
                              saveshift=fmap2src_shift,
                              unmaskshift=True,
                              dwell=src_echospacing,
                              unwarpdir=src_pedir[1])

        _: str = convertwarp(out=warp,
                             ref=ref,
                             postmat=affine,
                             shiftmap=fmap2src_shift,
                             shiftdir=src_pedir[1],
                             abs=True,
                             absout=True,
                             log=log)
        ## Inverse warp
        inv_warp: str = outputs.get('inv_warp')
        inv_warp: str = inv_warp(inwarp=warp,
                                 ref=src,
                                 outwarp=inv_warp,
                                 log=log)

        src2ref: str = outputs.get('src2ref')
        src2ref: str = applywarp(src=src,
                                 ref=ref,
                                 warp=warp,
                                 out=src2ref,
                                 interp='spline',
                                 log=log)
        tmp.rmdir()
    return affine, inv_affine, src2ref, warp, inv_warp


def nonlinear_reg(outdir: str,
                  src: str,
                  ref: str,
                  dim: int = 3,
                  quick: bool = False,
                  ref_brainmask: Optional[str] = None,
                  src_space: Optional[str] = None,
                  ref_space: Optional[str] = None,
                  basename: Optional[str] = None,
                  log: Optional[LogFile] = None,
                  **kwargs: Dict[str,str]
                 ) -> Tuple[str,str,str,str]:
    """Perform ``ANTs`` deformable registration."""
    if isinstance(src, list):
        pass
    elif isinstance(src, str):
        src: List[str] = [ src ]
    else:
        raise TypeError(f"Input for 'src' is not a list or string: {src}")
    
    if isinstance(ref, list):
        pass
    elif isinstance(ref, str):
        ref: List[str] = [ ref ]
    else:
        raise TypeError(f"Input for 'ref' is not a list or string: {ref}")
    
    if basename and src_space and ref_space:
        with WorkDir(src=outdir) as od:
            regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')
            basename: str = os.path.join(regdir,basename)
    elif basename:
        with WorkDir(src=outdir) as od:
            if not od.exists(): od.mkdir()
            regdir: str = od.abspath()
            basename: str = os.path.join(regdir,basename)
    else: 
        with NiiFile(src=src[0]) as sr:
            with NiiFile(src=ref[0]) as rf:
                with WorkDir(src=outdir) as od:

                    if not src_space: _, src_space, _ = sr.file_parts()
                    if not ref_space: _, ref_space, _ = rf.file_parts()
                    regdir: str = os.path.join(od.src,'reg',f'{src_space}_to_{ref_space}')

                    with WorkDir(src=regdir) as rd:
                        outdir: str = od.abspath()
                        regdir: str = rd.abspath()
                        basename: str = os.path.join(regdir,f'{src_space}_to_{ref_space}')
    
    # Define outputs
    outputs: Dict[str,str] = {
                                "antsout": f'{basename}.ants/output_',
                                "affine": f'{basename}_affine.mat',
                                "warp": f'{basename}_warp.nii.gz',
                                "inv_warp": f'{basename}_invwarp.nii.gz',
                                "inv_affine": f'{basename}_invaffine.mat',
                                "src2ref": f'{basename}_img.nii.gz'
                             }

    # Update outputs/templates from kwargs
    for k,v in kwargs.items():
        if k not in outputs:
            raise ValueError(f'Unknown output template {k}')
        outputs[k] = v
    
    # Perform ANTs registration
    antsRegistrationSyN(fixed=ref,
                        moving=src,
                        dim=dim,
                        out_prefix=outputs.get('antsout'),
                        quick=quick,
                        fixed_mask=ref_brainmask,
                        log=log)

    # Convert ANTs transformation matrix to FSL
    cmd: Command = Command("c3d_affine_tool")

    cmd.opt("-ref"); cmd.opt(f"{ref[0]}") 
    cmd.opt("-src"); cmd.opt(f"{src[0]}") 
    cmd.opt("-itk"); cmd.opt(f"{outputs.get('antsout')}0GenericAffine.mat")
    cmd.opt("-ras2fsl")
    cmd.opt("-o"); cmd.opt(f"{outputs.get('affine')}")
    cmd.run(log=log)

    affine: str = outputs.get('affine')

    # Convert ANTs warp to FSL
    with TmpDir(src=regdir) as tmp:
        cmd: Command = Command("c3d")

        cmd.opt("-mcs"); cmd.opt(f"{outputs.get('antsout')}1Warp.nii.gz")
        cmd.opt("-oo"); cmd.opt(f"{tmp.src}/wx.nii.gz")
        cmd.opt(f"{tmp.src}/wy.nii.gz"); cmd.opt(f"{tmp.src}/wz.nii.gz")
        cmd.run(log=log)

        tmp_invwarp: str = fslmaths(img=f"{tmp.src}/wy.nii.gz").mul(input=int(-1)).run(out=f"{tmp.src}/i_wy.nii.gz", log=log)

        warp: str = fslmerge(f"{outputs.get('warp')}",
                             "t", 
                             None, 
                             log, 
                             f"{tmp.src}/wx.nii.gz",  
                             tmp_invwarp,  
                             f"{tmp.src}/wz.nii.gz")
    
    warp: str = convertwarp(ref=ref[0],
                            premat=affine,
                            warp1=warp,
                            out=warp,
                            log=log)

    inv_warp: str = invwarp(inwarp=warp,
                            outwarp=outputs.get('inv_warp'),
                            ref=src[0],
                            log=log)

    src2ref: str = applywarp(src=src[0],
                             ref=ref[0],
                             warp=warp,
                             out=outputs.get('src2ref'),
                             interp='spline',
                             log=log)
    return (warp,
            inv_warp,
            affine,
            src2ref)


def antsRegistrationSyN(fixed: Union[List[str],str],
                        moving: str,
                        dim: int,
                        out_prefix: str,
                        fixed_mask: Optional[str] = None,
                        quick: bool = False,
                        log: Optional[LogFile] = None
                       ) -> None:
    """Python wrapper for ``ANTs`` ``antsRegistrationSyN``."""
    if isinstance(fixed, list):
        pass
    elif isinstance(fixed, str):
        fixed: List[str] = [ fixed ]
    else:
        raise TypeError(f"Input for 'fixed' is not a list or string: {fixed}")
    
    if isinstance(moving, list):
        pass
    elif isinstance(moving, str):
        moving: List[str] = [ moving ]
    else:
        raise TypeError(f"Input for 'moving' is not a list or string: {moving}")
    
    assert len(fixed) == len(moving), 'The same number of fixed and moving images must be provided'

    if quick:
        cmd: Command = Command("antsRegistrationSyNQuick.sh")
    else:
        cmd: Command = Command("antsRegistrationSyN.sh")
    
    # Get output directory name and create it
    with File(src=out_prefix) as otp:
        outdir, _, _ = otp.file_parts()
        with WorkDir(src=outdir) as _:
            pass
    
    cmd.opt("-d"); cmd.opt(f"{dim}")

    for fixed_img, moving_img in zip(fixed, moving):
        cmd.opt("-f"); cmd.opt(f"{fixed_img}")
        cmd.opt("-m"); cmd.opt(f"{moving_img}")
    
    cmd.opt("-o"), cmd.opt(f"{out_prefix}")
    cmd.opt("-t"); cmd.opt("s")
    cmd.opt("-j"); cmd.opt("1")

    if fixed_mask:
        cmd.opt("-x"); cmd.opt(f"{fixed_mask}")
    
    cmd.run(log=log)
    
    return None
