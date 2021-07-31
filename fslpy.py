# -*- coding: utf-8 -*-
"""Wrappers functions for FSL executable binaries.
"""
import glob
import os

from math import (
    e as EULER_CONST,
    log,
    sqrt
)

from typing import (
    List,
    Optional,
    Text,
    Tuple,
    Union
)

from util import (
    Command,
    DependencyError,
    File,
    NiiFile,
    LogFile
)

from enums import (
    ECModelFLM,
    ECModelSLM,
    MergeDim,
    ECInterp,
    ECresamp,
    ECOLType,
    FSLDataType,
    RegInterp
)

class FSLError(Exception):
    """Exception intended to be raised for FSL specific binaries and related wrapper functions."""
    pass

# TODO: Add verbose options to all wrapper functions

def fnirt():
    """work"""
    pass

def eddy(img: str,
         out: str,
         mask: str,
         idx: str,
         acqp: str,
         bvecs: str,
         bvals: str,
         slspec: Optional[str] = "",
         json_file: Optional[str] = "",
         mporder: int = 0,
         s2v_lambda: int = 1,
         topup: Optional[str] = "",
         field: Optional[str] = "",
         field_mat: Optional[str] = "",
         flm: str = "quadratic",
         slm: str = "none",
         fwhm: Union[int,str] = 0,
         s2v_fwhm: Union[int,str] = 0,
         niter: Union[int,str] = 5,
         s2v_niter: int = 5,
         cnr_maps: bool = False,
         residuals: bool = False,
         fep: bool = False,
         interp: str = "spline",
         s2v_interp: str = "trilinear",
         resamp: str = "jac",
         nvoxhp: int = 1000,
         initrand: int = 0,
         ff: float = 10.0,
         repol: bool = False,
         ol_nstd: int = 4,
         ol_nvox: int = 250,
         ol_type: str = "sw",
         ol_pos: bool = False,
         ol_sqr: bool = False,
         estimate_move_by_susceptibility: bool = False,
         mbs_niter: int = 10,
         mbs_lambda: int = 10,
         mbs_ksp: int = 10,
         dont_sep_offs_move: bool = False,
         dont_peas: bool = False,
         data_is_shelled: bool = False,
         b0_only: bool = False,
         dont_mask_output: bool = False,
         verbose: bool = False,
         very_verbose: bool = False,
         log: Optional[LogFile] = None
        ) -> Tuple[str]:
    """work"""
    # TODO:
    #   * Get the corresponding output files for 
    #       each eddy option.
    eddy_cmds: List[str] = [
        "eddy_cuda",
        "eddy_cuda7.5",
        "eddy_cuda8.0",
        "eddy_cuda9.1",
        "eddy_openmp",
        "eddy"
    ]

    # Select for most optimal implementation of eddy.
    #   Look for GPU implementations first, followed 
    #   by openmp, then eddy.
    for eddy_cmd in eddy_cmds:
        try:
            eddy_proc: Command = Command(eddy_cmd)
            
            (return_code, 
             _, 
             _) = eddy_proc.run(log=log)
            
            # Check if the dependency is met and if the return code
            #   of an empty command options list returns 1.
            if (eddy_proc.check_dependency()) and (return_code == 1):
                break
        except (DependencyError, FileNotFoundError):
            continue
    
    if out.endswith('.nii.gz'):
        out: str = out[:-7]
    elif out.endswith('.nii'):
        out: str = out[:-4]
    
    # Required eddy current correction options
    img:    NiiFile = NiiFile(file=img, assert_exists=True, validate_nifti=True)
    out:    NiiFile = NiiFile(file=out, assert_exists=False, validate_nifti=False)
    mask:   NiiFile = NiiFile(file=mask, assert_exists=True, validate_nifti=True)
    idx:    File = File(file=idx, ext=".idx", assert_exists=True)
    acqp:   File = File(file=acqp, ext=".acqp", assert_exists=True)
    bvecs:  File = File(file=bvecs, ext='.bvec', assert_exists=True)
    bvals:  File = File(file=bvals, ext='.bval', assert_exists=True)

    eddy_proc.opt(f"--imain={img.file}")
    eddy_proc.opt(f"--mask={mask.file}")
    eddy_proc.opt(f"--index={idx.file}")
    eddy_proc.opt(f"--acqp={acqp.file}")
    eddy_proc.opt(f"--bvecs={bvecs.file}")
    eddy_proc.opt(f"--bvals={bvals.file}")
    eddy_proc.opt(f"--out={out.file}")

    # Conventional eddy current correction options
    if topup:
        topup: NiiFile = NiiFile(file=topup, assert_exists=True, validate_nifti=True)
        eddy_proc.opt(f"--topup={topup.rm_ext()}")
    
    if field:
        field: NiiFile = NiiFile(file=field, assert_exists=True, validate_nifti=True)
        eddy_proc.opt(f"--field={field.file}")
    
    if field_mat:
        field_mat: File = File(file=field_mat, assert_exists=True)
        eddy_proc.opt(f"--field_mat={field_mat.file}")
    
    if flm:
        flm: str = ECModelFLM(flm).name
        eddy_proc.opt(f"--flm={flm}")
    
    if slm:
        slm: str = ECModelSLM(slm).name
        eddy_proc.opt(f"--slm={slm}")
    
    if fwhm:
        eddy_proc.opt(f"--fwhm={fwhm}")
    
    if niter:
        eddy_proc.opt(f"--niter={niter}")
    
    if cnr_maps:
        eddy_proc.opt("--cnr_maps")
    
    if residuals:
        eddy_proc.opt("--residuals")
    
    if fep:
        eddy_proc.opt("--fep")
    
    if interp:
        interp: str = ECInterp(interp).name
        eddy_proc.opt(f"--interp={interp}")
    
    if resamp:
        resamp: str = ECresamp(resamp).name
        eddy_proc.opt(f"--resamp={resamp}")
    
    if nvoxhp:
        eddy_proc.opt(f"--nvoxhp={nvoxhp}")
    
    if initrand:
        eddy_proc.opt(f"--initrand={initrand}")
    
    if ff:
        eddy_proc.opt(f"--ff={ff}")
    
    if repol:
        eddy_proc.opt("--repol")
    
    if ol_nstd:
        eddy_proc.opt(f"--ol_nstd={ol_nstd}")
    
    if ol_nvox:
        eddy_proc.opt(f"--ol_nvox={ol_nvox}")

    if ol_type:
        ol_type: str = ECOLType(ol_type).name
        eddy_proc.opt(f"--ol_type={ol_type}")

    if ol_pos:
        eddy_proc.opt("--ol_pos")
    
    if ol_sqr:
        eddy_proc.opt("--ol_sqr")
    
    if estimate_move_by_susceptibility:
        eddy_proc.opt("--estimate_move_by_susceptibility")
    
    if mbs_niter:
        eddy_proc.opt(f"--mbs_niter={mbs_niter}")
    
    if mbs_lambda:
        eddy_proc.opt(f"--mbs_lambda={mbs_lambda}")
    
    if mbs_ksp:
        eddy_proc.opt(f"--mbs_ksp={mbs_ksp}")
    
    if dont_sep_offs_move:
        eddy_proc.opt("--dont_sep_offs_move")
    
    if dont_peas:
        eddy_proc.opt("--dont_peas")
    
    if data_is_shelled:
        eddy_proc.opt("--data_is_shelled")
    
    if b0_only:
        eddy_proc.opt("--b0_only")
    
    if dont_mask_output:
        eddy_proc.opt("--dont_mask_output")
    
    if verbose:
        eddy_proc.opt("--verbose")
    
    if very_verbose:
        eddy_proc.opt("--very_verbose")
    
    # Eddy slice-to-volume (s2v) correction options
    if mporder:
        eddy_proc.opt(f"--mporder={mporder}")
    
        if slspec:
            slspec: File = File(file=slspec, assert_exists=True)
            eddy_proc.opt(f"--slspec={slspec.file}")
        elif json_file:
            json_file: File = File(file=json_file, assert_exists=True)
            eddy_proc.opt(f"--json={json_file.file}")
        else:
            raise FSLError("Either the 'slspec' or the 'json' option must be specified with the 'mporder' option.")
        
    if s2v_lambda:
        eddy_proc.opt(f"--s2v_lambda={s2v_lambda}")
    
    if s2v_fwhm:
        eddy_proc.opt(f"--s2v_fwhm={s2v_fwhm}")
    
    if s2v_niter:
        eddy_proc.opt(f"--s2v_niter={s2v_niter}")
    
    if s2v_interp: 
        s2v_interp: str = ECInterp(s2v_interp).name
        eddy_proc.opt(f"--s2v_interp={s2v_interp}")
    
    out: str = out.file + '.nii.gz'
    out: NiiFile = NiiFile(file=out, assert_exists=True, validate_nifti=True)

    eddy_proc.run(log=log)

    return (out.file)

def bet(img: str,
        out: str,
        mask: bool = False,
        robust: bool = False,
        seg: bool = True,
        frac_int: Optional[int] = None,
        verbose: bool = False,
        log: Optional[LogFile] = None
       ) -> Tuple[str]:
    """work"""
    img: NiiFile = NiiFile(file=img, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(file=out, assert_exists=False, validate_nifti=False)

    b: Command = Command("bet")
    b.opt(img.abs_path())
    b.opt(out.abs_path())

    if mask:
        b.opt("-m")
        mask_img: str = out.rm_ext() + "_mask.nii.gz"
    else:
        mask_img: str = ""
    
    if robust:
        b.opt("-R")

    if frac_int:
        b.opt("-f")
        b.opt(frac_int)
    
    if seg:
        pass
    else:
        b.opt("-n")
    
    if verbose:
        b.opt("-v")
    
    b.run(log=log)

    return (out.file, 
            mask_img)

def topup(img: str,
          param: str,
          out: str,
          config: Optional[str] = None,
          fout: bool = False,
          iout: bool = False,
          scale: int = 1,
          verbose: bool = False,
          log: Optional[LogFile] = None
         ) -> Tuple[str,str,str]:
    """work"""
    img:    NiiFile = NiiFile(file=img, assert_exists=True, validate_nifti=True)
    out:    NiiFile = NiiFile(file=out, assert_exists=False, validate_nifti=False)
    param:  File = File(file=param, assert_exists=True)

    top: Command = Command("topup")
    top.opt(f"--imain={img.abs_path()}")
    top.opt(f"--out={out.abs_path()}")
    top.opt(f"--datain={param.abs_path()}")
    top.opt(f"--scale={scale}")

    if config:
        config: File = File(file=config, assert_exists=True)
        top.opt(f"--config={config.abs_path()}")
    
    if fout:
        field_img: str = out.rm_ext() + "_field.nii.gz"
        field_img: File = File(file=field_img)
        top.opt(f"--fout={field_img.rm_ext()}")
    else:
        field_img: str = ""
        field_img: File = File(file=field_img)
    
    if iout:
        unwarped_img: str = out.rm_ext() + "_unwarped.nii.gz"
        unwarped_img: NiiFile = NiiFile(file=unwarped_img)
        top.opt(f"--iout={unwarped_img.rm_ext()}")
    else:
        unwarped_img: str = ""
        unwarped_img: File = File(file=unwarped_img)

    if verbose:
        top.opt("--verbose")
    
    top.run(log=log)

    return (out.file,
            field_img.file,
            unwarped_img.file)

def fslreorient2std(img: str,
                    out: Optional[str] = None,
                    out_mat: bool = False,
                    log: Optional[LogFile] = None
                   ) -> Tuple[str,str,str]:
    """work"""
    img: NiiFile = NiiFile(file=img, assert_exists=True, validate_nifti=True)
    reorient: Command = Command("fslreorient2std")

    if out_mat:
        out_mat: str = out.rm_ext() + "_native2std_reorient.mat"
        out_mat: File = File(file=out_mat)
        reorient.opt("-m")
        reorient.opt(f"{out_mat.file}")
    else:
        out_mat: str = ""
        out_mat: File = File(file=out_mat)
    
    reorient.opt(f"{img.file}")
    
    if out:
        out: NiiFile = NiiFile(file=out)
        reorient.opt(f"{out.file}")
    else:
        out: str = ""
        out: NiiFile = NiiFile(file=out)

    reorient.run(log=log)

    return (img.file, 
            out.file, 
            out_mat.file)

def fslroi(img: str,
           out: str,
           xmin: Optional[int] = None,
           xsize: Optional[int] = None,
           ymin: Optional[int] = None,
           ysize: Optional[int] = None,
           zmin: Optional[int] = None,
           zsize: Optional[int] = None,
           tmin: Optional[int] = None,
           tsize: Optional[int] = None,
           log: Optional[LogFile] = None
          ) -> str:
    """work"""
    img: NiiFile = NiiFile(file=img, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(file=out, assert_exists=False, validate_nifti=False)

    if (xmin and xsize and ymin and ysize and zmin and zsize) and (tmin and tsize):
        raise FSLError("Cannot specify both xyz and temporal dimensions.")
    elif (not (xmin and xsize and ymin and ysize and zmin and zsize) and 
            not (tmin and tsize)):
        raise FSLError("Neither xyz nor temporal dimensions were specified.")
    else:
        roi: Command = Command("fslroi")
        roi.opt(img.file)
        roi.opt(out.file)
    
    if (xmin and xsize and ymin and ysize and zmin and zsize):
        roi.opt(xmin); roi.opt(xsize)
        roi.opt(ymin); roi.opt(ysize)
        roi.opt(zmin); roi.opt(zsize)
    elif ((xmin is None) and
            (xsize is None) and
            (ymin is None) and
            (ysize is None) and
            (zmin is None) and
            (zsize is None)):
        pass
    else:
        raise FSLError("Either the xyz min or size was not specified.")

    if (tmin and tsize):
        roi.opt(tmin); roi.opt(tsize)
    elif (tmin is None) and \
        (tsize is None):
        pass
    else:
        raise FSLError("Either the temporal min or size was not specified.")
    
    roi.run(log=log)

    return out.file

def fslmerge(out: str,
             merge_opt: str = "t",
             tr: Optional[float] = None,
             log: Optional[LogFile] = None,
             *args
            ) -> str:
    """Merges a series of 3D NIFTI files."""
    out: NiiFile = NiiFile(file=out)
    merge_opt: str = MergeDim(merge_opt).name

    merge: Command = Command("fslmerge")
    merge.opt(f"-{merge_opt}")
    merge.opt(out.file)
    
    for arg in args:
        img: NiiFile = NiiFile(file=arg, assert_exists=True, validate_nifti=True)
        merge.opt(img.abs_path())
    
    if tr:
        merge.opt(tr)
    
    merge.run(log=log)

    return out.file

def catmats(matdir: str,
            out: str
           ) -> str:
    """Concatenate ``FSL`` linear trasformations files into a single file."""
    mats: List[str] = glob.glob(matdir, "MAT_????")
    with open(out,'w') as output_file:
        for mat_file in mats:
            with open(mat_file) as f:
                output_file.write(f.read())
    return out

def applywarp(src: str,
              ref: str,
              out: str,
              warp: Optional[str] = "",
              premat: Optional[str] = "",
              prematdir: Optional[str] = "",
              postmat: Optional[str] = "",
              postmatdir: Optional[str] = "",
              interp: str = "spline",
              paddingsize: Optional[int] = None,
              abs: bool = False,
              rel: bool = False,
              log: Optional[LogFile] = None
             ) -> str:
    """Applies ``FSL`` warps."""
    assert (warp or premat or postmat or prematdir or postmatdir), \
        "either a warp or mat (premat, postmat or prematdir) must be supplied"
    assert not (premat and prematdir), \
        "cannot use premat and prematdir arguments together"
    assert not (postmat and postmatdir), \
        "cannot use postmat and postmatdir arguments together"

    src: NiiFile = NiiFile(file=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(file=ref, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(file=out, assert_exists=False, validate_nifti=False)

    apply: Command = Command("applywarp")

    apply.opt(f"--in={src.abs_path()}")
    apply.opt(f"--ref={ref.abs_path()}")
    apply.opt(f"--out={out.file}")

    if interp:
        interp: str = RegInterp(interp).name
        apply.opt(f"--interp={interp}")
    
    if prematdir:
        premat: str = os.path.join(prematdir,'allmats.txt')
        premat: str = catmats(matdir=prematdir, out=premat)

    if postmatdir:
        postmat: str = os.path.join(postmatdir,'allmats.txt')
        postmat: str = catmats(matdir=postmatdir, out=postmat)

    if warp:
        warp: NiiFile = NiiFile(file=warp, assert_exists=True, validate_nifti=True)
        apply.opt(f"--warp={warp}")

    if premat:
        apply.opt(f"--premat={premat}")

    if postmat:
        apply.opt(f"--postmat={postmat}")

    if paddingsize and (isinstance(paddingsize,int)):
        apply.opt(f"--paddingsize={paddingsize}")

    if abs:
        apply.opt("--abs")

    if rel:
        apply.opt("--rel")

    apply.run(log=log)
    return out.file


def invxfm(inmat: str,
           outmat: str,
           log: Optional[LogFile] = None
          ) -> str:
    """Inverts ``FSL`` transformation matrices."""
    inmat:  File = File(file=inmat, assert_exists=True)
    outmat: File = File(file=outmat, assert_exists=False)

    inv: Command = Command("invxfm")

    inv.opt("-omat")
    inv.opt(outmat.file)

    inv.opt("-inverse")
    inv.opt(inmat.file)

    inv.run(log=log)
    return outmat

def applyxfm(src: str,
             ref: str,
             mat: str,
             out: str,
             interp: Optional[str] = "trilinear",
             log: Optional[LogFile] = None
            ) -> str:
    """Applies ``FSL`` transformation matrices."""
    src: NiiFile = NiiFile(file=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(file=ref, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(file=out, assert_exists=False, validate_nifti=False)
    mat: File = File(file=mat, assert_exists=True)

    xfm: Command = Command("flirt")

    xfm.opt("-init")
    xfm.opt(mat.file)

    xfm.opt("-in")
    xfm.opt(src.file)

    xfm.opt("-ref")
    xfm.opt(ref.file)

    xfm.opt("-applyxfm")

    xfm.opt("-out")
    xfm.opt(out.file)

    interp: str = RegInterp(interp).name
    xfm.opt("-interp")
    xfm.opt(interp)

    xfm.run(log=log)
    return out.file

def apply_isoxfm(src: str,
                 ref: str,
                 res: int,
                 out: str,
                 interp: Optional[str] = "interp",
                 log: Optional[str] = None
                ) -> str:
    """Resamples images to an isometric resolution."""
    src: NiiFile = NiiFile(file=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(file=ref, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(file=out, assert_exists=False, validate_nifti=False)

    fsldir: str = os.getenv('FSLDIR',None)
    assert fsldir is not None, 'FSLDIR environment must be set.'

    ident: str = os.path.join(fsldir,'etc/flirtsch/ident.mat')

    xfm: Command = Command("flirt")

    xfm.opt("-init")
    xfm.opt(ident)

    xfm.opt("-in")
    xfm.opt(src.file)

    xfm.opt("-ref")
    xfm.opt(ref.file)

    xfm.opt("-applyisoxfm")
    xfm.opt(f"{res}")

    xfm.opt("-out")
    xfm.opt(out.file)

    interp: str = RegInterp(interp).name
    xfm.opt("-interp")
    xfm.opt(interp)

    xfm.run(log=log)
    return out.file

def concatxfm(inmat1: str,
              inmat2: str,
              outmat: str,
              log: Optional[LogFile] = None
             ) -> str:
    """Concatenates two ``FSL`` transformation matrices."""
    inmat1: File = File(file=inmat1, assert_exists=True)
    inmat2: File = File(file=inmat2, assert_exists=True)

    concat: Command = Command("convert_xfm")

    concat.opt("-omat")
    concat.opt(outmat)

    concat.opt("-concat")
    concat.opt(inmat1.file)
    concat.opt(inmat2.file)

    concat.run(log=log)
    return outmat

def invwarp(inwarp: str,
            ref: str,
            outwarp:str,
            rel: bool = False,
            abs: bool = False,
            noconstraint: bool = False,
            jmin: float = 0.01,
            jmax: float = 100.0,
            verbose: bool = False,
            log: Optional[LogFile] = None
           ) -> str:
    """Invert existing warps."""
    inwarp:  NiiFile = NiiFile(file=inwarp, assert_exists=True, validate_nifti=True)
    ref:     NiiFile = NiiFile(file=ref, assert_exists=True, validate_nifti=True)
    outwarp: NiiFile = NiiFile(file=outwarp, assert_exists=False, validate_nifti=False)

    inv: Command = Command("invwarp")

    # Required arguments
    inv.opt(f"--warp={inwarp.abs_path()}")
    inv.opt(f"--ref={ref.abs_path()}")
    inv.opt(f"--out={outwarp.file}")

    # Optional arguments
    if rel:
        inv.opt("--rel")
    
    if abs:
        inv.opt("--abs")
    
    if noconstraint:
        inv.opt("--noconstraint")
    
    if jmin:
        inv.opt(f"--jmin={jmin}")
    
    if jmax:
        inv.opt(f"--jamx={jmax}")
    
    if verbose:
        inv.opt("--verbose")
    
    inv.run(log=log)
    return outwarp.file

def convertwarp(out: str,
                ref: str,
                warp1: Optional[str] = None,
                warp2: Optional[str] = None,
                premat: Optional[str] = None,
                midmat: Optional[str] = None,
                postmat: Optional[str] = None,
                shiftmap: Optional[str] = None,
                shiftdir: Optional[str] = None,
                abs: bool = False,
                absout: bool = False,
                rel: bool = False,
                relout: bool = False,
                log: Optional[LogFile] = None
               ) -> str:
    """Convert ``FSL`` non-linear warps."""
    assert (warp1 or warp2 or premat or midmat or postmat), \
        "either a warp (warp1 or warp2) or mat (premat, midmat, or " + \
        "postmat) must be supplied"
    
    ref: NiiFile = NiiFile(file=ref, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(file=out)

    warp1: NiiFile = NiiFile(file="", assert_exists=False, validate_nifti=False)
    warp2: NiiFile = NiiFile(file="", assert_exists=False, validate_nifti=False)
    shiftmap: NiiFile = NiiFile(file="", assert_exists=False, validate_nifti=False)

    postmat: File = File(file="", assert_exists=False)
    premat: File = File(file="", assert_exists=False)
    midmat: File = File(file="", assert_exists=False)

    cvwarp: Command = Command("convertwarp")
    
    cvwarp.opt(f"--ref={ref.file}")
    cvwarp.opt(f"--out={out.file}")

    if warp1:
        warp1: NiiFile = NiiFile(file=warp1, assert_exists=True, validate_nifti=True)
        cvwarp.opt(f"--warp1={warp1.file}")
    
    if warp2:
        warp2: NiiFile = NiiFile(file=warp2, assert_exists=True, validate_nifti=True)
        cvwarp.opt(f"--warp2={warp2.file}")
    
    if premat:
        premat: File = File(file=premat, assert_exists=True)
        cvwarp.opt(f"--premat={premat.file}")
    
    if midmat:
        midmat: File = File(file=midmat, assert_exists=True)
        cvwarp.opt(f"--midmat={midmat.file}")

    if postmat:
        postmat: File = File(file=postmat, assert_exists=True)
        cvwarp.opt(f"--postmat={postmat.file}")

    if shiftmap:
        shiftmap: NiiFile = NiiFile(file=shiftmap, assert_exists=True, validate_nifti=True)
        cvwarp.opt(f"--shiftmap={shiftmap.file}")

    if (shiftdir == 'x' or 
        shiftdir == 'x-' or 
        shiftdir == 'y' or 
        shiftdir == 'y-' or 
        shiftdir == 'z' or 
        shiftdir == 'z-'):
        cvwarp.opt(f"--shiftdir={shiftdir}")
    else:
        raise FSLError(f"Invalid input: {shiftdir}. Valid inputs for 'shiftdir' option includes: x,y,z,x-,y-,z-")

    if absout:
        cvwarp.opt(f"--absout")

    if relout:
        cvwarp.opt(f"--relout")

    if abs:
        cvwarp.opt(f"--abs")

    if rel:
        cvwarp.opt(f"--rel")

    cvwarp.run(log=log)

    return out.file


def fugue():
    """FMRIB's Utility for Geometric Unwarping of EPIs."""
    pass

def flirt():
    pass

def melodic(input: str,
            outdir: str,
            dim: Optional[int] = 0,
            tr: Optional[float] = None,
            mmthresh: Optional[float] = None,
            report: bool = False,
            prefix: Optional[str] = "",
            nomask: bool = False,
            updatemask: bool = False,
            nobet: bool = False,
            mask: Optional[str] = "",
            Oall: bool = False,
            verbose: bool = False,
            log: Optional[LogFile] = None
           ) -> str:
    """Multivariate Exploratory Linear Optimised ICA."""
    input: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)

    mel: Command = Command("melodic")

    mel.opt(f"--in={input.file}")
    mel.opt(f"--outdir={outdir}")

    if mmthresh:
        mel.opt(f"--mmthresh={mmthresh}")
    
    if dim:
        mel.opt(f"--dim={dim}")
    
    if tr:
        mel.opt(f"--tr={tr}")
    
    if report:
        mel.opt(f"--report")
    
    if Oall:
        mel.opt(f"--Oall")
    
    # if nomask:



def fsl_regfilt():
    pass

def mcflirt():
    pass

def slicer():
    pass

def cluster():
    pass

class fslmaths:
    """work"""

    def __init__(self,
                 img: str,
                 dt: Optional[str] = None,
                ):
        """Constructor"""
        img: NiiFile = NiiFile(file=img, assert_exists=True, validate_nifti=True)
        self._maths: Command = Command("fslmaths")

        if dt:
            dt: str = FSLDataType(dt).name
            self._maths.opt("-dt")
            self._maths.opt(dt)
        
        self._maths.opt(img.file)
    
    def abs(self):
        """work"""
        self._maths.opt("-abs")
        return self
    
    def bin(self):
        """work"""
        self._maths.opt("-bin")
        return self
    
    def binv(self):
        """work"""
        self._maths.opt("-binv")
        return self
    
    def recip(self):
        """work"""
        self._maths.opt("-recip")
        return self
    
    def Tmean(self):
        """work"""
        self._maths.opt("-Tmean")
        return self
    
    def Tstd(self):
        """work"""
        self._maths.opt("-Tstd")
        return self
    
    def Tmin(self):
        """work"""
        self._maths.opt("-Tmin")
        return self
    
    def Tmax(self):
        """work"""
        self._maths.opt("-Tmax")
        return self
    
    def sqrt(self):
        """work"""
        self._maths.opt("-sqrt")
        return self
    
    def sqr(self):
        """work"""
        self._maths.opt("-sqr")
        return self
    
    def fillh(self):
        """work"""
        self._maths.opt("-fillh")
        return self
    
    def ero(self,
            repeat: int = 1
           ):
        """work"""
        for _ in range(repeat):
            self._maths.opt("-ero")
        return self
    
    def dilM(self,
             repeat: int = 1
            ):
        """work"""
        for _ in range(repeat):
            self._maths.opt("-dilM")
        return self
    
    def dilF(self,
             repeat: int = 1
            ):
        """work"""
        for _ in range(repeat):
            self._maths.opt("-dilF")
        return self

    def add(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.opt("-add")

        if isinstance(input,int):
            self._maths.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._maths.opt(img.file)
        return self
    
    def sub(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.opt("-sub")

        if isinstance(input,int):
            self._maths.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._maths.opt(img.file)
        return self
    
    def mul(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.opt("-mul")

        if isinstance(input,int):
            self._maths.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._maths.opt(img.file)
        return self
    
    def div(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.opt("-div")

        if isinstance(input,int):
            self._maths.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._maths.opt(img.file)
        return self
    
    def mas(self,
            img: str
           ):
        """work"""
        img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
        self._maths.opt("-mas")
        self._maths.opt(img.file)
        return self

    def rem(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.opt("-rem")

        if isinstance(input,int):
            self._maths.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._maths.opt(img.file)
        return self
    
    def thr(self,
            num: Union[int,float]
           ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._maths.opt("-thr")
            self._maths.opt(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self
    
    def uthr(self,
             num: Union[int,float]
            ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._maths.opt("-uthr")
            self._maths.opt(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self

    def inm(self,
            num: Union[int,float]
           ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._maths.opt("-inm")
            self._maths.opt(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self
    
    def bptf(self,
             high_pass: Union[int,float],
             low_pass: Union[int,float],
             tr: Optional[Union[int,float]] = None,
             input_is_hz: bool = False,
             input_is_sec: bool = False
            ):
        """
        Input is assumed to be in sigma (in volume).
        """
        if input_is_hz and input_is_sec:
            raise RuntimeError("Both 'input_is_hz' and 'input_is_sec' were specified. ONLY ONE of these options may be specified.")
        elif (input_is_hz or input_is_sec) and (not tr):
            raise FSLError("The TR (Repetition Time) is required when either the 'input_is_hz' or 'input_is_sec' options are specified.")

        def _compute_sigma(tr: Union[int,float],
                           hz: Union[int,float],
                           compute_low_pass: bool = False
                          ):
            """
            relevant links: 
                * https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;fc5b33c5.1205
                * https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;f6fd75a6.1709
            """
            if ((isinstance(tr,int) or isinstance(tr,float)) and
                    (isinstance(hz,int) or isinstance(hz,float))):

                if compute_low_pass:
                    fwhm_kernel: float = 18.0
                else:
                    fwhm_kernel: float = sqrt(8*log(2,EULER_CONST))
                
                try:
                    sigma_vol:   float = 1/(fwhm_kernel * tr * hz)
                except ZeroDivisionError:
                    sigma_vol: float = 0.0
                return sigma_vol
            else:
                raise TypeError(f"Input TR: {tr} or Hz: {hz} is not a number.")
        
        if input_is_sec:
            try:
                high_pass: float = _compute_sigma(tr=tr, 
                                                  hz=(1/high_pass), 
                                                  compute_low_pass=False)
            except ZeroDivisionError:
                high_pass: float = 0.0
            
            try:
                low_pass:  float = _compute_sigma(tr=tr, 
                                                  hz=(1/low_pass), 
                                                  compute_low_pass=True)
            except ZeroDivisionError:
                low_pass: float = 0.0
            
        elif input_is_hz:
            high_pass: float = _compute_sigma(tr=tr, 
                                              hz=high_pass,
                                              compute_low_pass=False)
            low_pass:  float = _compute_sigma(tr=tr, 
                                              hz=low_pass,
                                              compute_low_pass=True)

        if ((isinstance(high_pass,int) or isinstance(high_pass,float)) and
                (isinstance(low_pass,int) or isinstance(low_pass,float))):
            self._maths.opt("-bptf")
            self._maths.opt(f"{high_pass}")
            self._maths.opt(f"{low_pass}")
        else:
            raise TypeError(f"Input high_pass: {high_pass} or low_pass: {low_pass} is not a number.")
        
        return self
    
    def run(self,
            out: str,
            odt: Optional[str] = None,
            log: Optional[LogFile] = None
           ) -> str:
        """work"""
        self._maths.opt(out)

        if odt:
            odt: str = FSLDataType(odt).name
            self._maths.opt("-odt")
            self._maths.opt(odt)
            
        self._maths.run(log=log)
        return out
