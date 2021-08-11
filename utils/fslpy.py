# -*- coding: utf-8 -*-
"""Wrappers functions for FSL executable binaries.
"""
# TODO: 
#   
#   * Write doc-strings for each wrapper function
#   * Add verbose options to all wrapper functions
#   * Make FSLDIR global variable
#       * Add the option prepend FSLDIR/bin to each command object

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

def fnirt(src, 
          ref, 
          aff: Optional[str] = None,
          imprefm: Optional[int] = None, 
          impinm: Optional[int] = None, 
          applyrefmask: Optional[str] = None,
          applyinmask: Optional[str] = None,
          subsamp: Optional[str] = None,
          miter: Optional[str] = None, 
          infwhm: Optional[str] = None,
          reffwhm: Optional[str] = None, 
          lmbda: Optional[str] = None, 
          estint: Optional[str] = None, 
          warpres: Optional[str] = None, 
          ssqlambda: Optional[str] = None,
          regmod: Optional[str] = None, 
          intmod: Optional[str] = None, 
          intorder: Optional[str] = None, 
          biasres: Optional[str] = None,
          biaslambda: Optional[str] = None, 
          refderiv: Optional[str] = None, 
          cout: Optional[str] = None, 
          intout: Optional[str] = None, 
          refout: Optional[str] = None,
          iout: Optional[str] = None, 
          interp: Optional[str] = None, 
          inwarp: Optional[str] = None, 
          minmet: Optional[str] = None, 
          verbose: bool = False,
          intin: Optional[str] = None, 
          jout: Optional[str] = None,
          log: Optional[LogFile] = None
         ) -> str:
    """Perform/compute non-linear image registrations."""
    src: NiiFile = NiiFile(file=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(file=ref, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("fnirt")

    # TODO: Verify that input options are either:
    #   * str
    #   * int or float

    if aff:
        aff: File = File(file=aff, assert_exists=True)
        cmd.opt(f"--aff={aff.file}")
    
    if imprefm:
        imprefm: File = File(file=imprefm, assert_exists=True)
        cmd.opt(f"--imprefm={imprefm.file}")

    if impinm:
        impinm: File = File(file=impinm, assert_exists=True)
        cmd.opt(f"--impinm={impinm.file}")

    if applyrefmask:
        applyrefmask: File = File(file=applyrefmask, assert_exists=True)
        cmd.opt(f"--applyrefmask={applyrefmask.file}")

    if applyinmask:
        applyinmask: File = File(file=applyinmask, assert_exists=True)
        cmd.opt(f"--applyinmask={applyinmask.file}")

    if subsamp:
        pass

    if miter:
        pass

    if infwhm:
        pass

    if reffwhm:
        pass

    if lmbda:
        pass

    if estint:
        pass

    if warpres:
        pass

    if ssqlambda:
        pass

    if regmod:
        pass

    if intmod:
        pass

    if intorder:
        pass

    if biasres:
        pass

    if biaslambda:
        pass

    if refderiv:
        pass

    if cout:
        pass

    if intout:
        pass

    if refout:
        pass

    if iout:
        pass

    if interp:
        pass

    if inwarp:
        pass

    if minmet:
        pass

    if verbose:
        pass

    if intin:
        pass

    if jout:
        pass

    cmd.run(log=log)

    return None


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
            cmd: Command = Command(eddy_cmd)
            
            (return_code, 
             _, 
             _) = cmd.run(log=log)
            
            # Check if the dependency is met and if the return code
            #   of an empty command options list returns 1.
            if (cmd.check_dependency()) and (return_code == 1):
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

    cmd.opt(f"--imain={img.file}")
    cmd.opt(f"--mask={mask.file}")
    cmd.opt(f"--index={idx.file}")
    cmd.opt(f"--acqp={acqp.file}")
    cmd.opt(f"--bvecs={bvecs.file}")
    cmd.opt(f"--bvals={bvals.file}")
    cmd.opt(f"--out={out.file}")

    # Conventional eddy current correction options
    if topup:
        topup: NiiFile = NiiFile(file=topup, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--topup={topup.rm_ext()}")
    
    if field:
        field: NiiFile = NiiFile(file=field, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--field={field.file}")
    
    if field_mat:
        field_mat: File = File(file=field_mat, assert_exists=True)
        cmd.opt(f"--field_mat={field_mat.file}")
    
    if flm:
        flm: str = ECModelFLM(flm).name
        cmd.opt(f"--flm={flm}")
    
    if slm:
        slm: str = ECModelSLM(slm).name
        cmd.opt(f"--slm={slm}")
    
    if fwhm:
        cmd.opt(f"--fwhm={fwhm}")
    
    if niter:
        cmd.opt(f"--niter={niter}")
    
    if cnr_maps:
        cmd.opt("--cnr_maps")
    
    if residuals:
        cmd.opt("--residuals")
    
    if fep:
        cmd.opt("--fep")
    
    if interp:
        interp: str = ECInterp(interp).name
        cmd.opt(f"--interp={interp}")
    
    if resamp:
        resamp: str = ECresamp(resamp).name
        cmd.opt(f"--resamp={resamp}")
    
    if nvoxhp:
        cmd.opt(f"--nvoxhp={nvoxhp}")
    
    if initrand:
        cmd.opt(f"--initrand={initrand}")
    
    if ff:
        cmd.opt(f"--ff={ff}")
    
    if repol:
        cmd.opt("--repol")
    
    if ol_nstd:
        cmd.opt(f"--ol_nstd={ol_nstd}")
    
    if ol_nvox:
        cmd.opt(f"--ol_nvox={ol_nvox}")

    if ol_type:
        ol_type: str = ECOLType(ol_type).name
        cmd.opt(f"--ol_type={ol_type}")

    if ol_pos:
        cmd.opt("--ol_pos")
    
    if ol_sqr:
        cmd.opt("--ol_sqr")
    
    if estimate_move_by_susceptibility:
        cmd.opt("--estimate_move_by_susceptibility")
    
    if mbs_niter:
        cmd.opt(f"--mbs_niter={mbs_niter}")
    
    if mbs_lambda:
        cmd.opt(f"--mbs_lambda={mbs_lambda}")
    
    if mbs_ksp:
        cmd.opt(f"--mbs_ksp={mbs_ksp}")
    
    if dont_sep_offs_move:
        cmd.opt("--dont_sep_offs_move")
    
    if dont_peas:
        cmd.opt("--dont_peas")
    
    if data_is_shelled:
        cmd.opt("--data_is_shelled")
    
    if b0_only:
        cmd.opt("--b0_only")
    
    if dont_mask_output:
        cmd.opt("--dont_mask_output")
    
    if verbose:
        cmd.opt("--verbose")
    
    if very_verbose:
        cmd.opt("--very_verbose")
    
    # Eddy slice-to-volume (s2v) correction options
    if mporder:
        cmd.opt(f"--mporder={mporder}")
    
        if slspec:
            slspec: File = File(file=slspec, assert_exists=True)
            cmd.opt(f"--slspec={slspec.file}")
        elif json_file:
            json_file: File = File(file=json_file, assert_exists=True)
            cmd.opt(f"--json={json_file.file}")
        else:
            raise FSLError("Either the 'slspec' or the 'json' option must be specified with the 'mporder' option.")
        
    if s2v_lambda:
        cmd.opt(f"--s2v_lambda={s2v_lambda}")
    
    if s2v_fwhm:
        cmd.opt(f"--s2v_fwhm={s2v_fwhm}")
    
    if s2v_niter:
        cmd.opt(f"--s2v_niter={s2v_niter}")
    
    if s2v_interp: 
        s2v_interp: str = ECInterp(s2v_interp).name
        cmd.opt(f"--s2v_interp={s2v_interp}")
    
    out: str = out.file + '.nii.gz'
    out: NiiFile = NiiFile(file=out, assert_exists=True, validate_nifti=True)

    cmd.run(log=log)

    # TODO:
    #   * Get the corresponding output files for 
    #       each eddy option.

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

    cmd: Command = Command("bet")
    cmd.opt(img.abs_path())
    cmd.opt(out.abs_path())

    if mask:
        cmd.opt("-m")
        mask_img: str = out.rm_ext() + "_mask.nii.gz"
    else:
        mask_img: str = ""
    
    if robust:
        cmd.opt("-R")

    if frac_int:
        cmd.opt("-f")
        cmd.opt(frac_int)
    
    if seg:
        pass
    else:
        cmd.opt("-n")
    
    if verbose:
        cmd.opt("-v")
    
    cmd.run(log=log)

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

    cmd: Command = Command("topup")
    cmd.opt(f"--imain={img.abs_path()}")
    cmd.opt(f"--out={out.abs_path()}")
    cmd.opt(f"--datain={param.abs_path()}")
    cmd.opt(f"--scale={scale}")

    if config:
        config: File = File(file=config, assert_exists=True)
        cmd.opt(f"--config={config.abs_path()}")
    
    if fout:
        field_img: str = out.rm_ext() + "_field.nii.gz"
        field_img: File = File(file=field_img)
        cmd.opt(f"--fout={field_img.rm_ext()}")
    else:
        field_img: str = ""
        field_img: File = File(file=field_img)
    
    if iout:
        unwarped_img: str = out.rm_ext() + "_unwarped.nii.gz"
        unwarped_img: NiiFile = NiiFile(file=unwarped_img)
        cmd.opt(f"--iout={unwarped_img.rm_ext()}")
    else:
        unwarped_img: str = ""
        unwarped_img: File = File(file=unwarped_img)

    if verbose:
        cmd.opt("--verbose")
    
    cmd.run(log=log)

    # TODO:
    #   * Get the corresponding output files for 
    #       each topup options.

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
    cmd: Command = Command("fslreorient2std")

    if out_mat:
        out_mat: str = out.rm_ext() + "_native2std_cmd.mat"
        out_mat: File = File(file=out_mat)
        cmd.opt("-m")
        cmd.opt(f"{out_mat.file}")
    else:
        out_mat: str = ""
        out_mat: File = File(file=out_mat)
    
    cmd.opt(f"{img.file}")
    
    if out:
        out: NiiFile = NiiFile(file=out)
        cmd.opt(f"{out.file}")
    else:
        out: str = ""
        out: NiiFile = NiiFile(file=out)

    cmd.run(log=log)

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
        cmd: Command = Command("fslroi")
        cmd.opt(img.file)
        cmd.opt(out.file)
    
    if (xmin and xsize and ymin and ysize and zmin and zsize):
        cmd.opt(xmin); cmd.opt(xsize)
        cmd.opt(ymin); cmd.opt(ysize)
        cmd.opt(zmin); cmd.opt(zsize)
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
        cmd.opt(tmin); cmd.opt(tsize)
    elif (tmin is None) and \
        (tsize is None):
        pass
    else:
        raise FSLError("Either the temporal min or size was not specified.")
    
    cmd.run(log=log)

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

    cmd: Command = Command("fslmerge")
    cmd.opt(f"-{merge_opt}")
    cmd.opt(out.file)
    
    for arg in args:
        img: NiiFile = NiiFile(file=arg, assert_exists=True, validate_nifti=True)
        cmd.opt(img.abs_path())
    
    if tr:
        cmd.opt(tr)
    
    cmd.run(log=log)

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

    cmd: Command = Command("applywarp")

    cmd.opt(f"--in={src.abs_path()}")
    cmd.opt(f"--ref={ref.abs_path()}")
    cmd.opt(f"--out={out.file}")

    if interp:
        interp: str = RegInterp(interp).name
        cmd.opt(f"--interp={interp}")
    
    if prematdir:
        premat: str = os.path.join(prematdir,'allmats.txt')
        premat: str = catmats(matdir=prematdir, out=premat)

    if postmatdir:
        postmat: str = os.path.join(postmatdir,'allmats.txt')
        postmat: str = catmats(matdir=postmatdir, out=postmat)

    if warp:
        warp: NiiFile = NiiFile(file=warp, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--warp={warp}")

    if premat:
        cmd.opt(f"--premat={premat}")

    if postmat:
        cmd.opt(f"--postmat={postmat}")

    if paddingsize and (isinstance(paddingsize,int)):
        cmd.opt(f"--paddingsize={paddingsize}")

    if abs:
        cmd.opt("--abs")

    if rel:
        cmd.opt("--rel")

    cmd.run(log=log)
    return out.file


def invxfm(inmat: str,
           outmat: str,
           log: Optional[LogFile] = None
          ) -> str:
    """Inverts ``FSL`` transformation matrices."""
    inmat:  File = File(file=inmat, assert_exists=True)
    outmat: File = File(file=outmat, assert_exists=False)

    cmd: Command = Command("invxfm")

    cmd.opt("-omat")
    cmd.opt(outmat.file)

    cmd.opt("-inverse")
    cmd.opt(inmat.file)

    cmd.run(log=log)
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

    cmd: Command = Command("flirt")

    cmd.opt("-init")
    cmd.opt(mat.file)

    cmd.opt("-in")
    cmd.opt(src.file)

    cmd.opt("-ref")
    cmd.opt(ref.file)

    cmd.opt("-applyxfm")

    cmd.opt("-out")
    cmd.opt(out.file)

    interp: str = RegInterp(interp).name
    cmd.opt("-interp")
    cmd.opt(interp)

    cmd.run(log=log)
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

    cmd: Command = Command("flirt")

    cmd.opt("-init")
    cmd.opt(ident)

    cmd.opt("-in")
    cmd.opt(src.file)

    cmd.opt("-ref")
    cmd.opt(ref.file)

    cmd.opt("-applyisoxfm")
    cmd.opt(f"{res}")

    cmd.opt("-out")
    cmd.opt(out.file)

    interp: str = RegInterp(interp).name
    cmd.opt("-interp")
    cmd.opt(interp)

    cmd.run(log=log)
    return out.file

def concatxfm(inmat1: str,
              inmat2: str,
              outmat: str,
              log: Optional[LogFile] = None
             ) -> str:
    """Concatenates two ``FSL`` transformation matrices."""
    inmat1: File = File(file=inmat1, assert_exists=True)
    inmat2: File = File(file=inmat2, assert_exists=True)

    cmd: Command = Command("convert_xfm")

    cmd.opt("-omat")
    cmd.opt(outmat)

    cmd.opt("-concat")
    cmd.opt(inmat1.file)
    cmd.opt(inmat2.file)

    cmd.run(log=log)
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

    cmd: Command = Command("invwarp")

    # Required arguments
    cmd.opt(f"--warp={inwarp.abs_path()}")
    cmd.opt(f"--ref={ref.abs_path()}")
    cmd.opt(f"--out={outwarp.file}")

    # Optional arguments
    if rel:
        cmd.opt("--rel")
    
    if abs:
        cmd.opt("--abs")
    
    if noconstraint:
        cmd.opt("--noconstraint")
    
    if jmin:
        cmd.opt(f"--jmin={jmin}")
    
    if jmax:
        cmd.opt(f"--jamx={jmax}")
    
    if verbose:
        cmd.opt("--verbose")
    
    cmd.run(log=log)
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

    cmd: Command = Command("convertwarp")
    
    cmd.opt(f"--ref={ref.file}")
    cmd.opt(f"--out={out.file}")

    if warp1:
        warp1: NiiFile = NiiFile(file=warp1, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--warp1={warp1.file}")
    
    if warp2:
        warp2: NiiFile = NiiFile(file=warp2, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--warp2={warp2.file}")
    
    if premat:
        premat: File = File(file=premat, assert_exists=True)
        cmd.opt(f"--premat={premat.file}")
    
    if midmat:
        midmat: File = File(file=midmat, assert_exists=True)
        cmd.opt(f"--midmat={midmat.file}")

    if postmat:
        postmat: File = File(file=postmat, assert_exists=True)
        cmd.opt(f"--postmat={postmat.file}")

    if shiftmap:
        shiftmap: NiiFile = NiiFile(file=shiftmap, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--shiftmap={shiftmap.file}")

    if (shiftdir == 'x' or 
        shiftdir == 'x-' or 
        shiftdir == 'y' or 
        shiftdir == 'y-' or 
        shiftdir == 'z' or 
        shiftdir == 'z-'):
        cmd.opt(f"--shiftdir={shiftdir}")
    else:
        raise FSLError(f"Invalid input: {shiftdir}. Valid inputs for 'shiftdir' option includes: x,y,z,x-,y-,z-")

    if absout:
        cmd.opt(f"--absout")

    if relout:
        cmd.opt(f"--relout")

    if abs:
        cmd.opt(f"--abs")

    if rel:
        cmd.opt(f"--rel")

    cmd.run(log=log)

    return out.file


def fugue():
    """FMRIB's Utility for Geometric Unwarping of EPIs."""
    pass

def flirt(src: str,
          ref: str,
          out: Optional[str] = None,
          omat: Optional[str] = None,
          dof: Optional[int] = None,
          cost: Optional[str] = None,
          wmseg: Optional[str] = None,
          init: Optional[str] = None,
          schedule: Optional[str] = None,
          echospacing: Optional[float] = None,
          pedir: Optional[str] = None,
          fieldmap: Optional[str] = None,
          fieldmapmask: Optional[str] = None,
          bbrslope: Optional[float] = None,
          bbrtype: Optional[str] = None,
          interp: Optional[str] = None,
          refweight: Optional[str] = None,
          applyisoxfm: Optional[float] = None,
          usesqform: bool = False,
          nosearch: bool = False, 
          verbose: Optional[int] = 0,
          searchrx: Tuple[int,int] = None, 
          searchry: Tuple[int,int] = None, 
          searchrz: Tuple[int,int] = None,
          log: Optional[LogFile] = None
         ) -> Tuple[str]:
    """work"""
    src: NiiFile = NiiFile(file=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(file=ref, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("flirt")
    
    cmd.opt("-in");  cmd.opt(f"{src.file}")
    cmd.opt("-ref"); cmd.opt(f"{ref.file}")

    if out:
        cmd.opt("-out"); cmd.opt(f"{out}")
    
    if omat:
        cmd.opt("-omat"); cmd.opt(f"{omat}")
    
    if dof:
        cmd.opt("-dof"); cmd.opt(f"{dof}")
    
    if cost:
        # TODO: set enum here
        cmd.opt("-cost"); cmd.opt(f"{cost}")
    
    if wmseg:
        wmseg: NiiFile = NiiFile(file=wmseg, assert_exists=True, validate_nifti=True)
        cmd.opt("-wmseg"); cmd.opt(f"{wmseg.file}")
    
    if init:
        init: File = File(file=init, assert_exists=True)
        cmd.opt("-init"); cmd.opt(f"{init.file}")
    
    if schedule:
        schedule: File = File(file=schedule, assert_exists=True)
        cmd.opt("-schedule"); cmd.opt(f"{schedule.file}")
    
    if echospacing:
        cmd.opt("-echospacing"); cmd.opt(f"{echospacing}")
    
    if pedir:
        cmd.opt("-pedir"); cmd.opt(f"{pedir}")
    
    if fieldmap:
        fieldmap: NiiFile = NiiFile(file=fieldmap, assert_exists=True, validate_nifti=True)
        cmd.opt("-fieldmap"); cmd.opt(f"{fieldmap.file}")
    
    if fieldmapmask:
        fieldmapmask: NiiFile = NiiFile(file=fieldmapmask, assert_exists=True, validate_nifti=True)
        cmd.opt("-fieldmapmask"); cmd.opt(f"{fieldmapmask.file}")
    
    if bbrslope:
        cmd.opt("-bbrslope"); cmd.opt(f"{bbrslope}")
    
    if bbrtype:
        # TODO: set enum here
        cmd.opt("-bbrtype"); cmd.opt(f"{bbrtype}")
    
    if interp:
        # TODO: set enum here
        cmd.opt("-interp"); cmd.opt(f"{interp}")
    
    if refweight:
        refweight: NiiFile = NiiFile(file=refweight, assert_exists=True, validate_nifti=True)
        cmd.opt("-refweight"); cmd.opt(f"{refweight.file}")
    
    if applyisoxfm:
        cmd.opt("-applyisoxfm"); cmd.opt(f"{applyisoxfm}")
    
    if usesqform:
        cmd.opt("-usesqform")
    
    if nosearch:
        cmd.opt("-nosearch")
    
    if verbose:
        cmd.opt("-verbose"); cmd.opt(f"{verbose}")
    
    if searchrx:
        cmd.opt("-searchrx")
        cmd.opt(f"{searchrx[0]}")
        cmd.opt(f"{searchrx[1]}")
    
    if searchry:
        cmd.opt("-searchry")
        cmd.opt(f"{searchry[0]}")
        cmd.opt(f"{searchry[1]}")
    
    if searchrz:
        cmd.opt("-searchrz")
        cmd.opt(f"{searchrz[0]}")
        cmd.opt(f"{searchrz[1]}")
    
    cmd.run(log=log)

    # TODO: set all possible output types here,
    #   set, their defaults above.

    return (out,
            omat)

def melodic(input: str,
            outdir: str,
            dim: Optional[int] = 0,
            tr: Optional[float] = None,
            mmthresh: Optional[float] = None,
            report: bool = False,
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

    cmd: Command = Command("melodic")

    cmd.opt(f"--in={input.file}")
    cmd.opt(f"--outdir={outdir}")

    if mask:
        mask: NiiFile = NiiFile(file=mask, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--mask={mask.file}")

    if mmthresh:
        cmd.opt(f"--mmthresh={mmthresh}")
    
    if dim:
        cmd.opt(f"--dim={dim}")
    
    if tr:
        cmd.opt(f"--tr={tr}")
    
    if report:
        cmd.opt("--report")
    
    if Oall:
        cmd.opt("--Oall")
    
    if nomask:
        cmd.opt("--nomask")
    
    if updatemask:
        cmd.opt("--update_mask")
    
    if nobet:
        cmd.opt("--nobet")
    
    if verbose:
        cmd.opt("--verbose")
    
    cmd.run(log=log)

    # TODO:
    #   * Get the corresponding output files for 
    #       each melodic options.

    return outdir
    
def fsl_regfilt(infile: str,
                outfile: str,
                mix: str,
                ics: List[int],
                log: Optional[LogFile] = None
               ) -> str:
    """Data de-noising by regression.
    """
    infile: NiiFile = NiiFile(file=infile, assert_exists=True, validate_nifti=True)
    outfile: NiiFile = NiiFile(file=outfile, assert_exists=False, validate_nifti=False)
    mix: File = File(file=mix, assert_exists=True)

    # NOTE: Semi-original code from dchp-fmri:
    #   link: https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/blob/master/dhcp/util/fslpy.py
    #   
    #   This approach shifted the IC number by 1.
    # 
    # icstr: str = '"'
    # for i in range(0, len(ics)-1):
    #     icstr: str = icstr + f"{ics[i]+1},"
    # icstr: str = icstr + f"{ics[-1]+1}"

    icstr: str = '"'
    for i in range(0, len(ics)-1):
        icstr: str = icstr + f"{ics[i]},"
    icstr: str = icstr + f'{ics[-1]}"'

    cmd: Command = Command("fsl_regfilt")

    cmd.opt(f"--in={infile.file}")
    cmd.opt(f"--out={outfile.file}")
    cmd.opt(f"--design={mix.file}")
    cmd.opt(f"--filter={icstr}")

    cmd.run(log=log)

    # TODO:
    #   * Get the corresponding output files for 
    #       each fsl_regfilt options.

    return outfile.file

def mcflirt(infile: str, 
            outfile: str, 
            reffile: Optional[str] = None, 
            spline_final: bool = True, 
            plots: bool = True,
            mats: bool = True, 
            refvol: Optional[str] = None,
            log: Optional[LogFile] = None
           ) -> str:
    """Rigid-body motion correction using ``mcflirt``."""
    infile: NiiFile = NiiFile(file=infile, assert_exists=True, validate_nifti=True)
    outfile: NiiFile = NiiFile(file=outfile)

    cmd: Command = Command("mcflirt")

    cmd.opt("-in"); cmd.opt(f"{infile.file}")
    cmd.opt("-out"); cmd.opt(f"{outfile.rm_ext()}")

    if reffile:
        reffile: File = File(file=reffile, assert_exists=True)
        cmd.opt("-reffile"); cmd.opt(f"{reffile.file}")
    
    if spline_final:
        cmd.opt("-spline_final")
    
    if plots:
        cmd.opt("-plots")
    
    if mats:
        cmd.opt("-mats")
    
    if refvol:
        refvol: NiiFile = NiiFile(file=refvol, assert_exists=True, validate_nifti=True)
        cmd.opt("-refvol"); cmd.opt(f"{refvol.file}")
    
    cmd.run(log=log)

    return None

def slicer(input: str,
           input2: Optional[str] = None,
           label: Optional[Union[int,str]] = None,
           lut: Optional[str] = None,
           intensity: Optional[Tuple[int]] = None,
           edgethreshold: Optional[int] = None,
           x: Optional[Tuple[int,str]] = None,
           y: Optional[Tuple[int,str]] = None,
           z: Optional[Tuple[int,str]] = None,
           log: Optional[LogFile] = None
          ) -> str:
    """Creates a combined NIFTI image using one or two NIFTI files.
    """
    # TODO: set output file(s) for this function
    input: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("slicer")
    cmd.opt(f"{input.file}")

    if input2:
        input2: NiiFile = NiiFile(file=input2, assert_exists=True, validate_nifti=True)
        cmd.opt(f"{input2.file}")
    
    if label:
        cmd.opt("-L")
        cmd.opt(f"{label}")
    
    if lut:
        lut: File = File(file=lut, assert_exists=True)
        cmd.opt("-l")
        cmd.opt(f"{lut.file}")
    
    if intensity:
        cmd.opt("-i")
        cmd.opt(f"{intensity[0]}")
        cmd.opt(f"{intensity[1]}")
    
    if edgethreshold:
        cmd.opt("-e")
        cmd.opt(f"{edgethreshold}")
    
    if x:
        cmd.opt("-x")
        cmd.opt(f"{x[0]}")
        cmd.opt(f"{x[1]}")

    if y:
        cmd.opt("-y")
        cmd.opt(f"{y[0]}")
        cmd.opt(f"{y[1]}")

    if z:
        cmd.opt("-z")
        cmd.opt(f"{z[0]}")
        cmd.opt(f"{z[1]}")
    
    cmd.run(log=log)

    return None

def cluster(infile: str,
            thresh: Optional[int] = None,
            oindex: Optional[str] = None,
            no_table: bool = False,
            log: Optional[LogFile] = None
           ):
    """Form clusters, report information about clusters and/or perform cluster-based inference.
    """
    infile: NiiFile = NiiFile(file=infile, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("cluster")
    cmd.opt(f"--in={infile.file}")
    
    if thresh:
        cmd.opt(f"--thresh={thresh}")
    
    if oindex:
        cmd.opt(f"--oindex={oindex}")
    
    if no_table:
        cmd.opt("--no_table")
    
    cmd.run(log=log)

    # TODO:
    #   * Get the corresponding output files for 
    #       each cluster options.
    #   * Get number and types of output files.

    return None


class fslmaths:
    """``FSL`` wrapper class for the ``fslmaths`` utility executable.
    Perform mathematical operations and/or manipulation of images.
    """

    def __init__(self,
                 img: str,
                 dt: Optional[str] = None,
                ):
        """Constructor"""
        img: NiiFile = NiiFile(file=img, assert_exists=True, validate_nifti=True)
        self._cmds: Command = Command("fslmaths")

        if dt:
            dt: str = FSLDataType(dt).name
            self._cmds.opt("-dt")
            self._cmds.opt(dt)
        
        self._cmds.opt(img.file)
    
    def abs(self):
        """work"""
        self._cmds.opt("-abs")
        return self
    
    def bin(self):
        """work"""
        self._cmds.opt("-bin")
        return self
    
    def binv(self):
        """work"""
        self._cmds.opt("-binv")
        return self
    
    def recip(self):
        """work"""
        self._cmds.opt("-recip")
        return self
    
    def Tmean(self):
        """work"""
        self._cmds.opt("-Tmean")
        return self
    
    def Tstd(self):
        """work"""
        self._cmds.opt("-Tstd")
        return self
    
    def Tmin(self):
        """work"""
        self._cmds.opt("-Tmin")
        return self
    
    def Tmax(self):
        """work"""
        self._cmds.opt("-Tmax")
        return self
    
    def sqrt(self):
        """work"""
        self._cmds.opt("-sqrt")
        return self
    
    def sqr(self):
        """work"""
        self._cmds.opt("-sqr")
        return self
    
    def fillh(self):
        """work"""
        self._cmds.opt("-fillh")
        return self
    
    def ero(self,
            repeat: int = 1
           ):
        """work"""
        for _ in range(repeat):
            self._cmds.opt("-ero")
        return self
    
    def dilM(self,
             repeat: int = 1
            ):
        """work"""
        for _ in range(repeat):
            self._cmds.opt("-dilM")
        return self
    
    def dilF(self,
             repeat: int = 1
            ):
        """work"""
        for _ in range(repeat):
            self._cmds.opt("-dilF")
        return self

    def add(self,
            input: Union[int,float,str]
           ):
        """work"""
        self._cmds.opt("-add")

        if isinstance(input,int) or isinstance(input,float):
            self._cmds.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.file)
        return self
    
    def sub(self,
            input: Union[int,float,str]
           ):
        """work"""
        self._cmds.opt("-sub")

        if isinstance(input,int) or isinstance(input,float):
            self._cmds.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.file)
        return self
    
    def mul(self,
            input: Union[int,float,str]
           ):
        """work"""
        self._cmds.opt("-mul")

        if isinstance(input,int) or isinstance(input,float):
            self._cmds.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.file)
        return self
    
    def div(self,
            input: Union[int,float,str]
           ):
        """work"""
        self._cmds.opt("-div")

        if isinstance(input,int) or isinstance(input,float):
            self._cmds.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.file)
        return self
    
    def mas(self,
            img: str
           ):
        """work"""
        img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
        self._cmds.opt("-mas")
        self._cmds.opt(img.file)
        return self

    def rem(self,
            input: Union[int,float,str]
           ):
        """work"""
        self._cmds.opt("-rem")

        if isinstance(input,int) or isinstance(input,float):
            self._cmds.opt(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.file)
        return self
    
    def thr(self,
            num: Union[int,float]
           ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._cmds.opt("-thr")
            self._cmds.opt(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self
    
    def uthr(self,
             num: Union[int,float]
            ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._cmds.opt("-uthr")
            self._cmds.opt(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self

    def inm(self,
            num: Union[int,float]
           ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._cmds.opt("-inm")
            self._cmds.opt(f"{num}")
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
        
        def _compute_hz(sec: float) -> float:
            """Computes frequency cutoff in hertz (Hz).
            """
            try:
                return 1/sec
            except ZeroDivisionError:
                return 0.0
            
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
            high_pass: float = _compute_sigma(tr=tr, 
                                              hz=_compute_hz(sec=high_pass), 
                                              compute_low_pass=False)
            low_pass:  float = _compute_sigma(tr=tr, 
                                              hz=_compute_hz(sec=low_pass),
                                              compute_low_pass=True)
        elif input_is_hz:
            high_pass: float = _compute_sigma(tr=tr, 
                                              hz=high_pass,
                                              compute_low_pass=False)
            low_pass:  float = _compute_sigma(tr=tr, 
                                              hz=low_pass,
                                              compute_low_pass=True)

        if ((isinstance(high_pass,int) or isinstance(high_pass,float)) and
                (isinstance(low_pass,int) or isinstance(low_pass,float))):
            self._cmds.opt("-bptf")
            self._cmds.opt(f"{high_pass}")
            self._cmds.opt(f"{low_pass}")
        else:
            raise TypeError(f"Input high_pass: {high_pass} or low_pass: {low_pass} is not a number.")
        
        return self
    
    def run(self,
            out: str,
            odt: Optional[str] = None,
            log: Optional[LogFile] = None
           ) -> str:
        """work"""
        self._cmds.opt(out)

        if odt:
            odt: str = FSLDataType(odt).name
            self._cmds.opt("-odt")
            self._cmds.opt(odt)
            
        self._cmds.run(log=log)
        return out