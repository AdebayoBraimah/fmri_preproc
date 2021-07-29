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
    
    if '.nii.gz' in out:
        out: str = out[:-7]
    elif '.nii' in out:
        out: str = out[:-4]
    
    # Required eddy current correction options
    img:    NiiFile = NiiFile(file=img)
    out:    NiiFile = NiiFile(file=out)
    mask:   NiiFile = NiiFile(file=mask)
    idx:    File = File(file=idx)
    acqp:   File = File(file=acqp)
    bvecs:  File = File(file=bvecs, ext='.bvec')
    bvals:  File = File(file=bvals, ext='.bval')

    eddy_proc.cmd_list.append(f"--imain={img.file}")
    eddy_proc.cmd_list.append(f"--mask={mask.file}")
    eddy_proc.cmd_list.append(f"--index={idx.file}")
    eddy_proc.cmd_list.append(f"--acqp={acqp.file}")
    eddy_proc.cmd_list.append(f"--bvecs={bvecs.file}")
    eddy_proc.cmd_list.append(f"--bvals={bvals.file}")
    eddy_proc.cmd_list.append(f"--out={out.file}")

    # Conventional eddy current correction options
    if topup:
        # Topup files' basenames
        topup: File = File(file=topup)
        eddy_proc.cmd_list.append(f"--topup={topup.rm_ext()}")
    
    if field:
        field: NiiFile = NiiFile(file=field)
        eddy_proc.cmd_list.append(f"--field={field.file}")
    
    if field_mat:
        field_mat: File = File(file=field_mat)
        eddy_proc.cmd_list.append(f"--field_mat={field_mat.file}")
    
    if flm:
        flm: str = ECModelFLM(flm).name
        eddy_proc.cmd_list.append(f"--flm={flm}")
    
    if slm:
        slm: str = ECModelSLM(slm).name
        eddy_proc.cmd_list.append(f"--slm={slm}")
    
    if fwhm:
        eddy_proc.cmd_list.append(f"--fwhm={fwhm}")
    
    if niter:
        eddy_proc.cmd_list.append(f"--niter={niter}")
    
    if cnr_maps:
        eddy_proc.cmd_list.append("--cnr_maps")
    
    if residuals:
        eddy_proc.cmd_list.append("--residuals")
    
    if fep:
        eddy_proc.cmd_list.append("--fep")
    
    if interp:
        interp: str = ECInterp(interp).name
        eddy_proc.cmd_list.append(f"--interp={interp}")
    
    if resamp:
        resamp: str = ECresamp(resamp).name
        eddy_proc.cmd_list.append(f"--resamp={resamp}")
    
    if nvoxhp:
        eddy_proc.cmd_list.append(f"--nvoxhp={nvoxhp}")
    
    if initrand:
        eddy_proc.cmd_list.append(f"--initrand={initrand}")
    
    if ff:
        eddy_proc.cmd_list.append(f"--ff={ff}")
    
    if repol:
        eddy_proc.cmd_list.append("--repol")
    
    if ol_nstd:
        eddy_proc.cmd_list.append(f"--ol_nstd={ol_nstd}")
    
    if ol_nvox:
        eddy_proc.cmd_list.append(f"--ol_nvox={ol_nvox}")

    if ol_type:
        ol_type: str = ECOLType(ol_type).name
        eddy_proc.cmd_list.append(f"--ol_type={ol_type}")

    if ol_pos:
        eddy_proc.cmd_list.append("--ol_pos")
    
    if ol_sqr:
        eddy_proc.cmd_list.append("--ol_sqr")
    
    if estimate_move_by_susceptibility:
        eddy_proc.cmd_list.append("--estimate_move_by_susceptibility")
    
    if mbs_niter:
        eddy_proc.cmd_list.append(f"--mbs_niter={mbs_niter}")
    
    if mbs_lambda:
        eddy_proc.cmd_list.append(f"--mbs_lambda={mbs_lambda}")
    
    if mbs_ksp:
        eddy_proc.cmd_list.append(f"--mbs_ksp={mbs_ksp}")
    
    if dont_sep_offs_move:
        eddy_proc.cmd_list.append("--dont_sep_offs_move")
    
    if dont_peas:
        eddy_proc.cmd_list.append("--dont_peas")
    
    if data_is_shelled:
        eddy_proc.cmd_list.append("--data_is_shelled")
    
    if b0_only:
        eddy_proc.cmd_list.append("--b0_only")
    
    if dont_mask_output:
        eddy_proc.cmd_list.append("--dont_mask_output")
    
    if verbose:
        eddy_proc.cmd_list.append("--verbose")
    
    if very_verbose:
        eddy_proc.cmd_list.append("--very_verbose")
    
    # Eddy slice-to-volume (s2v) correction options
    if mporder:
        eddy_proc.cmd_list.append(f"--mporder={mporder}")
    
        if slspec:
            slspec: File = File(file=slspec)
            eddy_proc.cmd_list.append(f"--slspec={slspec.file}")
        elif json_file:
            json_file: File = File(file=json_file)
            eddy_proc.cmd_list.append(f"--json={json_file.file}")
        else:
            raise FSLError("Either the 'slspec' or the 'json' option must be specified with the 'mporder' option.")
        
    if s2v_lambda:
        eddy_proc.cmd_list.append(f"--s2v_lambda={s2v_lambda}")
    
    if s2v_fwhm:
        eddy_proc.cmd_list.append(f"--s2v_fwhm={s2v_fwhm}")
    
    if s2v_niter:
        eddy_proc.cmd_list.append(f"--s2v_niter={s2v_niter}")
    
    if s2v_interp: 
        s2v_interp: str = ECInterp(s2v_interp).name
        eddy_proc.cmd_list.append(f"--s2v_interp={s2v_interp}")
    
    out: str = out.file + '.nii.gz'
    out: File = File(file=out)

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
    img: NiiFile = NiiFile(file=img)
    out: NiiFile = NiiFile(file=out)

    b: Command = Command("bet")
    b.cmd_list.append(img.abs_path())
    b.cmd_list.append(out.abs_path())

    if mask:
        b.cmd_list.append("-m")
        mask_img: str = out.rm_ext() + "_mask.nii.gz"
    else:
        mask_img: str = ""
    
    if robust:
        b.cmd_list.append("-R")

    if frac_int:
        b.cmd_list.append("-f")
        b.cmd_list.append(frac_int)
    
    if seg:
        pass
    else:
        b.cmd_list.append("-n")
    
    if verbose:
        b.cmd_list.append("-v")
    
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
    img: NiiFile = NiiFile(file=img)
    out: NiiFile = NiiFile(file=out)
    param: File = File(file=param)

    top: Command = Command("topup")
    top.cmd_list.append(f"--imain={img.abs_path()}")
    top.cmd_list.append(f"--out={out.abs_path()}")
    top.cmd_list.append(f"--datain={param.abs_path()}")
    top.cmd_list.append(f"--scale={scale}")

    if config:
        config: File = File(file=config)
        top.cmd_list.append(f"--config={config.abs_path()}")
    
    if fout:
        field_img: str = out.rm_ext() + "_field.nii.gz"
        field_img: File = File(file=field_img)
        top.cmd_list.append(f"--fout={field_img.rm_ext()}")
    else:
        field_img: str = ""
        field_img: File = File(file=field_img)
    
    if iout:
        unwarped_img: str = out.rm_ext() + "_unwarped.nii.gz"
        unwarped_img: File = File(file=unwarped_img)
        top.cmd_list.append(f"--iout={unwarped_img.rm_ext()}")
    else:
        unwarped_img: str = ""
        unwarped_img: File = File(file=unwarped_img)

    if verbose:
        top.cmd_list.append("--verbose")
    
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
    img: NiiFile = NiiFile(file=img)
    out: NiiFile = NiiFile(file=out)

    out_mat: str = out.rm_ext() + "_native2std_reorient.mat"
    out_mat: File = File(file=out_mat)

    reorient: Command = Command("fslreorient2std")

    if out_mat:
        reorient.cmd_list.append("-m")
        reorient.cmd_list.append(f"{out_mat.file}")
    else:
        out_mat: str = ""
        out_mat: File = File(file=out_mat)
    
    reorient.cmd_list.append(f"{img.file}")
    
    if out:
        reorient.cmd_list.append(f"{out.file}")
    else:
        out: str = ""
        out: File = File(file=out)

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
    img: NiiFile = NiiFile(file=img)
    out: NiiFile = NiiFile(file=out)

    if (xmin and xsize and ymin and ysize and zmin and zsize) and (tmin and tsize):
        raise FSLError("Cannot specify both xyz and temporal dimensions.")
    elif (not (xmin and xsize and ymin and ysize and zmin and zsize) and 
            not (tmin and tsize)):
        raise FSLError("Neither xyz nor temporal dimensions were specified.")
    else:
        roi: Command = Command("fslroi")
        roi.cmd_list.append(img.file)
        roi.cmd_list.append(out.file)
    
    if (xmin and xsize and ymin and ysize and zmin and zsize):
        roi.cmd_list.append(xmin); roi.cmd_list.append(xsize)
        roi.cmd_list.append(ymin); roi.cmd_list.append(ysize)
        roi.cmd_list.append(zmin); roi.cmd_list.append(zsize)
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
        roi.cmd_list.append(tmin); roi.cmd_list.append(tsize)
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
    merge.cmd_list.append(f"-{merge_opt}")
    merge.cmd_list.append(out.file)
    
    for arg in args:
        img: NiiFile = NiiFile(file=arg)
        merge.cmd_list.append(img.abs_path())
    
    if tr:
        merge.cmd_list.append(tr)
    
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

    src: NiiFile = NiiFile(file=src)
    ref: NiiFile = NiiFile(file=ref)
    out: NiiFile = NiiFile(file=out)

    apply: Command = Command("applywarp")

    apply.cmd_list.append(f"--in={src.abs_path()}")
    apply.cmd_list.append(f"--ref={ref.abs_path()}")
    apply.cmd_list.append(f"--out={out.file}")

    if interp:
        interp: str = RegInterp(interp).name
        apply.cmd_list.append(f"--interp={interp}")
    
    if prematdir:
        premat: str = os.path.join(prematdir,'allmats.txt')
        premat: str = catmats(matdir=prematdir, out=premat)

    if postmatdir:
        postmat: str = os.path.join(postmatdir,'allmats.txt')
        postmat: str = catmats(matdir=postmatdir, out=postmat)

    if warp:
        warp: NiiFile = NiiFile(file=warp)
        apply.cmd_list.append(f"--warp={warp}")

    if premat:
        apply.cmd_list.append(f"--premat={premat}")

    if postmat:
        apply.cmd_list.append(f"--postmat={postmat}")

    if paddingsize and (isinstance(paddingsize,int)):
        apply.cmd_list.append(f"--paddingsize={paddingsize}")

    if abs:
        apply.cmd_list.append("--abs")

    if rel:
        apply.cmd_list.append("--rel")

    apply.run(log=log)
    return out.file


def invxfm(inmat: str,
           outmat: str,
           log: Optional[LogFile] = None
          ) -> str:
    """Inverts ``FSL`` transformation matrices."""
    inv: Command = Command("invxfm")

    inv.cmd_list.append("-omat")
    inv.cmd_list.append(outmat)

    inv.cmd_list.append("-inverse")
    inv.cmd_list.append(inmat)

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
    src: NiiFile = NiiFile(file=src)
    ref: NiiFile = NiiFile(file=ref)
    out: NiiFile = NiiFile(file=out)
    mat: File = File(file=mat)

    xfm: Command = Command("flirt")

    xfm.cmd_list.append("-init")
    xfm.cmd_list.append(mat.file)

    xfm.cmd_list.append("-in")
    xfm.cmd_list.append(src.file)

    xfm.cmd_list.append("-ref")
    xfm.cmd_list.append(ref.file)

    xfm.cmd_list.append("-applyxfm")

    xfm.cmd_list.append("-out")
    xfm.cmd_list.append(out.file)

    interp: str = RegInterp(interp).name
    xfm.cmd_list.append("-interp")
    xfm.cmd_list.append(interp)

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
    src: NiiFile = NiiFile(file=src)
    ref: NiiFile = NiiFile(file=ref)
    out: NiiFile = NiiFile(file=out)

    fsldir: str = os.getenv('FSLDIR',None)
    assert fsldir is not None, 'FSLDIR environment must be set.'

    ident: str = os.path.join(fsldir,'etc/flirtsch/ident.mat')

    xfm: Command = Command("flirt")

    xfm.cmd_list.append("-init")
    xfm.cmd_list.append(ident)

    xfm.cmd_list.append("-in")
    xfm.cmd_list.append(src.file)

    xfm.cmd_list.append("-ref")
    xfm.cmd_list.append(ref.file)

    xfm.cmd_list.append("-applyisoxfm")
    xfm.cmd_list.append(f"{res}")

    xfm.cmd_list.append("-out")
    xfm.cmd_list.append(out.file)

    interp: str = RegInterp(interp).name
    xfm.cmd_list.append("-interp")
    xfm.cmd_list.append(interp)

    xfm.run(log=log)
    return out.file

def concatxfm(inmat1: str,
              inmat2: str,
              outmat: str,
              log: Optional[LogFile] = None
             ) -> str:
    """Concatenates two ``FSL`` transformation matrices."""
    inmat1: File = File(file=inmat1)
    inmat2: File = File(file=inmat2)

    concat: Command = Command("convert_xfm")

    concat.cmd_list.append("-omat")
    concat.cmd_list.append(outmat)

    concat.cmd_list.append("-concat")
    concat.cmd_list.append(inmat1.file)
    concat.cmd_list.append(inmat2.file)

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
    inwarp: NiiFile = NiiFile(file=inwarp)
    ref: NiiFile = NiiFile(file=ref)
    outwarp: NiiFile = NiiFile(file=outwarp)

    inv: Command = Command("invwarp")

    # Required arguments
    inv.cmd_list.append(f"--warp={inwarp.abs_path()}")
    inv.cmd_list.append(f"--ref={ref.abs_path()}")
    inv.cmd_list.append(f"--out={outwarp.file}")

    # Optional arguments
    if rel:
        inv.cmd_list.append("--rel")
    
    if abs:
        inv.cmd_list.append("--abs")
    
    if noconstraint:
        inv.cmd_list.append("--noconstraint")
    
    if jmin:
        inv.cmd_list.append(f"--jmin={jmin}")
    
    if jmax:
        inv.cmd_list.append(f"--jamx={jmax}")
    
    if verbose:
        inv.cmd_list.append("--verbose")
    
    inv.run(log=log)
    return outwarp.file

def convertwarp(out: str,
                ref: str,
                warp1: str,
                warp2: str,
                premat: Optional[str] = None,
                midmat: Optional[str] = None,
                postmat: Optional[str] = None,
                shiftmap: Optional[str] = None,
                shiftdir: Optional[str] = None,
                abs: bool = False,
                absout: bool = False,
                rel: bool = False,
                relout: bool = False
               ) -> None:
    """Convert ``FSL`` non-linear warps."""
    pass

def flirt():
    pass

def melodic():
    pass

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
        img: NiiFile = NiiFile(file=img)
        self._maths: Command = Command("fslmaths")

        if dt:
            dt: str = FSLDataType(dt).name
            self._maths.cmd_list.append("-dt")
            self._maths.cmd_list.append(dt)
        
        self._maths.cmd_list.append(img.file)
    
    def abs(self):
        """work"""
        self._maths.cmd_list.append("-abs")
        return self
    
    def bin(self):
        """work"""
        self._maths.cmd_list.append("-bin")
        return self
    
    def binv(self):
        """work"""
        self._maths.cmd_list.append("-binv")
        return self
    
    def recip(self):
        """work"""
        self._maths.cmd_list.append("-recip")
        return self
    
    def Tmean(self):
        """work"""
        self._maths.cmd_list.append("-Tmean")
        return self
    
    def Tstd(self):
        """work"""
        self._maths.cmd_list.append("-Tstd")
        return self
    
    def Tmin(self):
        """work"""
        self._maths.cmd_list.append("-Tmin")
        return self
    
    def Tmax(self):
        """work"""
        self._maths.cmd_list.append("-Tmax")
        return self
    
    def sqrt(self):
        """work"""
        self._maths.cmd_list.append("-sqrt")
        return self
    
    def sqr(self):
        """work"""
        self._maths.cmd_list.append("-sqr")
        return self
    
    def fillh(self):
        """work"""
        self._maths.cmd_list.append("-fillh")
        return self
    
    def ero(self,
            repeat: int = 1
           ):
        """work"""
        for _ in range(repeat):
            self._maths.cmd_list.append("-ero")
        return self
    
    def dilM(self,
             repeat: int = 1
            ):
        """work"""
        for _ in range(repeat):
            self._maths.cmd_list.append("-dilM")
        return self
    
    def dilF(self,
             repeat: int = 1
            ):
        """work"""
        for _ in range(repeat):
            self._maths.cmd_list.append("-dilF")
        return self

    def add(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.cmd_list.append("-add")

        if isinstance(input,int):
            self._maths.cmd_list.append(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input)
            self._maths.cmd_list.append(img.file)
        return self
    
    def sub(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.cmd_list.append("-sub")

        if isinstance(input,int):
            self._maths.cmd_list.append(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input)
            self._maths.cmd_list.append(img.file)
        return self
    
    def mul(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.cmd_list.append("-mul")

        if isinstance(input,int):
            self._maths.cmd_list.append(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input)
            self._maths.cmd_list.append(img.file)
        return self
    
    def div(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.cmd_list.append("-div")

        if isinstance(input,int):
            self._maths.cmd_list.append(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input)
            self._maths.cmd_list.append(img.file)
        return self
    
    def mas(self,
            img: str
           ):
        """work"""
        img: NiiFile = NiiFile(file=img)
        self._maths.cmd_list.append("-mas")
        self._maths.cmd_list.append(img.file)
        return self

    def rem(self,
            input: Union[int,str]
           ):
        """work"""
        self._maths.cmd_list.append("-rem")

        if isinstance(input,int):
            self._maths.cmd_list.append(f"{input}")
        elif isinstance(input,str):
            img: NiiFile = NiiFile(file=input)
            self._maths.cmd_list.append(img.file)
        return self
    
    def thr(self,
            num: Union[int,float]
           ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._maths.cmd_list.append("-thr")
            self._maths.cmd_list.append(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self
    
    def uthr(self,
             num: Union[int,float]
            ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._maths.cmd_list.append("-uthr")
            self._maths.cmd_list.append(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self

    def inm(self,
            num: Union[int,float]
           ):
        """work"""
        if isinstance(num,int) or isinstance(num,float):
            self._maths.cmd_list.append("-inm")
            self._maths.cmd_list.append(f"{num}")
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
        Input is assumed to in sigma (in volume).
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
            self._maths.cmd_list.append("-bptf")
            self._maths.cmd_list.append(f"{high_pass}")
            self._maths.cmd_list.append(f"{low_pass}")
        else:
            raise TypeError(f"Input high_pass: {high_pass} or low_pass: {low_pass} is not a number.")
        
        return self
    
    def run(self,
            out: str,
            odt: Optional[str] = None,
            log: Optional[LogFile] = None
           ) -> str:
        """work"""
        self._maths.cmd_list.append(out)

        if odt:
            odt: str = FSLDataType(odt).name
            self._maths.cmd_list.append("-odt")
            self._maths.cmd_list.append(odt)
            
        self._maths.run(log=log)
        return out


