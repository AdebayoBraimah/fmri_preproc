# -*- coding: utf-8 -*-
"""Wrappers functions for FSL executable binaries.
"""
from typing import (
    List,
    Optional,
    Tuple,
    Union
)

from numpy import log

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
    ECOLType
)

class FSLError(Exception):
    """Exception intended to be raised for FSL specific binaries and related wrapper functions."""
    pass

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
         flm: ECModelFLM = ECModelFLM.quadratic,
         slm: ECModelSLM = ECModelSLM.none,
         fwhm: Union[int,str] = 0,
         s2v_fwhm: Union[int,str] = 0,
         niter: Union[int,str] = 5,
         s2v_niter: int = 5,
         cnr_maps: bool = False,
         residuals: bool = False,
         fep: bool = False,
         interp: ECInterp = ECInterp.spline,
         s2v_interp: ECInterp = ECInterp.trilinear,
         resamp: ECresamp = ECresamp.jac,
         nvoxhp: int = 1000,
         initrand: int = 0,
         ff: float = 10.0,
         repol: bool = False,
         ol_nstd: int = 4,
         ol_nvox: int = 250,
         ol_type: ECOLType = ECOLType.sw,
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
        field: File = File(file=field)
        eddy_proc.cmd_list.append(f"--field={field.file}")
    
    if field_mat:
        field_mat: File = File(file=field_mat)
        eddy_proc.cmd_list.append(f"--field_mat={field_mat.file}")
    
    if flm == ECModelFLM.movement:
        eddy_proc.cmd_list.append("--flm=movement")
    elif flm == ECModelFLM.linear:
        eddy_proc.cmd_list.append("--flm=linear")
    elif flm == ECModelFLM.quadratic:
        eddy_proc.cmd_list.append("--flm=quadratic")
    elif flm == ECModelFLM.cubic:
        eddy_proc.cmd_list.append("--flm=cubic")
    
    if slm == ECModelSLM.none:
        eddy_proc.cmd_list.append("--slm=none")
    elif slm == ECModelSLM.linear:
        eddy_proc.cmd_list.append("--slm=linear")
    elif slm == ECModelSLM.quadratic:
        eddy_proc.cmd_list.append("--slm=quadratic")
    
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
    
    if interp == ECInterp.spline:
        eddy_proc.cmd_list.append("--interp=spline")
    elif interp == ECInterp.trilinear:
        eddy_proc.cmd_list.append("--interp=trilinear")
    
    if resamp == ECresamp.jac:
        eddy_proc.cmd_list.append("--resamp=jac")
    elif resamp == ECresamp.lsr:
        eddy_proc.cmd_list.append("--resamp=lsr")
    
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

    if ol_type == ECOLType.sw:
        eddy_proc.cmd_list.append("--ol_type=sw")
    elif ol_type == ECOLType.gw:
        eddy_proc.cmd_list.append("--ol_type=gw")
    elif ol_type == ECOLType.both:
        eddy_proc.cmd_list.append("--ol_type=both")

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
    
    if s2v_interp == ECInterp.spline:
        eddy_proc.cmd_list.append("--s2v_interp=spline")
    elif s2v_interp == ECInterp.trilinear:
        eddy_proc.cmd_list.append("--s2v_interp=trilinear")
    
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
    else:
        roi: Command = Command("fslroi")
        roi.cmd_list.append(img.file)
        roi.cmd_list.append(out.file)
    
    if (xmin and xsize and ymin and ysize and zmin and zsize):
        roi.cmd_list.append(xmin); roi.cmd_list.append(xsize)
        roi.cmd_list.append(ymin); roi.cmd_list.append(ysize)
        roi.cmd_list.append(zmin); roi.cmd_list.append(zsize)
    elif (xmin is None) and \
        (xsize is None) and \
        (ymin is None) and \
        (ysize is None) and \
        (zmin is None) and \
        (zsize is None):
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
             merge_opt: MergeDim = MergeDim.t,
             tr: Optional[float] = None,
             log: Optional[LogFile] = None,
             *args
            ) -> str:
    """work"""
    out: NiiFile = NiiFile(file=out)

    merge: Command = Command("fslmerge")

    if merge_opt == MergeDim.t:
        merge.cmd_list.append("-t")
    elif merge_opt == MergeDim.x:
        merge.cmd_list.append("-x")
    elif merge_opt == MergeDim.y:
        merge.cmd_list.append("-y")
    elif merge_opt == MergeDim.z:
        merge.cmd_list.append("-z")
    elif merge_opt == MergeDim.a:
        merge.cmd_list.append("-a")
    elif merge_opt == MergeDim.tr:
        merge.cmd_list.append("-tr")
    elif merge_opt == MergeDim.n:
        merge.cmd_list.append("-n")
    
    merge.cmd_list.append(out.file)
    
    for arg in args:
        img: NiiFile = NiiFile(file=arg)
        merge.cmd_list.append(img.abs_path())
    
    if tr:
        merge.cmd_list.append(tr)
    
    merge.run(log=log)

    return out.file

def catmats():
    """work"""
    pass

def applywarp():
    """work"""
    pass

def invxfm():
    """work"""
    pass

def applyxfm():
    """work"""
    pass

def apply_isoxfm():
    """work"""
    pass

def concatxfm():
    pass

def invwarp():
    pass

def convertwarp():
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
                 img: str):
        """Constructor"""
        # TODO: Add options for input data-type '-dt'
        #   and output data-type '-odt' options
        self.maths: Command = Command("fslmaths")
        self.maths.cmd_list.append(img)
    
    def abs(self):
        """work"""
        self.maths.cmd_list.append("-abs")
        return self
    
    def bin(self):
        """work"""
        self.maths.cmd_list.append("-bin")
        return self
    
    def binv(self):
        """work"""
        self.maths.cmd_list.append("-binv")
        return self
    
    def recip(self):
        """work"""
        self.maths.cmd_list.append("-recip")
        return self
    
    def Tmean(self):
        """work"""
        self.maths.cmd_list.append("-Tmean")
        return self
    
    def Tstd(self):
        """work"""
        self.maths.cmd_list.append("-Tstd")
        return self
    
    def Tmin(self):
        """work"""
        self.maths.cmd_list.append("-Tmin")
        return self
    
    def Tmax(self):
        """work"""
        self.maths.cmd_list.append("-Tmax")
        return self
    
    def fillh(self):
        """work"""
        self.maths.cmd_list.append("-fillh")
        return self
    
    def ero(self,
            repeat: int = 1
           ):
        """work"""
        for _ in range(repeat):
            self.maths.cmd_list.append("-ero")
        return self
    
    def dilM(self,
             repeat: int = 1
            ):
        """work"""
        for _ in range(repeat):
            self.maths.cmd_list.append("-dilM")
        return self
    
    def dilF(self,
             repeat: int = 1
            ):
        """work"""
        for _ in range(repeat):
            self.maths.cmd_list.append("-dilF")
        return self

    # Add more fslmaths options below this comment
    
    def run(self,
            out: str,
            log: Optional[LogFile] = None
           ) -> str:
        """work"""
        self.maths.cmd_list.append(out)
        self.maths.run(log=log)
        return out


