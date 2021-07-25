# -*- coding: utf-8 -*-
"""Wrappers functions for FSL executable binaries.
"""
from typing import (
    Dict,
    List,
    Optional,
    Tuple
)

from numpy import log

from util import (
    File,
    WorkDir,
    TmpDir,
    NiiFile,
    LogFile,
    Command
)

from enums import MergeDim

class FSLError(Exception):
    """Exception intended to be raised for FSL specific binaries and related wrapper functions."""
    pass

def fnirt():
    """work"""
    pass

def eddy():
    """work"""
    pass

def bet(nii: str,
        out: str,
        mask: bool = False,
        robust: bool = False,
        seg: bool = True,
        frac_int: Optional[int] = None,
        verbose: bool = False,
        log: Optional[LogFile] = None
       ) -> Tuple[str]:
    """work"""
    nii: NiiFile = NiiFile(file=nii)
    out: NiiFile = NiiFile(file=out)

    b: Command = Command("bet")
    b.cmd_list.append(nii.abs_path())
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

def fslmaths():
    pass

