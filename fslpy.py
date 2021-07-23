# -*- coding: utf-8 -*-
"""Wrappers functions for FSL executable binaries.
"""
from typing import (
    Dict,
    List,
    Optional,
    Tuple
)
from util import (
    File,
    WorkDir,
    TmpDir,
    NiiFile,
    LogFile,
    Command
)

def fnirt():
    """work"""
    pass

def eddy():
    """work"""
    pass

def bet():
    """work"""
    pass

def topup(img: str,
          param: str,
          out: str,
          config: Optional[str] = None,
          fout: Optional[str] = None,
          iout: Optional[str] = None,
          scale: int = 1,
          verbose: bool = False,
          log: Optional[LogFile] = None
         ) -> str:
    """work"""
    top: Command = Command("topup")
    top.cmd_list.append(f"--imain={img}")
    top.cmd_list.append(f"--out={out}")
    top.cmd_list.append(f"--datain={param}")
    top.cmd_list.append(f"--scale={scale}")

    if config:
        top.cmd_list.append(f"--config={config}")
    
    if fout:
        top.cmd_list.append(f"--fout={fout}")
    
    if iout:
        top.cmd_list.append(f"--iout={iout}")

    if verbose:
        top.cmd_list.append("--verbose")
    
    top.run(log=log)
    return out

def fslreorient2std(img: str,
                    out: Optional[str] = None,
                    out_mat: Optional[str] = None,
                    log: Optional[LogFile] = None
                   ) -> Tuple[str,str,str]:
    """work"""
    reorient: Command = Command("fslreorient2std")

    if out_mat:
        reorient.cmd_list.append("-m")
        reorient.cmd_list.append(f"{out_mat}")
    
    reorient.cmd_list.append(f"{img}")
    
    if out:
        reorient.cmd_list.append(f"{out}")

    reorient.run(log=log)

    return (img, 
            out, 
            out_mat)

def fslroi():
    """work"""
    pass

def fslmerge():
    """work"""
    pass

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

