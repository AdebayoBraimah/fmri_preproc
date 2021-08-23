# -*- coding: utf-8 -*-
"""Setup up acquisition parameters needed for EDDY, TOPUP, and FLIRT.
"""
import nibabel as nib
import numpy as np
import json

from typing import (
    Dict,
    List,
    Optional,
    Tuple,
    Union
)

from fmri_preproc.utils.enums import PhaseEncodeDirection

def write_func_params(epi: str,
                      echospacing: float,
                      pedir: Union[str,List[str]],
                      out: str,
                      epifactor: Optional[int] = None,
                      inplane_acc: Optional[float] = 1,
                      out_flirt: Optional[str] = None
                     ) -> Tuple[str,Union[str,None]]:
    """doc
    """
    if not isinstance(pedir, list):
        pedir: List[str] = pedir.split(",")

    # Verify input phase encoding directions
    pedir: List[str] = [ PhaseEncodeDirection(x.upper()).name for x in pedir ]

    eddyp: np.array = np.zeros((len(pedir),4))
    ax: List[str] = [''] * len(pedir)

    echospacing: float = float(echospacing)
    inplane_acc: float = float(inplane_acc)

    epi: nib.Nifti1Image = nib.load(epi)
    
    assert len(pedir) == epi.header.get('dim','')[4]

    # Loop through phase encoding directions and determine:
    #   * axis orientation
    #   * epi factor
    #   * timing for EDDY/TOPUP

    for i in range(len(pedir)):

        (ax[i], 
         axidx, 
         axdir) = _get_axis(epi=epi, pedir=pedir[i])
        
        eddyp[i, axidx] = axdir

        if epifactor:
            epif: int = int(epifactor)
        else:
            epif: int = epi.shape[axidx]
        
        eddyp[i, 3] = (epif - 1) * echospacing / inplane_acc
    
    # Write FLIRT parameters
    if out_flirt:
        param: Dict[str,str] = {'pedir': ax, 'echospacing': echospacing}
        out_flirt: str = _dict2json(param, out_flirt)
    else:
        out_flirt: str = None
    
    # Write EDDY/TOPUP parameters
    np.savetxt(out, eddyp, fmt='%i %i %i %f')

    return out, out_flirt

def _dict2json(dict: Dict, 
               jsonfile: str, 
               indent: int = 4
              ) -> str:
    """Write dictionary to json file."""
    with open(jsonfile, 'w') as outfile:
        json.dump(dict, outfile, indent=indent)
    return jsonfile

def _get_axis(epi: Union[str, nib.Nifti1Image],
              pedir: str,
             ) -> Tuple[str,str,str]:
    """doc
    """
    if not isinstance(epi, nib.Nifti1Image):
        epi: nib.Nifti1Image = nib.load(epi)
    ornt: nib.orientations = _get_orientation(epi=epi)

    if ornt[0] == pedir:
        ax = 'x'
        axidx = 0
        axdir = 1
    elif _flip(ornt[0]) == pedir:
        ax = 'x-'
        axidx = 0
        axdir = -1
    elif ornt[1] == pedir:
        ax = 'y'
        axidx = 1
        axdir = 1
    elif _flip(ornt[1]) == pedir:
        ax = 'y-'
        axidx = 1
        axdir = -1
    elif ornt[2] == pedir:
        ax = 'z'
        axidx = 2
        axdir = 1
    elif _flip(ornt[2]) == pedir:
        ax = 'z-'
        axidx = 2
        axdir = -1
    
    return ax, axidx, axdir

def _flip(d: str) -> str:
    """doc
    """
    return d[::-1]

def _get_orientation(epi: nib.Nifti1Image) -> nib.orientations:
    """doc
    """
    labels: Tuple[str,str] = (('RL', 'LR'), ('AP', 'PA'), ('SI', 'IS'))
    code: nib.orientations = nib.orientations.aff2axcodes(aff=epi.affine, labels=labels)
    return code

def phase_encode_dict() -> Dict[str,Tuple[int,str]]:
    """doc
    """
    pedict: Dict[str,Tuple[int,str]] = {
        'x': (1, 'x'),
        'y': (2, 'y'),
        'z': (3, 'z'),
        '-x': (-1, 'x-'),
        '-y': (-2, 'y-'),
        '-z': (-3, 'z-'),
        'x-': (-1, 'x-'),
        'y-': (-2, 'y-'),
        'z-': (-3, 'z-')
    }
    return pedict