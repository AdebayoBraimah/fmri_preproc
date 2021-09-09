# -*- coding: utf-8 -*-
"""``fmri_preproc`` neonatal structural masks setup.
"""
import numpy as np
import nibabel as nib

from typing import (
    Dict,
    List
)

from fmri_preproc.utils.fileio import NiiFile
from fmri_preproc.utils.enums import SegType


def get_dseg_labels(seg_type: str) -> Dict[str,List[int]]:
    """Gets image segmentation labels.
    """
    seg_type: str = SegType(seg_type.lower()).name

    if seg_type == 'drawem':
        labels: Dict[str,int] = {
                                    'csf': [1, 5],
                                    'gm': [2],
                                    'wm': [3],
                                    'sc': [7, 9],
                                    'bs': [8],
                                    'cb': [6]
                                }
    elif seg_type == 'freesurfer_aseg':
        labels: Dict[str,int] = {
                                    'gm': [3, 42],
                                    'wm': [2, 41, 77, 78, 79, 192, 250, 251, 252, 253, 254, 255],
                                    'sc': [9, 10, 11, 12, 13, 17, 18, 19, 20, 26, 27, 28, 31, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 59, 60, 63],
                                    'bs': [16],
                                    'cb': [7, 8, 46, 47]
                                }
    elif seg_type == 'fsl_fast':
        labels: Dict[str,int] = {
                                    'csf': [1],
                                    'gm': [2],
                                    'wm': [3],
                                    'sc': [2],
                                    'bs': [3],
                                    'cb': [3]
                                }
    else:
        raise RuntimeError('Unknown segmentation type.')
    return labels


def convert_to_fsl_fast(dseg: str,
                        seg_type: str,
                        out: str
                       ) -> str:
    """Converts segmentation types to FSL ``FAST`` segmentation files.
    """
    with NiiFile(src=dseg, assert_exists=True, validate_nifti=True) as dg:
        dseg: str = dg.abspath()

    dseg: nib.Nifti1Image = nib.load(dseg)
    d: np.array = dseg.get_data()
    d0: np.array = np.zeros(d.shape)

    label_map: Dict[str,List[int]] = get_dseg_labels(seg_type=seg_type)

    for name, idx in get_dseg_labels(seg_type="fsl_fast").items():
        if name not in label_map:
            continue
        else:
            for i in label_map[name]:
                d0[d == i] = idx
    
    nii: nib.Nifti1Image = nib.Nifti1Image(d0, dseg.affine, dseg.header)
    nii.to_filename(out)
    with NiiFile(src=out, assert_exists=True, validate_nifti=True) as ot:
        out: str = ot.abspath()
    return out


def create_mask(dseg: str,
                dseg_type: str,
                labels: List[str],
                out: str
               ) -> str:
    """Create segmentation mask.
    """
    dseg: nib.Nifti1Image = nib.load(dseg)
    d: np.array = dseg.get_data()
    d0: np.array = np.zeros(d.shape)

    label_map: Dict[str,List[int]] = get_dseg_labels(seg_type=dseg_type)

    for i in labels:
        for j in label_map[i]:
            d0[d == j] = 1
    
    nii: nib.Nifti1Image = nib.Nifti1Image(d0, dseg.affine, dseg.header)
    nii.to_filename(out)
    with NiiFile(src=out, assert_exists=True, validate_nifti=True) as ot:
        out: str = ot.abspath()
    return out
