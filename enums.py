# -*- coding: utf-8 -*-
"""Enums modules for ``fmri_preproc``.
"""
from enum import (
    Enum,
    unique
)

@unique
class NiiHeaderField(Enum):
    """NIFTI file header field options."""
    descrip: int = 1
    intent_name: int = 2

@unique
class MergeDim(Enum):
    """``fslmerge`` options."""
    t: int = 1
    x: int = 2
    y: int = 3
    z: int = 4
    a: int = 5
    tr: int = 6
    n: int = 7