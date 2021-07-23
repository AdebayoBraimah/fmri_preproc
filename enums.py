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