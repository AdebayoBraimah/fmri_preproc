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
    descrip:        int = 1
    intent_name:    int = 2

@unique
class MergeDim(Enum):
    """``fslmerge`` options."""
    t:  int = 1
    x:  int = 2
    y:  int = 3
    z:  int = 4
    a:  int = 5
    tr: int = 6
    n:  int = 7
    
@unique
class ECModelFLM(Enum):
    """``eddy`` first level EC model (``flm``) options."""
    movement:   int = 1
    linear:     int = 2
    quadratic:  int = 3     # default
    cubic:      int = 4

@unique
class ECModelSLM(Enum):
    """``eddy`` second level EC model (``slm``) options."""
    none:       int = 1     # default
    linear:     int = 2
    quadratic:  int = 3

@unique
class ECInterp(Enum):
    """``eddy`` interpolation model options."""
    spline:     int = 1     # default
    trilinear:  int = 2

@unique
class ECresamp(Enum):
    """``eddy`` image resampling methods."""
    jac: int = 1            # default
    lsr: int = 2

@unique
class ECOLType(Enum):
    """``eddy`` outlier type methods."""
    sw:     int = 1         # default
    gw:     int = 2
    both:   int = 3