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
    descrip:        str = "descrip"
    intent_name:    str = "intent_name"

@unique
class MergeDim(Enum):
    """``fslmerge`` options."""
    t:  str = "t"
    x:  str = "x"
    y:  str = "y"
    z:  str = "z"
    a:  str = "a"
    tr: str = "tr"
    n:  str = "n"
    
@unique
class ECModelFLM(Enum):
    """``eddy`` first level EC model (``flm``) options."""
    movement:   str = "movement"
    linear:     str = "linear"
    quadratic:  str = "quadratic"   # default
    cubic:      str = "cubic"

@unique
class ECModelSLM(Enum):
    """``eddy`` second level EC model (``slm``) options."""
    none:       str = "none"        # default
    linear:     str = "linear"
    quadratic:  str = "quadratic"

@unique
class ECInterp(Enum):
    """``eddy`` interpolation model options."""
    spline:     str = "spline"      # default
    trilinear:  str = "trilinear"

@unique
class ECresamp(Enum):
    """``eddy`` image resampling methods."""
    jac: str = "jac"                # default
    lsr: str = "lsr"

@unique
class ECOLType(Enum):
    """``eddy`` outlier type methods."""
    sw:     str = "sw"              # default
    gw:     str = "gw"
    both:   str = "both"

@unique
class FSLDataType(Enum):
    """``FSL`` input and output data types.

    NOTE:
        * The 'input' option will set the datatype to that of the original image.
    """
    char:   str = "char"
    short:  str = "short"
    int:    str = "int"
    float:  str = "float"
    double: str = "double"
    input:  str = "input"

@unique
class RegInterp(Enum):
    """``FSL`` registration interpolation methods."""
    nn:         str = "nn"                  # applywarp exclusive option
    trilinear:  str = "trilinear"
    sinc:       str = "sinc"
    spline:     str = "spline"
    nearest:    str = "nearestneighbour"

@unique
class SliceAcqOrder(Enum):
    """Slice acquisition order method."""
    interleaved:    str = "interleaved"
    default:        str = "default"
    single_shot:    str = "single-shot"
