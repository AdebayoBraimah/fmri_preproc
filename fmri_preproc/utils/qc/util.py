# -*- coding: utf-8 -*-
"""Quality control utility functions for the ``fmri_preproc`` neonatal rs-fMRI preprocessing pipeline.

.. autosummary::
    :nosignatures:

    dict2json
    json2dict
    run
    load_img
    bet
    split
    remove_ext
"""
import os
import json
import subprocess
import sys
import nibabel as nib
import numpy as np

from collections import OrderedDict
from tempfile import TemporaryDirectory

from typing import Any, Dict, Optional, Tuple


def dict2json(dict: Dict[Any, Any], jsonfile: str, indent: int = 4) -> None:
    """Write dictionary to json file."""
    for k, v in dict.items():

        if isinstance(v, np.ndarray):
            dict[k] = v.tolist()
        elif isinstance(v, nib.Nifti1Image):
            dict[k] = v.get_filename()

    with open(jsonfile, "w") as outfile:
        json.dump(dict, outfile, indent=indent)
    return None


def json2dict(jsonfile: str) -> Dict[Any, Any]:
    """Read dictionary from json file."""
    with open(jsonfile, "r") as infile:
        d = json.load(infile, object_pairs_hook=OrderedDict)
    return d


def run(cmd):
    """Executes command line command."""
    for i, v in enumerate(cmd):
        if isinstance(v, nib.Nifti1Image):
            cmd[i] = v.get_filename()
    print("RUNNING: " + " ".join(cmd))
    sys.stdout.flush()
    jobout = subprocess.check_output(cmd)
    return jobout.decode("utf-8").strip()


def load_img(d: str) -> nib.Nifti1Image:
    """Loads NIFTI-1 or 2 image as ``nibabel`` image object."""
    return nib.load(os.path.expanduser(d)) if isinstance(d, str) else d


def bet(
    func: nib.Nifti1Image, maskname: str, tmpdir: Optional[str] = None
) -> nib.Nifti1Image:
    """Performs brain extraction of the input ``nibabel`` image object.
    """
    with TemporaryDirectory(dir=tmpdir) as td:
        run(["fslmaths", func.get_filename(), "-Tmean", os.path.join(td, "mean")])
        run(["bet", os.path.join(td, "mean"), os.path.join(td, "brain"), "-n", "-m"])
        nib.load(os.path.join(td, "brain_mask.nii.gz")).to_filename(maskname)
        maskname: nib.Nifti1Image = nib.load(maskname)
    return maskname


def split(fname: str) -> Tuple[str, str, str]:
    """Split the filename into path, basename, and extension.

    Can handle zipped nifti (*.nii.gz)
    """
    if fname.endswith(".nii.gz"):
        return (
            os.path.dirname(fname),
            os.path.basename(os.path.splitext(os.path.splitext(fname)[0])[0]),
            ".nii.gz",
        )
    else:
        return (
            os.path.dirname(fname),
            os.path.basename(os.path.splitext(fname)[0]),
            os.path.splitext(fname)[1],
        )


def remove_ext(fname: str) -> str:
    """Remove the extension from the supplied filename.

    Can handle zipped nifti (*.nii.gz)
    """
    return os.path.join(split(fname)[0], split(fname)[1])
