# -*- coding: utf-8 -*-
"""Wrappers functions for FSL executable binaries.

.. autosummary::
    :nosignatures:

    fnirt
    eddy
    bet
    topup
    fslreorient2std
    fslroi
    sigloss
    fslmerge
    catmats
    applywarp
    invxfm
    applyxfm
    apply_isoxfm
    concatxfm
    invwarp
    convertwarp
    fugue
    flirt
    melodic
    fsl_regfilt
    mcflirt
    slicer
    cluster

.. autoclass:: fslmaths
    :members:

.. autoclass:: FSLError
    :members:
"""
# TODO:
#
#   * Write doc-strings for each wrapper function
#   * Add verbose options to all wrapper functions
#   * Make FSLDIR global variable
#       * Add the option prepend FSLDIR/bin to each command object

import glob
import os
import numpy as np

from typing import List, Optional, Tuple, Union

from fmri_preproc.utils.fileio import File, NiiFile

from fmri_preproc.utils.logutil import LogFile

from fmri_preproc.utils.command import Command, DependencyError

from fmri_preproc.utils.enums import (
    ECModelFLM,
    ECModelSLM,
    MergeDim,
    ECInterp,
    ECresamp,
    ECOLType,
    FSLDataType,
    RegInterp,
)


# Global FSLDIR variable
FSLDIR: str = os.getenv("FSLDIR", None)


class FSLError(Exception):
    """Exception intended to be raised for FSL specific binaries and related wrapper functions."""

    pass


def fnirt(
    src: str,
    ref: str,
    aff: Optional[str] = None,
    imprefm: Optional[int] = None,
    impinm: Optional[int] = None,
    applyrefmask: Optional[str] = None,
    applyinmask: Optional[str] = None,
    subsamp: Optional[str] = None,
    miter: Optional[str] = None,
    infwhm: Optional[str] = None,
    reffwhm: Optional[str] = None,
    lmbda: Optional[str] = None,
    estint: Optional[str] = None,
    warpres: Optional[str] = None,
    ssqlambda: Optional[str] = None,
    regmod: Optional[str] = None,
    intmod: Optional[str] = None,
    intorder: Optional[str] = None,
    biasres: Optional[str] = None,
    biaslambda: Optional[str] = None,
    refderiv: Optional[str] = None,
    cout: Optional[str] = None,
    intout: Optional[str] = None,
    refout: Optional[str] = None,
    iout: Optional[str] = None,
    interp: Optional[str] = None,
    inwarp: Optional[str] = None,
    minmet: Optional[str] = None,
    verbose: bool = False,
    intin: Optional[str] = None,
    jout: Optional[str] = None,
    log: Optional[LogFile] = None,
) -> str:
    """Perform/compute non-linear image registrations."""
    src: NiiFile = NiiFile(src=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(src=ref, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("fnirt")

    # TODO: Verify that input options are either:
    #   * str
    #   * int or float

    if aff:
        aff: File = File(src=aff, assert_exists=True)
        cmd.opt(f"--aff={aff.src}")

    if imprefm:
        imprefm: File = File(src=imprefm, assert_exists=True)
        cmd.opt(f"--imprefm={imprefm.src}")

    if impinm:
        impinm: File = File(src=impinm, assert_exists=True)
        cmd.opt(f"--impinm={impinm.src}")

    if applyrefmask:
        applyrefmask: File = File(src=applyrefmask, assert_exists=True)
        cmd.opt(f"--applyrefmask={applyrefmask.src}")

    if applyinmask:
        applyinmask: File = File(src=applyinmask, assert_exists=True)
        cmd.opt(f"--applyinmask={applyinmask.src}")

    if subsamp:
        pass

    if miter:
        pass

    if infwhm:
        pass

    if reffwhm:
        pass

    if lmbda:
        pass

    if estint:
        pass

    if warpres:
        pass

    if ssqlambda:
        pass

    if regmod:
        pass

    if intmod:
        pass

    if intorder:
        pass

    if biasres:
        pass

    if biaslambda:
        pass

    if refderiv:
        pass

    if cout:
        pass

    if intout:
        pass

    if refout:
        pass

    if iout:
        pass

    if interp:
        pass

    if inwarp:
        pass

    if minmet:
        pass

    if verbose:
        pass

    if intin:
        pass

    if jout:
        pass

    cmd.run(log=log)

    return None


def eddy(
    img: str,
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
    flm: str = "quadratic",
    slm: str = "none",
    fwhm: Union[int, str] = 0,
    s2v_fwhm: Union[int, str] = 0,
    niter: Union[int, str] = 5,
    s2v_niter: int = 5,
    cnr_maps: bool = False,
    residuals: bool = False,
    fep: bool = False,
    interp: str = "spline",
    s2v_interp: str = "trilinear",
    resamp: str = "jac",
    nvoxhp: int = 1000,
    initrand: int = 0,
    ff: float = 10.0,
    repol: bool = False,
    ol_nstd: int = 4,
    ol_nvox: int = 250,
    ol_type: str = "sw",
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
    log: Optional[LogFile] = None,
) -> Tuple[str, str, str]:
    """Performs eddy current and distortion correction of 4D MR images.

    NOTE: 
        The function return files are not exhaustive (i.e. not all of eddy's files are returned, just the necessary ones).
    """
    eddy_cmds: List[str] = [
        "eddy_cuda",
        "eddy_cuda9.1",
        "eddy_cuda8.0",
        "eddy_cuda7.5",
        "eddy_openmp",
        "eddy",
    ]

    # Select for most optimal implementation of eddy.
    #   Look for GPU implementations first, followed
    #   by openmp, then eddy.
    for eddy_cmd in eddy_cmds:
        try:
            if "cuda" in eddy_cmd:
                cuda_cmd: str = "nvcc"
            else:
                cuda_cmd: str = "echo"

            cmd1: Command = Command(eddy_cmd)
            cmd2: Command = Command(cuda_cmd)

            # Check if both of the dependencies are met.
            if (cmd1.check_dependency()) and (cmd2.check_dependency()):
                cmd: Command = cmd1
                break
        except (DependencyError, FileNotFoundError):
            continue

    # Required eddy current correction options
    img: NiiFile = NiiFile(src=img, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(src=out, assert_exists=False, validate_nifti=False)
    mask: NiiFile = NiiFile(src=mask, assert_exists=True, validate_nifti=True)
    idx: File = File(src=idx, ext=".idx", assert_exists=True)
    acqp: File = File(src=acqp, ext=".acqp", assert_exists=True)
    bvecs: File = File(src=bvecs, ext=".bvec", assert_exists=True)
    bvals: File = File(src=bvals, ext=".bval", assert_exists=True)

    cmd.opt(f"--imain={img.src}")
    cmd.opt(f"--mask={mask.src}")
    cmd.opt(f"--index={idx.src}")
    cmd.opt(f"--acqp={acqp.src}")
    cmd.opt(f"--bvecs={bvecs.src}")
    cmd.opt(f"--bvals={bvals.src}")
    cmd.opt(f"--out={out.rm_ext()}")

    # Conventional eddy current correction options
    if topup:
        topup: NiiFile = NiiFile(src=topup, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--topup={topup.rm_ext()}")

    if field:
        field: NiiFile = NiiFile(src=field, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--field={field.rm_ext()}")

    if field_mat:
        field_mat: File = File(src=field_mat, assert_exists=True)
        cmd.opt(f"--field_mat={field_mat.src}")

    if flm:
        flm: str = ECModelFLM(flm).name
        cmd.opt(f"--flm={flm}")

    if slm:
        slm: str = ECModelSLM(slm).name
        cmd.opt(f"--slm={slm}")

    if fwhm:
        cmd.opt(f"--fwhm={fwhm}")

    if niter:
        cmd.opt(f"--niter={niter}")

    if cnr_maps:
        cmd.opt("--cnr_maps")

    if residuals:
        cmd.opt("--residuals")

    if fep:
        cmd.opt("--fep")

    if interp:
        interp: str = ECInterp(interp).name
        cmd.opt(f"--interp={interp}")

    if resamp:
        resamp: str = ECresamp(resamp).name
        cmd.opt(f"--resamp={resamp}")

    if nvoxhp:
        cmd.opt(f"--nvoxhp={nvoxhp}")

    if initrand:
        cmd.opt(f"--initrand={initrand}")

    if ff:
        cmd.opt(f"--ff={ff}")

    if repol:
        cmd.opt("--repol")

    if ol_nstd:
        cmd.opt(f"--ol_nstd={ol_nstd}")

    if ol_nvox:
        cmd.opt(f"--ol_nvox={ol_nvox}")

    if ol_type:
        ol_type: str = ECOLType(ol_type).name
        cmd.opt(f"--ol_type={ol_type}")

    if ol_pos:
        cmd.opt("--ol_pos")

    if ol_sqr:
        cmd.opt("--ol_sqr")

    if estimate_move_by_susceptibility:
        cmd.opt("--estimate_move_by_susceptibility")

    if mbs_niter:
        cmd.opt(f"--mbs_niter={mbs_niter}")

    if mbs_lambda:
        cmd.opt(f"--mbs_lambda={mbs_lambda}")

    if mbs_ksp:
        cmd.opt(f"--mbs_ksp={mbs_ksp}")

    if dont_sep_offs_move:
        cmd.opt("--dont_sep_offs_move")

    if dont_peas:
        cmd.opt("--dont_peas")

    if data_is_shelled:
        cmd.opt("--data_is_shelled")

    if b0_only:
        cmd.opt("--b0_only")

    if dont_mask_output:
        cmd.opt("--dont_mask_output")

    if verbose:
        cmd.opt("--verbose")

    if very_verbose:
        cmd.opt("--very_verbose")

    # Eddy slice-to-volume (s2v) correction options
    if mporder:
        cmd.opt(f"--mporder={mporder}")

        if slspec:
            slspec: File = File(src=slspec, assert_exists=True)
            cmd.opt(f"--slspec={slspec.src}")
        elif json_file:
            json_file: File = File(src=json_file, assert_exists=True)
            cmd.opt(f"--json={json_file.src}")
        else:
            raise FSLError(
                "Either the 'slspec' or the 'json' option must be specified with the 'mporder' option."
            )

        if s2v_lambda:
            cmd.opt(f"--s2v_lambda={s2v_lambda}")

        if s2v_fwhm:
            cmd.opt(f"--s2v_fwhm={s2v_fwhm}")

        if s2v_niter:
            cmd.opt(f"--s2v_niter={s2v_niter}")

        if s2v_interp:
            s2v_interp: str = ECInterp(s2v_interp).name
            cmd.opt(f"--s2v_interp={s2v_interp}")

    cmd.run(log=log)

    with NiiFile(src=out.src, assert_exists=True, validate_nifti=True) as of:
        outdir, fname, ext = of.file_parts()
        eddy_mask: str = os.path.join(outdir, f"{fname}.eddy_output_mask{ext}")
        eddy_motion_par: str = os.path.join(outdir, f"{fname}.eddy_parameters")
        with NiiFile(src=eddy_mask, assert_exists=True, validate_nifti=True) as em:
            with File(src=eddy_motion_par, assert_exists=True) as emp:
                out: str = of.abspath()
                eddy_mask: str = em.abspath()
                eddy_motion_par: str = emp.abspath()

    return (out, eddy_motion_par, eddy_mask)


def bet(
    img: str,
    out: str,
    mask: bool = False,
    robust: bool = False,
    seg: bool = True,
    frac_int: Optional[float] = None,
    verbose: bool = False,
    log: Optional[LogFile] = None,
) -> Tuple[str, str]:
    """Performs brain extraction."""
    img: NiiFile = NiiFile(src=img, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(src=out, assert_exists=False, validate_nifti=False)

    cmd: Command = Command("bet")
    cmd.opt(img.abspath())
    cmd.opt(out.abspath())

    if mask:
        cmd.opt("-m")
        mask_img: str = out.rm_ext() + "_mask.nii.gz"
    else:
        mask_img: str = ""

    if robust:
        cmd.opt("-R")

    if frac_int:
        cmd.opt("-f")
        cmd.opt(f"{frac_int}")

    if seg:
        pass
    else:
        cmd.opt("-n")

    if verbose:
        cmd.opt("-v")

    cmd.run(log=log)

    return (out.src, mask_img)


def topup(
    img: str,
    param: str,
    out: str,
    config: Optional[str] = None,
    fout: bool = False,
    iout: bool = False,
    subsamp: Optional[str] = None,
    scale: int = 1,
    verbose: bool = False,
    log: Optional[LogFile] = None,
) -> Tuple[str, str, str]:
    """Performs field map estimation from reversed phase encoded fieldmaps or single-band references (e.g. b0s)."""
    img: NiiFile = NiiFile(src=img, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(src=out, assert_exists=False, validate_nifti=False)
    param: File = File(src=param, assert_exists=True)

    cmd: Command = Command("topup")
    cmd.opt(f"--imain={img.abspath()}")
    cmd.opt(f"--out={out.rm_ext()}")
    cmd.opt(f"--datain={param.abspath()}")
    cmd.opt(f"--scale={scale}")

    if config:
        config: File = File(src=config, assert_exists=True)
        cmd.opt(f"--config={config.abspath()}")

    if fout:
        field_img: str = out.rm_ext() + "_field.nii.gz"
        field_img: File = File(src=field_img)
        cmd.opt(f"--fout={field_img.rm_ext()}")
    else:
        field_img: str = ""
        field_img: File = File(src=field_img)

    if iout:
        unwarped_img: str = out.rm_ext() + "_unwarped.nii.gz"
        unwarped_img: NiiFile = NiiFile(src=unwarped_img)
        cmd.opt(f"--iout={unwarped_img.rm_ext()}")
    else:
        unwarped_img: str = ""
        unwarped_img: File = File(src=unwarped_img)

    if subsamp:
        cmd.opt(f"--subsamp={subsamp}")

    if verbose:
        cmd.opt("--verbose")

    cmd.run(log=log)

    return (out.src, field_img.src, unwarped_img.src)


def fslreorient2std(
    img: str,
    out: Optional[str] = None,
    out_mat: bool = False,
    log: Optional[LogFile] = None,
) -> Tuple[str, str, str]:
    """Reorients input image to MNI 152 orientation."""
    img: NiiFile = NiiFile(src=img, assert_exists=True, validate_nifti=True)
    cmd: Command = Command("fslreorient2std")

    if out:
        out: NiiFile = NiiFile(src=out)
    else:
        out: str = ""
        out: NiiFile = NiiFile(src=out)

    if out_mat:
        out_mat: str = out.rm_ext() + "_native2std_cmd.mat"
        out_mat: File = File(src=out_mat)
        cmd.opt("-m")
        cmd.opt(f"{out_mat.src}")
    else:
        out_mat: str = ""
        out_mat: File = File(src=out_mat)

    cmd.opt(f"{img.src}")

    if out:
        cmd.opt(f"{out.src}")

    cmd.run(log=log)

    return (img.src, out.src, out_mat.src)


def fslroi(
    img: str,
    out: str,
    xmin: Optional[int] = None,
    xsize: Optional[int] = None,
    ymin: Optional[int] = None,
    ysize: Optional[int] = None,
    zmin: Optional[int] = None,
    zsize: Optional[int] = None,
    tmin: Optional[int] = None,
    tsize: Optional[int] = None,
    log: Optional[LogFile] = None,
) -> str:
    """work"""
    img: NiiFile = NiiFile(src=img, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(src=out, assert_exists=False, validate_nifti=False)

    if (xmin and xsize and ymin and ysize and zmin and zsize) and (tmin and tsize):
        raise FSLError("Cannot specify both xyz and temporal dimensions.")
    elif (
        xmin is None
        and xsize is None
        and ymin is None
        and ysize is None
        and zmin is None
        and zsize is None
    ) and (tmin is None and tsize is None):
        raise FSLError("Neither xyz nor temporal dimensions were specified.")
    else:
        cmd: Command = Command("fslroi")
        cmd.opt(img.src)
        cmd.opt(out.src)

    if (
        xmin is not None
        and xsize is not None
        and ymin is not None
        and ysize is not None
        and zmin is not None
        and zsize is not None
    ):
        cmd.opt(f"{xmin}")
        cmd.opt(f"{xsize}")
        cmd.opt(f"{ymin}")
        cmd.opt(f"{ysize}")
        cmd.opt(f"{zmin}")
        cmd.opt(f"{zsize}")
    elif (
        (xmin is None)
        and (xsize is None)
        and (ymin is None)
        and (ysize is None)
        and (zmin is None)
        and (zsize is None)
    ):
        pass
    else:
        raise FSLError("Either the xyz min or size was not specified.")

    if tmin is not None and tsize is not None:
        cmd.opt(f"{tmin}")
        cmd.opt(f"{tsize}")
    elif (tmin is None) and (tsize is None):
        pass
    else:
        raise FSLError("Either the temporal min or size was not specified.")

    cmd.run(log=log)

    return out.src


def sigloss(
    img: str,
    out: str,
    te: Optional[float] = None,
    slicedir: Optional[str] = None,
    mask: Optional[str] = None,
    log: Optional[LogFile] = None,
) -> str:
    """Estimate signal loss from a field map (in rad/s).
    """
    with NiiFile(src=img, assert_exists=True, validate_nifti=True) as im:
        img: str = im.abspath()

    cmd: Command = Command("sigloss")
    cmd.opt("-i")
    cmd.opt(f"{img}")
    cmd.opt("-s")
    cmd.opt(f"{out}")

    if te:
        te: float = float(te)
        cmd.opt(f"--te={te}")
    if slicedir:
        cmd.opt(f"--slicedir={slicedir}")
    if mask:
        cmd.opt(f"--mask={mask}")

    cmd.run(log=log)
    return out


def fslmerge(
    out: str,
    merge_opt: str = "t",
    tr: Optional[float] = None,
    log: Optional[LogFile] = None,
    *args,
) -> str:
    """Merges a series of 3D NIFTI files."""
    out: NiiFile = NiiFile(src=out)
    merge_opt: str = MergeDim(merge_opt).name

    cmd: Command = Command("fslmerge")
    cmd.opt(f"-{merge_opt}")
    cmd.opt(out.src)

    for arg in args:
        img: NiiFile = NiiFile(src=arg, assert_exists=True, validate_nifti=True)
        cmd.opt(img.abspath())

    if tr:
        cmd.opt(f"{tr}")

    cmd.run(log=log)

    return out.src


def catmats(matdir: str, out: str) -> str:
    """Concatenate ``FSL`` linear trasformations files into a single file."""
    mats: List[str] = glob.glob(matdir, "MAT_????")
    with open(out, "w") as output_file:
        for mat_file in mats:
            with open(mat_file) as f:
                output_file.write(f.read())
    return out


def applywarp(
    src: str,
    ref: str,
    out: str,
    warp: Optional[str] = "",
    premat: Optional[str] = "",
    prematdir: Optional[str] = "",
    postmat: Optional[str] = "",
    postmatdir: Optional[str] = "",
    interp: str = "spline",
    paddingsize: Optional[int] = None,
    abs: bool = False,
    rel: bool = False,
    log: Optional[LogFile] = None,
) -> str:
    """Applies ``FSL`` warps."""
    assert (
        warp or premat or postmat or prematdir or postmatdir
    ), "either a warp or mat (premat, postmat or prematdir) must be supplied"
    assert not (
        premat and prematdir
    ), "cannot use premat and prematdir arguments together"
    assert not (
        postmat and postmatdir
    ), "cannot use postmat and postmatdir arguments together"

    src: NiiFile = NiiFile(src=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(src=ref, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(src=out, assert_exists=False, validate_nifti=False)

    cmd: Command = Command("applywarp")

    cmd.opt(f"--in={src.abspath()}")
    cmd.opt(f"--ref={ref.abspath()}")
    cmd.opt(f"--out={out.src}")

    if interp:
        interp: str = RegInterp(interp).name
        cmd.opt(f"--interp={interp}")

    if prematdir:
        premat: str = os.path.join(prematdir, "allmats.txt")
        premat: str = catmats(matdir=prematdir, out=premat)

    if postmatdir:
        postmat: str = os.path.join(postmatdir, "allmats.txt")
        postmat: str = catmats(matdir=postmatdir, out=postmat)

    if warp:
        warp: NiiFile = NiiFile(src=warp, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--warp={warp.src}")

    if premat:
        cmd.opt(f"--premat={premat}")

    if postmat:
        cmd.opt(f"--postmat={postmat}")

    if paddingsize and (isinstance(paddingsize, int)):
        cmd.opt(f"--paddingsize={paddingsize}")

    if abs:
        cmd.opt("--abs")

    if rel:
        cmd.opt("--rel")

    cmd.run(log=log)
    return out.src


def invxfm(inmat: str, outmat: str, log: Optional[LogFile] = None) -> str:
    """Inverts ``FSL`` transformation matrices."""
    inmat: File = File(src=inmat, assert_exists=True)
    outmat: File = File(src=outmat, assert_exists=False)

    cmd: Command = Command("convert_xfm")

    cmd.opt("-omat")
    cmd.opt(outmat.src)

    cmd.opt("-inverse")
    cmd.opt(inmat.src)

    cmd.run(log=log)
    return outmat.src


def applyxfm(
    src: str,
    ref: str,
    mat: str,
    out: str,
    interp: Optional[str] = "trilinear",
    log: Optional[LogFile] = None,
) -> str:
    """Applies ``FSL`` transformation matrices."""
    src: NiiFile = NiiFile(src=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(src=ref, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(src=out, assert_exists=False, validate_nifti=False)
    mat: File = File(src=mat, assert_exists=True)

    cmd: Command = Command("flirt")

    cmd.opt("-init")
    cmd.opt(mat.src)

    cmd.opt("-in")
    cmd.opt(src.src)

    cmd.opt("-ref")
    cmd.opt(ref.src)

    cmd.opt("-applyxfm")

    cmd.opt("-out")
    cmd.opt(out.src)

    interp: str = RegInterp(interp).name
    cmd.opt("-interp")
    cmd.opt(interp)

    cmd.run(log=log)
    return out.src


def apply_isoxfm(
    src: str,
    ref: str,
    res: int,
    out: str,
    interp: Optional[str] = "interp",
    log: Optional[str] = None,
) -> str:
    """Resamples images to an isometric resolution."""
    src: NiiFile = NiiFile(src=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(src=ref, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(src=out, assert_exists=False, validate_nifti=False)

    fsldir: str = os.getenv("FSLDIR", None)
    assert fsldir is not None, "FSLDIR environment must be set."

    ident: str = os.path.join(fsldir, "etc/flirtsch/ident.mat")

    cmd: Command = Command("flirt")

    cmd.opt("-init")
    cmd.opt(ident)

    cmd.opt("-in")
    cmd.opt(src.src)

    cmd.opt("-ref")
    cmd.opt(ref.src)

    cmd.opt("-applyisoxfm")
    cmd.opt(f"{res}")

    cmd.opt("-out")
    cmd.opt(out.src)

    interp: str = RegInterp(interp).name
    cmd.opt("-interp")
    cmd.opt(interp)

    cmd.run(log=log)
    return out.src


def concatxfm(
    inmat1: str, inmat2: str, outmat: str, log: Optional[LogFile] = None
) -> str:
    """Concatenates two ``FSL`` transformation matrices."""
    inmat1: File = File(src=inmat1, assert_exists=True)
    inmat2: File = File(src=inmat2, assert_exists=True)

    cmd: Command = Command("convert_xfm")

    cmd.opt("-omat")
    cmd.opt(outmat)

    cmd.opt("-concat")
    cmd.opt(inmat1.src)
    cmd.opt(inmat2.src)

    cmd.run(log=log)
    return outmat


def invwarp(
    inwarp: str,
    ref: str,
    outwarp: str,
    rel: bool = False,
    abs: bool = False,
    noconstraint: bool = False,
    jmin: float = 0.01,
    jmax: float = 100.0,
    verbose: bool = False,
    log: Optional[LogFile] = None,
) -> str:
    """Invert existing warps."""
    inwarp: NiiFile = NiiFile(src=inwarp, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(src=ref, assert_exists=True, validate_nifti=True)
    outwarp: NiiFile = NiiFile(src=outwarp, assert_exists=False, validate_nifti=False)

    cmd: Command = Command("invwarp")

    # Required arguments
    cmd.opt(f"--warp={inwarp.abspath()}")
    cmd.opt(f"--ref={ref.abspath()}")
    cmd.opt(f"--out={outwarp.src}")

    # Optional arguments
    if rel:
        cmd.opt("--rel")

    if abs:
        cmd.opt("--abs")

    if noconstraint:
        cmd.opt("--noconstraint")

    if jmin:
        cmd.opt(f"--jmin={jmin}")

    if jmax:
        cmd.opt(f"--jmax={jmax}")

    if verbose:
        cmd.opt("--verbose")

    cmd.run(log=log)
    return outwarp.src


def convertwarp(
    out: str,
    ref: str,
    warp1: Optional[str] = None,
    warp2: Optional[str] = None,
    premat: Optional[str] = None,
    midmat: Optional[str] = None,
    postmat: Optional[str] = None,
    shiftmap: Optional[str] = None,
    shiftdir: Optional[str] = None,
    abs: bool = False,
    absout: bool = False,
    rel: bool = False,
    relout: bool = False,
    log: Optional[LogFile] = None,
) -> str:
    """Convert ``FSL`` non-linear warps."""
    assert warp1 or warp2 or premat or midmat or postmat, (
        "either a warp (warp1 or warp2) or mat (premat, midmat, or "
        + "postmat) must be supplied"
    )

    ref: NiiFile = NiiFile(src=ref, assert_exists=True, validate_nifti=True)
    out: NiiFile = NiiFile(src=out)

    cmd: Command = Command("convertwarp")

    cmd.opt(f"--ref={ref.src}")
    cmd.opt(f"--out={out.src}")

    if warp1:
        warp1: NiiFile = NiiFile(src=warp1, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--warp1={warp1.src}")

    if warp2:
        warp2: NiiFile = NiiFile(src=warp2, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--warp2={warp2.src}")

    if premat:
        premat: File = File(src=premat, assert_exists=True)
        cmd.opt(f"--premat={premat.src}")

    if midmat:
        midmat: File = File(src=midmat, assert_exists=True)
        cmd.opt(f"--midmat={midmat.src}")

    if postmat:
        postmat: File = File(src=postmat, assert_exists=True)
        cmd.opt(f"--postmat={postmat.src}")

    if shiftmap:
        shiftmap: NiiFile = NiiFile(
            src=shiftmap, assert_exists=True, validate_nifti=True
        )
        cmd.opt(f"--shiftmap={shiftmap.src}")

    if shiftdir is None:
        pass
    elif (
        shiftdir == "x"
        or shiftdir == "x-"
        or shiftdir == "y"
        or shiftdir == "y-"
        or shiftdir == "z"
        or shiftdir == "z-"
    ):
        cmd.opt(f"--shiftdir={shiftdir}")
    else:
        raise FSLError(
            f"Invalid input: {shiftdir}. Valid inputs for 'shiftdir' option includes: x,y,z,x-,y-,z-"
        )

    if absout:
        cmd.opt(f"--absout")

    if relout:
        cmd.opt(f"--relout")

    if abs:
        cmd.opt(f"--abs")

    if rel:
        cmd.opt(f"--rel")

    cmd.run(log=log)

    return out.src


def fugue(
    unmaskshift: bool = False,
    despike: bool = False,
    unmaskfmap: bool = False,
    nocheck: bool = False,
    log: Optional[LogFile] = None,
    **kwargs,
) -> Tuple[str, str, str, str]:
    """FMRIB's Utility for Geometric Unwarping of EPIs."""
    cmd: Command = Command("fugue")

    if unmaskshift:
        cmd.opt("--unmaskshift")
    if despike:
        cmd.opt("--despike")
    if unmaskfmap:
        cmd.opt("--unmaskfmap")
    if nocheck:
        cmd.opt("--nocheck")

    unwarp: str = kwargs.get("unwarp", "")
    warp: str = kwargs.get("warp", "")
    savefmap: str = kwargs.get("savefmap", "")
    saveshift: str = kwargs.get("saveshift", "")

    for k, v in kwargs.items():
        cmd.opt(f"--{k}={v}")

    cmd.run(log=log)

    return (unwarp, warp, savefmap, saveshift)


def flirt(
    src: str,
    ref: str,
    out: Optional[str] = None,
    omat: Optional[str] = None,
    dof: Optional[int] = None,
    cost: Optional[str] = None,
    wmseg: Optional[str] = None,
    init: Optional[str] = None,
    schedule: Optional[str] = None,
    echospacing: Optional[float] = None,
    pedir: Optional[str] = None,
    fieldmap: Optional[str] = None,
    fieldmapmask: Optional[str] = None,
    bbrslope: Optional[float] = None,
    bbrtype: Optional[str] = None,
    interp: Optional[str] = None,
    refweight: Optional[str] = None,
    applyisoxfm: Optional[float] = None,
    usesqform: bool = False,
    nosearch: bool = False,
    verbose: Optional[int] = 0,
    searchrx: Tuple[int, int] = None,
    searchry: Tuple[int, int] = None,
    searchrz: Tuple[int, int] = None,
    log: Optional[LogFile] = None,
) -> Tuple[str]:
    """work"""
    src: NiiFile = NiiFile(src=src, assert_exists=True, validate_nifti=True)
    ref: NiiFile = NiiFile(src=ref, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("flirt")

    cmd.opt("-in")
    cmd.opt(f"{src.src}")
    cmd.opt("-ref")
    cmd.opt(f"{ref.src}")

    if out:
        cmd.opt("-out")
        cmd.opt(f"{out}")

    if omat:
        cmd.opt("-omat")
        cmd.opt(f"{omat}")

    if dof:
        cmd.opt("-dof")
        cmd.opt(f"{dof}")

    if cost:
        # TODO: set enum here
        cmd.opt("-cost")
        cmd.opt(f"{cost}")

    if wmseg:
        wmseg: NiiFile = NiiFile(src=wmseg, assert_exists=True, validate_nifti=True)
        cmd.opt("-wmseg")
        cmd.opt(f"{wmseg.src}")

    if init:
        init: File = File(src=init, assert_exists=True)
        cmd.opt("-init")
        cmd.opt(f"{init.src}")

    if schedule:
        schedule: File = File(src=schedule, assert_exists=True)
        cmd.opt("-schedule")
        cmd.opt(f"{schedule.src}")

    if echospacing:
        cmd.opt("-echospacing")
        cmd.opt(f"{echospacing}")

    if pedir:
        cmd.opt("-pedir")
        cmd.opt(f"{pedir}")

    if fieldmap:
        fieldmap: NiiFile = NiiFile(
            src=fieldmap, assert_exists=True, validate_nifti=True
        )
        cmd.opt("-fieldmap")
        cmd.opt(f"{fieldmap.src}")

    if fieldmapmask:
        fieldmapmask: NiiFile = NiiFile(
            src=fieldmapmask, assert_exists=True, validate_nifti=True
        )
        cmd.opt("-fieldmapmask")
        cmd.opt(f"{fieldmapmask.src}")

    if bbrslope:
        cmd.opt("-bbrslope")
        cmd.opt(f"{bbrslope}")

    if bbrtype:
        # TODO: set enum here
        cmd.opt("-bbrtype")
        cmd.opt(f"{bbrtype}")

    if interp:
        # TODO: set enum here
        cmd.opt("-interp")
        cmd.opt(f"{interp}")

    if refweight:
        refweight: NiiFile = NiiFile(
            src=refweight, assert_exists=True, validate_nifti=True
        )
        cmd.opt("-refweight")
        cmd.opt(f"{refweight.src}")

    if applyisoxfm:
        cmd.opt("-applyisoxfm")
        cmd.opt(f"{applyisoxfm}")

    if usesqform:
        cmd.opt("-usesqform")

    if nosearch:
        cmd.opt("-nosearch")

    if verbose:
        cmd.opt("-verbose")
        cmd.opt(f"{verbose}")

    if searchrx:
        cmd.opt("-searchrx")
        cmd.opt(f"{searchrx[0]}")
        cmd.opt(f"{searchrx[1]}")

    if searchry:
        cmd.opt("-searchry")
        cmd.opt(f"{searchry[0]}")
        cmd.opt(f"{searchry[1]}")

    if searchrz:
        cmd.opt("-searchrz")
        cmd.opt(f"{searchrz[0]}")
        cmd.opt(f"{searchrz[1]}")

    cmd.run(log=log)

    # TODO: set all possible output types here,
    #   set, their defaults above.

    return (out, omat)


def melodic(
    input: str,
    outdir: str,
    dim: Optional[int] = 0,
    tr: Optional[float] = None,
    mmthresh: Optional[float] = None,
    report: bool = False,
    nomask: bool = False,
    updatemask: bool = False,
    nobet: bool = False,
    mask: Optional[str] = "",
    Oall: bool = False,
    verbose: bool = False,
    log: Optional[LogFile] = None,
) -> str:
    """Multivariate Exploratory Linear Optimised ICA."""
    input: NiiFile = NiiFile(src=input, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("melodic")

    cmd.opt(f"--in={input.src}")
    cmd.opt(f"--outdir={outdir}")

    if mask:
        mask: NiiFile = NiiFile(src=mask, assert_exists=True, validate_nifti=True)
        cmd.opt(f"--mask={mask.src}")

    if mmthresh:
        cmd.opt(f"--mmthresh={mmthresh}")

    if dim:
        cmd.opt(f"--dim={dim}")

    if tr:
        cmd.opt(f"--tr={tr}")

    if report:
        cmd.opt("--report")

    if Oall:
        cmd.opt("--Oall")

    if nomask:
        cmd.opt("--nomask")

    if updatemask:
        cmd.opt("--update_mask")

    if nobet:
        cmd.opt("--nobet")

    if verbose:
        cmd.opt("--verbose")

    cmd.run(log=log)

    # TODO:
    #   * Get the corresponding output files for
    #       each melodic options.

    return outdir


def fsl_regfilt(
    infile: str, outfile: str, mix: str, ics: List[int], log: Optional[LogFile] = None
) -> str:
    """Data de-noising by regression.
    """
    infile: NiiFile = NiiFile(src=infile, assert_exists=True, validate_nifti=True)
    outfile: NiiFile = NiiFile(src=outfile, assert_exists=False, validate_nifti=False)
    mix: File = File(src=mix, assert_exists=True)

    # NOTE: Semi-original code from dchp-fmri:
    #   link: https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/blob/master/dhcp/util/fslpy.py
    #
    #   This approach shifted the IC number by 1.
    #
    # icstr: str = '"'
    # for i in range(0, len(ics)-1):
    #     icstr: str = icstr + f"{ics[i]+1},"
    # icstr: str = icstr + f"{ics[-1]+1}"

    icstr: str = '"'
    for i in range(0, len(ics) - 1):
        icstr: str = icstr + f"{ics[i]},"
    icstr: str = icstr + f'{ics[-1]}"'

    cmd: Command = Command("fsl_regfilt")

    cmd.opt(f"--in={infile.src}")
    cmd.opt(f"--out={outfile.src}")
    cmd.opt(f"--design={mix.src}")
    cmd.opt(f"--filter={icstr}")

    cmd.run(log=log)

    # TODO:
    #   * Get the corresponding output files for
    #       each fsl_regfilt options.

    return outfile.src


def mcflirt(
    infile: str,
    outfile: str,
    reffile: Optional[str] = None,
    spline_final: bool = True,
    plots: bool = True,
    mats: bool = True,
    refvol: Optional[int] = None,
    log: Optional[LogFile] = None,
) -> Tuple[str, str, str]:
    """Rigid-body motion correction using ``mcflirt``."""
    infile: NiiFile = NiiFile(src=infile, assert_exists=True, validate_nifti=True)
    outfile: File = File(src=outfile)
    matsdir: str = ""
    parfile: str = outfile.rm_ext() + ".par"

    cmd: Command = Command("mcflirt")

    cmd.opt("-in")
    cmd.opt(f"{infile.src}")
    cmd.opt("-out")
    cmd.opt(f"{outfile.rm_ext()}")

    if reffile:
        reffile: File = File(src=reffile, assert_exists=True)
        cmd.opt("-reffile")
        cmd.opt(f"{reffile.src}")

    if spline_final:
        cmd.opt("-spline_final")

    if plots:
        cmd.opt("-plots")

    if mats:
        matsdir: str = outfile.rm_ext() + ".mat"
        cmd.opt("-mats")

    if refvol:
        assert isinstance(
            refvol, int
        ), f"'refvol' option requires integers, refvol input: {refvol}"
        cmd.opt("-refvol")
        cmd.opt(f"{refvol}")

    cmd.run(log=log)

    return (outfile.src, parfile, matsdir)


def slicer(
    input: str,
    input2: Optional[str] = None,
    label: Optional[Union[int, str]] = None,
    lut: Optional[str] = None,
    intensity: Optional[Tuple[int]] = None,
    edgethreshold: Optional[int] = None,
    x: Optional[Tuple[int, str]] = None,
    y: Optional[Tuple[int, str]] = None,
    z: Optional[Tuple[int, str]] = None,
    log: Optional[LogFile] = None,
) -> str:
    """Creates a combined NIFTI image using one or two NIFTI files.
    """
    # TODO: set output file(s) for this function
    input: NiiFile = NiiFile(src=input, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("slicer")
    cmd.opt(f"{input.src}")

    if input2:
        input2: NiiFile = NiiFile(src=input2, assert_exists=True, validate_nifti=True)
        cmd.opt(f"{input2.src}")

    if label:
        cmd.opt("-L")
        cmd.opt(f"{label}")

    if lut:
        lut: File = File(src=lut, assert_exists=True)
        cmd.opt("-l")
        cmd.opt(f"{lut.src}")

    if intensity:
        cmd.opt("-i")
        cmd.opt(f"{intensity[0]}")
        cmd.opt(f"{intensity[1]}")

    if edgethreshold:
        cmd.opt("-e")
        cmd.opt(f"{edgethreshold}")

    if x:
        cmd.opt("-x")
        cmd.opt(f"{x[0]}")
        cmd.opt(f"{x[1]}")

    if y:
        cmd.opt("-y")
        cmd.opt(f"{y[0]}")
        cmd.opt(f"{y[1]}")

    if z:
        cmd.opt("-z")
        cmd.opt(f"{z[0]}")
        cmd.opt(f"{z[1]}")

    cmd.run(log=log)

    return None


def cluster(
    infile: str,
    thresh: Optional[int] = None,
    oindex: Optional[str] = None,
    no_table: bool = False,
    log: Optional[LogFile] = None,
):
    """Form clusters, report information about clusters and/or perform cluster-based inference.
    """
    infile: NiiFile = NiiFile(src=infile, assert_exists=True, validate_nifti=True)

    cmd: Command = Command("cluster")
    cmd.opt(f"--in={infile.src}")

    if thresh:
        cmd.opt(f"--thresh={thresh}")

    if oindex:
        cmd.opt(f"--oindex={oindex}")

    if no_table:
        cmd.opt("--no_table")

    cmd.run(log=log)

    # TODO:
    #   * Get the corresponding output files for
    #       each cluster options.
    #   * Get number and types of output files.

    return None


class fslmaths:
    """``FSL`` wrapper class for the ``fslmaths`` utility executable.
    Perform mathematical operations and/or manipulation of images.
    """

    def __init__(
        self, img: str, dt: Optional[str] = None,
    ):
        """Constructor"""
        img: NiiFile = NiiFile(src=img, assert_exists=True, validate_nifti=True)
        self._cmds: Command = Command("fslmaths")

        if dt:
            dt: str = FSLDataType(dt).name
            self._cmds.opt("-dt")
            self._cmds.opt(dt)

        self._cmds.opt(img.src)

    def abs(self):
        """work"""
        self._cmds.opt("-abs")
        return self

    def bin(self):
        """work"""
        self._cmds.opt("-bin")
        return self

    def binv(self):
        """work"""
        self._cmds.opt("-binv")
        return self

    def recip(self):
        """work"""
        self._cmds.opt("-recip")
        return self

    def Tmean(self):
        """work"""
        self._cmds.opt("-Tmean")
        return self

    def Tstd(self):
        """work"""
        self._cmds.opt("-Tstd")
        return self

    def Tmin(self):
        """work"""
        self._cmds.opt("-Tmin")
        return self

    def Tmax(self):
        """work"""
        self._cmds.opt("-Tmax")
        return self

    def sqrt(self):
        """work"""
        self._cmds.opt("-sqrt")
        return self

    def sqr(self):
        """work"""
        self._cmds.opt("-sqr")
        return self

    def fillh(self):
        """work"""
        self._cmds.opt("-fillh")
        return self

    def ero(self, repeat: int = 1):
        """work"""
        for _ in range(repeat):
            self._cmds.opt("-ero")
        return self

    def dilM(self, repeat: int = 1):
        """work"""
        for _ in range(repeat):
            self._cmds.opt("-dilM")
        return self

    def dilF(self, repeat: int = 1):
        """work"""
        for _ in range(repeat):
            self._cmds.opt("-dilF")
        return self

    def add(self, input: Union[int, float, str]):
        """work"""
        self._cmds.opt("-add")

        if isinstance(input, int) or isinstance(input, float):
            self._cmds.opt(f"{input}")
        elif isinstance(input, str):
            img: NiiFile = NiiFile(src=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.src)
        return self

    def sub(self, input: Union[int, float, str]):
        """work"""
        self._cmds.opt("-sub")

        if isinstance(input, int) or isinstance(input, float):
            self._cmds.opt(f"{input}")
        elif isinstance(input, str):
            img: NiiFile = NiiFile(src=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.src)
        return self

    def mul(self, input: Union[int, float, str]):
        """work"""
        self._cmds.opt("-mul")

        if isinstance(input, int) or isinstance(input, float):
            self._cmds.opt(f"{input}")
        elif isinstance(input, str):
            img: NiiFile = NiiFile(src=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.src)
        return self

    def div(self, input: Union[int, float, str]):
        """work"""
        self._cmds.opt("-div")

        if isinstance(input, int) or isinstance(input, float):
            self._cmds.opt(f"{input}")
        elif isinstance(input, str):
            img: NiiFile = NiiFile(src=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.src)
        return self

    def mas(self, img: str):
        """work"""
        img: NiiFile = NiiFile(src=img, assert_exists=True, validate_nifti=True)
        self._cmds.opt("-mas")
        self._cmds.opt(img.src)
        return self

    def rem(self, input: Union[int, float, str]):
        """work"""
        self._cmds.opt("-rem")

        if isinstance(input, int) or isinstance(input, float):
            self._cmds.opt(f"{input}")
        elif isinstance(input, str):
            img: NiiFile = NiiFile(src=input, assert_exists=True, validate_nifti=True)
            self._cmds.opt(img.src)
        return self

    def thr(self, num: Union[int, float]):
        """work"""
        if isinstance(num, int) or isinstance(num, float):
            self._cmds.opt("-thr")
            self._cmds.opt(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self

    def uthr(self, num: Union[int, float]):
        """work"""
        if isinstance(num, int) or isinstance(num, float):
            self._cmds.opt("-uthr")
            self._cmds.opt(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self

    def inm(self, num: Union[int, float]):
        """work"""
        if isinstance(num, int) or isinstance(num, float):
            self._cmds.opt("-inm")
            self._cmds.opt(f"{num}")
        else:
            raise TypeError(f"Input {num} is not a number.")
        return self

    def bptf(
        self,
        high_pass: Union[int, float],
        low_pass: Union[int, float],
        tr: Optional[Union[int, float]] = None,
        input_is_hz: bool = False,
        input_is_sec: bool = False,
    ):
        """
        Input is assumed to be in sigma (in volume).
        """
        if input_is_hz and input_is_sec:
            raise RuntimeError(
                "Both 'input_is_hz' and 'input_is_sec' were specified. ONLY ONE of these options may be specified."
            )
        elif (input_is_hz or input_is_sec) and (not tr):
            raise FSLError(
                "The TR (Repetition Time) is required when either the 'input_is_hz' or 'input_is_sec' options are specified."
            )

        def _compute_hz(sec: float) -> float:
            """Computes frequency cutoff in hertz (Hz).
            """
            try:
                return 1 / sec
            except ZeroDivisionError:
                return -1.0

        def _compute_sigma(
            tr: Union[int, float], hz: Union[int, float], compute_low_pass: bool = False
        ) -> float:
            """
            relevant links: 
                * https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;fc5b33c5.1205
                * https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;f6fd75a6.1709
            """
            if isinstance(float(tr), float) and isinstance(float(hz), float):

                if compute_low_pass:
                    fwhm_kernel: float = float(18)
                else:
                    fwhm_kernel: float = float(2)

                if hz < 0:
                    return -1.0

                try:
                    sigma_vol: float = np.round(1 / (fwhm_kernel * tr * hz))
                except ZeroDivisionError:
                    sigma_vol: float = -1.0
                return sigma_vol

            else:
                raise TypeError(f"Input TR: {tr} or Hz: {hz} is not a number.")

        if input_is_sec:
            high_pass: float = _compute_sigma(
                tr=tr, hz=_compute_hz(sec=high_pass), compute_low_pass=False
            )
            low_pass: float = _compute_sigma(
                tr=tr, hz=_compute_hz(sec=low_pass), compute_low_pass=True
            )
        elif input_is_hz:
            high_pass: float = _compute_sigma(
                tr=tr, hz=high_pass, compute_low_pass=False
            )
            low_pass: float = _compute_sigma(tr=tr, hz=low_pass, compute_low_pass=True)

        if isinstance(float(high_pass), float) and isinstance(float(low_pass), float):
            self._cmds.opt("-bptf")
            self._cmds.opt(f"{high_pass}")
            self._cmds.opt(f"{low_pass}")
        else:
            raise TypeError(
                f"Input high_pass: {high_pass} or low_pass: {low_pass} is not a number."
            )

        return self

    def run(
        self, out: str, odt: Optional[str] = None, log: Optional[LogFile] = None
    ) -> str:
        """work"""
        out: NiiFile = NiiFile(src=out)
        self._cmds.opt(out.src)

        if odt:
            odt: str = FSLDataType(odt).name
            self._cmds.opt("-odt")
            self._cmds.opt(odt)

        self._cmds.run(log=log)
        return out.src
