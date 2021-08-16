# -*- coding: utf-8 -*-
"""Performs motion correction and distortion correction of MR image data.

NOTE:
    * Check this resource pertaining to the odd number of slices: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;28f6c983.1806
"""
import os
import numpy as np
import nibabel as nib
import pandas as pd

from math import pi as PI
from shutil import copy

from typing import (
    List,
    Optional,
    Tuple,
    Union
)

from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.util import rel_sym_link

from fmri_preproc.utils.fslpy import (
    applywarp,
    eddy,
    FSLDIR,
    fslmaths,
    fslroi,
    mcflirt
)

from fmri_preproc.utils.io import (
    File,
    NiiFile
)

from fmri_preproc.utils.enums import (
    MotionMetric,
    SliceAcqOrder
)

def mcdc(func: str,
         outdir: str,
         func_echospacing: Optional[float] = 0.1,
         func_brainmask: Optional[str] = None,
         func_slorder: Optional[str] = None,
         fmap: Optional[str] = None,
         fmap2func_affine: Optional[str] = None,
         mb_factor: Optional[int] = None,
         dc: bool = False,
         s2v: bool = False,
         mbs: bool = False,
         use_mcflirt: bool = False,
         use_gpu: bool = False,
         ref: Union[str,int] = 0,
         log: Optional[LogFile] = None
        ) -> None:
    """Perform motion and distortion correction stage of ``fmri_preproc``.
    """

    # Logic tests
    _has_fmap = fmap is not None
    _has_slorder = (mb_factor is not None) or (func_slorder is not None)
    _has_acqp = func_echospacing is not None

    if dc:
        pass
    else:
        dc: bool = _has_fmap
    
    if mbs:
        pass
    else:
        mbs: bool = (_has_fmap & _has_acqp & dc)

    if s2v:
        pass
    else:
        s2v: bool = _has_slorder
    
    if use_mcflirt:
        pass
    else:
        use_mcflirt: bool = (not dc) & (not s2v) & (not mbs)

    # Check argument combinations
    if (dc & use_mcflirt) | (s2v & use_mcflirt) | (mbs & use_mcflirt):
        raise RuntimeError('Cannot use MCFLIRT with dc, s2v and/or mbs.')

    if dc and not (_has_fmap & _has_acqp):
        raise RuntimeError('fmap, fmap2func_affine, func_pedir, and func_echospacing are required for DC.')

    if mbs and not (_has_fmap & _has_acqp):
        raise RuntimeError('fmap, fmap2func_affine, func_pedir, and func_echospacing are required for mbs.')

    if mbs and not dc:
        raise RuntimeError('Cannot mbs without dc as they are both distortion corrections.')

    if (not use_mcflirt) & (func_brainmask is None):
        raise RuntimeError('func_brainmask is required to use EDDY')
    
    # Define output files
    mcdir: str = os.path.join(outdir,"mc")
    func_mcdc: str = os.path.join(mcdir,"prefiltered_func_data.nii.gz")
    motfile: str = os.path.join(mcdir,"prefiltered_func_data.par")
    func_mot: str = os.path.join(mcdir,"func_mcdc_motion.tsv")

    # Perform MCDC
    if use_mcflirt:
        (func_mcdc, 
         motfile, 
         matsdir) = mcflirt_mc(func=func,
                               func_mc=func_mcdc,
                               ref=ref,
                               log=log)

        mcf: np.array = np.loadtxt(motfile)
        mcf: pd.DataFrame = pd.DataFrame(mcf, columns=['RotX', 'RotY', 'RotZ', 'X', 'Y', 'Z'])
        mcf.to_csv(func_mot, sep='\t', index=None)
    else:
        (func_mcdc,
         motfile,
         eddy_mask) = eddy_mcdc(func=func,
                                func_brainmask=func_brainmask,
                                func_mcdc=func_mcdc,
                                func_sliceorder=func_slorder,
                                func_echospacing=func_echospacing,
                                fmap=fmap,
                                fmap2func_xfm=fmap2func_affine,
                                mb_factor=mb_factor,
                                mot_params=motfile,
                                mbs=mbs,
                                s2v_corr=s2v,
                                use_gpu=use_gpu,
                                log=log)

        mcf: np.array = np.loadtxt(motfile)
        mcf: pd.DataFrame = pd.DataFrame(mcf, columns=['RotX', 'RotY', 'RotZ', 'X', 'Y', 'Z'])
        mcf.to_csv(func_mot, sep='\t', index=None)
    
    # Calculate post-mc motion outliers
    # TODO: Pick up from here: https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/blob/master/dhcp/func/mcdc.py
    pass

def mcflirt_mc(func: str,
               func_mc: str,
               ref: Optional[Union[int, str]] = None,
               dc_warp: Optional[str] = None,
               log: Optional[LogFile] = None
              ) -> Tuple[str,str,str]:
    """Performs MCFLIRT-based motion and distortion correction.
    """
    func: NiiFile = NiiFile(file=func, assert_exists=True, validate_nifti=True)
    
    if log:
        log.log("Performing motion correction using MCFLIRT.")

    with File(file=func_mc, assert_exists=False) as f:
        outdir, _, _ = f.file_parts()
        with WorkDir(work_dir=outdir) as d:
            d.mkdir()

    func_mc: NiiFile = NiiFile(file=func_mc, assert_exists=False)

    if isinstance(ref, str):
        with NiiFile(file=ref, assert_exists=True, validate_nifti=True) as f:
            ref: str = f.abs_path()
    
    if dc_warp:
        with NiiFile(file=dc_warp, assert_exists=True, validate_nifti=True) as f:
            dc_warp: str = f.abs_path()
    
    ref_vol: Union[int, None] = ref if isinstance(ref, int) else None
    ref_file: Union[str, None] = ref if isinstance(ref, str) else None

    (func_mc,
     parfile, 
     matsdir) = mcflirt(infile=func.abs_path(),
                        outfile=func_mc.file,
                        reffile=ref_file,
                        refvol=ref_vol,
                        log=log)

    if dc_warp:
        with TmpDir(tmp_dir=outdir) as tmp:
            tmp.mkdir()
            func_ref: str = fslroi(img=func_mc,
                                   out=os.path.join(tmp.abs_path(),"func0_ref.nii.gz"),
                                   tmin=0,
                                   tsize=1)
            func_mc: str = applywarp(src=func_mc,
                                     ref=func_ref,
                                     out=func_mc,
                                     warp=dc_warp,
                                     interp='spline')
    return (func_mc, 
            parfile, 
            matsdir)

def eddy_mcdc(func: str,
              func_brainmask: str,
              func_mcdc: str,
              func_sliceorder: Optional[str] = None,
              func_echospacing: Optional[float] = 0.05,
              fmap: Optional[str] = None,
              fmap2func_xfm: Optional[str] = None,
              mb_factor: Optional[int] = 1,
              mot_params: Optional[str] = None,
              mbs: bool = False,
              s2v_corr: bool = False,
              use_gpu: bool = False,
              log: Optional[LogFile] = None
             ) -> Tuple[str,str,str]:
    """Perform EDDY-based motion and distortion correction.

    NOTE: Input ``fmri`` is **ASSUMED** to be in the PA phase encoding direction.
    """
    if log:
        log.log("Performing EDDY-based motion and distortion correction.")

    with NiiFile(file=func, assert_exists=True, validate_nifti=True) as fn:
        with NiiFile(file=func_brainmask, assert_exists=True, validate_nifti=True) as fb:
            with File(file=func_mcdc) as fmc:
                func: str = fn.abs_path()
                func_brainmask: str = fb.abs_path()

                outdir, _, _ = fmc.file_parts()
                eddy_dir: str = os.path.join(outdir,"eddy")
                eddy_basename: str = os.path.join(eddy_dir,"eddy_corr")

                with WorkDir(work_dir=eddy_dir) as d:
                    if log:
                        log.log("Creating eddy output directory.")
                    d.mkdir()
    
    # Logic tests
    _has_acqp = func_echospacing is not None
    _has_fmap = fmap is not None

    if _has_fmap and not _has_acqp:
        raise RuntimeError('Cannot do DC without func_echospacing.')

    if mbs and not _has_fmap:
        raise RuntimeError('Cannot do MBS without fmap.')

    # if s2v_corr and log:
    if s2v_corr and log and use_gpu:
        log.log(f'Performing slice-to-volume (S2V) motion correction')
    elif log:
        log.log(f'Performing rigid-body motion correction')

    if _has_acqp and _has_fmap:
        if mbs and log:
            log.log(f'Performing motion-by-susceptibility distortion correction')
        elif log:
            log.log(f'Performing static distortion correction')
    elif log:
        log.log(f'NO distortion correction')

    # Setup Eddy files
    num_vols: int = nib.load(filename=func).shape[3]
    slices: int = nib.load(filename=func).header.get('dim','')[3]

    idx: str = write_index(num_frames=num_vols, out_file=eddy_basename + "_fmri_pre-mcdc.idx")
    bvals: str = write_bvals(num_frames=num_vols, out_file=eddy_basename + "_fmri_pre-mcdc.bval")
    bvecs: str = write_bvecs(num_frames=num_vols, out_file=eddy_basename + "_fmri_pre-mcdc.bvec")
    acqp: str = write_func_acq_params(num_frames=num_vols, effective_echo_spacing=func_echospacing, out_prefix=eddy_basename + "_fmri_pre-mcdc.acqp")

    # Set default Eddy parameters
    niter: int = 10
    fwhm: str = "10,10,5,5,0,0,0,0,0,0"
    nvoxhp: int = 1000       # Not used for anything by eddy
    cnr_maps: bool = False   # Turned off because not relevant for B0's
    residuals: bool = False  # Turned off because not relevant for B0's

    # Set default Eddy s2v parameters
    s2v_niter:  Union[int, None] = None
    s2v_fwhm:   Union[int, None] = None
    s2v_lambda: Union[int, None] = None
    s2v_interp: Union[str, None] = None
    mporder:    Union[int, None] = None
    mbs_niter:  Union[int, None] = None
    mbs_lambda: Union[int, None] = None
    mbs_ksp:    Union[int, None] = None

    # if s2v_corr:
    if s2v_corr and use_gpu:
        if func_sliceorder:
            with File(file=func_sliceorder, assert_exists=True) as f:
                func_sliceorder: str = f.abs_path()
        else:
            func_sliceorder: str = write_slice_order(slices=slices, 
                                                    mb_factor=mb_factor, 
                                                    mode="single-shot",
                                                    out_file=eddy_basename + "_fmri_pre-mcdc.slice.order",
                                                    return_mat=False)

        s2v_niter:  int = 10
        s2v_fwhm:   int = 0
        s2v_lambda: int = 1
        s2v_interp: str = "trilinear"

        # Set mporder to 16 or smallest (advised by Jesper Anderson and Sean Fitzgibbon)
        mporder:int = np.loadtxt(func_sliceorder).shape[0] - 1
        if mporder > 16:
            mporder: int = 16
        
        if mbs:
            mbs_niter:  int = 20
            mbs_lambda: int = 5
            mbs_ksp:    int = 5

    # Prepare fieldmap transform
    if fmap2func_xfm:
        with File(file=fmap2func_xfm, assert_exists=True) as f:
            fmap2func_xfm: str = f.abs_path()
    else:
        fmap2func_xfm: str = os.path.join(FSLDIR,'etc', 'flirtsch', 'ident.mat')
    
    # Prepare fieldmap
    field_hz: Union[str, None] = None
    if fmap:
        with NiiFile(file=fmap, assert_exists=True, validate_nifti=True) as f:
            _, fname, _ = f.file_parts()
            field_hz: str = os.path.join(eddy_basename + "_" + fname + "_Hz_pre-mcdc")
            field_hz: str = fslmaths(img=fmap).div(2 * PI).run(out=field_hz, log=log)
        
    # Perform Eddy-based mcdc
    (eddy_corr,
     eddy_motion_par,
     eddy_mask) = eddy(img=func,
                       out=eddy_basename,
                       mask=func_brainmask,
                       use_gpu=use_gpu,
                       acqp=acqp,
                       bvecs=bvecs,
                       bvals=bvals,
                       idx=idx,
                       very_verbose=True,
                       niter=niter,
                       fwhm=fwhm,
                       s2v_niter=s2v_niter,
                       mporder=mporder,
                       nvoxhp=nvoxhp,
                       slspec=func_sliceorder,
                       b0_only=True,
                       field=field_hz,
                       field_mat=fmap2func_xfm,
                       s2v_lambda=s2v_lambda,
                       s2v_fwhm=s2v_fwhm,
                       dont_mask_output=True,
                       s2v_interp=s2v_interp,
                       data_is_shelled=True,
                       estimate_move_by_susceptibility=mbs,
                       mbs_niter=mbs_niter,
                       mbs_lambda=mbs_lambda,
                       mbs_ksp=mbs_ksp,
                       cnr_maps=cnr_maps,
                       residuals=residuals,
                       log=log)

    # NOTE: Eddy output files
    # 
    # eddy_corr.eddy_command_txt                    eddy_corr.eddy_values_of_all_input_parameters
    # eddy_corr.eddy_mbs_first_order_fields.nii.gz  eddy_corr_fmap_Hz_pre-mcdc.nii.gz
    # eddy_corr.eddy_movement_over_time             eddy_corr_fmri_pre-mcdc.acqp_functional.acqp
    # eddy_corr.eddy_movement_rms                   eddy_corr_fmri_pre-mcdc.bval
    # eddy_corr.eddy_output_mask.nii.gz             eddy_corr_fmri_pre-mcdc.bvec
    # eddy_corr.eddy_parameters                     eddy_corr_fmri_pre-mcdc.idx
    # eddy_corr.eddy_restricted_movement_rms        eddy_corr_fmri_pre-mcdc.slice.order
    # eddy_corr.eddy_rotated_bvecs                  eddy_corr.nii.gz

    # Necessary eddy output files
    func_mcdc: str = copy(eddy_corr, func_mcdc)

    if mot_params:
        rel_sym_link(target=eddy_motion_par, linkname=mot_params)
    
    # Re-write TR in output NIFTI file header
    func: nib.load(func)

    nib.Nifti1Image(nib.load(func_mcdc).get_fdata(),
                    header=func.header,
                    affine=func.affine
                    ).to_filename(func_mcdc)

    return (func_mcdc,
            eddy_mask,
            mot_params)

def write_bvals(num_frames: int,
                out_file: str = 'file.bval'
               ) -> str:
    """Creates/generates an arbitrary number b=0 b-values for a fMRI acquisition
    provided the number of dynamics/temporal frames.
    
    Usage example:
        >>> bval_file = write_bvals(num_frames=1200, out_file="fmri.bval")
        
    Arguments:
        num_frames: Number of temporal frames/dynamics.
        out_file: Output file name.
        
    Returns:
        Output file path as string that contains the corresponding b-values.
    
    Raises:
        TypeError: Exception that is raised when the input for ``num_frames`` is not an ``int``.
    """
    if not isinstance(num_frames,int):
        raise TypeError(f"Input for num_frames: {num_frames} is not an integer.")
        
    np.savetxt(out_file,
               np.zeros((1,num_frames),
                        dtype=np.int8,
                        order='C'),
               fmt='%i',
               delimiter=' ',
               encoding='utf-8')
    out_file: str = os.path.abspath(out_file)
    return out_file

def write_bvecs(num_frames: int,
                out_file: str = 'file.bvec'
               ) -> str:
    """Creates/generates an arbitrary number of x,y,z b-vectors for a fMRI acquisition
    provided the number of dynamics/temporal frames.
    
    Usage example:
        >>> bval_file = write_bvecs(num_frames=1200, out_file="fmri.bvec")
        
    Arguments:
        num_frames: Number of temporal frames/dynamics.
        out_file: Output file name.
        
    Returns:
        Output file path as string that contains the corresponding b-vectors.
    
    Raises:
        TypeError: Exception that is raised when the input for ``num_frames`` is not an ``int``.
    """
    if not isinstance(num_frames,int):
        raise TypeError(f"Input for num_frames: {num_frames} is not an integer.")
        
    np.savetxt(out_file,
               np.ones((1,num_frames),dtype=np.int8,order='C') * np.array([[1], [0], [0]]),
               fmt='%i',
               delimiter=' ',
               encoding='utf-8')
    out_file: str = os.path.abspath(out_file)
    return out_file

def write_func_acq_params(num_frames: int, 
                          effective_echo_spacing: Optional[float] = 0.05, 
                          out_prefix: Optional[str] = 'file'
                         ) -> str:
    """Creates acquisition parameters files for fMRI acquisition provided the number of frames and the 
    effective echo spacing.
    
    Usage example:
        >>> func_acqp = write_func_acq_params(num_frames=1200, 
        ...                                   effective_echo_spacing=0.05,
        ...                                   out_prefix="fmri")
        ...
        
    Arguments:
        num_frames: Number of temporal frames/dynamics.
        effective_echo_spacing: Effective echo spacing.
        out_prefix: Output file prefix.
        
    Returns:
        String that corresponds to the acquisition parameter file for eddy current correction.
    
    Raises:
        TypeError: Exception that is raised when either ``num_frames`` is not an ``int`` OR when ``effective_echo_spacing`` is not a ``float``.
    """
    if (not isinstance(effective_echo_spacing, int) and 
        not isinstance(effective_echo_spacing, float) or
        not isinstance(num_frames, int)):
        raise TypeError(f"Input for num_frames: {num_frames} is not an integer OR effective_echo_spacing: {effective_echo_spacing} is not a float.")
    
    out_func: str = out_prefix + "_functional.acqp"
    
    with open(out_func,'w') as f:
        for _ in range(0,num_frames):
            f.write(f"0 1 0 {effective_echo_spacing}\n")
        f.close()
    
    out_func: str = os.path.abspath(out_func)
    
    return out_func

def write_index(num_frames: int,
                out_file: str = 'file.idx'
               ) -> str:
    """Creates index files for use with FSL's ``eddy``.
    
    Usage example:
        >>> func_idx = write_index(num_frames=1200, out_file="fmri.idx")
        
    Arguments:
        num_frames: Number of temporal frames/dynamics.
        out_file: Output file name.
        
    Returns:
        Output file path as string that contains the corresponding frame indices.
    
    Raises:
        TypeError: Exception that is raised when the input for ``num_frames`` is not an ``int``.
    """
    if not isinstance(num_frames, int):
        raise TypeError(f"Input for num_frames: {num_frames} is not an integer.")
    np.savetxt(out_file, np.ones((1, num_frames)).T, fmt="%i")
    out_file: str = os.path.abspath(out_file)
    return out_file

def write_slice_order(slices: int,
                      mb_factor: int = 1,
                      mode: str = 'interleaved',
                      out_file: str = 'file.slice.order',
                      return_mat: Optional[bool] = False
                     ) -> Union[str,np.array]:
    """Generates the slice acquisition order file for use with ``eddy's`` slice-to-volume motion correction method.
    
    The file generated consists of an (N/m) x m matrix | N = number of slices in the acquisition direction (assumed to be
    the z-direction), and m = multi-band factor.
    
    NOTE:
        * Links for details:
            * Documentation: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#A--slspec
            * Implementation: https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/blob/master/dhcp/resources/default_func.slorder
            * Implementation: https://git.fmrib.ox.ac.uk/matteob/dHCP_neo_dMRI_pipeline_release/-/blob/master/slorder.txt
        * This wrapper function was written with help from Gregory Lee, PhD.
    
    Usage example:
        >>> sls_order = write_slice_order(slices=44,
        ...                               mb_factor=6,
        ...                               mode='single-shot',
        ...                               out_file='file.slice.order',
        ...                               return_mat=False)
        ...
        
    Arguments:
        slices: Number of slices in the acquisition direction.
        mb_factor: Multi-band factor.
        mode: Acquisition algorithm/method/scheme used to acquire the data. Valid options include:
            * ``interleaved``: Optimal slice acquisition technique in which non-adjacent slices are acquired (best for diffusion, and structural MR images).
            * ``single-shot``: Slice acquisition in which slices are acquired sequentially with an ascending slice order (best for functional MR images).
            * ``default``: Default acquisition order in which slices are acquired with an ascending slice order.
        out_file: Output file name.
        return_mat: Return a numpy multi-dimensional array (matrix) rather than a file.
        
    Returns:
        Output file path as string that contains the slice acquisition order OR a numpy multi-dimensional array if specified.
    
    Raises:
        TypeError: Exception that is raised in the case that either ``slices`` or ``mb_factor`` is not an ``int``.
    """
    # Check input types
    if not isinstance(slices, int) and not isinstance(mb_factor, int):
        raise TypeError(f"Input option slices: {slices} or mb_factor: {mb_factor} is not an integer.")
    
    # Locations (in the slices) divided by Multi-Band Factor
    locs: int = slices//mb_factor

    mode: str = SliceAcqOrder(mode.lower()).name
    
    if mode == 'interleaved':
        step: int = int(np.round(np.sqrt(locs)))
    elif mode == 'default':
        step: int = 2
    elif mode == 'single_shot':
        step: int = 1
    
    # Iterate through each MB acquisition to get slice ordering
    n: List[int] = []
    
    for s in range(step):
        for k in range(s, locs, step):
            if mb_factor != 1:
                a: List[int] = [ k + locs*j for j in range(mb_factor) ]
                n.append(a)
            else:
                a: int = k
                n.append(a)
    
    if return_mat:
        return np.array(n)
    
    slice_order: np.arrary = np.array(n)
    np.savetxt(out_file,
               slice_order, 
               fmt="%i")
    out_file: str = os.path.abspath(out_file)
    return out_file

def _dvars(img: nib.Nifti1Header) -> np.array:
    """Calculates DVARS.
    """
    dvars: np.array = np.square(np.diff(img, axis=3))
    dvars: np.array = np.nanmean(dvars, axis=(0, 1, 2))
    dvars: np.array = np.sqrt(dvars)
    dvars: np.array = 1000 * (dvars / np.nanmedian(img))
    dvars: np.array = np.concatenate(([0], dvars))
    return dvars


def _refrms(img: nib.Nifti1Header,
            ref: Optional[nib.Nifti1Header] = None
           ) -> np.array:
    """Calculates REFRMS.
    """
    if ref:
        pass
    else:
        ref: nib.Nifti1Header = img[:, :, :, int(np.round(img.shape[3] / 2))]

    rms: nib.Nifti1Header = img
    for i in range(0, rms.shape[3]):
        rms[:, :, :, i] = (rms[:, :, :, i] - ref)

    rms: np.array = rms / np.nanmedian(img)
    rms: np.array = np.square(rms)
    rms: np.array = np.nanmean(rms, axis=(0, 1, 2))
    rms: np.array = np.sqrt(rms)
    rms: np.array = np.abs(np.diff(rms))
    rms: np.array = np.concatenate(([0], rms))
    return rms


def motion_outlier(func: str,
                   metric_name: str,
                   plot_name: str,
                   metric: str = "dvars",
                   ref: Optional[str] = None,
                   thr: Optional[int] = None
                  ) -> str:
    """doc-string
    """
    func: NiiFile = NiiFile(file=func, assert_exists=True, validate_nifti=True)

    img0: nib.Nifti1Header = nib.load(func)
    image_data: np.array = img0.get_data().astype(float)

    # Threshold image
    thr2: float = np.nanpercentile(image_data, 2)
    thr98: float = np.nanpercentile(image_data, 98)
    robust_thr: float = thr2 + 0.1 * (thr98 - thr2)
    image_data[image_data < robust_thr] = np.nan

    # Calculate motion metric
    dvars: np.array = _dvars(img=image_data)
    refrms: np.array = _refrms(img=image_data, ref=ref)

    metric: str = MotionMetric(metric).name

    if metric == 'dvars':
        metric_data: np.array = dvars
    elif metric == 'refrms':
        metric_data: np.array = refrms
    
    # Calculate outlier
    if thr:
        pass
    else:
        q75, q25 = np.percentile(metric_data, [75, 25])
        iqr: Union[float,int] = q75 - q25
        thr: Union[float,int] = q75 + (1.5 * iqr)

    outlier: np.array = metric_data > thr

    # Save metric values to file
    pd.DataFrame(
        np.stack((dvars, refrms, outlier), axis=1),
        columns=['DVARS', 'RefRMS', 'Outlier' + metric.name.upper()],
    ).to_csv(metric_name, sep='\t', index=None)

    # Plot metric
    if plot_name:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        _, fname, _ = func.file_parts()

        hfig: plt = plt.figure()
        ax_hdl: plt = hfig.add_subplot(1, 1, 1)

        outlier[0] = False
        h: plt = ax_hdl.fill_between(
                        range(outlier.shape[0]),
                        outlier * np.amax(metric_data),
                        0,
                        color='red',
                        label="outliers",
                        interpolate=True)
        plt.setp(h, alpha=0.25)

        ax_hdl.plot(metric_data, label=metric.name)
        ax_hdl.set_xlim(1, len(metric_data))
        ax_hdl.set_ylim(0, np.nanmean(metric_data) + (4 * np.nanstd(metric_data)))
        ax_hdl.set_ylabel(metric.name)
        ax_hdl.set_title(fname)
        ax_hdl.set_xlabel("time (volumes)")
        plt.legend()
        plt.savefig(plot_name)

    return outlier, metric_data, thr

# This function should be used elsewhere, perhaps outside of this
#   package.
# 
# def _find_optimal_mb_factor_for_eddy_s2v(slices: int,
#                                          mb_factor: int,
#                                          attempts: Optional[int] = None,
#                                          reported_mb: Optional[int] = None,
#                                         ) -> int:
#     """Helper function that finds the optimal multi-band factor to best
#     create/reconstruct the slice order file.
#     """
#     if attempts:
#         pass
#     else:
#         attempts: int = 1
#         reported_mb: int = mb_factor
#     
#     if attempts == 10:
#         possible_mbs: List[int] = list(range(reported_mb-2,13))
# 
#         for mb in possible_mbs:
#             if (slices % mb) == 0:
#                 return mb
#             else:
#                 return reported_mb
#     
#     if (slices % mb_factor) == 0:
#         return mb_factor
#     elif (slices % mb_factor) >= 5:
#         return _find_optimal_mb_factor_for_eddy_s2v(slices, mb_factor-1, attempts+1, reported_mb)
#     elif (slices % mb_factor) <= 5:
#         return _find_optimal_mb_factor_for_eddy_s2v(slices, mb_factor+1, attempts+1, reported_mb)