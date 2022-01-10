# -*- coding: utf-8 -*-
"""Performs motion correction and distortion correction of MR image data.
"""
import os
import numpy as np
import nibabel as nib
import pandas as pd
from warnings import warn

from typing import (
    Dict,
    List,
    Optional,
    Tuple,
    Union
)

from fmri_preproc.utils.util import timeops
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.tempdir import TmpDir
from fmri_preproc.utils.acqparam import write_func_params
from fmri_preproc.utils.outputs.mcdc import MCDCFiles

from fmri_preproc.utils.fslpy import (
    applywarp,
    bet,
    eddy,
    FSLDIR,
    fslmaths,
    fslroi,
    mcflirt
)

from fmri_preproc.utils.fileio import (
    File,
    NiiFile
)

from fmri_preproc.utils.enums import (
    MotionMetric,
    PhaseEncodeDirection,
    SliceAcqOrder
)


# Globlally define (temporary) log file object
with TmpDir(src=os.getcwd()) as tmpd:
    with TmpDir.TmpFile(tmp_dir=tmpd.src, ext='.log') as tmpf:
        log: LogFile = LogFile(log_file=tmpf.src)


@timeops(log)
def mcdc(func: str,
         outdir: str,
         func_echospacing: Optional[float] = 0.1,
         func_pedir: Optional[Union[str,List[str]]] = None,
         epifactor: Optional[int] = None,
         inplane_acc: float = 1,
         func_brainmask: Optional[str] = None,
         func_slorder: Optional[str] = None,
         fmap: Optional[str] = None,
         fmap2func_affine: Optional[str] = None,
         mb_factor: Optional[int] = None,
         dc: bool = False,
         s2v: bool = False,
         mbs: bool = False,
         use_mcflirt: bool = False,
         ref: Union[str,int] = 0,
         log: Optional[LogFile] = None
        ) -> Tuple[str]:
    """Perform motion and distortion correction stage of ``fmri_preproc``.
    """
    if log:
        log.log("Performing Motion correction")

    # Logic tests
    _has_fmap: bool = fmap is not None
    _has_acqp: bool = (func_pedir is not None) & (func_echospacing is not None)

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
    
    if use_mcflirt:
        _has_fmap: bool = False
        _has_acqp: bool = False
        s2v: bool = False
        mbs: bool = False
        dc: bool = False
    else:
        use_mcflirt: bool = (not dc) & (not s2v) & (not mbs)

    # Check argument combinations
    if (dc & use_mcflirt) | (s2v & use_mcflirt) | (mbs & use_mcflirt):
        if log: log.warning("WARNING: MCFLIRT enabled. This has disabled dc, s2v and/or mbs.")
        warn("WARNING: MCFLIRT enabled. This has disabled dc, s2v and/or mbs.")

    if dc and not (_has_fmap & _has_acqp):
        if log: log.error("RuntimeError: fmap, fmap2func_affine, func_pedir, and func_echospacing are required for DC.")
        raise RuntimeError('fmap, fmap2func_affine, func_pedir, and func_echospacing are required for DC.')

    if mbs and not (_has_fmap & _has_acqp):
        if log: log.error("RuntimeError: fmap, fmap2func_affine, func_pedir, and func_echospacing are required for mbs.")
        raise RuntimeError('fmap, fmap2func_affine, func_pedir, and func_echospacing are required for mbs.')

    if mbs and not dc:
        if log: log.error("RuntimeError: Cannot mbs without dc as they are both distortion corrections.")
        raise RuntimeError('Cannot mbs without dc as they are both distortion corrections.')

    if (not use_mcflirt) & (func_brainmask is None):
        if log: log.error("RuntimeError: func_brainmask is required to use EDDY.")
        raise RuntimeError('func_brainmask is required to use EDDY.')
    
    # Define output files
    out: MCDCFiles = MCDCFiles(outdir=outdir)
    outputs: Dict[str,str] = out.outputs(dc=(not use_mcflirt))
    mcdir: str = outputs.get('mcdir')

    # Perform MCDC
    if use_mcflirt:
        (func_mcdc, 
         motfile, 
         _) = mcflirt_mc(func=func,
                         func_mc=outputs.get('func_mcdc'),
                         ref=ref,
                         log=log)
        eddy_output_mask: str = ""
    else:
        (func_mcdc,
         motfile,
         eddy_output_mask) = eddy_mcdc(func=func,
                                       func_brainmask=func_brainmask,
                                       func_mcdc=outputs.get('func_mcdc'),
                                       func_sliceorder=func_slorder,
                                       func_echospacing=func_echospacing,
                                       func_pedir=func_pedir,
                                       epifactor=epifactor,
                                       inplane_acc=inplane_acc,
                                       fmap=fmap,
                                       fmap2func_xfm=fmap2func_affine,
                                       mb_factor=mb_factor,
                                       mot_params=outputs.get('motparams'),
                                       mbs=mbs,
                                       s2v_corr=s2v,
                                       log=log)

    # Write motion regressors
    mcf: np.array = np.loadtxt(motfile)
    mcf: pd.DataFrame = pd.DataFrame(mcf, columns=['RotX', 'RotY', 'RotZ', 'X', 'Y', 'Z'])
    mcf.to_csv(outputs.get('motparams'), sep='\t', index=None)
    mcf: str = outputs.get('motparams')
    
    # Calculate post-mc motion outliers
    _: Tuple[str,int] = motion_outlier(func=func_mcdc,
                                       metric_name=outputs.get('func_metrics'),
                                       plot_name=outputs.get('func_out_plot'))
    
    with File(src=outputs.get('func_metrics'), assert_exists=True) as fnm:
        func_metrics: str = fnm.abspath()

    # Brain extract mcdc images
    mcdc: str = func_mcdc
    mcdc_mean: str = outputs.get('mcdc_mean')
    mcdc_std: str = outputs.get('mcdc_std')
    mcdc_tsnr: str = outputs.get('mcdc_tsnr')
    mcdc_brainmask: str = outputs.get('mcdc_brainmask')

    mcdc_mean: str = fslmaths(img=mcdc).Tmean().run(out=mcdc_mean, log=log)
    mcdc_std: str = fslmaths(img=mcdc).Tmean().run(out=mcdc_std, log=log)
    mcdc_tsnr: str = fslmaths(img=mcdc).Tmean().run(out=mcdc_tsnr, log=log)

    with TmpDir(src=mcdir) as tmp:
        brain: str = os.path.join(tmp.src,'brain.nii.gz')
        brain, _ = bet(img=mcdc_mean, out=brain, mask=False, frac_int=0.4, robust=True, log=log)
        mcdc_brainmask: str = fslmaths(img=brain).bin().run(out=mcdc_brainmask, log=log)
    
    # Create out-of-FOV masks (for EDDY)
    if os.path.exists(eddy_output_mask):
        fov_mask: str = outputs.get('func_mcdc_fovmask')
        fov_percent: str = outputs.get('func_mcdc_fovpercent')

        fov_mask: str = fslmaths(img=mcdc_brainmask).sub(eddy_output_mask).bin().run(out=fov_mask, log=log)
        fov_percent: str = fslmaths(img=fov_mask).Tmean().mul(100).run(out=fov_percent, log=log)
    else:
        fov_mask: str = None
        fov_percent: str = None
    
    return (func_mcdc,
            mcf,
            func_metrics,
            mcdc_mean,
            mcdc_std,
            mcdc_tsnr,
            mcdc_brainmask,
            fov_mask,
            fov_percent)


@timeops(log)
def mcflirt_mc(func: str,
               func_mc: str,
               ref: Optional[Union[int, str]] = None,
               dc_warp: Optional[str] = None,
               log: Optional[LogFile] = None
              ) -> Tuple[str,str,str]:
    """Performs MCFLIRT-based motion and distortion correction.
    """
    func: NiiFile = NiiFile(src=func, assert_exists=True, validate_nifti=True)
    
    if log:
        log.log("Performing motion correction using MCFLIRT.")

    with File(src=func_mc, assert_exists=False) as f:
        outdir, _, _ = f.file_parts()
        func_mc: str = f.abspath()
        with WorkDir(src=outdir) as _:
            pass

    if isinstance(ref, str):
        with NiiFile(src=ref, assert_exists=True, validate_nifti=True) as f:
            ref: str = f.abspath()
    
    if dc_warp:
        with NiiFile(src=dc_warp, assert_exists=True, validate_nifti=True) as f:
            dc_warp: str = f.abspath()
    
    ref_vol: Union[int, None] = ref if isinstance(ref, int) else None
    ref_file: Union[str, None] = ref if isinstance(ref, str) else None

    (func_mc,
     parfile, 
     matsdir) = mcflirt(infile=func.abspath(),
                        outfile=func_mc,
                        reffile=ref_file,
                        refvol=ref_vol,
                        log=log)

    if dc_warp:
        with TmpDir(src=outdir) as tmp:
            func_ref: str = fslroi(img=func_mc,
                                   out=os.path.join(tmp.abspath(),"func0_ref.nii.gz"),
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


@timeops(log)
def eddy_mcdc(func: str,
              func_brainmask: str,
              func_mcdc: str,
              func_sliceorder: Optional[str] = None,
              func_echospacing: Optional[float] = 0.1,
              func_pedir: Optional[Union[str,List[str]]] = None,
              epifactor: Optional[int] = None,
              inplane_acc: float = 1,
              fmap: Optional[str] = None,
              fmap2func_xfm: Optional[str] = None,
              mb_factor: Optional[int] = 1,
              mot_params: Optional[str] = None,
              mbs: bool = False,
              mporder: Optional[int] = None,
              s2v_corr: bool = False,
              log: Optional[LogFile] = None
             ) -> Tuple[str,str,str]:
    """Perform EDDY-based motion and distortion correction.

    NOTE: Input ``fmri`` acquisition is **ASSUMED** to be acquired in the same phase-encoding direction
        throughout each volume (e.g. the phase-encoding direction cannot change from PA -> AP to AP -> PA
        for several volumes).
    """
    if log: log.log("Performing EDDY-based motion and distortion correction.")

    with NiiFile(src=func, assert_exists=True, validate_nifti=True) as fn:
        with NiiFile(src=func_brainmask, assert_exists=True, validate_nifti=True) as fb:
            with File(src=func_mcdc) as fmc:
                func: str = fn.abspath()
                func_brainmask: str = fb.abspath()

                outdir, _, _ = fmc.file_parts()
                eddy_dir: str = os.path.join(outdir,"eddy")
                eddy_basename: str = os.path.join(eddy_dir,"eddy_corr")

                with WorkDir(src=eddy_dir) as _:
                    if log: log.log("Creating eddy output directory.")
    
    # Define output files
    outputs: Dict[str,str] = {
                                "idx": eddy_basename + "_fmri_pre-mcdc.idx",
                                "bvals": eddy_basename + "_fmri_pre-mcdc.bval",
                                "bvecs": eddy_basename + "_fmri_pre-mcdc.bvec",
                                "acqp": eddy_basename + "_fmri_pre-mcdc.acqp",
                                "slice_order": eddy_basename + "_fmri_pre-mcdc.slice.order",
                                "ident_matrix": os.path.join(FSLDIR,'etc', 'flirtsch', 'ident.mat')
                             }

    # Logic tests
    _has_acqp = (func_pedir is not None) & (func_echospacing is not None)
    _has_fmap = fmap is not None

    if _has_fmap and not _has_acqp:
        raise RuntimeError('Cannot do DC without func_echospacing.')

    if mbs and not _has_fmap:
        raise RuntimeError('Cannot do MBS without fmap.')

    if s2v_corr and log:
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

    idx: str = write_index(num_frames=num_vols, out_file=outputs.get('idx'))
    bvals: str = write_bvals(num_frames=num_vols, out_file=outputs.get('bvals'))
    bvecs: str = write_bvecs(num_frames=num_vols, out_file=outputs.get('bvecs'))
    acqp, _ = write_func_params(epi=func,
                                  echospacing=func_echospacing,
                                  pedir=PhaseEncodeDirection(func_pedir).name,
                                  out=outputs.get('acqp'),
                                  epifactor=epifactor,
                                  inplane_acc=inplane_acc)

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

    # TODO: mb_factor needs to be required in elif statement
    #   to ensure reliable performance.
    if s2v_corr:
        if func_sliceorder:
            with File(src=func_sliceorder, assert_exists=True) as f:
                # func_sliceorder: str = f.abspath()
                func_sliceorder: str = f.copy(dst=outputs.get('slice_order'))
        else:
            func_sliceorder: str = write_slice_order(slices=slices, 
                                                     mb_factor=mb_factor, 
                                                     mode="single-shot",
                                                     out_file=outputs.get('slice_order'),
                                                     return_mat=False)

        s2v_niter:  int = 10
        s2v_fwhm:   int = 0
        s2v_lambda: int = 1
        s2v_interp: str = "trilinear"

        # TODO: Add option to arbitrarily set the mporder.
        # 
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
        with File(src=fmap2func_xfm, assert_exists=True) as f:
            fmap2func_xfm: str = f.abspath()
    else:
        fmap2func_xfm: str = outputs.get('ident_matrix')
    
    # Prepare fieldmap by converting fieldmap from rad/s -> Hz
    field_hz: Union[str, None] = None
    if fmap:
        with NiiFile(src=fmap, assert_exists=True, validate_nifti=True) as f:
            _, fname, _ = f.file_parts()
            field_hz: str = os.path.join(eddy_basename + "_" + fname + "_Hz_pre-mcdc")
            field_hz: str = fslmaths(img=fmap).div(6.2832).run(out=field_hz, log=log)
        
    # Perform Eddy-based mcdc
    (eddy_corr,
     eddy_motion_par,
     eddy_mask) = eddy(img=func,
                       out=eddy_basename,
                       mask=func_brainmask,
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

    # Necessary eddy output files
    with NiiFile(src=eddy_corr, assert_exists=True, validate_nifti=True) as f:
        func_mcdc: str = f.copy(dst=func_mcdc)

    if mot_params:
        with File(src=eddy_motion_par) as f:
            mot_params: str = f.sym_link(mot_params)
    
    # Re-write TR in output NIFTI file header
    func: nib.Nifti1Image = nib.load(func)

    nib.Nifti1Image(nib.load(func_mcdc).get_fdata(),
                    header=func.header,
                    affine=func.affine
                    ).to_filename(func_mcdc)

    return (func_mcdc,
            mot_params,
            eddy_mask)


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


def _dvars(img: nib.Nifti1Image) -> np.array:
    """Calculates DVARS.
    """
    dvars: np.array = np.square(np.diff(img, axis=3))
    dvars: np.array = np.nanmean(dvars, axis=(0, 1, 2))
    dvars: np.array = np.sqrt(dvars)
    dvars: np.array = 1000 * (dvars / np.nanmedian(img))
    dvars: np.array = np.concatenate(([0], dvars))
    return dvars


def _refrms(img: nib.Nifti1Image,
            ref: Optional[nib.Nifti1Image] = None
           ) -> np.array:
    """Calculates REFRMS.
    """
    if ref:
        pass
    else:
        ref: nib.Nifti1Image = img[:, :, :, int(np.round(img.shape[3] / 2))]

    rms: nib.Nifti1Image = img
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
                   metric_name: Optional[str] = None,
                   plot_name: Optional[str] = None,
                   metric: str = "dvars",
                   ref: Optional[str] = None,
                   thr: Optional[int] = None
                  ) -> Tuple[str,int,Union[str,None]]:
    """ Estimate motion in 4D volume.
    """
    func: NiiFile = NiiFile(src=func, assert_exists=True, validate_nifti=True)

    img0: nib.Nifti1Image = nib.load(func.src)
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
    if metric_name:
        pd.DataFrame(
            np.stack((dvars, refrms, outlier), axis=1),
            columns=['DVARS', 'RefRMS', 'Outlier' + metric.upper()],
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

        ax_hdl.plot(metric_data, label=metric)
        ax_hdl.set_xlim(1, len(metric_data))
        ax_hdl.set_ylim(0, np.nanmean(metric_data) + (4 * np.nanstd(metric_data)))
        ax_hdl.set_ylabel(metric)
        ax_hdl.set_title(fname)
        ax_hdl.set_xlabel("time (volumes)")
        plt.legend()
        plt.savefig(plot_name)

    return outlier, metric_data, thr, metric_name, plot_name
