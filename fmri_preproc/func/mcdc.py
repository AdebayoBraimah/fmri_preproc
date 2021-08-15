# -*- coding: utf-8 -*-
"""Performs motion correction and distortion correction of MR image data.

NOTE:
    * Check this resource pertaining to the odd number of slices: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;28f6c983.1806
"""
from math import inf
import os
import numpy as np
import nibabel as nib
from numpy.lib.npyio import load
import pandas as pd

from typing import (
    List,
    Optional,
    Union
)

from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.func.fieldmap import _get_b0_conf

from fmri_preproc.utils.fslpy import (
    eddy,
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

def mcdc():
    """doc-string
    """
    pass

def mcflirt_mc(func: str,
               func_mc: str,
               ref: Optional[Union[int, str]] = None,
               dc_warp: Optional[str] = None,
               log: Optional[LogFile] = None
              ) -> None:
    """Performs MCFLIRT-based motion and distortion correction.
    """
    func: NiiFile = NiiFile(file=func, assert_exists=True, validate_nifti=True)
    
    with File(file=func_mc, assert_exists=False) as f:
        outdir, _, _ = f.file_parts()
        with WorkDir(work_dir=outdir) as d:
            d.mkdir()

    func_mc: File = File(file=func_mc, assert_exists=False)

    if isinstance(ref, str):
        with NiiFile(file=ref, assert_exists=True, validate_nifti=True) as f:
            ref: str = f.abs_path()
    
    if dc_warp:
        with NiiFile(file=dc_warp, assert_exists=True, validate_nifti=True) as f:
            dc_warp: str = f.abs_path()
    
    ref_vol: Union[int, None] = ref if isinstance(ref, int) else None
    ref_file: Union[str, None] = ref if isinstance(ref, str) else None

    par_file, matsdir = mcflirt(infile=func.abs_path(),
                                outfile=func_mc,
                                reffile=ref_file,
                                refvol=ref_vol,
                                log=log)

    return par_file, matsdir

def eddy_mcdc():
    """doc-string
    """
    pass

def generate_bvals(num_frames: int,
                   out_file: str = 'file.bval'
                  ) -> str:
    """Creates/generates an arbitrary number b=0 b-values for a fMRI acquisition
    provided the number of dynamics/temporal frames.
    
    Usage example:
        >>> bval_file = generate_bvals(num_frames=1200, out_file="fmri.bval")
        
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

def generate_bvecs(num_frames: int,
                   out_file: str = 'file.bvec'
                  ) -> str:
    """Creates/generates an arbitrary number of x,y,z b-vectors for a fMRI acquisition
    provided the number of dynamics/temporal frames.
    
    Usage example:
        >>> bval_file = generate_bvecs(num_frames=1200, out_file="fmri.bvec")
        
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

def generate_index(num_frames: int,
                   out_file: str = 'file.idx'
                  ) -> str:
    """Creates index files for use with FSL's ``eddy``.
    
    Usage example:
        >>> func_idx = generate_index(num_frames=1200, out_file="fmri.idx")
        
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
    return np.savetxt(out_file, np.ones((1, num_frames)).T, fmt="%i")

def generate_slice_order(slices: int,
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
        >>> sls_order = generate_slice_order(slices=44,
        ...                                  mb_factor=6,
        ...                                  mode='single-shot',
        ...                                  out_file='file.slice.order',
        ...                                  return_mat=False)
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
    if not isinstance(slices, int) or not isinstance(mb_factor, int):
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

    if metric is 'dvars':
        metric_data: np.array = dvars
    elif metric is 'refrms':
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
