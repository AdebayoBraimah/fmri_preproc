# -*- coding: utf-8 -*-
"""
"""
import os
import numpy as np

from typing import (
    List,
    Tuple,
    Optional,
    Union
)

from fmri_preproc.utils.enums import SliceAcqOrder

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

def generate_acq_params(num_frames: int,
                        effective_echo_spacing: float = 0.05,
                        out_prefix: str = 'file'
                       ) -> Tuple[str,str]:
    """Creates acquisition parameters files for fMRI acquisition provided the number of frames and the 
    effective echo spacing.
    
    Usage example:
        >>> dist_acqp, func_acqp = generate_acq_params(num_frames=1200, 
        ...                                            effective_echo_spacing=0.05,
        ...                                            out_prefix="fmri")
        ...
        
    Arguments:
        num_frames: Number of temporal frames/dynamics.
        effective_echo_spacing: Effective echo spacing.
        out_prefix: Output file prefix.
        
    Returns:
        Tuple of strings that correspond to the acquisition parameter files (for both distortion correction and
            eddy current correction).
    
    Raises:
        TypeError: Exception that is raised when either ``num_frames`` is not an ``int`` OR when ``effective_echo_spacing`` is not a ``float``.
    """
    if (not isinstance(effective_echo_spacing, int) or 
        not isinstance(effective_echo_spacing, float) or
        not isinstance(num_frames, int)):
        raise TypeError(f"Input for num_frames: {num_frames} is not an integer OR effective_echo_spacing: {effective_echo_spacing} is not a float.")
    
    # Generate distortion correction
    out_dist: str = out_prefix + "_distortion_correction.acqp"
    out_func: str = out_prefix + "_functional.acqp"
    
    # Write distortion correction acqp file
    with open(out_dist,'w') as f:
        f.write(f"0 1 0 {effective_echo_spacing}\n")
        f.write(f"0 -1 0 {effective_echo_spacing}\n")
        f.close()
    
    # Write functinoal acqp file
    with open(out_func,'w') as f:
        for i in range(0,num_frames):
            f.write(f"0 1 0 {effective_echo_spacing}\n")
        f.close()
    
    # Obtain absolute file paths
    out_dist: str = os.path.abspath(out_dist)
    out_func: str = os.path.abspath(out_func)
    
    return (out_dist, 
            out_func)

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
        Links for details:
            * Documentation: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#A--slspec
            * Implementation: https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/blob/master/dhcp/resources/default_func.slorder
            * Implementation: https://git.fmrib.ox.ac.uk/matteob/dHCP_neo_dMRI_pipeline_release/-/blob/master/slorder.txt
    
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
