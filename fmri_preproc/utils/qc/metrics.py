# -*- coding: utf-8 -*-
"""Quality control metrics for the ``fmri_preproc`` neonatal rs-fMRI preprocessing pipeline.
"""
import io
import os
import warnings
from tempfile import TemporaryDirectory as tempdir

from typing import (
    Dict,
    Optional, 
    Tuple, 
    Union
)

import nibabel as nib
import numpy as np
import pandas as pd
from nibabel.nifti1 import Nifti1Image as Image
from numpy import random
from statsmodels.tsa.ar_model import AR

from fmri_preproc.utils.qc.util import (
    run,
    split
)

def convert_to_z(d: np.ndarray, 
                 robust: bool = True, 
                 axis: Optional[Union[int,Tuple[int]]] = None
                ) -> np.ndarray:
    """doc-string
    """
    if robust:
        #   modified z-score
        #   more robust to outliers
        median = np.nanmedian(d, axis=axis, keepdims=True)
        diff = d - median
        # median absolute deviation (MAD)
        # https://en.wikipedia.org/wiki/Median_absolute_deviation)
        mad = np.nanmedian(np.abs(diff), axis=axis, keepdims=True)

        oldsettings = np.seterr(all='warn')
        np.seterr(divide='ignore', invalid='ignore')

        z = 0.6745 * diff / mad
        z[np.isnan(z)] = 0

        np.seterr(**oldsettings)

    else:
        mu = np.mean(d, axis=axis, keepdims=True)
        sigma = np.std(d, ddof=0, axis=axis, keepdims=True)
        z = (d - mu) / sigma

    return z


def dr(func: Image,
       func_brainmask: Image,
       spatial_template: Image,
       template2func_warp: Image,
       basename: str,
       standard: Image = None,
       tmpdir: Optional[str] = None,
       func2standard_warp: Image = None,
       **kwargs):
    """Performs dual regression.
    """
    # create output template dictionary
    ft: Dict[str,str] = dict(
        dr1='{basename}_dr1.txt',
        dr2='{basename}_dr2.nii.gz',
        dr2_standard='{basename}_standard_dr2.nii.gz',
        cnr='{basename}_cnr.nii.gz',
        cnr_standard='{basename}_standard_cnr.nii.gz',
    )

    # update template with kwargs
    for k, v in kwargs.items():
        if k not in ft:
            raise ValueError(f'unknown output template {k}')
        ft[k] = v
    
    dname, _, _ = split(fname=v)

    if not os.path.exists(dname):
        os.makedirs(name=dname, exist_ok=True)

    with tempdir(dir=tmpdir) as td:

        # warp spatial_maps to native space
        tmp_template_native = os.path.join(td, 'template_native.nii.gz')
        run([
            'applywarp', '-i', spatial_template, '-o',
            tmp_template_native, '-r', func_brainmask, '-w',
            template2func_warp, '--interp=spline'
        ])

        # DR stage 1: regress spatial maps onto functional timeseries
        # residual == noise
        # dr1 = os.path.join(basename + '_space-orig_dr1.txt')
        tmp_noise = os.path.join(td, 'dr_stage1_noise.nii.gz')
        run([
            'fsl_glm', '-i', func, '-d',
            tmp_template_native, '-o', ft.get('dr1'),
            '--demean', '-m', func_brainmask,
            '--out_res=' + tmp_noise
        ])

        # DR stage 2: regress stage 1 timeseries onto functional timeseries
        # dr2 = os.path.join(basename + '_space-orig_dr2.nii.gz')
        run([
            'fsl_glm', '-i', func, '-d', ft.get('dr1'), '-o', ft.get('dr2'),
            '--demean', '-m', func_brainmask, '--des_norm'
        ])
        dr2 = nib.load(ft.get('dr2'))

        dr_stats = img_descriptives(dr2, func_brainmask, prefix='dr2')

        # correlate func timeseries with dr2 timeseries
        dr_stats['spatial_corr'] = tscc(dr2, nib.load(tmp_template_native),
                                        func_brainmask)

        # calculate netmat
        nm, nmz = netmat(np.loadtxt(ft.get('dr1')).T)

        dr_stats['netmat'] = nm
        dr_stats['netmat_z'] = nmz

        # calculate the temporal stdev of the noise: Tstd(noise)
        tmp_noise_std = os.path.join(td, 'dr_stage1_noise_std.nii.gz')
        run(['fslmaths', tmp_noise, '-Tstd',
             '-mas', func_brainmask, tmp_noise_std])

        # calculate the temporal standard deviation of the contrast:
        # Tstd(func - noise)
        tmp_contrast_std = os.path.join(td, 'contrast_std.nii.gz')
        run([
            'fslmaths', func, '-sub', tmp_noise,
            '-Tstd', '-mas', func_brainmask, tmp_contrast_std
        ])

        # calculate the cnr: contrast_std / noise_std
        # cnr = basename + '_space-orig_cnr.nii.gz'
        run(['fslmaths', tmp_contrast_std, '-div', tmp_noise_std, '-nan', ft.get('cnr', make_dir=True)])
        cnr = nib.load(ft.get('cnr'))
        cnr_stats = img_descriptives(cnr, func_brainmask, prefix='cnr')

        # setup output images and stats
        images = {'dr1': ft.get('dr1'), 'dr2': dr2, 'cnr': cnr}
        stats = {**dr_stats, **cnr_stats}

        # warp to standard space
        if standard is not None and func2standard_warp is not None:
            for f in ['dr2', 'cnr']:
                # fname = images[f].get_filename().replace('_space-orig_',
                #                                          '_space-standard_')
                run([
                    'applywarp', '-i', images[f], '-o', ft.get(f + '_standard'),
                    '-r', standard,
                    '-w', func2standard_warp, '--interp=spline'
                ])
                images[f + '_standard'] = nib.load(ft.get(f + '_standard'))

    return images, stats


def tSNR(func: Image,
         func_brainmask: Image,
         basename: str,
         standard: Image = None,
         func2standard_warp: Image = None,
         **kwargs) -> Image:
    """Calculate temporal SNR."""

    # create output template dictionary
    ft: Dict[str,str] = dict(
        mean='{basename}_tmean.nii.gz',
        mean_standard='{basename}_standard_tmean.nii.gz',
        std='{basename}_tstd.nii.gz',
        std_standard='{basename}_standard_tstd.nii.gz',
        snr='{basename}_tsnr.nii.gz',
        snr_standard='{basename}_standard_tsnr.nii.gz',
    )

    # update template with kwargs
    for k, v in kwargs.items():
        if k not in ft:
            raise ValueError(f'unknown output template {k}')
        ft[k] = v

    dname, _, _ = split(fname=v)

    if not os.path.exists(dname):
        os.makedirs(name=dname, exist_ok=True)

    # with tempdir(dir=tmpdir) as td:

    fn = func.get_filename()

    mean = ft.get('mean')
    std = ft.get('std')

    run(['fslmaths', fn, '-Tmean', mean])
    run(['fslmaths', fn, '-Tstd', std])
    run(['fslmaths', mean, '-div', std, '-nan', ft.get('snr')])
    snr = nib.load(ft.get('snr'))

    stats = img_descriptives(snr, func_brainmask, prefix='snr')

    images = {'snr': snr}

    if standard is not None and func2standard_warp is not None:
        run([
            'applywarp', '-i', snr, '-o', ft.get('snr_standard', make_dir=True),
            '-r', standard,
            '-w', func2standard_warp, '--interp=spline'
        ])
        images['snr_standard'] = nib.load(ft.get('snr_standard'))

    return images, stats


# TODO: Add Q25, Q75, IQR,
def img_descriptives(img: Image, mask: Image, prefix: str = None):
    """Calculate descriptive statistics of image."""
    out = run([
        'fslstats', img, '-k', mask, '-M', '-S',
        '-P', '1', '-P', '5', '-P', '50', '-P', '95', '-P', '99'
    ]).split()

    keys = ['mean', 'std', 'p01', 'p05', 'p50', 'p95', 'p99']
    stats = {}
    for idx, k in enumerate(keys):
        stats[k if prefix is None else prefix + '_' + k] = out[idx]

    return stats


def ts_descriptives(ts: np.ndarray):
    ts_abs = np.abs(ts)

    stats = {}

    stats['rms'] = np.sqrt(np.nanmean(np.square(ts)))
    stats['max'] = np.amax(ts)
    stats['mean'] = np.mean(ts)
    stats['std'] = np.std(ts)

    stats['rms_abs'] = np.sqrt(np.nanmean(np.square(ts_abs)))
    stats['max_abs'] = np.amax(ts_abs)
    stats['mean_abs'] = np.mean(ts_abs)
    stats['std_abs'] = np.std(ts_abs)

    for p in [1, 5, 50, 95, 99]:
        stats['p{:02n}'.format(p)] = np.percentile(ts, p)
        stats['p{:02n}_abs'.format(p)] = np.percentile(ts_abs, p)

    return stats


def tscc(img1: nib.Nifti1Image,
         img2: nib.Nifti1Image,
         brainmask: nib.Nifti1Image):
    """Calculate timeseries cross-correlation from two images."""
    cc = run([
        'fslcc', '-t', '-10', '--noabs', '-m', brainmask, img2, img1
    ])
    cc = np.loadtxt(io.StringIO(cc))
    return cc[cc[:, 0] == cc[:, 1], 2]


def measurecost(src: nib.Nifti1Image,
                ref: nib.Nifti1Image,
                cost: str = 'normmi',
                boundarymask: nib.Nifti1Image = None,
                ref_brainmask: nib.Nifti1Image = None):
    valid_costfcn = ['mutualinfo', 'corratio', 'normcorr', 'normmi', 'leastsq', 'labeldiff', 'bbr']
    assert cost in valid_costfcn, (
        'Invalid costfcn {}.  Must be one of: {}'.format(cost, valid_costfcn)
    )

    fsldir = os.environ.get('FSLDIR')

    assert cost != 'bbr' or boundarymask is not None, (
        'boundary mask required if cost is bbr'
    )

    cmd = 'flirt -in {0} -ref {1} -schedule {2}/etc/flirtsch/measurecost1.sch -cost {3}'
    cmd = cmd.format(src.get_filename(), ref.get_filename(), fsldir, cost)

    if ref_brainmask is not None:
        cmd += ' -refweight {}'.format(ref_brainmask.get_filename())

    if cost == 'bbr' and boundarymask is not None:
        cmd += ' -wmseg {}'.format(boundarymask.get_filename())

    cst = float(run(cmd.split()).split()[0])

    return cst


def netmat(ts, varnorm: bool = True) -> Tuple:
    """
    Calculate netmats.

    Full correlation and (unregularised) partial correlation

    d           Time-series (nodes x time)
    varnorm     If True, temporal variance normalisation (overall stddev)

    """

    def _pc(dd):
        cv = np.cov(dd)
        pinvA = -np.linalg.inv(cv)
        iis = np.tile(np.sqrt(np.abs(pinvA.diagonal()))[:, np.newaxis], pinvA.shape[1])
        return pinvA / (iis * iis.T)

    def _fc(dd):
        return np.corrcoef(dd)

    def _ar1(dd):
        return AR(dd).fit(1).params[1]

    def _r_to_z(r):
        return 0.5 * np.log((1 + r) / (1 - r))

    # demean timeseries
    d = ts - np.mean(ts, axis=1, keepdims=True)

    # variance normalise
    if varnorm:
        d = d / np.std(d.ravel())

    # calculate full correlation
    fc = _fc(d)
    fc[np.diag_indices_from(fc)] = 0

    # calculate partial correlation
    pc = _pc(d)
    pc[np.diag_indices_from(pc)] = 0

    # calculate null distribution
    ar1 = np.median([_ar1(d[i, :]) for i in np.arange(d.shape[0])])

    # TODO: this should not be done as list appends (SLOW)
    null_dist = []
    null_dist.append(random.randn(25, 1))
    for i in np.arange(d.shape[1] - 1):
        null_dist.append(null_dist[i] * ar1 + random.randn(25, 1))

    null_dist = np.concatenate(null_dist, axis=1)

    # calculate fc of null dist
    null_fc = _fc(null_dist)
    null_fc = null_fc[~np.eye(null_fc.shape[0], dtype=bool)]
    null_fc_z = _r_to_z(null_fc)

    # calculate pc of the null dist
    null_pc = _pc(null_dist)
    null_pc = null_pc[~np.eye(null_pc.shape[0], dtype=bool)]
    null_pc_z = _r_to_z(null_pc)

    rtoz_cor = 1 / np.std(null_pc_z)
    # print(rtoz_cor)
    pc_z = _r_to_z(pc) * rtoz_cor

    rtoz_cor = 1 / np.std(null_fc_z)
    # print(rtoz_cor)
    fc_z = _r_to_z(fc) * rtoz_cor

    # copy partial correlation to lower triangle
    li = np.tril_indices(fc.shape[0], k=-1)

    netmat = fc.copy()
    netmat[li] = pc[li].copy()

    netmat_z = fc_z.copy()
    netmat_z[li] = pc_z[li].copy()

    return netmat, netmat_z


def motparams(func: nib.Nifti1Image,
              motparams: pd.DataFrame = None,
              tmpdir: str = None):
    """doc-string.
    """
    if motparams is None:
        with tempdir(dir=tmpdir) as td:
            outname = os.path.join(td, 'mc')
            run(['mcflirt', '-in', func, '-meanvol', '-plots', '-o', outname])
            motparams = np.loadtxt(outname + '.par')
            motparams = pd.DataFrame(motparams, columns=['RotX', 'RotY', 'RotZ', 'X', 'Y', 'Z'])

    stats = {}
    for key in motparams.columns:
        stats[key] = motparams[key].values

    tr_stats = ts_descriptives(motparams[['X', 'Y', 'Z']].values)
    for k, v in tr_stats.items():
        stats['tr_' + k] = v

    rot_stats = ts_descriptives(motparams[['RotX', 'RotY', 'RotZ']].values)
    for k, v in rot_stats.items():
        stats['rot_' + k] = v

    if 'framewise_displacement' in motparams.columns:
        fd_stats = ts_descriptives(motparams['framewise_displacement'].values)
    else:
        fd_stats = fd(motparams[['RotX', 'RotY', 'RotZ', 'X', 'Y', 'Z']].values, order='mcflirt')

    for k, v in fd_stats.items():
        stats[k] = v

    return stats


def dvars(func: nib.Nifti1Image,
          mask: nib.Nifti1Image = None,
          thr: float = None):
    """Calculate DVARS."""

    func = func.get_data().astype(np.float32)

    Nt = func.shape[-1]
    func = func.reshape((-1, Nt), order='F')

    if mask is not None:
        mask = mask.get_data()
        mask = mask.ravel(order='F')
    else:
        th2, th98 = np.percentile(func, [2, 98])
        robthr = th2 + 0.1 * (th98 - th2)
        mask = np.mean(func, axis=1) >= robthr

    func = func[mask == 1, :]

    dvars = np.square(np.diff(func, axis=1))
    dvars = np.mean(dvars, axis=0)
    dvars = np.sqrt(dvars)
    dvars = 1000 * (dvars / np.median(func))
    dvars = np.concatenate(([0], dvars))

    # calculate outlier
    if thr is None:
        q75, q25 = np.percentile(dvars, [75, 25])
        iqr = q75 - q25
        thr = q75 + (1.5 * iqr)
    outlier = dvars > thr

    stats = {'dvars': dvars, 'outlier': outlier}
    for k, v in ts_descriptives(dvars).items():
        stats['dvars_' + k] = v

    return stats


def fd(motparams: np.ndarray, order: str = 'mcflirt'):
    """Calculate framewise displacement (Power et al, 2012)."""

    if order not in ['mcflirt', 'eddy']:
        raise (RuntimeError(f'Invalid order: {order}'))

    if order == 'eddy':
        motparams = np.hstack([motparams[:, 3:], motparams[:, :3]])

    rot = motparams[:, :3]
    rot = np.diff(rot, axis=0)
    rot = np.abs(rot) * 50

    tr = motparams[:, 3:]
    tr = np.diff(tr, axis=0)
    tr = np.abs(tr)

    fd = np.sum(np.concatenate((rot, tr), axis=1), axis=1)

    stats = {'fd': fd}
    for k, v in ts_descriptives(fd).items():
        stats['fd_' + k] = v

    return stats


def img_cumstats(img, ddof=0):
    N = len(img)

    d0 = img[0].get_data()
    img[0].uncache()

    for i in np.arange(1, N):
        d0 += img[i].get_data()
        img[i].uncache()

    mn = d0 / N

    diffsq = (img[0].get_data() - mn) ** 2
    img[0].uncache()

    for i in np.arange(1, N):
        diffsq += (img[i].get_data() - mn) ** 2
        img[i].uncache()

    std = np.sqrt((diffsq / (N - ddof)))

    stats = {
        'mean': nib.Nifti1Image(mn, img[0].affine, img[0].header),
        'std': nib.Nifti1Image(std, img[0].affine, img[0].header)
    }

    return stats


def z_smoothness(img, brain_thr=20):
    """Estimate smoothness in IS direction."""
    # ax = acqp._get_axis(img, 'IS')[0]
    #
    # assert ax is 'z', ("z_smoothness currently requires IS to be in the",
    #                    " z dimension")

    # load img
    d = nib.load(img)
    d0 = np.copy(d.get_data())

    # threshold at brain_thr to get rid of non-brain
    # TODO: replace with EM gaussian mixture model
    d0[d0 < brain_thr] = np.nan

    # calculate mean intensity xy plane for each volume
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        smth = np.nanmean(np.nanmean(d0, axis=0), axis=0)

    # estimate smoothness across z dimension
    smth = np.nanstd(np.diff(smth, axis=0), axis=0)

    return smth.astype(float)
