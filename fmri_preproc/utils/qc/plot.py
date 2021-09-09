# -*- coding: utf-8 -*-
"""Quality control plots for the ``fmri_preproc`` neonatal rs-fMRI preprocessing pipeline.
"""
from token import OP
import matplotlib as mpl

mpl.use('Agg')
import matplotlib.pyplot as plt

import matplotlib.ticker as ticker
import nibabel as nib
import numpy as np
import base64
import seaborn as sns
import nilearn.plotting as plotting
import ptitprince as pt

from tempfile import TemporaryFile
from matplotlib.colors import Colormap

from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
    Union
)


def colorbar(ax: plt.subplot, 
             cmap: Union[Colormap,str], 
             vmin: int = 0, 
             vmax: int = 1, 
             color: str = 'k', 
             figsize: Tuple[int,int] = (1, 5), 
             orientation: str = 'vertical', 
             label: str = ''
            ) -> None:
    """doc-string
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=figsize)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = ax
    ax.set_visible(False)
    cb = plt.colorbar(sm, orientation=orientation, fraction=1)
    cb.set_label(label, color=color)
    cb.ax.yaxis.set_tick_params(color=color)
    cb.ax.xaxis.set_tick_params(color=color)
    cb.outline.set_edgecolor(color)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color=color)
    plt.setp(plt.getp(cb.ax.axes, 'xticklabels'), color=color)
    return None


def raincloudplot(d0: np.ndarray, 
                  ax: Optional[plt.subplots] = None, 
                  xticklabels: Optional[str] = None, 
                  xtickrot: Optional[int] = None, 
                  xlabel: Optional[str] = None, 
                  ylabel: Optional[str] = None, 
                  title: Optional[str] = None, 
                  style: str = 'whitegrid', 
                  figsize: Tuple[int,int] = (12, 12), 
                  marker: Optional[str] = None, 
                  z: bool = False,
                  markercolor: str = 'r', 
                  markersize: int = 100, 
                  orient: str = 'v', 
                  markerstyle: str = 'x', 
                  markerlegend: Optional[str] = None, 
                  ylim: Optional[int] = None,
                  **kwargs
                 ) -> None:
    """doc-string
    """
    if z:
        mu = np.mean(d0)
        sigma = np.std(d0, ddof=0)

        d0 = (d0 - mu) / sigma

        marker = None if marker is None else (marker - mu) / sigma

    with sns.axes_style(style=style):

        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=figsize)

        pt.RainCloud(
            data=d0, ax=ax, palette="pastel", orient=orient, **kwargs
        )
        # sns.despine(left=True, bottom=True, ax=ax)

        if marker is not None:

            sc = []

            for idx, c in enumerate(d0.columns):
                y = marker[c].values
                x = np.ones(np.array(y).shape) * idx

                if orient == 'h':
                    tmp = x
                    x = y
                    y = tmp

                sc += [ax.scatter(x, y, s=markersize, marker=markerstyle, c=markercolor, linewidths=1, zorder=100)]

            if markerlegend is not None:
                ax.legend(sc, markerlegend)

        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)
        if xtickrot is not None:
            for l in ax.get_xticklabels():
                l.set_rotation(xtickrot)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if title is not None:
            ax.set_title(title)
        if ylim is not None:
            ax.set_ylim(ylim)
    return None


def distplot(d0: np.ndarray,
             ax: Optional[plt.subplots] = None, 
             xlabel: Optional[str] = None, 
             ylabel: Optional[str] = None, 
             title: Optional[str] = None, 
             style: str = 'darkgrid',
             figsize: Tuple[int,int] = (12, 12),
             **kwargs
            ) -> None:
    """doc-string
    """
    with sns.axes_style(style=style):

        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=figsize)
        sns.distplot(d0, ax=ax, kde=False, norm_hist=False, **kwargs)

        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if title is not None:
            ax.set_title(title)
    return None


def barplot(d0: np.ndarray,
            ax: Optional[plt.subplots] = None, 
            xticklabels: Optional[str] = None, 
            xtickrot: Optional[int] = None, 
            ylabel: Optional[str] = None, 
            title: Optional[str] = None, 
            xlabel: Optional[str] = None, 
            style: str = 'darkgrid',
            figsize: Tuple[int,int] = (12, 12),
            **kwargs
           ) -> None:
    """doc-string
    """
    with sns.axes_style(style=style):

        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=figsize)
        sns.barplot(data=d0, ax=ax, **kwargs)

        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)
        if xtickrot is not None:
            for l in ax.get_xticklabels():
                l.set_rotation(xtickrot)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if title is not None:
            ax.set_title(title)
    return None


def violinplot(d0: np.ndarray, 
               ax: Optional[plt.subplots] = None, 
               xticklabels: Optional[str] = None, 
               xtickrot: Optional[str] = None, 
               xlabel: Optional[str] = None, 
               ylabel: Optional[str] = None, 
               title: Optional[str] = None, 
               style: str = 'darkgrid', 
               figsize: Tuple[int,int] = (12, 12),
               marker: Optional[str] = None, 
               z: bool = False,
               markercolor: str = '#FFFF00', 
               markersize: int = 250, 
               **kwargs
              ) -> None:
    """doc-string
    """
    if z:
        mu = np.mean(d0)
        sigma = np.std(d0, ddof=0)

        d0 = (d0 - mu) / sigma

        marker = None if marker is None else (marker - mu) / sigma

    with sns.axes_style(style=style):

        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=figsize)

        sns.violinplot(
            data=d0, ax=ax, width=.5, palette="pastel",
            linewidth=1, inner="point", **kwargs
        )
        sns.despine(left=True, bottom=True, ax=ax)

        if marker is not None:

            for idx, c in enumerate(d0.columns):
                y = marker[c].values
                x = np.ones(np.array(y).shape) * idx
                ax.scatter(x, y, s=markersize, marker='*', c=markercolor, linewidths=1)

        if xticklabels is not None:
            ax.set_xticklabels(xticklabels)
        if xtickrot is not None:
            for l in ax.get_xticklabels():
                l.set_rotation(xtickrot)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if title is not None:
            ax.set_title(title)
    return None


def motion_parameters(mp: np.ndarray,
                      figsize: Tuple[int,int] = (22, 8),
                      ax: Optional[plt.subplots] = None, 
                      title: str = 'Motion Parameters',
                      xlabel: str = 'TRs'
                     ) -> None:
    """Plot motion parameter timeseries.
    """
    mp = np.array(mp)

    rot = mp[:, :3]
    tr = mp[:, 3:]

    nTR = mp.shape[0]

    if ax is None:
        _, (ax0, ax1) = plt.subplots(2, 1, figsize=figsize)

    # plot rotations
    ln0 = ax0.plot(rot)
    ax0.legend(['X', 'Y', 'Z'], loc=1)
    ax0.grid()
    ax0.set_ylabel('Rotation (radians)')
    ax0.set_xlim(0, nTR)
    ax0.set_xticklabels([])

    # plot translations
    ln1 = ax1.plot(tr)
    ax1.legend(['X', 'Y', 'Z'], loc=1)
    ax1.grid()
    ax1.set_ylabel('Translation (mm)')
    ax1.set_xlim(0, nTR)

    if xlabel is not None:
        ax1.set_xlabel(xlabel)
    if title is not None:
        ax0.set_title(title)
    return None


def framewise_displacement(fd: np.ndarray,
                           figsize: Tuple[int,int] = (22, 4),
                           ax: Optional[plt.subplots] = None, 
                           title: str = 'Framewise Displacement',
                           ylabel: str = 'mm',
                           xlabel: str = 'TRs'
                          ) -> None:
    """Plot framewise displacement of the timeseries.
    """
    fd = np.array(fd)

    nTR = fd.shape[0]

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=figsize)

    # plot fd
    ln = ax.plot(fd)
    ax.grid()
    ax.set_xlim(0, nTR)

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if title is not None:
        ax.set_title(title)
    return None


def dvars(dvars: np.ndarray,
          figsize: Tuple[int,int] = (22, 4),
          ax: Optional[plt.subplots] = None, 
          title: str = 'DVARs',
          ylabel: str = 'DVARs',
          xlabel: str = 'TRs',
          legend: Optional[str] = None
         ) -> None:
    """Plot DVARS of the timeseries.
    """
    dvars = np.array(dvars)
    nTR = dvars.shape[0]

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=figsize)

    # plot dvars
    ln = ax.plot(dvars)
    ax.grid()
    ax.set_xlim(0, nTR)

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if title is not None:
        ax.set_title(title)
    if legend is not None:
        ax.legend(legend)
    return None


def netmat(netmat: np.ndarray,
           figsize: Tuple[int,int] = (12, 10),
           ax: Optional[plt.subplots] = None, 
           cmap: Union[Colormap,str] = 'bwr',
           title: str = 'netmat',
           ylabel: Optional[str] = None,
           xlabel: Optional[str] = None,
           colorbar: bool = True,
           colorbarlabel: Optional[str] = None,
           cbar_shrink: int = 0.75,
           clim: Tuple[int,int] = (-3, 3)
          ) -> None:
    """Plot network matrix (netmat).
    """
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=figsize)

    mp = ax.matshow(netmat)
    mp.set_clim(clim[0], clim[1])
    mp.set_cmap(cmap)

    # Add colorbar
    if colorbar:
        cbar = plt.colorbar(mp, ax=ax, shrink=cbar_shrink)
        if colorbarlabel is not None:
            cbar.set_label(colorbarlabel)

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if title is not None:
        ax.set_title(title)
    return None


def spatialcor(sc: np.ndarray,
               figsize: Tuple[int,int] = (12, 10),
               ax: Optional[plt.subplots] = None, 
               title: str = 'Correlation with spatial template',
               ylabel: str = 'Correlation',
               xlabel: Optional[str] = None,
               xticklabels: Optional[str] = None,
               xtickrot: Optional[int] = None,
               ymin: int = 0,
               ymax: int = 1
              ) -> None:
    """Plot spatial correlation as barplot.
    """
    sc = np.array(sc)

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=figsize)

    ax.set_axisbelow(True)
    ax.grid()

    bar = ax.bar(np.arange(np.max(sc.shape)), sc)
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(np.arange(sc.shape[0]))

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if title is not None:
        ax.set_title(title)

    if xticklabels is not None:
        ax.set_xticklabels(xticklabels)

    if xtickrot is not None:
        for l in ax.get_xticklabels():
            l.set_rotation(xtickrot)
    return None


def voxplot(img: nib.Nifti1Image,
            dseg: Optional[nib.Nifti1Image] = None,
            dseg_labels: Optional[Dict[str,int]] = None,
            brain_mask: Optional[nib.Nifti1Image] = None,
            ax: Optional[plt.subplots] = None, 
            figsize: Tuple[int,int] = (12, 6),
            title: Optional[str] = None,
            xlabel: str = 'TRs',
            ylabel: str = 'Voxels',
            colorbar: bool = True,
            colorbarlabel: Optional[str] = None,
            cmap: str = 'bwr',
            clim: Optional[float] = None,
            cbar_shrink: float = 0.75,
            zscore: bool = False,
            grid_color: str = 'w',
            grid_style: str = ':',
            grid_width: int = 2,
            img_mean: Optional[nib.Nifti1Image] = None,
            img_std: Optional[nib.Nifti1Image] = None,
            interp: str = 'nearest'
           ) -> None:
    """doc-string
    """
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    img = img.get_data()
    img = np.reshape(img, (-1, img.shape[-1]))

    if dseg is not None and dseg_labels is not None:
        idx = []
        count = [0]
        labels = []

        dseg = dseg.get_data().ravel()
        for k, m in dseg_labels.items():
            labels += [k]
            idx0 = list(np.where(np.isin(dseg, m))[0])
            idx += idx0
            count += [int(len(idx0))]
        img = img[idx, :]
    elif brain_mask is not None:
        brain_mask = brain_mask.get_data().ravel().astype(bool)
        img = img[brain_mask, :]

    # Regress the temporal mean and linear trend as per Power
    r0 = np.ones((2, img.shape[1]))
    r0[1, :] = np.linspace(0, 1, img.shape[1])

    beta = np.dot(np.linalg.pinv(r0.T), img.T)
    pred = np.dot(r0.T, beta)
    img = img - pred.T

    if zscore:
        img_mean = np.mean(img) if img_mean is None else img_mean
        img_std = np.std(img) if img_std is None else img_std
        img = (img - img_mean) / img_std

    aximg = ax.imshow(img, aspect='auto', interpolation=interp)
    aximg.set_cmap(cmap)

    if clim is not None:
        aximg.set_clim(clim)

    if dseg is not None and dseg_labels is not None:
        count = np.cumsum(np.array(count))[:-1]
        ax.yaxis.set_major_locator(ticker.FixedLocator(count))
        ax.set_yticklabels(dseg_labels.keys())
    else:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    ax.xaxis.grid(False)
    ax.grid(linewidth=grid_width, axis='y', color=grid_color, linestyle=grid_style)

    # Add colorbar
    if colorbar:
        cbar = plt.colorbar(aximg, ax=ax, shrink=cbar_shrink)
        if colorbarlabel is not None:
            cbar.set_label(colorbarlabel)

    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if title is not None:
        ax.set_title(title)
    return None


def plot_overlay(base: Union[str,nib.Nifti1Image,Any], 
                 overlay: Optional[Union[str,nib.Nifti1Image,Any]] = None,
                 levels: List[int] = [0.5, 1.5, 2.5, 3.5], 
                 linewidth: int = 1, 
                 contrast: int = 0, 
                 colors: List[str] =['yellow', 'red', 'aqua', 'lime'], 
                 **kwargs
                ) -> None:
    """doc-string
    """
    display: plotting.plot_anat = plotting.plot_anat(base, dim=contrast * -1, **kwargs)
    if overlay is not None:
        display.add_contours(overlay, linewidths=linewidth, levels=levels, colors=colors)
    return None


def plot_overlay_contour(base: Union[str,nib.Nifti1Image,Any], 
                         overlay: Optional[Union[str,nib.Nifti1Image,Any]] = None,
                         levels: List[int] = [0.5, 1.5, 2.5, 3.5], 
                         linewidth: int = 1, 
                         contrast: int = 0, 
                         colors: List[str] = ['yellow', 'red', 'aqua', 'lime'], 
                         **kwargs
                        ) -> None:
    """doc-string
    """
    display = plotting.plot_anat(base, dim=contrast * -1, **kwargs)
    if overlay is not None:
        display.add_contours(overlay, linewidths=linewidth, levels=levels, colors=colors)
    return None


def plot_overlay_edges(base: Union[str,nib.Nifti1Image,Any], 
                       overlay: Optional[Union[str,nib.Nifti1Image,Any]] = None,
                       contrast: int = 0, 
                       color: Optional[str] = None, 
                       **kwargs
                      ) -> None:
    """doc-string
    """
    display = plotting.plot_anat(base, dim=contrast * -1, **kwargs)
    if overlay is not None:
        display.add_edges(overlay)
    return None


def to_file(fcn: callable, 
            fname: str, 
            *args, 
            **kwargs
           ) -> None:
    """doc-string.
    """
    dpi = kwargs.pop('dpi', 100)
    fcn(*args, **kwargs)
    fig = plt.gcf()
    # fig.tight_layout()
    plt.savefig(fname, dpi=dpi, bbox_inches='tight')
    # fig.savefig(fname, dpi=dpi)
    plt.close()
    return None


def to_base64(fcn: callable, 
              *args, 
              **kwargs
             ) -> str:
    """doc-string.
    """
    ext = kwargs.pop('ext', '.svg')
    with TemporaryFile(suffix=ext) as tmp:
        to_file(fcn, tmp, *args, **kwargs)
        tmp.seek(0)
        s: base64.b64encode = base64.b64encode(tmp.read()).decode("utf-8")
    return f'data:image/{ext};base64,{s}'
    # return 'data:image/{};base64,'.format(ext) + s
