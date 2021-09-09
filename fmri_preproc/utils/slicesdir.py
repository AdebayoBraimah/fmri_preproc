#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Generate slice-wise images for quality control purposes for the ``fmri_proc`` rs-fMRI preprocessing pipeline.
"""
from fmri_preproc.utils.qc import plot
from fmri_preproc.utils.qc import util
from fmri_preproc.utils.util import dict2json, json2dict
from nilearn import plotting as nilearn_plot, image as nilearn_img
import os
import jinja2
import tqdm

from fmri_preproc import HTMLDIR
from fmri_preproc.utils.workdir import WorkDir


def slicesdir(primary,
              secondary=None,
              tertiary=None,
              overlay=None,
              overlay_type='contour',
              mask=None,
              workdir='slicesdir',
              title='slicesdir.py',
              view=('axial', 'sagittal'),
              contrast=0,
              vmin=None,
              vmax=None,
              primary_label=None,
              secondary_label=None,
              tertiary_label=None):

    if secondary is not None and len(primary) != len(secondary):
        raise RuntimeError('number of secondary images must equal the number of primary images')

    if tertiary is not None and len(primary) != len(tertiary):
        raise RuntimeError('number of tertiary images must equal the number of primary images')

    if overlay is not None and len(primary) != len(overlay):
        raise RuntimeError('number of overlays must equal the number of primary images')

    if mask is not None and len(primary) != len(mask):
        raise RuntimeError('number of mask must equal the number of primary images')

    if not os.path.exists(workdir):
        os.makedirs(workdir)

    # env = jinja2.Environment(
    #     loader=jinja2.PackageLoader('dhcp', 'resources'),
    #     autoescape=jinja2.select_autoescape(['html', 'xml'])
    # )

    views_dict = {
        "sagittal": "x",
        "coronal": "y",
        "axial": "z",
    }

    view = view.split(',') if isinstance(view, str) else view
    # view = [views_dict[v] for v in view]

    # data_json = Path(workdir).join(f'data.json')
    # if data_json.exists():
    #     data = json2dict(data_json)
    # else:
    data = []

    pbar = tqdm.trange(len(primary))

    for idx in pbar:

        primary0 = primary[idx]

        # pbar.set_description(f"{Path(primary0).filename}")
        _, fname, _ = util.split(fname=primary0)
        pbar.set_description(f"{fname}")

        secondary0 = secondary if secondary is None else secondary[idx]
        tertiary0 = tertiary if tertiary is None else tertiary[idx]
        overlay0 = overlay if overlay is None else overlay[idx]
        mask0 = mask if mask is None else mask[idx]

        # print(idx, primary0, secondary0, overlay0)

        primary_label = 'Primary' if primary_label is None else primary_label
        secondary_label = 'Secondary' if secondary_label is None else secondary_label
        tertiary_label = 'Tertiary' if tertiary_label is None else tertiary_label

        data0 = {
            'primary_name': primary0,
            'secondary_name': secondary0,
            'tertiary_name': tertiary0,
            'primary_label': primary_label,
            'secondary_label': secondary_label,
            'tertiary_label': tertiary_label,
            'overlay_name': overlay0,
            'mask_name': mask0,
            'primary_img': [],
            'secondary_img': [],
            'tertiary_img': [],
        }

        # plot primary image

        if overlay_type == 'contour':
            plt_fcn = plot.plot_overlay_contour
        elif overlay_type == 'edges':
            plt_fcn = plot.plot_overlay_edges
        else:
            raise RuntimeError(f'Unknown overlay type: {overlay_type}')

        if mask0 is not None:
            primary0 = nilearn_img.math_img("img1 * img2", img1=primary0, img2=mask[idx])

            if overlay0 is not None:
                overlay0 = nilearn_img.math_img("img1 * img2", img1=overlay0, img2=mask[idx])

        for v in view:

            coords = nilearn_plot.find_cut_slices(primary0, direction=views_dict[v], n_cuts=9)

            # fname = Path(workdir).join(f'primary{idx}_{v}.png')

            with WorkDir(src=workdir) as wd:
                fname: str = wd.join(f'primary{idx}_{v}.png')

            plot.to_file(
                plt_fcn,
                fname,
                primary0,
                overlay=overlay0,
                annotate=False,
                black_bg=True,
                display_mode=views_dict[v],
                cut_coords=coords,
                colorbar=True,
                contrast=contrast,
                vmin=vmin,
                vmax=vmax,
                title=primary_label if secondary is not None else None,
            )

            data0['primary_img'] += [fname.replace(workdir, '.')]

        # plot secondary image

        if secondary0 is not None:

            for v in view:

                if mask0 is not None:
                    secondary0 = nilearn_img.math_img("img1 * img2", img1=secondary0, img2=mask[idx])

                # fname = Path(workdir).join(f'secondary{idx}_{v}.png')

                with WorkDir(src=workdir) as wd:
                    fname: str = wd.join(f'secondary{idx}_{v}.png')

                plot.to_file(
                    plt_fcn,
                    fname,
                    secondary0,
                    overlay=overlay0,
                    annotate=False,
                    black_bg=True,
                    display_mode=views_dict[v],
                    cut_coords=coords,
                    colorbar=True,
                    contrast=contrast,
                    vmin=vmin,
                    vmax=vmax,
                    title=secondary_label,
                )
                data0['secondary_img'] += [fname.replace(workdir, '.')]

        if tertiary0 is not None:

            for v in view:

                if mask0 is not None:
                    tertiary0 = nilearn_img.math_img("img1 * img2", img1=tertiary0, img2=mask[idx])

                # fname = Path(workdir).join(f'tertiary{idx}_{v}.png')

                with WorkDir(src=workdir) as wd:
                    fname: str = wd.join(f'tertiary{idx}_{v}.png')

                plot.to_file(
                    plt_fcn,
                    fname,
                    tertiary0,
                    overlay=overlay0,
                    annotate=False,
                    black_bg=True,
                    display_mode=views_dict[v],
                    cut_coords=coords,
                    colorbar=True,
                    contrast=contrast,
                    vmin=vmin,
                    vmax=vmax,
                    title=tertiary_label,
                )
                data0['tertiary_img'] += [fname.replace(workdir, '.')]


        data += [data0]

        with WorkDir(src=workdir) as wd:
            dict2json(data, wd.join('data.json'))

    parse(workdir, title=title)


def parse(workdir, title='slicesdir.py'):

    with WorkDir(src=workdir) as wd:
        data = json2dict(wd.join('data.json'))
    
    env = jinja2.Environment(
        loader=jinja2.PackageLoader(HTMLDIR),
        autoescape=jinja2.select_autoescape(['html', 'xml'])
    )

    template = env.get_template('slicesdir.html')
    html = template.render(
        title=title,
        workdir=os.path.abspath(os.path.expanduser(workdir)),
        data=data,
    )

    with open(os.path.join(workdir, 'index.html'), 'w') as outfile:
        outfile.write(html)

