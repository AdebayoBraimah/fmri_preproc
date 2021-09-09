#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Quality control for the ``fmri_proc`` rs-fMRI preprocessing pipeline.
"""
import os
import nibabel as nb
import numpy as np
import jinja2
import pandas as pd
from collections import OrderedDict
import nilearn.plotting as nilearn_plot

import fmri_preproc.utils.qc.util as util
import fmri_preproc.utils.qc.metrics as metrics
import fmri_preproc.utils.qc.plot as plot
import warnings

from fmri_preproc.utils.qc.util import load_img

from tempfile import TemporaryDirectory

class Subject(object):
    """class doc-string
    """

    def __init__(self, workdir: str):
        """class constructor.
        """
        self.properties = OrderedDict()
        self.workdir = workdir

        if not os.path.exists(workdir):
            os.makedirs(workdir)

        # self._json = os.path.join(workdir, 'qc.json')
        self._json = 'qc.json'


    @property
    def workdir(self):
        return self._workdir

    @workdir.setter
    def workdir(self, d):
        if not os.path.exists(d):
            os.makedirs(d)
        self._workdir = d

    def parse(self, template: str, outname: str, group_json: str=None):

        env = jinja2.Environment(
            loader=jinja2.PackageLoader('dhcp', 'resources'),
            autoescape=jinja2.select_autoescape(['html', 'xml']),
            extensions=['jinja2.ext.do']
        )

        args = dict(
            ind_obj=self,
            numpy=np,
            pandas=pd,
            util=util,
            plot=plot,
            nl_plot=nilearn_plot
        )

        if group_json is not None:
            group_json = Group(group_json)
            args['grp_obj'] = group_json

        template = env.get_template(template)
        html = template.render(**args)

        with open(outname, 'w') as outfile:
            outfile.write(html)

    def has(self, key):
        return key in self.properties.keys()

    def add(self, key, val):
        self.add_property(key, val)

    def add_property(self, key, val):

        if self.has(key):
            warnings.warn(f'key ({key}) is not unique.  Will overwrite previous value.')

        if isinstance(val, str) and os.path.isfile(val):
            val = os.path.relpath(val, self.workdir)
        elif isinstance(val, nb.Nifti1Image):
            val = os.path.relpath(val.get_filename(), self.workdir)

        self.properties[key] = val

    def add_property_dict(self, dict, prefix=None):
        for k, v in dict.items():
            if prefix is not None:
                k = prefix+'_'+k
            self.add_property(k, v)

    def get(self, key):
        val = self.properties[key]

        if isinstance(val, str) and os.path.isfile(os.path.join(self.workdir, val)):
            val = os.path.join(self.workdir, val)

        if isinstance(val, str) and (val.endswith('.nii.gz') or val.endswith('.nii')):
            val = nb.load(val)

        return val

    def get_filename(self, key):
        fname = os.path.join(self.workdir, self.properties[key])
        assert os.path.isfile(fname), 'Invalid file: {}'.format(fname)
        return fname

    def get_image(self, key):
        fname = os.path.join(self.workdir, self.properties[key])
        assert os.path.isfile(fname), 'Invalid file: {}'.format(fname)
        return nb.load(fname)

    def get_nparray(self, key):
        return np.array(self.properties[key])

    def to_dataframe(self):
        return pd.DataFrame(pd.Series(self.properties)).T

    def to_json(self, outname):

        if not outname.endswith('.json'):
            outname = outname + '.json'
        # self.properties['key_order'] = list(self.properties.keys())
        util.dict2json(self.properties, os.path.join(self.workdir, outname))

    @classmethod
    def from_json(clss, jsonfile):
        js = util.json2dict(jsonfile)
        workdir = os.path.dirname(jsonfile)
        d = clss(workdir=workdir)
        d.add_property_dict(js)
        d.fname = jsonfile

        return d

    @classmethod
    def from_qcdir(clss, workdir):
        d = clss(workdir=workdir)

        jsonfile = os.path.join(workdir, 'qc.json')
        if os.path.exists(jsonfile):
            # js = pd.read_json(jsonfile, typ='series')
            js = util.json2dict(jsonfile)
            d.add_property_dict(js)

        return d

    def add_subjectinfo(self, subject_id=None, session_id=None, scan_age=None, birth_age=None):
        self.add('subid', subject_id)
        self.add('sesid', session_id)
        self.add('scan_age', scan_age)
        self.add('birth_age', birth_age)
        self.to_json(self._json)

    def add_reg(self, label, source, ref, ref_brainmask, ref_boundarymask=None, ref_dseg=None):

        # assert not self.has(label + '_source_fname'), "Label {} is not unique".format(label)

        source = load_img(source)
        ref = load_img(ref)
        ref_brainmask = load_img(ref_brainmask)
        ref_boundarymask = load_img(ref_boundarymask)
        ref_dseg = load_img(ref_dseg)

        self.add(label + '_source_fname', source.get_filename())
        self.add(label + '_ref_fname', ref.get_filename())
        self.add(label + '_qctype', 'reg')

        normmi = metrics.measurecost(source, ref, cost='normmi', ref_brainmask=ref_brainmask)
        normcorr = metrics.measurecost(source, ref, cost='normcorr', ref_brainmask=ref_brainmask)

        self.add(label + '_normmi', normmi)
        self.add(label + '_normcorr', normcorr)

        if ref_boundarymask is not None:
            self.add(label + '_ref_boundary_fname', ref_boundarymask.get_filename())
            bbr = metrics.measurecost(
                source, ref, cost='bbr', ref_brainmask=ref_brainmask, boundarymask=ref_boundarymask
            )
            self.add(label + '_bbr', bbr)

        if ref_dseg is not None:
            self.add(label + '_ref_dseg_fname', ref_dseg.get_filename())

        # --- write to json ---

        self.to_json(self._json)

    def add_motparam(self, label, motparams=None, func=None):

        # assert any([motparams, func]), 'Either motparams or func arg must be specified'

        # check if key/label already exists
        # assert label not in self.properties['labels_motp_qc'], "Label {} is not unique".format(label)
        # assert not self.has(label + '_fname'), "Label {} is not unique".format(label)
        # self.properties['labels_motp_qc'].append(label)

        # --- FILE INFO ---
        fname = motparams if motparams is not None else func
        self.add(label + '_fname', fname)
        self.add(label + '_qctype', 'motparam')

        # --- MOTION STATS ---

        if motparams is not None:
            motparams = pd.read_csv(motparams, delimiter='\t', index_col=None)
            # motparams = np.loadtxt(motparams)

        mp_stats = metrics.motparams(func, motparams, tmpdir=self.workdir)
        self.add_property_dict(mp_stats, prefix=label)

        # framewise displacement
        #
        # fd = metrics.fd(mp_stats['mp'])
        # self.add_property_dict(fd, prefix=label)

        # --- write to json ---

        self.to_json(self._json)

    def add_func(self,
                 func,
                 label,
                 brainmask=None,
                 standard=None,
                 func2standard_warp=None,
                 template=None,
                 template2func_warp=None,
                 template_dseg=None,
                 template_dseg_labels=None):

        # check if key/label already exists
        # assert label not in self.properties['labels_func_qc'], "Label {} is not unique".format(label)
        # assert not self.has(label + '_fname'), "Label {} is not unique".format(label)
        # self.properties['labels_func_qc'].append(label)

        func = load_img(func)
        brainmask = load_img(brainmask)
        standard = load_img(standard)
        func2standard_warp = load_img(func2standard_warp)
        template = load_img(template)
        template2func_warp = load_img(template2func_warp)
        template_dseg = load_img(template_dseg)

        # TODO: add template_dseg and template_dseg_labels

        basename = os.path.join(self.workdir, label)

        # create brainmask

        if brainmask is None:
            brainmask = basename + '_brainmask.nii.gz'
            with TemporaryDirectory(dir=self.workdir) as td:
                util.run(['fslmaths', func.get_filename(), '-Tmean', os.path.join(td, 'mean')])
                util.run([
                    'bet', os.path.join(td, 'mean'), os.path.join(td, 'bet'), '-R',
                    '-n', '-m'])
                util.run(['fslmaths', os.path.join(td, 'bet_mask'), '-dilM', brainmask])
            brainmask = nb.load(brainmask)
            self.add(label + '_brainmask', brainmask)
        else:
            self.add(label + '_brainmask', brainmask)

        # transform template_dseg to func space

        if template_dseg is not None and template2func_warp is not None:
            dseg = basename + '_dseg.nii.gz'
            util.run([
                'applywarp', '-i', template_dseg.get_filename(),
                '-r', brainmask.get_filename(),
                '-w', template2func_warp,
                '-o', dseg,
                '--interp=nn'
            ])
            dseg = nb.load(dseg)
            self.add(label + '_func_dseg', dseg)
            self.add(label + '_func_dseg_labels', template_dseg_labels)

        # --- FILE INFO ---

        self.add(label + '_fname', func)
        self.add(label + '_qctype', 'func')

        # --- SNR ---

        # TODO: fold snr and snr_stats into one dictionary
        snr, snr_stats = metrics.tSNR(
            func=func,
            func_brainmask=brainmask,
            basename=basename,
            standard=standard,
            func2standard_warp=func2standard_warp,
            tmpdir=self.workdir)
        self.add_property_dict({**snr, **snr_stats}, prefix=label)

        # --- DUAL REGRESSION, CNR, SPATIAL_CORR ---

        if template is not None and template2func_warp is not None:
            # TODO: fold dr and dr_stats into one dictionary
            # TODO: do this in template space (rather than native)
            # TODO: add spatial and netmat group similarity
            dr, dr_stats = metrics.dr(
                func=func,
                func_brainmask=brainmask,
                spatial_template=template,
                template2func_warp=template2func_warp,
                basename=basename,
                func2standard_warp=func2standard_warp,
                standard=standard,
                tmpdir=self.workdir)
            self.add_property_dict({**dr, **dr_stats}, prefix=label)

        # --- DVARS ---

        dvars = metrics.dvars(func, brainmask)
        self.add_property_dict(dvars, prefix=label)

        func.uncache()

        # --- write to json ---

        self.to_json(self._json)

    def add_fmap(self, label, fmap, fmap_mag, fmap_brainmask, spinecho=None):

        # assert not self.has(f'{label}_ph_fname'), f"Label {label} is not unique"

        fmap = load_img(fmap)
        fmap_mag = load_img(fmap_mag)
        spinecho = load_img(spinecho)
        fmap_brainmask = load_img(fmap_brainmask)

        self.add(label + '_ph_fname', fmap.get_filename())
        self.add(label + '_mag_fname', fmap_mag.get_filename())
        self.add(label + '_brainmask_fname', fmap_brainmask.get_filename())
        self.add(label + '_qctype', 'fmap')

        if spinecho is not None:
            self.add(label + '_spinecho_fname', spinecho.get_filename())
            zsmth = metrics.z_smoothness(spinecho.get_filename())
            self.add(f'{label}_spinecho_zsmoothness', zsmth)
            # TODO: add min zsmoothness per pedir

        # --- write to json ---

        self.to_json(self._json)

    def add_fix(self, label, ic, mix, labels=None):

        # assert not self.has(f'{label}_ic_fname'), f"Label {label} is not unique"

        ic = load_img(ic)

        self.add(label + '_ic_fname', ic.get_filename())
        self.add(label + '_mix_fname', mix)

        if labels is not None:
            self.add(label + '_labels_fname', labels)

            with open(labels, 'r') as f:
                tmp = f.readlines()

            unknown = 0
            signal = 0
            noise = 0

            for l in tmp[1:-1]:
                if 'unknown' in l.lower():
                    unknown += 1
                if 'signal' in l.lower():
                    signal += 1
                if 'noise' in l.lower():
                    noise += 1

            self.add(label + '_unknown_N', unknown)
            self.add(label + '_signal_N', signal)
            self.add(label + '_noise_N', noise)

        # --- write to json ---

        self.to_json(self._json)

    def __str__(self):
        return "[workdir: " + self.workdir + "] " + str(self.properties)


class Group(object):
    """class doc-string
    """

    def __init__(self, json, unique_key=None, drop_na=False):
        self.fname = json
        self.group_scores = pd.read_json(json, typ='frame', orient='split')
        
        # This FIX is for compatibility with old qc files that contained the key 'age' instead of 'scan_age'
        self.group_scores = self.group_scores.rename(columns={'age': 'scan_age'})
        
        if unique_key is not None:
            self.group_scores.set_index(['subid', 'sesid'])
        if drop_na:
            self.group_scores = self.group_scores.dropna(axis=0)
        # self.group_scores = self.group_scores.reindex(columns=key_order)

    def parse(self, template: str, outname: str):

        env = jinja2.Environment(
            loader=jinja2.PackageLoader('report', 'templates'),
            autoescape=jinja2.select_autoescape(['html', 'xml'])
        )

        template = env.get_template(template)
        html = template.render(
            group_obj=self,
            numpy=np,
            util=util,
            plot=plot,
            metrics=metrics
        )

        with open(os.path.expanduser(outname), 'w') as outfile:
            outfile.write(html)

    def get_images(self, key):
        fnames = self.group_scores[key].values.tolist()
        wdirs = self.group_scores['workdir'].values.tolist()
        img = []
        for d, f in zip(wdirs, fnames):
            if f is not None and d is not None:
                img += [nb.load(os.path.join(d, f))]
        return img

    @classmethod
    def from_json(clss, outfile, unique_key, *args):
        json = []
        for f in args:
            f = os.path.expanduser(f)
            tmp = util.json2dict(f)
            tmp['workdir'] = os.path.dirname(f)
            key_order = list(tmp.keys())
            json += [pd.Series(tmp)]
        json = pd.concat(json, axis=1)

        grp_scores = pd.DataFrame(json.T)
        grp_scores = grp_scores.reindex(columns=key_order)

        if unique_key is not None:
            grp_scores = grp_scores.set_index(unique_key)

        grp_scores.reset_index().to_json(os.path.expanduser(outfile), orient='split')

        d = clss(outfile, unique_key=unique_key)
        return d

