# -*- coding: utf-8 -*-
"""Utility module for ``fmri_preproc`` resting-state fMRI pre-processing pipeline.
"""
import os
import json
import tarfile
import urllib.request as urllib

from time import time
from tqdm import tqdm
from tempfile import TemporaryDirectory

from typing import (
    Any,
    Dict,
    Optional,
)

from fmri_preproc.utils.fileio import File
from fmri_preproc.utils.logutil import LogFile
from fmri_preproc.utils.tempdir import TmpDir

from fmri_preproc import (
    ATLASDIR,
    GROUP_MAP_DIR,
    GROUP_QC_DIR
)


# Globlally define (temporary) log file object
with TmpDir(src=os.getcwd()) as tmpd:
    with TmpDir.TmpFile(tmp_dir=tmpd.src, ext='.log') as tmpf:
        log: LogFile = LogFile(log_file=tmpf.src)


def timeops(log: Optional[LogFile] = None) -> callable:
    """Decorator function that times some operation and writes that time to
    a log file object.

    Usage example:
        >>> from fmri_preproc.utils.logutil import LogFile
        >>> log = LogFile('my_log_file.log')
        >>>
        >>> @timeops(log)
        >>> def my_func(args*, log):
        ...     for i in args:
        ...         log.log(f"This is an arg: {i}")
        ...     return None
        ...
        >>> # The length of time to complete the operation 
        >>> # should be written to the log file.
        >>> myfunc(args*, log)  

    Arguments:
        log: Log file object to be written to.
    """
    def decor(func: callable) -> callable:
        """Inner decorated function that accepts functions."""
        def timed(*args,**kwargs) -> callable:
            """Nested decorator function the performs timing of an operation.
            """
            start: float = time()
            if log: log.log(f"BEGIN: {func.__name__}", use_header=True)
            result: callable = func(*args,**kwargs)
            end: float = time()
            if log: log.log(f"END: {func.__name__}  |  Time elapsed: {(end - start):2f} sec.", use_header=True)
            return result
        return timed
    return decor


def json2dict(jsonfile: str) -> Dict[Any,Any]:
    """Read JSON file to dictionary.
    """
    with open(jsonfile, 'r') as file:
        d: Dict[Any,Any] = json.load(file)
    return d


def dict2json(dict: Dict[Any,Any],
              jsonfile: str,
              indent: int = 4
             ) -> str:
    """Write dictionary to JSON file.
    """
    with open(jsonfile, 'w') as out:
        json.dump(dict, out, indent=indent)
    return out


def update_sidecar(file: str, **kwargs) -> str:
    """Updates a JSON sidecar/file.
    """
    with File(src=file, assert_exists=False) as f:
        dirname, basename, _ = f.file_parts()
        jsonfile: str = os.path.join(dirname,basename + '.json')
        with File(src=jsonfile) as jf:
            jsonfile: str = jf.abspath()
    
    d: Dict[Any,Any] = load_sidecar(file=jsonfile)
    d.update(**kwargs)
    jsonfile: str = dict2json(dict=d, jsonfile=jsonfile, indent=4)

    return jsonfile


def load_sidecar(file: str) -> Dict[Any,Any]:
    """Reads in a JSON sidecar/file.
    """
    d: Dict[Any,Any] = {}
    with File(src=file, assert_exists=False) as f:
        dirname, basename, _ = f.file_parts()
        jsonfile: str = os.path.join(dirname,basename + '.json')
        with File(src=jsonfile) as jf:
            if os.path.exists(jf.abspath()):
                with open(jf.abspath(),'r') as j:
                    d.update(json.load(j))
    return d


def get_fsl_version() -> str:
    """Returns a string that represents the version of ``FSL`` in the system path.
    """
    fsl_version_file: str = os.path.join(os.environ['FSLDIR'], 'etc/fslversion')
    with open(fsl_version_file, 'r') as file:
        ver: str = file.read().split(':')[0]
    return ver


def fetch_dhcp_volumetric_atlas(path: Optional[str] = None, 
                                extended: bool = True
                               ) -> None:
    """doc-string
    """
    path: str = ATLASDIR if path is None else os.path.realpath(os.path.expanduser(path))

    if extended:
        atlas_url: str = 'https://users.fmrib.ox.ac.uk/~seanf/dhcp-augmented-volumetric-atlas-extended.tar.gz'
        path: str = os.path.join(path,'dhcp_volumetric_atlas_extended')
    else:
        atlas_url: str = 'https://users.fmrib.ox.ac.uk/~seanf/dhcp-augmented-volumetric-atlas.tar.gz'
        path: str = os.path.join(path,'dhcp_volumetric_atlas')

    if not os.path.exists(path):
        os.makedirs(path)
    
    print('Download atlas:')

    with TemporaryDirectory(dir=path) as tmp:

        with DownloadBar(unit='B', unit_scale=True, miniters=1,
                         desc=atlas_url.split('/')[-1]) as t:  # all optional kwargs
            file_tmp = urllib.urlretrieve(atlas_url, filename=os.path.join(tmp, 'atlas.tar.gz'), reporthook=t.update_to)[0]
        
        print('Unpack atlas:')
        with tarfile.open(name=file_tmp) as tar:
            for m in tqdm(iterable=tar.getmembers(), total=len(tar.getmembers())):
                tar.extract(member=m, path=path)
        
        unc_data: str = os.path.join(ATLASDIR,'UNC.tar.gz')
        path: str = os.path.join(ATLASDIR,'UNCInfantAtlas_2020')
        if os.path.exists(unc_data):
            print('Unpack UNC atlas:')
            with tarfile.open(name=unc_data) as tar:
                for m in tqdm(iterable=tar.getmembers(), total=len(tar.getmembers())):
                    tar.extract(member=m, path=path)
    return None


def fetch_dhcp_group_qc(path: Optional[str] = None) -> None:
    rdata_url: str = 'https://users.fmrib.ox.ac.uk/~seanf/grp_qc_512_anon.json'

    path = GROUP_QC_DIR if path is None else os.path.realpath(os.path.expanduser(path))

    if not os.path.exists(path):
        os.makedirs(path)

    path: str = os.path.join(path, os.path.basename(rdata_url))

    print('Download group QC:')

    with DownloadBar(unit='B', unit_scale=True, miniters=1,
                     desc=rdata_url.split('/')[-1]) as t:  # all optional kwargs
        urllib.urlretrieve(rdata_url, filename=path, reporthook=t.update_to)
    return None


def fetch_dhcp_group_maps(path: Optional[str] = None) -> None:
    rdata_url = 'https://users.fmrib.ox.ac.uk/~seanf/group_maps.nii.gz'

    path = GROUP_MAP_DIR if path is None else os.path.realpath(os.path.expanduser(path))

    if not os.path.exists(path):
        os.makedirs(path)

    path: str = os.path.join(path, os.path.basename(rdata_url))

    print('Download group maps:')

    with DownloadBar(unit='B', unit_scale=True, miniters=1,
                     desc=rdata_url.split('/')[-1]) as t:  # all optional kwargs
        urllib.urlretrieve(rdata_url, filename=path, reporthook=t.update_to)
    return None



# FROM TQDM DOCS https://github.com/tqdm/tqdm#table-of-contents
class DownloadBar(tqdm):
    """Provides ``update_to(n)`` which uses ``tqdm.update(delta_n)``."""

    def update_to(self, b=1, bsize=1, tsize=None):
        """
        b  : int, optional
            Number of blocks transferred so far [default: 1].
        bsize  : int, optional
            Size of each block (in tqdm units) [default: 1].
        tsize  : int, optional
            Total size (in tqdm units). If [default: None] remains unchanged.
        """
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)  # will also set self.n = b * bsize
