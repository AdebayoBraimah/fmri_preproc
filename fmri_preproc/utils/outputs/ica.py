# -*- coding: utf-8 -*-
"""Output data filepaths dictionaries for the ``ica`` module.

.. autosummary::
    :nosignatures:

    ICAFiles

.. autoclass:: ICAFiles
    :members:
"""
import os

from typing import Dict, Union

from fmri_preproc.utils.outputs.outdict import OutDict
from fmri_preproc.utils.workdir import WorkDir


class ICAFiles(OutDict):
    """class doc-string.
    """

    def __init__(self, outdir: Union[WorkDir, str]) -> None:
        """Class constructor method.
        NOTE: ``outdir`` should be subject working directory.
        """
        self.outdir: Union[WorkDir, str] = outdir
        super(ICAFiles, self).__init__()

    def outputs(self) -> Dict[str, str]:
        """doc-string.
        """
        if isinstance(self.outdir, WorkDir):
            outdir: str = self.outdir.src
        else:
            outdir: str = self.outdir

        icadir: str = os.path.join(outdir, "ica")

        if os.path.exists(icadir):
            remove_dir: bool = False
        else:
            remove_dir: bool = True

        with WorkDir(src=icadir) as od:
            self.output: Dict[str, str] = {
                "icadir": od.abspath(),
                "func_filt": od.join("func_filtered.nii.gz"),
                "meldir": od.join("func_filtered.ica"),
            }

            if remove_dir:
                od.rmdir()

        return self.output
