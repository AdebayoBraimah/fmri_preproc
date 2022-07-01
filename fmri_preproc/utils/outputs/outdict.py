# -*- coding: utf-8 -*-
"""Output data template dictionaries for the ``fmri_preproc`` modules.

.. autosummary::
    :nosignatures:

    OutDict

.. autoclass:: OutDict
    :members:
"""
import os
from abc import ABC

from typing import Any, Dict, KeysView


class OutDict(ABC):
    """Abstract base class doc-string.
    """

    __slots__ = "output"

    def __init__(self) -> None:
        """Abstract base class constructor.
        """
        self.output: Dict = {}
        return None

    def check_exists(self, *args,) -> bool:
        """Check if input argument(s) exist. 
        If all of the input arguments exists, then True is returned, or False otherwise.
        If no input arguments are provided, then all dictionary keys are checked.
        """
        if len(args) == 0:
            args: KeysView = self.output.keys()

        for arg in args:
            file: str = self.output.get(arg, "")
            if not os.path.exists(file):
                return False
        return True

    def set_key_to_none(self) -> Dict:
        """Sets values of keys to ``None`` IF the value is a 
        file or directory that does not exist.
        """
        return key_to_none(d=self.output)


def key_to_none(d: Dict[Any, Any]) -> Dict[Any, Any]:
    """Iterates through some input dictionary's keys, determines if the associated
    key-mapped value is a directory or file. If the value is a directory or file that
    DOES NOT EXIST, then it is set to ``None``.
    """
    # Only looking for files that do not exist.
    # Leave other non-filepath dictionary values as is.
    for key, val in d.items():
        if (
            (key == "subid")
            or (key == "sesid")
            or (key == "scan_pma")
            or (key == "birth_ga")
        ):
            continue
        try:
            val: str = os.path.abspath(val)
            if not os.path.exists(val):
                d[key] = None
        except TypeError:
            pass
    return d
