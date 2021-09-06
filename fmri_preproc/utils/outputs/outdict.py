# -*- coding: utf-8 -*-
"""Output data template dictionaries for the ``fmri_preproc`` modules.
"""
import os
from abc import ABC

from typing import (
    Dict,
    KeysView
)

class OutDict(ABC):
    """Abstract base class doc-string.
    """
    __slots__ = [ "output" ]

    def __init__(self) -> None:
        """Abstract base class constructor.
        """
        self.output: Dict = {}
        return None
    
    def check_exists(self,
                     *args,
                    ) -> bool:
        """Check if input argument(s) exist. 
        If all of the input arguments exists, then True is returned, or False otherwise.
        If no input arguments are provided, then all dictionary keys are checked.
        """
        if len(args) == 0:
            args: KeysView = self.output.keys()
        
        for arg in args:
            file: str = self.output.get(arg,'')
            if not os.path.exists(file):
                return False
        return True
