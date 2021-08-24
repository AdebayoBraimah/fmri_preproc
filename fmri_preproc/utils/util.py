# -*- coding: utf-8 -*-
"""Utility module for ``fmri_preproc`` resting-state fMRI pre-processing pipeline.
"""
import os
import subprocess
import shutil
from time import time 

from typing import (
    Any,
    Dict,
    List,
    Optional,
    Tuple,
    Union
)

from fmri_preproc.utils.logutil import LogFile

def timeops(log: Optional[LogFile] = None) -> Any:
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
    def decor(func: Any) -> Any:
        """Inner decorated function that accepts functions."""
        def timed(*args,**kwargs) -> Any:
            """Nested decorator function the performs timing of an operation.
            """
            start: float = time()
            if log: log.log(f"BEGIN {func.__name__}", use_header=True)
            result: Any = func(*args,**kwargs)
            end: float = time()
            if log: log.log(f"END: {func.__name__}  |  Time elapsed: {(end - start):2f} sec.", use_header=True)
            return result
        return timed
    return decor
