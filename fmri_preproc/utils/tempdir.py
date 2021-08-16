# -*- coding: utf-8 -*-
"""Temporary working directory and file module.
"""
import os
import random

from typing import Optional

from fmri_preproc.utils.workdir import WorkDir
from fmri_preproc.utils.io import File

class TmpDir(WorkDir):
    """Temporary directory class that creates (random) temporary directories and files given a parent directory. This class inherits methods from the ``WorkDir`` base class.
    
    Attributes:
        work_dir: Input temproary working directory.
        parent_dir: Parent directory of the specified temproary directory.
    
    Usage example:
            >>> with TmpDir("/path/to/temporary_directory",False) as tmp_dir:
            ...     tmp_dir.mk_tmp_dir()
            ...     # do more stuff
            ...     tmp_dir.rm_tmp_dir(rm_parent=False)
            ...
            >>> # or
            >>>
            >>> tmp_dir = TmpDir("/path/to/temporary_directory")
            >>> tmp_dir
            "/path/to/temporary_directory"
            >>> tmp_dir.rm_tmp_dir(rm_parent=False)
        
    Arguments:
        tmp_dir: Temporary parent directory name/path.
        use_cwd: Use current working directory as working direcory.
    """
    
    def __init__(self,
                 tmp_dir: str,
                 use_cwd: bool = False
                ) -> None:
        """Initialization method for the TmpDir child class.
        
        Usage example:
            >>> with TmpDir("/path/to/temporary_directory",False) as tmp_dir:
            ...     tmp_dir.mk_tmp_dir()
            ...     # do more stuff
            ...     tmp_dir.rm_tmp_dir(rm_parent=False)
            ...
            >>> # or
            >>> tmp_dir = TmpDir("/path/to/temporary_directory")
            >>> tmp_dir
            "/path/to/temporary_directory"
            >>> tmp_dir.rm_tmp_dir(rm_parent=False)
        
        Arguments:
            tmp_dir: Temporary parent directory name/path.
            use_cwd: Use current working directory as working direcory.
        """
        _n: int = 10000
        tmp_dir: str = os.path.join(tmp_dir,'tmp_dir_' + 
                                    str(random.randint(0,_n)))
        
        if use_cwd:
            _cwd: str = os.getcwd()
            tmp_dir = os.path.join(_cwd,tmp_dir)
        
        super(TmpDir, self).__init__(tmp_dir,
                                     use_cwd)
    
    class TmpFile(File):
        """Sub-class of ``TmpDir`` class, which creates and manipulates temporary files via inheritance from the ``File`` object base class.
        
        Attributes:
            file: Temporary file name.
            ext: File extension of input file. If no extension is provided, it is inferred.
        
        Usage example:
            >>> tmp_directory = TmpDir("/path/to/temporary_directory")
            >>>
            >>> temp_file = TmpDir.TmpFile(tmp_directory.tmp_dir,
            ...                             ext="txt")
            ...
            >>> temp_file
            "/path/to/temporary_directory/temporary_file.txt"
        
        Arguments:
            tmp_dir: Temporary directory name.
            tmp_file: Temporary file name.
            ext: Temporary directory file extension.
        """

        def __init__(self,
                     tmp_dir: str,
                     tmp_file: Optional[str] = "",
                     ext: Optional[str] = "",
                    ) -> None:
            """Initialization method for the TmpFile sub-class that inherits from the ``File`` base class, allowing for the creation and maninuplation of temporary files.
            
            Usage example:
                >>> tmp_directory = TmpDir("/path/to/temporary_directory")
                >>>
                >>> temp_file = TmpDir.TmpFile(tmp_directory.tmp_dir,
                ...                             ext="txt")
                ...
                >>> temp_file
                "/path/to/temporary_directory/temporary_file.txt"
            
            Arguments:
                tmp_dir: Temporary directory name.
                tmp_file: Temporary file name.
                ext: File extension.
            """          
            if tmp_file:
                pass
            else:
                _n: int = 10000 
                tmp_file: str = "tmp_file_" + str(random.randint(0,_n))
            
            if ext:
                tmp_file: str = tmp_file + f".{ext}"

            tmp_file: str = os.path.join(tmp_dir, 
                                         tmp_file)

            super(TmpDir.TmpFile, self).__init__(tmp_file, 
                                                 ext)