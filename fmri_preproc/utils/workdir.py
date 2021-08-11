# -*- coding: utf-8 -*-
"""Working directory module.
"""
import os
import shutil

class WorkDir:
    """Working directory base class that instantiates ``WorkDir`` objects that creates and manipulates working directories.

    Attributes:
        work_dir: Input working directory.
        parent_dir: Parent directory.

    Usage example:
            >>> # Using class object as context manager
            >>> ## Create work directory , then clean-up (remove it)
            >>> with WorkDir(work_dir="/path/to/working_directory", use_cwd=False) as work:
            ...     work.mkdir()
            ...     work.rmdir(rm_parent=False)
            ...
            >>> # or
            >>>
            >>> work = TmpDir(work_dir="/path/to/working_directory", 
            ...               use_cwd=False)
            >>> work.mkdir()
            >>> work
            "/path/to/working_directory"
            >>> work.rmdir(rm_parent=False)

    Arguments:
        work_dir: Working directory name/path. This directory need not exist at runtime.
        use_cwd: Use current working directory as the parent directory.
    """
    __slots__ = [ 
                    "work_dir", 
                    "parent_dir"
                ]

    def __init__(self, 
                 work_dir: str,
                 use_cwd: bool = False
                ) -> None:
        """Initialization method for the ``WorkDir`` base class.

        Usage example:
            >>> # Using class object as context manager
            >>> ## Create work directory , then clean-up (remove it)
            >>> with WorkDir(work_dir="/path/to/working_directory", use_cwd=False) as work:
            ...     work.mkdir()
            ...     work.rmdir(rm_parent=False)
            ...
            >>> # or
            >>>
            >>> work = TmpDir(work_dir="/path/to/working_directory", 
            ...               use_cwd=False)
            >>> work.mkdir()
            >>> work
            "/path/to/working_directory"
            >>> work.rmdir(rm_parent=False)
        
        Arguments:
            work_dir: Working directory name/path. This directory need not exist at runtime.
            use_cwd: Use current working directory as the parent directory.
        """
        self.work_dir: str = work_dir
        self.parent_dir: str = os.path.dirname(self.work_dir)

        if use_cwd:
            _cwd: str = os.getcwd()
            self.work_dir: str = os.path.join(_cwd,self.work_dir)
            self.parent_dir: str = os.path.dirname(self.work_dir)

    def __enter__(self):
        """Context manager entrance method."""
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        """Context manager exit method."""
        return False
    
    def __repr__(self):
        """Representation request method."""
        return self.work_dir
    
    def mkdir(self) -> None:
        """Makes/creates the working directory.

        This class method is analogous to UNIX's ``mkdir -p`` command and option combination.

        Usage example:
            >>> # Using class object as context manager
            >>> with WorkDir(work_dir="/path/to/working_directory", use_cwd=False) as work:
            ...     work.mkdir()
            ...
            >>> # or
            >>>
            >>> work = TmpDir(work_dir="/path/to/working_directory", 
            ...               use_cwd=False)
            >>> work.mkdir()
            >>> work
            "/path/to/working_directory"
        """
        if not os.path.exists(self.work_dir):
            return os.makedirs(self.work_dir)
        else:
            print("Working directory already exists.")
            return None
    
    def rmdir(self, 
              rm_parent: bool = False
             ) -> None:
        """Removes working directory, and the parent directory if indicated to do so.

        This class method is analogous to UNIX's ``rm -rf`` command and option combination.

        Usage example:
            >>> # Using class object as context manager
            >>> with WorkDir(work_dir="/path/to/working_directory", use_cwd=False) as work:
            ...     work.mkdir()
            ...     work.rmdir(rm_parent=False)
            ...
            >>> # or
            >>>
            >>> work = TmpDir(work_dir="/path/to/working_directory", 
            ...               use_cwd=False)
            >>> work.mkdir()
            >>> work.rmdir(rm_parent=False)

        Arguments:
            rm_parent: Removes parent directory as well.
        """
        if rm_parent and os.path.exists(self.parent_dir):
            return shutil.rmtree(self.parent_dir,ignore_errors=True)
        elif os.path.exists(self.work_dir):
            return shutil.rmtree(self.work_dir,ignore_errors=True)
        else:
            print("Working directory does not exist.")
            return None