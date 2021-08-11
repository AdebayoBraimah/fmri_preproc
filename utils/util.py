# -*- coding: utf-8 -*-
"""Utility module for classes and functions that pertain to wrapping command line executables.
"""
import os
import subprocess
import shutil

from typing import (
    Dict,
    List,
    Optional,
    Tuple,
    Union
)

from fileio import File
from logutil import LogFile

# import subprocess
# import logging
# import os
# import random
# import shutil
# import nibabel as nib
# 
# from enums import NiiHeaderField
# 
# from typing import (
#     Dict, 
#     List, 
#     Optional, 
#     Tuple,
#     Union
# )

class DependencyError(Exception):
    """Exception intended for unment dependencies"""
    pass

# class InvalidNiftiFileError(Exception):
#     """Exception intended for invalid NIFTI files."""
#     pass
# 
# class NiftiFileIOWarning(Warning):
#     """Warning that is raised in the case of char byte overflow written to NIFTI file headers."""
#     pass

# class File:
#     """File object base class. This class creates a ``File`` object that encapsulates a number of methods and properites for file and filename handling, and file manipulation.
#     
#     Attributes:
#         file: Input file.
#         ext: File extension of input file. If no extension is provided, it is inferred.
# 
#     Usage example:
#         >>> # Using class object as context manager
#         >>> with File("file_name.txt") as file:
#         ...     file.touch()
#         ...     file.write_txt("some text")
#         ...
#         >>> # or
#         >>> 
#         >>> file = File("file_name.txt")
#         >>> file
#         "file_name.txt"
# 
#     Arguments:
#         file: Input file (need not exist at runtime/instantiated).
#         ext: File extension of input file. If no extension is provided, it is inferred.
#     """
#     __slots__ = [ 
#                     "file", 
#                     "ext"
#                 ]
#     
#     def __init__(self,
#                  file: str,
#                  ext: Optional[str] = "",
#                  assert_exists: bool = False
#                 ) -> None:
#         """Initialization method for the File base class.
#         
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with File("file_name.txt") as file:
#             ...     file.touch()
#             ...     file.write_txt("some text")
#             ...
#             >>> # or
#             >>> 
#             >>> file = File("file_name.txt")
#             >>> file
#             "file_name.txt"
# 
#         Arguments:
#             file: Input file (need not exist at runtime/instantiated).
#             ext: File extension of input file. If no extension is provided, it is inferred.
#             assert_exists: Asserts that the specified input file must exist. 
#         """
#         self.file: str = file
# 
#         # if ext:
#         #     self.ext: str = ext
#         # elif self.file.endswith('.gz'):
#         #     self.ext: str = self.file[-(7):]
#         # else:
#         #     self.ext: str = self.file[-(4):]
#         
#         if ext:
#             self.ext: str = ext
#         else:
#             _, self.ext = os.path.splitext(self.file)
#         
#         if assert_exists:
#             assert os.path.exists(self.file), f"Input file {self.file} does not exist."
#     
#     def __enter__(self):
#         """Context manager entrance method."""
#         return self
# 
#     def __exit__(self, exc_type, exc_val, traceback):
#         """Context manager exit method."""
#         return False
# 
#     def __repr__(self):
#         """Representation request method."""
#         return self.file
#         
#     def touch(self) -> None:
#         """Creates empty file.
# 
#         This class mehtod is analagous to UNIX's ``touch`` command.
#         
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with File("file_name.txt") as file:
#             ...     file.touch()
#             ...
#             >>> # or
#             >>> 
#             >>> file = File("file_name.txt")
#             >>> file.touch()
#         """
#         if os.path.exists(self.file):
#             print(f"The file: {self.file} already exists.")
#         else:
#             with open(self.file,'w') as _:
#                 pass
#         return None
#     
#     def abs_path(self,
#                  follow_sym_links: bool = False
#                 ) -> str:
#         """Returns absolute path of file.
#         
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with File("file_name.txt") as file:
#             ...     file.touch()
#             ...     print(file.abs_path())
#             ...
#             "abspath/to/file_namt.txt"
#             >>> 
#             >>> # or
#             >>> 
#             >>> file = File("file_name.txt")
#             >>> file.abs_path()
#             "abspath/to/file_namt.txt"
#         
#         Arguments:
#             follow_sym_links: If set to true, the absolute path of the symlinked file is returned.
#         """
#         if follow_sym_links and os.path.exists(self.file):
#             return os.path.abspath(os.path.realpath(self.file))
#         elif os.path.exists(self.file):
#             return os.path.abspath(self.file)
#         else:
#             self.touch()
#             file_path: str = os.path.abspath(self.file)
#             os.remove(self.file)
#             return file_path
#     
#     def rm_ext(self,
#                ext: str = "") -> str:
#         """Removes file extension from the file.
#         
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with File("file_name.txt") as file:
#             ...     file.touch()
#             ...     print(file.rm_ext())
#             ...
#             "file_name"
#             >>> 
#             >>> # or
#             >>> 
#             >>> file = File("file_name.txt")
#             >>> file.rm_ext()
#             "file_name"
#         
#         Arguments:
#             ext: File extension.
#         
#         Returns:
#             Filename as string with no extension.
#         """
#         if ext:
#             ext_len: int = len(ext)
#             return self.file[:-(ext_len)]
#         elif self.ext:
#             ext_len = len(self.ext)
#             return self.file[:-(ext_len)]
#         else:
#             return self.file[:-(4)]
#         
#     def write_txt(self,
#                   txt: str = ""
#                  ) -> None:
#         """Writes/appends text to file.
# 
#         NOTE:
#             Text written to file is ALWAYS utf-8 encoded.
#         
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with File("file_name.txt") as file:
#             ...     file.write_txt("<Text to be written>")
#             ...
#             >>> # or
#             >>> 
#             >>> file = File("file_name.txt")
#             >>> file.write_txt("<Text to be written>")
#         
#         Arguments:
#             txt: Text/string to be written to file.
#         """
#         with open(self.file, mode="a", encoding='utf-8') as tmp_file:
#             tmp_file.write(txt)
#             tmp_file.close()
#         return None
# 
#     def file_parts(self,
#                    ext: str = ""
#                   ) -> Tuple[str,str,str]:
#         """Similar to MATLAB's ``fileparts``, this function splits a file and its path into its constituent parts:
# 
#             * file path
#             * filename
#             * extension
#         
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with File("file_name.txt") as file:
#             ...     print(file.file_parts())
#             ...
#             ("path/to/file", "filename", ".txt")
#             >>> 
#             >>> # or
#             >>> 
#             >>> file = File("file_name.txt")
#             >>> file.file.file_parts()
#             ("path/to/file", "filename", ".txt")
#         
#         Arguments:
#             ext: File extension, needed if the file extension of file object is longer than 4 characters.
#         
#         Returns:
#             Tuple: 
#                 * Absolute file path, excluding filename.
#                 * Filename, excluding extension.
#                 * File extension.
#         """
#         file: File = self.file
#         file: str = self.abs_path()
#         
#         path, _filename = os.path.split(file)
#         
#         if ext:
#             ext_num: int = len(ext)
#             _filename: str = _filename[:-(ext_num)]
#             [filename, _ext] = os.path.splitext(_filename)
#         elif self.ext:
#             ext: str = self.ext
#             ext_num: int = len(ext)
#             _filename: str = _filename[:-(ext_num)]
#             [filename, _ext] = os.path.splitext(_filename)
#         else:
#             [filename, ext] = os.path.splitext(_filename)
#         
#         return (path, 
#                 filename, 
#                 ext)

# class WorkDir:
#     """Working directory base class that instantiates ``WorkDir`` objects that creates and manipulates working directories.
# 
#     Attributes:
#         work_dir: Input working directory.
#         parent_dir: Parent directory.
# 
#     Usage example:
#             >>> # Using class object as context manager
#             >>> ## Create work directory , then clean-up (remove it)
#             >>> with WorkDir(work_dir="/path/to/working_directory", use_cwd=False) as work:
#             ...     work.mkdir()
#             ...     work.rmdir(rm_parent=False)
#             ...
#             >>> # or
#             >>>
#             >>> work = TmpDir(work_dir="/path/to/working_directory", 
#             ...               use_cwd=False)
#             >>> work.mkdir()
#             >>> work
#             "/path/to/working_directory"
#             >>> work.rmdir(rm_parent=False)
# 
#     Arguments:
#         work_dir: Working directory name/path. This directory need not exist at runtime.
#         use_cwd: Use current working directory as the parent directory.
#     """
#     __slots__ = [ 
#                     "work_dir", 
#                     "parent_dir"
#                 ]
# 
#     def __init__(self, 
#                  work_dir: str,
#                  use_cwd: bool = False
#                 ) -> None:
#         """Initialization method for the ``WorkDir`` base class.
# 
#         Usage example:
#             >>> # Using class object as context manager
#             >>> ## Create work directory , then clean-up (remove it)
#             >>> with WorkDir(work_dir="/path/to/working_directory", use_cwd=False) as work:
#             ...     work.mkdir()
#             ...     work.rmdir(rm_parent=False)
#             ...
#             >>> # or
#             >>>
#             >>> work = TmpDir(work_dir="/path/to/working_directory", 
#             ...               use_cwd=False)
#             >>> work.mkdir()
#             >>> work
#             "/path/to/working_directory"
#             >>> work.rmdir(rm_parent=False)
#         
#         Arguments:
#             work_dir: Working directory name/path. This directory need not exist at runtime.
#             use_cwd: Use current working directory as the parent directory.
#         """
#         self.work_dir: str = work_dir
#         self.parent_dir: str = os.path.dirname(self.work_dir)
# 
#         if use_cwd:
#             _cwd: str = os.getcwd()
#             self.work_dir: str = os.path.join(_cwd,self.work_dir)
#             self.parent_dir: str = os.path.dirname(self.work_dir)
# 
#     def __enter__(self):
#         """Context manager entrance method."""
#         return self
# 
#     def __exit__(self, exc_type, exc_val, traceback):
#         """Context manager exit method."""
#         return False
#     
#     def __repr__(self):
#         """Representation request method."""
#         return self.work_dir
#     
#     def mkdir(self) -> None:
#         """Makes/creates the working directory.
# 
#         This class method is analogous to UNIX's ``mkdir -p`` command and option combination.
# 
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with WorkDir(work_dir="/path/to/working_directory", use_cwd=False) as work:
#             ...     work.mkdir()
#             ...
#             >>> # or
#             >>>
#             >>> work = TmpDir(work_dir="/path/to/working_directory", 
#             ...               use_cwd=False)
#             >>> work.mkdir()
#             >>> work
#             "/path/to/working_directory"
#         """
#         if not os.path.exists(self.work_dir):
#             return os.makedirs(self.work_dir)
#         else:
#             print("Working directory already exists.")
#             return None
#     
#     def rmdir(self, 
#               rm_parent: bool = False
#              ) -> None:
#         """Removes working directory, and the parent directory if indicated to do so.
# 
#         This class method is analogous to UNIX's ``rm -rf`` command and option combination.
# 
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with WorkDir(work_dir="/path/to/working_directory", use_cwd=False) as work:
#             ...     work.mkdir()
#             ...     work.rmdir(rm_parent=False)
#             ...
#             >>> # or
#             >>>
#             >>> work = TmpDir(work_dir="/path/to/working_directory", 
#             ...               use_cwd=False)
#             >>> work.mkdir()
#             >>> work.rmdir(rm_parent=False)
# 
#         Arguments:
#             rm_parent: Removes parent directory as well.
#         """
#         if rm_parent and os.path.exists(self.parent_dir):
#             return shutil.rmtree(self.parent_dir,ignore_errors=True)
#         elif os.path.exists(self.work_dir):
#             return shutil.rmtree(self.work_dir,ignore_errors=True)
#         else:
#             print("Working directory does not exist.")
#             return None

# class TmpDir(WorkDir):
#     """Temporary directory class that creates (random) temporary directories and files given a parent directory. This class inherits methods from the ``WorkDir`` base class.
#     
#     Attributes:
#         work_dir: Input temproary working directory.
#         parent_dir: Parent directory of the specified temproary directory.
#     
#     Usage example:
#             >>> with TmpDir("/path/to/temporary_directory",False) as tmp_dir:
#             ...     tmp_dir.mk_tmp_dir()
#             ...     # do more stuff
#             ...     tmp_dir.rm_tmp_dir(rm_parent=False)
#             ...
#             >>> # or
#             >>>
#             >>> tmp_dir = TmpDir("/path/to/temporary_directory")
#             >>> tmp_dir
#             "/path/to/temporary_directory"
#             >>> tmp_dir.rm_tmp_dir(rm_parent=False)
#         
#     Arguments:
#         tmp_dir: Temporary parent directory name/path.
#         use_cwd: Use current working directory as working direcory.
#     """
#     
#     def __init__(self,
#                  tmp_dir: str,
#                  use_cwd: bool = False
#                 ) -> None:
#         """Initialization method for the TmpDir child class.
#         
#         Usage example:
#             >>> with TmpDir("/path/to/temporary_directory",False) as tmp_dir:
#             ...     tmp_dir.mk_tmp_dir()
#             ...     # do more stuff
#             ...     tmp_dir.rm_tmp_dir(rm_parent=False)
#             ...
#             >>> # or
#             >>> tmp_dir = TmpDir("/path/to/temporary_directory")
#             >>> tmp_dir
#             "/path/to/temporary_directory"
#             >>> tmp_dir.rm_tmp_dir(rm_parent=False)
#         
#         Arguments:
#             tmp_dir: Temporary parent directory name/path.
#             use_cwd: Use current working directory as working direcory.
#         """
#         _n: int = 10000
#         tmp_dir: str = os.path.join(tmp_dir,'tmp_dir_' + 
#                                     str(random.randint(0,_n)))
#         
#         if use_cwd:
#             _cwd: str = os.getcwd()
#             tmp_dir = os.path.join(_cwd,tmp_dir)
#         
#         super(TmpDir, self).__init__(tmp_dir,
#                                      use_cwd)
#     
#     class TmpFile(File):
#         """Sub-class of ``TmpDir`` class, which creates and manipulates temporary files via inheritance from the ``File`` object base class.
#         
#         Attributes:
#             file: Temporary file name.
#             ext: File extension of input file. If no extension is provided, it is inferred.
#         
#         Usage example:
#             >>> tmp_directory = TmpDir("/path/to/temporary_directory")
#             >>>
#             >>> temp_file = TmpDir.TmpFile(tmp_directory.tmp_dir,
#             ...                             ext="txt")
#             ...
#             >>> temp_file
#             "/path/to/temporary_directory/temporary_file.txt"
#         
#         Arguments:
#             tmp_dir: Temporary directory name.
#             tmp_file: Temporary file name.
#             ext: Temporary directory file extension.
#         """
# 
#         def __init__(self,
#                      tmp_dir: str,
#                      tmp_file: Optional[str] = "",
#                      ext: Optional[str] = "",
#                     ) -> None:
#             """Initialization method for the TmpFile sub-class that inherits from the ``File`` base class, allowing for the creation and maninuplation of temporary files.
#             
#             Usage example:
#                 >>> tmp_directory = TmpDir("/path/to/temporary_directory")
#                 >>>
#                 >>> temp_file = TmpDir.TmpFile(tmp_directory.tmp_dir,
#                 ...                             ext="txt")
#                 ...
#                 >>> temp_file
#                 "/path/to/temporary_directory/temporary_file.txt"
#             
#             Arguments:
#                 tmp_dir: Temporary directory name.
#                 tmp_file: Temporary file name.
#                 ext: File extension.
#             """          
#             if tmp_file:
#                 pass
#             else:
#                 _n: int = 10000 
#                 tmp_file: str = "tmp_file_" + str(random.randint(0,_n))
#             
#             if ext:
#                 tmp_file: str = tmp_file + f".{ext}"
# 
#             tmp_file: str = os.path.join(tmp_dir, 
#                                          tmp_file)
# 
#             super(TmpDir.TmpFile, self).__init__(tmp_file, 
#                                                  ext)
        
# class NiiFile(File):
#     """NIFTI file class specific for NIFTI files which inherits class methods from the ``File`` base class.
#     
#     Usage example:
#         >>> # Using class object as context manager
#         >>> with NiiFile("file.nii") as nii:
#         ...     print(nii.file_parts())
#         ...
#         ("path/to/file", "file", ".nii")
#         >>> 
#         >>> # or
#         >>> 
#         >>> nii = NiiFile("file.nii")
#         >>> nii
#         "file.nii"
#         >>> nii.abs_path()
#         "abspath/to/file.nii"
#         >>> 
#         >>> nii.rm_ext()
#         "file"
#         >>>
#         >>> nii.file_parts()
#         ("path/to/file", "file", ".nii")
#     
#     Arguments:
#         file: Path to NIFTI file.
#         
#     Raises:
#         InvalidNiftiFileError: Exception that is raised in the case **IF** the specified NIFTI file exists, but is an invalid NIFTI file.
#     """
# 
#     def __init__(self,
#                  file: str,
#                  assert_exists: bool = False,
#                  validate_nifti: bool = False
#                 ) -> None:
#         """Initialization method for the NiiFile class.
#         
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with NiiFile("file.nii") as nii:
#             ...     print(nii.file_parts())
#             ...
#             "abspath/to/file.nii"
#             "file"
#             ("path/to/file", "file", ".nii")
#             >>> 
#             >>> # or
#             >>> 
#             >>> nii = NiiFile("file.nii")
#             >>> nii
#             "file.nii"
#             >>> nii.abs_path()
#             "abspath/to/file.nii"
#             >>> 
#             >>> nii.rm_ext()
#             "file"
#             >>>
#             >>> nii.file_parts()
#             ("path/to/file", "file", ".nii")
#         
#         Arguments:
#             file: Path to NIFTI file.
#             assert_exists: Asserts that the specified input file must exist. 
#             validate_nifti: Validates the input NIFTI file if it exists.
#         
#         Raises:
#             InvalidNiftiFileError: Exception that is raised in the case **IF** the specified NIFTI file exists, but is an invalid NIFTI file.
#         """
#         self.file: str = file
#         
#         if self.file.endswith(".nii.gz"):
#             self.ext: str = ".nii.gz"
#         elif self.file.endswith(".nii"):
#             self.ext: str = ".nii"
#         else:
#             self.ext: str = ".nii.gz"
#             self.file: str = self.file + self.ext
#         
#         super(NiiFile, self).__init__(self.file)
# 
#         if assert_exists:
#             assert os.path.exists(self.file), f"Input file {self.file} does not exist."
#         
#         if validate_nifti and os.path.exists(self.file):
#             try:
#                 _: nib.Nifti1Header = nib.load(filename=self.file)
#             except Exception as e:
#                 print(e)
#                 raise InvalidNiftiFileError(f"The NIFTI file {self.file} is not a valid NIFTI file.")
#         
#     # Overwrite several File base class methods
#     def touch(self) -> None:
#         """This class method is not implemented and will simply return None, and is not relevant/needed for NIFTI files.
#         """
#         return None
# 
#     def write_txt(self,
#                   txt: str = "",
#                   header_field: str = "intent_name"
#                  ) -> None:
#         """This class method writes relevant information to the NIFTI file header.
#         This is done by writing text to either the ``descrip`` or ``intent_name``
#         field of the NIFTI header.
# 
#         NOTE:
#             * The ``descrip`` NIFTI header field has limitation of 24 bytes - meaning that only a string of 24 characters can be written without truncation.
#             * The ``intent_name`` NIFTI header field has limitation of 16 bytes - meaning that only a string of 16 characters can be written without truncation.
#         
#         Usage example:
#             >>> # Using class object as context manager
#             >>> with NiiFile("file.nii") as nii:
#             ...     nii.write_txt(txt='Source NIFTI',
#             ...                   header_field='intent_name')
#             ...
#             >>> # or
#             >>> 
#             >>> nii = NiiFile("file.nii")
#             >>> nii.write_txt(txt='Source NIFTI',
#             ...               header_field='intent_name')
# 
#         Arguments:
#             txt: Input text to be added to the NIFTI file header.
#             header_field: Header field to have text added to.
# 
#         Raises:
#             NiftiFileIOWarning: Warning that is raised if the byte character limit is surpassed for the specified header field.
#         """
#         img: nib.Nifti1Header = nib.load(self.file)
#         header_field: str = NiiHeaderField(header_field).name
# 
#         if header_field == 'descrip':
#             if len(txt) >= 24:
#                 raise NiftiFileIOWarning(f"The input string is longer than the allowed limit of 24 bytes/characters for the '{header_field}' header field.")
#             img.header['descrip'] = txt
#         elif header_field == 'intent_name':
#             if len(txt) >= 16:
#                 raise NiftiFileIOWarning(f"The input string is longer than the allowed limit of 16 bytes/characters for the '{header_field}' header field.")
#             img.header['intent_name'] = txt
#         return None

# class LogFile(File):
#     """Convenience class that creates a log file object for logging purposes. Due to how this class is constructed - its 
#     intended use case requires that this class is instantiated/called once and ONLY once.
#     
#     Once a class instance has been instantiated, then it and its associated methods can be used.
#     
#     Attributes:
#         log_file: Log filename.
#     
#     Usage examples:
#         >>> log = LogFile("file.log",False)
#         >>> log
#         "file.log"
# 
#     Arguments:
#         file: Log filename (need not exist at runtime).
#         print_to_screen: If true, prints output to standard output (stdout) as well.
#     """
#     
#     def __init__(self,
#                  log_file: str = "",
#                  print_to_screen: bool = False
#                 ) -> None:
#         """Initialization method for the LogFile class. Initiates logging and its associated methods (from the ``logging`` module).
#         
#         Usage examples:
#             >>> log = LogFile("file.log",False)
#             >>> log
#             "file.log"
#         
#         Arguments:
#             file: Log filename (need not exist at runtime).
#             print_to_screen: If true, prints output to standard output (stdout) as well.
#         """
#         self.log_file: str = log_file
#         
#         # Set-up logging to file
#         logging.basicConfig(level=logging.INFO,
#                             format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
#                             datefmt='%d-%m-%y %H:%M:%S',
#                             filename=self.log_file,
#                             filemode='a')
#         
#         # Define a Handler which writes INFO messages or higher to the sys.stderr
#         if print_to_screen:
#             self.console = logging.StreamHandler()
#             self.console.setLevel(logging.INFO)
#             logging.getLogger().addHandler(self.console)
#             
#         # Define logging
#         self.logger = logging.getLogger(__name__)
#         super(LogFile, self).__init__(self.log_file)
#         
#     def info(self,
#             msg: str = "") -> None:
#         """Writes information to log file.
#         
#         Usage examples:
#             >>> log = LogFile("file.log")
#             >>> log.info("<str>")
#         
#         Arguments:
#             msg: String to be printed to log file.
#         """
#         self.logger.info(msg)
#         
#     def debug(self,
#             msg: str = "") -> None:
#         """Writes debug information to file.
#         
#         Usage examples:
#             >>> log = LogFile("file.log")
#             >>> log.debug("<str>")
#         
#         Arguments:
#             msg: String to be printed to log file.
#         """
#         self.logger.debug(msg)
#         
#     def error(self,
#             msg: str = "") -> None:
#         """Writes error information to file.
#         
#         Usage examples:
#             >>> log = LogFile("file.log")
#             >>> log.error("<str>")
#         
#         Arguments:
#             msg: String to be printed to log file.
#         """
#         self.logger.error(msg)
#         
#     def warning(self,
#             msg: str = "") -> None:
#         """Writes warnings to file.
#         
#         Usage examples:
#             >>> log = LogFile("file.log")
#             >>> log.warning("<str>")
#         
#         Arguments:
#             msg: String to be printed to log file.
#         """
#         self.logger.warning(msg)
#     
#     def log(self,
#             log_cmd: str = "") -> None:
#         """Log function for logging commands and messages to some log file.
#         
#         Usage examples:
#             >>> log = LogFile("file.log")
#             >>> log.log("<str>")
#             
#         Arguments:
#             log_cmd: Message to be written to log file
#         """
#         self.info(log_cmd)

class Command:
    """Creates a command and an empty command list for UNIX command line programs/applications. Primary use and
    use-cases are intended for the subprocess module and its associated classes (i.e. Popen/call/run).

    The input argument is a command (string), and a mutable list is returned (, that can later be appended to).

    NOTE: 
        The specified command used must be in system path.
    
    Attributes:
        command: Command to be performed on the command line.
        cmd_list: Mutable list that can be appended to.
    
    Usage example:
        >>> echo = Command("echo")
        >>> echo.opt("Hi!")
        >>> echo.opt("I have arrived!")
        >>>
        >>> echo
        echo Hi! I have arrived!
        >>>
        >>> echo.run()
        Hi! I have arrived!
    
    Arguments:
        command: Command to be used. 
    
    Returns:
        Mutable list that can be appended to.
    """
    __slots__ = [ 
                    "command", 
                    "cmd_list"
                ]

    def __init__(self,
                 command: str
                ) -> None:
        """Initialization method for the Command class. Initializes a command to be used on UNIX command line.
        The input argument is a command (string), and a mutable list is returned (, that can later
        be appended to).
        
        Usage example:
            >>> echo = Command("echo")
            >>> echo.opt("Hi!")
            >>> echo.opt("I have arrived!")
            >>>
            >>> echo
            echo 'Hi! I have arrived!'
            >>>
            >>> echo.run()
            Hi! I have arrived!
        
        Arguments:
            command: Command to be used. 
                NOTE: command used must be in system path
        
        Returns:
            Mutable list that can be appended to.
        """
        self.command: str = command
        self.cmd_list: List[str] = [ f"{self.command}" ]
    
    def __repr__(self):
        """NOTE: This returns a string represnted as a joined list of strings."""
        return ' '.join(self.cmd_list)
    
    def opt(self,
            cmd: str
           ) -> None:
        """Appends a string or command-line option to an existing command list.

        Usage example:
            >>> echo = Command("echo")
            >>> echo.opt(cmd="Hi!")
            >>> echo.opt(cmd="I have arrived!")
            >>>
            >>> echo
            echo 'Hi! I have arrived!'
            >>>
            >>> echo.run()
            Hi! I have arrived!

        Arguments:
            cmd: The command or option to be appended to the command list.
        """
        return self.cmd_list.append(cmd)
        
    def check_dependency(self,
                         err_msg: Optional[str] = None,
                         path_envs: Optional[List[str]] = []
                        ) -> bool:
        """Checks dependency of some command line executable. Should the 
        dependency not be met, then an exception is raised. Check the 
        system path should problems arise and ensure that the executable
        of interest is installed.
        
        Usage example:
            >>> figlet = Command("figlet")
            >>> figlet.check_dependency()   # Raises exception if not in system path
        
        Arguments:
            err_msg: Error message to print to screen.
            path_envs: List of directory paths to append to the system's 'PATH' variable.

        Returns:
            Returns True if dependency is met, raises exception otherwise.
        
        Raises:
            DependencyError: Dependency error exception is raised if the dependency is not met.
        """
        # Append to PATH environmental variable
        mod_path_env: str = os.environ['PATH']
        if path_envs:
            for path in path_envs:
                mod_path_env += os.pathsep + path

        os.environ['PATH'] = mod_path_env
        
        if not shutil.which(self.command):
            if err_msg:
                print(f"\n \t {err_msg} \n")
            else:
                print(f"\n \t The required dependency {self.command} is not installed or in the system path. \n")
            raise DependencyError(f"Command executable not found in system path: {self.command}.")
        else:
            return True
        
    def run(self,
            log: Optional[LogFile] = None,
            debug: bool = False,
            dryrun: bool = False,
            path_envs: List[str] = [],
            env: Dict = {},
            stdout: str = "",
            shell: bool = False
           ) -> Tuple[int,Union[File,None]]:
        """Uses python's built-in subprocess class to execute (run) a command from an input command list.
        The standard output and error can optionally be written to file.
        
        NOTE: 
            * The contents of the ``stdout`` output file will be empty if ``shell`` is set to True.
            * **IF** ``check_dependency`` was used with the ``path_envs`` argument, then the default system ``PATH`` variable has been updated to include the list of specified paths.
        
        Usage example:
            >>> # Create command and cmd_list
            >>> echo = Command("echo")
            >>> echo.opt("Hi!")
            >>> echo.opt("I have arrived!")
            >>>
            >>> # Run/execute command
            >>> echo.run()
            (0, '', '')
        
        Arguments:
            log: LogFile object
            debug: Sets logging function verbosity to DEBUG level
            dryrun: Dry run -- does not run task. Command is recorded to log file.
            path_envs: List of directory paths to append to the system's 'PATH' variable.
            env: Dictionary of environment variables to add to subshell.
            stdout: Output file to write standard output to.
            shell: Use shell to execute command.
            
        Returns:
            Tuple:
                * Return code for command execution.
                * Standard output writtent to file should the 'stdout' option be used.
                * Standard error writtent to file should the 'stdout' option be used.
        """
        # Create command str for log
        cmd: str = ' '.join(self.cmd_list)      # Join list for logging purposes
        
        if log:
            if debug:
                log.debug(f"Running:\t\t {cmd}")
            else:
                log.info(f"Running:\t\t {cmd}")
        
        if log:
            if dryrun:
                log.info("Performing command as dryrun")
                return (0,'','')
        
        # Append to PATH environmental variable
        mod_path_env: str = os.environ['PATH']
        if path_envs:
            for path in path_envs:
                mod_path_env += os.pathsep + path

        # Define environment variables
        merged_env: Dict = os.environ
        if env and path_envs:
            merged_env.update(env)
            merged_env['PATH'] = mod_path_env
        elif env:
            merged_env.update(env)
        elif path_envs:
            merged_env['PATH'] = mod_path_env
        
        # Execute/run command
        p: subprocess.Popen = subprocess.Popen(self.cmd_list,
                                               shell=shell,
                                               env=merged_env,
                                               stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)

        # Write log files
        out,err = p.communicate()
        out = out.decode('utf-8')
        err = err.decode('utf-8')

        # Write std output/error files
        if stdout:
            stderr: str = os.path.splitext(stdout)[0] + ".err"
                
            stdout: File = File(stdout)
            stderr: File = File(stderr)
            
            stdout.write_txt(out)
            stderr.write_txt(err)
        else:
            stdout: None = None
            stderr: None = None

        if p.returncode != 0:
            if log:
                log.error(f"command: {cmd} \n Failed with returncode {p.returncode}")
            else:
                print(f"command: {cmd} \n Failed with returncode {p.returncode}")

        if len(out) > 0:
            if log:
                if debug:
                    log.debug(out)
                else:
                    log.info(out)

        if len(err) > 0:
            if log:
                if debug:
                    log.info(err)
                else:
                    log.warning(err)
            else:
                print(f"ERROR: {err}")
        return (p.returncode, 
                stdout, 
                stderr)
