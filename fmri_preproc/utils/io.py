# -*- coding: utf-8 -*-
"""File IO methods, functions and operations.
"""
import os
import nibabel as nib

from fmri_preproc.utils.enums import NiiHeaderField

from typing import (
    Optional,
    Tuple
)

class InvalidNiftiFileError(Exception):
    """Exception intended for invalid NIFTI files."""
    pass

class NiftiFileIOWarning(Warning):
    """Warning that is raised in the case of char byte overflow written to NIFTI file headers."""
    pass

class File:
    """File object base class. This class creates a ``File`` object that encapsulates a number of methods and properites for file and filename handling, and file manipulation.
    
    Attributes:
        file: Input file.
        ext: File extension of input file. If no extension is provided, it is inferred.

    Usage example:
        >>> # Using class object as context manager
        >>> with File("file_name.txt") as file:
        ...     file.touch()
        ...     file.write_txt("some text")
        ...
        >>> # or
        >>> 
        >>> file = File("file_name.txt")
        >>> file
        "file_name.txt"

    Arguments:
        file: Input file (need not exist at runtime/instantiated).
        ext: File extension of input file. If no extension is provided, it is inferred.
    """
    __slots__ = [ 
                    "file", 
                    "ext"
                ]
    
    def __init__(self,
                 file: str,
                 ext: Optional[str] = "",
                 assert_exists: bool = False
                ) -> None:
        """Initialization method for the File base class.
        
        Usage example:
            >>> # Using class object as context manager
            >>> with File("file_name.txt") as file:
            ...     file.touch()
            ...     file.write_txt("some text")
            ...
            >>> # or
            >>> 
            >>> file = File("file_name.txt")
            >>> file
            "file_name.txt"

        Arguments:
            file: Input file (need not exist at runtime/instantiated).
            ext: File extension of input file. If no extension is provided, it is inferred.
            assert_exists: Asserts that the specified input file must exist. 
        """
        self.file: str = file
        
        if ext:
            self.ext: str = ext
        elif self.file.endswith('.gz'):
            self.ext: str = self.file[-7:]
        else:
            _, self.ext = os.path.splitext(self.file)
        
        if assert_exists:
            assert os.path.exists(self.file), f"Input file {self.file} does not exist."
    
    def __enter__(self):
        """Context manager entrance method."""
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        """Context manager exit method."""
        return False

    def __repr__(self):
        """Representation request method."""
        return self.file
        
    def touch(self) -> None:
        """Creates empty file.

        This class mehtod is analagous to UNIX's ``touch`` command.
        
        Usage example:
            >>> # Using class object as context manager
            >>> with File("file_name.txt") as file:
            ...     file.touch()
            ...
            >>> # or
            >>> 
            >>> file = File("file_name.txt")
            >>> file.touch()
        """
        if os.path.exists(self.file):
            print(f"The file: {self.file} already exists.")
        else:
            with open(self.file,'w') as _:
                pass
        return None
    
    def abs_path(self,
                 follow_sym_links: bool = False
                ) -> str:
        """Returns absolute path of file.
        
        Usage example:
            >>> # Using class object as context manager
            >>> with File("file_name.txt") as file:
            ...     file.touch()
            ...     print(file.abs_path())
            ...
            "abspath/to/file_namt.txt"
            >>> 
            >>> # or
            >>> 
            >>> file = File("file_name.txt")
            >>> file.abs_path()
            "abspath/to/file_namt.txt"
        
        Arguments:
            follow_sym_links: If set to true, the absolute path of the symlinked file is returned.
        """
        if follow_sym_links and os.path.exists(self.file):
            return os.path.abspath(os.path.realpath(self.file))
        else:
            return os.path.abspath(self.file)
    
    def rm_ext(self,
               ext: str = "") -> str:
        """Removes file extension from the file.
        
        Usage example:
            >>> # Using class object as context manager
            >>> with File("file_name.txt") as file:
            ...     file.touch()
            ...     print(file.rm_ext())
            ...
            "file_name"
            >>> 
            >>> # or
            >>> 
            >>> file = File("file_name.txt")
            >>> file.rm_ext()
            "file_name"
        
        Arguments:
            ext: File extension.
        
        Returns:
            Filename as string with no extension.
        """
        if ext:
            ext_len: int = len(ext)
            return self.file[:-(ext_len)]
        elif self.ext:
            ext_len = len(self.ext)
            return self.file[:-(ext_len)]
        else:
            return self.file[:-(4)]
        
    def write_txt(self,
                  txt: str = ""
                 ) -> None:
        """Writes/appends text to file.

        NOTE:
            Text written to file is ALWAYS utf-8 encoded.
        
        Usage example:
            >>> # Using class object as context manager
            >>> with File("file_name.txt") as file:
            ...     file.write_txt("<Text to be written>")
            ...
            >>> # or
            >>> 
            >>> file = File("file_name.txt")
            >>> file.write_txt("<Text to be written>")
        
        Arguments:
            txt: Text/string to be written to file.
        """
        with open(self.file, mode="a", encoding='utf-8') as tmp_file:
            tmp_file.write(txt)
            tmp_file.close()
        return None

    def file_parts(self,
                   ext: str = ""
                  ) -> Tuple[str,str,str]:
        """Similar to MATLAB's ``fileparts``, this function splits a file and its path into its constituent parts:

            * file path
            * filename
            * extension
        
        Usage example:
            >>> # Using class object as context manager
            >>> with File("file_name.txt") as file:
            ...     print(file.file_parts())
            ...
            ("path/to/file", "filename", ".txt")
            >>> 
            >>> # or
            >>> 
            >>> file = File("file_name.txt")
            >>> file.file.file_parts()
            ("path/to/file", "filename", ".txt")
        
        Arguments:
            ext: File extension, needed if the file extension of file object is longer than 4 characters.
        
        Returns:
            Tuple: 
                * Absolute file path, excluding filename.
                * Filename, excluding extension.
                * File extension.
        """
        file: File = self.file
        file: str = os.path.abspath(file)
        
        path, _filename = os.path.split(file)
        
        if ext:
            ext_num: int = len(ext)
            _filename: str = _filename[:-(ext_num)]
            [filename, _] = os.path.splitext(_filename)
        elif self.ext:
            ext: str = self.ext
            ext_num: int = len(ext)
            _filename: str = _filename[:-(ext_num)]
            [filename, _] = os.path.splitext(_filename)
        else:
            [filename, ext] = os.path.splitext(_filename)
        
        return (path, 
                filename, 
                ext)

class NiiFile(File):
    """NIFTI file class specific for NIFTI files which inherits class methods from the ``File`` base class.
    
    Usage example:
        >>> # Using class object as context manager
        >>> with NiiFile("file.nii") as nii:
        ...     print(nii.file_parts())
        ...
        ("path/to/file", "file", ".nii")
        >>> 
        >>> # or
        >>> 
        >>> nii = NiiFile("file.nii")
        >>> nii
        "file.nii"
        >>> nii.abs_path()
        "abspath/to/file.nii"
        >>> 
        >>> nii.rm_ext()
        "file"
        >>>
        >>> nii.file_parts()
        ("path/to/file", "file", ".nii")
    
    Arguments:
        file: Path to NIFTI file.
        
    Raises:
        InvalidNiftiFileError: Exception that is raised in the case **IF** the specified NIFTI file exists, but is an invalid NIFTI file.
    """

    def __init__(self,
                 file: str,
                 assert_exists: bool = False,
                 validate_nifti: bool = False
                ) -> None:
        """Initialization method for the NiiFile class.
        
        Usage example:
            >>> # Using class object as context manager
            >>> with NiiFile("file.nii") as nii:
            ...     print(nii.file_parts())
            ...
            "abspath/to/file.nii"
            "file"
            ("path/to/file", "file", ".nii")
            >>> 
            >>> # or
            >>> 
            >>> nii = NiiFile("file.nii")
            >>> nii
            "file.nii"
            >>> nii.abs_path()
            "abspath/to/file.nii"
            >>> 
            >>> nii.rm_ext()
            "file"
            >>>
            >>> nii.file_parts()
            ("path/to/file", "file", ".nii")
        
        Arguments:
            file: Path to NIFTI file.
            assert_exists: Asserts that the specified input file must exist. 
            validate_nifti: Validates the input NIFTI file if it exists.
        
        Raises:
            InvalidNiftiFileError: Exception that is raised in the case **IF** the specified NIFTI file exists, but is an invalid NIFTI file.
        """
        self.file: str = file
        super(NiiFile, self).__init__(self.file)

        if self.file.endswith(".nii.gz"):
            self.ext: str = ".nii.gz"
        elif self.file.endswith(".nii"):
            self.ext: str = ".nii"
        else:
            self.ext: str = ".nii.gz"
            self.file: str = self.file + self.ext

        if assert_exists:
            assert os.path.exists(self.file), f"Input file {self.file} does not exist."
        
        if validate_nifti and os.path.exists(self.file):
            try:
                _: nib.Nifti1Header = nib.load(filename=self.file)
            except Exception as e:
                # print(e)
                raise InvalidNiftiFileError(f"The NIFTI file {self.file} is not a valid NIFTI file and raised the error {e}.")
        
    # Overwrite several File base class methods
    def touch(self) -> None:
        """This class method is not implemented and will simply return None, and is not relevant/needed for NIFTI files.
        """
        return None

    def write_txt(self,
                  txt: str = "",
                  header_field: str = "intent_name"
                 ) -> None:
        """This class method writes relevant information to the NIFTI file header.
        This is done by writing text to either the ``descrip`` or ``intent_name``
        field of the NIFTI header.

        NOTE:
            * The ``descrip`` NIFTI header field has limitation of 24 bytes - meaning that only a string of 24 characters can be written without truncation.
            * The ``intent_name`` NIFTI header field has limitation of 16 bytes - meaning that only a string of 16 characters can be written without truncation.
        
        Usage example:
            >>> # Using class object as context manager
            >>> with NiiFile("file.nii") as nii:
            ...     nii.write_txt(txt='Source NIFTI',
            ...                   header_field='intent_name')
            ...
            >>> # or
            >>> 
            >>> nii = NiiFile("file.nii")
            >>> nii.write_txt(txt='Source NIFTI',
            ...               header_field='intent_name')

        Arguments:
            txt: Input text to be added to the NIFTI file header.
            header_field: Header field to have text added to.

        Raises:
            NiftiFileIOWarning: Warning that is raised if the byte character limit is surpassed for the specified header field.
        """
        img: nib.Nifti1Header = nib.load(self.file)
        header_field: str = NiiHeaderField(header_field).name

        if header_field == 'descrip':
            if len(txt) >= 24:
                raise NiftiFileIOWarning(f"The input string is longer than the allowed limit of 24 bytes/characters for the '{header_field}' header field.")
            img.header['descrip'] = txt
        elif header_field == 'intent_name':
            if len(txt) >= 16:
                raise NiftiFileIOWarning(f"The input string is longer than the allowed limit of 16 bytes/characters for the '{header_field}' header field.")
            img.header['intent_name'] = txt
        return None