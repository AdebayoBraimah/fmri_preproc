# -*- coding: utf-8 -*-
"""Abstract base class file IO methods, functions and operations.
"""
import os
import subprocess
from warnings import warn

from shutil import (
    copy,
    copytree
)

from typing import (
    Any,
    List, 
    Tuple,
    Union
)

from abc import (
    ABC,
    abstractmethod
)

class IOBaseObj(ABC):
    """IO abstract base class (``ABC``) object that encapsulates methods related to file and directory manipulation. 

    This ``ABC`` cannot be directly instantiated, and **MUST** used by a child/sub-class that inherits from this class. 
    Additionally, several class methods (shown in abstract methods) **MUST** be overwritten when inheriting from this class. 

    Attributes:
        src: Input string that represents a file or directory.
    
    Abstract methods:
        abspath: Returns the absolute file path.
        sym_link: Creates a symbolic link with an absolute or relative file path.
        copy: Copies a file or recursively copies a directory.
    
    Usage example:
        >>> # Initialize child class and inherit 
        >>> #   from IOBaseObj ABC
        >>> class SomeFileClass(IOBaseObj):
        ...     def __init__(self, src):
        ...         super().__init__(src)
        ...     
        ...     # Overwrite IOBaseObj ABC methods
        ...     def abspath(self, follow_sym_links):
        ...         return super().abspath(follow_sym_links)
        ...
        ...     def sym_link(self, dst, relative):
        ...         return super().sym_link(dst, relative)
        ...     
        ...     def copy(self, dst):
        ...         return super().copy(dst)
        ...         

    Arguments:
        src: Input string that represents a file or directory.
    """
    __slots__ = [ "src" ]

    def __init__(self,
                 src: str) -> None:
        """Constructor that initializes ``IOBaseObj`` abstract base class."""
        self.src: str = src
        super(IOBaseObj, self).__init__()
    
    def __enter__(self):
        """Context manager entrance method."""
        return self
    
    def __exit__(self, exc_type, exc_val, traceback):
        """Context manager exit method."""
        return False
    
    def __repr__(self):
        """Representation request method."""
        return self.src
    
    @abstractmethod
    def abspath(self,
                follow_sym_links: bool = False
               ) -> Union[str,None]:
        """Returns the absolute file path.
        
        Usage example:
            >>> # Initialize child class and inherit 
            >>> #   from IOBaseObj ABC
            >>> class SomeFileClass(IOBaseObj):
            ...     def __init__(self, src):
            ...         super().__init__(src)
            ...     
            ...     # Overwrite IOBaseObj ABC methods
            ...     def abspath(self, follow_sym_links):
            ...         return super().abspath(follow_sym_links)
            ...
            ...     def sym_link(self, dst, relative):
            ...         return super().sym_link(dst, relative)
            ...     
            ...     def copy(self, dst):
            ...         return super().copy(dst)
            ...         
            >>> # Using class object as context manager
            >>> with SomeFileClass("file_name.txt") as file:
            ...     print(file.abspath())
            ...
            "abspath/to/file_namt.txt"
            >>>
            >>> # OR
            >>> file = SomeFileClass("file_name.txt")
            >>> file.abspath()
            "abspath/to/file_namt.txt"
        
        Arguments:
            follow_sym_links: If set to true, the absolute path of the symlinked file is returned.
        
        Returns:
            String that represents the absolute file path if it exists, otherwise ``None`` is returned.
        """
        if follow_sym_links and os.path.exists(self.src):
            return os.path.abspath(os.path.realpath(self.src))
        else:
            return os.path.abspath(self.src)
    
    @abstractmethod
    def sym_link(self, 
                 dst: str, 
                 relative: bool = False
                ) -> str:
        """Creates a symbolic link with an absolute or relative file path.

        Usage example:
            >>> # Initialize child class and inherit 
            >>> #   from IOBaseObj ABC
            >>> class SomeFileClass(IOBaseObj):
            ...     def __init__(self, src):
            ...         super().__init__(src)
            ...     
            ...     # Overwrite IOBaseObj ABC methods
            ...     def abspath(self, follow_sym_links):
            ...         return super().abspath(follow_sym_links)
            ...
            ...     def sym_link(self, dst, relative):
            ...         return super().sym_link(dst, relative)
            ...     
            ...     def copy(self, dst):
            ...         return super().copy(dst)
            ...         
            >>> # Using class object as context manager
            >>> with SomeFileClass("file_name.txt") as file:
            ...     linked_file: str = file.sym_link("file2.txt")
            ...     print(linked_file)
            ...
            "file2.txt"
            >>>
            >>> # OR
            >>> file = SomeFileClass("file_name.txt")
            >>> file.sym_link("file2.txt")
            "file2.txt"

        Arguments:
            dst: Destination file path.
            relative: Symbolically link the file using a relative path.

        Returns:
            String that reprents the sym linked file path.
        """
        src: str = self.src
        src: str = self.abspath(follow_sym_links=True)

        if os.path.exists(dst):
            os.remove(dst)
            warn(f"WARNING: Symlinked file of the name {dst} already exists. It is being replaced.")
        
        if relative:
            src: str = os.path.relpath(src, dst)
        
        # Create command list
        cmd: List[str] = [ "ln", 
                           "-s",
                           f"{src}",
                           f"{dst}"
                         ]
        
        # Execute command
        p: subprocess.Popen = subprocess.Popen(cmd)
        _: Tuple[Any] = p.communicate()
        dst: str = os.path.abspath(dst)
        return dst

    @abstractmethod
    def copy(self, 
             dst: str
            ) -> str:
        """Copies file or recursively copies a directory to some destination.

        Usage example:
            >>> # Initialize child class and inherit 
            >>> #   from IOBaseObj ABC
            >>> class SomeFileClass(IOBaseObj):
            ...     def __init__(self, src):
            ...         super().__init__(src)
            ...     
            ...     # Overwrite IOBaseObj ABC methods
            ...     def abspath(self, follow_sym_links):
            ...         return super().abspath(follow_sym_links)
            ...
            ...     def sym_link(self, dst, relative):
            ...         return super().sym_link(dst, relative)
            ...     
            ...     def copy(self, dst):
            ...         return super().copy(dst)
            ...         
            >>> # Using class object as context manager
            >>> with SomeFileClass("file_name.txt") as file:
            ...     new_file: str = file.copy("file2.txt")
            ...     print(new_file)
            ...
            "/abs/path/to/file2.txt"
            >>>
            >>> # OR
            >>> file = SomeFileClass("file_name.txt")
            >>> file.copy("file2.txt")
            "/abs/path/to/file2.txt"

        Arguments:
            dst: Destination file path.

        Return:
            String that corresponds to the copied file or directory.
        """
        src: str = self.abspath(follow_sym_links=True)
        if os.path.isfile(src):
            return os.path.abspath(copy(src=src, dst=dst))
        elif os.path.isdir(src):
            return os.path.abspath(copytree(src=src, dst=dst))
