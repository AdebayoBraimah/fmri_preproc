# -*- coding: utf-8 -*-
"""Command module for executing/running shell and command line tasks and operations.
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

from fmri_preproc.utils.fileio import File
from fmri_preproc.utils.logutil import LogFile

class DependencyError(Exception):
    """Exception intended for unment dependencies"""
    pass

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
    __slots__ = (
                    "command", 
                    "cmd_list"
                )

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
        return f"<{self.__class__.__name__} {' '.join(self.cmd_list)}>"
    
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
            shell: bool = False,
            raise_exc: bool = True
           ) -> Tuple[int,Union[str,None]]:
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
            log: ``LogFile`` object.
            debug: Sets logging function verbosity to DEBUG level.
            dryrun: Dry run -- does not run task. Command is recorded to log file.
            path_envs: List of directory paths to append to the system's 'PATH' variable.
            env: Dictionary of environment variables to add to subshell.
            stdout: Output file to write standard output to.
            shell: Use shell to execute command.
            raise_exec: If true, raises ``RuntimeError`` exception if the return code of the command is not 0.
            
        Returns:
            Tuple:
                * str: return code for command execution.
                * str: standard output writtent to file should the 'stdout' option be used.
                * str: standard error writtent to file should the 'stdout' option be used.
        
        Raises:
            RuntimeError: Exception that is raised if the return code of the command is not 0 and the ``raise_exc`` argument is set to ``True``.
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

            with File(src=stdout) as stout:
                with File(src=stderr) as sterr:
                    stout.write(out)
                    sterr.write(err)
                    stdout: str = stout.abspath()
                    stderr: str = sterr.abspath()
        else:
            stdout: str = None
            stderr: str = None

        if p.returncode != 0:
            if log:
                log.error(f"command: {cmd} \n Failed with returncode {p.returncode}")
            else:
                print(f"command: {cmd} \n Failed with returncode {p.returncode}")
            
            if raise_exc:
                raise RuntimeError(f"command: {cmd} \n\nFailed with returncode {p.returncode}")

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
    