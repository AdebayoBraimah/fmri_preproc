# -*- coding: utf-8 -*-
"""Logging utility.
"""
import logging

from fmri_preproc.utils.io import File

class LogFile(File):
    """Convenience class that creates a log file object for logging purposes. Due to how this class is constructed - its 
    intended use case requires that this class is instantiated/called once and ONLY once.
    
    Once a class instance has been instantiated, then it and its associated methods can be used.
    
    Attributes:
        log_file: Log filename.
    
    Usage examples:
        >>> log = LogFile("file.log",False)
        >>> log
        "file.log"

    Arguments:
        file: Log filename (need not exist at runtime).
        print_to_screen: If true, prints output to standard output (stdout) as well.
    """
    
    def __init__(self,
                 log_file: str = "",
                 print_to_screen: bool = False
                ) -> None:
        """Initialization method for the LogFile class. Initiates logging and its associated methods (from the ``logging`` module).
        
        Usage examples:
            >>> log = LogFile("file.log",False)
            >>> log
            "file.log"
        
        Arguments:
            file: Log filename (need not exist at runtime).
            print_to_screen: If true, prints output to standard output (stdout) as well.
        """
        self.log_file: str = log_file
        
        # Set-up logging to file
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%d-%m-%y %H:%M:%S',
                            filename=self.log_file,
                            filemode='a')
        
        # Define a Handler which writes INFO messages or higher to the sys.stderr
        if print_to_screen:
            self.console = logging.StreamHandler()
            self.console.setLevel(logging.INFO)
            logging.getLogger().addHandler(self.console)
            
        # Define logging
        self.logger = logging.getLogger(__name__)
        super(LogFile, self).__init__(self.log_file)
        
    def info(self,
            msg: str = "") -> None:
        """Writes information to log file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.info("<str>")
        
        Arguments:
            msg: String to be printed to log file.
        """
        self.logger.info(msg)
        
    def debug(self,
            msg: str = "") -> None:
        """Writes debug information to file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.debug("<str>")
        
        Arguments:
            msg: String to be printed to log file.
        """
        self.logger.debug(msg)
        
    def error(self,
            msg: str = "") -> None:
        """Writes error information to file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.error("<str>")
        
        Arguments:
            msg: String to be printed to log file.
        """
        self.logger.error(msg)
        
    def warning(self,
            msg: str = "") -> None:
        """Writes warnings to file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.warning("<str>")
        
        Arguments:
            msg: String to be printed to log file.
        """
        self.logger.warning(msg)
    
    def log(self,
            log_cmd: str = "") -> None:
        """Log function for logging commands and messages to some log file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.log("<str>")
            
        Arguments:
            log_cmd: Message to be written to log file
        """
        self.info(log_cmd)