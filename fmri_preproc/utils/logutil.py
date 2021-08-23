# -*- coding: utf-8 -*-
"""Logging utility module.
"""
import logging
from datetime import datetime
from fmri_preproc.utils.fileio import File
from fmri_preproc.utils.enums import LogLevel

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
                 print_to_screen: bool = False,
                 level: str = "info" # TODO: Need enum here
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

        level: str = LogLevel(level.lower()).name
        
        if level == "info":
            level: logging.INFO = logging.INFO
        elif level == "debug":
            level: logging.DEBUG = logging.DEBUG
        
        logging.basicConfig(level=level,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%m-%d-%y %H:%M:%S',
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
             msg: str = "",
             use_header: bool = False) -> None:
        """Writes information to log file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.info("<str>")
        
        Arguments:
            msg: String to be printed to log file.
            use_header: Give log message a section header.
        """
        if use_header:
            self.logger.info(self._section_header(msg))
        else:
            self.logger.info(msg)
        
    def debug(self,
              msg: str = "",
              use_header: bool = False) -> None:
        """Writes debug information to file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.debug("<str>")
        
        Arguments:
            msg: String to be printed to log file.
            use_header: Give log message a section header.
        """
        if use_header:
            self.logger.debug(self._section_header(msg))
        else:
            self.logger.debug(msg)
    
    def critical(self,
                 msg: str,
                 use_header: bool = False) -> None:
        """Write critical messages/information to file.

        Usage example:
            >>> log = LogFile("file.log")
            >>> log.critical("<str>")

        Arguments:
            msg: String to be printed to log file.
            use_header: Give log message a section header. 
        """
        if use_header:
            self.logger.critical(self._section_header(msg))
        else:
            self.logger.critical(msg)
        
    def error(self,
              msg: str = "",
              use_header: bool = False) -> None:
        """Writes error information to file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.error("<str>")
        
        Arguments:
            msg: String to be printed to log file.
            use_header: Give log message a section header.
        """
        if use_header:
            self.logger.error(self._section_header(msg))
        else:
            self.logger.error(msg)
        
    def warning(self,
                msg: str = "",
                use_header: bool = False) -> None:
        """Writes warnings to file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.warning("<str>")
        
        Arguments:
            msg: String to be printed to log file.
            use_header: Give log message a section header.
        """
        if use_header:
            self.logger.warning(self._section_header(msg))
        else:
            self.logger.warning(msg)
    
    def log(self,
            log_cmd: str = "",
            use_header: bool = False) -> None:
        """Log function for logging commands and messages to some log file.
        
        Usage examples:
            >>> log = LogFile("file.log")
            >>> log.log("<str>")
            
        Arguments:
            log_cmd: Message to be written to log file.
            use_header: Give log message a section header.
        """
        if use_header:
            self.info(self._section_header(log_cmd))
        else:
            self.info(log_cmd)
    
    def _section_header(self,
                        msg: str) -> str:
        """Helper function that adds a section header that consists of
        a line break, with the date, time, and message string.

        Usage example:
            >>> _section_header("INFO: This is a test")

            -------------------------------------------------------------
            Mon Aug 23 13:34:21 2021 : INFO: This is a test
            -------------------------------------------------------------

        Arguments:
            msg: Message string to have section header.

        Returns:
            String that represents the message with header.
        """
        header: str = f"""\n
-------------------------------------------------------------
{datetime.now().ctime()} : {msg}
-------------------------------------------------------------
        """
        return header