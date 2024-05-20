"""
Module for logging PSI4 output to both, console and a file. 
Ending will be replaced with .log and a timestamp will be added.

Usage:
    logger = Logger("log")
    logger.print_log("This is a log message")
    
To minimize changing code, you can add the following function to your code:

def print(*args, **kwargs):
    # Redirect standard print function
    logger.print_log(*args, **kwargs) 

"""

import os
import sys
from datetime import datetime


class Logger:
    """
    A class for logging information to both, console and a file

    Args:
        filename (str): The name of the log file.

    Attributes:
        filename (str): The name of the log file.

    """

    def __init__(self, filename):
        self.now = datetime.now()
        self.filename = os.path.join(
            os.path.splitext(filename)[0]
            + "_"
            + self.now.strftime("%y%m%d_%H%M%S")
            + ".log",
        )

    def print_log(self, *args):
        """
        the acutal function to print the log
        """

        print(*args)

        with open(
            os.path.join(
                os.path.dirname(os.path.abspath(sys.argv[0])),
                self.filename,
            ),
            "a",
            encoding="utf-8",
        ) as file:
            print(*args, file=file)

    def set_filename(self, new_filename):
        """
        Set a new filename for the log file; necesary? perhaps not used in the end.
        """
        self.filename = new_filename
