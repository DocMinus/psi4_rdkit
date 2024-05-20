"""
Module for logging PSI output (and only that) to a file with given filename.
Fileending will be replaced with log and a timestamp will be added.
"""

import logging
import os
import sys
from datetime import datetime

import psi4


class Psi4Logger:
    """
    A class for setting up logging and configuring the Psi4 output file.

    Args:
        filename (str): The name of the log file.

    Attributes:
        filename (str): The name of the log file.
        logfile (str): The absolute path to the log file.
        logger (logging.Logger): The logger object for logging messages.

    """

    def __init__(self, filename):
        self.now = datetime.now()
        self.filename = os.path.join(
            os.path.splitext(filename)[0]
            + "_"
            + self.now.strftime("%y%m%d_%H%M%S")
            + ".log",
        )
        self.logfile = os.path.join(
            os.path.dirname(os.path.abspath(sys.argv[0])),
            self.filename,
        )

        # Create logger
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

        # Create console handler and set level to info
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)
        self.logger.addHandler(handler)

        # Create file handler and set level to info
        handler = logging.FileHandler(self.logfile)
        handler.setLevel(logging.INFO)
        self.logger.addHandler(handler)

        # Set Psi4 output file
        psi4.core.set_output_file(self.logfile, False)
