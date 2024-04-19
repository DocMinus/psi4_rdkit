# example logger for Psi4; not used per se here but could be.
# could make a function out of this....
# using a different log function in the each script for now.
import logging
import sys

import psi4

# Create logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create console handler and set level to info
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
logger.addHandler(handler)

# Create file handler and set level to info
handler = logging.FileHandler("output.dat")
handler.setLevel(logging.INFO)
logger.addHandler(handler)

# Set Psi4 output file
psi4.core.set_output_file("output.dat", False)

# Now output from Psi4 will be written to 'output.dat' and printed to console
