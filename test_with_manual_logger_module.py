"""
just a test script to test the logger
"""

from mymodules.manual_logger import Logger


def print(*args, **kwargs):
    # not pylint compliant, but works to redirect print to the Logger without changing code
    logger.print_log(*args, **kwargs)


logger = Logger("testfile.txt")

# Calculate the heat of formation
print("\n------------------------------------------ \n\n\n")
print("test stuff")
print("\n----------------------------- \n\n\n")
