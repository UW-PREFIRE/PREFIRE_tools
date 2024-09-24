"""
Utility routines for selected file system operations.

This program requires Python version 3.6 or later, and is importable as a 
Python module.
"""

# From the Python standard library:
import os
import errno

  # From other external Python packages:

  # Custom utilities:


#--------------------------------------------------------------------------
def mkdir_p(path):
    """Emulates 'mkdir -p' functionality"""
    try:
        os.makedirs(path)
    except OSError as this_exception:
        if this_exception.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
