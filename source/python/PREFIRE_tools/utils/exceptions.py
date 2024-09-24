"""
Define some Python custom exceptions

This program requires Python version 3.6 or later, and is importable as a 
Python module.
"""


class NoValidInputData(Exception):
    "Raised when there are no valid input data found."
    pass


class NoValidNewInputData(Exception):
    "Raised when there are no new valid input data found."
    pass
