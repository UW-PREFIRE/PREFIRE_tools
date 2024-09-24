"""
Routines for dealing with common operations with numeric objects.

This program requires python version 3.6 or later, and is importable as a
python module.
"""

  # From the Python standard library:
import numpy as np
import copy

  # From other external Python packages:

  # Custom utilities:


def contains_NaN_or_Inf(obj):
    """Data object contains an unmasked NaN or Inf value (if not str-like)?"""
    answer = False
    if not (isinstance(obj, str) or isinstance(obj, bytes)
            or isinstance(obj, bytearray)):
        if not np.any(np.isfinite(obj)):  # Failed NaN and/or Inf check
            answer = True
    return answer


def replace_NaN_or_Inf(obj, fill_value=None, warn_with_moniker=None):
    """Attempt to 'repair' NaN or Inf values (in a numeric scalar or numpy
       array) using an appropriate mask or fill value. (Optionally, emit a
       warning message.)"""

    if warn_with_moniker is not None:
        print("WARNING: {} contains NaN or Inf (or equivalent) -- will attempt "
              "to replace those values with the appropriate mask or fill "
              "value.".format(warn_with_moniker))

    tmp_o = copy.deepcopy(obj)
    try:
        try:  # Some sort of numpy array?
            if isinstance(tmp_o, np.ma.core.MaskedArray):
                np.ma.masked_invalid(tmp_o)
            elif fill_value is not None:
                tmp_o[~np.isfinite(tmp_o)] = fill_value
            else:
                raise ValueError
        except:   # Scalar obj?
            if fill_value is not None:
                tmp_o = fill_value
            else:
                raise ValueError
    except:
        if warn_with_moniker is not None:
            emsg = warn_with_moniker
        else:
            emsg = "object"
        raise ValueError("{} contains NaN or Inf (or equivalent), and is not "
                         "auto-repairable".format(emsg))
    return tmp_o
