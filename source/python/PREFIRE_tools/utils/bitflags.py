"""
Module containing methods to work with "bitflags", where the individual binary
 bits in an integer datatype contain boolean info.

This program requires python version 3.6 or later, and is importable as a 
python module.
"""

  # From the Python standard library:
from itertools import chain

  # From other external Python packages:
import numpy as np

  # Custom utilities:


#--------------------------------------------------------------------------
def bit_meaning(int_v, offset):
   """Returns True if the bit at 'offset' is one, and False otherwise."""
   x = int_v.dtype.type(1)
   mask = x << offset
   if (int_v & mask) != 0:  # nonzero result (2**offset) if the bit is one
      return True
   else:
      return False


#--------------------------------------------------------------------------
def set_bit(int_v, offset):
   """Returns an integer with the bit at 'offset' set to one."""
   x = int_v.dtype.type(1)
   mask = x << offset
   return (int_v | mask)


#--------------------------------------------------------------------------
def clear_bit(int_v, offset):
    """Returns an integer with the bit at 'offset' cleared (set to zero)."""
    x = int_v.dtype.type(1)
    mask = ~(x << offset)
    return (int_v & mask)


#--------------------------------------------------------------------------
def toggle_bit(int_v, offset):
    """Returns an integer with the bit at 'offset' inverted, 0 -> 1, 1 -> 0."""
    x = int_v.dtype.type(1)
    mask = x << offset
    return (int_v ^ mask)


#--------------------------------------------------------------------------
def boolrep_to_bit(bool_v, int_v, offset):
   """Returns an integer with the bit at 'offset' set to one if 'bool_v' is
      True, and set to zero otherwise."""
   if bool_v:
       return set_bit(int_v, offset)
   else:
       return clear_bit(int_v, offset)


#--------------------------------------------------------------------------
def get_bit_from_bitflags_v(int_v, offset):
   """Given a scalar, list, or NumPy array of integers ('int_v'; each element of
      which contains multiple bitflags) and a bitflag ID ('offset'), returns a
      corresponding boolean variable (scalar, list, or NumPy array). If the
      input was a NumPy array, the output array retains its shape.  However, if
      nested lists were the input, the output list is a flattened version of
      that."""
   try:
       len_iv = len(int_v)  # Check whether this is a scalar
       try:
           shape_iv = int_v.shape
           tmp_iv0 = int_v.flatten()
           tmp_iv1 = np.array([bit_meaning(df, offset) for df in tmp_iv0])
           return np.reshape(tmp_iv1, shape_iv)
       except:
           # Input is likely a list or nested list:
           nested = any(isinstance(elem, list) for elem in int_v)
           if nested:
               tmp_iv0 = list(chain.from_iterable(int_v))
           else:
               tmp_iv0 = int_v
           return [bit_meaning(df, offset) for df in tmp_iv0]
   except TypeError:
       # Input is likely a scalar
       return bit_meaning(int_v, offset)


#--------------------------------------------------------------------------
def apply_bit_to_bitflags_v(offset, bool_v, int_v):
   """Given a bitflag ID ('offset'), a boolean variable ('bool_v'; scalar, list,
      or NumPy array; each element of which indicates the bit value/state), and
      an integer variable ('int_v'; scalar, list, or NumPy array; each element
      of which contains multiple bitflags), returns a modified copy of the
      integer variable (scalar, list, or NumPy array; a nested list will be
      flattened) with the specified bit set according to the boolean input.
      Note that 'bool_v' and 'int_v' must be the same size when flattened."""
   try:
       bvl = len(bool_v)  # Check whether this is a scalar
       bv_is_scalar = False
   except TypeError:
       bv_is_scalar = True
   try:
       ivl = len(int_v)  # Check whether this is a scalar
       iv_is_scalar = False
   except TypeError:
       iv_is_scalar = True

   if bv_is_scalar != iv_is_scalar:
       raise ValueError("(only one is scalar) The input boolean field must be "
                        "the same size/shape as the input integer bitflags "
                        "field.")

   if bv_is_scalar:
       # Input is likely a scalar:
       return boolrep_to_bit(bool_v, int_v, offset)
   else:
       # Input is likely a list or NumPy array:
       try:
           shape_bv = bool_v.shape
           tmp_bv = bool_v.flatten()
           bv_is_list = False
       except AttributeError:
           # Likely a list or nested list:
           nested = any(isinstance(elem, list) for elem in bool_v)
           if nested:
               tmp_bv = list(chain.from_iterable(bool_v))
           else:
               tmp_bv = bool_v.copy()
           bv_is_list = True

       try:
           shape_iv = int_v.shape
           tmp_iv = int_v.flatten()
           iv_is_list = False
       except AttributeError:
           # Likely a list or nested list:
           nested = any(isinstance(elem, list) for elem in int_v)
           if nested:
               tmp_iv = list(chain.from_iterable(int_v))
           else:
               tmp_iv = int_v.copy()
           iv_is_list = True

       tmp_iv1 = [boolrep_to_bit(bv, iv, offset) for bv, iv in
                  zip(tmp_bv, tmp_iv)]

       if iv_is_list:
           return tmp_iv1
       else:
           tmp_iv2 = np.array(tmp_iv1)
           return np.reshape(tmp_iv2, shape_iv)
