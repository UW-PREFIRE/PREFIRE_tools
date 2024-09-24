"""
Utility routines for dealing with time/dates.

This program requires Python version 3.6 or later, and is importable as a 
Python module.
"""

# From the Python standard library:
import os
import sys
import datetime

  # From other external Python packages:
import numpy as np

  # Custom utilities:
import PREFIRE_tools.utils.unit_conv as uc
import PREFIRE_tools.filepaths as filepaths


one_s_td = datetime.timedelta(seconds=1)


#--------------------------------------------------------------------------
def _parse_IERS_bulletinC():
    """Reads/parses a local copy of IERS Bulletin C, obtained from
       https://data.iana.org/time-zones/data/leap-seconds.list
       and returns a dictionary with useful fields.
       *** That file must be kept up-to-date (new one every ~6 months) ***"""

    leap_s_info = {}

    # Continuous_time (similar to TAI) datetime of the epoch
    #  1900-01-01T00:00:00 "faux-UTC" -- but the leap second correction
    #  is 0 for this date, so the naive version is equal to the aware version:
    ctime_of_epoch1900_DT = datetime.datetime(1900, 1, 1, 0, 0, 0,
                                              tzinfo=datetime.timezone.utc)
    leap_s_info["ctime_of_epoch1900_DT"] = ctime_of_epoch1900_DT

    # Filepath of IERS Bulletin C:
    ls_fpath = os.path.join(filepaths.package_ancillary_data_dir,
                            "leap-seconds.list")

    # Determine and store reference datetimes at each leap second ("Ls")
    #  modification:
    ref_ctime1900_atLs_s = np.loadtxt(ls_fpath, dtype="float64", comments='#',
                   usecols=0)  # [seconds since 1900-01-01T00:00:00 "faux-UTC"]
    ref_ctime1900_atLs_td = ref_ctime1900_atLs_s*one_s_td
    leap_s_info["refDT_atLs"] = ctime_of_epoch1900_DT+ref_ctime1900_atLs_td

    # Determine and store continuous_time vs UTC offset at each reference leap
    #  second ("Ls") modification:
    tmp = np.loadtxt(ls_fpath, dtype="int8", comments='#', usecols=1)  # [s]
    leap_s_info["refDT_atLs_ctimeOffsetFromUTC_td"] = tmp*one_s_td

    # Note that dictionary keys "refDT_atLs" and
    #  "refDT_atLs_ctimeOffsetFromUTC_td" together represent a reference table
    #  of cumulative (ctime - UTC) offsets due to leap-seconds.
    return leap_s_info


#--------------------------------------------------------------------------
def init_leap_s_for_ctimeRefEpoch(ref_epoch_for_ctime,
                                  epoch_for_ctime_is_actual_UTC):
    """Initializes a leap-second reference dictionary at a given
       continuous_time reference epoch."""

    # *IMPORTANT*: 1) This routine will likely not work correctly for a
    #                 continuous_time reference epoch that is within ~10 seconds
    #                 of a leap-second modification.
    #              2) This routine should only be used for epochs in 1972 or
    #                 after, due to how UTC is defined.

    #-- Input parameters:
    # ref_epoch_for_ctime : (list/array of 6 ints) Reference epoch where
    #                                               ctime = 0,
    #                   [YYYY, MM, DD, hh, mm, ss], e.g., [1993, 1, 1, 0, 0, 0]
    # epoch_for_ctime_is_actual_UTC : (boolean) Is the ctime reference epoch
    #                                     actual-UTC (i.e., leap-second aware),
    #                                     versus "faux-UTC"?
    #                             If True, the practical consequences are:
    #                               * can use less-ambiguous units such as
    #                                 "seconds since 1993-01-01T00:00:00 UTC"
    #                               * the (ctime - UTC) offset will only contain
    #                                 any leap-seconds after the ctime epoch
    #                                 (not the full leap-second cumulative sum)

    # In general, the mathematical relationships between leap-second-aware UTC
    #  (i.e., actual-UTC) and ctime (epoch on/before ~1958) are:
    #      ctime = UTC + cumulative_sum_of_leap_seconds
    #      UTC = ctime - cumulative_sum_of_leap_seconds

    # Obtain leap-second reference table fields:
    leap_s_info = _parse_IERS_bulletinC()
    leap_s_info["epoch_for_ctime_is_actual_UTC"] = epoch_for_ctime_is_actual_UTC

    # Determine ctime epoch ("Ep") datetime and the full (ctime - UTC) offset
    #  at that epoch:
    refDT_atEp = datetime.datetime(*ref_epoch_for_ctime,
                                   tzinfo=datetime.timezone.utc)
    itmp = np.nonzero(leap_s_info["refDT_atLs"] <= refDT_atEp)
    i_e = itmp[0][-1]
    leap_s_info["ref_ctimeOffsetFromUTC_atEp_td"] = (
           leap_s_info["refDT_atLs_ctimeOffsetFromUTC_td"][i_e])  # full offset
    leap_s_info["ref_ctimeOffsetFromUTC_atEp_s"] = (
          leap_s_info["ref_ctimeOffsetFromUTC_atEp_td"].total_seconds())  # [s]
    leap_s_info["ref_DT_atEp"] = refDT_atEp

    # Continuous_time ([seconds from the given epoch]) of each leap-second
    #  modification, as well as the ctime (at the given epoch) offset from UTC
    #  of each leap-second modification:
    if epoch_for_ctime_is_actual_UTC:
        tmp0_td = ([x-leap_s_info["ref_ctimeOffsetFromUTC_atEp_td"] for x in
                              leap_s_info["refDT_atLs_ctimeOffsetFromUTC_td"]])
        tmp_s = [x.total_seconds() for x in tmp0_td]  # [s]
        leap_s_info["ctime_atEp_offsetFromUTC_atLs_td"] = tmp0_td
        leap_s_info["ctime_atEp_offsetFromUTC_atLs_s"] = tmp_s  # [s]

        tmp1_td = ([x-leap_s_info["ref_DT_atEp"] for x in
                                                    leap_s_info["refDT_atLs"]])
        ref_ctime_atEp_and_atLs_td = ([x+y-one_s_td for x, y in
                                                        zip(tmp1_td, tmp0_td)])
        leap_s_info["ref_ctime_atEp_and_atLs_s"] = np.array(
                [x.total_seconds() for x in ref_ctime_atEp_and_atLs_td])  # [s]
    else:
        tmp0_td = leap_s_info["refDT_atLs_ctimeoffsetFromUTC_td"]
        tmp_s = [x.total_seconds() for x in tmp0_td]  # [s]
        leap_s_info["ctime_atEp_offsetFromUTC_atLs_td"] = tmp0_td
        leap_s_info["ctime_atEp_offsetFromUTC_atLs_s"] = tmp_s  # [s]

        tmp1_td = ([x-leap_s_info["ref_DT_atEp"] for x in
                                                    leap_s_info["refDT_atLs"]])
        ref_ctime_atEp_and_atLs_td = ([x-one_s_td for x in tmp1_td])
        leap_s_info["ref_ctime_atEp_and_atLs_s"] = np.array(
                [x.total_seconds() for x in ref_ctime_atEp_and_atLs_td])  # [s]

#^
#    for i in range(len(leap_s_info["refDT_atLs"])):
#        print('rrr',i,leap_s_info["ref_ctime_atEp_and_atLs_s"][i],
#              leap_s_info["ctime_atEp_offsetFromUTC_atLs_s"][i])
#^

    return leap_s_info


#--------------------------------------------------------------------------
def _change_ctime_epoch_type(ctime, units, leap_s_info, input_is_naive):
    """Converts continuous_time value(s) between those referenced to an
       actual-UTC epoch and those referenced to a "faux-UTC" epoch,
       where 'epoch' is otherwise the same."""
    # === INPUTS ===
    # ctime :  continuous_time scalar value or NumPy array of values
    # units :  units of input/output ctime values, valid values:
    #                                      's', 'seconds', 'ms', 'milliseconds'
    # leap_s_info :  dictionary with useful leap-second-related fields
    # input_is_naive :  Is the input ctime "naive" (i.e., referenced to a
    #                    "faux-UTC" epoch)?
    #
    # === OUTPUT ===
    # returns a ctime scalar value or NumPy array of values

    # Convert value to input/output units, if needed:
    x = leap_s_info["ref_ctimeOffsetFromUTC_atEp_s"]  # [s]
    if units == "ms" or units == "milliseconds":
        x *= uc.s_to_ms  # [s] -> [ms]
    elif not (units == 's' or units == "seconds"):
        print(("ERROR: Given ctime units ({}) are unsupported.".format(units)))
        sys.exit(1)

    if input_is_naive:
        return ctime-x  # [units]
    else:
        return ctime+x  # [units]


#--------------------------------------------------------------------------
def ctimeN_to_ctime(ctimeN, units, leap_s_info):
    """Converts ctimeN value(s) (referenced to a "faux-UTC" epoch) to
       ctime value(s) (referenced to an actual-UTC epoch), where 'epoch' is
       otherwise the same."""
    # === OUTPUT ===
    # returns a ctime scalar value or NumPy array of values

    return _change_ctime_epoch_type(ctimeN, units, leap_s_info,
                                    input_is_naive=True)  # [units]


#--------------------------------------------------------------------------
def ctime_to_ctimeN(ctime, units, leap_s_info):
    """Converts ctime value(s) (referenced to an actual-UTC epoch) to
       ctimeN value(s) (referenced to a "faux-UTC" epoch), where 'epoch' is
       otherwise the same."""
    # === OUTPUT ===
    # returns a ctimeN scalar value or NumPy array of values

    return _change_ctime_epoch_type(ctime, units, leap_s_info,
                                    input_is_naive=False)  # [units]


#--------------------------------------------------------------------------
def ctime_to_UTC_DT(ctime, units, leap_s_info):
    """Converts ctime (referenced to the epoch given in 'leap_s_info') to
        actual-UTC (which takes leap-seconds into account) datetime
        object(s)."""
    # === INPUTS ===
    # ctime :  ctime value(s) (scalar or array-like), **in ascending order**
    # units :  units of input ctime values, valid values:
    #                                      's', 'seconds', 'ms', 'milliseconds'
    # leap_s_info :  dictionary with useful leap-second-related fields
    #
    # === OUTPUTS ===
    # Tuple (UTC_DT, ctime_minus_UTC), where
    #   UTC_DT :  A single (or NumPy array of) actual-UTC datetime object(s)
    #              (leap-second-aware)
    #   ctime_minus_UTC :  [s] scalar value or NumPy array of (ctime - UTC)
    #                        difference

    # Convert ctime value(s) to units of seconds, if needed:
    if units == "ms" or units == "milliseconds":
        x0 = ctime*uc.ms_to_s  # [ms] -> [s]
    elif not (units == 's' or units == "seconds"):
        print(("ERROR: The given ctime units ({}) are "
               "unsupported.".format(units)))
        sys.exit(1)
    else:
        x0 = ctime  # [s]

    # Operate a bit differently if input ctime is a scalar value:
    if isinstance(x0, (list, tuple, np.ndarray)):  # "array-like"
        x = np.array(x0)
        originally_a_scalar = False
    else:  # Assumed to be a scalar value
        x = np.array([x0])
        originally_a_scalar = True

    n_arr = len(x)
    n_lse = len(leap_s_info["ctime_atEp_offsetFromUTC_atLs_td"])

    itmp = np.nonzero(leap_s_info["ref_ctime_atEp_and_atLs_s"] <= x[0])
    i_lse_b = itmp[0][-1]  # Leap-second-reference index for first ctime value
    ib = 0
    ind = []
    if i_lse_b == n_lse-1:  # Last leap-second-reference index
        ind = [(0, n_arr-1, i_lse_b)]
    else:
        for i_lse in range(i_lse_b, n_lse-1):
            itmp = (np.nonzero(x[ib:n_arr] <
                            leap_s_info["ref_ctime_atEp_and_atLs_s"][i_lse+1]))

            if len(itmp[0]) >= 1:
                ee = itmp[0][-1]
                ie = ib+ee
                ind.append((ib, ie, i_lse))
            elif i_lse == n_lse-2:  # must use last ref index
                ind.append((ib, n_arr-1, n_lse-1))

            ib = ie+1
            if ib > n_arr-1:
                break

    UTC_DT = np.zeros((n_arr,)).astype(datetime.datetime)
    ctime_minus_UTC = np.zeros((n_arr,))
    for i in range(0, len(ind)):
        ib, ie, i_lse = ind[i]
        UTC_naive_DT = np.array([leap_s_info["ref_DT_atEp"]+
                          datetime.timedelta(seconds=xx) for xx in x[ib:ie+1]])
        UTC_DT[ib:ie+1] = (UTC_naive_DT-
                        leap_s_info["ctime_atEp_offsetFromUTC_atLs_td"][i_lse])
        tmp_td = leap_s_info["ctime_atEp_offsetFromUTC_atLs_td"][i_lse]
        ctime_minus_UTC[ib:ie+1] = tmp_td.total_seconds()  # [s]

    if originally_a_scalar:
        return (UTC_DT[0], ctime_minus_UTC[0])
    else:
        return (UTC_DT, ctime_minus_UTC)


#--------------------------------------------------------------------------
def UTC_DT_to_ctime(UTC_DT, units, leap_s_info):
    """Converts actual-UTC (which takes leap-seconds into account) datetime
        object(s) to ctime (referenced to the epoch given in 'leap_s_info')."""
    # === INPUTS ===
    # UTC_DT :  A single (or NumPy array of) actual-UTC datetime object(s)
    #            (leap-second-aware, and UTC_DT.tzinfo must be "UTC"),
    #            **in ascending time order**
    # units :  units of output ctime values, valid values:
    #                                      's', 'seconds', 'ms', 'milliseconds'
    # leap_s_info :  dictionary with useful leap-second-related fields
    #
    # === OUTPUTS ===
    # Tuple (ctime, ctime_minus_UTC), where
    #   ctime :  [units as specified] single (or NumPy array of) ctime value(s)
    #   ctime_minus_UTC :  [s] scalar value or NumPy array of (ctime - UTC)
    #                        difference

    # Operate a bit differently if input UTC_DT is a scalar value:
    if isinstance(UTC_DT, (list, tuple, np.ndarray)):  # "array-like"
        x = np.array(UTC_DT)
        originally_a_scalar = False
    else:  # Assumed to be a scalar value
        x = np.array([UTC_DT])
        originally_a_scalar = True

    n_arr = len(x)
    n_lse = len(leap_s_info["refDT_atLs"])

    itmp = np.nonzero(leap_s_info["refDT_atLs"] <= x[0])
    i_lse_b = itmp[0][-1]  # Leap-second-reference index for first UTC_DT obj
    ib = 0
    ind = []
    if i_lse_b == n_lse-1:  # Last leap-second-reference index
        ind = [(0, n_arr-1, i_lse_b)]
    else:
        for i_lse in range(i_lse_b, n_lse-1):
            itmp = np.nonzero(x[ib:n_arr] < leap_s_info["refDT_atLs"][i_lse+1])
            if len(itmp[0]) >= 1:
                ee = itmp[0][-1]
                ie = ib+ee
                ind.append((ib, ie, i_lse))
            elif i_lse == n_lse-2:  # must use last ref index
                ind.append((ib, n_arr-1, n_lse-1))

            ib = ie+1
            if ib > n_arr-1:
                break

    ctime = np.zeros((n_arr,))
    ctime_minus_UTC = np.zeros((n_arr,))
    for i in range(0, len(ind)):
        ib, ie, i_lse = ind[i]
        tmp_td = leap_s_info["ctime_atEp_offsetFromUTC_atLs_td"][i_lse]
        ctime_minus_UTC[ib:ie+1] = tmp_td.total_seconds()  # [s]

        ctime0 = np.array([(xx-leap_s_info["ref_DT_atEp"]).total_seconds()
                                                   for xx in x[ib:ie+1]]) # [s]

        ctime[ib:ie+1] = ctime0+ctime_minus_UTC[ib:ie+1]  # [s]

    # Convert ctime value(s) to desired units, if needed:
    if units == "ms" or units == "milliseconds":
        ctime *= uc.s_to_ms  # [s] -> [ms]
    elif not (units == 's' or units == "seconds"):
        print(("ERROR: The given ctime units ({}) are "
               "unsupported.".format(units)))
        sys.exit(1)

    if originally_a_scalar:
        return (ctime[0], ctime_minus_UTC[0])
    else:
        return (ctime, ctime_minus_UTC)
