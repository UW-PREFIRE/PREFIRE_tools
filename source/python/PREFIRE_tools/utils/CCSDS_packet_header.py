"""
Contains routines that provide useful information about and/or operate on
CCSDS headers.

This program requires Python version 3.6 or later, and is importable as a 
Python module.

Binary format information from 'D-103781, PREFIRE Instrument Command and
 Telemetry Dictionary (CTD)', dated 29 Sep 2022 (RevB).
"""

  # From the Python standard library:
from collections import OrderedDict
import struct

  # From other external Python packages:

  # Custom utilities:


CCSDS_sync_marker = int("0x352EF853", 16)
bytesize_of_CCSDS_sync_marker = 4

payload_CCSDS_hdr = OrderedDict(CCSDS_Sync_Marker = (0, 3, "0x352ef853"),
                                CCSDS_Info_A = (4, 4, "0x03"),
                                APID_LSB = (5, 5, None),
                                CCSDS_Info_B = (6, 7, None),
                                N_Pktd_Bytes_m_1 = (8, 9, None))

#%struct_fmt_ref = ('', 'B', 'H', '', 'I')


#--------------------------------------------------------------------------
def CCSDS_header_specs(verbose=False):
    """Determine some useful information about a CCSDS header (in general)."""

    bytesize_of_CCSDS_hdr = 0
#%    CCSDS_hdr_unpack_fmtstr = '>'  # Assume big-endian byte order
    for key in payload_CCSDS_hdr:
        B_start, B_stop, nominal_hexvalue = payload_CCSDS_hdr[key]
        nB = B_stop-B_start+1
        bytesize_of_CCSDS_hdr += nB
#%        CCSDS_hdr_unpack_fmtstr += struct_fmt_ref[nB]

    if verbose:
        print("General CCSDS header info: bytesize = {:d}".format(
                                                        bytesize_of_CCSDS_hdr))

    return bytesize_of_CCSDS_hdr


#--------------------------------------------------------------------------
def seek_CCSDS_sync_marker(in_f):
    """Find next CCSDS sync marker and reposition the file just before it."""

    # Assuming that the 4-byte CCSDS sync marker binary sequence is very
    # improbable to find outside a real CCSDS header, this routine allows
    # packets with malformed or missing CCSDS headers to be gracefully skipped.
    #
    # Returns the number of bytes traversed to find the next CCSDS sync marker,
    #  or -1 if an EOF is encountered

    # Search for the next CCSDS sync marker using a 1-byte stride:
    search_Bcount = 0
    while True:
        sync_marker_buffer = in_f.read(bytesize_of_CCSDS_sync_marker)
        if len(sync_marker_buffer) != bytesize_of_CCSDS_sync_marker:  # EOF
            return -1
        sm_tmp = struct.unpack(">I", sync_marker_buffer)[0]
        if sm_tmp == CCSDS_sync_marker:
            in_f.seek(-bytesize_of_CCSDS_sync_marker,
                      1)  # Reposition at the start of the CCSDS header
            if search_Bcount > 0:
                print(("Extended CCSDS sync marker search; skipped {:d} "
                       "bytes.".format(search_Bcount)))
            return search_Bcount
            break
        else:
            in_f.seek(-bytesize_of_CCSDS_sync_marker+1,
                      1)  # Reposition one byte from last read position
        search_Bcount += 1


#--------------------------------------------------------------------------
def about_CCSDS_header(hdr_buffer):
    """Determine some useful information about a given CCSDS header."""


    ib, ie, _ = payload_CCSDS_hdr["CCSDS_Sync_Marker"]
    sync_marker = int.from_bytes(hdr_buffer[ib:ie+1], "big", signed=False)

    ib, ie, _ = payload_CCSDS_hdr["CCSDS_Info_A"]
    APID_MSB = int.from_bytes(hdr_buffer[ib:ie+1], "big", signed=False)
    APID_MSB &= ~(1 << 3)  # Extract appropriate (3) bits

    ib, ie, _ = payload_CCSDS_hdr["APID_LSB"]
    APID_LSB = int.from_bytes(hdr_buffer[ib:ie+1], "big", signed=False)
   
    APID_str_repr = "{0:08b}".format(APID_MSB)+"{0:08b}".format(APID_LSB)
    APID = int(APID_str_repr, 2)

    ib, ie, _ = payload_CCSDS_hdr["N_Pktd_Bytes_m_1"]
    remaining_packet_bytes = 1+int.from_bytes(hdr_buffer[ib:ie+1], "big",
                                              signed=False)

    return (sync_marker, APID, remaining_packet_bytes)
