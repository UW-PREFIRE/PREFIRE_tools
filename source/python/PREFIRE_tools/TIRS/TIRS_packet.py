"""
Contains routines that provide information about and/or operate on flat-binary
 TIRS science packets.

This program requires Python version 3.6 or later, and is importable as a 
Python module.

Binary format information from 'D-103781, PREFIRE Instrument Command and
 Telemetry Dictionary (CTD)', dated 29 Sep 2022 (RevB).
"""

  # From the Python standard library:
from collections import OrderedDict
import struct
import binascii

  # From other external Python packages:
import numpy as np

  # Custom utilities:
import PREFIRE_tools.utils.unit_conv as uc
import PREFIRE_tools.utils.CCSDS_packet_header as Cph
from PREFIRE_tools.utils.time import ctimeN_to_ctime

  # offset DN, copied from (PREFIRE_L1)prefire_TIRS.m
_TIRS_offset_DN = 24000

payload_scipkt_APID = int("0x300", 16)
payload_scipkt_data = OrderedDict( \
                             Engineering_Data = [10, 89, None],
                             CnDH_Packet_ACKs = [90, 91, None],
                             CnDH_Packet_NACKs = [92, 93, None],
                             MoE_Packet_ACKs = [94, 95, None],
                             MoE_Packet_NACKs = [96, 97, None],
                             MoE_Packet_Sequence = [98, 99, None],
                             MoE_NVM_Faults = [100, 101, None],
                             Timestamp_Sec_Since_Epoch = [102, 105, None],
                             Timestamp_MilliSec_Since_PPS = [106, 109, (30, 0)],
                             Encoder_Position = [110, 111, None],
                             CRC_16b_A = [112, 113, (10, 111)],
                             ROIC_Image = [114, 1137, None],
                             CRC_16b_B = [1138, 1139, (114, 1137)])
payload_full_CRC_scipkt = OrderedDict(CRC_16b = [1140, 1141, (4, 1139)])

CRC_seed = int("0xFFFF", 16)

# Encoder position range (at end of each integration period) for each type of
#  nominal look:  [TIRS01, TIRS02]
# (Note that as of November 2023, this no longer works for TIRS01.)
space_enc_pos_range = [(4864-4, 4864+4), (4869-4, 4869+4)]
caltgt_enc_pos_range = [(3960-4, 3960+4), (3962-4, 3962+4)]
obstgt_enc_pos_range = [(3146-4, 3146+4), (3148-4, 3148+4)]

# (0)nominal obstgt, (1)possible obstgt, (10)nominal caltgt,
# (11)possible caltgt, (20)nominal space, (21)possible space,
# (30)skew based on encoder only, (31)skew based on block frame count 1,
# (32)skew based on block frame count 2, (33)skew based on encoder + ROIC DN,
# (34)skew assumed via disposition
OBSTGT_NOM, OBSTGT_POSS, CALTGT_NOM, CALTGT_POSS, SPACE_NOM, SPACE_POSS, \
    SKEW_VIA_ENCODER, SKEW_VIA_BLKCNT1, SKEW_VIA_BLKCNT2, SKEW_VIA_DN, \
    SKEW_VIA_DISP = (0, 1, 10, 11, 20, 21, 30, 31, 32, 33, 34)

N_XTRACK = 8  # Number of cross-orbital-track scenes
N_SPECTRAL = 64  # Number of spectral channels per scene

ROIC_tau = 0.7007  # [s] TIRS ROIC integration period

struct_fmt_ref = ('', 'B', 'H', '', 'I')

ib, ie, _ = payload_scipkt_data["Engineering_Data"]
n_engd = int((ie-ib+1)/2)
engd_unpack_fmtstr = '>'
engd_unpack_fmtstr += 'H'*n_engd

ib = payload_scipkt_data["CnDH_Packet_ACKs"][0]
ie = payload_scipkt_data["MoE_NVM_Faults"][1]
n_misc = int((ie-ib+1)/2)
misc_unpack_fmtstr = '>'
misc_unpack_fmtstr += 'H'*n_misc

n_CDH_therm = 2
n_ROIC_therm = 4
n_TIRS_therm = 8
n_MoE_therm = 1

B_ib, B_ie, b = payload_scipkt_data["Timestamp_MilliSec_Since_PPS"]
Ts_ms_i = (B_ie-B_ib+1)*8-b[1]-1  # Bit index to use

ROIC_unpack_fmtstr = '>'
ROIC_unpack_fmtstr += 'H'*(N_XTRACK*N_SPECTRAL)

T_cfg = {}
T_cfg["CDH"] = {"Eq_type": 0,
                "R_pullup": 442.,  # [kohm]
                "R_0": 440.,  # [kohm]
                "T_0_inv": 1./298.15,  # [1/K]
                "beta_inv": 1./4138.}  # [?]
T_cfg["ROIC"] = {"Eq_type": 1,
                 "R_pullup": 10.e3,  # [ohm]
                 "T_0": -599.6,  # [K]
                 "A": -3.617e-6,  # [C/(ohm)^2]
                 "B": 1.0112e-1}  # [C/ohm]
T_cfg["TIRS"] = {"Eq_type": 0,
                 "R_pullup": 10.,  # [kohm]
                 "R_0": 5.,  # [kohm]
                 "T_0_inv": 1./298.15,  # [1/K]
                 "beta_inv": 1./3891.}  # [1/K]
T_cfg["MoE"] = {"Eq_type": 0,
                "R_pullup": 442.,  # [kohm]
                "R_0": 440.,  # [kohm]
                "T_0_inv": 1./298.15,  # [1/K]
                "beta_inv": 1./4138.}  # [?]

I_cfg = {}
I_cfg["PDU_12V"] = {"V_factor": 3.3/4096.,  # [V]
                    "gain_amp": 5.006}  # [V/A]
I_cfg["PDU_SCBATT"] = {"V_factor": 3.3/4096.,  # [V]
                       "gain_amp": 2.503}  # [V/A]
I_cfg["MoE_12V"] = {"V_factor": 3.3/4096.,  # [V]
                    "gain_amp": 8.742}  # [V/A]


#--------------------------------------------------------------------------
def get_scipkt_bytebuffer_indices(bytesize_of_CCSDS_hdr):
    """Begin/end NumPy indices of 'payload_scipkt_data' fields in sci_buffer."""

    modified_odict = OrderedDict()
    for key in payload_scipkt_data:
        this_list = payload_scipkt_data[key]
        modified_odict[key] = (this_list[0]-bytesize_of_CCSDS_hdr,
                               this_list[1]-bytesize_of_CCSDS_hdr+1,
                               this_list[2])
   
    return (modified_odict)


#--------------------------------------------------------------------------
def get_T_value_from_DN(DN, DN_type):
    """Get engineering temperature value(s) from raw DN(s), given their type."""

    c = T_cfg[DN_type]  # Select cfg relevant for this DN_type

    if c["Eq_type"] == 0:
        DN_term = np.log((c["R_pullup"]/c["R_0"])*(DN/(4096.-DN)))  # [-]
        T_inv = c["T_0_inv"]+c["beta_inv"]*DN_term  # [1/K]
        T_val = 1./T_inv  # [K]
    elif c["Eq_type"] == 1:
        R_RTP = c["R_pullup"]*(DN/(4096.-DN))
        T_val = c["A"]*R_RTP*R_RTP+c["B"]*R_RTP+c["T_0"]+273.15  # [K]

    return T_val


#--------------------------------------------------------------------------
def get_I_value_from_DN(DN, DN_type):
    """Get engineering electrical current value(s) from raw DN(s), given
       their type."""

    c = I_cfg[DN_type]  # Select cfg relevant for this DN_type

    DN_term = DN*c["V_factor"]  # [V]
    I_val = DN_term/c["gain_amp"]  # [A]

    return I_val


def srd_to_sci_DN(srd_DN):
    """Helper to apply transformations to raw DN (in a "SRD" convention),
       resulting in DN in a "SCI"ence convention."""

    # The former is shaped (N_SPECTRAL*N_XTRACK,nframe), and the latter is shaped
    # ((N_SPECTRAL,N_XTRACK,nframe), and has 3 transformations applied (in this
    # order): (1) apply polarity
    #         (2) reshape (N_SPECTRAL*N_XTRACK,) to (N_SPECTRAL,N_XTRACK)
    #         (3) apply ROIC map
    #
    # Note that we do one extra step which is to add 2*offset to samples
    # with polarity flips; this should 're-align' samples with opposing
    # polarity and make it easier to plot / inspect. This is done as part
    # of step (1).
    #
    # Otherwise this is a nearly line by line copy of the MATLAB.
    #
    # Parameters
    # ----------
    # srd_DN: ndarray
    #     raw data shaped (nframe, N_SPECTRAL*N_XTRACK), unsigned 16-bit integer.
    #
    # Outputs
    # -------
    # sci_DN: ndarray
    #     raw data in shape (N_SPECTRAL,N_XTRACK,nframe), with the above
    #     transformations applied. Data is converted to signed 32-bit integers.

    nspectral, nscene, nframe = (N_SPECTRAL, N_XTRACK, srd_DN.shape[0])
    srd_tmp = np.zeros((nspectral, nscene, nframe), dtype="int32")

    for i,j in np.ndindex((nscene, nspectral)):
        k = N_XTRACK*j+i
        o = _TIRS_offset_DN * 2 * ((j+1) % 2)
        srd_tmp[j,i,:] = (-1)**(j+1) * srd_DN[:,k] + o

    sci_DN = np.zeros((nspectral, nscene, nframe), dtype="int32")

    # science order      <>  detector order
    sci_DN[ 0:32,3,:] =         srd_tmp[ 0:32,0,:]     #DS-A0a
    sci_DN[ 0:32,2,:] =         srd_tmp[32:64,0,:]     #CS-A0b
    sci_DN[ 0:32,5,:] = np.flip(srd_tmp[ 0:32,1,:], 0) #FS-A1a r
    sci_DN[ 0:32,4,:] = np.flip(srd_tmp[32:64,1,:], 0) #ES-A1b r
    sci_DN[32:64,0,:] =         srd_tmp[ 0:32,2,:]     #AL-B0a
    sci_DN[32:64,1,:] =         srd_tmp[32:64,2,:]     #BL-B0b
    sci_DN[ 0:32,1,:] =         srd_tmp[ 0:32,3,:]     #BS-B1a
    sci_DN[ 0:32,0,:] =         srd_tmp[32:64,3,:]     #AS-B1b
    sci_DN[32:64,4,:] = np.flip(srd_tmp[ 0:32,4,:], 0) #EL-C0a r
    sci_DN[32:64,5,:] = np.flip(srd_tmp[32:64,4,:], 0) #FL-C0b r
    sci_DN[32:64,2,:] =         srd_tmp[ 0:32,5,:]     #CL-C1a
    sci_DN[32:64,3,:] =         srd_tmp[32:64,5,:]     #DL-C1b
    sci_DN[ 0:32,7,:] = np.flip(srd_tmp[ 0:32,6,:], 0) #HS-D0a r
    sci_DN[ 0:32,6,:] = np.flip(srd_tmp[32:64,6,:], 0) #GS-D0b r
    sci_DN[32:64,6,:] = np.flip(srd_tmp[ 0:32,7,:], 0) #GL-D1a r
    sci_DN[32:64,7,:] = np.flip(srd_tmp[32:64,7,:], 0) #HL-D1b r

    return sci_DN


#--------------------------------------------------------------------------
def read_one_scipkt(in_f, bytesize_of_CCSDS_hdr, psp_buffinfo, nonscipkt_log,
                    leap_s_info, check_CRCs=True):
    """Read a TIRS science packet (no other packet types), and return its data
       along with some other info."""

    sm_bytes_traversed = Cph.seek_CCSDS_sync_marker(in_f)
    if sm_bytes_traversed < 0:  # EOF
        return (None, None, None, None, None, None)

    CCSDS_header_buffer = in_f.read(bytesize_of_CCSDS_hdr)

    if len(CCSDS_header_buffer) != bytesize_of_CCSDS_hdr:  # EOF
        return (None, None, None, None, None, None)

    while len(CCSDS_header_buffer) == bytesize_of_CCSDS_hdr:  # Not EOF
        sync_marker, APID, remaining_packet_bytes = \
                                    Cph.about_CCSDS_header(CCSDS_header_buffer)

        if sync_marker != Cph.CCSDS_sync_marker:
            msg = "should be CCSDS sync marker {} {}".format(sync_marker,
                                                         Cph.CCSDS_sync_marker)
            raise RuntimeError(msg)

        if APID != payload_scipkt_APID:
            #== Not a sci packet; find the next packet:
            print(("non-sci packet APID ({}), remaining packet bytes = "
                   "{:d}".format(hex(APID), remaining_packet_bytes)))
            info = [APID, in_f.tell()-bytesize_of_CCSDS_hdr,
                    remaining_packet_bytes+bytesize_of_CCSDS_hdr]
            sm_bytes_traversed = Cph.seek_CCSDS_sync_marker(in_f)
            info.append(int(remaining_packet_bytes == sm_bytes_traversed))
            nonscipkt_log.append(info)
            if sm_bytes_traversed < 0:  # EOF
                return (None, None, None, None, None, None)
            elif sm_bytes_traversed != remaining_packet_bytes:
                print ("WARNING: this was not a valid packet.")
            CCSDS_header_buffer = in_f.read(bytesize_of_CCSDS_hdr)  # Next hdr
            continue  # Continue to next packet
        else:
            #== This is a sci packet, but need to check the CRCs:
            sci_buffer = in_f.read(remaining_packet_bytes)
            if check_CRCs:
                fullpkt_buffer = CCSDS_header_buffer+sci_buffer

                  # Check topmost-level CRC first:
                ib, ie, rb_t = tuple(payload_full_CRC_scipkt["CRC_16b"])
                CRC_16b = int.from_bytes(fullpkt_buffer[ib:ie+1], "big",
                                         signed=False)
                this_CRC = binascii.crc_hqx(fullpkt_buffer[rb_t[0]:rb_t[1]+1],
                                            CRC_seed)
                top_CRC_check_failed = ( this_CRC != CRC_16b )

                  # Check ROIC CRC:
                ib, ie, rb_t = tuple(payload_scipkt_data["CRC_16b_B"])
                CRC_16b = int.from_bytes(fullpkt_buffer[ib:ie+1], "big",
                                         signed=False)
                this_CRC = binascii.crc_hqx(fullpkt_buffer[rb_t[0]:rb_t[1]+1],
                                            CRC_seed)
                ROIC_CRC_check_failed = ( this_CRC != CRC_16b )

                  # Check ENG CRC:
                ib, ie, rb_t = tuple(payload_scipkt_data["CRC_16b_A"])
                CRC_16b = int.from_bytes(fullpkt_buffer[ib:ie+1], "big",
                                         signed=False)
                this_CRC = binascii.crc_hqx(fullpkt_buffer[rb_t[0]:rb_t[1]+1],
                                            CRC_seed)
                ENG_CRC_check_failed = ( this_CRC != CRC_16b )

                if any([top_CRC_check_failed, ROIC_CRC_check_failed,
                        ENG_CRC_check_failed]):
                    print("sci packet failed CRC checks")
                    info = [APID, in_f.tell()-bytesize_of_CCSDS_hdr,
                            remaining_packet_bytes+bytesize_of_CCSDS_hdr]
                    sm_bytes_traversed = Cph.seek_CCSDS_sync_marker(in_f)
                    if top_CRC_check_failed:
                        ichk = 2
                    elif ROIC_CRC_check_failed and ENG_CRC_check_failed:
                        ichk = 3
                    elif ENG_CRC_check_failed:
                        ichk = 4
                    else:
                        ichk = 5
                    info.append(ichk)
                    nonscipkt_log.append(info)
                    CCSDS_header_buffer = in_f.read(bytesize_of_CCSDS_hdr)
                    continue  # Continue to next packet

        break

    # Timestamp (at end of each integration period):
    ib, ie, _ = psp_buffinfo["Timestamp_Sec_Since_Epoch"]
    ctimeN_ms = float(int.from_bytes(sci_buffer[ib:ie], "big",
                                        signed=False)*uc.int_s_to_ms)  # [ms]
    ib, ie, _ = psp_buffinfo["Timestamp_MilliSec_Since_PPS"]
    n = int.from_bytes(sci_buffer[ib:ie], "big", signed=False)
    PPS_ctime_mismatch = ((n >> Ts_ms_i) & 1)
    ms = n
    ms &= ~(1 << Ts_ms_i)  # [ms] Extract 'Ts_ms_i' bits
    ctimeN_ms += float(ms)  # [ms]

    ctime_s = ctimeN_to_ctime(ctimeN_ms*uc.ms_to_s, 's', leap_s_info)  # [s]

    # Encoder position (at end of each integration period):
    ib, ie, _ = psp_buffinfo["Encoder_Position"]
    encoder_pos = int.from_bytes(sci_buffer[ib:ie], "big", signed=False)

    return (CCSDS_header_buffer, sci_buffer, encoder_pos, ctime_s,
            PPS_ctime_mismatch)


#--------------------------------------------------------------------------
def extract_from_sci_buffer(sci_buffer_l, psp_buffinfo):
    """Extract (and return) information from a list of TIRS science data
       buffers."""

    n_sb = len(sci_buffer_l)

    EngData = np.empty((n_sb, n_engd), dtype="uint16")
    srd_DN = np.empty((n_sb, N_SPECTRAL*N_XTRACK),
                      dtype="uint16")
    misc_tlm = np.empty((n_sb, n_misc), dtype="uint16")

   #== First, extract larger chunks of data and convert them to NumPy arrays:

    ibA, ieA, _ = psp_buffinfo["Engineering_Data"]
    ibB, ieB, _ = psp_buffinfo["ROIC_Image"]
    ibC, ieC = (psp_buffinfo["CnDH_Packet_ACKs"][0],
                psp_buffinfo["MoE_NVM_Faults"][1])

    for i in range(n_sb):
        EngData[i,:] = np.array(struct.unpack(engd_unpack_fmtstr,
                                              sci_buffer_l[i][ibA:ieA]))

          # Leave in raw 1D array order:
        srd_DN[i,:] = np.array(struct.unpack(ROIC_unpack_fmtstr,
                                      sci_buffer_l[i][ibB:ieB]), dtype="uint16")

        misc_tlm[i,:] = np.array(struct.unpack(misc_unpack_fmtstr,
                                      sci_buffer_l[i][ibC:ieC]), dtype="uint16")

   #== Then expand/retrieve some of those values for output:

    CDH_T = get_T_value_from_DN(EngData[:,2:4], "CDH").astype("float32")
    ROIC_T = get_T_value_from_DN(EngData[:,4:8], "ROIC").astype("float32")
    TIRS_T = get_T_value_from_DN(EngData[:,8:16], "TIRS").astype("float32")
    MoE_T = get_T_value_from_DN(EngData[:,37], "MoE").astype("float32")

    PDU_12V_I = get_I_value_from_DN(EngData[:,18], "PDU_12V").astype("float32")
    PDU_SCBATT_I = get_I_value_from_DN(EngData[:,26], "PDU_SCBATT").astype(
                                                                      "float32")
    MoE_12V_I = get_I_value_from_DN(EngData[:,32], "MoE_12V").astype("float32")

    sci_DN = srd_to_sci_DN(srd_DN)

    return (EngData, misc_tlm, sci_DN, CDH_T, ROIC_T, TIRS_T, MoE_T, PDU_12V_I,
            PDU_SCBATT_I, MoE_12V_I)
