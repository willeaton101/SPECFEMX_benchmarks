# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Functions for loading and processing data from YSPEC, SPECFEMX, AxiSEM3D or NMSYN
# ______________________________________________________________________________________________________________________
import numpy as np
import obspy

def get_data(type, filepath, meta, measure="D", label=None, yspec="", grav_yspec="", nograv_fp=""):
    # ==================================================================================================================
    # DESCRIPTION:
    # Load data from YSPEC, NMSYNG or SPECFEMX simulation and return to user as obspy stream
    # INPUTS:
    #    type      [str]        - Type of simulation data. Can be "spfmx", "nmsyn" or "yspec"
    #    filepath  [str]        - Directory filepath holding the dataset
    #    meta      [metadata]   - PlotMetadata object (see plot_meta_data_class.py)
    #    measure   [str][opt.]  - Displacement (D), Velocity (V) or Acceleration (A)
    #    label     [str][opt.]  - Label the dataset for plotting etc.
    #    yspec     [str][opt.]  - Required for YSPEC data - general prefix in output files
    # OUTPUTS:
    #    st        [Obspy st.]  - Obspy stream object with all data
    # ==================================================================================================================

    st = obspy.Stream() # Create stream for output

    channel_loop_no = range(meta.no_chnls)


    for stn in meta.stn_list:
        # Loop through each channel

        for i in channel_loop_no:
            chl = meta.channels[i] # Get channel letter

            if chl == "G":
                gravity=True
            else:
                gravity=False


            if meta.custom_networks:
                network = meta.network_list[stn]
            else:
                network = meta.network


            # Load the required data:
            data = np.loadtxt(get_loadstr(type, filepath, stn, measure, chl, gravity, yspec, network=network)) # format: [t, trace]

            # Subtract for gravity:
            if meta.gravity_subtract==True:
                # Load the non-gravity data:
                nograv = np.loadtxt(get_loadstr(type, nograv_fp, stn, measure, chl, gravity, grav_yspec, network=meta.network))
                d = data[:,1] - nograv[:,1]

            else:
                d = data[:,1]

            tr = _set_tr_stats(d, type=type, stn=stn, chl=chl, dt=data[1,0]-data[0,0], label=label,
                               network=network, attach_coords=meta.attach_coords, meta=meta)

            st += tr # Add trace to stream
    print(f"Loaded {type} data")
    return st




def _set_tr_stats(data, type, stn, chl, dt, label, network, meta, attach_coords=True):
    # ==================================================================================================================
    # DESCRIPTION:
    # Sets the metadata or 'stats' for the Obspy traces
    # INPUTS:
    #    data   [1D arr.]          - Trace data
    #    type   [str]              - type of simulation data - e.g. "nmsyn"
    #    stn    [str]              - station e.g. "X151"
    #    chl    [str]              - Channel - e.g. N, E, Z or G
    #    dt     [float]            - timestep
    #    label  [str]              - label for simulation data
    # OUTPUTS:
    #    tr     [Obspy trace]      - Trace with metadata
    # ==================================================================================================================

    # Set data and simulation_metadata
    tr = obspy.Trace()


    tr.stats.network = network

    tr.stats.station = stn
    tr.stats.channel = chl
    tr.stats.delta = dt

    tr.stats["type"] = type # Add simulation type to trace dict

    if label == None:
        tr.stats["label"] = type # Add simulation type to trace dict
    else:
        tr.stats["label"] = label


    tr.stats["coordinates"] = {}  # add the coordinates to your dictionary,
    tr.stats["coordinates"]["latitude"] = None
    tr.stats["coordinates"]["longitude"] = None

    if attach_coords == True:
        if meta.Xstations:
            tr.stats["coordinates"]["latitude"]  = float(stn[1:]) - 91
            tr.stats["coordinates"]["longitude"] = 0.0
        else:
            tr.stats["coordinates"]["latitude"]  = float(meta.lat_list[stn])
            tr.stats["coordinates"]["longitude"] = float(meta.lon_list[stn])




    ## Adding data:
    tr.data = data

    return tr




def get_loadstr(type, filepath, stn, measure, chl, gravity, yspec, network="Y5"):
    # ==================================================================================================================
    # DESCRIPTION:
    # Generates string with the filepath and file name for loading of data
    # INPUTS:
    #    type     [str]    - type of simulation data - e.g. "nmsyn"
    #    filepath [str]    - Directory filepath holding the dataset
    #    stn      [str]    - station e.g. "X151"
    #    measure  [str]    - Displacement (D), Velocity (V) or Acceleration (A)
    #    chl      [str]    - Channel - e.g. N, E, Z or G
    #    gravity  [bool]   - True or False for whether loading gravity data
    #    yspec    [str]    - Feneral prefix in output files
    #    network  [str]    - Network for stations
    # OUTPUTS:
    #    load_str [str]    - String holding file path and name
    # ==================================================================================================================



    # Strings require alteration for gravity style data
    if gravity==True:
        spfmx_chl = "MXG"
        nmsyn_chl = "GRV"
    else:
        #spfmx_chl = f"MX{chl}{measure}"
        spfmx_chl = f"MX{chl}{measure}" #remove the D at the end
        nmsyn_chl = f"LH{chl}{measure}"

    # Generate load string:
    # YOU MAY NEED TO EDIT THESE FOR YOUR DATA - e.g. removing the conv_
    t = type.upper()
    if   t == "YSPEC":
        # Example: conv_HN_test1_run5.yspec.X1.DN
        # load_str = f"{filepath}/conv_HN_test1_run{run}.yspec.{stn}.{measure}{chl}"
        load_str = f"{filepath}/{yspec}.{stn}.{measure}{chl}"

    elif t == "SPFMX":
        # Example: conv_X1.Y5.MXEV.sem.ascii
        #load_str = f"{filepath}/conv_{network}.{stn}.{spfmx_chl}.sem.ascii"
        if yspec==0:
            load_str = f"{filepath}/conv_{network}.{stn}.{spfmx_chl}.sem.ascii"
        else:
            load_str = f"{filepath}/conv_{stn}.{network}.{spfmx_chl}.sem.ascii"

    elif t == "NMSYN":
        # Example: X21.Y5.LHE.nmsyn
        load_str = f"{filepath}/conv_{stn}.{network}.{nmsyn_chl}.nmsyn"
        #load_str = f"{filepath}/{stn}.{network}.{nmsyn_chl}.nmsyn"

    elif t == "AXISEM":
        # Example: Y5.X2.N.ascii
        load_str = f"{filepath}/conv_{network}.{stn}.ascii.D{chl}"
    else:
        raise ValueError("Must be spfmx/yspec/nmsyng...")

    print('loading :', load_str)

    return load_str


def process_stream(st_stream, meta, reverse_P=False, reverse_all=False, reverse_G=False, reverse_Z=False, reverse_T=False):
    # ==================================================================================================================
    # DESCRIPTION:
    # Filters and slices streams. Also sets early parts of some traces to 0 to remove convolution errors.
    # INPUTS:
    #    st_stream   [Obspy Stream] - Stream for processing
    #    meta        [metadata]     - PlotMetadata object (see plot_meta_data_class.py)
    # OUTPUTS:
    #    stream      [Obspy stream] - Processed stream
    # ==================================================================================================================

    # Select channel, filter, slice
    stream = st_stream
    stream = stream.filter(type='bandpass', freqmin=meta.fmin, freqmax=meta.fmax)
    stream = stream.slice(starttime=stream[0].stats.starttime + meta.tmin,
                          endtime=stream[0].stats.starttime + meta.tmax)

    for i in range(len(stream)):
        tr = stream[i]


        tr.stats.starttime = tr.stats.starttime - meta.time_offset[tr.stats.type]


        if np.logical_and(reverse_P, tr.stats.channel == "P"):
            tr.data = tr.data*(-1)
        if np.logical_and(reverse_T, tr.stats.channel == "T"):
            tr.data = tr.data * (-1)
        if np.logical_and(reverse_Z, tr.stats.channel == "Z"):
            tr.data = tr.data*(-1)
        if np.logical_and(reverse_G,  tr.stats.channel == "G"):
            tr.data = tr.data*(-1)
        if reverse_all:
                tr.data = tr.data * (-1)
    return stream