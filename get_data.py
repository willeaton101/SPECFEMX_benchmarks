# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Functions for loading and processing data from YSPEC, SPECFEMX, AxiSEM3D or NMSYN
# ______________________________________________________________________________________________________________________
import os

import matplotlib.pyplot as plt
import numpy as np
import obspy


def get_data(type, filepath, meta, measure="D", label=None, yspec="", grav_yspec="", nograv_fp="", rename=None):
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
            if label=='REAL':
                d = obspy.read(get_loadstr(type, filepath, stn, measure, chl, gravity, yspec, network=network))[0]
            else:

                data = np.loadtxt(get_loadstr(type, filepath, stn, measure, chl, gravity, yspec, network=network)) # format: [t, trace]
                d = data[:,1]


            if label=='REAL':
                tr = _set_real_stats(d, stn=stn, attach_coords=meta.attach_coords, type=type, meta=meta, rename=rename)
            else:
                tr = _set_tr_stats(d, type=type, stn=stn, chl=chl, dt=data[1, 0] - data[0, 0], label=label,
                                   network=network, attach_coords=meta.attach_coords, meta=meta, rename=rename)

            st += tr # Add trace to stream
    print(f"Loaded {type} data")
    return st


def _set_real_stats(tr, stn, attach_coords, meta, type, label=None, rename=None):


    chl = tr.stats.channel

    if rename!=None:
        try:
            chl = rename[chl]
        except:
            'Didnt rename channel: ', chl

    if chl[:2]=='BH':
        tr.stats.channel = chl[2]



    tr.stats["type"] = 'real'
    tr.stats["coordinates"] = {}  # add the coordinates to your dictionary,
    tr.stats["coordinates"]["latitude"] = None
    tr.stats["coordinates"]["longitude"] = None

    if label == None:
        tr.stats["label"] = type # Add simulation type to trace dict
    else:
        tr.stats["label"] = label

    if attach_coords == True:
        if meta.Xstations:
            tr.stats["coordinates"]["latitude"]  = float(stn[1:]) - 91
            tr.stats["coordinates"]["longitude"] = 0.0
        else:
            tr.stats["coordinates"]["latitude"]  = float(meta.lat_list[stn])
            tr.stats["coordinates"]["longitude"] = float(meta.lon_list[stn])

    return tr


def _set_tr_stats(data, type, stn, chl, dt, label, network, meta, attach_coords=True, rename=None):
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


    if rename!=None:
        try:
            chl = rename[chl]
            print('Renamed channel to ', chl)

        except:
            print('Didnt rename channel: ', chl)


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

    if chl[1:3] == 'C.':
        if type!='nmsyn':
            measure= chl[2:]
            chl= chl[:2]

    if 'GRAVIMETER' in chl:
        measure = ''
        chl = '.'+chl

    if 'STRAIN' in chl:
        measure= ''
        #chl = chl[:2]


    # Strings require alteration for gravity style data
    if gravity==True:
        spfmx_chl = "MXG"
        nmsyn_chl = "GRV"
    else:
        #spfmx_chl = f"MX{chl}{measure}"
        spfmx_chl = f"MX{chl}{measure}" #remove the D at the end
        nmsyn_chl = f"LH{chl}"#{measure}"

        globe_chl = f"MX{chl}"

    # Generate load string:
    # YOU MAY NEED TO EDIT THESE FOR YOUR DATA - e.g. removing the conv_
    t = type.upper()
    if   t == "YSPEC":
        # Example: conv_HN_test1_run5.yspec.X1.DN
        # load_str = f"{filepath}/conv_HN_test1_run{run}.yspec.{stn}.{measure}{chl}"
        load_str = f"{filepath}/{yspec}.{stn}.{measure}{chl}"

    elif t == "SPFMX":

        if network=="SY":
            network = "Y5"



        # Example: conv_X1.Y5.MXEV.sem.ascii
        #load_str = f"{filepath}/conv_{network}.{stn}.{spfmx_chl}.sem.ascii"
        if yspec==0:
            load_str = f"{filepath}/conv_{network}.{stn}.{spfmx_chl}.sem.ascii"
        elif yspec==2:
            load_str = f"{filepath}/{stn}.{network}.{spfmx_chl}.sem.ascii"
        else:
            load_str = f"{filepath}/conv_{stn}.{network}.{spfmx_chl}.sem.ascii"

    elif t == "3DGLOBE":

        if measure!='D':
            raise Warning("Specfem3D globe automatically outputs disp only. Assuming displacement.")

        if network=="SY":
            network = "Y5"

        # Example: conv_X1.Y5.MXEV.sem.ascii
        #load_str = f"{filepath}/conv_{network}.{stn}.{spfmx_chl}.sem.ascii"
        if yspec==0:
            load_str = f"{filepath}/conv_{network}.{stn}.{globe_chl}.sem.ascii"
        elif yspec==2:
            load_str = f"{filepath}/{network}.{stn}.{globe_chl}.sem.ascii"
        else:
            load_str = f"{filepath}/conv_{stn}.{network}.{spfmx_chl}.sem.ascii"


    elif t == "NMSYN":
        # Example: X21.Y5.LHE.nmsyn
        if yspec==2:
            load_str = f"{filepath}/{stn}.{network}.{nmsyn_chl}.nmsyn"
        else:
            load_str = f"{filepath}/conv_{stn}.{network}.{nmsyn_chl}.nmsyn"
        #load_str = f"{filepath}/{stn}.{network}.{nmsyn_chl}.nmsyn"

    elif t=="REAL":
        # Example: IU.INCN.BHZ.D.mseed
        load_str = f"{filepath}/{network}.{stn}.BH{chl}.{measure}.mseed"

    elif t == "AXISEM":
        # Example: Y5.X2.N.ascii
        load_str = f"{filepath}/conv_{network}.{stn}.ascii.D{chl}"
    elif t == "QSSP":
        # Example: label_X11.Y5.DN
        if network=="SY":
            network = "Y5"
        elif network =='GG':
            network = 'G'
        load_str = f"{filepath}/{yspec}_{stn}.{network}.{measure}{chl}"
    else:
        raise ValueError("Must be spfmx/yspec/nmsyng...")

    print('loading :', load_str)

    return load_str


def process_stream(st_stream, meta, fixyspecamp=False, reverse=[], slice_by_picks=None,
                   filters=None):
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


    if filters==None:
        stream = stream.filter(type='bandpass', freqmin=meta.fmin, freqmax=meta.fmax)
    else:
        for IF in filters:
            stream = stream.filter(type=IF[0], freq=IF[1], corners=IF[2], zerophase=IF[3])



    for i in range(len(stream)):

        endtime = stream[0].stats.starttime + meta.tmax + meta.tmin
        if slice_by_picks != None:
            endtime = stream[0].stats.starttime + slice_by_picks[stream[i].stats.station]
            print(slice_by_picks[stream[i].stats.station])

        stream[i] = stream[i].slice(starttime=stream[0].stats.starttime + meta.tmin, endtime=endtime)


        if np.logical_and(fixyspecamp, stream[i].stats.label=='YSPEC'):
            stream[i].data = stream[i].data * 1e-7



        #stream[i].stats.starttime = stream[i].stats.starttime - meta.time_offset[stream[i].stats.type]


        for c in reverse:
            if c==stream[i].stats.channel:
                stream[i].data = stream[i].data * (-1)
        '''if np.logical_and(reverse_P, tr.stats.channel == "P"):
            
        if np.logical_and(reverse_T, tr.stats.channel == "T"):
            tr.data = tr.data * (-1)
        if np.logical_and(reverse_Z, tr.stats.channel == "Z"):
            tr.data = tr.data*(-1)
        if np.logical_and(reverse_G, tr.stats.channel == "G"):
            tr.data = tr.data*(-1)
        if reverse_all:
                tr.data = tr.data * (-1)'''
    return stream




def calc_internal_gravity(st, chl, rho, meta):
    Gconstant = 6.6723e-11

    # Check that we have both channels present:
    ZZ = st.select(channel=f'{chl}')
    PP = st.select(channel=f'{chl}C.PGRAV')

    # Loop stations:
    assert (len(ZZ)==len(PP))
    for istn in range(len(ZZ)):

        Z = ZZ.select(station=meta.stn_list[istn])
        P = PP.select(station=meta.stn_list[istn])

        if np.logical_or(len(Z)!=1, len(P)!=1):
            raise ValueError(f'Error: Wrong amount of channel data in stream: Z has {len(Z)}, P has {len(P)} traces')

        from copy import copy
        from obspy.core.trace import Trace
        tr = Trace()

        tr.data = P[0].data - 4 * np.pi * Gconstant * rho * Z[0].data

        # Add stats:
        tr.stats.type           = P[0].stats.type
        tr.stats.sampling_rate  = copy(P[0].stats.sampling_rate)
        tr.stats.coordinates    = copy(P[0].stats.coordinates)
        tr.stats.network        = copy(P[0].stats.network)
        tr.stats.label          = copy(P[0].stats.label)
        tr.stats.station        = copy(P[0].stats.station)
        tr.stats.channel        = f'{chl}C.PGRAV_INSIDE'

        st += tr
    return st



def calc_external_gravity(st, chl, rho, meta):
    Gconstant = 6.6723e-11

    # Check that we have both channels present:
    ZZ         = st.select(channel=f'{chl}')
    PP         = st.select(channel=f'{chl}C.PGRAV_INSIDE')

    # Loop stations:
    assert (len(ZZ)==len(PP))
    for istn in range(len(ZZ)):

        Z = ZZ.select(station=meta.stn_list[istn])
        P = PP.select(station=meta.stn_list[istn])

        if np.logical_or(len(Z)!=1, len(P)!=1):
            raise ValueError(f'Error: Wrong amount of channel data in stream: Z has {len(Z)}, P has {len(P)} traces')

        from copy import copy
        from obspy.core.trace import Trace
        tr = Trace()

        tr.data = P[0].data + 4 * np.pi * Gconstant * rho * Z[0].data

        # Add stats:
        tr.stats.type           = P[0].stats.type
        tr.stats.sampling_rate  = copy(P[0].stats.sampling_rate)
        tr.stats.coordinates    = copy(P[0].stats.coordinates)
        tr.stats.network        = copy(P[0].stats.network)
        tr.stats.label          = copy(P[0].stats.label)
        tr.stats.station        = copy(P[0].stats.station)
        tr.stats.channel        = f'{chl}C.PGRAV'

        st += tr
    return st


def calculate_QSSP_free_air(st, chl, meta):
    Gconstant = 6.6723e-11
    print('Warning: QSSP Free air - ensure using vertical ACCELERATION not DISPLACEMENT ')

    # Check that we have both channels present:
    ZZ = st.select(channel=f'{chl}')
    PP = st.select(channel=f'{chl}C.PGRAV_INSIDE')
    GRAVI = st.select(channel=f'GRAVIMETER{chl}')


    # Loop stations:
    assert (len(ZZ)==len(PP)==len(GRAVI))
    for istn in range(len(ZZ)):

        Z = ZZ.select(station=meta.stn_list[istn])
        P = PP.select(station=meta.stn_list[istn])
        G = GRAVI.select(station=meta.stn_list[istn])

        if np.logical_or(len(Z)!=1, np.logical_or(len(P)!=1, len(G)!=1) ):
            raise ValueError(f'Error: Wrong amount of channel data in stream: Z has {len(Z)}, P has {len(P)}, GRAVIMETER has {len(GRAVI)} traces')

        from copy import copy
        from obspy.core.trace import Trace
        tr = Trace()

        # gravimeter = acc + change in grav + free air
        tr.data = G[0].data - Z[0].data + P[0].data


        # Add stats:
        tr.stats.type           = P[0].stats.type
        tr.stats.sampling_rate  = copy(P[0].stats.sampling_rate)
        tr.stats.coordinates    = copy(P[0].stats.coordinates)
        tr.stats.network        = copy(P[0].stats.network)
        tr.stats.label          = copy(P[0].stats.label)
        tr.stats.station        = copy(P[0].stats.station)
        tr.stats.channel        = f'{chl}C.FA'


        # Mask the aliasing:
        """from wetools import obspy_gen_mpl
        t,y = obspy_gen_mpl(tr)
        tmask = t<0
        y[tmask] = 0
        tr.data = y"""

        st += tr

    return st



def calculate_free_air_from_displacement(st, chl, rho, meta):
    Gconstant = 6.6723e-11
    print('Warning: Free air from displacement - ensure using vertical DISPLACEMENT ')

    # Check that we have both channels present:
    ZZ = st.select(channel=f'{chl}')

    assert(chl=='Z')

    # Loop stations:
    for istn in range(len(ZZ)):

        Z = ZZ.select(station=meta.stn_list[istn])

        if len(Z)!=1:
            raise ValueError(f'Error: Wrong amount of channel data in stream: Z has {len(Z)} traces')

        from copy import copy
        from obspy.core.trace import Trace
        tr = Trace()

        # gravimeter = acc + change in grav + free air
        tr.data = Z[0].data * (4*np.pi*Gconstant*rho - 2*9.81/6371000)
        print((4*np.pi*Gconstant*rho - 2*9.81/6371000))


        # Add stats:
        tr.stats.type           = Z[0].stats.type
        tr.stats.sampling_rate  = copy(Z[0].stats.sampling_rate)
        tr.stats.coordinates    = copy(Z[0].stats.coordinates)
        tr.stats.network        = copy(Z[0].stats.network)
        tr.stats.label          = copy(Z[0].stats.label)
        tr.stats.station        = copy(Z[0].stats.station)
        tr.stats.channel        = f'{chl}C.FA'

        st += tr

    return st



