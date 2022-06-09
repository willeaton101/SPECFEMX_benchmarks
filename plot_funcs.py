# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Functions for plotting data from YSPEC, SPECFEMX, AxiSEM3D or NMSYN
# ______________________________________________________________________________________________________________________
import matplotlib.pyplot as plt
import numpy as np

LW = 1.5    # Global linewidth

def normalise(array):
    norm = array/np.amax(np.abs(array))
    return norm

def plot_station(st, station, chls, meta, style='dark'):
    # ==================================================================================================================
    # DESCRIPTION:
    # Plots single station for multiple different simulation types (e.g. nmsyng, yspec)
    # INPUTS:
    #    st      [Obspy stream]     - Stream holding all traces
    #    station [str]              - station e.g. "X151"
    #    chls    [list]             - channels to be plotted
    #    meta    [metadata]         - metadata for plotting - see plot_meta_data_class.py
    # OUTPUTS:
    #    Returns matplotlib figure for that station
    # ==================================================================================================================
    # Set dark mode:
    if style == 'dark':
        plt.style.use('dark_background')

    n = len(chls)                                                               # number of channels to loop through
    fig, ax = plt.subplots(n, 1, figsize=(12, 7.5), sharex=True, sharey=True)   # create figure.

    # Loop through each channel:
    for k in range(n):
        st_plot = st.select(station=station, channel=chls[k])
        leg_lines = []
        leg_string = []

        # Check if there is a line to plot for each type
        for type in meta.types:
            # If there is then for each line of that type, plot it
            for i in range(len(st_plot)):
                tr = st_plot[i]
                if tr.stats.type == type:
                    time = np.linspace(0, tr.stats.npts * tr.stats.delta, tr.stats.npts) + meta.tmin

                    # Adjust start time for specific traces:
                    time = time - meta.time_offset[tr.stats.type]

                    offset = meta.calc_offset(stn_lat=tr.stats.coordinates.latitude, stn_lon=tr.stats.coordinates.longitude)

                    c = meta.colours[tr.stats.label]
                    if n==1:
                        line, = ax.plot(time, normalise(tr.data), '-', linewidth=LW, color=c)
                    else:
                        if tr.stats.type=="spfmx":
                            line, = ax[k].plot(time, normalise(tr.data), ':', linewidth=LW, color=c)
                        else:
                            line, = ax[k].plot(time, normalise(tr.data), '-', linewidth=LW, color=c)

                    leg_lines.append(line)
                    leg_string.append(str(tr.stats.label))

        # Add bells and whistles:
        plt.suptitle(f"Station {station}: Latitude = {tr.stats.coordinates.latitude} ; Δ = {np.round(offset, 2)} degrees ; {meta.fmin_mHz}-{meta.fmax_mHz} mHz")
        title_str = f"Channel: {meta.chl_options[chls[k]]}"
        x_label   = "Time [s]"
        if tr.stats.channel == "G":
            y_label   = "Norm. Gravitational potential perturbation"
        else:
            y_label   = "Norm. displacement"

        if n==1:
            ax.set_title(title_str)
            ax.legend(leg_lines, leg_string)
            ax.set_xlabel(x_label)
            #ax.set_ylim([-1, 1])
            ax.set_xlim([meta.tmin, meta.tmax])
            ax.set_ylabel(y_label)
        else:
            ax[k].set_title(title_str)
            ax[k].set_xlim([meta.tmin, meta.tmax])
            ax[0].legend(leg_lines, leg_string)
            ax[-1].set_xlabel(x_label)
            #ax[k].set_ylim([-1, 1])
            ax[k].set_ylabel(y_label)

    fig.set_tight_layout(True)

    return fig



def plot_record(stream, meta, amplitude=10, orientation="HORIZONTAL", style='dark'):
    # ==================================================================================================================
    # DESCRIPTION:
    # Plots record section multiple different simulation types (e.g. nmsyng, yspec)
    # INPUTS:
    #    stream         [Obspy stream]  - Stream holding all traces
    #    meta           [metadata]      - Metadata for plot
    #    amplitude      [float]         - Amplitude of each trace in terms of azimuthal distance (y axis)
    #    orientation    [str]           - Can be horizontal or vertical
    # OUTPUTS:
    #    Returns Matplotlib figure
    # ==================================================================================================================
    if style == 'dark':
        plt.style.use('dark_background')

    # Counters/lists needed for flexible legend
    spfmx_leg_ctr = 0
    yspec_leg_ctr = 0
    nmsyng_leg_ctr = 0
    axi_leg_ctr = 0
    leg_list = []
    leg_names = []
    offsets = []

    fig_record, ax_r = plt.subplots(figsize=(9,7))

    for tr in stream:
        # Get type of simulation for data
        type = tr.stats.type

        # Get time and data of trace
        time = np.linspace(0, tr.stats.delta * tr.stats.npts, tr.stats.npts)
        if type == 'spfmx':
            time = time - 1.79550004

        data = amplitude * normalise(tr.data) # Normalised data

        # Calculate offset and store copy for ylims
        offset = meta.calc_offset(stn_lat=tr.stats.coordinates.latitude, stn_lon=tr.stats.coordinates.longitude)
        offsets.append(offset)

        # Sort out Legend colours
        c = meta.colours[tr.stats.label]

        # Prepare plotting based on orientation:
        if orientation.upper()=="VERTICAL":
            x = offset + data
            y = time
        else:
            x = time
            y = offset + data

        # loop required because of issues with Fig legend.
        if type == 'spfmx':
            if tr.stats.label=='APRIL17':
                line_spfmx, = ax_r.plot(x, y, ':', color=c, linewidth=LW)
            else:
                line_spfmx, = ax_r.plot(x, y, ':', color=c, linewidth=LW)
            spfmx_leg_ctr, leg_list, leg_names = _check_ctr(spfmx_leg_ctr, line_spfmx, "SPECFEM-X", leg_list, leg_names)
        elif type == 'yspec':
            line_yspec, = ax_r.plot(x, y, color=c, linewidth=LW)
            yspec_leg_ctr, leg_list, leg_names = _check_ctr(yspec_leg_ctr, line_yspec, "YSPEC", leg_list, leg_names)
        elif type == 'nmsyn':
            line_nmsyn, = ax_r.plot(x, y, color=c, linewidth=LW)
            nmsyng_leg_ctr, leg_list, leg_names = _check_ctr(nmsyng_leg_ctr, line_nmsyn, "NMSYNG", leg_list, leg_names)
        elif type == 'axisem':
            line_axi, = ax_r.plot(x, y, color=c, linewidth=LW)
            axi_leg_ctr, leg_list, leg_names = _check_ctr(axi_leg_ctr, line_axi, "AXISEM", leg_list, leg_names)

    # Set figure metadetails
    if orientation.upper() == "VERTICAL":
        x_low = np.min(np.array(offsets)) - amplitude
        x_high = np.max(np.array(offsets)) + amplitude
        xlabel = "Distance, Δ [degrees]"

        y_low = meta.tmin
        y_high = meta.tmax
        ylabel = "Time [s]"
    else:
        x_low = meta.tmin
        x_high = meta.tmax
        xlabel = "Time [s]"

        y_low = np.amin(np.array(offsets)) - amplitude
        y_high = np.amax(np.array(offsets)) + amplitude
        #ylabel = "Distance, Δ [degrees]"
        ylabel = "Station Latitude [degrees]"

    fs = 12
    ax_r.set_xlabel(xlabel, fontsize=fs)
    ax_r.set_ylabel(ylabel, fontsize=fs)
    ax_r.set_xlim([x_low, x_high])
    ax_r.set_ylim([y_low, y_high])

    ax_r.set_title(f"Record Section: Channel {stream[0].stats.channel}")
    #ax_r.legend(leg_list, leg_names, loc="upper left")
    fig_record.set_tight_layout(True)

    return fig_record

def _check_ctr(counter, line, name, leg_lines, leg_names):
    # ==================================================================================================================
    # DESCRIPTION:
    # Internal function for plot_record
    # INPUTS:
    #    counter   [int]    - Counter determining if this if statement has been triggered before
    #    line      [dict]   - Matplotlib 2D line object that will be used in legend
    #    name      [float]  - Associated name for this line object in the legend
    #    leg_lines [dict]   - List holding all 2D lines for legend. Line will be added to this.
    #    leg_names [float]  - List holding all name strings for legend. Name will be added to this.
    # OUTPUTS:
    #    Returns updated versions of counter, leg_lines and leg_names
    # ==================================================================================================================

    if counter == 0:
        leg_lines.append(line)
        leg_names.append(name)
        counter += 1

    return counter, leg_lines, leg_names



def plot_spectra(st, station, meta, chls=["Z", "T", "P"], style='dark'):
    # ==================================================================================================================
    # DESCRIPTION:
    # Plots single station for multiple different simulation types (e.g. nmsyng, yspec)
    # INPUTS:
    #    st      [Obspy stream]     - Stream holding all traces
    #    station [str]              - station e.g. "X151"
    #    meta    [dict]             - Metadata for plot
    #    chls    [list]             - list of channels to plot
    # OUTPUTS:
    #    Returns matplotlib figure
    # ==================================================================================================================
    if style=='dark':
        plt.style.use('dark_background')

    st_stn = st.select(station=station)
    N = len(chls)

    fig, ax = plt.subplots(N, 1, sharex=True, figsize=(12, 7.5))

    # Loop through channels:
    for k in range(N):
        chl = chls[k]
        leg_lines = []
        leg_string = []

        for type in meta.types:
            st_plot = st_stn.select(channel=chl)

            for i in range(len(st_plot)):
                tr = st_plot[i]
                offset = meta.calc_offset(stn_lat=tr.stats.coordinates.latitude, stn_lon=tr.stats.coordinates.longitude)

                if tr.stats.type == type:

                    data = tr.data

                    freq  = np.fft.fftfreq(n=tr.stats.npts, d=tr.stats.delta)

                    # Avoid freq of 0
                    mask  = [freq>0]
                    period = np.reciprocal(freq[tuple(mask)])
                    power = np.abs(np.fft.fft(data))[tuple(mask)]


                    label = tr.stats.label
                    c = meta.colours[label]

                    line, = ax[k].plot(period, normalise(power), linewidth=LW, color=c)
                    leg_lines.append(line)
                    leg_string.append(str(type))

        plt.suptitle(f"Station {station}: Latitude = {tr.stats.coordinates.latitude} ; Δ = {np.round(offset, 2)} degrees ; {meta.fmin_mHz}-{meta.fmax_mHz} mHz")
        ax[k].axvline(x=meta.period_max, ymin=0, ymax=1, linestyle="--")
        ax[k].axvline(x=meta.period_min, ymin=0, ymax=1, linestyle="--")
        ax[k].set_title(f"Channel: {chl}")
        ax[k].legend(leg_lines, leg_string)
        ax[-1].set_xlabel("Period [s]")
        ax[k].set_xscale("log")

        freqbuffer = 0.05*(meta.period_max - meta.period_min )

        ax[k].set_xlim([meta.period_min-freqbuffer, meta.period_max+freqbuffer])
        ax[k].set_ylabel("Power")

    return fig 