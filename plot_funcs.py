# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Functions for plotting data from YSPEC, SPECFEMX, AxiSEM3D or NMSYN
# ______________________________________________________________________________________________________________________
import matplotlib.pyplot as plt
import numpy as np
from wetools import *

LW = 1.3   # Global linewidth

def normalise(array):
    norm = array/np.amax(np.abs(array))
    return norm


def plot_station(st, station, chls, meta, style='dark', norm_amp=False, min_norm_time=0,
                 return_ax=False, plot_picks=False, custom_ylim=None):
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
    fig, ax = plt.subplots(n, 1, figsize=(12, 7.5), sharex=True, sharey=False)   # create figure.


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
                    time = time + meta.time_offset[tr.stats.type]

                    # Get offset of station from src:
                    offset = meta.calc_offset(stn_lat=tr.stats.coordinates.latitude,
                                              stn_lon=tr.stats.coordinates.longitude)
                    # Get axis for plotting:
                    if n==1:
                        axtmp = ax
                    else:
                        axtmp = ax[k]

                    # Get line colour:
                    c = meta.colours[tr.stats.label]

                    # Amplitude:
                    if norm_amp:
                        t, x = obspy_gen_mpl(tr)
                        dtmp = normalise_after_time(t, tr.data, min_norm_time)
                    else:
                        dtmp = tr.data

                    # Plot!

                    if plot_picks:
                        axtmp.axvline(float(meta.picks[tr.stats.station][0]) , alpha=0.5, color='k')





                    line, = axtmp.plot(time, dtmp, '-', linewidth=LW, color=c)
                    leg_lines.append(line)
                    leg_string.append(str(tr.stats.label))

        # Add bells and whistles:
        plt.suptitle(f"Station {station}: Latitude = {tr.stats.coordinates.latitude} ; Δ = {np.round(offset, 2)} degrees ; {meta.fmin_mHz}-{meta.fmax_mHz} mHz")
        title_str = f"Channel: {meta.chl_options[chls[k]]}"
        x_label   = "Time [s]"
        #if tr.stats.channel == "G":
        #   # y_label   = "Gravitational potential perturbation"
        #else:
           # y_label   = "Displacement"

        if norm_amp:
            #y_label = 'Normalised ' + y_label
            axtmp.set_ylim([-1, 1])

        if custom_ylim!=None:
            axtmp.set_ylim(custom_ylim)
        # Bells and whistles:
        axtmp.set_title(title_str)
        axtmp.legend(leg_lines, leg_string, loc="upper right")
        axtmp.set_xlabel(x_label)


        axtmp.set_xlim([meta.tmin, meta.tmax])
        #axtmp.set_xlim([meta.tmin, meta.picks[station]])


        #axtmp.set_ylabel(y_label)


    fig.set_tight_layout(True)

    if return_ax:
        return fig, ax
    else:
        return fig



def plot_record(stream, meta, amplitude=10, orientation="HORIZONTAL", style='dark',
                min_norm_time=0, custom_norm=None):
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
    globe_leg_ctr = 0

    yspec_leg_ctr = 0
    nmsyng_leg_ctr = 0
    qssp_leg_ctr = 0
    axi_leg_ctr = 0
    real_leg_ctr = 0
    leg_list = []
    leg_names = []
    offsets = []

    fig_record, ax_r = plt.subplots(figsize=(9,7))

    for tr in stream:
        # Get type of simulation for data
        type = tr.stats.type
        # Get time and data of trace
        time = np.linspace(0, tr.stats.delta * tr.stats.npts, tr.stats.npts)
        # Apply the offset
        time = time + meta.time_offset[type]

        if custom_norm!=None:
            data = amplitude * tr.data/custom_norm
        else:
            data = amplitude * normalise_after_time(time, tr.data, min_norm_time) # Normalised data



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
            line_spfmx, = ax_r.plot(x, y, '-', color=c, linewidth=LW)
            spfmx_leg_ctr, leg_list, leg_names = _check_ctr(spfmx_leg_ctr, line_spfmx, "SPECFEM-X", leg_list, leg_names)
        if type.upper() == '3DGLOBE':
            line_globe, = ax_r.plot(x, y, '-', color=c, linewidth=LW)
            globe_leg_ctr, leg_list, leg_names = _check_ctr(globe_leg_ctr, line_globe, "SPECFEM3D Globe", leg_list, leg_names)
        elif type == 'yspec':
            line_yspec, = ax_r.plot(x, y, color=c, linewidth=LW)
            yspec_leg_ctr, leg_list, leg_names = _check_ctr(yspec_leg_ctr, line_yspec, "YSPEC", leg_list, leg_names)
        elif type == 'nmsyn':
            line_nmsyn, = ax_r.plot(x, y, color=c, linewidth=LW)
            nmsyng_leg_ctr, leg_list, leg_names = _check_ctr(nmsyng_leg_ctr, line_nmsyn, "NMSYNG", leg_list, leg_names)
        elif type == 'axisem':
            line_axi, = ax_r.plot(x, y, color=c, linewidth=LW)
            axi_leg_ctr, leg_list, leg_names = _check_ctr(axi_leg_ctr, line_axi, "AXISEM", leg_list, leg_names)
        elif type == 'qssp':
            line_axi, = ax_r.plot(x, y, color=c, linewidth=LW)
            qssp_leg_ctr, leg_list, leg_names = _check_ctr(qssp_leg_ctr, line_axi, "QSSP", leg_list, leg_names)
        elif type == 'real':
            line_real, = ax_r.plot(x, y, color=c, linewidth=LW)
            real_leg_ctr, leg_list, leg_names = _check_ctr(real_leg_ctr, line_real, "REAL", leg_list, leg_names)
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
        label = "Distance, Δ [degrees]"
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



def plot_spectra(st, station, meta, chls=["Z", "T", "P"], style='light'):
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

                    line, = ax[k].plot(period, power, linewidth=LW, color=c)
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



def plot_synthetic_L2misfit_OLD(d, stn, label1, label2, meta, resample=None, mismax=None, normalise_seis=False, slicer=None):
    fig, ax = plt.subplots(len(meta.plot_channels)+2, figsize=(12, 7.5), sharex=True)   # create figure.
    fig.set_tight_layout(True)

    maxvals = []
    stndata = d.select(station=stn)

    plt.suptitle(f"Station {stn}: Lat. = {stndata[0].stats.coordinates.latitude}, Lon. = {stndata[0].stats.coordinates.longitude}; {meta.fmin_mHz}-{meta.fmax_mHz} mHz")

    # Initialise L2 misfit time-series
    L2norm = np.zeros(len(stndata[0].data))

    chlctr = 0
    for chl in meta.plot_channels:
        chldata  = stndata.select(channel=chl)

        found1 = False
        found2 = False
        for i in chldata.traces:
            if i.stats.label == label1:
                l1 = i
                found1 = True
            elif i.stats.label == label2:
                l2 = i
                found2 = True

        # Check both are found:
        if np.logical_and(found1, found2):
           # Found both
            pass
        else:
            raise ValueError(f"Error: Found: Label {label1} {found1}  - Label {label2} {found2} ")

        if slicer!=None:
            for slice in slicer:
                slicerlabel = slice[0]
                slicer_st = slice[1] # start time
                slicer_et = slice[2] # end time:


                if l1.stats.label==slicerlabel:
                    l1 = l1.slice(starttime=l1.stats.starttime + slicer_st,
                                              endtime=l1.stats.endtime + slicer_et)
                elif l2.stats.label==slicerlabel:
                    l2 = l2.slice(starttime=l2.stats.starttime + slicer_st, endtime=l2.stats.endtime + slicer_et)


        if resample!=None:
            #print(f'Resampling for L2 misfit: samprate {resample} Hz')
            l1  = l1.resample(sampling_rate=resample)
            l2  = l2.resample(sampling_rate=resample)

        l1x, l1y = obspy_gen_mpl(l1)
        l2x, l2y = obspy_gen_mpl(l2)

        # Get shortest length
        len1 = len(l1x)
        len2 = len(l2x)
        lmin = np.min([len1, len2])
        #cutoff
        l1x = l1x[:lmin]
        l2x = l2x[:lmin]
        l1y = l1y[:lmin]
        l2y = l2y[:lmin]

        if normalise_seis == True:
            l1y = normalise(l1y)
            l2y = normalise(l2y)

        L2norm = L2norm[:lmin]

        maxvals.append(np.amax([np.amax(np.abs(l1y)), np.amax(np.abs(l2y))]))

        square_misfit = np.square(l2y - l1y)
        L2norm += square_misfit


        ax[chlctr].plot(l1x, l1y)
        ax[chlctr].plot(l2x, l2y)
        ax[chlctr].set_title(chl)

        #ax[2].plot(l1x, square_misfit)

        chlctr += 1


    maxval = np.amax(np.abs(np.array(maxvals)))

    L2norm = L2norm**0.5
    ax[chlctr+1].plot(l1x, L2norm, 'k')
    ax[chlctr].plot(l1x, L2norm/maxval * 100, 'k')


    mean = np.mean(L2norm/maxval * 100)
    print(mean)
    #normsum = np.sum(L2norm)
    #print("normsum", normsum)
    #ax[-1].text(x=(meta.tmax-meta.tmin)*0.1 + meta.tmin, y=0.9*maxval, s=normsum)


    ax[chlctr+1].set_title('Root sum of squared misfits')
    ax[chlctr].set_title('Normalised to maximum amplitude of seismogram')

    ax[0].legend([label1, label2])


    ax[-1].set_xlabel('Time [s]')
    ax[0].set_xlim([meta.tmin, meta.tmax])
    ax[1].set_xlim([meta.tmin, meta.tmax])
    ax[2].set_xlim([meta.tmin, meta.tmax])
    ax[3].set_xlim([meta.tmin, meta.tmax])


    if mismax !=None:
        ax[-2].set_ylim([0, mismax])
        ax[-1].set_ylim([0, mismax/100])


    return fig# , ax, normsum



def plot_cross_correlations(stream, chls, meta, norm_traces=False, correlations=None):

    need_plot = {}
    no_chls = len(chls)

    fig, ax = plt.subplots(no_chls, figsize=(12, 7.5), sharex=True)
    fig.set_tight_layout(True)

    ax[-1].set_xlim([meta.tmin, meta.tmax])
    # Plot the data; initialise list to hold figures for outputting to user:
    chlctr = 0
    for chl in chls:

        ax[chlctr].set_title(f'Channel = {chl}')
        legend_str = []
        st = stream.select(channel=chl)

        # Generate cross-correlation permutations:
        sims = []
        freqs = []
        for i in range(len(st)):
            sims.append(st[i].stats.label)
            freqs.append(st[i].stats.sampling_rate)

            need_plot[st[i].stats.label] = True


        minfreq = np.min(np.array(freqs))

        # Downsample and slice to smallest array size of the options
        st.resample(sampling_rate=minfreq)
        arrlens = []
        for i in range(len(st)):
            arrlens.append(len(st[i].data))
        minlen = np.min(np.array(arrlens))
        for i in range(len(st)):
            st[i].data = st[i].data[:minlen]

        # Generate permutations:
        lensims = len(sims)
        perms = []
        for j in range(len(sims)):
            for k in range(j + 1, lensims):
                perms.append([sims[j], sims[k]])

        yspecctr = 0
        no_perms = len(perms)
        for m in range(no_perms):
            l1 = perms[m][0]
            l2 = perms[m][1]

            # Assign the two traces:
            for j in range(lensims):

                if st[j].stats.label == 'YSPEC':
                    if yspecctr == 0:
                        amp = 1e-7
                        yspecctr += 1
                    else:
                        amp = 1
                else:
                    amp = 1

                if st[j].stats.label == l1:
                    d1 = st[j]

                    d1.data = d1.data * amp
                elif st[j].stats.label == l2:
                    d2 = st[j]
                    d2.data = d2.data * amp


            d1x, d1y = obspy_gen_mpl(d1)
            d2x, d2y = obspy_gen_mpl(d2)

            if norm_traces:
                d1y = normalise(d1y)
                d2y = normalise(d2y)

            cc = calc_correlation(d1y, d2y)
            ccc = np.around(cc[0] * 100, 4)

            if correlations!=None:
                correlations[chl].append(ccc)

            if need_plot[l1]:
                ax[chlctr].plot(d1x, d1y, color=meta.colours[l1])
                need_plot[l1]=False
                legend_str.append(l1)

            if need_plot[l2]:
                ax[chlctr].plot(d2x, d2y, color=meta.colours[l2])
                need_plot[l2]=False
                legend_str.append(l2)


            if norm_traces:
                ax[chlctr].set_ylim([-1,1])

            ax[chlctr].text(x=10000, y=0.5 + m*0.2 * np.amax(d1y), s=f"{l1} - {l2} :  {ccc} %")
            #ax[m].text(x=10000, y=0.6 * np.amax(d1y), s=f"normalised: {np.around(normcc[0] * 100, 4)} %")
            ax[-1].set_xlabel('Time [s]')

            ax[chlctr].legend(legend_str)

        # Update for each channel for plotting index:
        chlctr+=1

    if correlations != None:
        return fig, correlations
    else:
        return fig



def calc_correlation(y1, y2):
    cc  = np.correlate(y1, y2)
    ac1 = np.correlate(y1, y1)**0.5
    ac2 = np.correlate(y2, y2)**0.5
    normcc = cc/(ac1*ac2)
    return normcc





def plot_synthetic_L2misfit(d, stn, label1, label2, meta, resample=None,
                            mismax=None, normalise_seis=False, return_meanL2=False):
    no_chls = len(meta.plot_channels)

    means = []

    fig, ax = plt.subplots(no_chls*2, figsize=(12, 7.5), sharex=True)   # create figure.
    fig.set_tight_layout(True)

    stndata = d.select(station=stn)

    plt.suptitle(f"Station {stn}: Lat. = {stndata[0].stats.coordinates.latitude}, Lon. = {stndata[0].stats.coordinates.longitude}; {meta.fmin_mHz}-{meta.fmax_mHz} mHz")

    # Initialise L2 misfit time-series


    chlctr = 0
    for chl in meta.plot_channels:


        chldata  = stndata.select(channel=chl)

        found1 = False
        found2 = False
        for i in chldata.traces:
            if i.stats.label == label1:
                l1 = i
                found1 = True
            elif i.stats.label == label2:
                l2 = i
                found2 = True

        # Check both are found:
        if np.logical_and(found1, found2):
           # Found both
            pass
        else:
            raise ValueError(f"Error: Found: Label {label1} {found1}  - Label {label2} {found2} ")



        # Plot original
        l1x, l1y = obspy_gen_mpl_withtime(l1, starttime=-meta.time_offset[l1.stats.label.lower()])
        l2x, l2y = obspy_gen_mpl_withtime(l2, starttime=-meta.time_offset[l2.stats.label.lower()])
        if normalise_seis == True:
            l1y = normalise(l1y)
            l2y = normalise(l2y)

        ax[chlctr].plot(l1x, l1y, 'k')
        ax[chlctr].plot(l2x, l2y, 'green')


        L2norm, mean, l1x, l1y, l2x, l2y = _calc_L2(l1,
                                                    l2,
                                                    resample=resample,
                                                    normalise_seis=normalise_seis,
                                                    calc_L2_function=True,
                                                    meta=meta,
                                                    pick=meta.picks[stn])
        means.append(mean)

        ax[chlctr].plot(l1x, l1y)
        ax[chlctr].plot(l2x, l2y)

        if normalise_seis:
            ax[chlctr].set_ylim([-1,1])

        #ax[chlctr].plot(l2x[mask], l2ymasked, ':')
        ax[chlctr].set_title(chl)


        ax[chlctr+no_chls].set_title(f"mean: {np.around(mean,4)*100} %")

        ax[chlctr+no_chls].plot(l2x, L2norm)


        ax[chlctr+no_chls].set_ylim([0, 0.4])
        ax[chlctr+no_chls].set_xlim([0, 17000])


        chlctr += 1

    print(stn, means)

    if return_meanL2:
        return fig, means
    else:
        return fig# , ax, normsum


def _calc_L2(l1,l2, resample, calc_L2_function, meta, normalise_seis=False, pick=0):

    if resample != None:
        # print(f'Resampling for L2 misfit: samprate {resample} Hz')
        l1 = l1.resample(sampling_rate=resample)
        l2 = l2.resample(sampling_rate=resample)

    l1x, l1y = obspy_gen_mpl_withtime(l1, starttime=-meta.time_offset[l1.stats.label.lower()])
    l2x, l2y = obspy_gen_mpl_withtime(l2, starttime=-meta.time_offset[l2.stats.label.lower()])

    # Start at the begininng of the latest time:
    mint = np.max([l1x[0], l2x[0]])
    l1xmask = l1x >= mint
    l2xmask = l2x >= mint
    l1x = l1x[l1xmask]
    l1y = l1y[l1xmask]
    l2x = l2x[l2xmask]
    l2y = l2y[l2xmask]

    # NOW slice by pick time - Consider only the parts after the SPECFEMX Z pick
    l1mask = l1x>pick
    l2mask = l2x>pick

    l1x = l1x[l1mask]
    l1y = l1y[l1mask]
    l2x = l2x[l2mask]
    l2y = l2y[l2mask]

    # Get shortest length (slicing at the top end)
    len1 = len(l1x)
    len2 = len(l2x)
    lmin = np.min([len1, len2])

    # cutoff
    l1x = l1x[:lmin]
    l2x = l2x[:lmin]
    l1y = l1y[:lmin]
    l2y = l2y[:lmin]

    if normalise_seis == True:
        l1y = normalise(l1y)
        l2y = normalise(l2y)

    # Get maxval for this channel:
    maxval = np.amax([np.amax(np.abs(l1y)), np.amax(np.abs(l2y))])

    square_misfit = np.square(l2y - l1y) / (maxval ** 2)
    L2norm = np.sqrt(square_misfit)

    mean = np.mean(L2norm)

    if calc_L2_function:
        return L2norm, mean, l1x, l1y, l2x, l2y
    else:
        return L2norm, mean