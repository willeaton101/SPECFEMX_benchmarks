# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   This is the script I use day-to-day for comparing results.
# ______________________________________________________________________________________________________________________

from get_data import get_data, process_stream
from plot_funcs import plot_record, plot_station, plot_spectra
import obspy
from rotate_stream import rotate_stream_data
from metadata import PlotMetadata
import matplotlib.pyplot as plt
import plotly.graph_objects as go



# ______________________________________________________________________________________________________________________
# SPECIFY METADATA AND PARAMETERS FOR PROCESSING/PLOTTING:
# Feel free to edit these - name needs to correspond to the label you have given the data when it loads
legend_colours = {
  "SPECFEMX": "tomato", #"SPECFEMX": "sandybrown"
  "NMSYN": "black", #"NMSYN": "royalblue",
  "YSPEC": "dimgrey", #"YSPEC": "purple",
  "No gravitation": "purple",
  "NEX160": "salmon",
  "pure": "y"  ,
  "test2": "r",
  "NEX192": "purple",
  "NOATT": "orange",
  "Self-gravitating" : "royalblue"
}
meta = PlotMetadata(chls= "ZTPG",                      # Channel names to load
                    out_chls = ["G"],       # Channels for plotting
                    rotate = False,                   # I dont think this actually does anything!
                    gravity_subtract=False,
                    geoco = 1.0,                        # Ellipticity
                    fmin = 0.002,
                    fmax = 0.02,
                    ev_coords = (-13.82, -67.25),     # Geographic source coordinates
                    tmin=0,                           # For slicing time array
                    tmax = 17950,                     # For slicing time array
                    colours = legend_colours,         # As above
                    t_offset_spfmx = 8.999999999e-02,               # Time offset for SPFMX data
                    t_offset_nmsyn = 1.439962400e-02,   # Time offset for NMSYN data
                    t_offset_yspec = 1.439962400e-02,   # Time offset for YSPEC data
                    t_offset_axisem=9.6300000000e-02 ,  # Time offset for AXISEM data
                    no_stns=17)                          # You will need to edit gen_stn_list() in plot_meta_data_class.py
                                                        # Depending on what the stations are called in data file names


# ______________________________________________________________________________________________________________________
# LOAD DATA
# BEST CURRENT RESULTS:
st_yspec  = get_data(type="yspec", filepath="./test1/self_gravitation/yspec/conv", meta=meta, yspec="conv_selfgrav_tilt.yspec", label="YSPEC")
st_nmsyn  = get_data(type="nmsyn", filepath="./test1/self_gravitation/nmsyn/conv", meta=meta,  label="NMSYN")


meta.set_channels("ZNEG")
st_spfmx  = get_data(type="spfmx", filepath="./test1/self_gravitation/specfemx/NEX192_May26/conv", meta=meta,  label="SPECFEMX")
rotate_stream_data(stream=st_spfmx, method="NE->TP", meta=meta)


# Process the data:
data_list = [st_yspec, st_nmsyn, st_spfmx]               # List of streams you want to put into some master, processed stream
pdf_name = "JT_even.pdf"


st = obspy.Stream()                            # This will be the master trace with all data
for stream in data_list:
    st +=  process_stream(stream, meta=meta)   # Check this as it might need editing for your simulation params.


# Plot the data; initialise list to hold figures for outputting to user:
fig_list = []






# Record sections
for c in meta.plot_channels:
    fig_list.append(plot_record(st.select(channel=c), meta, amplitude=5, orientation="horizontal", style='light'))


# Station figures
#for station in meta.stn_list:
    #fig_list.append(plot_station(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light'))
    #fig_list.append(plot_spectra(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light'))

    #plot_station(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light')
    #plot_spectra(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light')


# Write output to a pdf
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name)
for fig in fig_list:
    pdf.savefig(fig)
pdf.close()