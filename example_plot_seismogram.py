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

# ______________________________________________________________________________________________________________________
# SPECIFY METADATA AND PARAMETERS FOR PROCESSING/PLOTTING:
# Feel free to edit these - name needs to correspond to the label you have given the data when it loads
legend_colours = {
  "spfmx": "royalblue",
  "nmsyn": "g",
  "selfgrav": "c",
  "high-res": "salmon",
  "ZTP": "y"  ,
  "conv": "r",
  "yspec": "orchid",
  "table": "purple",
  "": "c",
  "specfemx" : "royalblue"
}
meta = PlotMetadata(chls= "ZTPG",                     # Channel names to load
                    out_chls = ["Z", "T", "P", "G"],  # Channels for plotting
                    rotate = False,                   # I dont think this actually does anything!
                    geoco = 1,                        # Ellipticity
                    fmin = 0.002,
                    fmax = 0.02,
                    ev_coords = (-13.82, -67.25),     # Geographic source coordinates
                    tmin=0,                           # For slicing time array
                    tmax = 17500,                     # For slicing time array
                    colours = legend_colours,         # As above
                    t_offset_spfmx = 0,               # Time offset for SPFMX data compared (times series doesnt start
                                                      # at 0
                    no_stns=17)                       # You will need to edit gen_stn_list() in plot_meta_data_class.py
                                                      # Depending on what the stations are called in data file names

# ______________________________________________________________________________________________________________________
# LOAD DATA
# BEST CURRENT RESULTS:
#st_yspec  = get_data(type="yspec", filepath="./Restart_Jan14/9_mode_breakdown/yspec/tilt/conv", meta=meta, yspec="conv_selfgrav_tilt.yspec", label="yspec")
st_nmsyn  = get_data(type="nmsyn", filepath="./self_gravitation/nmsyn/conv", meta=meta,  label="nmsyn")

# Process the data:
data_list = [st_nmsyn]               # List of streams you want to put into some master, processed stream
st = obspy.Stream()                            # This will be the master trace with all data
for stream in data_list:
    #rotate_stream_data(stream, method="NE->TP", meta=meta)
    st +=  process_stream(stream, meta=meta)   # Check this as it might need editing for your simulation params.


# Plot the data; initialise list to hold figures for outputting to user:
fig_list = []

# Record sections
for c in meta.plot_channels:
    fig_list.append(plot_record(st.select(channel=c), meta, amplitude=2, orientation="horizontal"))

# Station figures
for station in meta.stn_list:
    fig_list.append(plot_station(st=st, station=station, chls=meta.plot_channels, meta=meta))
    fig_list.append(plot_spectra(st=st, station=station, chls=meta.plot_channels, meta=meta))

# Write output to a pdf
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("test_results.pdf")
for fig in fig_list:
    pdf.savefig(fig)
pdf.close()