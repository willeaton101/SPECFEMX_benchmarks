# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   This is the script I use day-to-day for comparing results.
# ______________________________________________________________________________________________________________________
import obspy
from get_data import get_data, process_stream
from plot_funcs import plot_record, plot_station, plot_spectra, plot_synthetic_L2misfit
from rotate_stream import  rotate_stream_with_obspy, rotate_stream_data
from metadata import PlotMetadata


fig_list = []
# ______________________________________________________________________________________________________________________
# SPECIFY METADATA AND PARAMETERS FOR PROCESSING/PLOTTING:
# Feel free to edit these - name needs to correspond to the label you have given the data when it loads
legend_colours = {
  "SPECFEMX": "tomato", #"SPECFEMX": "sandybrown"
  "Cowling": "tomato", #"SPECFEMX": "sandybrown"
  "NMSYN": "black", #"NMSYN": "royalblue",
  "YSPEC": "purple", #"YSPEC": "purple",
  "Self Gravitation": "rebeccapurple",
  "YSPEC_ell": "salmon",
  "pure": "y"  ,
  "test2": "r",
  "NEX96": "orange",
  "S40RTS" : "royalblue",
  "NMSYN_alt": "plum",
  "WGS84" : "cyan"
}

fmin = 0.002
fmax = 0.02
meta = PlotMetadata(chls= "ZTPG",                      # Channel names to load
                    out_chls = ["Z","T","P", "G"],   # Channels for plotting,
                    types=["yspec", "nmsyn", "spfmx"],
                    rotate = False,                   # I dont think this actually does anything!
                    gravity_subtract=False,
                    geoco = 1.0,                        # Ellipticity
                    fmin = fmin,
                    fmax = fmax,
                    ev_coords = (-13.82, -67.25),     # Geographic source coordinates
                    tmin=0,                           # For slicing time array
                    tmax = 17000,                     # For slicing time array
                    colours = legend_colours,         # As above
                    t_offset_spfmx  = 8.999999999e-02,               # Time offset for SPFMX data
                    t_offset_nmsyn  = 1.439962400e-02,   # Time offset for NMSYN data
                    t_offset_yspec  = 1.439962400e-02,   # Time offset for YSPEC data
                    t_offset_axisem = 9.6300000000e-02,  # Time offset for AXISEM data
                    no_stns=2)                          # You will need to edit gen_stn_list() in plot_meta_data_class.py
                                                         # Depending on what the stations are called in data file names

# ______________________________________________________________________________________________________________________
# LOAD DATA
# BEST CURRENT RESULTS:
#st_yspec  = get_data(type="yspec", filepath="./test2/self_gravitation/yspec/conv", meta=meta, yspec="conv_selfgrav_att.yspec", label="YSPEC")
#st_nmsyn  = get_data(type="nmsyn", filepath="./test2/self_gravitation/nmsyn/conv", meta=meta,  label="NMSYN")

# Test 1 SG:
"""st_nmsyn  = get_data(type="nmsyn", filepath="./test1/self_gravitation/nmsyn/conv", meta=meta,  label="NMSYN")
meta.set_channels("ZNEG")
st_spfmx_sg      = get_data(type="spfmx", filepath="./test1/self_gravitation/specfemx/NEX192_May26/conv", yspec=1, meta=meta,  label="SPECFEMX")
rotate_stream_data(stream=st_spfmx_sg, method="NE->TP", meta=meta)
# Process the data:
data_list = [ st_nmsyn, st_spfmx_sg]               # List of streams you want to put into some master, processed stream
pdf_name = f"check_grav_ampSG_{fmin}_{fmax}.pdf"
"""

# Test 2 SG:
st_nmsyn  = get_data(type="nmsyn", filepath="./test2/self_gravitation/nmsyn/conv", meta=meta,  label="NMSYN")
st_yspec  = get_data(type="yspec", filepath="./test2/self_gravitation/yspec/alt/conv",
                     meta=meta, yspec="conv_TEST2ALT.yspec", label="YSPEC")

meta.set_channels("ZNEG")
#st_spfmx_sg      = get_data(type="spfmx", filepath="./test2/self_gravitation/specfemx/conv", yspec=1, meta=meta,  label="SPECFEMX")
#rotate_stream_data(stream=st_spfmx_sg, method="NE->TP", meta=meta)

st_nmsynalt  = get_data(type="nmsyn", filepath="./test2/self_gravitation/nmsyn/alt/conv", meta=meta,  label="NMSYN_alt")
rotate_stream_data(stream=st_nmsynalt, method="NE->TP", meta=meta)

# Process the data:
data_list = [ st_nmsyn, st_yspec, st_nmsynalt]               # List of streams you want to put into some master, processed stream
pdf_name = f"test2_check_{fmin}_{fmax}.pdf"




st = obspy.Stream()                            # This will be the master trace with all data
for stream in data_list:
    st +=  process_stream(stream, meta=meta, reverse_P=False)   # Check this as it might need editing for your simulation params.


# Plot the data; initialise list to hold figures for outputting to user:


# Record sections
for c in meta.plot_channels:
    fig_list.append(plot_record(st.select(channel=c), meta, amplitude=2, orientation="horizontal", style='light'))

# Station figures
for station in meta.stn_list:

    """fig_list.append(plot_synthetic_L2misfit(d=st, stn=station,
                                            label1='NMSYN',
                                            label2='SPECFEMX',
                                            meta=meta,
                                            resample=4.0,
                                            mismax=15,
                                            normalise_seis=True))"""

    fig_list.append(plot_station(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light', norm_amp=True))
    #plot_station(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light')
    #plot_station(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light')
    #plot_spectra(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light')


# Write output to a pdf
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name)
for fig in fig_list:
    pdf.savefig(fig)
pdf.close()
