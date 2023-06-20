# This script is basicaclly reloading/processing/comparing SPECFEMX results for different NEX at increasing frequencies
# until there are discrepencies
import matplotlib.pyplot as plt
from get_data import get_data, process_stream
from plot_funcs import plot_record, plot_station, plot_spectra
import obspy
import numpy as np
from rotate_stream import rotate_stream_data
from metadata import PlotMetadata


legend_colours = {
  "NEX192": "tomato",
  "NEX128": "rebeccapurple",
  "NEX96" : "royalblue"
}
pdf_name = "CONVERGENCE_TEST.pdf"

freqs = [[0.001, 0.01],  [0.001, 0.016], [0.001, 0.02], [0.001, 0.03], [0.001, 0.04]]

chl = 'G'


fig_list = []
for f in freqs:

    meta = PlotMetadata(chls= "ZNEG",                      # Channel names to load
                        out_chls = ["G"],       # Channels for plotting
                        rotate = False,                   # I dont think this actually does anything!
                        gravity_subtract=False,
                        geoco = 1.0,                        # Ellipticity
                        fmin = f[0],
                        fmax = f[1],
                        ev_coords = (-13.82, -67.25),     # Geographic source coordinates
                        tmin=0,                           # For slicing time array
                        tmax = 10000,                     # For slicing time array
                        colours = legend_colours,         # As above
                        t_offset_spfmx = 8.999999999e-02,               # Time offset for SPFMX data
                        t_offset_nmsyn = 1.439962400e-02,   # Time offset for NMSYN data
                        t_offset_yspec = 1.439962400e-02,   # Time offset for YSPEC data
                        t_offset_axisem=9.6300000000e-02,   # Time offset for AXISEM data
                        no_stns=1)                          # You will need to edit gen_stn_list() in plot_meta_data_class.py
                                                            # Depending on what the stations are called in data file names


    st_spfmx_96  = get_data(type="spfmx", filepath="./test1/self_gravitation/specfemx/NEX96/conv", meta=meta,  label="NEX96")
    rotate_stream_data(stream=st_spfmx_96, method="NE->TP", meta=meta)
    st_spfmx128  = get_data(type="spfmx", filepath="./test1/self_gravitation/specfemx/NEX128/conv", meta=meta,  label="NEX128")
    rotate_stream_data(stream=st_spfmx128, method="NE->TP", meta=meta)
    st_spfmx192  = get_data(type="spfmx", filepath="./test1/self_gravitation/specfemx/NEX192_May26/conv", meta=meta,  label="NEX192")
    rotate_stream_data(stream=st_spfmx192, method="NE->TP", meta=meta)


    data_list = [st_spfmx_96, st_spfmx128, st_spfmx192]
    st = obspy.Stream()
    for stream in data_list:
        st +=  process_stream(stream, meta=meta)


    # Station figures
    for station in meta.stn_list:
        fig_list.append(plot_station(st=st, station=station, chls=meta.plot_channels, meta=meta, style='light'))



# Write output to a pdf
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_name)
for fig in fig_list:
    pdf.savefig(fig)
pdf.close()