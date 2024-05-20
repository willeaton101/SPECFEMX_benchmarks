# Function to rotate the specfemx data in the hope that I wont have to do this in the future
# USE THIS WITH CAUTION - THERE ARE A LOT OF HARD-CODED PARTS ESPECIALLY IN FILE PATHS
import numpy as np
from get_data import get_data
from metadata import PlotMetadata
from rotate_stream import rotate_stream_data
# Load raw data:

meta = PlotMetadata(chls= "ZNE",                     # Channel names to load
                    out_chls = None,      # Channels for plotting
                    rotate = False,                   # I dont think this actually does anything!
                    geoco = 1,                        # Ellipticity
                    fmin = 0.00001,
                    fmax = 100,
                    ev_coords = (-13.82, -67.25),     # Geographic source coordinates
                    tmin=0,                           # For slicing time array
                    tmax = 20000,                     # For slicing time array
                    colours = None,                   # As above
                    t_offset_spfmx = 0)                        # Time offset for SPFMX data



c = "ZTP"   # Output channels
# Get a time array for output later:
t = np.loadtxt("./test1/cowling/specfemx/raw/X1.Y5.MXND.sem.ascii")[:,0]


st = get_data(type="spfmx", filepath="./test1/cowling/specfemx/raw", meta=meta,  label="spfmx")

rotate_stream_data(stream=st, method="NE->TP", meta=meta)

for stn in meta.stn_list:
    for ch in c:
        tr = st.select(station=stn, channel=ch)[0]

        x = np.transpose(np.array([t, tr.data]))
        np.savetxt(fname=f"./test1/cowling/specfemx/raw_rotated/{stn}.Y5.MX{ch}D.sem.ascii", X=x, fmt="%f")