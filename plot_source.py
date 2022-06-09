import numpy as np
from obspy.imaging.beachball import beachball

"""
BOLIVIA
PDE 1994  6  9  0 33 16.40 -13.8300  -67.5600 637.0 6.9 6.8 NORTHERN BOLIVIA
event name:     060994A
time shift:     29.0000
half duration:  20.0000
latitude:      -13.8200
longitude:     -67.2500
depth:         647.1000
Mrr:      -7.590000e+27
Mtt:       7.750000e+27
Mpp:      -1.600000e+26
Mrt:      -2.503000e+28
Mrp:       4.200000e+26
Mtp:      -2.480000e+27"""

mt = [-7.590000e+27, 7.750000e+27, -1.600000e+26, -2.503000e+28, 4.200000e+26, -2.480000e+27]
#mt = [0, 0, 0, 1, 0, 1]
beachball(mt, size=800, linewidth=1.5, facecolor='k')

import plotly.express as px
import plotly.graph_objects as go
import geopandas as gpd
import pandas as pd


lat = list(np.linspace(-80, 80, 17))
lat.append(-13.82)
lon = list(np.zeros(17))
lon.append(-67.25)
ms = list(np.zeros(18)+1)


stations = pd.DataFrame(
    {'Stns': ['X11', 'X21', 'X31', 'X41', 'X51', 'X61', 'X71', 'X81', 'X91', 'X101',
              'X111', 'X121', 'X131', 'X141', 'X151', 'X161', 'X171', 'src'],
     'Latitude': lat,
     'Longitude': lon,
     'msize': ms})



fig = px.scatter_geo(stations, lat=stations.Latitude, lon=stations.Longitude,
                symbol_sequence=['triangle-down', 'triangle-down', 'triangle-down', 'triangle-down', 'triangle-down', 'triangle-down', 'triangle-down', 'triangle-down', 'triangle-down', 'triangle-down',
                                 'triangle-down', 'triangle-down', 'triangle-down', 'triangle-down', 'triangle-down','triangle-down','circle'],
                size=stations.msize,
                size_max=3)

fig.update_geos(projection_type="orthographic")

fig.show()
#fig.write_image("fig1.pdf")
