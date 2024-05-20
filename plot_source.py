import numpy as np
from obspy.imaging.beachball import beachball
import plotly.express as px
import pandas as pd

# Plot MT with obspy for Bolivia:
mt = [-7.590000e+27, 7.750000e+27, -1.600000e+26, -2.503000e+28, 4.200000e+26, -2.480000e+27]

# Plot MT with obspy for Tohoku:
#mt = [1.730000e+29, -2.810000e+28, -1.450000e+29, 2.120000e+29, 4.550000e+29, -6.570000e+28]
beachball(mt, size=800, linewidth=1.5, facecolor='k')





# Plot stations on map:
lat = list(np.linspace(-80, 80, 17))
lon = list(np.zeros(17))

# Bolivia
#lat.append(-13.82)
#lon.append(-67.25)

# Tohoku
lat.append(37.5200)
lon.append(143.0500)


ms = list(np.zeros(18)+50)

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
