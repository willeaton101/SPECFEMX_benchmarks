# _______________________________________________________________________________________________________________________
# Author:       W Eaton, Princeton Uni. 2022
# Contact:      weaton@princeton.edu
# Last edit:    20th Jan 2022
# Notes:
#   Script for rotating a stream. This is the same as in ObsPy but I write it out for clarity
# ______________________________________________________________________________________________________________________
import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.signal import rotate

# Following the method used by JT in NMSYNG (DT98 eqn 10.3)
def rotate_stream_data(stream, method, meta):
    for stn in meta.stn_list:
        # Get station coordinates:
        stn_stream = stream.select(station=stn)
        lat_stn = float(stn_stream[0].stats.coordinates.latitude)
        lon_stn = float(stn_stream[0].stats.coordinates.longitude)

        pi = np.pi
        lat_src = meta.src_lat
        lon_src = meta.src_lon


        # Convert geographic decimal --> geocentric radian values (if geoco==1 then geographic == geocentric):
        # Note here that theta is the CO-latitude
        theta_src = (pi/2) - np.arctan(meta.geoco * np.tan(lat_src * pi/180 ))  # Source colatitude
        theta_stn = (pi/2) - np.arctan(meta.geoco * np.tan(lat_stn * pi/180 ))  # Station colatitude

        phi_src   = lon_src*pi/180                                               # Source longitude
        phi_stn   = lon_stn*pi/180                                               # Station longitude

        # Calculate epicentral distance $\Theta$ (scalar):
        dist = np.arccos(np.cos(theta_stn)*np.cos(theta_src) + np.sin(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

        rot1 = (1/np.sin(dist)) * \
               (np.sin(theta_stn)*np.cos(theta_src)  -  np.cos(theta_stn)*np.sin(theta_src)*np.cos(phi_stn - phi_src))

        rot2 = (1/np.sin(dist)) * (np.sin(theta_src)*np.sin(phi_stn - phi_src))

        # Conversion from RTP --> ZNE (where R=Z, 2D rotation) appears to use the following matrix:
        #   [N, E]' = [-rot1, -rot2; -rot2, rot1][T, P]' where T and P are theta, Phi
        #   Below we shall name the rotation matrix Q:
        # Hence to get the T and P matrix we should be multiplying [N,E] by the inverse of Q:
        Q    = np.array([[-rot1, -rot2], [rot2, -rot1]])
        Qinv = np.linalg.inv(Q)


        if method == "NE->TP":
            N = stn_stream.select(channel="N")[0].data
            E = stn_stream.select(channel="E")[0].data
            data_NE = np.array([N,E])
            data_TP = np.matmul(Qinv, data_NE)

            # Now writing back to stream:
            old_chls = ["N", "E"]
            new_chls = ["T", "P"]
            for i in range(2):
                if new_chls[i] == "P":
                    data_TP[i, :] = data_TP[i,:]*(-1)
                stream.select(station=stn, channel=old_chls[i])[0].data = data_TP[i,:]
                stream.select(station=stn, channel=old_chls[i])[0].stats.channel = new_chls[i]
        else:
            raise ValueError("Currently method must be NE->TP")




def rotate_stream_with_obspy(stream, method, meta, radius=6378137.0, geoco=None):

    # Sort out the geoco factor
    if geoco==None:
        print('Rotating with Geoco from meta:', meta.geoco)
        geoco = meta.geoco
    if geoco=='WGS84':
        geoco = 1 - 0.0033528106647474805
    f = 1 - geoco
    print('Rotating with f = 1 - geoco: f =', f)


    for stn in meta.stn_list:
        # Get station coordinates:
        stn_stream = stream.select(station=stn)
        lat_stn = float(stn_stream[0].stats.coordinates.latitude)
        lon_stn = float(stn_stream[0].stats.coordinates.longitude)

        # Source coordinates
        lat_src = meta.src_lat
        lon_src = meta.src_lon

        # Data for traces
        N = stn_stream.select(channel="N")[0].data
        E = stn_stream.select(channel="E")[0].data

        # get the backazimuth used for rotation #
        # returns (Great circle distance in m, azimuth A->B in degrees, azimuth B->A in degrees)
        distance, azimuth, backazimuth = gps2dist_azimuth(lat_src, lon_src, lat_stn, lon_stn, a=radius, f=f)

        if method=="NE->TP":
            radial, transverse = rotate.rotate_ne_rt(n=N, e=E, ba=backazimuth)

            stream.select(station=stn, channel='N')[0].data = radial
            stream.select(station=stn, channel='N')[0].stats.channel = 'T'

            stream.select(station=stn, channel='E')[0].data = transverse*(-1)
            stream.select(station=stn, channel='E')[0].stats.channel = 'P'

        else:
            raise ValueError('Only method "NE->TP" is available right now')


