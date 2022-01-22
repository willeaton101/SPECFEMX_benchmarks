## Test 1
Last updated Jan 22nd 2022, Will Eaton, Princeton University

### Current test results:

##### Self gravitating:
- Both NMSYN and YSPEC results are given for the self-gravitating example. They are in very good agreement apart from small differences on the Phi (transverse) channel. 
- Results are shown in selfgrav_results.pdf 

##### Cowling approximation: 
- Based on the agreement between NMSYN and YSPEC for self-gravitating, and the fact that the cowling approximation is not currently implemented in the mode catalogues for NMSYN, I only use YSPEC for this. 
- In this case, there are clear discrepencies between the NEX256 SPECFEM-X data and the YSPEC data, particularly on the vertical and theta (radial) channels. Good agreement is seen on the Phi (transverse) channel.


### Test description and Inputs:
We compare synthetic seismograms from two different programs with SPECFEMX. For a simulation with
self-gravitation I use YSPEC and NMSYN. For the simulation using a Cowling approximation, only YSPEC results
are used. Other simulation parameters are as follows:

Earth model:
- Completely isotropic PREM with ocean replaced by crust. One-layered crust.
- No attenuation
- Cowling or Self-gravitation depending on simulation
- No ellipticity (geoco = 1.0)

Modes:
- All modes within bandwidth
- 0.2 - 25 mHz (40 - 5000 seconds)
- Results then filtered between 50-500 seconds

Source:
See CMTSOLUTION in nmsyng inputs

Stations:
See STATIONS_18 in nmsyng inputs