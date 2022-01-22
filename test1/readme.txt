## Test 1

We compare synthetic seismograms from two different programs with SPECFEMX. For a simulation with
self-gravitation I use YSPEC and NMSYN. For the simulation using a Cowling approximation, only YSPEC results
are used. Other simulation parameters are as follows:

Earth model:
- Completely isotropic PREM with ocean replaced by crust. One-layered crust.
- No attenuation
- Cowling or Self-gravitation depending on simulation
- No elliticity (geoco = 1.0)

Modes:
- All modes within bandwidth
- 0.2 - 25 mHz (40 - 5000 seconds)
- Results then filtered between 50-500 seconds

Source:
See CMTSOLUTION in nmsyng inputs

Stations:
See STATIONS_18 in nmsyng inputs