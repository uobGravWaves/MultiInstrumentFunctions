**get_limbsounders()**

Unified loader for AIRS, ACE, GNSS, HIRDLS, MIPAS,MLS, SABER and SOFIE using the formats stored at the University of Bath:

ACE, GNSS: custom Matlab save format, contact for details
AIRS: Lars Hoffman's 3D netCDf retrieval: https://datapub.fz-juelich.de/slcs/airs/gravity_waves/data/retrieval/
HIRDLS: standard HDF5 product from NASA: https://disc.gsfc.nasa.gov/
MIPAS: the Oxford product, available from CEDA as netCDF: https://catalogue.ceda.ac.uk/uuid/2d14cb23b3318b12c0aea374998e4a4d
MLS: standard HDF5 product from NASA: https://disc.gsfc.nasa.gov/
SABER: monthly T-O3-H2O files from GATS: https://data.gats-inc.com/saber/custom/Temp_O3_H2O/v2.0/
SOFIE: GATS level 2 netCDF format: https://data.gats-inc.com/sofie/



gwanalyse_limb


Several functions used in get_limbsounders() and gwanalyse_limb() can be found in https://github.com/corwin365/MatlabFunctions . Over time these will be better-integrated.
