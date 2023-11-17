# **get_limbsounders()**

Unified loader for AIRS, ACE, GNSS, HIRDLS, MIPAS,MLS, SABER and SOFIE, loading data from the following data formats:
 
**AIRS:** <br />
    Lars Hoffman's 3D  retrieval netCDF files: https://datapub.fz-juelich.de/slcs/airs/gravity_waves/data/retrieval/  
    
**HIRDLS:** <br />
    standard HDF5 product from NASA: https://disc.gsfc.nasa.gov/  
    
**MIPAS:** <br />
    the Oxford product, available from CEDA as netCDF: https://catalogue.ceda.ac.uk/uuid/2d14cb23b3318b12c0aea374998e4a4d  

**MLS:** <br />
    standard HDF5 product from NASA: https://disc.gsfc.nasa.gov/  

**SABER:** <br />
    monthly T-O3-H2O files from GATS: https://data.gats-inc.com/saber/custom/Temp_O3_H2O/v2.0/  

**SOFIE:** <br />
    GATS level 2 netCDF format: https://data.gats-inc.com/sofie/  

**ACE, GNSS:** <br />
    custom Matlab save format produced from the original data files - contact us for details  

See file headers for options.

Several functions used in this routine can currently be found in https://github.com/corwin365/MatlabFunctions; over time these will become better-integrated into this repository.

# **gwanalyse_limb()**

Unified function to compute gravity waves using data outputted in the format produced by _get_limbsounders()_. 

See file headers for options.  

Several functions used in this routine can currently be found in https://github.com/corwin365/MatlabFunctions; over time these will become better-integrated into this repository.



