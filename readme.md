# MultiInstrumentLoaders

## **get_limbsounders()**

Unified loader for AIRS, ACE, GNSS, HIRDLS, MIPAS,MLS, SABER and SOFIE, loading data from the below-listed data formats. See file headers for options.

Several child functions used in this routine can currently be found in https://github.com/corwin365/MatlabFunctions; over time these will become better-integrated into this repository.

    
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
    
**AIRS:** <br />
    Lars Hoffman's 3D  retrieval netCDF files: https://datapub.fz-juelich.de/slcs/airs/gravity_waves/data/retrieval/  
    _Important note: this function returns the AIRS data as a 2D heights vs profiles matrix, which is not standard for AIRS as it takes no account of the granular structure, and is thus only useful in some cases. It also thins out the data by default to reduce volume- how much it does this can be specified in the function._

**ACE, GNSS, Misc:** <br />
    custom Matlab save format produced from the original ACE and GNSS files - contact us for details. With additional  options set to locate the files, any other files in this data format can also be loaded.

<br /><br />

## **gwanalyse_limb()**

Unified function to compute gravity waves using data outputted in the format produced by _get_limbsounders()_.  See file headers for options.  

Several functions used in this routine can currently be found in https://github.com/corwin365/MatlabFunctions; over time these will become better-integrated into this repository.



<br /><br /><br /><br />
# Structs

Functions for handling Matlab structs more efficiently, plus some functions used internally in these.

## Struct-handling functions

**cat_struct()** <br />
    For a pair of structs containing fields which share a length in one (or more) dimensions, concatenate all fields along the common specified dimension. Fields which are exceptions can be set to ignore.

**reduce_struct()** <br />
    For a structure containing a series of fields which share a length in one (or more) dimensions, index all fields along the common specific dimension. Fields which are exceptions can be set to ignore.
    
**spawn_uniform_struct()** <br />
    Produce a struct containing a series of identically-sized fields.   

## Other

Primarily intended for internal use in the above, but can be useful independently.

**expose_dim()** <br />
    Take an n-dimensional Matlab array and reshape such that the array becomes 2D with the chosen dimension as the first and all other dimensions merged in the second.  Useful for e.g. feeding multiple independent series in higher dimensions into functions that can vectorise do multiple 1D operations like interp1().

**index_dim()** <br />
    Select elements of an array along a selected dimension.  


**list_non_modal_size()** <br />
    List the fields in a structure which differ from the most common size.
    



