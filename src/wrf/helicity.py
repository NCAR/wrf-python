from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)

import numpy as np

from .constants import Constants
from .extension import _srhel, _udhel
from .destag import destagger
from .util import extract_vars, extract_global_attrs, either
from .metadecorators import copy_and_set_metadata

@copy_and_set_metadata(copy_varname="HGT", name="srh", 
                       description="storm relative helicity",
                       units="m2 s-2")
def get_srh(wrfin, timeidx=0, method="cat", squeeze=True, 
            cache=None, meta=True, _key=None, top=3000.0):
    """Return the storm relative helicity.
    
    The *top* argument specifies the top of the integration in [m].
    
    This functions extracts the necessary variables from the NetCDF file 
    object in order to perform the calculation.
    
    Args:
    
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF 
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile` 
            or an iterable sequence of the aforementioned types.
        
        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The 
            desired time index. This value can be a positive integer, 
            negative integer, or 
            :data:`wrf.ALL_TIMES` (an alias for None) to return 
            all times in the file or sequence. The default is 0.
        
        method (:obj:`str`, optional): The aggregation method to use for 
            sequences.  Must be either 'cat' or 'join'.  
            'cat' combines the data along the Time dimension.  
            'join' creates a new dimension for the file index.  
            The default is 'cat'.
        
        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions 
            with a size of 1 from being automatically removed from the shape 
            of the output. Default is True.
        
        cache (:obj:`dict`, optional): A dictionary of (varname, ndarray) 
            that can be used to supply pre-extracted NetCDF variables to the 
            computational routines.  It is primarily used for internal 
            purposes, but can also be used to improve performance by 
            eliminating the need to repeatedly extract the same variables 
            used in multiple diagnostics calculations, particularly when using 
            large sequences of files. 
            Default is None.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        _key (:obj:`int`, optional): A caching key. This is used for internal 
            purposes only.  Default is None.
            
        top (:obj:`float`, optional): The top of the integration in [m].  
            Default is 3000.0.
   
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        storm relative helicity.
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    # Top can either be 3000 or 1000 (for 0-1 srh or 0-3 srh)
    
    ncvars = extract_vars(wrfin, timeidx, ("HGT", "PH", "PHB"),
                          method, squeeze, cache, meta=False,
                          _key=_key)
    
    ter = ncvars["HGT"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    
    # As coded in NCL, but not sure this is possible
    varname = either("U", "UU")(wrfin)
    u_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    u = destagger(u_vars[varname], -1) 
    
    varname = either("V", "VV")(wrfin)
    v_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    v = destagger(v_vars[varname], -2)

    geopt = ph + phb
    geopt_unstag = destagger(geopt, -3)
    
    z = geopt_unstag / Constants.G
    
    # Re-ordering from high to low
    u1 = np.ascontiguousarray(u[...,::-1,:,:])
    v1 = np.ascontiguousarray(v[...,::-1,:,:])
    z1 = np.ascontiguousarray(z[...,::-1,:,:])
    
    srh = _srhel(u1, v1, z1, ter, top)
    
    return srh

@copy_and_set_metadata(copy_varname="MAPFAC_M", name="updraft_helicity", 
                       description="updraft helicity",
                       units="m2 s-2")
def get_uh(wrfin, timeidx=0, method="cat", squeeze=True, 
           cache=None, meta=True, _key=None,
           bottom=2000.0, top=5000.0):
    
    """Return the updraft helicity.
    
    The *bottom* and *top* arguments specify the bottom and top limits 
    for the integration in [m].
    
    This functions extracts the necessary variables from the NetCDF file 
    object in order to perform the calculation.
    
    Args:
    
        wrfin (:class:`netCDF4.Dataset`, :class:`Nio.NioFile`, or an \
            iterable): WRF-ARW NetCDF 
            data as a :class:`netCDF4.Dataset`, :class:`Nio.NioFile` 
            or an iterable sequence of the aforementioned types.
        
        timeidx (:obj:`int` or :data:`wrf.ALL_TIMES`, optional): The 
            desired time index. This value can be a positive integer, 
            negative integer, or 
            :data:`wrf.ALL_TIMES` (an alias for None) to return 
            all times in the file or sequence. The default is 0.
        
        method (:obj:`str`, optional): The aggregation method to use for 
            sequences.  Must be either 'cat' or 'join'.  
            'cat' combines the data along the Time dimension.  
            'join' creates a new dimension for the file index.  
            The default is 'cat'.
        
        squeeze (:obj:`bool`, optional): Set to False to prevent dimensions 
            with a size of 1 from being automatically removed from the shape 
            of the output. Default is True.
        
        cache (:obj:`dict`, optional): A dictionary of (varname, ndarray) 
            that can be used to supply pre-extracted NetCDF variables to the 
            computational routines.  It is primarily used for internal 
            purposes, but can also be used to improve performance by 
            eliminating the need to repeatedly extract the same variables 
            used in multiple diagnostics calculations, particularly when using 
            large sequences of files. 
            Default is None.
        
        meta (:obj:`bool`, optional): Set to False to disable metadata and 
            return :class:`numpy.ndarray` instead of 
            :class:`xarray.DataArray`.  Default is True.
            
        _key (:obj:`int`, optional): A caching key. This is used for internal 
            purposes only.  Default is None.
        
        bottom (:obj:`float`, optional): The bottom limit for the integration 
            in [m]. Default is 2000.0.
           
        top (:obj:`float`, optional): The top limit for the integration in [m].  
            Default is 5000.0.
   
    Returns:
        :class:`xarray.DataArray` or :class:`numpy.ndarray`: The 
        updraft helicity.
        If xarray is enabled and the *meta* parameter is True, then the result 
        will be a :class:`xarray.DataArray` object.  Otherwise, the result will 
        be a :class:`numpy.ndarray` object with no metadata.
    
    """
    
    ncvars = extract_vars(wrfin, timeidx, ("W", "PH", "PHB", "MAPFAC_M"),
                          method, squeeze, cache, meta=False, _key=None)
    
    wstag = ncvars["W"]
    ph = ncvars["PH"]
    phb = ncvars["PHB"]
    mapfct = ncvars["MAPFAC_M"]
    
    attrs  = extract_global_attrs(wrfin, attrs=("DX", "DY"))
    dx = attrs["DX"]
    dy = attrs["DY"]
    
    # As coded in NCL, but not sure this is possible
    varname = either("U", "UU")(wrfin)
    u_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    u = destagger(u_vars[varname], -1) 
    
    varname = either("V", "VV")(wrfin)
    v_vars = extract_vars(wrfin, timeidx, varname, method, squeeze, cache,
                          meta=False, _key=_key)
    v = destagger(v_vars[varname], -2) 
    
    zp = ph + phb
    
    uh = _udhel(zp, mapfct, u, v, wstag, dx, dy, bottom, top)
    
    return uh

    
    
    
    