Plotting Examples
=================

The examples below show how wrf-python can be used to make plots with 
matplotlib (with basemap and cartopy) and PyNGL.  None of these examples 
make use of xarray's builtin plotting functions, since additional work is most
likely needed to extend xarray in order to work correctly.  This is planned 
for a future release.

Matplotlib With Cartopy
-------------------------

Cartopy is becoming the main tool for base mapping with matplotlib, but you should 
be aware of a few shortcomings when working with WRF data.

- The builtin tranformations of coordinates when calling the contouring functions
  do not work correctly with the rotated pole projection.  The 
  transform_points method needs to be called manually on the latitude and 
  longitude arrays.
  
- The rotated pole projection requires the x and y limits to be set manually
  using set_xlim and set_ylim.

- You can't place latitude and longitude labels on the axes when using 
  any projection other than Mercator or LatLon.


Plotting a Two-dimensional Field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: _static/images/cartopy_slp.png    
   :scale: 100%
   :align: center
   
.. code-block:: python

    from __future__ import (absolute_import, division, print_function, unicode_literals)
    
    from netCDF4 import Dataset   
    import matplotlib.pyplot as plt
    from matplotlib.cm import get_cmap
    import cartopy.crs as crs
    from cartopy.feature import NaturalEarthFeature
    
    from wrf import npvalues, getvar, smooth2d
    
    ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")
    
    # Get the sea level pressure
    slp = getvar(ncfile, "slp")
    
    # Smooth the sea level pressure since it tends to be noisy near the mountains
    smooth_slp = smooth2d(slp, 3)
    
    # Get the numpy array from the XLAT and XLONG coordinates
    lats = npvalues(slp.coords["XLAT"])
    lons = npvalues(slp.coords["XLONG"])
    
    # Get the wrf.WrfProj object
    wrf_proj = slp.attrs["projection"]
    
    # The WrfProj.cartopy() method returns a cartopy.crs projection object
    cart_proj = wrf_proj.cartopy()
    
    # Create a figure that's 10x10
    fig = plt.figure(figsize=(10,10))
    
    # Get the GeoAxes set to the projection used by WRF
    ax = plt.axes(projection=cart_proj)
    
    # Download and add the states and coastlines
    states = NaturalEarthFeature(category='cultural', scale='50m', facecolor='none',
                                 name='admin_1_states_provinces_shp')
    ax.add_feature(states, linewidth=.5)
    ax.coastlines('50m', linewidth=0.8)
    
    # Make the contour outlines and filled contours for the smoothed sea level pressure.
    # The transform keyword indicates that the lats and lons arrays are lat/lon coordinates and tells 
    # cartopy to transform the points in to the WRF projection set for the GeoAxes.
    plt.contour(lons, lats, npvalues(smooth_slp), 10, colors="black", transform=crs.PlateCarree())
    plt.contourf(lons, lats, npvalues(smooth_slp), 10, transform=crs.PlateCarree())
    
    # Add a color bar
    plt.colorbar(ax=ax, shrink=.47)
    
    # Set the map limits
    ax.set_xlim(wrf_proj.cartopy_xlim())
    ax.set_ylim(wrf_proj.cartopy_ylim())
    
    # Add the gridlines
    ax.gridlines()
