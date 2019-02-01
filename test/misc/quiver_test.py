from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

from wrf import (getvar, interplevel, to_np, latlon_coords, get_cartopy,
                 cartopy_xlim, cartopy_ylim)

# Open the NetCDF file
ncfile = Dataset("wrfout_d01_2016-10-07_00_00_00")

# Extract the pressure and wind variables

p = getvar(ncfile, "pressure")
# Note: This is m/s.
ua = getvar(ncfile, "ua")
va = getvar(ncfile, "va")

# Interpolate u, and v winds to 950 hPa
u_950 = interplevel(ua, p, 950)
v_950 = interplevel(va, p, 950)

# Get the lat/lon coordinates
lats, lons = latlon_coords(u_950)

# Get the map projection information
cart_proj = get_cartopy(u_950)

# Create the figure
fig = plt.figure(figsize=(12, 9))
ax = plt.axes(projection=cart_proj)

# Download and add the states and coastlines
states = NaturalEarthFeature(category='cultural', scale='50m',
                             facecolor='none',
                             name='admin_1_states_provinces_shp')
ax.add_feature(states, linewidth=0.5)
ax.coastlines('50m', linewidth=0.8)

# Add the 950 hPa wind barbs, only plotting every 'thin'ed barb. Adjust thin
# as needed. Also, no scaling is done for the arrows, so you might need to
# mess with the scale argument.
thin = 50
plt.quiver(to_np(lons[::thin, ::thin]), to_np(lats[::thin, ::thin]),
           to_np(u_950[::thin, ::thin]), to_np(v_950[::thin, ::thin]),
           transform=crs.PlateCarree())

# Set the map bounds
ax.set_xlim(cartopy_xlim(u_950))
ax.set_ylim(cartopy_ylim(v_950))

ax.gridlines()

plt.title("Arrows")

plt.show()
