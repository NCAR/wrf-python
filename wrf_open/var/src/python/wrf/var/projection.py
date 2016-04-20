from __future__ import (absolute_import, division, print_function, 
                        unicode_literals)
import numpy as np
import math

from .config import basemap_enabled, cartopy_enabled, pyngl_enabled
from .constants import Constants

if cartopy_enabled():
    from cartopy import crs
    
if basemap_enabled():
    from mpl_toolkits.basemap import Basemap
    
if pyngl_enabled():
    from Ngl import Resources

__all__ = ["WrfProj", "LambertConformal", "Mercator",
           "PolarStereographic", "LatLon", "RotatedLatLon",
           "getproj"]


if cartopy_enabled():
    class MercatorWithLatTS(crs.Mercator):
        def __init__(self, central_longitude=0.0,
                     latitude_true_scale=0.0,
                     min_latitude=-80.0, 
                     max_latitude=84.0,
                     globe=None):
            proj4_params = [("proj", "merc"),
                ("lon_0", central_longitude),
                ("lat_ts", latitude_true_scale),
                ("k", 1),
                ("units", "m")]
            super(crs.Mercator, self).__init__(proj4_params, globe=globe)

            # Calculate limits.
            limits = self.transform_points(crs.Geodetic(),
                               np.array([-180, 180]) + central_longitude,
                               np.array([min_latitude, max_latitude]))
            
            # When using a latitude of true scale, the min/max x-limits get set 
            # to the same value, so make sure the left one is negative
            xlimits = limits[..., 0]
            
            if xlimits[0] == xlimits[1]:
                if xlimits[0] < 0:
                    xlimits[1] = -xlimits[1]
                else:
                    xlimits[0] = -xlimits[0]
            
            self._xlimits = tuple(xlimits)
            self._ylimits = tuple(limits[..., 1])
            
            self._threshold = np.diff(self.x_limits)[0] / 720

def _ismissing(val):
    return val is None or val > 90. or val < -90.

class WrfProj(object):
    def __init__(self, bottom_left=None, top_right=None, 
                 lats=None, lons=None, **proj_params):
        if bottom_left is not None and top_right is not None:
            self.ll_lat = bottom_left[0]
            self.ll_lon = bottom_left[1]
            self.ur_lat = top_right[0]
            self.ur_lon = top_right[1]
            self.bottom_left = bottom_left
            self.top_right = top_right
        elif lats is not None and lons is not None:
            self.ll_lat = lats[0,0]
            self.ur_lat = lats[-1,-1]
            self.ll_lon = lons[0,0]
            self.ur_lon = lons[-1,-1]
            self.bottom_left = [self.ll_lat, self.ll_lon]
            self.top_right = [self.ur_lat, self.ur_lon]
        else:
            raise ValueError("invalid corner point arguments")
        
        # These indicate the center of the nest/domain, not necessarily the 
        # center of the projection
        self._cen_lat = proj_params.get("CEN_LAT", None)
        self._cen_lon = proj_params.get("CEN_LON", None)
        
        self.truelat1 = proj_params.get("TRUELAT1", None)
        self.truelat2 = (proj_params.get("TRUELAT2", None)
                         if not _ismissing(proj_params.get("TRUELAT2", None)) 
                         else None)
        self.moad_cen_lat = proj_params.get("MOAD_CEN_LAT", None)
        self.stand_lon = proj_params.get("STAND_LON", None)
        self.pole_lat = proj_params.get("POLE_LAT", None)
        self.pole_lon = proj_params.get("POLE_LON", None)
        
        # Just in case...
        if self.moad_cen_lat is None:
            self.moad_cen_lat = self._cen_lat
        
        if self.stand_lon is None:
            self.stand_lon = self._cen_lon
            
    @property
    def _basemap(self):
        return None
    
    @property
    def _cf_params(self):
        return None
    
    @property
    def _cartopy(self):
        return None
    
    @property
    def _cart_extents(self):
        return ([self.ll_lon, self.ur_lon], [self.ll_lat, self.ur_lat]) 
    
    @property
    def _pyngl(self):
        return None
    
    @property
    def _proj4(self):
        return None
    
    @property
    def _globe(self):
        return (None if not cartopy_enabled() 
                else crs.Globe(ellipse=None,
                               semimajor_axis=Constants.WRF_EARTH_RADIUS,
                               semiminor_axis=Constants.WRF_EARTH_RADIUS))
     
    def cartopy_xlim(self):
        """Return the x extents in projected coordinates (for cartopy)"""
        return self._cart_extents[0]
    
    def cartopy_ylim(self):
        """Return the y extents in projected coordinates (for cartopy)"""
        return self._cart_extents[1]
    
    def __repr__(self):
        args = ("bottom_left={}, top_right={}, "
                "stand_lon={}, moad_cen_lat={}, "
                "pole_lat={}, pole_lon={}".format((self.ll_lat, self.ll_lon),
                                                  (self.ur_lat, self.ur_lon),
                                                  self.stand_lon, 
                                                  self.moad_cen_lat,
                                                  self.pole_lat,
                                                  self.pole_lon))
        return "{}({})".format(self.__class__.__name__, args)
    
    def basemap(self):
        """Return a mpl_toolkits.basemap.Basemap instance for the
        projection"""
        if not basemap_enabled():
            raise RuntimeError("'mpl_toolkits.basemap' is not "
                               "installed or is disabled")
        return self._basemap
    
    def cartopy(self):
        """Return a cartopy.crs.Projection subclass for the
        projection"""
        if not cartopy_enabled():
            raise RuntimeError("'cartopy' is not "
                               "installed or is disabled")
        return self._cartopy
    
    def pyngl(self):
        """Return the PyNGL resources for the projection"""
        if not pyngl_enabled():
            raise RuntimeError("'pyngl' is not "
                               "installed or is disabled")
        return self._pyngl
    
    def proj4(self):
        """Return the proj4 string for the map projection"""
        return self._proj4
    
    def cf(self):
        """Return a dictionary with the NetCDF CF parameters for the 
        projection"""
        return self._cf_params
    
    
class LambertConformal(WrfProj):
    def __init__(self, bottom_left=None, top_right=None, 
                 lats=None, lons=None, **proj_params):
        super(LambertConformal, self).__init__(bottom_left, 
                    top_right, lats, lons, **proj_params)
        
        self._std_parallels = [self.truelat1]
        if self.truelat2 is not None:
            self._std_parallels.append(self.truelat2)
            
    @property
    def _cf_params(self):
        _cf_params = {}
        _cf_params["grid_mapping_name"] = "lambert_conformal_conic";
        _cf_params["standard_parallel"] = self._std_parallels
        _cf_params["longitude_of_central_meridian"] = self.stand_lon
        _cf_params["latitude_of_projection_origin"] = self.moad_cen_lat
        _cf_params["semi_major_axis"] = Constants.WRF_EARTH_RADIUS
        
        return _cf_params
    
    @property
    def _pyngl(self):
        if not pyngl_enabled():
            return None
        
        truelat2 = (self.truelat1 
                if _ismissing(self.truelat2) 
                else self.truelat2)
                                                
        _pyngl = Resources()
        _pyngl.mpProjection = "LambertConformal"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = self.ll_lon
        _pyngl.mpLeftCornerLatF = self.ll_lat
        _pyngl.mpRightCornerLonF = self.ur_lon
        _pyngl.mpRightCornerLatF = self.ur_lat
        _pyngl.mpLambertMeridianF = self.stand_lon
        _pyngl.mpLambertParallel1F = self.truelat1
        _pyngl.mpLambertParallel2F = truelat2
        
        return _pyngl
    
    @property
    def _basemap(self):
        if not basemap_enabled():
            return None
        
        _basemap = Basemap(projection = "lcc",
            lon_0 = self.stand_lon,
            lat_0 = self.moad_cen_lat,
            lat_1 = self.truelat1,
            lat_2 = self.truelat2,
            llcrnrlat = self.ll_lat,
            urcrnrlat = self.ur_lat,
            llcrnrlon = self.ll_lon,
            urcrnrlon = self.ur_lon,
            rsphere = Constants.WRF_EARTH_RADIUS,
            resolution = 'l')
        
        return _basemap
    
    @property
    def _cartopy(self):
        if not cartopy_enabled():
            return None
            
        _cartopy = crs.LambertConformal(
            central_longitude = self.stand_lon,
            central_latitude = self.moad_cen_lat,
            standard_parallels = self._std_parallels,
            globe = self._globe)
        
        return _cartopy
            
    @property
    def _cart_extents(self):
        # Need to modify the extents for the new projection
        pc = crs.PlateCarree()
        xs, ys, zs  = self._cartopy.transform_points(pc,
                             np.array([self.ll_lon, self.ur_lon]),
                             np.array([self.ll_lat, self.ur_lat])).T

                           
        _xlimits = xs.tolist()
        _ylimits = ys.tolist()
        
        return (_xlimits, _ylimits)
    
    @property
    def _proj4(self):
        truelat2 = (self.truelat1 
                    if _ismissing(self.truelat2) 
                    else self.truelat2)
        
        _proj4 = ("+proj=lcc +units=meters +a={} +b={} +lat_1={} "
                       "+lat_2={} +lat_0={} +lon_0={}".format(
                                            Constants.WRF_EARTH_RADIUS,
                                            Constants.WRF_EARTH_RADIUS,
                                            self.truelat1,
                                            truelat2, 
                                            self.moad_cen_lat,
                                            self.stand_lon))
        return _proj4
            
class Mercator(WrfProj):
    def __init__(self, bottom_left=None, top_right=None, 
                 lats=None, lons=None, **proj_params):
        super(Mercator, self).__init__(bottom_left, top_right, 
                                           lats, lons, **proj_params)
        
        self._lat_ts = (None 
            if self.truelat1 == 0. or _ismissing(self.truelat1) 
            else self.truelat1) 
        
    @property
    def _cf_params(self):
        
        _cf_params = {}
        _cf_params["grid_mapping_name"] = "mercator"
        _cf_params["longitude_of_projection_origin"] = self.stand_lon
        _cf_params["standard_parallel"] = self.truelat1
        
        return _cf_params
    
    @property
    def _pyngl(self):
        if not pyngl_enabled():
            return None
        
        _pyngl = Resources()
        _pyngl.mpProjection = "Mercator"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = self.ll_lon
        _pyngl.mpLeftCornerLatF = self.ll_lat
        _pyngl.mpRightCornerLonF = self.ur_lon
        _pyngl.mpRightCornerLatF = self.ur_lat
        _pyngl.mpCenterLatF = 0.0
        _pyngl.mpCenterLonF = self.stand_lon
        
        return _pyngl
    
    @property
    def _basemap(self):
        if not basemap_enabled():
            return None
                  
        _basemap = Basemap(projection = "merc",
                lon_0 = self.stand_lon,
                lat_0 = self.moad_cen_lat,
                lat_ts = self._lat_ts,
                llcrnrlat = self.ll_lat,
                urcrnrlat = self.ur_lat,
                llcrnrlon = self.ll_lon,
                urcrnrlon = self.ur_lon,
                rsphere = Constants.WRF_EARTH_RADIUS,
                resolution = 'l')
        
        return _basemap
    
    @property
    def _cartopy(self):
        if not cartopy_enabled():
            return None
        
        if self._lat_ts == 0.0:
            _cartopy = crs.Mercator(
                                central_longitude = self.stand_lon,
                                globe = self._globe)
        
        else:
            _cartopy = MercatorWithLatTS(
                central_longitude = self.stand_lon,
                latitude_true_scale = self._lat_ts,
                globe = self._globe)
            
        return _cartopy
    
    @property
    def _cart_extents(self):
                
        # Need to modify the extents for the new projection
        pc = crs.PlateCarree()
        xs, ys, zs  = self._cartopy.transform_points(pc,
                             np.array([self.ll_lon, self.ur_lon]),
                             np.array([self.ll_lat, self.ur_lat])).T
                            
        _xlimits = xs.tolist()
        _ylimits = ys.tolist()
        
        return (_xlimits, _ylimits)
    
    @property
    def _proj4(self):
        
        _proj4 = ("+proj=merc +units=meters +a={} +b={} "
                       "+lon_0={} +lat_ts={}".format(
                                            Constants.WRF_EARTH_RADIUS,
                                            Constants.WRF_EARTH_RADIUS,
                                            self.stand_lon,
                                            self._lat_ts))
        
        return _proj4
        
class PolarStereographic(WrfProj):
    def __init__(self, bottom_left=None, top_right=None, 
                 lats=None, lons=None, **proj_params):
        super(PolarStereographic, self).__init__(bottom_left, 
                        top_right, lats, lons, **proj_params)
        self._hemi = -90. if self.truelat1 < 0 else 90.
        self._lat_ts = (None 
                  if _ismissing(self.truelat1) 
                  else self.truelat1)
    
    @property
    def _cf_params(self):
        _cf_params = {}
        _cf_params["grid_mapping_name"] = "polar_stereographic"
        _cf_params["straight_vertical_longitude_from_pole"] = (
                                                               self.stand_lon)
        _cf_params["standard_parallel"] = self.truelat1
        _cf_params["latitude_of_projection_origin"] = self._hemi
        
        return _cf_params
    
    @property
    def _pyngl(self):
        if not pyngl_enabled():
            return None
        
        _pyngl = Resources()
        _pyngl.mpProjection = "Stereographic"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = self.ll_lon
        _pyngl.mpLeftCornerLatF = self.ll_lat
        _pyngl.mpRightCornerLonF = self.ur_lon
        _pyngl.mpRightCornerLatF = self.ur_lat
        
        _pyngl.mpCenterLonF = self.stand_lon
        if self._hemi > 0:
            _pyngl.mpCenterLatF = 90.0
        else:
            _pyngl.mpCenterLatF = -90.0
        
        return _pyngl
    
    @property
    def _basemap(self):
        if not basemap_enabled():
            return None

        _basemap = Basemap(projection = "stere",
            lon_0 = self.stand_lon,
            lat_0 = self._hemi,
            lat_ts = self._lat_ts,
            llcrnrlat = self.ll_lat,
            urcrnrlat = self.ur_lat,
            llcrnrlon = self.ll_lon,
            urcrnrlon = self.ur_lon,
            rsphere = Constants.WRF_EARTH_RADIUS,
            resolution = 'l')
        
        return _basemap
    
    @property
    def _cartopy(self):
        if not cartopy_enabled():
            return None
        
        _cartopy = crs.Stereographic(central_latitude=self._hemi, 
                                          central_longitude=self.stand_lon, 
                                          true_scale_latitude=self._lat_ts, 
                                          globe=self._globe)
        return _cartopy
    
    @property
    def _cart_extents(self):
        # Need to modify the extents for the new projection
        pc = crs.PlateCarree()
        xs, ys, zs  = self._cartopy.transform_points(pc,
                             np.array([self.ll_lon, self.ur_lon]),
                             np.array([self.ll_lat, self.ur_lat])).T
                            
        _xlimits = xs.tolist()
        _ylimits = ys.tolist()
        
        return (_xlimits, _ylimits)
    
    @property
    def _proj4(self):
        _proj4 = ("+proj=stere +units=meters +a={} +b={} "
                       "+lat0={} +lon_0={} +lat_ts={}".format(
                                            Constants.WRF_EARTH_RADIUS,
                                            Constants.WRF_EARTH_RADIUS,
                                            self._hemi,
                                            self.stand_lon,
                                            self._lat_ts))
        
        return _proj4
            
                  

class LatLon(WrfProj):
    def __init__(self, bottom_left=None, top_right=None, 
                 lats=None, lons=None, **proj_params):
        super(LatLon, self).__init__(bottom_left, top_right, 
                                         lats, lons, **proj_params)
    
    @property
    def _cf_params(self):
        _cf_params = {}
        _cf_params["grid_mapping_name"] = "latitude_longitude"
        return _cf_params
    
    @property
    def _pyngl(self):
        if not pyngl_enabled():
            return None
        
        _pyngl = Resources()
        _pyngl.mpProjection = "CylindricalEquidistant"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = self.ll_lon
        _pyngl.mpLeftCornerLatF = self.ll_lat
        _pyngl.mpRightCornerLonF = self.ur_lon
        _pyngl.mpRightCornerLatF = self.ur_lat
        _pyngl.mpCenterLonF = self.stand_lon
        _pyngl.mpCenterLatF = self.moad_cen_lat
        
        return _pyngl
    
    @property
    def _basemap(self):
        if not basemap_enabled():
            return None
        
        _basemap = Basemap(projection = "cyl",
            lon_0 = self.stand_lon,
            lat_0 = self.moad_cen_lat,
            llcrnrlat = self.ll_lat,
            urcrnrlat = self.ur_lat,
            llcrnrlon = self.ll_lon,
            urcrnrlon = self.ur_lon,
            rsphere = Constants.WRF_EARTH_RADIUS,
            resolution = 'l')
        
        return _basemap
    
    @property
    def _cartopy(self):
        if not cartopy_enabled():
            return None
        
        _cartopy = crs.PlateCarree(central_longitude=self.stand_lon,
                                            globe=self._globe)
        
        return _cartopy
    
    @property
    def _cart_extents(self):
        return ([self.ll_lon, self.ur_lon], [self.ll_lat, self.ur_lat])
    
    @property
    def _proj4(self):
        _proj4 = ("+proj=eqc +units=meters +a={} +b={} "
                       "+lon_0={}".format(Constants.WRF_EARTH_RADIUS,
                                          Constants.WRF_EARTH_RADIUS,
                                          self.stand_lon))
        return _proj4

# Notes (may not be correct since this projection confuses me):
# Each projection system handles this differently.
# 1) In WRF, if following the WPS instructions, POLE_LON is mainly used to 
#    determine north or south hemisphere.  In other words, it determines if 
#    the globe is tipped toward or away from you.  If a non-0 or non-180 
#    value is chosen, PyNGL cannot plot it.  
# 2) In WRF, POLE_LAT is always positive, but should be negative in the 
#    proj4 based systems when using the southern hemisphere projections.
# 3) In cartopy, pole_longitude is used to describe the dateline, which 
#    is 180 degrees away from the normal central (standard) longitude 
#    (e.g. center of the projection), according to the cartopy developer.  
# 4) In basemap, lon_0 should be set to the central (standard) longitude.
# 5) In either cartopy, basemap or pyngl, I'm not sure that projections with
#    a pole_lon not equal to 0 or 180 can be plotted.  Hopefully people 
#    follow the WPS instructions, otherwise I need to see a sample file and 
#    a lot of rum.
# 6) For items in 3 - 4, the "longitude" (lon_0 or pole_longitude) is 
#    determined by WRF's 
#    STAND_LON values, with the following calculations based on hemisphere:
#    BASEMAP:  NH:  -STAND_LON;  SH:  180.0 - STAND_LON
#    CARTOPY:  NH:  -STAND_LON - 180.; SH:  -STAND_LON
# 9) For PYNGL/NCL, you only need to set the center lat and center lon,
#    Center lat is the offset of the pole from +/- 90 degrees.  Center 
#    lon is -STAND_LON in NH and 180.0 - STAND_LON in SH.  
# 10) It also appears that NetCDF CF has no clear documentation on what 
#    each parameter means.  Going to assume it is the same as basemap, since
#    basemap appears to mirror the WMO way of doing things (tilt earth, then 
#    spin globe).
# 11) Basemap and cartopy produce projections that differ in their extent 
#     calculations by either using negative values or 0-360 (basemap).  For
#     this reason, the proj4 string for this class will use cartopy's values
#     to keep things in the -180 to 180, -90 to 90 range.
# 12) This projection makes me sad.           
class RotatedLatLon(WrfProj):
    def __init__(self, bottom_left=None, top_right=None, 
                 lats=None, lons=None, **proj_params):
        super(RotatedLatLon, self).__init__(bottom_left, top_right, 
                                    lats, lons, **proj_params)
        
        # Need to determine hemisphere, typically pole_lon is 0 for southern
        # hemisphere, 180 for northern hemisphere.  If not, going to have 
        # to guess based on other parameters, but hopefully people follow
        # the WPS instructions and this never happens.
        self._north = True
        if self.pole_lon is not None:
            if self.pole_lon == 0.:
                self._north = False
            elif self.pole_lon != 180.:
                if self.moad_cen_lat is not None and self.moad_cen_lat < 0.0:
                    # Only probably true
                    self._north = False
        else:
            if self.moad_cen_lat is not None and self.moad_cen_lat < 0.0:
                # Only probably true
                self._north = False 
        
        if self.pole_lat is not None and self.stand_lon is not None:
            self._pyngl_cen_lat = (90. - self.pole_lat if self._north 
                             else self.pole_lat - 90.0)
            self._pyngl_cen_lon = (-self.stand_lon if self._north 
                             else 180.0 - self.stand_lon)
            self._bm_lon_0 = (-self.stand_lon if self._north 
                              else 180.0 - self.stand_lon)
            self._bm_cart_pole_lat = (self.pole_lat if self._north 
                                      else -self.pole_lat )
            # The important point is that pole longitude is the position 
            # of the dateline of the new projection, not its central 
            # longitude (per the creator of cartopy).  This is based on 
            # how it's handled by agencies like WMO, but not proj4.
            self._cart_pole_lon = (-self.stand_lon - 180.0 if self._north 
                             else -self.stand_lon)
        else:
            self._pyngl_cen_lat = self.moad_cen_lat
            self._pyngl_cen_lon = self.stand_lon
            self._bm_cart_pole_lat = (90.0 - self.moad_cen_lat if self._north 
                                else -90.0 - self.moad_cen_lat)
            self._bm_lon_0 = (-self.stand_lon if self._north 
                              else 180.0 - self.stand_lon)
            self._cart_pole_lon = (-self.stand_lon - 180.0 if self._north 
                             else -self.stand_lon)
        
    @property
    def _cf_params(self):   
        _cf_params = {}
        # Assuming this follows the same guidelines as cartopy
        _cf_params["grid_mapping_name"] = "rotated_latitude_longitude"
        _cf_params["grid_north_pole_latitude"] = self._bm_cart_pole_lat
        _cf_params["grid_north_pole_longitude"] = self.pole_lon
        _cf_params["north_pole_grid_longitude"] = self._bm_lon_0
        
        return _cf_params
    
    @property
    def _pyngl(self):
        if not pyngl_enabled():
            return None
        
        _pyngl = Resources()
        _pyngl.mpProjection = "CylindricalEquidistant"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = self.ll_lon
        _pyngl.mpLeftCornerLatF = self.ll_lat
        _pyngl.mpRightCornerLonF = self.ur_lon
        _pyngl.mpRightCornerLatF = self.ur_lat
        _pyngl.mpCenterLatF = self._pyngl_cen_lat
        _pyngl.mpCenterLonF = self._pyngl_cen_lon
        
        return _pyngl
    
    @property
    def _basemap(self):
        if not basemap_enabled():
            return None
        
        _basemap = Basemap(projection = "rotpole",
                                o_lat_p = self._bm_cart_pole_lat,
                                o_lon_p = self.pole_lon,
                                llcrnrlat = self.ll_lat,
                                urcrnrlat = self.ur_lat,
                                llcrnrlon = self.ll_lon,
                                urcrnrlon = self.ur_lon,
                                lon_0 = self._bm_lon_0,
                                rsphere = Constants.WRF_EARTH_RADIUS,
                                resolution = 'l')
        return _basemap
    
    @property
    def _cartopy(self):
        if not cartopy_enabled():
            return None
        
        _cartopy = crs.RotatedPole(
                                pole_longitude=self._cart_pole_lon, 
                                pole_latitude=self._bm_cart_pole_lat, 
                                central_rotated_longitude=(
                                        180.0 - self.pole_lon), # Probably
                                globe = self._globe)
        return _cartopy
    
    @property
    def _cart_extents(self):
        # Need to modify the extents for the new projection
        pc = crs.PlateCarree()
        xs, ys, zs  = self._cartopy.transform_points(pc,
                             np.array([self.ll_lon, self.ur_lon]),
                             np.array([self.ll_lat, self.ur_lat])).T
                            
        _xlimits = xs.tolist()
        _ylimits = ys.tolist()
        
        return (_xlimits, _ylimits)
    
    @property
    def _proj4(self):
        _proj4 = ("+proj=ob_tran +o_proj=latlon "
                       "+a={} +b={} +to_meter={} +o_lon_p={} +o_lat_p={} "
                       "+lon_0={}".format(Constants.WRF_EARTH_RADIUS,
                                          Constants.WRF_EARTH_RADIUS,
                                          math.radians(1),
                                          180.0 - self.pole_lon,
                                          self._bm_cart_pole_lat,
                                          180.0 + self._cart_pole_lon))
        
        return _proj4
        
def getproj(bottom_left=None, top_right=None, 
                 lats=None, lons=None, **proj_params):
    
    proj_type = proj_params.get("MAP_PROJ", 0)
    if proj_type == 1:
        return LambertConformal(bottom_left, top_right, 
                                    lats, lons, **proj_params)
    elif proj_type == 2:
        return PolarStereographic(bottom_left, top_right, 
                                      lats, lons, **proj_params)
    elif proj_type == 3:
        return Mercator(bottom_left, top_right, 
                            lats, lons, **proj_params)
    elif proj_type == 0 or proj_type == 6:
        if (proj_params.get("POLE_LAT", None) == 90. 
            and proj_params.get("POLE_LON", None) == 0.):
            return LatLon(bottom_left, top_right, 
                              lats, lons, **proj_params)
        else:
            return RotatedLatLon(bottom_left, top_right, 
                                 lats, lons, **proj_params)
    else:
        # Unknown projection
        return WrfProj(bottom_left, top_right, 
                       lats, lons, **proj_params)
    
        