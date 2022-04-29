from __future__ import (absolute_import, division, print_function)
import numpy as np
import math
from decimal import Decimal, Context, ROUND_HALF_UP

from .config import basemap_enabled, cartopy_enabled, pyngl_enabled
from .constants import Constants, ProjectionTypes
from .projutils import dict_keys_to_upper
from .py3compat import viewitems

if cartopy_enabled():
    from cartopy import crs

if basemap_enabled():
    from mpl_toolkits.basemap import Basemap

if pyngl_enabled():
    from Ngl import Resources



if cartopy_enabled():
    class MercatorWithLatTS(crs.Mercator):
        """A :class:`cartopy.crs.Mercator` subclass that adds support for
        a latitude of true scale parameter.

        See Also:

            :class:`cartopy.crs.Mercator`

        """
        def __init__(self, central_longitude=0.0,
                     latitude_true_scale=0.0,
                     min_latitude=-80.0,
                     max_latitude=84.0,
                     globe=None):
            """Initialize a :class:`wrf.MercatorWithLatTS` object.

            Args:

                central_longitude (:obj:`float`, optional): The central
                    longitude.  Default is 0.0.

                latitude_true_scale (:obj:`float`, optional): The latitude
                    of true scale.  Default is 0.0.

                min_latitude (:obj:`float`, optional): The maximum southerly
                    extent of the projection.  Default is -80.0.

                max_latitude (:obj:`float`, optional): The maximum northerly
                    extent of the projection.  Default is 84.0.

                globe (:class:`cartopy.crs.Globe`, optional): A globe object.
                    If omitted, a default globe is created.

            """
            proj4_params = [("proj", "merc"),
                            ("lon_0", central_longitude),
                            ("lat_ts", latitude_true_scale),
                            ("k", 1),
                            ("units", "m")]
            super(crs.Mercator, self).__init__(proj4_params, globe=globe)

            # Need to have x/y limits defined for the initial hash which
            # gets used within transform_points for caching
            self._x_limits = self._y_limits = None

            # Calculate limits.
            limits = self.transform_points(
                crs.Geodetic(),
                np.array([-180, 180]) + central_longitude,
                np.array([min_latitude, max_latitude]))

            # When using a latitude of true scale, the min/max x-limits get set
            # to the same value, so make sure the left one is negative
            xlimits = limits[..., 0]

            if math.fabs(xlimits[0] - xlimits[1]) < 1e-6:
                if xlimits[0] < 0:
                    xlimits[1] = -xlimits[1]
                else:
                    xlimits[0] = -xlimits[0]

            # Compatibility with cartopy >= 0.17
            self._x_limits = tuple(xlimits)
            self._y_limits =  tuple(limits[..., 1])

            self._threshold = min(np.diff(self.x_limits)[0] / 720,
                                  np.diff(self.y_limits)[0] / 360)


def _ismissing(val, islat=True):
    """Return True if a value is None or out of bounds.

    This function is used to check for invalid latitude/longitude values.

    Args:

        val (numeric): A numeric value.

        islat (:obj:`bool`): Set to False if checking for longitude values.

    Returns:

        :obj:`bool`: True if the value is None, or an out of bounds value.

    """
    if islat:
        if val is None:
            return True

        if math.fabs(val) > 90.:
            return True
    else:
        if val is None:
            return True

        if math.fabs(val) > 360.:
            return True

    return False


class WrfProj(object):
    """A base class for storing map projection information from WRF data.

    Subclasses of this type will be stored in the 'projection' attribute
    entry within a :attr:`xarray.DataArray.attrs` dictionary.  This base class
    contains the methods required to extract the appropriate projection class
    for PyNGL, matplotlib basemap, and cartopy.

    Attributes:

        map_proj (:obj:`int`): Model projection integer id.

        truelat1 (:obj:`float`): True latitude 1.

        truelat2 (:obj:`float`): True latitude 2.

        moad_cen_lat (:obj:`float`): Mother of all domains center latitude.

        stand_lon (:obj:`float`): Standard longitude.

        pole_lat (:obj:`float`): The pole latitude.

        pole_lon (:obj:`float`): The pole longitude.

        dx (:obj:`float`): The x grid spacing.

        dy (:obj:`float`): The y grid spacing.


    """
    def __init__(self, **proj_params):
        """Initialize a :class:`wrf.WrfProj` object.

        Args:

            **proj_params:  Map projection optional keyword arguments, that
                have the same names as found in WRF output NetCDF global
                attributes (case insensitive):

                - 'MAP_PROJ': The map projection type as an integer.
                - 'TRUELAT1': True latitude 1.
                - 'TRUELAT2': True latitude 2.
                - 'MOAD_CEN_LAT': Mother of all domains center latitude.
                - 'STAND_LON': Standard longitude.
                - 'POLE_LAT': Pole latitude.
                - 'POLE_LON': Pole longitude.

        """

        up_proj_params = dict_keys_to_upper(proj_params)

        self.map_proj = up_proj_params.get("MAP_PROJ", None)

        # These indicate the center of the nest/domain, not necessarily the
        # center of the projection
        self._cen_lat = up_proj_params.get("CEN_LAT", None)
        self._cen_lon = up_proj_params.get("CEN_LON", None)

        self.truelat1 = up_proj_params.get("TRUELAT1", None)
        self.truelat2 = (up_proj_params.get("TRUELAT2", None)
                         if not _ismissing(up_proj_params.get("TRUELAT2",
                                                              None))
                         else None)
        self.moad_cen_lat = up_proj_params.get("MOAD_CEN_LAT", None)
        self.stand_lon = up_proj_params.get("STAND_LON", None)
        self.pole_lat = up_proj_params.get("POLE_LAT", None)
        self.pole_lon = up_proj_params.get("POLE_LON", None)

        self.dx = up_proj_params.get("DX", None)
        self.dy = up_proj_params.get("DY", None)

        # Just in case...
        if self.moad_cen_lat is None:
            self.moad_cen_lat = self._cen_lat

        if self.stand_lon is None:
            self.stand_lon = self._cen_lon

    @staticmethod
    def _context_equal(x, y, ctx):
        """Return True if both objects are equal based on the provided context.

        Args:

            x (numeric): A numeric value.

            y (numeric): A numeric value.

            ctx (:class:`decimal.Context`): A decimal Context object.

        Returns:

            :obj:`bool`: True if the values are equal based on the provided
            context, False otherwise.

        """
        if x is not None:
            if y is None:
                return False

            # Note:  The float conversion is because these may come in as
            # numpy.float32 or numpy.float64, which Decimal does not know
            # how to handle.
            if (Decimal(float(x)).normalize(ctx) !=
                    Decimal(float(y)).normalize(ctx)):
                return False
        else:
            if y is not None:
                return False

        return True

    def __eq__(self, other):
        """Return True if this projection object is the same as *other*.

        Note:  WRF can use either floats or doubles, so only going to
        guarantee floating point precision equivalence, in case the different
        types are ever compared (or a projection is hand constructed). For WRF
        domains, 7-digit equivalence should be good enough.

        """

        if self.map_proj is not None:
            if self.map_proj != other.map_proj:
                return False
        else:
            if other.map_proj is not None:
                return False

        # Only checking for floating point equal.
        ctx = Context(prec=7, rounding=ROUND_HALF_UP)

        return (WrfProj._context_equal(self.truelat1, other.truelat1, ctx) and
                WrfProj._context_equal(self.truelat2, other.truelat2, ctx) and
                WrfProj._context_equal(self.moad_cen_lat, other.moad_cen_lat,
                                       ctx) and
                WrfProj._context_equal(self.stand_lon, other.stand_lon,
                                       ctx) and
                WrfProj._context_equal(self.pole_lat, other.pole_lat, ctx) and
                WrfProj._context_equal(self.pole_lon, other.pole_lon, ctx) and
                WrfProj._context_equal(self.dx, other.dx, ctx) and
                WrfProj._context_equal(self.dy, other.dy, ctx))

    def _basemap(self, geobounds, **kwargs):
        return None

    def _cf_params(self):
        return None

    def _cartopy(self):
        return None

    def _calc_extents(self, geobounds):
        # Need to modify the extents for the new projection
        pc = crs.PlateCarree()
        xs, ys, _ = self._cartopy().transform_points(
            pc,
            np.array([geobounds.bottom_left.lon, geobounds.top_right.lon]),
            np.array([geobounds.bottom_left.lat, geobounds.top_right.lat])).T
        _xlimits = xs.tolist()
        _ylimits = ys.tolist()

        return (_xlimits, _ylimits)

    def _cart_extents(self, geobounds):
        try:
            _ = len(geobounds)
        except TypeError:  # Only a single object
            extents = self._calc_extents(geobounds)
        else:
            extents = np.empty(geobounds.shape, np.object)

            for idxs, geobnd_val in np.ndenumerate(geobounds):
                extents[idxs] = self._calc_extents(geobnd_val)

        return extents

    def _pyngl(self, geobounds):
        return None

    def _proj4(self):
        return None

    def _globe(self):
        return (None if not cartopy_enabled()
                else crs.Globe(ellipse=None,
                               semimajor_axis=Constants.WRF_EARTH_RADIUS,
                               semiminor_axis=Constants.WRF_EARTH_RADIUS,
                               nadgrids="@null"))

    def cartopy_xlim(self, geobounds):
        """Return the x extents in projected coordinates for cartopy.

        Returns:

            :obj:`list`: A pair of [xmin, xmax].

        See Also:

            :mod:`cartopy`, :mod:`matplotlib`

        """
        try:
            _ = len(geobounds)
        except TypeError:
            x_extents = self._cart_extents(geobounds)[0]
        else:
            extents = self._cart_extents(geobounds)
            x_extents = np.empty(extents.shape, np.object)

            for idxs, extent in np.ndenumerate(extents):
                x_extents[idxs] = extent[0]

        return x_extents

    def cartopy_ylim(self, geobounds):
        """Return the y extents in projected coordinates for cartopy.

        Returns:

            :obj:`list`: A pair of [ymin, ymax].

        See Also:

            :mod:`cartopy`, :mod:`matplotlib`

        """
        try:
            _ = len(geobounds)
        except TypeError:
            y_extents = self._cart_extents(geobounds)[1]
        else:
            extents = self._cart_extents(geobounds)
            y_extents = np.empty(extents.shape, np.object)

            for idxs, extent in np.ndenumerate(extents):
                y_extents[idxs] = extent[1]

        return y_extents

    def __repr__(self):
        args = ("stand_lon={}, moad_cen_lat={}, "
                "truelat1={}, truelat2={}, "
                "pole_lat={}, pole_lon={}".format(self.stand_lon,
                                                  self.moad_cen_lat,
                                                  self.truelat1,
                                                  self.truelat2,
                                                  self.pole_lat,
                                                  self.pole_lon))
        return "{}({})".format(self.__class__.__name__, args)

    def basemap(self, geobounds, **kwargs):
        """Return a :class:`matplotlib.mpl_toolkits.basemap.Basemap` object
        for the map projection.

        Args:

            geobounds (:class:`wrf.GeoBounds`, optional):  The geobounds to
                get the extents.  If set to None and using the *var* parameter,
                the geobounds will be taken from the variable.  If using a
                file, then the geobounds will be taken from the native grid.

            **kwargs: Keyword arguments for creating a
                :class:`matplotlib.mpl_toolkits.basemap.Basemap`.  By default,
                the domain bounds will be set to the native projection, the
                resolution will be set to 'l', and the other projection
                parameters will be set by the information in the file.

        Returns:

            :class:`matplotlib.mpl_toolkits.basemap.Basemap`: A Basemap
            object for the projection.

        See Also:

            :class:`matplotlib.mpl_toolkits.basemap.Basemap`

        """
        if not basemap_enabled():
            raise RuntimeError("'mpl_toolkits.basemap' is not "
                               "installed or is disabled")

        return self._basemap(geobounds, **kwargs)

    def cartopy(self):
        """Return a :class:`cartopy.crs.Projection` subclass for the
        map projection.

        Returns:

            :class:`cartopy.crs.Projection`: A Projection subclass for the
            map projection.

        See Also:

            :class:`cartopy.crs.Projection`

        """
        if not cartopy_enabled():
            raise RuntimeError("'cartopy' is not "
                               "installed or is disabled")
        return self._cartopy()

    def pyngl(self, geobounds, **kwargs):
        """Return a :class:`Ngl.Resources` object for the map projection.

        Args:

            geobounds (:class:`wrf.GeoBounds`, optional):  The geobounds to
                get the extents.  If set to None and using the *var* parameter,
                the geobounds will be taken from the variable.  If using a
                file, then the geobounds will be taken from the native grid.

            **kwargs:  Additional PyNGL resources to set while creating the
                :class:`Ngl.Resources` object.

        Returns:

            :class:`Ngl.Resources`: A dict-like object that contains the
            PyNGL resources for the map projection.

        See Also:

            `PyNGL <https://www.pyngl.ucar.edu/>`_

        """
        if not pyngl_enabled():
            raise RuntimeError("'pyngl' is not "
                               "installed or is disabled")
        return self._pyngl(geobounds, **kwargs)

    def proj4(self):
        """Return the PROJ.4 string for the map projection.

        Returns:

            :obj:`str`: A string suitable for use with the PROJ.4 library.

        See Also:

            PROJ.4 <https://trac.osgeo.org/proj/>`_

        """
        return self._proj4()

    def cf(self):
        """Return a dictionary with the NetCDF CF parameters for the
        projection.

        Returns:

        :obj:`dict`: A dictionary with the NetCDF CF parameter names and
        projection parameter values.

        """
        return self._cf_params()


# Used for 'missing' projection values during the 'join' method
class NullProjection(WrfProj):
    """A :class:`wrf.WrfProj` subclass for empty projections.

    The :class:`NullProjection` is primarily used for creating missing
    projections when using the 'join' method.

    """
    def __init__(self):
        """Initialize a :class:`wrf.NullProjection` object."""
        pass

    def __repr__(self):
        return "{}()".format(self.__class__.__name__)


class LambertConformal(WrfProj):
    """A :class:`wrf.WrfProj` subclass for Lambert Conformal Conic projections.

    See Also:

        :class:`wrf.WrfProj`, :class:`wrf.LatLon`,
        :class:`wrf.PolarStereographic`,
        :class:`Mercator`, :class:`RotatedLatLon`

    """
    def __init__(self, **proj_params):
        """Initialize a :class:`wrf.LambertConformal` object.

        Args:

            **proj_params:  Map projection optional keyword arguments, that
                have the same names as found in WRF output NetCDF global
                attributes:

                - 'TRUELAT1': True latitude 1.
                - 'TRUELAT2': True latitude 2.
                - 'MOAD_CEN_LAT': Mother of all domains center latitude.
                - 'STAND_LON': Standard longitude.
                - 'POLE_LAT': Pole latitude.
                - 'POLE_LON': Pole longitude.

        """
        super(LambertConformal, self).__init__(**proj_params)

        self._std_parallels = [self.truelat1]
        if self.truelat2 is not None:
            self._std_parallels.append(self.truelat2)

    def _cf_params(self):
        _cf_params = {}
        _cf_params["grid_mapping_name"] = "lambert_conformal_conic"
        _cf_params["standard_parallel"] = self._std_parallels
        _cf_params["longitude_of_central_meridian"] = self.stand_lon
        _cf_params["latitude_of_projection_origin"] = self.moad_cen_lat
        _cf_params["earth_radius"] = Constants.WRF_EARTH_RADIUS

        return _cf_params

    def _pyngl(self, geobounds, **kwargs):
        if not pyngl_enabled():
            return None

        truelat2 = (self.truelat1 if _ismissing(self.truelat2)
                    else self.truelat2)

        _pyngl = Resources()
        _pyngl.mpProjection = "LambertConformal"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpLambertMeridianF = self.stand_lon
        _pyngl.mpLambertParallel1F = self.truelat1
        _pyngl.mpLambertParallel2F = truelat2

        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = geobounds.bottom_left.lon
        _pyngl.mpLeftCornerLatF = geobounds.bottom_left.lat
        _pyngl.mpRightCornerLonF = geobounds.top_right.lon
        _pyngl.mpRightCornerLatF = geobounds.top_right.lat

        for key, val in viewitems(kwargs):
            setattr(_pyngl, key, val)

        return _pyngl

    def _basemap(self, geobounds, **kwargs):
        if not basemap_enabled():
            return None

        local_kwargs = dict(projection="lcc",
                            lon_0=self.stand_lon,
                            lat_0=self.moad_cen_lat,
                            lat_1=self.truelat1,
                            lat_2=self.truelat2,
                            llcrnrlat=geobounds.bottom_left.lat,
                            urcrnrlat=geobounds.top_right.lat,
                            llcrnrlon=geobounds.bottom_left.lon,
                            urcrnrlon=geobounds.top_right.lon,
                            rsphere=Constants.WRF_EARTH_RADIUS,
                            resolution='l')
        local_kwargs.update(kwargs)

        _basemap = Basemap(**local_kwargs)

        return _basemap

    def _cartopy(self):
        if not cartopy_enabled():
            return None

        # Set cutoff to -30 for NH, +30.0 for SH.
        cutoff = -30.0 if self.moad_cen_lat >= 0 else 30.0

        _cartopy = crs.LambertConformal(
            central_longitude=self.stand_lon,
            central_latitude=self.moad_cen_lat,
            standard_parallels=self._std_parallels,
            globe=self._globe(),
            cutoff=cutoff)

        return _cartopy

    def _proj4(self):
        truelat2 = (self.truelat1
                    if _ismissing(self.truelat2)
                    else self.truelat2)

        _proj4 = ("+proj=lcc +units=m +a={} +b={} +lat_1={} "
                  "+lat_2={} +lat_0={} +lon_0={} +nadgrids=@null".format(
                      Constants.WRF_EARTH_RADIUS,
                      Constants.WRF_EARTH_RADIUS,
                      self.truelat1,
                      truelat2,
                      self.moad_cen_lat,
                      self.stand_lon))
        return _proj4


class Mercator(WrfProj):
    """A :class:`wrf.WrfProj` subclass for Mercator projections.

    See Also:

        :class:`wrf.WrfProj`, :class:`wrf.LatLon`,
        :class:`wrf.PolarStereographic`,
        :class:`RotatedLatLon`, :class:`LambertConformal`

    """
    def __init__(self, **proj_params):
        """Initialize a :class:`wrf.Mercator` object.

        Args:

            **proj_params:  Map projection optional keyword arguments, that
                have the same names as found in WRF output NetCDF global
                attributes:

                - 'TRUELAT1': True latitude 1.
                - 'TRUELAT2': True latitude 2.
                - 'MOAD_CEN_LAT': Mother of all domains center latitude.
                - 'STAND_LON': Standard longitude.
                - 'POLE_LAT': Pole latitude.
                - 'POLE_LON': Pole longitude.

        """
        super(Mercator, self).__init__(**proj_params)

        self._lat_ts = (
            None if self.truelat1 == 0. or _ismissing(self.truelat1)
            else self.truelat1)

        self._stand_lon = (0. if _ismissing(self.stand_lon, islat=False)
                           else self.stand_lon)

    def _cf_params(self):
        _cf_params = {}
        _cf_params["grid_mapping_name"] = "mercator"
        _cf_params["longitude_of_projection_origin"] = self.stand_lon
        _cf_params["standard_parallel"] = self.truelat1
        _cf_params["earth_radius"] = Constants.WRF_EARTH_RADIUS

        return _cf_params

    def _pyngl(self, geobounds, **kwargs):
        if not pyngl_enabled():
            return None

        _pyngl = Resources()
        _pyngl.mpProjection = "Mercator"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpCenterLatF = 0.0
        _pyngl.mpCenterLonF = self._stand_lon

        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = geobounds.bottom_left.lon
        _pyngl.mpLeftCornerLatF = geobounds.bottom_left.lat
        _pyngl.mpRightCornerLonF = geobounds.top_right.lon
        _pyngl.mpRightCornerLatF = geobounds.top_right.lat

        for key, val in viewitems(kwargs):
            setattr(_pyngl, key, val)

        return _pyngl

    def _basemap(self, geobounds, **kwargs):
        if not basemap_enabled():
            return None

        local_kwargs = dict(projection="merc",
                            lon_0=self._stand_lon,
                            lat_0=self.moad_cen_lat,
                            lat_ts=self._lat_ts,
                            llcrnrlat=geobounds.bottom_left.lat,
                            urcrnrlat=geobounds.top_right.lat,
                            llcrnrlon=geobounds.bottom_left.lon,
                            urcrnrlon=geobounds.top_right.lon,
                            rsphere=Constants.WRF_EARTH_RADIUS,
                            resolution="l")
        local_kwargs.update(kwargs)

        _basemap = Basemap(**local_kwargs)

        return _basemap

    def _cartopy(self):
        if not cartopy_enabled():
            return None

        if self._lat_ts == 0.0:
            _cartopy = crs.Mercator(central_longitude=self._stand_lon,
                                    globe=self._globe())
        else:
            _cartopy = MercatorWithLatTS(central_longitude=self._stand_lon,
                                         latitude_true_scale=self._lat_ts,
                                         globe=self._globe())

        return _cartopy

    def _proj4(self):

        _proj4 = ("+proj=merc +units=m +a={} +b={} "
                  "+lon_0={} +lat_ts={} +nadgrids=@null".format(
                                            Constants.WRF_EARTH_RADIUS,
                                            Constants.WRF_EARTH_RADIUS,
                                            self._stand_lon,
                                            self._lat_ts))

        return _proj4


class PolarStereographic(WrfProj):
    """A :class:`wrf.WrfProj` subclass for Polar Stereographic projections.

    See Also:

        :class:`wrf.WrfProj`, :class:`wrf.LatLon`,
        :class:`wrf.RotatedLatLon`,
        :class:`Mercator`, :class:`LambertConformal`

    """
    def __init__(self, **proj_params):
        """Initialize a :class:`wrf.PolarStereographic` object.

        Args:

            **proj_params:  Map projection optional keyword arguments, that
                have the same names as found in WRF output NetCDF global
                attributes:

                - 'TRUELAT1': True latitude 1.
                - 'TRUELAT2': True latitude 2.
                - 'MOAD_CEN_LAT': Mother of all domains center latitude.
                - 'STAND_LON': Standard longitude.
                - 'POLE_LAT': Pole latitude.
                - 'POLE_LON': Pole longitude.

        """
        super(PolarStereographic, self).__init__(**proj_params)
        self._hemi = -90. if self.truelat1 < 0 else 90.
        self._lat_ts = (None if _ismissing(self.truelat1) else self.truelat1)

    def _cf_params(self):
        _cf_params = {}
        _cf_params["grid_mapping_name"] = "polar_stereographic"
        _cf_params["straight_vertical_longitude_from_pole"] = (
                                                               self.stand_lon)
        _cf_params["standard_parallel"] = self.truelat1
        _cf_params["latitude_of_projection_origin"] = self._hemi
        _cf_params["earth_radius"] = Constants.WRF_EARTH_RADIUS

        return _cf_params

    def _pyngl(self, geobounds, **kwargs):
        if not pyngl_enabled():
            return None

        _pyngl = Resources()
        _pyngl.mpProjection = "Stereographic"
        _pyngl.mpDataBaseVersion = "MediumRes"

        _pyngl.mpCenterLonF = self.stand_lon
        if self._hemi > 0:
            _pyngl.mpCenterLatF = 90.0
        else:
            _pyngl.mpCenterLatF = -90.0

        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = geobounds.bottom_left.lon
        _pyngl.mpLeftCornerLatF = geobounds.bottom_left.lat
        _pyngl.mpRightCornerLonF = geobounds.top_right.lon
        _pyngl.mpRightCornerLatF = geobounds.top_right.lat

        for key, val in viewitems(kwargs):
            setattr(_pyngl, key, val)

        return _pyngl

    def _basemap(self, geobounds, **kwargs):
        if not basemap_enabled():
            return None

        local_kwargs = dict(projection="stere",
                            lon_0=self.stand_lon,
                            lat_0=self._hemi,
                            lat_ts=self._lat_ts,
                            llcrnrlat=geobounds.bottom_left.lat,
                            urcrnrlat=geobounds.top_right.lat,
                            llcrnrlon=geobounds.bottom_left.lon,
                            urcrnrlon=geobounds.top_right.lon,
                            rsphere=Constants.WRF_EARTH_RADIUS,
                            resolution="l")
        local_kwargs.update(kwargs)

        _basemap = Basemap(**local_kwargs)

        return _basemap

    def _cartopy(self):
        if not cartopy_enabled():
            return None

        _cartopy = crs.Stereographic(central_latitude=self._hemi,
                                     central_longitude=self.stand_lon,
                                     true_scale_latitude=self._lat_ts,
                                     globe=self._globe())
        return _cartopy

    def _proj4(self):
        _proj4 = ("+proj=stere +units=m +a={} +b={} "
                  "+lat_0={} +lon_0={} +lat_ts={} +nadgrids=@null".format(
                                            Constants.WRF_EARTH_RADIUS,
                                            Constants.WRF_EARTH_RADIUS,
                                            self._hemi,
                                            self.stand_lon,
                                            self._lat_ts))

        return _proj4


class LatLon(WrfProj):
    """A :class:`wrf.WrfProj` subclass for Lat Lon projections.

    See Also:

        :class:`wrf.WrfProj`, :class:`wrf.RotatedLatLon`,
        :class:`wrf.PolarStereographic`,
        :class:`Mercator`, :class:`LambertConformal`

    """
    def __init__(self, **proj_params):
        """Initialize a :class:`wrf.LatLon` object.

        Args:

            **proj_params:  Map projection optional keyword arguments, that
                have the same names as found in WRF output NetCDF global
                attributes:

                - 'TRUELAT1': True latitude 1.
                - 'TRUELAT2': True latitude 2.
                - 'MOAD_CEN_LAT': Mother of all domains center latitude.
                - 'STAND_LON': Standard longitude.
                - 'POLE_LAT': Pole latitude.
                - 'POLE_LON': Pole longitude.

        """
        super(LatLon, self).__init__(**proj_params)

    def _cf_params(self):
        _cf_params = {}
        _cf_params["grid_mapping_name"] = "latitude_longitude"
        _cf_params["earth_radius"] = Constants.WRF_EARTH_RADIUS
        return _cf_params

    def _pyngl(self, geobounds, **kwargs):
        if not pyngl_enabled():
            return None

        _pyngl = Resources()
        _pyngl.mpProjection = "CylindricalEquidistant"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpCenterLonF = self.stand_lon
        _pyngl.mpCenterLatF = self.moad_cen_lat

        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = geobounds.bottom_left.lon
        _pyngl.mpLeftCornerLatF = geobounds.bottom_left.lat
        _pyngl.mpRightCornerLonF = geobounds.top_right.lon
        _pyngl.mpRightCornerLatF = geobounds.top_right.lat

        for key, val in viewitems(kwargs):
            setattr(_pyngl, key, val)

        return _pyngl

    def _basemap(self, geobounds, **kwargs):
        if not basemap_enabled():
            return None

        local_kwargs = dict(projection="cyl",
                            lon_0=self.stand_lon,
                            lat_0=self.moad_cen_lat,
                            llcrnrlat=geobounds.bottom_left.lat,
                            urcrnrlat=geobounds.top_right.lat,
                            llcrnrlon=geobounds.bottom_left.lon,
                            urcrnrlon=geobounds.top_right.lon,
                            rsphere=Constants.WRF_EARTH_RADIUS,
                            resolution="l")

        local_kwargs.update(kwargs)
        _basemap = Basemap(**local_kwargs)

        return _basemap

    def _cartopy(self):
        if not cartopy_enabled():
            return None

        _cartopy = crs.PlateCarree(central_longitude=self.stand_lon,
                                   globe=self._globe())

        return _cartopy

    def _cart_extents(self, geobounds):
        return ([geobounds.bottom_left.lon, geobounds.top_right.lon],
                [geobounds.bottom_left.lat, geobounds.top_right.lat])

    def _proj4(self):
        _proj4 = ("+proj=eqc +units=m +a={} +b={} "
                  "+lon_0={} +nadgrids=@null".format(
                      Constants.WRF_EARTH_RADIUS,
                      Constants.WRF_EARTH_RADIUS,
                      self.stand_lon))
        return _proj4


# Notes (may not be correct since this projection confuses me):
# Each projection system handles this differently.
# 1) In WRF, if following the WPS instructions, POLE_LON is mainly used to
#    determine north or south hemisphere.  In other words, it determines if
#    the globe is tipped toward or away from you.
# 2) In WRF, POLE_LAT is always positive, but should be negative in the
#    proj4 based systems when using the southern hemisphere projections.
# 3) In cartopy, pole_longitude is used to describe the dateline, which
#    is 180 degrees away from the normal central (standard) longitude
#    (e.g. center of the projection), according to the cartopy developer.
# 4) In basemap, lon_0 should be set to the central (standard) longitude.
# 5) In either cartopy, basemap or pyngl, I'm not sure that projections with
#    a pole_lon not equal to 0 or 180 can be plotted.  Hopefully people
#    follow the WPS instructions, otherwise I need to see a sample file.
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
    """A :class:`wrf.WrfProj` subclass for Rotated Lat Lon projections.

    See Also:

        :class:`wrf.WrfProj`, :class:`wrf.LatLon`,
        :class:`wrf.PolarStereographic`,
        :class:`Mercator`, :class:`LambertConformal`

    """
    def __init__(self, **proj_params):
        """Initialize a :class:`wrf.RotatedLatLon` object.

        Args:

            **proj_params:  Map projection optional keyword arguments, that
                have the same names as found in WRF output NetCDF global
                attributes:

                - 'TRUELAT1': True latitude 1.
                - 'TRUELAT2': True latitude 2.
                - 'MOAD_CEN_LAT': Mother of all domains center latitude.
                - 'STAND_LON': Standard longitude.
                - 'POLE_LAT': Pole latitude.
                - 'POLE_LON': Pole longitude.

        """
        super(RotatedLatLon, self).__init__(**proj_params)

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
                                      else -self.pole_lat)
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

    def _cf_params(self):
        _cf_params = {}
        # Assuming this follows the same guidelines as cartopy
        _cf_params["grid_mapping_name"] = "rotated_latitude_longitude"
        _cf_params["grid_north_pole_latitude"] = self._bm_cart_pole_lat
        _cf_params["grid_north_pole_longitude"] = self.pole_lon
        _cf_params["north_pole_grid_longitude"] = self._bm_lon_0
        _cf_params["earth_radius"] = Constants.WRF_EARTH_RADIUS

        return _cf_params

    def _pyngl(self, geobounds, **kwargs):
        if not pyngl_enabled():
            return None

        _pyngl = Resources()
        _pyngl.mpProjection = "CylindricalEquidistant"
        _pyngl.mpDataBaseVersion = "MediumRes"
        _pyngl.mpCenterLatF = self._pyngl_cen_lat
        _pyngl.mpCenterLonF = self._pyngl_cen_lon

        _pyngl.mpLimitMode = "Corners"
        _pyngl.mpLeftCornerLonF = geobounds.bottom_left.lon
        _pyngl.mpLeftCornerLatF = geobounds.bottom_left.lat
        _pyngl.mpRightCornerLonF = geobounds.top_right.lon
        _pyngl.mpRightCornerLatF = geobounds.top_right.lat

        for key, val in viewitems(kwargs):
            setattr(_pyngl, key, val)

        return _pyngl

    def _basemap(self, geobounds, **kwargs):
        if not basemap_enabled():
            return None

        local_kwargs = dict(projection="rotpole",
                            o_lat_p=self._bm_cart_pole_lat,
                            o_lon_p=self.pole_lon,
                            llcrnrlat=geobounds.bottom_left.lat,
                            urcrnrlat=geobounds.top_right.lat,
                            llcrnrlon=geobounds.bottom_left.lon,
                            urcrnrlon=geobounds.top_right.lon,
                            lon_0=self._bm_lon_0,
                            rsphere=Constants.WRF_EARTH_RADIUS,
                            resolution="l")

        local_kwargs.update(kwargs)
        _basemap = Basemap(**local_kwargs)

        return _basemap

    def _cartopy(self):
        if not cartopy_enabled():
            return None

        _cartopy = crs.RotatedPole(pole_longitude=self._cart_pole_lon,
                                   pole_latitude=self._bm_cart_pole_lat,
                                   central_rotated_longitude=(
                                       180.0 - self.pole_lon),  # Probably
                                   globe=self._globe())
        return _cartopy

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


def getproj(**proj_params):
    """Return a :class:`wrf.WrfProj` subclass.

    This functions serves as a factory function for returning a
    :class:`wrf.WrfProj` subclass from the specified map projection parameters.

    Args:

        **proj_params:  Map projection optional keyword arguments, that
            have the same names as found in WRF output NetCDF global
            attributes:

            - 'MAP_PROJ': The map projection type as an integer.
            - 'TRUELAT1': True latitude 1.
            - 'TRUELAT2': True latitude 2.
            - 'MOAD_CEN_LAT': Mother of all domains center latitude.
            - 'STAND_LON': Standard longitude.
            - 'POLE_LAT': Pole latitude.
            - 'POLE_LON': Pole longitude.

    Returns:

        :class:`wrf.WrfProj`: A :class:`wrf.WrfProj` subclass for the
        specified map projection parameters.

    """
    up_proj_params = dict_keys_to_upper(proj_params)

    proj_type = up_proj_params.get("MAP_PROJ", 0)
    if proj_type == ProjectionTypes.LAMBERT_CONFORMAL:
        return LambertConformal(**proj_params)
    elif proj_type == ProjectionTypes.POLAR_STEREOGRAPHIC:
        return PolarStereographic(**proj_params)
    elif proj_type == ProjectionTypes.MERCATOR:
        return Mercator(**proj_params)
    elif (proj_type == ProjectionTypes.ZERO or
          proj_type == ProjectionTypes.LAT_LON):
        if (up_proj_params.get("POLE_LAT", None) == 90.
                and up_proj_params.get("POLE_LON", None) == 0.):
            return LatLon(**proj_params)
        else:
            return RotatedLatLon(**proj_params)
    else:
        # Unknown projection
        return WrfProj(**proj_params)
