import unittest as ut
import numpy.testing as nt
import numpy as np
import numpy.ma as ma
import os
import sys
import subprocess

from wrf import xy_to_ll_proj, ll_to_xy_proj, to_np


class WRFLatLonProjTest(ut.TestCase):
    longMessage = True


def make_test(xy_or_ll_out):
    def test(self):

        # Python 3
        if sys.version_info > (3, ):
            assert_raises_regex = self.assertRaisesRegex
            xrange = range
        else:
            assert_raises_regex = self.assertRaisesRegexp

        if xy_or_ll_out == "xy":
            # Test the required failures
            with assert_raises_regex(ValueError, ".*map_proj.*"):
                ll_to_xy_proj(30, -110)

            with assert_raises_regex(ValueError, ".*ref_lat.*"):
                ll_to_xy_proj(30, -110, map_proj=1)

            with assert_raises_regex(ValueError, ".*ref_lon.*"):
                ll_to_xy_proj(30, -110, map_proj=1, ref_lat=45)

            with assert_raises_regex(ValueError, ".*known_x.*"):
                ll_to_xy_proj(30, -110, map_proj=1, ref_lat=45.0,
                              ref_lon=-120.)

            with assert_raises_regex(ValueError, ".*known_y.*"):
                ll_to_xy_proj(30, -110, map_proj=1, ref_lat=45.0,
                              ref_lon=-120., known_x=1)

            with assert_raises_regex(ValueError, ".*dx.*"):
                ll_to_xy_proj(30, -110, map_proj=1, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0)

            # Now test the projections

            # Lambert Conformal - truelat1, stand_lon required
            with assert_raises_regex(ValueError, ".*truelat1.*"):
                ll_to_xy_proj(30, -110, map_proj=1, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.)

            with assert_raises_regex(ValueError, ".*stand_lon.*"):
                ll_to_xy_proj(30, -110, map_proj=1, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.,
                              truelat1=60.)

            # Make sure it runs with all params set vs not
            p_all = ll_to_xy_proj(28., -89., map_proj=1, ref_lat=17.803,
                                  ref_lon=-100.7747, known_x=0, known_y=0,
                                  dx=30000., truelat1=30., truelat2=30.,
                                  stand_lon=-89., pole_lat=90., pole_lon=0.)

            p_def = ll_to_xy_proj(28., -89., map_proj=1, ref_lat=17.803,
                                  ref_lon=-100.7747, known_x=0, known_y=0,
                                  dx=30000., truelat1=30., stand_lon=-89.)

            nt.assert_allclose(to_np(p_all), to_np(p_def))

            # Polar Stereographic - truelat1, stand_lon

            with assert_raises_regex(ValueError, ".*truelat1.*"):
                ll_to_xy_proj(30, -110, map_proj=2, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.)

            with assert_raises_regex(ValueError, ".*stand_lon.*"):
                ll_to_xy_proj(30, -110, map_proj=2, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.,
                              truelat1=60.)

            p_all = ll_to_xy_proj(28., -89., map_proj=2, ref_lat=17.933,
                                  ref_lon=-100.0735, known_x=0, known_y=0,
                                  dx=30000., truelat1=30., truelat2=30.,
                                  stand_lon=-89., pole_lat=90., pole_lon=0.)

            p_def = ll_to_xy_proj(28., -89., map_proj=2, ref_lat=17.933,
                                  ref_lon=-100.0735, known_x=0, known_y=0,
                                  dx=30000., truelat1=30., stand_lon=-89.)

            nt.assert_allclose(to_np(p_all), to_np(p_def))

            # Mercator - truelat1

            with assert_raises_regex(ValueError, ".*truelat1.*"):
                ll_to_xy_proj(30, -110, map_proj=2, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.)

            p_all = ll_to_xy_proj(28., -89., map_proj=3, ref_lat=19.1075,
                                  ref_lon=-101.008, known_x=0, known_y=0,
                                  dx=30000., truelat1=30., truelat2=30.,
                                  stand_lon=-89., pole_lat=90., pole_lon=0.)

            p_def = ll_to_xy_proj(28., -89., map_proj=3, ref_lat=19.1075,
                                  ref_lon=-101.008, known_x=0, known_y=0,
                                  dx=30000., truelat1=30.)

            nt.assert_allclose(to_np(p_all), to_np(p_def))

            # Lat/lon - stand_lon, dy, pole_lat, pole_lon

            with assert_raises_regex(ValueError, ".*stand_lon.*"):
                ll_to_xy_proj(30, -110, map_proj=6, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=.2698388)

            with assert_raises_regex(ValueError, ".*dy.*"):
                ll_to_xy_proj(30, -110, map_proj=6, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              stand_lon=89.0,
                              dx=.2698388)

            with assert_raises_regex(ValueError, ".*pole_lat.*"):
                ll_to_xy_proj(30, -110, map_proj=6, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              stand_lon=89.0,
                              dx=.2698388, dy=.2698388)

            with assert_raises_regex(ValueError, ".*pole_lon.*"):
                ll_to_xy_proj(30, -110, map_proj=6, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              stand_lon=89.0,
                              dx=.2698388, dy=.2698388,
                              pole_lat=62.0)

            p_all = ll_to_xy_proj(28., -89., map_proj=6, ref_lat=17.6759,
                                  ref_lon=-101.4286, known_x=0, known_y=0,
                                  dx=30000, dy=30000,
                                  truelat1=30., truelat2=30.,
                                  stand_lon=-89., pole_lat=62.0,
                                  pole_lon=180.)

            p_def = ll_to_xy_proj(28., -89., map_proj=6, ref_lat=17.6759,
                                  ref_lon=-101.4286, known_x=0, known_y=0,
                                  stand_lon=-89.,
                                  dx=30000, dy=30000, pole_lat=62.0,
                                  pole_lon=180.)

            nt.assert_allclose(to_np(p_all), to_np(p_def))

        if xy_or_ll_out == "ll":

            # Test the required failures
            with assert_raises_regex(ValueError, ".*map_proj.*"):
                xy_to_ll_proj(45, 50)

            with assert_raises_regex(ValueError, ".*ref_lat.*"):
                xy_to_ll_proj(45, 50, map_proj=1)

            with assert_raises_regex(ValueError, ".*ref_lon.*"):
                xy_to_ll_proj(45, 50, map_proj=1, ref_lat=45)

            with assert_raises_regex(ValueError, ".*known_x.*"):
                xy_to_ll_proj(45, 50, map_proj=1, ref_lat=45.0,
                              ref_lon=-120.)

            with assert_raises_regex(ValueError, ".*known_y.*"):
                xy_to_ll_proj(45, 50, map_proj=1, ref_lat=45.0,
                              ref_lon=-120., known_x=1)

            with assert_raises_regex(ValueError, ".*dx.*"):
                xy_to_ll_proj(45, 50, map_proj=1, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0)

            # Now test the projections

            # Lambert Conformal - truelat1, stand_lon required
            with assert_raises_regex(ValueError, ".*truelat1.*"):
                xy_to_ll_proj(45, 50, map_proj=1, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.)

            with assert_raises_regex(ValueError, ".*stand_lon.*"):
                xy_to_ll_proj(45, 50, map_proj=1, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.,
                              truelat1=60.)

            # Make sure it runs with all params set vs not
            p_all = xy_to_ll_proj(45, 50, map_proj=1, ref_lat=17.803,
                                  ref_lon=-100.7747, known_x=0, known_y=0,
                                  dx=30000., truelat1=30., truelat2=30.,
                                  stand_lon=-89., pole_lat=90., pole_lon=0.)

            p_def = xy_to_ll_proj(45, 50, map_proj=1, ref_lat=17.803,
                                  ref_lon=-100.7747, known_x=0, known_y=0,
                                  dx=30000., truelat1=30.,
                                  stand_lon=-89.)

            nt.assert_allclose(to_np(p_all), to_np(p_def))

            # Polar Stereographic - truelat1, stand_lon

            with assert_raises_regex(ValueError, ".*truelat1.*"):
                xy_to_ll_proj(45, 50, map_proj=2, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.)

            with assert_raises_regex(ValueError, ".*stand_lon.*"):
                xy_to_ll_proj(45, 50, map_proj=2, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.,
                              truelat1=60.)

            p_all = xy_to_ll_proj(45, 50, map_proj=2, ref_lat=17.933,
                                  ref_lon=-100.0735, known_x=0, known_y=0,
                                  dx=30000., truelat1=30., truelat2=30.,
                                  stand_lon=-89., pole_lat=90., pole_lon=0.)

            p_def = xy_to_ll_proj(45, 50, map_proj=2, ref_lat=17.933,
                                  ref_lon=-100.0735, known_x=0, known_y=0,
                                  dx=30000., truelat1=30.,
                                  stand_lon=-89.)

            nt.assert_allclose(to_np(p_all), to_np(p_def))

            # Mercator - truelat1

            with assert_raises_regex(ValueError, ".*truelat1.*"):
                xy_to_ll_proj(45, 50, map_proj=2, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=3000.)

            p_all = xy_to_ll_proj(45, 50, map_proj=3, ref_lat=19.1075,
                                  ref_lon=-101.008, known_x=0, known_y=0,
                                  dx=30000., truelat1=30., truelat2=30.,
                                  stand_lon=-89., pole_lat=90., pole_lon=0.)

            p_def = xy_to_ll_proj(45, 50, map_proj=3, ref_lat=19.1075,
                                  ref_lon=-101.008, known_x=0, known_y=0,
                                  dx=30000., truelat1=30.)

            nt.assert_allclose(to_np(p_all), to_np(p_def))

            # Lat/lon - stand_lon, dy, pole_lat, pole_lon

            with assert_raises_regex(ValueError, ".*stand_lon.*"):
                xy_to_ll_proj(45, 50, map_proj=6, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              dx=.2698388)

            with assert_raises_regex(ValueError, ".*dy.*"):
                xy_to_ll_proj(45, 50, map_proj=6, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              stand_lon=89.0,
                              dx=.2698388)

            with assert_raises_regex(ValueError, ".*pole_lat.*"):
                xy_to_ll_proj(45, 50, map_proj=6, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              stand_lon=89.0,
                              dx=.2698388, dy=.2698388)

            with assert_raises_regex(ValueError, ".*pole_lon.*"):
                xy_to_ll_proj(45, 50, map_proj=6, ref_lat=45.0,
                              ref_lon=-120., known_x=0, known_y=0,
                              stand_lon=89.0,
                              dx=.2698388, dy=.2698388,
                              pole_lat=62.0)

            p_all = xy_to_ll_proj(64, 40, map_proj=6, ref_lat=17.6759,
                                  ref_lon=-101.4286, known_x=0, known_y=0,
                                  dx=30000, dy=30000,
                                  truelat1=30., truelat2=30.,
                                  stand_lon=-89., pole_lat=62.0,
                                  pole_lon=180.)

            p_def = xy_to_ll_proj(64, 40, map_proj=6, ref_lat=17.6759,
                                  ref_lon=-101.4286, known_x=0, known_y=0,
                                  stand_lon=-89.,
                                  dx=30000, dy=30000, pole_lat=62.0,
                                  pole_lon=180.)

            nt.assert_allclose(to_np(p_all), to_np(p_def))

    return test


if __name__ == "__main__":

    for v in ("xy", "ll"):
        test_func = make_test(v)
        setattr(WRFLatLonProjTest, 'test_{0}'.format(v), test_func)

    ut.main()
