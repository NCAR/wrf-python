import unittest as ut
from os.path import join, basename

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from netCDF4 import Dataset as NetCDF


PYNGL = True
try:
    import Ngl
except ImportError:
    PYNGL = False
    
BASEMAP = True
try:
    import mpl_toolkits.basemap
except ImportError:
    BASEMAP = False

CARTOPY = True
try:
    from cartopy import crs, feature
except ImportError:
    CARTOPY = False
    

from wrf import get_proj_params
from wrf.projection import getproj, RotatedLatLon, PolarStereographic

FILE_DIR = "/Users/ladwig/Documents/wrf_files/"
WRF_FILES = [
            join(FILE_DIR, "rotated_pole", "EAS_geo_em.d01.nc"),
            join(FILE_DIR, "rotated_pole", "EUR_geo_em.d01.nc"),
            join(FILE_DIR,"wrfout_d01_2016-02-25_18_00_00"),
            join(FILE_DIR, "wrfout_d01_2008-09-29_23-30-00"),
            join(FILE_DIR, "wrfout_d01_2010-06-13_21:00:00")]


def nz_proj():
    lats = np.array([[-47.824014, -47.824014],
                     [-32.669853, -32.669853]])
    lons = np.array([[163.839595, -179.693502],
                     [163.839595, -179.693502]])
    
    params = {"MAP_PROJ" : 6,
              "CEN_LAT" : -41.814869,
              "CEN_LON" : 179.693502,
              "TRUELAT1" : 0,
              "TRUELAT2": 0,
              "MOAD_CEN_LAT" : -41.814869,
              "STAND_LON" : 180.0 - 179.693502,
              "POLE_LAT" : 48.185131,
              "POLE_LON" : 0.0}
    
    return lats, lons, RotLatLonProj(lats=lats, lons=lons, **params)

def argentina_proj():
    lats = np.array([[-57.144064, -57.144064],
                     [-21.154470, -21.154470]])
    lons = np.array([[-86.893797, -37.089724],
                     [-86.893797, -37.089724]])
    
    params = {"MAP_PROJ" : 6,
              "CEN_LAT" : -39.222954,
              "CEN_LON" : -65.980109,
              "TRUELAT1" : 0,
              "TRUELAT2": 0,
              "MOAD_CEN_LAT" : -39.222954,
              "STAND_LON" : 180.0 - -65.980109,
              "POLE_LAT" : 90 + -39.222954,
              "POLE_LON" : 0.0}
    
    return lats, lons, RotLatLonProj(lats=lats, lons=lons, **params)

def south_polar_proj():
    lats = np.array([[-30.0,-30.0],
                     [-30.0,-30.0]])
    lons = np.array([[-120, 60],
                     [-120, 60]])
    
    params = {"MAP_PROJ" : 2,
              "CEN_LAT" : -90.0,
              "CEN_LON" : 0,
              "TRUELAT1" : -10.0,
              "MOAD_CEN_LAT" : -90.0,
              "STAND_LON" : 0}
    
    return lats, lons, PolarStereographicProj(lats=lats, lons=lons, **params)

def north_polar_proj():
    lats = np.array([[30.0,30.0],
                     [30.0,30.0]])
    lons = np.array([[-45, 140],
                     [-45, 140]])
    
    params = {"MAP_PROJ" : 2,
              "CEN_LAT" : 90.0,
              "CEN_LON" : 10,
              "TRUELAT1" : 10.0,
              "MOAD_CEN_LAT" : 90.0,
              "STAND_LON" : 10}
    
    return lats, lons, PolarStereographicProj(lats=lats, lons=lons, **params)


def dateline_rot_proj():
    lats = np.array([[60.627974, 60.627974],
                     [71.717521, 71.717521]])
    lons = np.array([[170.332771, -153.456292],
                     [170.332771, -153.456292]])
    
    params = {"MAP_PROJ" : 6,
              "CEN_LAT" : 66.335764,
              "CEN_LON" : -173.143792,
              "TRUELAT1" : 0,
              "TRUELAT2": 0,
              "MOAD_CEN_LAT" : 66.335764,
              "STAND_LON" :  173.143792,
              "POLE_LAT" : 90.0 - 66.335764,
              "POLE_LON" : 180.0}
    return lats, lons, RotLatLonProj(lats=lats, lons=lons, **params)

class WRFProjTest(ut.TestCase):
    longMessage = True

def make_test(wrf_file=None, fixed_case=None):
    if wrf_file is not None:
        ncfile = NetCDF(wrf_file)
        lats, lons, proj_params = get_proj_params(ncfile)
        proj = getproj(lats=lats, lons=lons, **proj_params)
        name_suffix = basename(wrf_file)
    elif fixed_case is not None:
        name_suffix = fixed_case
        if fixed_case == "south_rot":
            lats, lons, proj = nz_proj()
        elif fixed_case == "arg_rot":
            lats, lons, proj = argentina_proj()
        elif fixed_case == "south_polar":
            lats, lons, proj = south_polar_proj()
        elif fixed_case == "north_polar":
            lats, lons, proj = north_polar_proj()
        elif fixed_case == "dateline_rot":
            lats,lons,proj = dateline_rot_proj()
    
    print ("wrf proj4: {}".format(proj.proj4()))
    if PYNGL:
        # PyNGL plotting
        wks_type = "png"
        wks = Ngl.open_wks(wks_type,"pyngl_{}".format(name_suffix))
        mpres = proj.pyngl()
        map = Ngl.map(wks,mpres)
        
        Ngl.delete_wks(wks)
    
    if BASEMAP: 
        # Basemap plotting
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        
        # Define and plot the meridians and parallels
        min_lat = np.amin(lats)
        max_lat = np.amax(lats)
        min_lon = np.amin(lons)
        max_lon = np.amax(lons)
        
        parallels = np.arange(np.floor(min_lat), np.ceil(max_lat), 
                              (max_lat - min_lat)/5.0)
        meridians = np.arange(np.floor(min_lon), np.ceil(max_lon), 
                              (max_lon - min_lon)/5.0)
    
        bm = proj.basemap()
        bm.drawcoastlines(linewidth=.5)
        #bm.drawparallels(parallels,labels=[1,1,1,1],fontsize=10)
        #bm.drawmeridians(meridians,labels=[1,1,1,1],fontsize=10)
        print ("basemap proj4: {}".format(bm.proj4string)) 
        plt.savefig("basemap_{}.png".format(name_suffix))
        plt.close(fig)
    
    if CARTOPY:
        # Cartopy plotting
        fig = plt.figure(figsize=(10,10))
        ax = plt.axes([0.1,0.1,0.8,0.8], projection=proj.cartopy())
        print ("cartopy proj4: {}".format(proj.cartopy().proj4_params))
        
        ax.coastlines('50m', linewidth=0.8)
        #print proj.x_extents()
        #print proj.y_extents()
        ax.set_xlim(proj.cartopy_xlim())
        ax.set_ylim(proj.cartopy_ylim())
        ax.gridlines()
        plt.savefig("cartopy_{}.png".format(name_suffix))
        plt.close(fig)

if __name__ == "__main__":
    for wrf_file in WRF_FILES:
        test_func = make_test(wrf_file=wrf_file)
        setattr(WRFProjTest, "test_proj", test_func)
    
    test_func2 = make_test(fixed_case="south_rot")
    setattr(WRFProjTest, "test_south_rot", test_func2)
     
    test_func3 = make_test(fixed_case="arg_rot")
    setattr(WRFProjTest, "test_arg_rot", test_func3)
    
    test_func4 = make_test(fixed_case="south_polar")
    setattr(WRFProjTest, "test_south_polar", test_func4)
     
    test_func5 = make_test(fixed_case="north_polar")
    setattr(WRFProjTest, "test_north_polar", test_func5)
     
    test_func6 = make_test(fixed_case="dateline_rot")
    setattr(WRFProjTest, "test_dateline_rot", test_func6)
    
    ut.main()