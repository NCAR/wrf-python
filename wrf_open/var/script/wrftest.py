
import wrf.var as w
import numpy as n

from netCDF4 import Dataset as NetCDF 

def main():
    wrfnc = NetCDF("/Users/bladwig/wrfout_d03_2003-05-07_09:00:00")
    
    # Cape NO RESULTS FOR LCL OR LFC
    cape, cin, lcl, lfc = w.getvar(wrfnc, "cape2d")
    #cape, cin = w.getvar(wrfnc, "cape3d")
    print n.amax(cape)
    print n.amax(cin)
    print n.amax(lcl)
    print n.amax(lfc)
    
    
    # DBZ 
    dbz = w.getvar(wrfnc, "dbz")
    print n.amax(dbz)
    
    # DP
    dp = w.getvar(wrfnc, "dp", units="f")
    print n.amax(dp)
    
    dp2 = w.getvar(wrfnc, "dp2m", units="f")
    print n.amax(dp2)
    
    # Height
    ht = w.getvar(wrfnc, "height", msl=False, units="m")
    print n.amax(ht)
    
    geopt = w.getvar(wrfnc, "geopt")
    print n.amax(geopt)
    
    # Helicity
    srh = w.getvar(wrfnc, "srh")
    print n.amax(srh)
    
    uhel = w.getvar(wrfnc, "uhel")
    print n.amax(uhel)
    
    # Omega (Not sure if this is correct, and units aren't C)
    omega = w.getvar(wrfnc, "omega")
    print n.amax(omega)
    
    # Precip Water (NOT SURE)
    pw = w.getvar(wrfnc, "pw")
    print n.amax(pw)
    
    # RH
    rh = w.getvar(wrfnc, "rh")
    print n.amax(rh)
    
    rh2 = w.getvar(wrfnc, "rh2m")
    print n.amax(rh2)
    
    # SLP
    slp = w.getvar(wrfnc, "slp", units="hpa")
    print n.amax(slp)
    
    # TEMP
    t = w.getvar(wrfnc, "temp", units="f")
    print n.amax(t)
    
    # ETH VALUES SEEM HIGH....
    eth = w.getvar(wrfnc, "theta_e", units="k")
    print n.amax(eth)
    
    tv = w.getvar(wrfnc, "tv", units="k")
    print n.amax(tv)
    
    # Note: NCL says this is in 'C', but appears to be 'K'
    tw = w.getvar(wrfnc, "tw", units="f")
    print n.amax(tw)
    
    # WIND
    umet,vmet = w.getvar(wrfnc, "uvmet", units="kts")
    print n.amax(umet)
    print n.amax(vmet)
    
    umet10,vmet10 = w.getvar(wrfnc, "uvmet10", units="kts")
    print n.amax(umet10)
    print n.amax(vmet10)
    
    
    
    # TERRAIN
    ter = w.getvar(wrfnc, "terrain", units="dm")
    print n.amax(ter)
    
    # VORTICITY
    avo = w.getvar(wrfnc, "avo")
    print n.amax(avo)
    
    pvo = w.getvar(wrfnc, "pvo")
    print n.amax(pvo)
    
    # LAT/LON
    lat = w.getvar(wrfnc, "lat")
    print n.amax(lat)
    print n.amin(lat)
    
    lon = w.getvar(wrfnc, "lon")
    print n.amax(lon)
    print n.amin(lon)
    
    i,j = w.get_ij(wrfnc, -97.516540, 35.467787)
    print i,j
    
    lon, lat = w.get_ll(wrfnc, 33.5, 33.5)
    print lon, lat
    
    #ETA -- Result somewhat different than geopt
    z = w.convert_eta(wrfnc, msl=False, units="m")
    print n.amax(z)
    
    diff = n.abs(z - ht)/ht * 100.0
    print n.amin(diff), n.amax(diff)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

if __name__ == "__main__":
    main()
