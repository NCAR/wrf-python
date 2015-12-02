from wrf.var.extension import computetd
from wrf.var.decorators import convert_units

__all__ = ["get_dp", "get_dp_2m"]

@convert_units("temp", "c")
def get_dp(wrfnc, units="c", timeidx=0):
    
    p = wrfnc.variables["P"][timeidx,:,:,:]
    pb = wrfnc.variables["PB"][timeidx,:,:,:]
    qvapor = wrfnc.variables["QVAPOR"][timeidx,:,:,:]
    
    # Algorithm requires hPa
    full_p = .01*(p + pb)
    qvapor[qvapor < 0] = 0
    
    td = computetd(full_p, qvapor)
    return td
    
@convert_units("temp", "c")
def get_dp_2m(wrfnc, units="c", timeidx=0):

    # Algorithm requires hPa
    psfc = .01*(wrfnc.variables["PSFC"][timeidx,:,:])
    q2 = wrfnc.variables["Q2"][timeidx,:,:]
    q2[q2 < 0] = 0
    
    td = computetd(psfc, q2)
    
    return td

