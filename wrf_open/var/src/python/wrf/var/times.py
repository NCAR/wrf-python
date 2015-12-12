
import datetime as dt

__all__ = ["get_times"]

def _make_time(timearr):
    return dt.datetime.strptime("".join(timearr[:]), "%Y-%m-%d_%H:%M:%S")

# TODO:  handle list of files and multiple times
def get_times(wrfnc):
    times = wrfnc.variables["Times"][:,:]
    return [_make_time(times[i,:]) for i in xrange(times.shape[0])]
