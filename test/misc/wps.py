# Hastily made script to read WPS intermediate files
from __future__ import print_function
import struct

import numpy as np

# Number of bytes used at the start and end of a fortran record to
# indicate the record size
SIZE_BYTES = 4


class WPSData(object):
    def __init__(self, ifv=None, hdate=None, xfcst=None, map_source=None,
                 field=None, units=None, desc=None, xlvl=None, nx=None,
                 ny=None, iproj=None, startloc=None, startlat=None,
                 startlon=None, deltalat=None, deltalon=None, nlats=None,
                 dx=None, dy=None, xlonc=None, truelat1=None, truelat2=None,
                 earth_radius=None, is_wind_earth_rel=None, slab=None):

        self.ifv = ifv
        self.hdate = hdate
        self.xfcst = xfcst
        self.map_source = map_source
        self.field = field
        self.units = units
        self.desc = desc
        self.xlvl = xlvl
        self.nx = nx
        self.ny = ny
        self.iproj = iproj
        self.startloc = startloc
        self.startlat = startlat
        self.startlon = startlon
        self.deltalat = deltalat
        self.deltalon = deltalon
        self.nlats = nlats
        self.dx = dx
        self.dy = dy
        self.xlonc = xlonc
        self.truelat1 = truelat1
        self.truelat2 = truelat2
        self.earth_radius = earth_radius
        self.is_wind_earth_rel = is_wind_earth_rel
        self.slab = slab


def _parse_record1(data):
    result = {}
    result["ifv"] = struct.unpack(">i", data)

    return result


def _parse_record2(data):
    result = {}
    parsed = struct.unpack(">24sf32s9s25s46sfiii", data)
    result["hdate"] = parsed[0]
    result["xfcst"] = parsed[1]
    result["map_source"] = parsed[2]
    result["field"] = parsed[3]
    result["units"] = parsed[4]
    result["desc"] = parsed[5]
    result["xlvl"] = parsed[6]
    result["nx"] = parsed[7]
    result["ny"] = parsed[8]
    result["iproj"] = parsed[9]

    return result


def _parse_record3(data, iproj):
    result = {}
    if iproj == 0:
        fmt = ">8sfffff"
        parsed = struct.unpack(fmt, data)
        result["startloc"] = parsed[0]
        result["startlat"] = parsed[1]
        result["startlon"] = parsed[2]
        result["deltalat"] = parsed[3]
        result["deltalon"] = parsed[4]
        result["earth_radius"] = parsed[5]
    elif iproj == 1:
        fmt = ">8sffffff"
        parsed = struct.unpack(fmt, data)
        result["startloc"] = parsed[0]
        result["startlat"] = parsed[1]
        result["startlon"] = parsed[2]
        result["dx"] = parsed[3]
        result["dy"] = parsed[4]
        result["truelat1"] = parsed[5]
        result["earth_radius"] = parsed[6]
    elif iproj == 3:
        fmt = ">8sffffffff"
        parsed = struct.unpack(fmt, data)
        result["startloc"] = parsed[0]
        result["startlat"] = parsed[1]
        result["startlon"] = parsed[2]
        result["dx"] = parsed[3]
        result["dy"] = parsed[4]
        result["xlonc"] = parsed[5]
        result["truelat1"] = parsed[6]
        result["truelat2"] = parsed[7]
        result["earth_radius"] = parsed[8]
    elif iproj == 4:
        fmt = ">8sfffff"
        parsed = struct.unpack(fmt, data)
        result["startloc"] = parsed[0]
        result["startlat"] = parsed[1]
        result["startlon"] = parsed[2]
        result["nlats"] = parsed[3]
        result["deltalon"] = parsed[4]
        result["earth_radius"] = parsed[5]
    elif iproj == 5:
        fmt = ">8sfffffff"
        parsed = struct.unpack(fmt, data)
        result["startloc"] = parsed[0]
        result["startlat"] = parsed[1]
        result["startlon"] = parsed[2]
        result["dx"] = parsed[3]
        result["dy"] = parsed[4]
        result["xlonc"] = parsed[5]
        result["truelat1"] = parsed[6]
        result["earth_radius"] = parsed[7]

    return result


def _parse_record4(data):
    result = {}
    result["is_wind_earth_rel"] = struct.unpack(">i", data)

    return result


def _parse_record5(data, nx, ny):
    result = {}

    size = nx * ny
    fmt = ">{}f".format(size)
    parsed = struct.unpack(fmt, data)

    arr = np.array(parsed, dtype=np.float32)
    result["slab"] = arr.reshape((nx, ny), order="F")

    return result


_PARSE_MAP = {0: _parse_record1,
              1: _parse_record2,
              2: _parse_record3,
              3: _parse_record4,
              4: _parse_record5}


def parse_record(record_idx, data, iproj=None, nx=None, ny=None):

    if record_idx == 0 or record_idx == 1 or record_idx == 3:
        return _PARSE_MAP[record_idx](data)
    elif record_idx == 2:
        return _PARSE_MAP[record_idx](data, iproj)
    elif record_idx == 4:
        return _PARSE_MAP[record_idx](data, nx, ny)


def read_wps(wps_file, field=None):
    with open(wps_file, "rb") as f:
        wps_params = {}
        record_list = []
        raw_data = f.read()

        record_size_idx = 0
        end_of_file_idx = len(raw_data) - 1

        while True:
            iproj = None
            nx = None
            ny = None
            keep_record = True
            for record_idx in range(5):
                # Each record has the size (in SIZE_BYTES bytes) at the
                # start and end of each record.  This might be compiler
                # dependent though, so this might need to be modified.
                # Also, the WPS files are stored big endian.

                record_size = struct.unpack(
                    ">i",
                    raw_data[record_size_idx: record_size_idx + SIZE_BYTES])
                record_start = record_size_idx + SIZE_BYTES
                record_end = record_start + record_size[0]
                record_data = raw_data[record_start:record_end]

                parsed_record = parse_record(record_idx, record_data, iproj,
                                             nx, ny)

                try:
                    field_name = parsed_record["field"].strip()
                except KeyError:
                    pass
                else:
                    if field is not None:
                        if field_name != field:
                            keep_record = False

                try:
                    iproj = parsed_record["iproj"]
                except KeyError:
                    pass

                try:
                    nx = parsed_record["nx"]
                except KeyError:
                    pass

                try:
                    ny = parsed_record["ny"]
                except KeyError:
                    pass

                wps_params.update(parsed_record)

                record_size_idx = record_end + SIZE_BYTES

            if keep_record:
                record_list.append(WPSData(**wps_params))

            # Repeat for all record slabs
            if record_end + SIZE_BYTES > end_of_file_idx:
                break

    return record_list
