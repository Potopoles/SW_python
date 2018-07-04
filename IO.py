import numpy as np
from netCDF4 import Dataset


def output_to_NC(GR, outCounter, HGHT, UWIND, VWIND, WIND):
    filename = 'output/out'+str(outCounter).zfill(4)+'.nc'

    ncf = Dataset(filename, 'w', format='NETCDF4')
    ncf.close()

    ncf = Dataset(filename, 'a', format='NETCDF4')

    # DIMENSIONS
    time_dim = ncf.createDimension('time', None)
    lon_dim = ncf.createDimension('lon', GR.nx)
    lons_dim = ncf.createDimension('lons', GR.nxs)
    lat_dim = ncf.createDimension('lat', GR.ny)
    lats_dim = ncf.createDimension('lats', GR.nys)

    # DIMENSION VARIABLES
    time = ncf.createVariable('time', 'f8', ('time',) )
    lon = ncf.createVariable('lon', 'f4', ('lon',) )
    lons = ncf.createVariable('lons', 'f4', ('lons',) )
    lat = ncf.createVariable('lat', 'f4', ('lat',) )
    lats = ncf.createVariable('lats', 'f4', ('lats',) )


    time[:] = 0
    lon[:] = GR.lon_rad[GR.ii,GR.nb+1]
    lons[:] = GR.lonis_rad[GR.iis,GR.nb+1]
    lat[:] = GR.lat_rad[GR.nb+1,GR.jj]
    lats[:] = GR.latjs_rad[GR.nb+1,GR.jjs]

    # VARIABLES
    HGHT_out = ncf.createVariable('HGHT', 'f4', ('time', 'lat', 'lon',) )
    UWIND_out = ncf.createVariable('UWIND', 'f4', ('time', 'lat', 'lons',) )
    VWIND_out = ncf.createVariable('VWIND', 'f4', ('time', 'lats', 'lon',) )
    WIND_out = ncf.createVariable('WIND', 'f4', ('time', 'lat', 'lon',) )

    HGHT_out[-1,:,:] = HGHT[GR.iijj].T
    UWIND_out[-1,:,:] = UWIND[GR.iisjj].T
    VWIND_out[-1,:,:] = VWIND[GR.iijjs].T
    WIND_out[-1,:,:] = WIND[GR.iijj].T

    ncf.close()
