import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import interp2d
from boundaries import exchange_BC_rigid_y, exchange_BC_periodic_x

def load_topo(GR):
    HSURF = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    filename = '../elevation/elev.1-deg.nc'
    ncf = Dataset(filename, 'r', format='NETCDF4')
    lon_inp = ncf['lon'][:]
    lat_inp = ncf['lat'][:]
    hsurf_inp = ncf['data'][0,:,:]
    interp = interp2d(lon_inp, lat_inp, hsurf_inp)
    HSURF[GR.iijj] = interp(GR.lon_deg[GR.ii,GR.nb+1], GR.lat_deg[GR.nb+1,GR.jj]).T
    HSURF[HSURF < 0] = 0
    HSURF = exchange_BC_periodic_x(GR, HSURF)
    HSURF = exchange_BC_rigid_y(GR, HSURF)

    n_smooth = 20
    tau = 0.2
    for i in range(0,n_smooth):
        HSURF[GR.iijj] = HSURF[GR.iijj] + tau*(HSURF[GR.iijj_im1] + HSURF[GR.iijj_ip1] + \
                                                HSURF[GR.iijj_jm1] + HSURF[GR.iijj_jp1] - \
                                                4*HSURF[GR.iijj]) 
        HSURF = exchange_BC_periodic_x(GR, HSURF)
        HSURF = exchange_BC_rigid_y(GR, HSURF)

    return(HSURF)


def output_to_NC(GR, outCounter, HGHT, HTOP, UWIND, VWIND, WIND,
                HSURF, TRACER,
                mean_ekin):
    filename = '../output/out'+str(outCounter).zfill(4)+'.nc'

    ncf = Dataset(filename, 'w', format='NETCDF4')
    ncf.close()

    ncf = Dataset(filename, 'a', format='NETCDF4')

    # DIMENSIONS
    time_dim = ncf.createDimension('time', None)
    bnds_dim = ncf.createDimension('bnds', 2)
    lon_dim = ncf.createDimension('lon', GR.nx)
    lons_dim = ncf.createDimension('lons', GR.nxs)
    lat_dim = ncf.createDimension('lat', GR.ny)
    lats_dim = ncf.createDimension('lats', GR.nys)

    # DIMENSION VARIABLES
    time = ncf.createVariable('time', 'f8', ('time',) )
    bnds = ncf.createVariable('bnds', 'f8', ('bnds',) )
    lon = ncf.createVariable('lon', 'f4', ('lon',) )
    lons = ncf.createVariable('lons', 'f4', ('lons',) )
    lat = ncf.createVariable('lat', 'f4', ('lat',) )
    lats = ncf.createVariable('lats', 'f4', ('lats',) )

    time[:] = outCounter
    bnds[:] = [0,1]
    lon[:] = GR.lon_rad[GR.ii,GR.nb+1]
    lons[:] = GR.lonis_rad[GR.iis,GR.nb+1]
    lat[:] = GR.lat_rad[GR.nb+1,GR.jj]
    lats[:] = GR.latjs_rad[GR.nb+1,GR.jjs]

    # VARIABLES
    HGHT_out = ncf.createVariable('HGHT', 'f4', ('time', 'lat', 'lon',) )
    HTOP_out = ncf.createVariable('HTOP', 'f4', ('time', 'lat', 'lon',) )
    UWIND_out = ncf.createVariable('UWIND', 'f4', ('time', 'lat', 'lons',) )
    VWIND_out = ncf.createVariable('VWIND', 'f4', ('time', 'lats', 'lon',) )
    WIND_out = ncf.createVariable('WIND', 'f4', ('time', 'lat', 'lon',) )
    HSURF_out = ncf.createVariable('HSURF', 'f4', ('bnds', 'lat', 'lon',) )
    TRACER_out = ncf.createVariable('TRACER', 'f4', ('time', 'lat', 'lon',) )
    mean_ekin_out = ncf.createVariable('mean_ekin', 'f4', ('time', 'bnds',) )

    HGHT_out[-1,:,:] = HGHT[GR.iijj].T
    HTOP_out[-1,:,:] = HTOP[GR.iijj].T
    UWIND_out[-1,:,:] = UWIND[GR.iisjj].T
    VWIND_out[-1,:,:] = VWIND[GR.iijjs].T
    WIND_out[-1,:,:] = WIND[GR.iijj].T
    HSURF_out[0,:,:] = HSURF[GR.iijj].T
    HSURF_out[1,:,:] = HSURF[GR.iijj].T
    TRACER_out[-1,:,:] = TRACER[GR.iijj].T
    mean_ekin_out[-1,:] = mean_ekin

    ncf.close()
