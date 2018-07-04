import numpy as np
from namelist import *
from boundaries import exchange_BC_all

def initialize_fields(GR):
    # CREATE ARRAYS
    # height
    HGHT = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    # wind velocities
    UWIND = np.full( (GR.nxs+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VWIND = np.full( (GR.nx+2*GR.nb,GR.nys+2*GR.nb), np.nan)
    WIND = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    # mass fluxes at velocity points
    UFLX = np.full( (GR.nxs+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VFLX = np.full( (GR.nx+2*GR.nb,GR.nys+2*GR.nb), np.nan)
    # mass fluxes at mass points
    UFLX = np.full( (GR.nxs+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    UFLX_mp = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VFLX_mp = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    # momentum fluxes at velocity points
    UUFLX = np.full( (GR.nxs+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VVFLX = np.full( (GR.nx+2*GR.nb,GR.nys+2*GR.nb), np.nan)

    # INITIAL CONDITIONS
    HGHT[GR.iijj] = h0
    UWIND[GR.iisjj] = u0   
    VWIND[GR.iijjs] = 0.

    HGHT = gaussian2D(GR, HGHT, hpert, np.pi*3/4, 0, np.pi/10, np.pi/10)

    # BOUNDARY CONDITIONS
    HGHT, UWIND, VWIND = exchange_BC_all(GR, HGHT, UWIND, VWIND)

    return(HGHT, UWIND, VWIND, WIND, UFLX, VFLX, UFLX_mp, VFLX_mp, UUFLX, VVFLX)




def gaussian2D(GR, FIELD, pert, lon0_rad, lat0_rad, lonSig_rad, latSig_rad):

    dimx,dimy = FIELD.shape

    if (dimy == GR.nys+2*GR.nb): # staggered in y 
        selinds = GR.iijjs
    elif (dimx == GR.nxs+2*GR.nb): # staggered in x 
        selinds = GR.iisjj
    else: # unstaggered in y and x 
        selinds = GR.iijj

    perturb = pert*np.exp( \
            - np.power(GR.lonjs_rad[selinds] - lon0_rad, 2)/(2*lonSig_rad**2) \
            - np.power(GR.latjs_rad[selinds] - lat0_rad, 2)/(2*latSig_rad**2) )
    FIELD[selinds] = FIELD[selinds] + perturb

    return(FIELD)
