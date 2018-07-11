import numpy as np
from namelist import *
from boundaries import exchange_BC_all, exchange_BC
from IO import load_topo

def initialize_fields(GR):
    # CREATE ARRAYS
    # height
    HGHT = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    HTOP = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    # wind velocities
    UWIND = np.full( (GR.nxs+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VWIND = np.full( (GR.nx+2*GR.nb,GR.nys+2*GR.nb), np.nan)
    WIND = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    # mass fluxes at velocity points
    UFLX = np.full( (GR.nxs+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VFLX = np.full( (GR.nx+2*GR.nb,GR.nys+2*GR.nb), np.nan)

    # FOR MASSPOINT_FLUX_TENDENCY_UPSTREAM:
    # mass fluxes at mass points
    UFLXMP = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VFLXMP = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    # momentum fluxes at velocity points
    UUFLX = np.full( (GR.nxs+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VUFLX = np.full( (GR.nx+2*GR.nb,GR.nys+2*GR.nb), np.nan)
    UVFLX = np.full( (GR.nxs+2*GR.nb,GR.ny+2*GR.nb), np.nan)
    VVFLX = np.full( (GR.nx+2*GR.nb,GR.nys+2*GR.nb), np.nan)

    # FOR WIND TENDENCY JACOBSON
    # surface height
    HSURF = load_topo(GR) 
    #HSURF[GR.iijj] = 0.
    #HSURF = exchange_BC(GR, HSURF)
    # tracer
    TRACER = np.full( (GR.nx+2*GR.nb,GR.ny+2*GR.nb), np.nan)

    
    # INITIAL CONDITIONS
    HGHT[GR.iijj] = h0 - HSURF[GR.iijj]
    HGHT = gaussian2D(GR, HGHT, hpert, np.pi*3/4, 0, np.pi/10, np.pi/10)
    HGHT = random2D(GR, HGHT, h_random_pert)
    HTOP[GR.iijj] = HSURF[GR.iijj] + HGHT[GR.iijj]

    UWIND[GR.iisjj] = u0   
    UWIND = gaussian2D(GR, UWIND, upert, np.pi*3/4, 0, np.pi/10, np.pi/10)
    VWIND[GR.iijjs] = v0
    VWIND = gaussian2D(GR, VWIND, vpert, np.pi*3/4, 0, np.pi/10, np.pi/10)

    TRACER[GR.iijj] = 0.
    TRACER = gaussian2D(GR, TRACER, 10, np.pi*3/4, 0, np.pi/10, np.pi/10)

    # BOUNDARY CONDITIONS
    HGHT, UWIND, VWIND, TRACER = exchange_BC_all(GR, HGHT, UWIND, VWIND, TRACER)

    return(HGHT, HTOP, UWIND, VWIND, WIND,
            UFLX, VFLX, UFLXMP, VFLXMP,
            UUFLX, VUFLX, UVFLX, VVFLX,
            HSURF, TRACER)



def random2D(GR, FIELD, pert):
    FIELD = FIELD + pert*np.random.rand(FIELD.shape[0], FIELD.shape[1])
    return(FIELD)

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
