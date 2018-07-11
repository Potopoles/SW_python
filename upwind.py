import numpy as np
from height import height_tendency_upstream
from wind import masspoint_flux_tendency_upstream
from tracer import tracer_tendency_upstream
from boundaries import exchange_BC


def tendencies_upwind(GR, HGHT, TRACER, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX):

    # PROGNOSE HGHT
    dHGHTdt = height_tendency_upstream(GR, HGHT, UWIND, VWIND, UFLX, VFLX)

    # PROGNOSE WIND
    dUFLXMPdt, dVFLXMPdt = masspoint_flux_tendency_upstream(GR, UFLXMP, VFLXMP, HGHT,
                                                    UWIND, VWIND,
                                                    UUFLX, VUFLX, UVFLX, VVFLX,
                                                    HSURF)

    # PROGNOSE TRACER
    dTRACERdt = tracer_tendency_upstream(GR, TRACER, HGHT, UWIND, VWIND, UFLX, VFLX, WIND)


    return(dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt)


def proceed_timestep_upwind(GR, UFLXMP, VFLXMP, HGHT, TRACER,
                    dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt):

    # TIME STEPPING
    UFLXMP[GR.iijj] = UFLXMP[GR.iijj] + GR.dt*dUFLXMPdt
    VFLXMP[GR.iijj] = VFLXMP[GR.iijj] + GR.dt*dVFLXMPdt
    HGHT[GR.iijj] = HGHT[GR.iijj] + GR.dt*dHGHTdt
    TRACER[GR.iijj] = TRACER[GR.iijj] + GR.dt*dTRACERdt

    UFLXMP = exchange_BC(GR, UFLXMP)
    VFLXMP = exchange_BC(GR, VFLXMP)
    HGHT = exchange_BC(GR, HGHT)
    TRACER = exchange_BC(GR, TRACER)

    return(UFLXMP, VFLXMP, HGHT, TRACER)


def diagnose_fields_upwind(GR, HGHT, TRACER,
                    UWIND, VWIND, UFLXMP, VFLXMP, HSURF):

    # DIAGNOSTICS 
    UWIND[GR.iisjj] = ( UFLXMP[GR.iisjj_im1] + UFLXMP[GR.iisjj] ) \
                    / ( HGHT[GR.iisjj_im1] + HGHT[GR.iisjj] )
    VWIND[GR.iijjs] = ( VFLXMP[GR.iijjs_jm1] + VFLXMP[GR.iijjs] ) \
                    / ( HGHT[GR.iijjs_jm1] + HGHT[GR.iijjs] )
    UWIND = exchange_BC(GR, UWIND)
    VWIND = exchange_BC(GR, VWIND)

    return(UWIND, VWIND)
