import copy
import numpy as np
from height import height_tendency_jacobson
from wind import wind_tendency_jacobson
from tracer import tracer_tendency_jacobson
from boundaries import exchange_BC

def tendencies_jacobson(GR, HGHT, TRACER, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX):

    # PROGNOSE HGHT
    dHGHTdt, UFLX, VFLX = height_tendency_jacobson(GR, HGHT, UWIND, VWIND, UFLX, VFLX)

    # PROGNOSE WIND
    dUFLXdt, dVFLXdt = wind_tendency_jacobson(GR, UWIND, VWIND, UFLX, VFLX, 
                                                    HGHT, HSURF)

    # PROGNOSE TRACER
    dTRACERdt = tracer_tendency_jacobson(GR, TRACER, HGHT, UWIND, VWIND, UFLX, VFLX, WIND)


    return(dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt)


def proceed_timestep_jacobson(GR, UWIND, VWIND, HGHT,
                    TRACER,
                    dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt):

    # TIME STEPPING
    HGHT_OLD = copy.deepcopy(HGHT)
    HGHTA_is_OLD, HGHTA_js_OLD = interp_HGHTA(GR, HGHT_OLD)

    HGHT[GR.iijj] = HGHT[GR.iijj] + GR.dt*dHGHTdt
    HGHT = exchange_BC(GR, HGHT)
    HGHTA_is_NEW, HGHTA_js_NEW = interp_HGHTA(GR, HGHT)

    UWIND[GR.iisjj] = UWIND[GR.iisjj] * HGHTA_is_OLD/HGHTA_is_NEW \
                        + GR.dt*dUFLXdt/HGHTA_is_NEW
    VWIND[GR.iijjs] = VWIND[GR.iijjs] * HGHTA_js_OLD/HGHTA_js_NEW \
                        + GR.dt*dVFLXdt/HGHTA_js_NEW
    UWIND = exchange_BC(GR, UWIND)
    VWIND = exchange_BC(GR, VWIND)

    TRACER[GR.iijj] = TRACER[GR.iijj] + GR.dt*dTRACERdt
    TRACER = exchange_BC(GR, TRACER)

    return(UWIND, VWIND, HGHT, TRACER)


#def diagnose_fields_jacobson(GR, HGHT, TRACER,
#                    UWIND, VWIND, UFLX, VFLX, HSURF):
#
#    # DIAGNOSTICS 
#    UWIND[GR.iisjj] = UFLX[GR.iisjj] / (GR.A[GR.iisjj_im1]*HGHT[GR.iisjj_im1] + GR.A[GR.iisjj]*HGHT[GR.iisjj])/2
#    VWIND[GR.iijjs] = VFLX[GR.iijjs] / (GR.A[GR.iijjs_jm1]*HGHT[GR.iijjs_jm1] + GR.A[GR.iijjs]*HGHT[GR.iijjs])/2
#    UWIND = exchange_BC(GR, UWIND)
#    VWIND = exchange_BC(GR, VWIND)
#
#    return(UWIND, VWIND)




def interp_HGHTA(GR, HGHT):

    HGHTA_is = 1/8*(    HGHT[GR.iisjj_im1_jp1] * GR.A[GR.iisjj_im1_jp1] + \
                        HGHT[GR.iisjj_jp1]     * GR.A[GR.iisjj_jp1] + \
                    2 * HGHT[GR.iisjj_im1]     * GR.A[GR.iisjj_im1] + \
                    2 * HGHT[GR.iisjj]         * GR.A[GR.iisjj] + \
                        HGHT[GR.iisjj_im1_jm1] * GR.A[GR.iisjj_im1_jm1] + \
                        HGHT[GR.iisjj_jm1]     * GR.A[GR.iisjj_jm1]     )

    HGHTA_js = 1/8*(    HGHT[GR.iijjs_ip1_jm1] * GR.A[GR.iijjs_ip1_jm1] + \
                        HGHT[GR.iijjs_ip1]     * GR.A[GR.iijjs_ip1] + \
                    2 * HGHT[GR.iijjs_im1_jm1]     * GR.A[GR.iijjs_im1_jm1] + \
                    2 * HGHT[GR.iijjs]         * GR.A[GR.iijjs] + \
                        HGHT[GR.iijjs_im1_jm1] * GR.A[GR.iijjs_im1_jm1] + \
                        HGHT[GR.iijjs_im1]     * GR.A[GR.iijjs_im1]     )

    return(HGHTA_is, HGHTA_js)




