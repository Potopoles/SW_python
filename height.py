import numpy as np
from namelist import inpRate, outRate, i_pseudo_radiation

def height_tendency_upstream(GR, HGHT, UWIND, VWIND, UFLX, VFLX):
    UFLX[GR.iisjj] = \
            GR.dy * (np.maximum(UWIND[GR.iisjj],0) * HGHT[GR.iisjj_im1] + \
                        np.minimum(UWIND[GR.iisjj],0) * HGHT[GR.iisjj])

    VFLX[GR.iijjs] = \
            GR.dxjs[GR.iijjs] * ( np.maximum(VWIND[GR.iijjs],0) * HGHT[GR.iijjs_jm1] + \
                                    np.minimum(VWIND[GR.iijjs],0) * HGHT[GR.iijjs] )

    fluxdiv = ( - (UFLX[GR.iijj_ip1] - UFLX[GR.iijj]) - (VFLX[GR.iijj_jp1] - VFLX[GR.iijj]) ) / GR.A[GR.iijj]

    dHGHTdt = fluxdiv

    if i_pseudo_radiation:
        radiation = - outRate*HGHT[GR.iijj] + inpRate*np.cos(GR.lat_rad[GR.iijj])
        dHGHTdt += radiation

    return(dHGHTdt)




def height_tendency_jacobson(GR, HGHT, UWIND, VWIND, UFLX, VFLX):

    UFLX[GR.iisjj] = \
            (HGHT[GR.iisjj_im1] + HGHT[GR.iisjj])/2 * UWIND[GR.iisjj] * GR.dy \

    VFLX[GR.iijjs] = \
            (HGHT[GR.iijjs_jm1] + HGHT[GR.iijjs])/2 * VWIND[GR.iijjs] * GR.dxjs[GR.iijjs] \

    fluxdiv = ( - (UFLX[GR.iijj_ip1] - UFLX[GR.iijj]) - (VFLX[GR.iijj_jp1] - VFLX[GR.iijj]) ) / GR.A[GR.iijj]

    dHGHTdt = fluxdiv

    if i_pseudo_radiation:
        radiation = - outRate*HGHT[GR.iijj] + inpRate*np.cos(GR.lat_rad[GR.iijj])
        dHGHTdt += radiation

    return(dHGHTdt, UFLX, VFLX)
