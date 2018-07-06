import numpy as np
from namelist import inpRate, outRate, i_pseudo_radiation

def height_tendency_upstream(GR, HGHT, UWIND, VWIND, UFLX, VFLX):
    UFLX[GR.iisjj] = \
            np.maximum(UWIND[GR.iisjj],0) * HGHT[GR.iisjj_im1] * GR.dy + \
            np.minimum(UWIND[GR.iisjj],0) * HGHT[GR.iisjj] * GR.dy

    VFLX[GR.iijjs] = \
            np.maximum(VWIND[GR.iijjs],0) * HGHT[GR.iijjs_jm1] * GR.dx[GR.iijjs_jm1] + \
            np.minimum(VWIND[GR.iijjs],0) * HGHT[GR.iijjs] * GR.dx[GR.iijjs]


    if i_pseudo_radiation:
        radiation = - outRate*HGHT[GR.iijj] + inpRate*np.cos(GR.lat_rad[GR.iijj])

    if i_pseudo_radiation:
        dHGHTdt = ( - (UFLX[GR.iijj_ip1] - UFLX[GR.iijj]) - (VFLX[GR.iijj_jp1] - VFLX[GR.iijj]) ) / GR.A[GR.iijj] \
                 + radiation
    else:
        dHGHTdt = ( - (UFLX[GR.iijj_ip1] - UFLX[GR.iijj]) - (VFLX[GR.iijj_jp1] - VFLX[GR.iijj]) ) / GR.A[GR.iijj]

    return(dHGHTdt)
