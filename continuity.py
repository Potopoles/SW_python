import numpy as np

def colp_tendency_upstream(GR, COLP, UWIND, VWIND, UFLX, VFLX):
    UFLX[GR.iisjj] = \
            np.maximum(UWIND[GR.iisjj],0) * COLP[GR.iisjj_im1] * GR.dy + \
            np.minimum(UWIND[GR.iisjj],0) * COLP[GR.iisjj] * GR.dy

    VFLX[GR.iijjs] = \
            np.maximum(VWIND[GR.iijjs],0) * COLP[GR.iijjs_jm1] * GR.dx[GR.iijjs_jm1] + \
            np.minimum(VWIND[GR.iijjs],0) * COLP[GR.iijjs] * GR.dx[GR.iijjs]

    dCOLPdt = ( - (UFLX[GR.iijj_ip1] - UFLX[GR.iijj]) - (VFLX[GR.iijj_jp1] - VFLX[GR.iijj]) ) / GR.A[GR.iijj]

    return(dCOLPdt)
