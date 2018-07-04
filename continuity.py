import numpy as np

def diag_mass_flux(GR, UFLX, VFLX, UWIND, VWIND, HGHT):

    UFLX[GR.iisjj] = np.maximum(UWIND[GR.iisjj],0)*HGHT[GR.iisjj_im1] + \
            np.minimum(UWIND[GR.iisjj],0)*HGHT[GR.iisjj]

    VFLX[GR.iijjs] = np.maximum(VWIND[GR.iijjs],0)*HGHT[GR.iijjs_jm1] + \
            np.minimum(VWIND[GR.iijjs],0)*HGHT[GR.iijjs]

    return(UFLX, VFLX)
