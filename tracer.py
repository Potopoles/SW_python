import numpy as np

def tracer_tendency_upstream(GR, TRACER, HGHT, UWIND, VWIND, UFLX, VFLX, WIND):
    UFLX[GR.iisjj] = \
            GR.dy * (np.maximum(UWIND[GR.iisjj],0) * TRACER[GR.iisjj_im1] * \
                    HGHT[GR.iisjj_im1] + \
                    np.minimum(UWIND[GR.iisjj],0) * TRACER[GR.iisjj] * \
                    HGHT[GR.iisjj])

    VFLX[GR.iijjs] = \
            GR.dxjs[GR.iijjs] * ( np.maximum(VWIND[GR.iijjs],0) * TRACER[GR.iijjs_jm1] * \
                                    HGHT[GR.iijjs_jm1] + \
                                    np.minimum(VWIND[GR.iijjs],0) * TRACER[GR.iijjs] * \
                                    HGHT[GR.iijjs] )

    dTRACERdt = ( - (UFLX[GR.iijj_ip1] - UFLX[GR.iijj]) - \
                    (VFLX[GR.iijj_jp1] - VFLX[GR.iijj]) ) / \
            (GR.A[GR.iijj] * HGHT[GR.iijj]) + \
            np.maximum(0., 0.00001*WIND[GR.iijj]*(5.0 - TRACER[GR.iijj]))

    return(dTRACERdt)





def tracer_tendency_jacobson(GR, TRACER, HGHT, UWIND, VWIND, UFLX, VFLX, WIND):

    dTRACERdt = ( + UFLX[GR.iijj    ] *\
                     (TRACER[GR.iijj_im1] + TRACER[GR.iijj    ])/2 \
                  - UFLX[GR.iijj_ip1] *\
                     (TRACER[GR.iijj    ] + TRACER[GR.iijj_ip1])/2 \
                  + VFLX[GR.iijj    ] *\
                     (TRACER[GR.iijj_jm1] + TRACER[GR.iijj    ])/2 \
                  - VFLX[GR.iijj_jp1] *\
                     (TRACER[GR.iijj    ] + TRACER[GR.iijj_jp1])/2 \
                ) / (GR.A[GR.iijj] * HGHT[GR.iijj]) #\
                #+ np.maximum(0., 0.00001*WIND[GR.iijj]*(5.0 - TRACER[GR.iijj]))

    return(dTRACERdt)

