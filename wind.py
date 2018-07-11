import numpy as np
from boundaries import exchange_BC
from constants import con_g, con_rE


def masspoint_flux_tendency_upstream(GR, UFLXMP, VFLXMP, HGHT,
                            UWIND, VWIND,
                            UUFLX, VUFLX, UVFLX, VVFLX,
                            HSURF):

    UFLXMP[GR.iijj] = HGHT[GR.iijj]*(UWIND[GR.iijj] + UWIND[GR.iijj_ip1])/2
    VFLXMP[GR.iijj] = HGHT[GR.iijj]*(VWIND[GR.iijj] + VWIND[GR.iijj_jp1])/2
    UFLXMP = exchange_BC(GR, UFLXMP)
    VFLXMP = exchange_BC(GR, VFLXMP)

    UUFLX[GR.iisjj] = GR.dy * ( \
            np.maximum(UWIND[GR.iisjj],0) * UFLXMP[GR.iisjj_im1] + \
            np.minimum(UWIND[GR.iisjj],0) * UFLXMP[GR.iisjj] )
    VUFLX[GR.iijjs] = GR.dxjs[GR.iijjs] * ( \
            np.maximum(VWIND[GR.iijjs],0) * UFLXMP[GR.iijjs_jm1] + \
            np.minimum(VWIND[GR.iijjs],0) * UFLXMP[GR.iijjs] )

    UVFLX[GR.iisjj] = GR.dy * ( \
            np.maximum(UWIND[GR.iisjj],0) * VFLXMP[GR.iisjj_im1] + \
            np.minimum(UWIND[GR.iisjj],0) * VFLXMP[GR.iisjj] )
    VVFLX[GR.iijjs] = GR.dxjs[GR.iijjs] * ( \
            np.maximum(VWIND[GR.iijjs],0) * VFLXMP[GR.iijjs_jm1] + \
            np.minimum(VWIND[GR.iijjs],0) * VFLXMP[GR.iijjs] )  

    corx = GR.corf[GR.iijj] * HGHT[GR.iijj] * (VWIND[GR.iijj] + VWIND[GR.iijj_jp1])/2
    cory = GR.corf[GR.iijj] * HGHT[GR.iijj] * (UWIND[GR.iijj] + UWIND[GR.iijj_ip1])/2

    dUFLXMPdt = - ( UUFLX[GR.iijj_ip1] - UUFLX[GR.iijj] + \
                    VUFLX[GR.iijj_jp1] - VUFLX[GR.iijj]) / GR.A[GR.iijj] + \
            - con_g*HGHT[GR.iijj]*( HSURF[GR.iijj_ip1] - HSURF[GR.iijj_im1] + \
                                    HGHT[GR.iijj_ip1] - HGHT[GR.iijj_im1] ) / \
                                    (2*GR.dx[GR.iijj]) + corx
    dVFLXMPdt = - ( UVFLX[GR.iijj_ip1] - UVFLX[GR.iijj] + \
                    VVFLX[GR.iijj_jp1] - VVFLX[GR.iijj]) / GR.A[GR.iijj] + \
            - con_g*HGHT[GR.iijj]*( HSURF[GR.iijj_jp1] - HSURF[GR.iijj_jm1] + \
                                    HGHT[GR.iijj_jp1] - HGHT[GR.iijj_jm1] ) / \
                                    (2*GR.dy) - cory

    return(dUFLXMPdt, dVFLXMPdt)




def wind_tendency_jacobson(GR, UWIND, VWIND, UFLX, VFLX, 
                                HGHT, HSURF):

    # HORIZONTAL ADVECTION
    horAdv_UWIND = horizontal_advection_UWIND(GR, UWIND, UFLX, VFLX)
    horAdv_VWIND = horizontal_advection_VWIND(GR, VWIND, UFLX, VFLX)

    # CORIOLIS AND SPHERICAL GRID CONVERSION (TODO)
    coriolis_UWIND, coriolis_VWIND = coriolis_term(GR, UWIND, VWIND, HGHT)

    # PRESSURE GRADIENT
    preGrad_UWIND, preGrad_VWIND = pressure_gradient_term(GR, HSURF, HGHT)

    dUFLXdt = horAdv_UWIND + coriolis_UWIND + preGrad_UWIND
    dVFLXdt = horAdv_VWIND + coriolis_VWIND + preGrad_VWIND

    return(dUFLXdt, dVFLXdt)



def coriolis_term(GR, UWIND, VWIND, HGHT):

    coriolis_UWIND = con_rE*GR.dlon_rad*GR.dlon_rad/2*(\

                          HGHT[GR.iisjj_im1] * \
                        ( VWIND[GR.iisjj_im1    ] + VWIND[GR.iisjj_im1_jp1] )/2 * \
                        ( GR.corf_is[GR.iisjj] * con_rE *\
                          np.cos(GR.latis_rad[GR.iisjj]) + \
                          ( UWIND[GR.iisjj_im1    ] + UWIND[GR.iisjj        ] )/2 * \
                          np.sin(GR.latis_rad[GR.iisjj]) )\

                        + HGHT[GR.iisjj    ] * \
                        ( VWIND[GR.iisjj        ] + VWIND[GR.iisjj_jp1    ] )/2 * \
                        ( GR.corf_is[GR.iisjj] * con_rE * \
                          np.cos(GR.latis_rad[GR.iisjj]) + \
                          ( UWIND[GR.iisjj        ] + UWIND[GR.iisjj_ip1    ] )/2 * \
                          np.sin(GR.latis_rad[GR.iisjj]) )\
                        )

    coriolis_VWIND = - con_rE*GR.dlon_rad*GR.dlon_rad/2*(\

                          HGHT[GR.iijjs_jm1] * \
                        ( UWIND[GR.iijjs_jm1    ] + UWIND[GR.iijjs_ip1_jm1] )/2 * \
                        ( GR.corf_js[GR.iijjs_jm1] * con_rE *\
                          np.cos(GR.latjs_rad[GR.iijjs_jm1]) +\
                          ( UWIND[GR.iijjs_jm1    ] + UWIND[GR.iijjs_ip1_jm1] )/2 * \
                          np.sin(GR.latjs_rad[GR.iijjs_jm1]) )\

                        + HGHT[GR.iijjs    ] * \
                        ( UWIND[GR.iijjs        ] + UWIND[GR.iijjs_ip1    ] )/2 * \
                        ( GR.corf_js[GR.iijjs] * con_rE *\
                          np.cos(GR.latjs_rad[GR.iijjs]) +\
                          ( UWIND[GR.iijjs        ] + UWIND[GR.iijjs_ip1    ] )/2 * \
                          np.sin(GR.latjs_rad[GR.iijjs]) )\
                        )

    return(coriolis_UWIND, coriolis_VWIND)

def pressure_gradient_term(GR, HSURF, HGHT):

    preGrad_UWIND = - con_g*GR.dy*\
                        ( HGHT[GR.iisjj]  + HGHT[GR.iisjj_im1] )/2*\
                        ( HSURF[GR.iisjj] - HSURF[GR.iisjj_im1] + \
                          HGHT[GR.iisjj]  - HGHT[GR.iisjj_im1]  )

    preGrad_VWIND = - con_g*GR.dx[GR.iijjs]*\
                        ( HGHT[GR.iijjs]  + HGHT[GR.iijjs_jm1] )/2*\
                        ( HSURF[GR.iijjs] - HSURF[GR.iijjs_jm1] + \
                          HGHT[GR.iijjs]  - HGHT[GR.iijjs_jm1]  )

    return(preGrad_UWIND, preGrad_VWIND)

def horizontal_advection_UWIND(GR, UWIND, UFLX, VFLX):
    BFLX = np.full( (GR.nx +2*GR.nb,GR.ny +2*GR.nb), np.nan )
    CFLX = np.full( (GR.nxs+2*GR.nb,GR.nys+2*GR.nb), np.nan )
    DFLX = np.full( (GR.nx +2*GR.nb,GR.nys+2*GR.nb), np.nan )
    EFLX = np.full( (GR.nx +2*GR.nb,GR.nys+2*GR.nb), np.nan )

    BFLX[GR.iijj]   = 1/12 * (  UFLX[GR.iijj_jm1    ]   +   UFLX[GR.iijj_ip1_jm1]   +\
                            2*( UFLX[GR.iijj        ]   +   UFLX[GR.iijj_ip1    ] ) +\
                                UFLX[GR.iijj_jp1    ]   +   UFLX[GR.iijj_ip1_jp1]    )
    BFLX = exchange_BC(GR, BFLX)

    CFLX[GR.iisjjs] = 1/12 * (  VFLX[GR.iisjjs_im1_jm1]   +   VFLX[GR.iisjjs_jm1    ]   +\
                            2*( VFLX[GR.iisjjs_im1    ]   +   VFLX[GR.iisjjs        ] ) +\
                                VFLX[GR.iisjjs_im1_jp1]   +   VFLX[GR.iisjjs_jp1    ]    )
    CFLX = exchange_BC(GR, CFLX)

    DFLX[GR.iijjs]  = 1/24 * (  VFLX[GR.iijjs_jm1    ]   + 2*VFLX[GR.iijjs        ]   +\
                                VFLX[GR.iijjs_jp1    ]   +   UFLX[GR.iijjs_jm1    ]   +\
                                UFLX[GR.iijjs        ]   +   UFLX[GR.iijjs_ip1_jm1]   +\
                                UFLX[GR.iijjs_ip1    ]                                  )
    DFLX = exchange_BC(GR, DFLX)

    EFLX[GR.iijjs]  = 1/24 * (  VFLX[GR.iijjs_jm1    ]   + 2*VFLX[GR.iijjs        ]   +\
                                VFLX[GR.iijjs_jp1    ]   -   UFLX[GR.iijjs_jm1    ]   +\
                              - UFLX[GR.iijjs        ]   -   UFLX[GR.iijjs_ip1_jm1]   +\
                              - UFLX[GR.iijjs_ip1    ]                                  )
    EFLX = exchange_BC(GR, EFLX)
    #print(DFLX)
    #quit()


    horAdv_UWIND =  + BFLX [GR.iisjj_im1    ] * \
                    ( UWIND[GR.iisjj_im1    ] + UWIND[GR.iisjj        ] )/2 \
                    - BFLX [GR.iisjj        ] * \
                    ( UWIND[GR.iisjj        ] + UWIND[GR.iisjj_ip1    ] )/2 \
                    \
                    + CFLX [GR.iisjj        ] * \
                    ( UWIND[GR.iisjj_jm1    ] + UWIND[GR.iisjj        ] )/2 \
                    - CFLX [GR.iisjj_jp1    ] * \
                    ( UWIND[GR.iisjj        ] + UWIND[GR.iisjj_jp1    ] )/2 \
                    \
                    + DFLX [GR.iisjj_im1    ] * \
                    ( UWIND[GR.iisjj_im1_jm1] + UWIND[GR.iisjj        ] )/2 \
                    - DFLX [GR.iisjj_jp1    ] * \
                    ( UWIND[GR.iisjj        ] + UWIND[GR.iisjj_ip1_jp1] )/2 \
                    \
                    + EFLX [GR.iisjj        ] * \
                    ( UWIND[GR.iisjj_ip1_jm1] + UWIND[GR.iisjj        ] )/2 \
                    - EFLX [GR.iisjj_im1_jp1] * \
                    ( UWIND[GR.iisjj        ] + UWIND[GR.iisjj_im1_jp1] )/2 

    #print(horAdv_UWIND)
    #quit()
    return( horAdv_UWIND )



def horizontal_advection_VWIND(GR, VWIND, UFLX, VFLX):
    RFLX = np.zeros( (GR.nx +2*GR.nb,GR.ny +2*GR.nb) )
    QFLX = np.zeros( (GR.nxs+2*GR.nb,GR.nys+2*GR.nb) )
    SFLX = np.zeros( (GR.nxs+2*GR.nb,GR.ny +2*GR.nb) )
    TFLX = np.zeros( (GR.nxs+2*GR.nb,GR.ny +2*GR.nb) )

    RFLX[GR.iijj] = 1/12 * (    VFLX[GR.iijj_im1    ]   +   VFLX[GR.iijj_im1_jp1]   +\
                            2*( VFLX[GR.iijj        ]   +   VFLX[GR.iijj_jp1    ] ) +\
                                VFLX[GR.iijj_ip1    ]   +   VFLX[GR.iijj_ip1_jp1]    )
    RFLX = exchange_BC(GR, RFLX)

    QFLX[GR.iisjjs] = 1/12 * (  UFLX[GR.iisjjs_im1_jm1]   +   UFLX[GR.iisjjs_im1    ]   +\
                            2*( UFLX[GR.iisjjs_jm1    ]   +   UFLX[GR.iisjjs        ] ) +\
                                UFLX[GR.iisjjs_ip1_jm1]   +   UFLX[GR.iisjjs_ip1    ]    )
    QFLX = exchange_BC(GR, QFLX)

    SFLX[GR.iisjj]  = 1/24 * (  VFLX[GR.iisjj_im1     ]   +   VFLX[GR.iisjj_im1_jp1]   +\
                                VFLX[GR.iisjj         ]   +   VFLX[GR.iisjj_jp1    ]   +\
                                UFLX[GR.iisjj_im1     ]   + 2*UFLX[GR.iisjj        ]   +\
                                UFLX[GR.iisjj_ip1     ]                                  )
    SFLX = exchange_BC(GR, SFLX)

    TFLX[GR.iisjj]  = 1/24 * (  VFLX[GR.iisjj_im1    ]   +   VFLX[GR.iisjj_im1_jp1]   +\
                                VFLX[GR.iisjj        ]   +   VFLX[GR.iisjj_jp1    ]   +\
                              - UFLX[GR.iisjj_im1    ]   - 2*UFLX[GR.iisjj        ]   +\
                              - UFLX[GR.iisjj_ip1    ]                                  )
    TFLX = exchange_BC(GR, TFLX)


    horAdv_VWIND =  + RFLX [GR.iijjs_jm1    ] * \
                    ( VWIND[GR.iijjs_jm1    ] + VWIND[GR.iijjs        ] )/2 \
                    - RFLX [GR.iijjs        ] * \
                    ( VWIND[GR.iijjs        ] + VWIND[GR.iijjs_jp1    ] )/2 \
                    \
                    + QFLX [GR.iijjs        ] * \
                    ( VWIND[GR.iijjs_im1    ] + VWIND[GR.iijjs        ] )/2 \
                    - QFLX [GR.iijjs_ip1    ] * \
                    ( VWIND[GR.iijjs        ] + VWIND[GR.iijjs_ip1    ] )/2 \
                    \
                    + SFLX [GR.iijjs_jm1    ] * \
                    ( VWIND[GR.iijjs_im1_jm1] + VWIND[GR.iijjs        ] )/2 \
                    - SFLX [GR.iijjs_ip1    ] * \
                    ( VWIND[GR.iijjs        ] + VWIND[GR.iijjs_ip1_jp1] )/2 \
                    \
                    + TFLX [GR.iijjs_ip1_jm1] * \
                    ( VWIND[GR.iijjs_ip1_jm1] + VWIND[GR.iijjs        ] )/2 \
                    - TFLX [GR.iijjs        ] * \
                    ( VWIND[GR.iijjs        ] + VWIND[GR.iijjs_im1_jp1] )/2 
    #print(horAdv_VWIND)
    #quit()

    return( horAdv_VWIND )
