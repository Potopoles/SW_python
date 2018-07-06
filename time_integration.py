import copy
import numpy as np
from height import height_tendency_upstream, height_tendency_jacobson
from wind import masspoint_flux_tendency_upstream, wind_tendency_jacobson
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


def tendencies_jacobson(GR, HGHT, TRACER, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX):

    # PROGNOSE HGHT
    dHGHTdt, UFLX, VFLX = height_tendency_jacobson(GR, HGHT, UWIND, VWIND, UFLX, VFLX)

    # PROGNOSE WIND
    dUFLXMPdt, dVFLXMPdt = wind_tendency_jacobson(GR, UFLXMP, VFLXMP, HGHT,
                                                    UWIND, VWIND,
                                                    UUFLX, VUFLX, UVFLX, VVFLX,
                                                    HSURF)

    # PROGNOSE TRACER
    dTRACERdt = tracer_tendency_upstream(GR, TRACER, HGHT, UWIND, VWIND, UFLX, VFLX, WIND)


    return(dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt)


def proceed_timestep(GR, UFLXMP, VFLXMP, HGHT, TRACER,
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


def diagnose_fields(GR, HGHT, TRACER,
                    UWIND, VWIND, UFLXMP, VFLXMP, HSURF):

    # DIAGNOSTICS 
    UWIND[GR.iisjj] = ( UFLXMP[GR.iisjj_im1] + UFLXMP[GR.iisjj] ) \
                    / ( HGHT[GR.iisjj_im1] + HGHT[GR.iisjj] )
    VWIND[GR.iijjs] = ( VFLXMP[GR.iijjs_jm1] + VFLXMP[GR.iijjs] ) \
                    / ( HGHT[GR.iijjs_jm1] + HGHT[GR.iijjs] )
    UWIND = exchange_BC(GR, UWIND)
    VWIND = exchange_BC(GR, VWIND)

    return(UWIND, VWIND)






def euler_forward(GR, HGHT, TRACER,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX,
                    HSURF):

    #dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
    dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_jacobson(GR, 
                    HGHT, TRACER, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX)

    UFLXMP, VFLXMP, HGHT, TRACER = proceed_timestep(GR, UFLXMP, VFLXMP, HGHT, TRACER,
                                                dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt)

    UWIND, VWIND = diagnose_fields(GR, HGHT, TRACER,
                    UWIND, VWIND, UFLXMP, VFLXMP, HSURF)

    return(HGHT, TRACER,
            UWIND, VWIND,
            UFLX, VFLX, UFLXMP, VFLXMP,
            UUFLX, UVFLX, VUFLX, VVFLX,
            HSURF)





def matsuno(GR, HGHT, TRACER,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX,
                    HSURF):

    ########## ESTIMATE
    #dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
    dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_jacobson(GR, 
                    HGHT, TRACER, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX)

    # has to happen after masspoint_flux_tendency function
    UFLXMP_OLD = copy.deepcopy(UFLXMP)
    VFLXMP_OLD = copy.deepcopy(VFLXMP)
    HGHT_OLD = copy.deepcopy(HGHT)
    TRACER_OLD = copy.deepcopy(TRACER)

    UFLXMP, VFLXMP, HGHT, TRACER = proceed_timestep(GR, UFLXMP, VFLXMP, HGHT, TRACER,
                                                dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt)

    UWIND, VWIND = diagnose_fields(GR, HGHT, TRACER,
                    UWIND, VWIND, UFLXMP, VFLXMP, HSURF)

    ########## FINAL
    #dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
    dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_jacobson(GR, 
                    HGHT, TRACER, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX)

    UFLXMP, VFLXMP, HGHT, TRACER = proceed_timestep(GR, UFLXMP_OLD, VFLXMP_OLD, HGHT_OLD, TRACER_OLD,
                                                dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt)

    UWIND, VWIND = diagnose_fields(GR, HGHT, TRACER,
                    UWIND, VWIND, UFLXMP, VFLXMP, HSURF)

    return(HGHT, TRACER,
            UWIND, VWIND,
            UFLX, VFLX, UFLXMP, VFLXMP,
            UUFLX, UVFLX, VUFLX, VVFLX,
            HSURF)







def RK4(GR, HGHT, TRACER,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX,
                    HSURF):

    ########## level 1
    dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies(GR, 
                    HGHT, TRACER, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX)

    # has to happen after masspoint_flux_tendency function
    UFLXMP_OLD = copy.deepcopy(UFLXMP)
    VFLXMP_OLD = copy.deepcopy(VFLXMP)
    HGHT_OLD = copy.deepcopy(HGHT)
    TRACER_OLD = copy.deepcopy(TRACER)

    UFLXMP_INT = copy.deepcopy(UFLXMP)
    VFLXMP_INT = copy.deepcopy(VFLXMP)
    HGHT_INT = copy.deepcopy(HGHT)
    TRACER_INT = copy.deepcopy(TRACER)

    dUFLXMP = GR.dt*dUFLXMPdt
    dVFLXMP = GR.dt*dVFLXMPdt
    dHGHT = GR.dt*dHGHTdt
    dTRACER = GR.dt*dTRACERdt

    UFLXMP_INT[GR.iijj] = UFLXMP_OLD[GR.iijj] + dUFLXMP/2
    VFLXMP_INT[GR.iijj] = VFLXMP_OLD[GR.iijj] + dVFLXMP/2
    HGHT_INT[GR.iijj] = HGHT_OLD[GR.iijj] + dHGHT/2
    TRACER_INT[GR.iijj] = TRACER_OLD[GR.iijj] + dTRACER/2
    UFLXMP_INT = exchange_BC(GR, UFLXMP_INT)
    VFLXMP_INT = exchange_BC(GR, VFLXMP_INT)
    HGHT_INT = exchange_BC(GR, HGHT_INT)
    TRACER_INT = exchange_BC(GR, TRACER_INT)

    UFLXMP[GR.iijj] = UFLXMP_OLD[GR.iijj] + dUFLXMP/6
    VFLXMP[GR.iijj] = VFLXMP_OLD[GR.iijj] + dVFLXMP/6
    HGHT[GR.iijj] = HGHT_OLD[GR.iijj] + dHGHT/6
    TRACER[GR.iijj] = TRACER_OLD[GR.iijj] + dTRACER/6
    
    UWIND, VWIND = diagnose_fields(GR, HGHT_INT, TRACER_INT,
                    UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

    ########## level 2
    dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies(GR, 
                    HGHT_INT, TRACER_INT, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP_INT, VFLXMP_INT,
                    UUFLX, UVFLX, VUFLX, VVFLX)

    dUFLXMP = GR.dt*dUFLXMPdt
    dVFLXMP = GR.dt*dVFLXMPdt
    dHGHT = GR.dt*dHGHTdt
    dTRACER = GR.dt*dTRACERdt

    UFLXMP_INT[GR.iijj] = UFLXMP_OLD[GR.iijj] + dUFLXMP/2
    VFLXMP_INT[GR.iijj] = VFLXMP_OLD[GR.iijj] + dVFLXMP/2
    HGHT_INT[GR.iijj] = HGHT_OLD[GR.iijj] + dHGHT/2
    TRACER_INT[GR.iijj] = TRACER_OLD[GR.iijj] + dTRACER/2
    UFLXMP_INT = exchange_BC(GR, UFLXMP_INT)
    VFLXMP_INT = exchange_BC(GR, VFLXMP_INT)
    HGHT_INT = exchange_BC(GR, HGHT_INT)
    TRACER_INT = exchange_BC(GR, TRACER_INT)

    UFLXMP[GR.iijj] = UFLXMP[GR.iijj] + dUFLXMP/3
    VFLXMP[GR.iijj] = VFLXMP[GR.iijj] + dVFLXMP/3
    HGHT[GR.iijj] = HGHT[GR.iijj] + dHGHT/3
    TRACER[GR.iijj] = TRACER[GR.iijj] + dTRACER/3
    
    UWIND, VWIND = diagnose_fields(GR, HGHT_INT, TRACER_INT,
                    UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

    ########## level 3
    dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies(GR, 
                    HGHT_INT, TRACER_INT, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP_INT, VFLXMP_INT,
                    UUFLX, UVFLX, VUFLX, VVFLX)

    dUFLXMP = GR.dt*dUFLXMPdt
    dVFLXMP = GR.dt*dVFLXMPdt
    dHGHT = GR.dt*dHGHTdt
    dTRACER = GR.dt*dTRACERdt

    UFLXMP_INT[GR.iijj] = UFLXMP_OLD[GR.iijj] + dUFLXMP
    VFLXMP_INT[GR.iijj] = VFLXMP_OLD[GR.iijj] + dVFLXMP
    HGHT_INT[GR.iijj] = HGHT_OLD[GR.iijj] + dHGHT
    TRACER_INT[GR.iijj] = TRACER_OLD[GR.iijj] + dTRACER
    UFLXMP_INT = exchange_BC(GR, UFLXMP_INT)
    VFLXMP_INT = exchange_BC(GR, VFLXMP_INT)
    HGHT_INT = exchange_BC(GR, HGHT_INT)
    TRACER_INT = exchange_BC(GR, TRACER_INT)

    UFLXMP[GR.iijj] = UFLXMP[GR.iijj] + dUFLXMP/3
    VFLXMP[GR.iijj] = VFLXMP[GR.iijj] + dVFLXMP/3
    HGHT[GR.iijj] = HGHT[GR.iijj] + dHGHT/3
    TRACER[GR.iijj] = TRACER[GR.iijj] + dTRACER/3
    
    UWIND, VWIND = diagnose_fields(GR, HGHT_INT, TRACER_INT,
                    UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

    ########## level 4
    dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies(GR, 
                    HGHT_INT, TRACER_INT, HSURF,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP_INT, VFLXMP_INT,
                    UUFLX, UVFLX, VUFLX, VVFLX)

    dUFLXMP = GR.dt*dUFLXMPdt
    dVFLXMP = GR.dt*dVFLXMPdt
    dHGHT = GR.dt*dHGHTdt
    dTRACER = GR.dt*dTRACERdt

    UFLXMP[GR.iijj] = UFLXMP[GR.iijj] + dUFLXMP/6
    VFLXMP[GR.iijj] = VFLXMP[GR.iijj] + dVFLXMP/6
    HGHT[GR.iijj] = HGHT[GR.iijj] + dHGHT/6
    TRACER[GR.iijj] = TRACER[GR.iijj] + dTRACER/6
    UFLXMP = exchange_BC(GR, UFLXMP)
    VFLXMP = exchange_BC(GR, VFLXMP)
    HGHT = exchange_BC(GR, HGHT)
    TRACER = exchange_BC(GR, TRACER)
    
    UWIND, VWIND = diagnose_fields(GR, HGHT, TRACER,
                    UWIND, VWIND, UFLXMP, VFLXMP, HSURF)

    return(HGHT, TRACER,
            UWIND, VWIND,
            UFLX, VFLX, UFLXMP, VFLXMP,
            UUFLX, UVFLX, VUFLX, VVFLX,
            HSURF)
