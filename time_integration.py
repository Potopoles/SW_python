import copy
import numpy as np
from boundaries import exchange_BC
from upwind import tendencies_upwind, proceed_timestep_upwind, diagnose_fields_upwind
from jacobson import tendencies_jacobson, proceed_timestep_jacobson#, diagnose_fields_jacobson
from jacobson import interp_HGHTA

######################################################################################
######################################################################################
######################################################################################

def euler_forward(GR, HGHT, TRACER,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX,
                    HSURF, i_spatial_discretization):

    if i_spatial_discretization == 'UPWIND':
        dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
                        HGHT, TRACER, HSURF,
                        UWIND, VWIND, WIND,
                        UFLX, VFLX, UFLXMP, VFLXMP,
                        UUFLX, UVFLX, VUFLX, VVFLX)

        UFLXMP, VFLXMP, HGHT, TRACER = proceed_timestep_upwind(GR, UFLXMP, VFLXMP, 
                                            HGHT, TRACER, dHGHTdt, dUFLXMPdt, 
                                            dVFLXMPdt, dTRACERdt)

        UWIND, VWIND = diagnose_fields_upwind(GR, HGHT, TRACER,
                        UWIND, VWIND, UFLXMP, VFLXMP, HSURF)

    elif i_spatial_discretization == 'JACOBSON':
        dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt = tendencies_jacobson(GR, HGHT, TRACER, HSURF,
                                                                    UWIND, VWIND, WIND,
                                                                    UFLX, VFLX)

        UWIND, VWIND, HGHT, TRACER = proceed_timestep_jacobson(GR, UWIND, VWIND, HGHT,
                                                    TRACER,
                                                    dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt)

        #UWIND, VWIND = diagnose_fields_jacobson(GR, HGHT, TRACER,
        #                UWIND, VWIND, UFLX, VFLX, HSURF)


    return(HGHT, TRACER,
            UWIND, VWIND,
            UFLX, VFLX, UFLXMP, VFLXMP,
            UUFLX, UVFLX, VUFLX, VVFLX,
            HSURF)




######################################################################################
######################################################################################
######################################################################################

def matsuno(GR, HGHT, TRACER,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX,
                    HSURF, i_spatial_discretization):

    if i_spatial_discretization == 'UPWIND':
        ########## ESTIMATE
        dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
                        HGHT, TRACER, HSURF,
                        UWIND, VWIND, WIND,
                        UFLX, VFLX, UFLXMP, VFLXMP,
                        UUFLX, UVFLX, VUFLX, VVFLX)

        # has to happen after masspoint_flux_tendency function
        UFLXMP_OLD = copy.deepcopy(UFLXMP)
        VFLXMP_OLD = copy.deepcopy(VFLXMP)
        HGHT_OLD = copy.deepcopy(HGHT)
        TRACER_OLD = copy.deepcopy(TRACER)


        UFLXMP, VFLXMP, HGHT, TRACER = proceed_timestep_upwind(GR, UFLXMP, VFLXMP, HGHT,
                                            TRACER, dHGHTdt, dUFLXMPdt, dVFLXMPdt,
                                            dTRACERdt)

        UWIND, VWIND = diagnose_fields_upwind(GR, HGHT, TRACER,
                        UWIND, VWIND, UFLXMP, VFLXMP, HSURF)

        ########## FINAL
        dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
                        HGHT, TRACER, HSURF,
                        UWIND, VWIND, WIND,
                        UFLX, VFLX, UFLXMP, VFLXMP,
                        UUFLX, UVFLX, VUFLX, VVFLX)

        UFLXMP, VFLXMP, HGHT, TRACER = proceed_timestep_upwind(GR, UFLXMP_OLD,
                                            VFLXMP_OLD, HGHT_OLD, TRACER_OLD,
                                            dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt)

        UWIND, VWIND = diagnose_fields_upwind(GR, HGHT, TRACER,
                        UWIND, VWIND, UFLXMP, VFLXMP, HSURF)


    elif i_spatial_discretization == 'JACOBSON':
        ########## ESTIMATE
        dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt = tendencies_jacobson(GR, HGHT, TRACER, HSURF,
                                                                    UWIND, VWIND, WIND,
                                                                    UFLX, VFLX)

        # has to happen after masspoint_flux_tendency function
        UWIND_OLD = copy.deepcopy(UWIND)
        VWIND_OLD = copy.deepcopy(VWIND)
        HGHT_OLD = copy.deepcopy(HGHT)
        TRACER_OLD = copy.deepcopy(TRACER)

        UWIND, VWIND, HGHT, TRACER = proceed_timestep_jacobson(GR, UWIND, VWIND, HGHT,
                                                    TRACER,
                                                    dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt)

        #UWIND, VWIND = diagnose_fields_jacobson(GR, HGHT, TRACER,
        #                UWIND, VWIND, UFLX, VFLX, HSURF)

        ########## FINAL
        dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt = tendencies_jacobson(GR, HGHT, TRACER, HSURF,
                                                                    UWIND, VWIND, WIND,
                                                                    UFLX, VFLX)

        UWIND, VWIND, HGHT, TRACER = proceed_timestep_jacobson(GR, UWIND_OLD, VWIND_OLD,
                                                    HGHT_OLD, TRACER_OLD,
                                                    dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt)

        #UWIND, VWIND = diagnose_fields_jacobson(GR, HGHT, TRACER,
        #                UWIND, VWIND, UFLX, VFLX, HSURF)

    return(HGHT, TRACER,
            UWIND, VWIND,
            UFLX, VFLX, UFLXMP, VFLXMP,
            UUFLX, UVFLX, VUFLX, VVFLX,
            HSURF)






######################################################################################
######################################################################################
######################################################################################

def RK4(GR, HGHT, TRACER,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX,
                    HSURF, i_spatial_discretization):

    if i_spatial_discretization == 'UPWIND':
        ########## level 1
        dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
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
        
        UWIND, VWIND = diagnose_fields_upwind(GR, HGHT_INT, TRACER_INT,
                        UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

        ########## level 2
        dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
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
        
        UWIND, VWIND = diagnose_fields_upwind(GR, HGHT_INT, TRACER_INT,
                        UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

        ########## level 3
        dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
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
        
        UWIND, VWIND = diagnose_fields_upwind(GR, HGHT_INT, TRACER_INT,
                        UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

        ########## level 4
        dHGHTdt, dUFLXMPdt, dVFLXMPdt, dTRACERdt = tendencies_upwind(GR, 
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
        
        UWIND, VWIND = diagnose_fields_upwind(GR, HGHT, TRACER,
                        UWIND, VWIND, UFLXMP, VFLXMP, HSURF)


    elif i_spatial_discretization == 'JACOBSON':
        ########## level 1
        dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt = tendencies_jacobson(GR, HGHT, TRACER, HSURF,
                                                                    UWIND, VWIND, WIND,
                                                                    UFLX, VFLX)

        # has to happen after masspoint_flux_tendency function
        UWIND_START = copy.deepcopy(UWIND)
        VWIND_START = copy.deepcopy(VWIND)
        HGHT_START = copy.deepcopy(HGHT)
        TRACER_START = copy.deepcopy(TRACER)

        UWIND_INT = copy.deepcopy(UWIND)
        VWIND_INT = copy.deepcopy(VWIND)
        HGHT_INT = copy.deepcopy(HGHT)
        TRACER_INT = copy.deepcopy(TRACER)

        dUFLX = GR.dt*dUFLXdt
        dVFLX = GR.dt*dVFLXdt
        dHGHT = GR.dt*dHGHTdt
        dTRACER = GR.dt*dTRACERdt

        # TIME STEPPING
        HGHTA_is_START, HGHTA_js_START = interp_HGHTA(GR, HGHT_START)
        HGHT_INT[GR.iijj] = HGHT_START[GR.iijj] + dHGHT/2
        HGHT_INT = exchange_BC(GR, HGHT_INT)
        HGHTA_is_NEW, HGHTA_js_NEW = interp_HGHTA(GR, HGHT_INT)
        UWIND_INT[GR.iisjj] = UWIND_START[GR.iisjj] * HGHTA_is_START/HGHTA_is_NEW \
                            + dUFLX/2/HGHTA_is_NEW
        VWIND_INT[GR.iijjs] = VWIND_START[GR.iijjs] * HGHTA_js_START/HGHTA_js_NEW \
                            + dVFLX/2/HGHTA_js_NEW
        UWIND_INT = exchange_BC(GR, UWIND_INT)
        VWIND_INT = exchange_BC(GR, VWIND_INT)
        TRACER_INT[GR.iijj] = TRACER_START[GR.iijj] + dTRACER/2
        TRACER_INT = exchange_BC(GR, TRACER_INT)

        HGHTA_is_START, HGHTA_js_START = interp_HGHTA(GR, HGHT_START)
        HGHT[GR.iijj] = HGHT_START[GR.iijj] + dHGHT/6
        HGHT = exchange_BC(GR, HGHT)
        HGHTA_is_NEW, HGHTA_js_NEW = interp_HGHTA(GR, HGHT)
        UWIND[GR.iisjj] = UWIND_START[GR.iisjj] * HGHTA_is_START/HGHTA_is_NEW \
                        + dUFLX/6/HGHTA_is_NEW
        VWIND[GR.iijjs] = VWIND_START[GR.iijjs] * HGHTA_js_START/HGHTA_js_NEW \
                            + dVFLX/6/HGHTA_js_NEW
        UWIND = exchange_BC(GR, UWIND)
        VWIND = exchange_BC(GR, VWIND)
        TRACER[GR.iijj] = TRACER_START[GR.iijj] + dTRACER/6
        TRACER = exchange_BC(GR, TRACER)

        #UWIND, VWIND = diagnose_fields_jacobson(GR, HGHT_INT, TRACER_INT,
        #                UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

        ########## level 2
        dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt = tendencies_jacobson(GR, HGHT_INT, TRACER_INT,
                                                    HSURF, UWIND_INT, VWIND_INT, WIND,
                                                    UFLX, VFLX)

        dUFLX = GR.dt*dUFLXdt
        dVFLX = GR.dt*dVFLXdt
        dHGHT = GR.dt*dHGHTdt
        dTRACER = GR.dt*dTRACERdt

        HGHT_INT[GR.iijj] = HGHT_START[GR.iijj] + dHGHT/2
        HGHT_INT = exchange_BC(GR, HGHT_INT)
        HGHTA_is_NEW, HGHTA_js_NEW = interp_HGHTA(GR, HGHT_INT)
        UWIND_INT[GR.iisjj] = UWIND_START[GR.iisjj] * HGHTA_is_START/HGHTA_is_NEW \
                            + dUFLX/2/HGHTA_is_NEW
        VWIND_INT[GR.iijjs] = VWIND_START[GR.iijjs] * HGHTA_js_START/HGHTA_js_NEW \
                            + dVFLX/2/HGHTA_js_NEW
        UWIND_INT = exchange_BC(GR, UWIND_INT)
        VWIND_INT = exchange_BC(GR, VWIND_INT)
        TRACER_INT[GR.iijj] = TRACER_START[GR.iijj] + dTRACER/2
        TRACER_INT = exchange_BC(GR, TRACER_INT)

        HGHT_OLD = copy.deepcopy(HGHT)
        HGHTA_is_OLD, HGHTA_js_OLD = interp_HGHTA(GR, HGHT_OLD)
        HGHT[GR.iijj] = HGHT[GR.iijj] + dHGHT/3
        HGHT = exchange_BC(GR, HGHT)
        HGHTA_is_NEW, HGHTA_js_NEW = interp_HGHTA(GR, HGHT)
        UWIND[GR.iisjj] = UWIND[GR.iisjj] * HGHTA_is_OLD/HGHTA_is_NEW \
                        + dUFLX/3/HGHTA_is_NEW
        VWIND[GR.iijjs] = VWIND[GR.iijjs] * HGHTA_js_OLD/HGHTA_js_NEW \
                            + dVFLX/3/HGHTA_js_NEW
        UWIND = exchange_BC(GR, UWIND)
        VWIND = exchange_BC(GR, VWIND)
        TRACER[GR.iijj] = TRACER[GR.iijj] + dTRACER/3
        TRACER = exchange_BC(GR, TRACER)
        
        #UWIND, VWIND = diagnose_fields_jacobson(GR, HGHT_INT, TRACER_INT,
        #                UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

        ########## level 3
        dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt = tendencies_jacobson(GR, HGHT_INT, TRACER_INT,
                                                    HSURF, UWIND_INT, VWIND_INT, WIND,
                                                    UFLX, VFLX)

        dUFLX = GR.dt*dUFLXdt
        dVFLX = GR.dt*dVFLXdt
        dHGHT = GR.dt*dHGHTdt
        dTRACER = GR.dt*dTRACERdt

        HGHT_INT[GR.iijj] = HGHT_START[GR.iijj] + dHGHT
        HGHT_INT = exchange_BC(GR, HGHT_INT)
        HGHTA_is_NEW, HGHTA_js_NEW = interp_HGHTA(GR, HGHT_INT)
        UWIND_INT[GR.iisjj] = UWIND_START[GR.iisjj] * HGHTA_is_START/HGHTA_is_NEW \
                            + dUFLX/HGHTA_is_NEW
        VWIND_INT[GR.iijjs] = VWIND_START[GR.iijjs] * HGHTA_js_START/HGHTA_js_NEW \
                            + dVFLX/HGHTA_js_NEW
        UWIND_INT = exchange_BC(GR, UWIND_INT)
        VWIND_INT = exchange_BC(GR, VWIND_INT)
        TRACER_INT[GR.iijj] = TRACER_START[GR.iijj] + dTRACER
        TRACER_INT = exchange_BC(GR, TRACER_INT)

        HGHT_OLD = copy.deepcopy(HGHT)
        HGHTA_is_OLD, HGHTA_js_OLD = interp_HGHTA(GR, HGHT_OLD)
        HGHT[GR.iijj] = HGHT[GR.iijj] + dHGHT/3
        HGHT = exchange_BC(GR, HGHT)
        HGHTA_is_NEW, HGHTA_js_NEW = interp_HGHTA(GR, HGHT)
        UWIND[GR.iisjj] = UWIND[GR.iisjj] * HGHTA_is_OLD/HGHTA_is_NEW \
                        + dUFLX/3/HGHTA_is_NEW
        VWIND[GR.iijjs] = VWIND[GR.iijjs] * HGHTA_js_OLD/HGHTA_js_NEW \
                            + dVFLX/3/HGHTA_js_NEW
        UWIND = exchange_BC(GR, UWIND)
        VWIND = exchange_BC(GR, VWIND)
        TRACER[GR.iijj] = TRACER[GR.iijj] + dTRACER/3
        TRACER = exchange_BC(GR, TRACER)
        
        #UWIND, VWIND = diagnose_fields_jacobson(GR, HGHT_INT, TRACER_INT,
        #                UWIND, VWIND, UFLXMP_INT, VFLXMP_INT, HSURF)

        ########## level 4
        dHGHTdt, dUFLXdt, dVFLXdt, dTRACERdt = tendencies_jacobson(GR, HGHT_INT, TRACER_INT,
                                                    HSURF, UWIND_INT, VWIND_INT, WIND,
                                                    UFLX, VFLX)

        dUFLX = GR.dt*dUFLXdt
        dVFLX = GR.dt*dVFLXdt
        dHGHT = GR.dt*dHGHTdt
        dTRACER = GR.dt*dTRACERdt

        HGHT_OLD = copy.deepcopy(HGHT)
        HGHTA_is_OLD, HGHTA_js_OLD = interp_HGHTA(GR, HGHT_OLD)
        HGHT[GR.iijj] = HGHT[GR.iijj] + dHGHT/6
        HGHT = exchange_BC(GR, HGHT)
        HGHTA_is_NEW, HGHTA_js_NEW = interp_HGHTA(GR, HGHT)
        UWIND[GR.iisjj] = UWIND[GR.iisjj] * HGHTA_is_OLD/HGHTA_is_NEW \
                        + dUFLX/6/HGHTA_is_NEW
        VWIND[GR.iijjs] = VWIND[GR.iijjs] * HGHTA_js_OLD/HGHTA_js_NEW \
                            + dVFLX/6/HGHTA_js_NEW
        UWIND = exchange_BC(GR, UWIND)
        VWIND = exchange_BC(GR, VWIND)
        TRACER[GR.iijj] = TRACER[GR.iijj] + dTRACER/6
        TRACER = exchange_BC(GR, TRACER)

        #UWIND, VWIND = diagnose_fields_jacobson(GR, HGHT, TRACER,
        #                UWIND, VWIND, UFLXMP, VFLXMP, HSURF)

    return(HGHT, TRACER,
            UWIND, VWIND,
            UFLX, VFLX, UFLXMP, VFLXMP,
            UUFLX, UVFLX, VUFLX, VVFLX,
            HSURF)
