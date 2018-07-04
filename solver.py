import numpy as np

from grid import Grid
from fields import initialize_fields
from IO import output_to_NC
from boundaries import exchange_BC_fluxes_mp, exchange_BC_height, \
                        exchange_BC_winds

from continuity import diag_mass_flux


GR = Grid()

HGHT, UWIND, VWIND, WIND, UFLX, VFLX, UFLX_mp, VFLX_mp, \
UUFLX, VVFLX = initialize_fields(GR)

outCounter = 0
WIND[GR.iijj] = np.sqrt( ((UWIND[GR.iijj] + UWIND[GR.iijj_ip1])/2)**2 + \
                ((VWIND[GR.iijj] + VWIND[GR.iijj_jp1])/2)**2 )
output_to_NC(GR, outCounter, HGHT, UWIND, VWIND, WIND)

while GR.ts < GR.nts:
    GR.ts += 1
    GR.sim_time_sec = GR.ts*GR.dt

    vmax = max(np.max(UWIND[GR.iisjj]), np.max(VWIND[GR.iijjs]))
    sum_hght = np.sum(HGHT[GR.iijj]*GR.A[GR.iijj])/np.sum(GR.A[GR.iijj])
    print('#### ' + str(GR.ts) + '  ' + str(np.round(GR.sim_time_sec/3600/24,2)) + \
            ' days  vmax: ' + str(np.round(vmax,1)) + ' m/s  hght: ' + \
            str(np.round(sum_hght,3)) + ' m')

    # PROGNOSE HGHT
    UFLX[GR.iisjj] = np.maximum(UWIND[GR.iisjj],0)*HGHT[GR.iisjj_im1] + \
            np.minimum(UWIND[GR.iisjj],0)*HGHT[GR.iisjj]

    VFLX[GR.iijjs] = np.maximum(VWIND[GR.iijjs],0)*HGHT[GR.iijjs_jm1] + \
            np.minimum(VWIND[GR.iijjs],0)*HGHT[GR.iijjs]

    dHGHTdt = - ( UFLX[GR.iijj_ip1] - UFLX[GR.iijj] ) / GR.dx[GR.iijj] \
            - ( VFLX[GR.iijj_jp1] - VFLX[GR.iijj] ) / GR.dy

    print(np.sum(dHGHTdt))
    #print(dHGHTdt)


    # PROGNOSE WIND
    UFLX_mp[GR.iijj] = HGHT[GR.iijj]*(UWIND[GR.iijj] + UWIND[GR.iijj_ip1])/2
    VFLX_mp[GR.iijj] = HGHT[GR.iijj]*(VWIND[GR.iijj] + VWIND[GR.iijj_jp1])/2
    UFLX_mp, VFLX_mp = exchange_BC_fluxes_mp(GR, UFLX_mp, VFLX_mp)

    UUFLX[GR.iisjj] = np.maximum(UWIND[GR.iisjj],0)*UFLX_mp[GR.iisjj_im1] + \
            np.minimum(UWIND[GR.iisjj],0)*UFLX_mp[GR.iisjj]

    VVFLX[GR.iijjs] = np.maximum(VWIND[GR.iijjs],0)*VFLX_mp[GR.iijjs_jm1] + \
            np.minimum(VWIND[GR.iijjs],0)*VFLX_mp[GR.iijjs]

    dUFLX_mpdt = - ( UUFLX[GR.iijj_ip1] - UUFLX[GR.iijj] ) / GR.dx[GR.iijj] + \
                - HGHT[GR.iijj]*( HGHT[GR.iijj_ip1] - HGHT[GR.iijj_im1] ) / (2*GR.dx[GR.iijj])

    dVFLX_mpdt = - ( VVFLX[GR.iijj_jp1] - VVFLX[GR.iijj] ) / GR.dy + \
                - HGHT[GR.iijj]*( HGHT[GR.iijj_jp1] - HGHT[GR.iijj_jm1] ) / (2*GR.dy)

    # TIME STEPPING
    UFLX_mp[GR.iijj] = UFLX_mp[GR.iijj] + GR.dt*dUFLX_mpdt
    VFLX_mp[GR.iijj] = VFLX_mp[GR.iijj] + GR.dt*dVFLX_mpdt
    UFLX_mp, VFLX_mp = exchange_BC_fluxes_mp(GR, UFLX_mp, VFLX_mp)
    HGHT[GR.iijj] = HGHT[GR.iijj] + GR.dt*dHGHTdt
    HGHT = exchange_BC_height(GR, HGHT)

    UWIND[GR.iisjj] = ( UFLX_mp[GR.iisjj_im1] + UFLX_mp[GR.iisjj] ) \
                    / ( HGHT[GR.iisjj_im1] + HGHT[GR.iisjj] )
    #UWIND[GR.iisjj] = 0
    VWIND[GR.iijjs] = ( VFLX_mp[GR.iijjs_jm1] + VFLX_mp[GR.iijjs] ) \
                    / ( HGHT[GR.iijjs_jm1] + HGHT[GR.iijjs] )
    #print(VWIND)
    UWIND,VWIND = exchange_BC_winds(GR, UWIND, VWIND)
    #print(VWIND)
    #quit()
    #VWIND[GR.iijjs] = 0


    if GR.ts % GR.i_out_nth_ts == 0:
        outCounter = outCounter + 1
        WIND[GR.iijj] = np.sqrt( ((UWIND[GR.iijj] + UWIND[GR.iijj_ip1])/2)**2 + \
                        ((VWIND[GR.iijj] + VWIND[GR.iijj_jp1])/2)**2 )
        output_to_NC(GR, outCounter, HGHT, UWIND, VWIND, WIND)

