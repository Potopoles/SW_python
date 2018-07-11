import numpy as np

from grid import Grid
from fields import initialize_fields
from IO import output_to_NC
from namelist import i_time_stepping, i_spatial_discretization
if i_time_stepping == 'EULER_FORWARD':
    from time_integration import euler_forward as time_stepper
elif i_time_stepping == 'MATSUNO':
    from time_integration import matsuno as time_stepper
elif i_time_stepping == 'RK4':
    from time_integration import RK4 as time_stepper

GR = Grid()

HGHT, HTOP, UWIND, VWIND, WIND, \
UFLX, VFLX, UFLXMP, VFLXMP, \
UUFLX, VUFLX, UVFLX, VVFLX, \
HSURF, TRACER = initialize_fields(GR)

outCounter = 0
WIND[GR.iijj] = np.sqrt( ((UWIND[GR.iijj] + UWIND[GR.iijj_ip1])/2)**2 + \
                ((VWIND[GR.iijj] + VWIND[GR.iijj_jp1])/2)**2 )
mean_ekin = np.sum( 0.5*1*WIND[GR.iijj]**2*HGHT[GR.iijj]*GR.A[GR.iijj] ) / \
                np.sum( HGHT[GR.iijj]*GR.A[GR.iijj] )
output_to_NC(GR, outCounter, HGHT, HTOP, UWIND, VWIND, WIND,
            HSURF, TRACER,
            mean_ekin)

while GR.ts < GR.nts:
    GR.ts += 1
    GR.sim_time_sec = GR.ts*GR.dt

    vmax = max(np.max(UWIND[GR.iisjj]), np.max(VWIND[GR.iijjs]))
    mean_hght = np.sum(HGHT[GR.iijj]*GR.A[GR.iijj])/np.sum(GR.A[GR.iijj])
    mean_tracer = np.sum(TRACER[GR.iijj]*GR.A[GR.iijj]*HGHT[GR.iijj])/ \
            np.sum(GR.A[GR.iijj]*HGHT[GR.iijj])
    WIND[GR.iijj] = np.sqrt( ((UWIND[GR.iijj] + UWIND[GR.iijj_ip1])/2)**2 + \
                    ((VWIND[GR.iijj] + VWIND[GR.iijj_jp1])/2)**2 )
    mean_ekin = np.sum( 0.5*1*WIND[GR.iijj]**2*HGHT[GR.iijj]*GR.A[GR.iijj] ) / \
            np.sum( HGHT[GR.iijj]*GR.A[GR.iijj] )
    print('#### ' + str(GR.ts) + '  ' + str(np.round(GR.sim_time_sec/3600/24,2)) + \
            '  days  vmax: ' + str(np.round(vmax,1)) + '  m/s  hght: ' + \
            str(np.round(mean_hght,2)) + '  m ekin: ' + \
            str(np.round(mean_ekin,3)) + '  tracer: ' + \
            str(np.round(mean_tracer,7)))

    HGHT, TRACER, \
    UWIND, VWIND, \
    UFLX, VFLX, UFLXMP, VFLXMP, \
    UUFLX, UVFLX, VUFLX, VVFLX, \
    HSURF = time_stepper(GR, HGHT, TRACER,
                    UWIND, VWIND, WIND,
                    UFLX, VFLX, UFLXMP, VFLXMP,
                    UUFLX, UVFLX, VUFLX, VVFLX,
                    HSURF, i_spatial_discretization)

    # OUTPUT DIAGNOSTICS
    HTOP[GR.iijj] = HSURF[GR.iijj] + HGHT[GR.iijj] 


    if GR.ts % GR.i_out_nth_ts == 0:
        outCounter = outCounter + 1
        print('write fields')
        WIND[GR.iijj] = np.sqrt( ((UWIND[GR.iijj] + UWIND[GR.iijj_ip1])/2)**2 + \
                        ((VWIND[GR.iijj] + VWIND[GR.iijj_jp1])/2)**2 )
        mean_ekin = np.sum( 0.5*1*WIND[GR.iijj]**2*HGHT[GR.iijj]*GR.A[GR.iijj] ) / \
                        np.sum( HGHT[GR.iijj]*GR.A[GR.iijj] )
        output_to_NC(GR, outCounter, HGHT, HTOP, UWIND, VWIND, WIND,
                    HSURF, TRACER,
                    mean_ekin)

