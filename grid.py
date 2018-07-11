import numpy as np
from namelist import *
from constants import con_rE, con_omega
from boundaries import exchange_BC_rigid_y, exchange_BC

class Grid:

    def __init__(self):

        # GRID DEFINITION IN DEGREES 
        self.lon0_deg = lon0_deg
        self.lon1_deg = lon1_deg
        self.lat0_deg = lat0_deg
        self.lat1_deg = lat1_deg
        self.dlon_deg = dlon_deg
        self.dlat_deg = dlat_deg

        # GRID DEFINITION IN RADIANS
        self.lon0_rad = self.lon0_deg/180*np.pi
        self.lon1_rad = self.lon1_deg/180*np.pi
        self.lat0_rad = self.lat0_deg/180*np.pi
        self.lat1_rad = self.lat1_deg/180*np.pi
        self.dlon_rad = self.dlon_deg/180*np.pi
        self.dlat_rad = self.dlat_deg/180*np.pi

        # NUMBER OF GRID POINTS IN EACH DIMENSION
        self.nx = int((self.lon1_deg - self.lon0_deg)/self.dlon_deg)
        self.nxs = self.nx + 1
        self.ny = int((self.lat1_deg - self.lat0_deg)/self.dlat_deg)
        self.nys = self.ny + 1
        self.nb = nb

        # INDEX ARRAYS
        self.ii = np.arange((self.nb),(self.nx+self.nb)) 
        self.jj = np.arange((self.nb),(self.ny+self.nb)) 
        self.iis = np.arange((self.nb),(self.nxs+self.nb)) 
        self.jjs = np.arange((self.nb),(self.nys+self.nb)) 

        self.iijj           = np.ix_(self.ii   ,self.jj  )
        self.iijj_im1       = np.ix_(self.ii-1 ,self.jj  )
        self.iijj_im1_jp1   = np.ix_(self.ii-1 ,self.jj+1)
        self.iijj_ip1       = np.ix_(self.ii+1 ,self.jj  )
        self.iijj_ip1_jm1   = np.ix_(self.ii+1 ,self.jj-1)
        self.iijj_ip1_jp1   = np.ix_(self.ii+1 ,self.jj+1)
        self.iijj_jm1       = np.ix_(self.ii   ,self.jj-1)
        self.iijj_jp1       = np.ix_(self.ii   ,self.jj+1)

        self.iisjj          = np.ix_(self.iis  ,self.jj  )
        self.iisjj_jm1      = np.ix_(self.iis  ,self.jj-1)
        self.iisjj_jp1      = np.ix_(self.iis  ,self.jj+1)
        self.iisjj_im1      = np.ix_(self.iis-1,self.jj  )
        self.iisjj_im1_jm1  = np.ix_(self.iis-1,self.jj-1)
        self.iisjj_im1_jp1  = np.ix_(self.iis-1,self.jj+1)
        self.iisjj_ip1      = np.ix_(self.iis+1,self.jj  )
        self.iisjj_ip1_jm1  = np.ix_(self.iis+1,self.jj-1)
        self.iisjj_ip1_jp1  = np.ix_(self.iis+1,self.jj+1)

        self.iijjs          = np.ix_(self.ii  ,self.jjs  )
        self.iijjs_im1      = np.ix_(self.ii-1,self.jjs  )
        self.iijjs_ip1      = np.ix_(self.ii+1,self.jjs  )
        self.iijjs_jm1      = np.ix_(self.ii  ,self.jjs-1)
        self.iijjs_im1_jm1  = np.ix_(self.ii-1,self.jjs-1)
        self.iijjs_im1_jp1  = np.ix_(self.ii-1,self.jjs+1)
        self.iijjs_ip1_jm1  = np.ix_(self.ii+1,self.jjs-1)
        self.iijjs_ip1_jp1  = np.ix_(self.ii+1,self.jjs+1)
        self.iijjs_jp1      = np.ix_(self.ii  ,self.jjs+1)

        self.iisjjs         = np.ix_(self.iis  ,self.jjs  )
        self.iisjjs_im1     = np.ix_(self.iis-1,self.jjs  )
        self.iisjjs_im1_jm1 = np.ix_(self.iis-1,self.jjs-1)
        self.iisjjs_im1_jp1 = np.ix_(self.iis-1,self.jjs+1)
        self.iisjjs_ip1     = np.ix_(self.iis+1,self.jjs  )
        self.iisjjs_ip1_jm1 = np.ix_(self.iis+1,self.jjs-1)
        self.iisjjs_jm1     = np.ix_(self.iis  ,self.jjs-1)
        self.iisjjs_jp1     = np.ix_(self.iis  ,self.jjs+1)

        # 2D MATRIX OF LONGITUDES AND LATITUDES IN DEGREES
        self.lon_deg = np.full( (self.nx+2*self.nb,self.ny+2*self.nb), np.nan)
        self.lat_deg = np.full( (self.nx+2*self.nb,self.ny+2*self.nb), np.nan)
        self.lonis_deg = np.full( (self.nxs+2*self.nb,self.ny+2*self.nb), np.nan)
        self.latis_deg = np.full( (self.nxs+2*self.nb,self.ny+2*self.nb), np.nan)
        self.lonjs_deg = np.full( (self.nx+2*self.nb,self.nys+2*self.nb), np.nan)
        self.latjs_deg = np.full( (self.nx+2*self.nb,self.nys+2*self.nb), np.nan)

        for j in range(self.nb, self.ny+self.nb):
            self.lon_deg[self.ii,j] = self.lon0_deg + \
                                    (self.ii-self.nb+0.5)*self.dlon_deg
            self.lonis_deg[self.iis,j] = self.lon0_deg + \
                                    (self.iis-self.nb)*self.dlon_deg
        for j_s in range(self.nb, self.nys+self.nb):
            self.lonjs_deg[self.ii,j_s] = self.lon0_deg + \
                                    (self.ii-self.nb+0.5)*self.dlon_deg

        for i in range(self.nb, self.nx+self.nb):
            self.lat_deg[i,self.jj] = self.lat0_deg + \
                                    (self.jj-self.nb+0.5)*self.dlat_deg
            self.latjs_deg[i,self.jjs] = self.lat0_deg + \
                                    (self.jjs-self.nb)*self.dlat_deg
        for i_s in range(self.nb, self.nxs+self.nb):
            self.latis_deg[i_s,self.jj] = self.lat0_deg + \
                                    (self.jj-self.nb+0.5)*self.dlat_deg

        # 2D MATRIX OF LONGITUDES AND LATITUDES IN RADIANS
        self.lon_rad = self.lon_deg/180*np.pi
        self.lat_rad = self.lat_deg/180*np.pi
        self.lonis_rad = self.lonis_deg/180*np.pi
        self.latis_rad = self.latis_deg/180*np.pi
        self.lonjs_rad = self.lonjs_deg/180*np.pi
        self.latjs_rad = self.latjs_deg/180*np.pi

        # 2D MATRIX OF GRID SPACING IN METERS
        self.dx = np.full( (self.nx+2*self.nb,self.ny+2*self.nb), np.nan)
        self.dxjs = np.full( (self.nx+2*self.nb,self.nys+2*self.nb), np.nan)
        self.dxis = np.full( (self.nxs+2*self.nb,self.ny+2*self.nb), np.nan)

        self.dx[self.iijj] = np.cos( self.lat_rad[self.iijj] )*self.dlon_rad*con_rE 
        self.dxjs[self.iijjs] = np.cos( self.latjs_rad[self.iijjs] )*self.dlon_rad*con_rE 
        self.dx = exchange_BC_rigid_y(self, self.dx)
        self.dxjs = exchange_BC_rigid_y(self, self.dxjs)
        #self.dxis[self.iisjj] = np.cos( self.latis_rad[self.iisjj] )*self.dlon_rad*con_rE 
        self.dy = self.dlat_rad*con_rE

        if not i_curved_earth:
            maxdx = np.max(self.dx[self.iijj])
            self.dx[self.iijj] = maxdx
            self.dxjs[self.iijjs] = maxdx
            self.dxis[self.iisjj] = maxdx

        self.A = np.full( (self.nx+2*self.nb,self.ny+2*self.nb), np.nan)
        for i in self.ii:
            for j in self.jj:
                #lon0 = self.lonis_rad[i,j]
                #lon1 = self.lonis_rad[i+1,j]
                #lat0 = self.latjs_rad[i,j]
                #lat1 = self.latjs_rad[i,j+1]
                #self.A[i,j] = lat_lon_recangle_area(lon0,lon1,lat0,lat1, i_curved_earth)
                self.A[i,j] = lat_lon_recangle_area(self.lat_rad[i,j],
                                    self.dlon_rad, self.dlat_rad, i_curved_earth)
        self.A = exchange_BC(self, self.A)

        print('fraction of earth covered: ' + str(np.round(np.sum(self.A[self.iijj])/(4*np.pi*con_rE**2),2)))
        print('fraction of cylinder covered: ' + str(np.round(np.sum(self.A[self.iijj])/(2*np.pi**2*con_rE**2),2)))

        # CORIOLIS FORCE
        self.corf = np.full( (self.nx+2*self.nb,self.ny+2*self.nb), np.nan)
        self.corf_is = np.full( (self.nxs+2*self.nb,self.ny+2*self.nb), np.nan)
        self.corf_js = np.full( (self.nx+2*self.nb,self.nys+2*self.nb), np.nan)
        self.corf[self.iijj] = 2*con_omega*np.sin(self.lat_rad[self.iijj])
        self.corf_is[self.iisjj] = 2*con_omega*np.sin(self.latis_rad[self.iisjj])
        self.corf_js[self.iijjs] = 2*con_omega*np.sin(self.latjs_rad[self.iijjs])

        # TIME STEP
        mindx = np.nanmin(self.dx)
        self.CFL = CFL
        self.i_out_nth_hour = i_out_nth_hour
        self.i_sim_n_days = i_sim_n_days
        self.dt = int(self.CFL*mindx/340)
        while i_out_nth_hour*3600 % self.dt > 0:
            self.dt -= 1
        self.nts = i_sim_n_days*3600*24/self.dt
        self.ts = 0
        self.i_out_nth_ts = int(self.i_out_nth_hour*3600 / self.dt)





#def lat_lon_recangle_area(lon0,lon1,lat0,lat1, i_curved_earth):
#    if i_curved_earth:
#        A = 2*np.pi * np.abs(np.sin(lat0) - np.sin(lat1)) * \
#                np.abs(lon0 - lon1)/(2*np.pi) * con_rE**2
#    else:
#        A = np.abs(lat0 - lat1) * np.abs(lon0 - lon1) * con_rE**2
#    return(A)
def lat_lon_recangle_area(lat,dlon,dlat, i_curved_earth):
    if i_curved_earth:
        A = np.cos(lat) * \
                dlon * dlat * con_rE**2
    else:
        A = dlon * dlat * con_rE**2
    return(A)
