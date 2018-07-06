# GRID PARAMS
nb = 1

lat0_deg = -80
lat1_deg = 80
lon0_deg = 0
lon1_deg = 360

dlat_deg = 5
dlon_deg = 5

i_curved_earth = 0

# SIMULATION
i_sim_n_days = 5
i_out_nth_hour = 24
CFL = 0.5

# SPATIAL DISCRETIZATION

# TIME DISCRETIZATION
i_time_stepping = 'EULER_FORWARD'
i_time_stepping = 'MATSUNO'
#i_time_stepping = 'RK4'

# INITIAL CONDITIONS
u0 = 10
h0 = 10000
hpert = 1000
h_random_pert = 0

# PSEUDO RADIATION
i_pseudo_radiation = 0
inpRate = 0.002
outRate = 0.001*1E-3

