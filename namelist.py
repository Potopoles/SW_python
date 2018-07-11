# GRID PARAMS
nb = 1

lat0_deg = -80
lat1_deg = 80
lon0_deg = 0
lon1_deg = 360

dlat_deg = 5
dlon_deg = 5

i_curved_earth = 1

# SIMULATION
i_sim_n_days = 0.5
i_out_nth_hour = 2
CFL = 0.5

# SPATIAL DISCRETIZATION
i_spatial_discretization = 'UPWIND'
i_spatial_discretization = 'JACOBSON'

# TIME DISCRETIZATION
i_time_stepping = 'EULER_FORWARD'
i_time_stepping = 'MATSUNO'
i_time_stepping = 'RK4'

# INITIAL CONDITIONS
u0 = 0
v0 = 0
h0 = 10000
hpert = 1000
upert = 0
vpert = 0
h_random_pert = 0

# PSEUDO RADIATION
i_pseudo_radiation = 0
inpRate = 0.002
outRate = 0.001*1E-3

