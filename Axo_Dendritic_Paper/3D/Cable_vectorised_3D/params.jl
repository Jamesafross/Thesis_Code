


max_T = 20
#mesh paramaters
min_R1 = 0;
max_R1 = 25;
min_R1 = 0;
max_R2 = 25;
min_X = -0.05
max_X = 0.05
dt = 0.1
dx = 0.01;
dr1 = 0.1
dr2 = 0.1
tau = 0.1
#finds the grid number for origin on the cable
#PARAMETERS
v = 18 #connectivity speed
sigma = 1.0; #
h = 0.15; #threshold
α = 1 / sigma;
β = 1 / v
τ = 1  #time constant for cable equation
D = (0.1)^2#diffusion coefficient cable equation
k = 0.02#diffusion coefficient for x spatial interaction for wave PDE
alfa = 1.0
W = 1.0
type = 1 #type of firing rate function (0 for heaviside, 1 for sigmoid)
beta = 100 #sigmoid steepness (if using...)
pos = 0.02 #position on cable (multiples of dx)

V_plus = 70
