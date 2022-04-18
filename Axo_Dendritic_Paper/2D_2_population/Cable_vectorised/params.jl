# HOPF Pattern

#PARAMETERS
#PARAMETERS
Wc = 1
 W_EE =Wc*2.5
 W_EI = Wc*2.2
 W_IE = Wc*3.0
 W_II = Wc*2.0
 tau_E = 1  #milliseconds #x10
 tau_I = 1
#finds the grid number for origin on the cable
  v = 1
  v_EE = v #mm/s  #connectivity speed
  v_EI = v #mm/s  #connectivity speed
  v_IE = v #mm/s  #connectivity speed
  v_II = v  #mm/s  #connectivity speed
  sigma_I = 1
  sigma_E = 1
  σ_EE = sigma_E*1.8;   #x10^-1 mm #
  σ_EI = sigma_I*0.8;   #x10^-1 mm #
  σ_IE = sigma_E*1.8;   #x10^-1 mm #
  σ_II = sigma_I*0.8;   #x10^-1 mm #
  h = 0; #threshold
  β_EE = 1 / v_EE
  β_EI = 1 / v_EI
  β_IE = 1 / v_IE
  β_II = 1 / v_II
  τ_EE = 1   #time constant for cable equation
  D_E = (0.5)#diffusion coefficient cable equation
  D_I = (0.5) #diffusion coefficient cable equation

  κ_EE = 0.0
  κ_EI = 0.0
  κ_IE = 0.0
  κ_II = 0.0
  α_EE = 0.5
  α_EI = 0.5
  α_IE = 0.7
  α_II = 0.2
 #milliseconds #x10^-3 seconds
  type = "tanh" #type of firing rate function (0 for heaviside, 1 for sigmoid)
  beta = 5#Sigmoid steepness (if using...
  d = 0.1
  d_EE = d #position on cable (multiples of dx)
  d_EI = d
  d_IE = d
  d_II = d
  g_EXT = 300
  g_EEc = g_EXT
  g_EIc = g_EXT
  g_IEc = g_EXT
  g_IIc = g_EXT
  V_plus = 0.3
  V_neg = -V_plus


#Turing Hopf

#PARAMETERS
Wc = 1
 W_EE =Wc*1.0
 W_EI = Wc*1.2
 W_IE = Wc*2.0
 W_II = Wc*1.7
 tau_E = 1  #milliseconds #x10
 tau_I = 5
#finds the grid number for origin on the cable
  v = 0.2
  v_EE = v #mm/s  #connectivity speed
  v_EI = v #mm/s  #connectivity speed
  v_IE = v #mm/s  #connectivity speed
  v_II = v  #mm/s  #connectivity speed
  sigma_I = 1
  sigma_E = 1
  σ_EE = sigma_E*1.4;   #x10^-1 mm #
  σ_EI = sigma_I*0.2;   #x10^-1 mm #
  σ_IE = sigma_E*1.5;   #x10^-1 mm #
  σ_II = sigma_I*0.2;   #x10^-1 mm #
  h = 0; #threshold
  β_EE = 1 / v_EE
  β_EI = 1 / v_EI
  β_IE = 1 / v_IE
  β_II = 1 / v_II
  τ_EE = 1   #time constant for cable equation
  D_E = (0.1)#diffusion coefficient cable equation
  D_I = (0.1) #diffusion coefficient cable equation

  κ_EE = 0.0
  κ_EI = 0.0
  κ_IE = 0.0
  κ_II = 0.0
  α_EE = 0.5
  α_EI = 0.5
  α_IE = 0.7
  α_II = 0.2
 #milliseconds #x10^-3 seconds
  type = "tanh" #type of firing rate function (0 for heaviside, 1 for sigmoid)
  beta = 5#Sigmoid steepness (if using...
  d = 0.1
  d_EE = d #position on cable (multiples of dx)
  d_EI = d
  d_IE = d
  d_II = d
  g_EXT = 300
  g_EEc = g_EXT
  g_EIc = g_EXT
  g_IEc = g_EXT
  g_IIc = g_EXT
  V_plus = 0.1
  V_neg = -V_plus
