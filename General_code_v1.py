# This code solves:
# a. 12 coupled differential equations viz. 5 for number concentration of the particles; 
#    5 for the volume concentration of the particles; 2 for ions
# b. Can take four different temperature inputs viz. 300, 600, 1200, 2400
# c. beta (for charging) has been evaluated at these 4 temperatures using Fuchs Model (Makela paper): 
#    The code for the Fuchs Model is written by Yang Wang in Matlab


import os
import math as ma
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode

ke = 8.99e9                                     # in SI units Nm2C-2  1/(4*pi*eps0) 
kb = 1.38e-23                                   # in SI units
Na = 6.022e23                                   # Avagadro's Number
e  = 1.602e-19                                  # Charge on one electron

dp_zero_0 = 0.4e-9
T  = 1200 # can take values 300, 600, 1200, 2400                                  # in Kelvin
CHAR = 1
COAG = 1
DISP = 1
QYN  = 0  # 0 = no change in ion concentration, 1 = solve differential equation for ion balance to see the changing ion concentration


# initial conditions
ion = 16
par =  18
n_m_0      = 1.0*10**(ion)
n_p_0      = 1.0*10**(ion)

N_0_0  = 1.0*10**(par)

if CHAR == 1 and COAG ==1:
    file_add = '_coag_char'
elif CHAR == 1 and COAG == 0:
    file_add = '_char'
elif CHAR == 0 and COAG == 1:
    file_add = '_coag'
else:
    pass

if DISP == 1:
    file_add = file_add + '_disp'
else:
    pass

if QYN == 1:
    file_add = file_add + '_qyn1' 
else:
    file_add = file_add + '_qyn0' 

qm = 4 # maximum charge
N = 30 # data points in beta_char

T = float(T)
P  = 101325.0                                   # in N/m2
M = 80.0*1e-3/Na

Mg = 28.97                                      # Molecular Weight for air
mu = 1.716e-5 * (T/273)**(2.0/3)                # viscosity of air   
lam = mu/P * ma.sqrt(np.pi*kb*T/(2.0*Mg/Na))    # mean free path in air(m)

q_neg = 0.0                                  # rate of production of negative ion
q_pos = 0.0                                  # rate of production of positive ion
alpha = 1e-13


# for integration
t_final = 5.0 # time in seconds
t_start = 1e-5
delta_t = 1.1


N_m4_0 = 0.0
N_m3_0 = 0.0
N_m2_0 = 0.0
N_m1_0 = 0.0

N_p1_0 = 0.0
N_p2_0 = 0.0
N_p3_0 = 0.0
N_p4_0 = 0.0


v_zero_0     = ma.pi/6 * (dp_zero_0)**3


V_m4_0 = 0.0
V_m3_0 = 0.0
V_m2_0 = 0.0
V_m1_0 = 0.0
V_0_0  = v_zero_0*N_0_0
V_p1_0 = 0.0
V_p2_0 = 0.0
V_p3_0 = 0.0
V_p4_0 = 0.0


Initial = [ N_m4_0, N_m3_0, N_m2_0, N_m1_0, N_0_0, N_p1_0, N_p2_0, N_p3_0, N_p4_0, \
                        V_m4_0, V_m3_0, V_m2_0, V_m1_0, V_0_0, V_p1_0, V_p2_0, V_p3_0, V_p4_0, n_m_0, n_p_0]
# raw_input('ii')

# ############################################## FOR CHARGING ##############################################################
def beta_pos(dp):
    x_data = x_data_pos
    y_data = y_data_pos
    
    beta = np.zeros((2*qm+1))
    if dp < x_data[0]:
        pass
    else:
        for i in range(0,N):
            if dp >= x_data[i] and dp < x_data[i+1]:
                beta = ((x_data[i+1] - dp) * y_data[i] + (dp - x_data[i])  * y_data[i+1])/(x_data[i+1] - x_data[i])

    return beta



def beta_neg(dp):
    x_data = x_data_neg
    y_data = y_data_neg
    
    beta = np.zeros((2*qm+1))
    if dp < x_data[0]:
        pass
    else:
        for i in range(0,N):
            if dp >= x_data[i] and dp < x_data[i+1]:
                beta = ((x_data[i+1] - dp) * y_data[i] + (dp - x_data[i])  * y_data[i+1])/(x_data[i+1] - x_data[i])

    return beta




def xydata_pos():
    filename = 'beta_pos_' + str(int(T)) + '.txt'
    file = open(filename, 'r')
    i = 0
    x_data = np.zeros((N))
    y_data = np.zeros((N,2*qm+1))
    for line in file:
        line = line.split("\t")
        x_data[i] = 1e-9*float(line[0])
        for j in range (0,2*qm+1):
            y_data[i][j] = float(line[j+1])
        i = i+1
    
    return x_data, y_data



def xydata_neg():
    filename = 'beta_neg_' + str(int(T)) + '.txt'
    file = open(filename, 'r')
    i = 0
    x_data = np.zeros((N))
    y_data = np.zeros((N,2*qm+1))
    for line in file:
        line = line.split("\t")
        x_data[i] = 1e-9*float(line[0])
        for j in range (0,2*qm+1):
            y_data[i][j] = float(line[j+1])
        i = i+1
    
    return x_data, y_data



############################################## FOR COAGULATION ############################################################
# DIFFUSIVITY
def D(dp):
#     m = pho * v
#     c = ma.sqrt(8*kb*T/(ma.pi*m))
    Kn = 2.0*lam/dp
    D = kb*T/(3.0*ma.pi*mu*dp) * ((5.0 + 4.0*Kn + 6.0*Kn*Kn + 18.0*Kn*Kn*Kn)/(5.0 - Kn + (8.0 + ma.pi)*Kn*Kn))
    return D



# COAGULATION COEFFICIENT USING FUCHS SUTIGUN
def beta(q1,q2,dp1,dp2):
    q1 = float(q1)
    q2 = float(q2)
    c1 = ma.sqrt(8.0*kb*T/(ma.pi*M))
    c2 = ma.sqrt(8.0*kb*T/(ma.pi*M))    
    
    l1 = 8.0*D(dp1)/(ma.pi*c1)
    l2 = 8.0*D(dp2)/(ma.pi*c2)
    
    g1 = 1.0/(3.0*dp1*l1) * ((dp1+l1)**3-(dp1*dp1+l1*l1)**(3.0/2)) - dp1
    g2 = 1.0/(3.0*dp2*l2) * ((dp2+l2)**3-(dp2*dp2+l2*l2)**(3.0/2)) - dp2
    

    coeff1 = 2.0*ma.pi*(dp1 + dp2) * (D(dp1) + D(dp2))
    coeff2 = (dp1 + dp2)/(dp1 + dp2 + 2.0*ma.sqrt(g1*g1 + g2*g2)) + (8.0*(D(dp1)+D(dp2)))/((dp1+dp2)*ma.sqrt(c1*c1+c2*c2))
    
    if dp1 < dp_zero_0 or dp2 < dp_zero_0:
        beta = 0.0
    else:
        if q1 == 0 or q2 == 0:
            W = 1.0
        else:
            y = ke * e* e * q1 * q2/ (kb * T * (dp1 + dp2)/2)
            W = 1.0/y * (ma.exp(y) - 1)
        beta = (1.0/W) * (coeff1/coeff2)
    
    return beta


# COEFFICIENT ACCOUNTING FOR ELECTROSTATIC DISPERSION
def B(dp): 
    A1 = 1.257 
    A2 = 0.400 
    A3 = 0.55
    if dp > dp_zero_0:
        C = 1.0 + 2.0*lam/dp * (A1 + A2 * ma.exp (-A3 * dp/lam))
        B =  C/(3.0 * ma.pi * mu * dp)
    else:
        B = 0.0
    return B


# RIGHT HAND SIDE FOR THE SOLVER
def rhs(t,y):
    ######################################## INITILAIZING ####################################
    # find the important quantities from y
    N = y[0:2*qm + 1]
    V = y[2*qm + 1:2*(2*qm +1)]
    n_neg = y[2*(2*qm +1)]
    n_pos = y[2*(2*qm +1)+1]
    
    # initializing 
    v = np.zeros((2*qm +1)) # v = V/N
    dp= np.zeros((2*qm +1)) # v = pi/6 * dp*dp*dp
    
#     print v
    for i in range(0,2*qm +1):
        if N[i] >= 1e-90:
            v[i] = V[i]/N[i]
            dp[i]= (6.0/np.pi * v[i])**(1.0/3)
    
    
    dN = np.zeros((2*qm +1))
    dN_coag = np.zeros((2*qm +1))
    dN_char = np.zeros((2*qm +1))
    dN_disp = np.zeros((2*qm +1))
    dV = np.zeros((2*qm +1))
    dV_coag = np.zeros((2*qm +1))
    dV_char = np.zeros((2*qm +1))
    dV_disp = np.zeros((2*qm +1))
    sumN1 = np.zeros((2*qm +1))
    sumN2 = np.zeros((2*qm +1))
    
    sumV1 = np.zeros((2*qm +1))
    sumV2 = np.zeros((2*qm +1))
    
    
    ############################################ COAGULATION #######################################
    for k in range(-qm,qm+1):
        
        for i in range(-qm,qm+1):
            
            sumN2[k+qm] = sumN2[k+qm] + beta(i,k,dp[i+qm],dp[k+qm]) * N[i+qm] * N[k+qm]
            sumV2[k+qm] = sumV2[k+qm] + beta(i,k,dp[i+qm],dp[k+qm]) * N[i+qm] * N[k+qm] * v[k+qm]
            
            for j in range(-qm,qm+1):
                
                if i + j == k:
                    sumN1[k+qm] = sumN1[k+qm] + beta(i,j,dp[i+qm],dp[j+qm]) * N[i+qm] * N[j+qm]
                    sumV1[k+qm] = sumV1[k+qm] + beta(i,j,dp[i+qm],dp[j+qm]) * N[i+qm] * N[j+qm] * (v[i+qm] + v[j+qm])
#                     print i,j,k
#                     print sumN1[k+qm], beta(i,j,dp[i+qm],dp[j+qm]), N[i+qm], N[j+qm]
                else:
                    pass
        
        dN_coag[k+qm] = 1.0/2 * sumN1[k+qm] - sumN2[k+qm]
        dV_coag[k+qm] = 1.0/2 * sumV1[k+qm] - sumV2[k+qm]


        
        
    ############################################# CHARGING #########################################
    #### For Number and Volume Concentration:
    
    dN_char[-qm+qm] = beta_neg(dp[-(qm-1)+qm])[-(qm-1)+qm] * n_neg * N[-(qm-1)+qm] \
                    - beta_pos(dp[-qm+qm])[-qm+qm] * n_pos * N[-qm+qm]
    dV_char[-qm+qm] = beta_neg(dp[-(qm-1)+qm])[-(qm-1)+qm] * n_neg * V[-(qm-1)+qm] \
                    - beta_pos(dp[-qm+qm])[-qm+qm] * n_pos * V[-qm+qm]
    for i in range(-(qm-1),(qm)):
        dN_char[i+qm] = beta_neg(dp[i+1+qm])[i+1+qm] * n_neg * N[i+1+qm] \
                     +  beta_pos(dp[i-1+qm])[i-1+qm] * n_pos * N[i-1+qm] \
                     -  beta_neg(dp[i+qm])[i+qm] * n_neg * N[i+qm] \
                     -  beta_pos(dp[i+qm])[i+qm] * n_pos * N[i+qm]
        dV_char[i+qm] = beta_neg(dp[i+1+qm])[i+1+qm] * n_neg * V[i+1+qm] \
                     +  beta_pos(dp[i-1+qm])[i-1+qm] * n_pos * V[i-1+qm] \
                     -  beta_neg(dp[i+qm])[i+qm] * n_neg * V[i+qm] \
                     -  beta_pos(dp[i+qm])[i+qm] * n_pos * V[i+qm]    
                    
    dN_char[+qm+qm] = beta_pos(dp[+(qm-1)+qm])[+(qm-1)+qm] * n_pos * N[+(qm-1)+qm] \
                    - beta_neg(dp[+qm+qm])[+qm+qm] * n_neg * N[+qm+qm] 
    dV_char[+qm+qm] = beta_pos(dp[+(qm-1)+qm])[+(qm-1)+qm] * n_pos * V[+(qm-1)+qm] \
                    - beta_neg(dp[+qm+qm])[+qm+qm] * n_neg * V[+qm+qm] 
        
 
    
    ######################################### DISPERSION ################################################
    temp_sum = 0.0
    for i in range(-qm,qm+1):
        temp_sum = temp_sum + float(i)*N[i+qm]
    average_sum = temp_sum
    for i in range(-qm,qm+1):
        dN_disp[i+qm] = - 4.0 * ma.pi*B(dp[i+qm]) * e*e* N[i+qm] * (i) * average_sum
        dV_disp[i+qm] = - 4.0 * ma.pi*B(dp[i+qm]) * e*e* V[i+qm] * (i) * average_sum 
        
    
    ##### For Ion Number Concentration:     
    temp_sum_neg = 0.0
    temp_sum_pos = 0.0
    for i in range(-qm+1,qm+1):
        temp_sum_neg = temp_sum_neg + beta_neg(dp[i+qm])[i+qm] * N[i+qm]
    for i in range(-qm,qm):
        temp_sum_pos = temp_sum_pos + beta_pos(dp[i+qm])[i+qm] * N[i+qm]
        
    if QYN == 1: 
        dn_neg = q_neg - alpha*n_neg*n_pos -  n_neg * temp_sum_neg
        dn_pos = q_pos - alpha*n_neg*n_pos -  n_pos * temp_sum_pos
    else:
        dn_neg = 0.0
        dn_pos = 0.0
        
        
        
    ############################### COMBINING CHARGING AND COAGULATION ##################################
    if CHAR == 1 and COAG == 1:
        dN = dN_coag + dN_char
        dV = dV_coag + dV_char

    elif CHAR == 1 and COAG == 0:
        dN = dN_char
        dV = dV_char

    elif CHAR == 0 and COAG == 1:
        dN = dN_coag
        dV = dV_coag
        
    else:
        pass
    
    ########################## COMBINING DISPERSION #######################################################
    if DISP == 1:
        dN = dN + dN_disp
        dV = dV + dV_disp
    else:
        pass
    
    dn = np.asarray([dn_neg, dn_pos])
    dy = np.concatenate([dN,dV,dn])
    return dy




################################################ MAIN PROGRAM ################################################
def main():
    
    # Start by specifying the integrator:
    # use ``vode`` with "backward differentiation formula"
    r = ode(rhs).set_integrator('vode', method='bdf', nsteps='1000000000')
 

    # Number of time steps: 1 extra for initial condition
#     num_steps = np.floor((t_final - t_start)/delta_t) + 1
    num_steps = np.floor(ma.log10(t_final/t_start)/ma.log10(delta_t))+2
    

    # Set initial condition(s): for integrating variable and time!
    r.set_initial_value(Initial, t_start)
    
 
    # Additional Python step: create vectors to store variables
    t = np.zeros((num_steps, 1))
    N_m4_t = np.zeros((num_steps, 1))
    N_m3_t = np.zeros((num_steps, 1))
    N_m2_t = np.zeros((num_steps, 1))
    N_m1_t = np.zeros((num_steps, 1))
    N_0_t = np.zeros((num_steps, 1))
    N_p1_t = np.zeros((num_steps, 1))
    N_p2_t = np.zeros((num_steps, 1))
    N_p3_t = np.zeros((num_steps, 1))
    N_p4_t = np.zeros((num_steps, 1))
    

    V_m4_t = np.zeros((num_steps, 1))
    V_m3_t = np.zeros((num_steps, 1))
    V_m2_t = np.zeros((num_steps, 1))
    V_m1_t = np.zeros((num_steps, 1))
    V_0_t = np.zeros((num_steps, 1))
    V_p1_t = np.zeros((num_steps, 1))
    V_p2_t = np.zeros((num_steps, 1))
    V_p3_t = np.zeros((num_steps, 1))
    V_p4_t = np.zeros((num_steps, 1))
    
    n_m_t      = np.zeros((num_steps, 1)) 
    n_p_t      = np.zeros((num_steps, 1)) 

    dp_m4_t = np.zeros((num_steps, 1))
    dp_m3_t = np.zeros((num_steps, 1))
    dp_m2_t = np.zeros((num_steps, 1))
    dp_m1_t = np.zeros((num_steps, 1))
    dp_0_t = np.zeros((num_steps, 1))
    dp_p1_t = np.zeros((num_steps, 1))
    dp_p2_t = np.zeros((num_steps, 1))
    dp_p3_t = np.zeros((num_steps, 1))
    dp_p4_t = np.zeros((num_steps, 1))
    
    N_tot_t = np.zeros((num_steps, 1))
    dp_avg_t= np.zeros((num_steps, 1))
    average_charge = np.zeros((num_steps,1))
    
    t[0] = t_start
    N_m4_t[0] = N_m4_0 
    N_m3_t[0] = N_m3_0 
    N_m2_t[0] = N_m2_0 
    N_m1_t[0] = N_m1_0 
    N_0_t[0]  = N_0_0 
    N_p1_t[0] = N_p1_0 
    N_p2_t[0] = N_p2_0 
    N_p3_t[0] = N_p3_0 
    N_p4_t[0] = N_p4_0 

    V_m4_t[0] = V_m4_0 
    V_m3_t[0] = V_m3_0 
    V_m2_t[0] = V_m2_0 
    V_m1_t[0] = V_m1_0 
    V_0_t[0]  = V_0_0 
    V_p1_t[0] = V_p1_0 
    V_p2_t[0] = V_p2_0 
    V_p3_t[0] = V_p3_0 
    V_p4_t[0] = V_p4_0 
    
    n_m_t[0]  = n_m_0
    n_p_t[0]  = n_p_0

 
    # Integrate the ODE(s) across each delta_t timestep
    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t*delta_t)
 
        # Store the results to plot later
        t[k] = r.t
        
        N_m4_t[k] = r.y[0] 
        N_m3_t[k] = r.y[1] 
        N_m2_t[k] = r.y[2] 
        N_m1_t[k] = r.y[3] 
        N_0_t[k]  = r.y[4] 
        N_p1_t[k] = r.y[5] 
        N_p2_t[k] = r.y[6] 
        N_p3_t[k] = r.y[7] 
        N_p4_t[k] = r.y[8] 

        V_m4_t[k] = r.y[9] 
        V_m3_t[k] = r.y[10] 
        V_m2_t[k] = r.y[11] 
        V_m1_t[k] = r.y[12] 
        V_0_t[k]  = r.y[13] 
        V_p1_t[k] = r.y[14] 
        V_p2_t[k] = r.y[15] 
        V_p3_t[k] = r.y[16] 
        V_p4_t[k] = r.y[17] 

        n_m_t[k]  = r.y[18]
        n_p_t[k]  = r.y[19]
 
        k += 1
    
#     file1 = 'Raw_data_' + str(int(T))
#     with open(file1, 'w') as f:
#         # Print & save the solution.
#         for i in range(0,k):
#             print >> f, t[i], n_p_t[i] # N_minustwo_t[i], N_minusone_t[i], N_zero_t[i], N_plusone_t[i], N_plustwo_t[i], \
#                              # V_minustwo_t[i], V_minusone_t[i], V_zero_t[i], V_plusone_t[i], V_plustwo_t[i], \
#                              # n_neg_t[i],
    
    
    
    #########################   CASE A ##############################
    # creates a temp file to store all the data for CASE A temporarily
    with open('temp', 'w') as f:
        # Print & save the solution.
        for i in range(0,k):

            if N_m4_t[i] < 1e-10:
                dp_m4_t[i] = 0
            else:
                dp_m4_t[i] = (6.0/ma.pi * V_m4_t[i]/N_m4_t[i])**(1.0/3)
            
            if N_m3_t[i] < 1e-10:
                dp_m3_t[i] = 0
            else:
                dp_m3_t[i] = (6.0/ma.pi * V_m3_t[i]/N_m3_t[i])**(1.0/3)
            
            if N_m2_t[i] < 1e-10:
                dp_m2_t[i] = 0
            else:
                dp_m2_t[i] = (6.0/ma.pi * V_m2_t[i]/N_m2_t[i])**(1.0/3)
            
            if N_m1_t[i] < 1e-10:
                dp_m1_t[i] = 0
            else:
                dp_m1_t[i] = (6.0/ma.pi * V_m1_t[i]/N_m1_t[i])**(1.0/3)
            
            if N_0_t[i] < 1e-10:
                dp_0_t[i] = 0
            else:
                dp_0_t[i] = (6.0/ma.pi * V_0_t[i]/N_0_t[i])**(1.0/3)
            
            if N_p1_t[i] < 1e-10:
                dp_p1_t[i] = 0
            else:
                dp_p1_t[i] = (6.0/ma.pi * V_p1_t[i]/N_p1_t[i])**(1.0/3)
                
            if N_p2_t[i] < 1e-10:
                dp_p2_t[i] = 0
            else:
                dp_p2_t[i] = (6.0/ma.pi * V_p2_t[i]/N_p2_t[i])**(1.0/3)
                
            if N_p3_t[i] < 1e-10:
                dp_p3_t[i] = 0
            else:
                dp_p3_t[i] = (6.0/ma.pi * V_p3_t[i]/N_p3_t[i])**(1.0/3)
                
            if N_p4_t[i] < 1e-10:
                dp_p4_t[i] = 0
            else:
                dp_p4_t[i] = (6.0/ma.pi * V_p4_t[i]/N_p4_t[i])**(1.0/3)

            
            
            temp =  dp_m4_t[i]*N_m4_t[i] + dp_m3_t[i]*N_m3_t[i]+ dp_m2_t[i]*N_m2_t[i]\
                    + dp_m1_t[i]*N_m1_t[i] + dp_0_t[i]*N_0_t[i] + dp_p1_t[i]*N_p1_t[i] + dp_p2_t[i]*N_p2_t[i]\
                     + dp_p3_t[i]*N_p3_t[i]  + dp_p4_t[i]*N_p4_t[i]  
            
            
            N_tot_t[i] =  N_m4_t[i] + N_m3_t[i] + N_m2_t[i] + N_m1_t[i] + N_0_t[i] + N_p1_t[i] \
                          + N_p2_t[i] + N_p3_t[i]  + N_p4_t[i]  
            dp_avg_t[i] = temp/N_tot_t[i]
            average_charge[i] = (-4.0*N_m4_t[i] + -3.0*N_m3_t[i] + -2.0*N_m2_t[i] + -1.0*N_m1_t[i] + 0.0*N_0_t[i] + 1.0*N_p1_t[i] \
                          + 2.0*N_p2_t[i] + 3.0*N_p3_t[i]  + 4.0*N_p4_t[i] )/N_tot_t[i]
            print >> f, t[i], N_tot_t[i], dp_avg_t[i], n_p_t[i], n_m_t[i]
    
    # Remove the brackets from temp file [ and ]
    with open('temp', 'r') as my_file:
        text = my_file.read()
        text = text.replace("[", "")
        text = text.replace("]", "")     
   
    # Use proper name for proper sub-case for CASE A
    file2 = 'caseA_'  +   'ion_' + str(int(ion)) +   '_par_' + str(int(par)) + file_add +'.txt'
    with open(file2, 'w') as my_file:
        my_file.write(text)
    
    # Remove the temp file from the folder
    os.remove('temp')


    ##########################   CASE B ##############################
#     with open('temp', 'w') as f:
#         # Print & save the solution.
#         for i in range(0,k):
#             if N_minustwo_t[i] < 1e-10:
#                 dp_minustwo_t[i] = 0
#             else:
#                 dp_minustwo_t[i] = (6.0/ma.pi * V_minustwo_t[i]/N_minustwo_t[i])**(1.0/3)
                
#             if N_minusone_t[i] < 1e-10:
#                 dp_minusone_t[i] = 0
#             else:            
#                 dp_minusone_t[i] = (6.0/ma.pi * V_minusone_t[i]/N_minusone_t[i])**(1.0/3)
            
#             if N_zero_t[i] < 1e-10:
#                 dp_zero_t[i] = 0
#             else:            
#                 dp_zero_t[i]     = (6.0/ma.pi * V_zero_t[i]    /N_zero_t[i]    )**(1.0/3)
            
#             if N_plusone_t[i] < 1e-10:
#                 dp_plusone_t[i] = 0
#             else:
#                 dp_plusone_t[i]  = (6.0/ma.pi * V_plusone_t[i] /N_plusone_t[i] )**(1.0/3)
            
#             if N_plustwo_t[i] < 1e-10:
#                 dp_plustwo_t[i] = 0
#             else:
#                 dp_plustwo_t[i]  = (6.0/ma.pi * V_plustwo_t[i] / N_plustwo_t[i])**(1.0/3)
            
#             temp = dp_minustwo_t[i]*N_minustwo_t[i] + dp_minusone_t[i]*N_minusone_t[i] + dp_zero_t[i]*N_zero_t[i] \
#                    + dp_plusone_t[i]*N_plusone_t[i] + dp_plustwo_t[i]*N_plustwo_t[i]
#             N_tot_t[i] = N_minustwo_t[i] + N_minusone_t[i] + N_zero_t[i] + N_plusone_t[i] +  N_plustwo_t[i]
#             dp_avg_t[i] = temp/N_tot_t[i]
#             print >> f, t[i], N_minusone_t[i]/N_tot_t[i], N_zero_t[i]/N_tot_t[i], N_plusone_t[i]/N_tot_t[i], dp_avg_t[i]
    
#     # Remove the brackets from temp file [ and ]
#     with open('temp', 'r') as my_file:
#         text = my_file.read()
#         text = text.replace("[", "")
#         text = text.replace("]", "")     
   
#     # Use proper name for proper sub-case for CASE B
#     file2 = 'caseB_' + 'f_' + str(int(factor)) + '_dp_' + str(int(dp_zero_0*1e9)) + '_T_' + str(int(T))
#     with open(file2, 'w') as my_file:
#         my_file.write(text)
    
#     # Remove the temp file from the folder
#     os.remove('temp')    
    
#     plt.figure(1)
#     plt.grid(True)
#     plt.xlabel('Time (in s)')
#     plt.ylabel('Particle Number Concentrations (#/m3)')
#     line_plus, = plt.plot(t, N_plusone_t,'b+', label='+')
#     line_minus, = plt.plot(t, N_minusone_t,'b-', label='-')
#     line_zero, = plt.plot(t, N_zero_t, 'g*', label = '0')
#     plt.legend(handles=[line_plus, line_minus, line_zero])
#     plt.title('Comparison of positive, negative and neutral particles concentrations')
#     plt.show()     
    
#     plt.figure(2)
#     plt.grid(True)
#     plt.xlabel('Time (in s)')
#     plt.ylabel('Ion Concentrations (#/m3)')
#     line_plus, = plt.plot(t, n_plus_t ,'g+', label='+')
#     line_minus, = plt.plot(t, n_minus_t ,'b-', label='-')
#     plt.legend(handles=[line_plus, line_minus])
#     plt.title('Comparison of positive and negative ions concentrations')
#     plt.show()    


if __name__ == '__main__':
    global x_data_pos
    global y_data_pos
    global x_data_neg
    global y_data_neg
    x_data_pos, y_data_pos = xydata_pos()
    x_data_neg, y_data_neg = xydata_neg()
#     s = raw_input('bass')  
    main()
  

    
    
    
    
    
    
    
    
    