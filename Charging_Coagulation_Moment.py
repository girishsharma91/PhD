#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Written by Girish Sharma under the guidance of Prof. Pratim Biswas
# This code solves:
# 29 coupled differential equations viz.  3*9 for -4,-3,-2,-1,0,+1,+2,+3,+4 charged particles; 2 for ions


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
sigmag0 = 1.0
sigmag_ss = 1.35
T  = 2400 # can take values 300, 600, 1200, 2400                                  # in Kelvin
CHAR = 1
COAG = 1
QYN  = 0  # 0 = no change in ion concentration, 1 = solve differential equation for ion balance to see the changing ion concentration


# initial conditions
ion = 16
par = 20
n_m_0      = 1.0*10**(ion)
n_p_0      = 1.0*10**(ion)

N_0_0  = 1.0*10**(par)
qm = 4 # maximum charge

T = float(T)
P  = 101325.0                                   # in N/m2
M = 80.0*1e-3/Na
rho = 4230.0 

Mg = 28.97                                      # Molecular Weight for air
mu = 1.716e-5 * (T/273)**(2.0/3)                # viscosity of air   
lam = mu/P * ma.sqrt(np.pi*kb*T/(2.0*Mg/Na))    # mean free path in air(m)
print lam
############## FOR NON-DIMENSIONAL ######################
if CHAR == 1 and COAG == 1:
    file_add = '_coag_char'
    NT = N_0_0*n_p_0/(N_0_0+n_p_0)
elif CHAR == 1 and COAG == 0:
    file_add = '_char'
    NT = n_p_0
elif CHAR == 0 and COAG == 1:
    file_add = '_coag'
    NT = N_0_0
else:
    pass
KK = 8.0*kb*T/(3.0*mu)
tau = 1.0/(KK*NT)

# NT = 1.0
# KK = 1.0
# tau = 1.0


if QYN == 1:
    file_add = file_add + '_qyn1' 
else:
    file_add = file_add + '_qyn0' 


q_neg = 0.0/(NT*NT)                                     # rate of production of negative ion
q_pos = 0.0/(NT*NT)                                     # rate of production of positive ion
alpha = 1e-13


# for integration
t_final = 1.0e-2/tau # time in seconds
t_start = 1e-5/tau
delta_t = 1.1


M0_m4_0 = 0.0 # Moment (zero) of charge (minus 4) at time (zero)
M0_m3_0 = 0.0
M0_m2_0 = 0.0
M0_m1_0 = 0.0
M0_0_0 = N_0_0
M0_p1_0 = 0.0
M0_p2_0 = 0.0
M0_p3_0 = 0.0
M0_p4_0 = 0.0


v_zero_0     = ma.pi/6 * (dp_zero_0)**3
v1 = v_zero_0

M1_m4_0 = 0.0
M1_m3_0 = 0.0
M1_m2_0 = 0.0
M1_m1_0 = 0.0
M1_0_0  = v1*N_0_0
M1_p1_0 = 0.0
M1_p2_0 = 0.0
M1_p3_0 = 0.0
M1_p4_0 = 0.0

M2_m4_0 = 0.0
M2_m3_0 = 0.0
M2_m2_0 = 0.0
M2_m1_0 = 0.0
M2_0_0  = v1*v1*N_0_0
M2_p1_0 = 0.0
M2_p2_0 = 0.0
M2_p3_0 = 0.0
M2_p4_0 = 0.0

# Initial = [ M0_m4_0, M0_m3_0, M0_m2_0, M0_m1_0, M0_0_0, M0_p1_0, M0_p2_0, M0_p3_0, M0_p4_0, \
#             M1_m4_0, M1_m3_0, M1_m2_0, M1_m1_0, M1_0_0, \
#             M1_p1_0, M1_p2_0, M1_p3_0, M1_p4_0, M2_m4_0, \
#             M2_m3_0, M2_m2_0, M2_m1_0,\
#             M2_0_0, M2_p1_0, M2_p2_0, \
#             M2_p3_0, M2_p4_0, n_m_0, n_p_0]

Initial = [ M0_m4_0/NT, M0_m3_0/NT, M0_m2_0/NT, M0_m1_0/NT, M0_0_0/NT, M0_p1_0/NT, M0_p2_0/NT, M0_p3_0/NT, M0_p4_0/NT, \
            M1_m4_0/(NT*v_zero_0), M1_m3_0/(NT*v_zero_0), M1_m2_0/(NT*v_zero_0), M1_m1_0/(NT*v_zero_0), M1_0_0/(NT*v_zero_0), \
            M1_p1_0/(NT*v_zero_0), M1_p2_0/(NT*v_zero_0), M1_p3_0/(NT*v_zero_0), M1_p4_0/(NT*v_zero_0), M2_m4_0/(NT*v_zero_0*v_zero_0), \
            M2_m3_0/(NT*v_zero_0*v_zero_0), M2_m2_0/(NT*v_zero_0*v_zero_0), M2_m1_0/(NT*v_zero_0*v_zero_0),\
            M2_0_0/(NT*v_zero_0*v_zero_0), M2_p1_0/(NT*v_zero_0*v_zero_0), M2_p2_0/(NT*v_zero_0*v_zero_0), \
            M2_p3_0/(NT*v_zero_0*v_zero_0), M2_p4_0/(NT*v_zero_0*v_zero_0), n_m_0/NT, n_p_0/NT]
# # raw_input('ii')

b = -0.42/3
# Inputs:
# M0,M1,M2 vectors with all the charges, q charge, kth moment
def mom(M0_in,M1_in,M2_in,q_in,k_in):
    q = int(q_in)
    k = float(k_in)
    M0 = M0_in[q+qm]
    M1 = M1_in[q+qm]
    M2 = M2_in[q+qm]
    
    if M0*NT < 1.0:
        Moment = 0.0
    else:
        Moment =   M0**(1.0-3.0*k/2+k*k/2)*M1**(2.0*k-k*k)*M2**(-k/2+k*k/2) # 
    
    return Moment



############################# HELP FOR COAGULATION PART ##############################################
coag_coeff1 = (3.0/(4.0*ma.pi))**(1.0/6) * (6.0*kb*T/rho)**(1.0/2)


def Q1(i,j):
    if i*j > 0:
        ans =  (6.0/ma.pi)**(1.0/3) * 0.5968 * 1.0/(i*j) * ((kb*T)/(2.0*ke*e*e)) * (v_zero_0**(1.0/2))
    else:
        ans = 0.0
    return ans
    
def Q2(i,j):
    if i*j < 0:
        ans = (6.0/ma.pi)**(2.0*b) * 0.6205 * 1.0/abs(i*j) * ((kb*T)/(2.0*ke*e*e))**(6.0*b) * (v_zero_0**(1.0/6 + 2.0*b))
    else:
        ans = 0.0
    return ans    

def Q0(i,j):
    if i*j == 0:
        ans =  1.0 * (v_zero_0**(1.0/6))
    else:
        ans = 0.0
    return ans


###################################### HELP FOR CHARGING PART #############################################

def AB():
    
    # Assumption is that positive ion will not attach on negatively charged particles and vice versa. FOR SIMPLFICATION,
    # can be modified later.
    p_cons = ke*e*e/(kb*T)
    K = 1.0
    P = 1.01e5
    T0 = 298.0
    P0 = 1.01e5
    Mp = 150.0
    Mn = 120.0
    mp = Mp*1e-3/6.022e23
    mn = Mn*1e-3/6.022e23
    zp0 = ma.exp(-0.0347*(ma.log(Mp))**2-0.0376*ma.log(Mp)+1.4662)*1e-4    #Kilpatrick
    zn0 = ma.exp(-0.0347*(ma.log(Mn))**2-0.0376*ma.log(Mn)+1.4662)*1e-4    #Kilpatrick
    zp = zp0*P0*T/P/T0       #change with temperature and pressure
    zn = zn0*P0*T/P/T0       #change with temperature and pressure
    M = 28.9*1e-3/6.02e23    #air molecule mass
    diffp = kb*T*zp/e        #positive ion diffusion coefficient
    diffn = kb*T*zn/e        #negative ion diffusion coefficient
    
    cp = (8.0*kb*T/3.14/mp)**0.5    #positive ion velocity;
    cn = (8.0*kb*T/3.14/mn)**0.5    #negative ion velocity;

    lambdap = 16.0*2.0**0.5*diffp/3.0/3.14/cp*(M/(M+mp))**0.5   #positive ion mfp, or 32*diffp/3/3.14/cp*M/(M+mp)
    lambdan = 16.0*2.0**0.5*diffn/3.0/3.14/cn*(M/(M+mn))**0.5 
    
    
    Ap = np.zeros((2*qm +1))
    Bp = np.zeros((2*qm +1))
    An = np.zeros((2*qm +1))
    Bn = np.zeros((2*qm +1)) 
    x  = np.zeros((2*qm +1)) 
    
    for i in range(-qm,1):
        if (i == 0):
            Ap[i+qm] = (3.0/(4.0*ma.pi))**(2.0/3) * ma.pi * cp/(1+p_cons*i/lambdap + 1.0/2* (p_cons*i/lambdap)**2) 
            An[i+qm] = (3.0/(4.0*ma.pi))**(2.0/3) * ma.pi * cn/(1+p_cons*i/lambdan + 1.0/2* (p_cons*i/lambdan)**2) 
                
            Bp[i+qm] = (3.0/(4.0*ma.pi))**(1.0/2) * ma.pi * cp/(1+p_cons*i/lambdap + 1.0/2* (p_cons*i/lambdap)**2) \
                       * (2.0*ma.sqrt(p_cons*K/3))
            Bn[i+qm] = (3.0/(4.0*ma.pi))**(1.0/2) * ma.pi * cn/(1+p_cons*i/lambdan + 1.0/2* (p_cons*i/lambdan)**2) \
                       * (2.0*ma.sqrt(p_cons*K/3))
            
            Ap[i+qm] = (v_zero_0)**(2.0/3) * Ap[i+qm]
            An[i+qm] = (v_zero_0)**(2.0/3) * An[i+qm]
                
            Bp[i+qm] = (v_zero_0)**(1.0/2) * Bp[i+qm]
            Bn[i+qm] = (v_zero_0)**(1.0/2) * Bn[i+qm]
                
            
        else:
            x[i+qm]  = 1.0 + 1.0/2*K/(ma.sqrt(abs(i)))
            Ap[i+qm] = (3.0/(4.0*ma.pi))**(2.0/3) * ma.pi * cp/(1+p_cons*i/lambdap + 1.0/2* (p_cons*i/lambdap)**2) \
                        *(x[i+qm]**2*(1+2.0/3*p_cons*i/lambdap))
            An[i+qm] = (3.0/(4.0*ma.pi))**(2.0/3) * ma.pi * cn/(1+p_cons*i/lambdan + 1.0/2* (p_cons*i/lambdan)**2) \
                        *(x[i+qm]**2*(1+2.0/3*p_cons*i/lambdan))
                
            Bp[i+qm] = (3.0/(4.0*ma.pi))**(1.0/3) * ma.pi * cp/(1+p_cons*i/lambdap + 1.0/2* (p_cons*i/lambdap)**2) \
                       *(1.0/3*p_cons*K/(x[i+qm]*x[i+qm]-1) - 2.0/3 * p_cons*i * x[i+qm] )
            Bn[i+qm] = (3.0/(4.0*ma.pi))**(1.0/3) * ma.pi * cn/(1+p_cons*i/lambdan + 1.0/2* (p_cons*i/lambdan)**2) \
                       *(1.0/3*p_cons*K/(x[i+qm]*x[i+qm]-1) - 2.0/3 * p_cons*i * x[i+qm] ) 
                
            Ap[i+qm] = (v_zero_0)**(2.0/3) * Ap[i+qm]
            An[i+qm] = (v_zero_0)**(2.0/3) * An[i+qm]
                
            Bp[i+qm] = (v_zero_0)**(1.0/3) * Bp[i+qm]
            Bn[i+qm] = (v_zero_0)**(1.0/3) * Bn[i+qm]
        
#             if (Ap[i+qm]*vg[i+qm]**(2.0/3) + Bp[i+qm]*vg[i+qm]**(1.0/3)) < 0:
#                 Ap[i+qm] = 0.0
#                 An[i+qm] = 0.0
#                 Bp[i+qm] = 0.0
#                 Bn[i+qm] = 0.0    
    
           
    return Ap,An,Bp,Bn      

# RIGHT HAND SIDE FOR THE SOLVER
def rhs(t,y):
    ######################################## INITILAIZING ####################################
    # find the important quantities from y
    M0 = y[0:2*qm + 1]
    M1 = y[2*qm + 1:2*(2*qm +1)]
    M2 = y[2*(2*qm +1):3*(2*qm +1)]
    n_neg = y[3*(2*qm +1)]
    n_pos = y[3*(2*qm +1)+1]
    
    
    sigmag_rhs = np.zeros((2*qm+1))
    temp_sigmag = 0.0
    
    for i in range(0,2*qm+1):
        if M0[i] * NT < 1 or M0[i]*M2[i]/M1[i]**2.0 < 1:
            sigmag_rhs[i] = 1.0
        else:
            sigmag_rhs[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0[i]*M2[i]/M1[i]**2.0)))

        temp_sigmag =  temp_sigmag + (sigmag_rhs[i]*M0[i])
    
    sigmag_avg = temp_sigmag/np.sum(M0)
    coag_coeff2 = 0.39 + 0.5 * sigmag_avg - 0.214 * sigmag_avg**(2.0) + 0.029 * (sigmag_avg)**(3.0)
   

            
            
            
            
            
            
            
            
    dM0 = np.zeros((2*qm +1))
    gain_M0 = np.zeros((2*qm +1))
    loss_M0 = np.zeros((2*qm +1))
    dM0_coag = np.zeros((2*qm +1))
    dM0_char = np.zeros((2*qm +1))
    
    dM1 = np.zeros((2*qm +1))
    gain_M1 = np.zeros((2*qm +1))
    loss_M1 = np.zeros((2*qm +1))
    dM1_coag = np.zeros((2*qm +1))
    dM1_char = np.zeros((2*qm +1))
    
    dM2 = np.zeros((2*qm +1))
    gain_M2 = np.zeros((2*qm +1))
    loss_M2 = np.zeros((2*qm +1))
    dM2_coag = np.zeros((2*qm +1))
    dM2_char = np.zeros((2*qm +1))
    
#     sumN1 = np.zeros((2*qm +1))
#     sumN2 = np.zeros((2*qm +1))
    
#     sumV1 = np.zeros((2*qm +1))
#     sumV2 = np.zeros((2*qm +1))
    
    
#     ############################################ COAGULATION #######################################
#     for k in range(-qm,qm+1):
        
#         for i in range(-qm,qm+1):
            
#             sumN2[k+qm] = sumN2[k+qm] + beta(i,k,dp[i+qm],dp[k+qm]) * N[i+qm] * N[k+qm]
#             sumV2[k+qm] = sumV2[k+qm] + beta(i,k,dp[i+qm],dp[k+qm]) * N[i+qm] * N[k+qm] * v[k+qm]
            
#             for j in range(-qm,qm+1):
                
#                 if i + j == k:
#                     sumN1[k+qm] = sumN1[k+qm] + beta(i,j,dp[i+qm],dp[j+qm]) * N[i+qm] * N[j+qm]
#                     sumV1[k+qm] = sumV1[k+qm] + beta(i,j,dp[i+qm],dp[j+qm]) * N[i+qm] * N[j+qm] * (v[i+qm] + v[j+qm])
# #                     print i,j,k
# #                     print sumN1[k+qm], beta(i,j,dp[i+qm],dp[j+qm]), N[i+qm], N[j+qm]
#                 else:
#                     pass
        
#         dN_coag[k+qm] = 1.0/2 * sumN1[k+qm] - sumN2[k+qm]
#         dV_coag[k+qm] = 1.0/2 * sumV1[k+qm] - sumV2[k+qm]


    for k in range(-qm,qm+1):
        q=k
        gain_M0[k+qm] = 0.0
        gain_M1[k+qm] = 0.0
        gain_M2[k+qm] = 0.0
    
        loss_M0[k+qm] = 0.0
        loss_M1[k+qm] = 0.0
        loss_M2[k+qm] = 0.0
        for i in range(-qm,qm+1):
            for j in range(-qm,qm+1):
                if i + j == k:
                    
                    gain_M0[k+qm] = gain_M0[k+qm] + \
                                    (Q1(i,j)*( mom(M0,M1,M2,i,1.0/2)*mom(M0,M1,M2,j,0.0) + 3.0*mom(M0,M1,M2,i,-1.0/6)*mom(M0,M1,M2,j,2.0/3) + 
                                    3.0*mom(M0,M1,M2,i,1.0/6)*mom(M0,M1,M2,j,1.0/3) + mom(M0,M1,M2,i,1.0)*mom(M0,M1,M2,j,-1.0/2) ) + 
                                    Q2(i,j)*( mom(M0,M1,M2,i,1.0/6+b)*mom(M0,M1,M2,j,b) + mom(M0,M1,M2,i,-1.0/2+b)*mom(M0,M1,M2,j,2.0/3+b) 
                                    + 2.0*mom(M0,M1,M2,i,-1.0/6+b)*mom(M0,M1,M2,j,1.0/3+b) ) 
                                    + Q0(i,j)*( mom(M0,M1,M2,i,1.0/6)*mom(M0,M1,M2,j,0.0) + mom(M0,M1,M2,i,-1.0/2)*mom(M0,M1,M2,j,2.0/3+0.0) 
                                    + 2.0*mom(M0,M1,M2,i,-1.0/6)*mom(M0,M1,M2,j,1.0/3) ))

                    gain_M1[k+qm] = gain_M1[k+qm] + \
                                    (Q1(i,j)*( mom(M0,M1,M2,i,3.0/2)*mom(M0,M1,M2,j,0.0) + 3.0*mom(M0,M1,M2,i,5.0/6)*mom(M0,M1,M2,j,2.0/3) 
                                    + 3.0*mom(M0,M1,M2,i,7.0/6)*mom(M0,M1,M2,j,1.0/3) + mom(M0,M1,M2,i,2.0)*mom(M0,M1,M2,j,-1.0/2)
                                    + 3.0*mom(M0,M1,M2,i,4.0/3)*mom(M0,M1,M2,j,1.0/6) + 3.0*mom(M0,M1,M2,i,5.0/3)*mom(M0,M1,M2,j,-1.0/6)
                                    + 2.0*mom(M0,M1,M2,i,1.0/2)*mom(M0,M1,M2,j,1.0)  )
                                    + Q2(i,j)* ( mom(M0,M1,M2,i,7.0/6+b)*mom(M0,M1,M2,j,b) + mom(M0,M1,M2,i,1.0/2+b)*mom(M0,M1,M2,j,2.0/3+b)
                                    + 2.0*mom(M0,M1,M2,i,5.0/6+b)*mom(M0,M1,M2,j,1.0/3+b) + mom(M0,M1,M2,i,5.0/3+b)*mom(M0,M1,M2,j,-1.0/2+b)
                                    + mom(M0,M1,M2,i,1.0+b)*mom(M0,M1,M2,j,1.0/6+b) + 2.0*mom(M0,M1,M2,i,4.0/3+b)*mom(M0,M1,M2,j,-1.0/6+b) )
                                    + Q0(i,j)* ( mom(M0,M1,M2,i,7.0/6)*mom(M0,M1,M2,j,0) + mom(M0,M1,M2,i,1.0/2)*mom(M0,M1,M2,j,2.0/3)
                                    + 2.0*mom(M0,M1,M2,i,5.0/6)*mom(M0,M1,M2,j,1.0/3) + mom(M0,M1,M2,i,5.0/3)*mom(M0,M1,M2,j,-1.0/2)
                                    + mom(M0,M1,M2,i,1.0)*mom(M0,M1,M2,j,1.0/6) + 2.0*mom(M0,M1,M2,i,4.0/3)*mom(M0,M1,M2,j,-1.0/6) ))

                    gain_M2[k+qm] = gain_M2[k+qm] + \
                                    (Q1(i,j)*( mom(M0,M1,M2,i,5.0/2)*mom(M0,M1,M2,j,0.0) + 3.0*mom(M0,M1,M2,i,11.0/6)*mom(M0,M1,M2,j,2.0/3) 
                                    + 3.0*mom(M0,M1,M2,i,13.0/6)*mom(M0,M1,M2,j,1.0/3) + mom(M0,M1,M2,i,3.0)*mom(M0,M1,M2,j,-1.0/2)
                                    + 3.0*mom(M0,M1,M2,i,7.0/3)*mom(M0,M1,M2,j,1.0/6) + 3.0*mom(M0,M1,M2,i,8.0/3)*mom(M0,M1,M2,j,-1.0/6)
                                    + 3.0*mom(M0,M1,M2,i,3.0/2)*mom(M0,M1,M2,j,1.0) + 6.0*mom(M0,M1,M2,i,5.0/6)*mom(M0,M1,M2,j,5.0/3)
                                    + 6.0*mom(M0,M1,M2,i,7.0/6)*mom(M0,M1,M2,j,4.0/3) + 3.0*mom(M0,M1,M2,i,2.0)*mom(M0,M1,M2,j,1.0/2))
                                    + Q2(i,j)*( mom(M0,M1,M2,i,13.0/6+b)*mom(M0,M1,M2,j,b) + mom(M0,M1,M2,i,3.0/2+b)*mom(M0,M1,M2,j,2.0/3+b)
                                    + 2.0*mom(M0,M1,M2,i,11.0/6+b)*mom(M0,M1,M2,j,1.0/3+b) + 2.0*mom(M0,M1,M2,i,7.0/6+b)*mom(M0,M1,M2,j,1.0+b)
                                    + 2.0*mom(M0,M1,M2,i,1.0/2+b)*mom(M0,M1,M2,j,5.0/3+b) + 4.0*mom(M0,M1,M2,i,5.0/6+b)*mom(M0,M1,M2,j,4.0/3+b)
                                    + mom(M0,M1,M2,i,1.0/6+b)*mom(M0,M1,M2,j,2.0+b) + mom(M0,M1,M2,i,-1.0/2+b)*mom(M0,M1,M2,j,8.0/3+b)
                                    + 2.0*mom(M0,M1,M2,i,-1.0/6+b)*mom(M0,M1,M2,j,7.0/3+b) )
                                    + Q0(i,j)*( mom(M0,M1,M2,i,13.0/6)*mom(M0,M1,M2,j,0) + mom(M0,M1,M2,i,3.0/2)*mom(M0,M1,M2,j,2.0/3)
                                    + 2.0*mom(M0,M1,M2,i,11.0/6)*mom(M0,M1,M2,j,1.0/3) + 2.0*mom(M0,M1,M2,i,7.0/6)*mom(M0,M1,M2,j,1.0)
                                    + 2.0*mom(M0,M1,M2,i,1.0/2)*mom(M0,M1,M2,j,5.0/3) + 4.0*mom(M0,M1,M2,i,5.0/6)*mom(M0,M1,M2,j,4.0/3)
                                    + mom(M0,M1,M2,i,1.0/6)*mom(M0,M1,M2,j,2.0) + mom(M0,M1,M2,i,-1.0/2)*mom(M0,M1,M2,j,8.0/3)
                                    + 2.0*mom(M0,M1,M2,i,-1.0/6)*mom(M0,M1,M2,j,7.0/3) ))

            loss_M0[k+qm] = loss_M0[k+qm] + \
                            (Q1(q,i)* ( mom(M0,M1,M2,q,0.0)*mom(M0,M1,M2,i,1.0/2) + 3.0*mom(M0,M1,M2,q,2.0/3)*mom(M0,M1,M2,i,-1.0/6) + 
                            3.0*mom(M0,M1,M2,q,1.0/3)*mom(M0,M1,M2,i,1.0/6) + mom(M0,M1,M2,q,-1.0/2)*mom(M0,M1,M2,i,1.0) +
                            mom(M0,M1,M2,q,1.0/2)*mom(M0,M1,M2,i,0.0) + 3.0*mom(M0,M1,M2,q,-1.0/6)*mom(M0,M1,M2,i,2.0/3) +
                            3.0*mom(M0,M1,M2,q,1.0/6)*mom(M0,M1,M2,i,1.0/3) + mom(M0,M1,M2,q,1.0)*mom(M0,M1,M2,i,-1.0/2))
                            + Q2(q,i)* ( mom(M0,M1,M2,q,1.0/6+b)*mom(M0,M1,M2,i,b) + mom(M0,M1,M2,q,-1.0/2+b)*mom(M0,M1,M2,i,2.0/3+b)
                            + 2.0*mom(M0,M1,M2,q,-1.0/6+b)*mom(M0,M1,M2,i,1.0/3+b) + mom(M0,M1,M2,q,1.0/6+b)*mom(M0,M1,M2,i,b) +
                            mom(M0,M1,M2,q,-1.0/2+b)*mom(M0,M1,M2,i,2.0/3+b) + 2.0*mom(M0,M1,M2,q,-1.0/6+b)*mom(M0,M1,M2,i,1.0/3+b))
                            + Q0(q,i)* ( mom(M0,M1,M2,q,1.0/6)*mom(M0,M1,M2,i,0) + mom(M0,M1,M2,q,-1.0/2)*mom(M0,M1,M2,i,2.0/3)
                            + 2.0*mom(M0,M1,M2,q,-1.0/6)*mom(M0,M1,M2,i,1.0/3) + mom(M0,M1,M2,q,1.0/6)*mom(M0,M1,M2,i,0) +
                            mom(M0,M1,M2,q,-1.0/2)*mom(M0,M1,M2,i,2.0/3) + 2.0*mom(M0,M1,M2,q,-1.0/6)*mom(M0,M1,M2,i,1.0/3)))
            loss_M1[k+qm] = loss_M1[k+qm] +  \
                            (Q1(q,i)* ( mom(M0,M1,M2,q,3.0/2)*mom(M0,M1,M2,i,0.0) + 3.0*mom(M0,M1,M2,q,5.0/6)*mom(M0,M1,M2,i,2.0/3) +
                            3.0*mom(M0,M1,M2,q,7.0/6)*mom(M0,M1,M2,i,1.0/3) + mom(M0,M1,M2,q,2.0)*mom(M0,M1,M2,i,-1.0/2) +
                            3.0*mom(M0,M1,M2,q,4.0/3)*mom(M0,M1,M2,i,1.0/6) + 3.0*mom(M0,M1,M2,q,5.0/3)*mom(M0,M1,M2,i,-1.0/6) +
                            mom(M0,M1,M2,q,1.0/2)*mom(M0,M1,M2,i,1.0) + mom(M0,M1,M2,q,1.0)*mom(M0,M1,M2,i,1.0/2) )
                            + Q2(q,i)* ( mom(M0,M1,M2,q,7.0/6+b)*mom(M0,M1,M2,i,b) + mom(M0,M1,M2,q,1.0/2+b)*mom(M0,M1,M2,i,2.0/3+b) +
                            2.0*mom(M0,M1,M2,q,5.0/6+b)*mom(M0,M1,M2,i,1.0/3+ b) +  mom(M0,M1,M2,q,5.0/3+b)*mom(M0,M1,M2,i,-1.0/2+ b) +
                            mom(M0,M1,M2,q,1.0+b)*mom(M0,M1,M2,i,1.0/6+ b) +  2.0*mom(M0,M1,M2,q,4.0/3+b)*mom(M0,M1,M2,i,-1.0/6+ b))
                            + Q0(q,i)* ( mom(M0,M1,M2,q,7.0/6)*mom(M0,M1,M2,i,0) + mom(M0,M1,M2,q,1.0/2)*mom(M0,M1,M2,i,2.0/3) +
                            2.0*mom(M0,M1,M2,q,5.0/6)*mom(M0,M1,M2,i,1.0/3) +  mom(M0,M1,M2,q,5.0/3)*mom(M0,M1,M2,i,-1.0/2) +
                            mom(M0,M1,M2,q,1.0)*mom(M0,M1,M2,i,1.0/6) +  2.0*mom(M0,M1,M2,q,4.0/3)*mom(M0,M1,M2,i,-1.0/6)))
            loss_M2[k+qm] = loss_M2[k+qm] +  \
                            (Q1(q,i)* ( mom(M0,M1,M2,q,5.0/2)*mom(M0,M1,M2,i,0.0) + 3.0*mom(M0,M1,M2,q,11.0/6)*mom(M0,M1,M2,i,2.0/3) +
                            3.0*mom(M0,M1,M2,q,13.0/6)*mom(M0,M1,M2,i,1.0/3) + mom(M0,M1,M2,q,3.0)*mom(M0,M1,M2,i,-1.0/2) +
                            3.0* mom(M0,M1,M2,q,7.0/3)*mom(M0,M1,M2,i,1.0/6) + 3.0*mom(M0,M1,M2,q,8.0/3)*mom(M0,M1,M2,i,-1.0/6) + 
                            mom(M0,M1,M2,q,3.0/2)*mom(M0,M1,M2,i,1.0) + mom(M0,M1,M2,q,2.0)*mom(M0,M1,M2,i,1.0/2) )
                            + Q2(q,i)* ( mom(M0,M1,M2,q,13.0/6+b)*mom(M0,M1,M2,i,b) + mom(M0,M1,M2,q,3.0/2+b)*mom(M0,M1,M2,i,2.0/3+b) +
                            2.0*mom(M0,M1,M2,q,11.0/6+b)*mom(M0,M1,M2,i,1.0/3+b) + mom(M0,M1,M2,q,8.0/3+b)*mom(M0,M1,M2,i,-1.0/2+b) +
                            mom(M0,M1,M2,q,2.0+b)*mom(M0,M1,M2,i,1.0/6+b) + 2.0*mom(M0,M1,M2,q,7.0/3+b)*mom(M0,M1,M2,i,-1.0/6+b))
                            + Q0(q,i)* ( mom(M0,M1,M2,q,13.0/6)*mom(M0,M1,M2,i,0) + mom(M0,M1,M2,q,3.0/2)*mom(M0,M1,M2,i,2.0/3) +
                            2.0*mom(M0,M1,M2,q,11.0/6)*mom(M0,M1,M2,i,1.0/3) + mom(M0,M1,M2,q,8.0/3)*mom(M0,M1,M2,i,-1.0/2) +
                            mom(M0,M1,M2,q,2.0)*mom(M0,M1,M2,i,1.0/6) + 2.0*mom(M0,M1,M2,q,7.0/3)*mom(M0,M1,M2,i,-1.0/6)))
        dM0_coag[k+qm] =  coag_coeff1 * coag_coeff2 *(gain_M0[k+qm] - loss_M0[k+qm]) 
        dM1_coag[k+qm] =  coag_coeff1 * coag_coeff2 *(gain_M1[k+qm] - loss_M1[k+qm])
        dM2_coag[k+qm] =  coag_coeff1 * coag_coeff2 *(gain_M2[k+qm] - loss_M2[k+qm])
                
        
        
    ############################################# CHARGING #########################################
    #### For Number and Volume Concentration:
    
    
    [Ap,An,Bp,Bn]=   AB() 
    #print 'Ap, Bp' , Ap, Bp
    
    dM0_char[-qm+qm] = n_neg * (An[-(qm-1)+qm] * mom(M0,M1,M2,-(qm-1),2.0/3) + Bn[-(qm-1)+qm] *mom(M0,M1,M2,-(qm-1),1.0/3)) \
                       -n_pos * (Ap[-qm+qm] * mom(M0,M1,M2,-qm,2.0/3) + Bp[-qm+qm] * mom(M0,M1,M2,-qm,1.0/3))
    dM1_char[-qm+qm] = n_neg * (An[-(qm-1)+qm] * mom(M0,M1,M2,-(qm-1),2.0/3+1) + Bn[-(qm-1)+qm] *mom(M0,M1,M2,-(qm-1),1.0/3+1)) \
                       -n_pos * (Ap[-qm+qm] * mom(M0,M1,M2,-qm,2.0/3+1) + Bp[-qm+qm] * mom(M0,M1,M2,-qm,1.0/3+1))
    dM2_char[-qm+qm] = n_neg * (An[-(qm-1)+qm] * mom(M0,M1,M2,-(qm-1),2.0/3+2) + Bn[-(qm-1)+qm] *mom(M0,M1,M2,-(qm-1),1.0/3+2)) \
                       -n_pos * (Ap[-qm+qm] * mom(M0,M1,M2,-qm,2.0/3+2) + Bp[-qm+qm] * mom(M0,M1,M2,-qm,1.0/3+2))
    
    
    for i in range(-(qm-1),(qm)):
        if (i == -1):
            dM0_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/2)) + \
                             n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3)) - \
                             n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3)) - \
                             n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3))
            dM1_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3+1) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/2+1)) + \
                             n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3+1) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3+1)) - \
                             n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3+1) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3+1)) - \
                             n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3+1) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3+1))
            dM2_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3+2) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/2+2)) + \
                             n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3+2) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3+2)) - \
                             n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3+2) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3+2)) - \
                             n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3+2) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3+2))
            
            
        elif (i == +1):
            dM0_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3)) \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/2))\
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3))\
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3))
            dM1_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3+1) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3+1)) \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3+1) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/2+1))\
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3+1) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3+1))\
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3+1) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3+1))
            dM2_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3+2) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3+2)) \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3+2) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/2+2))\
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3+2) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3+2))\
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3+2) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3+2))
            
            
        elif (i == 0):
            dM0_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3)) \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3)) \
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/2)) \
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/2))
            dM1_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3+1) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3+1))  \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3+1) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3+1))  \
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3+1) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/2+1))  \
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3+1) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/2+1))
            dM2_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3+2) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3+2))  \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3+2) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3+2))  \
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3+2) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/2+2))   \
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3+2) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/2+2))
        
        
        else:
            dM0_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3))  \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3))  \
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3))  \
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3))
            dM1_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3+1) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3+1))  \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3+1) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3+1))  \
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3+1) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3+1))  \
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3+1) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3+1))
            dM2_char[i+qm] = n_neg * (An[i+1+qm] * mom(M0,M1,M2,i+1,2.0/3+2) + Bn[i+1+qm] *mom(M0,M1,M2,i+1,1.0/3+2))  \
                          +  n_pos * (Ap[i-1+qm] * mom(M0,M1,M2,i-1,2.0/3+2) + Bp[i-1+qm] *mom(M0,M1,M2,i-1,1.0/3+2))  \
                          -  n_neg * (An[i+qm] * mom(M0,M1,M2,i,2.0/3+2) + Bn[i+qm] *mom(M0,M1,M2,i,1.0/3+2))   \
                          -  n_pos * (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3+2) + Bp[i+qm] *mom(M0,M1,M2,i,1.0/3+2))


    dM0_char[+qm+qm] = n_pos * (An[+(qm-1)+qm] * mom(M0,M1,M2,+(qm-1),2.0/3) + Bn[+(qm-1)+qm] *mom(M0,M1,M2,+(qm-1),1.0/3)) \
                      -n_neg * (Ap[+qm+qm] * mom(M0,M1,M2,+qm,2.0/3) + Bp[+qm+qm] * mom(M0,M1,M2,+qm,1.0/3))
    dM1_char[+qm+qm] = n_pos * (An[+(qm-1)+qm] * mom(M0,M1,M2,+(qm-1),2.0/3+1) + Bn[+(qm-1)+qm] *mom(M0,M1,M2,+(qm-1),1.0/3+1)) \
                      -n_neg * (Ap[+qm+qm] * mom(M0,M1,M2,+qm,2.0/3+1) + Bp[+qm+qm] * mom(M0,M1,M2,+qm,1.0/3+1))
    dM2_char[+qm+qm] = n_pos * (An[+(qm-1)+qm] * mom(M0,M1,M2,+(qm-1),2.0/3+2) + Bn[+(qm-1)+qm] *mom(M0,M1,M2,+(qm-1),1.0/3+2)) \
                      -n_neg * (Ap[+qm+qm] * mom(M0,M1,M2,+qm,2.0/3+2) + Bp[+qm+qm] * mom(M0,M1,M2,+qm,1.0/3+2)) 
        
 
    #print 'dM0_char[0] = ', dM0_char[qm] 
    #print 'dM0_coag[0] = ', dM0_coag[qm]
    
    ##### For Ion Number Concentration:     
    temp_sum_neg = 0.0
    temp_sum_pos = 0.0
    for i in range(-qm+1,qm):
        if (i == 0):
            temp_sum_neg = temp_sum_neg + (An[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bn[0+qm] *mom(M0,M1,M2,i,1.0/2))
        else:
            temp_sum_neg = temp_sum_neg + (An[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bn[i+qm] *mom(M0,M1,M2,i,2.0/3)) 
            
    for i in range(-qm,qm):
        if (i == 0):
            temp_sum_pos = temp_sum_pos + (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bp[0+qm] *mom(M0,M1,M2,i,1.0/2))
        else:
            temp_sum_pos = temp_sum_pos + (Ap[i+qm] * mom(M0,M1,M2,i,2.0/3) + Bp[i+qm] *mom(M0,M1,M2,i,2.0/3)) 
        
    if QYN == 1: 
        dn_neg = tau*NT * (q_neg - alpha*n_neg*n_pos -  n_neg * temp_sum_neg)
        dn_pos = tau*NT * (q_pos - alpha*n_neg*n_pos -  n_pos * temp_sum_pos) 
    else:
        dn_neg = 0.0
        dn_pos = 0.0


    ############################### COMBINING CHARGING AND COAGULATION ##################################
    if n_pos * NT < 100:
        CHARn = 0
        dn_neg = 0.0
        dn_pos = 0.0
    else:
        CHARn = CHAR
    
    #print CHARn
    if CHARn == 1 and COAG == 1:
        dM0 = tau*NT * (dM0_coag + dM0_char) 
        dM1 = tau*NT * (dM1_coag + dM1_char) 
        dM2 = tau*NT * (dM2_coag + dM2_char) 

    elif CHARn == 1 and COAG == 0:
        dM0 = tau*NT   *  (dM0_char) 
        dM1 = tau*NT   *  (dM1_char) 
        dM2 = tau*NT   *  (dM2_char) 

    elif CHARn == 0 and COAG == 1:
        dM0 = tau*NT * (dM0_coag) 
        dM1 = tau*NT * (dM1_coag) 
        dM2 = tau*NT * (dM2_coag)
#         if t * tau > .130:
#            # print dN
    else:
        pass
    
   
    ############################### STACKING ##################################    
    dn = np.asarray([dn_neg, dn_pos])
    dy = np.concatenate([dM0,dM1,dM2,dn])
    # print t*tau, dn_neg, n_neg
    
    return dy




################################################ MAIN PROGRAM ################################################
def main():
    
    # Start by specifying the integrator:
    # use ``vode`` with "backward differentiation formula"
    r = ode(rhs).set_integrator('vode', method='bdf', nsteps='1000000000', atol = '1e-4', rtol = '1e-7')   # atol = '1e-14', rtol = '1e-21'
 

    # Number of time steps: 1 extra for initial condition
#     num_steps = np.floor((t_final - t_start)/delta_t) + 1
    num_steps = np.floor(ma.log10(t_final/t_start)/ma.log10(delta_t))+2
    

    # Set initial condition(s): for integrating variable and time!
    r.set_initial_value(Initial, t_start)
    
 
    # Additional Python step: create vectors to store variables
    t = np.zeros((num_steps, 1))
    M0_m4_t = np.zeros((num_steps, 1))
    M0_m3_t = np.zeros((num_steps, 1))
    M0_m2_t = np.zeros((num_steps, 1))
    M0_m1_t = np.zeros((num_steps, 1))
    M0_0_t = np.zeros((num_steps, 1))
    M0_p1_t = np.zeros((num_steps, 1))
    M0_p2_t = np.zeros((num_steps, 1))
    M0_p3_t = np.zeros((num_steps, 1))
    M0_p4_t = np.zeros((num_steps, 1))
    

    M1_m4_t = np.zeros((num_steps, 1))
    M1_m3_t = np.zeros((num_steps, 1))
    M1_m2_t = np.zeros((num_steps, 1))
    M1_m1_t = np.zeros((num_steps, 1))
    M1_0_t = np.zeros((num_steps, 1))
    M1_p1_t = np.zeros((num_steps, 1))
    M1_p2_t = np.zeros((num_steps, 1))
    M1_p3_t = np.zeros((num_steps, 1))
    M1_p4_t = np.zeros((num_steps, 1))

    M2_m4_t = np.zeros((num_steps, 1))
    M2_m3_t = np.zeros((num_steps, 1))
    M2_m2_t = np.zeros((num_steps, 1))
    M2_m1_t = np.zeros((num_steps, 1))
    M2_0_t = np.zeros((num_steps, 1))
    M2_p1_t = np.zeros((num_steps, 1))
    M2_p2_t = np.zeros((num_steps, 1))
    M2_p3_t = np.zeros((num_steps, 1))
    M2_p4_t = np.zeros((num_steps, 1))
    
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
    
    sigmag_m4_t = np.zeros((num_steps, 1))
    sigmag_m3_t = np.zeros((num_steps, 1))
    sigmag_m2_t = np.zeros((num_steps, 1))
    sigmag_m1_t = np.zeros((num_steps, 1))
    sigmag_0_t = np.zeros((num_steps, 1))
    sigmag_p1_t = np.zeros((num_steps, 1))
    sigmag_p2_t = np.zeros((num_steps, 1))
    sigmag_p3_t = np.zeros((num_steps, 1))
    sigmag_p4_t = np.zeros((num_steps, 1))
    
    M0_tot_t = np.zeros((num_steps, 1))
    dp_avg_t= np.zeros((num_steps, 1))
    sigmag_avg_t= np.zeros((num_steps, 1))
    average_charge = np.zeros((num_steps,1))
    

    
    t[0] = t_start
    M0_m4_t[0] = Initial[0] 
    M0_m3_t[0] = Initial[1]  
    M0_m2_t[0] = Initial[2]  
    M0_m1_t[0] = Initial[3] 
    M0_0_t[0]  = Initial[4]  
    M0_p1_t[0] = Initial[5]  
    M0_p2_t[0] = Initial[6]  
    M0_p3_t[0] = Initial[7]  
    M0_p4_t[0] = Initial[8]  

    M1_m4_t[0] = Initial[9]  
    M1_m3_t[0] = Initial[10]  
    M1_m2_t[0] = Initial[11]  
    M1_m1_t[0] = Initial[12]  
    M1_0_t[0]  = Initial[13]  
    M1_p1_t[0] = Initial[14] 
    M1_p2_t[0] = Initial[15]  
    M1_p3_t[0] = Initial[16] 
    M1_p4_t[0] = Initial[17] 
    
    M2_m4_t[0] = Initial[18]  
    M2_m3_t[0] = Initial[19] 
    M2_m2_t[0] = Initial[20]  
    M2_m1_t[0] = Initial[21]  
    M2_0_t[0]  = Initial[22] 
    M2_p1_t[0] = Initial[23]  
    M2_p2_t[0] = Initial[24] 
    M2_p3_t[0] = Initial[25]  
    M2_p4_t[0] = Initial[26] 
    
    n_m_t[0]  = Initial[27] 
    n_p_t[0]  = Initial[28] 

 
    # Integrate the ODE(s) across each delta_t timestep
    k = 1
    while r.successful() and k < num_steps:
        r.integrate(r.t*delta_t)
 
        # Store the results to plot later
        t[k] = r.t
        print r.t*tau
        M0_m4_t[k] = r.y[0] 
        M0_m3_t[k] = r.y[1] 
        M0_m2_t[k] = r.y[2] 
        M0_m1_t[k] = r.y[3] 
        M0_0_t[k]  = r.y[4] 
        M0_p1_t[k] = r.y[5] 
        M0_p2_t[k] = r.y[6] 
        M0_p3_t[k] = r.y[7] 
        M0_p4_t[k] = r.y[8] 

        M1_m4_t[k] = r.y[9] 
        M1_m3_t[k] = r.y[10] 
        M1_m2_t[k] = r.y[11] 
        M1_m1_t[k] = r.y[12] 
        M1_0_t[k]  = r.y[13] 
        M1_p1_t[k] = r.y[14] 
        M1_p2_t[k] = r.y[15] 
        M1_p3_t[k] = r.y[16] 
        M1_p4_t[k] = r.y[17] 

        M2_m4_t[k] = r.y[18] 
        M2_m3_t[k] = r.y[19] 
        M2_m2_t[k] = r.y[20] 
        M2_m1_t[k] = r.y[21] 
        M2_0_t[k]  = r.y[22] 
        M2_p1_t[k] = r.y[23] 
        M2_p2_t[k] = r.y[24] 
        M2_p3_t[k] = r.y[25] 
        M2_p4_t[k] = r.y[26] 

        n_m_t[k]  = r.y[27]
        n_p_t[k]  = r.y[28]
 
        k += 1
    
    
#     file1 = 'Raw_data_' + str(int(T))
#     with open(file1, 'w') as f:
#         # Print & save the solution.
#         for i in range(0,k):
#             print >> f, t[i], n_p_t[i] N_minustwo_t[i], N_minusone_t[i], N_zero_t[i], N_plusone_t[i], N_plustwo_t[i], \
#                              # V_minustwo_t[i], V_minusone_t[i], V_zero_t[i], V_plusone_t[i], V_plustwo_t[i], \
#                              # n_neg_t[i],
    
    
    with open('rawtemp', 'w') as f:
        for i in range(0,k):
            print >> f, t[i]*tau, M0_m4_t[i], M0_m3_t[i], M0_m2_t[i], M0_m1_t[i], M0_0_t[i],\
            M0_p1_t[i], M0_p2_t[i], M0_p3_t[i], M0_p4_t[i], M1_m4_t[i], M1_m3_t[i], M1_m2_t[i], M1_m1_t[i], M1_0_t[i],\
            M1_p1_t[i], M1_p2_t[i], M1_p3_t[i], M1_p4_t[i], M2_m4_t[i], M2_m3_t[i], M2_m2_t[i], M2_m1_t[i], M2_0_t[i],\
            M2_p1_t[i], M2_p2_t[i], M2_p3_t[i], M2_p4_t[i], n_m_t[i], n_p_t[i]
    
    with open('rawtemp', 'r') as my_file:
        text = my_file.read()
        text = text.replace("[", "")
        text = text.replace("]", "")     
   
    
    file4 = 'raw_'   +   'ion_' + str(int(ion)) +   '_par_' + str(int(par)) + '_' + str(int(T)) + '_' + file_add +'.txt'
    with open(file4, 'w') as my_file:
        my_file.write(text)
    
    # Remove the temp file from the folder
    os.remove('rawtemp')    
    
    
    
    #########################   CASE A ##############################
    # creates a temp file to store all the data for CASE A temporarily
    with open('temp', 'w') as f:
        # Print & save the solution.
        for i in range(0,k):

            if M0_m4_t[i] * NT < 1:
                dp_m4_t[i] = 0
            else:
                dp_m4_t[i] = (6.0/ma.pi * (M1_m4_t[i]**(2.0)/(M0_m4_t[i]**(3.0/2) * M2_m4_t[i]**(1.0/2))) * v_zero_0)**(1.0/3)
            
            if M0_m3_t[i]* NT < 1:
                dp_m3_t[i] = 0
            else:
                dp_m3_t[i] = (6.0/ma.pi * (M1_m3_t[i]**(2.0)/(M0_m3_t[i]**(3.0/2) * M2_m3_t[i]**(1.0/2)))  * v_zero_0)**(1.0/3)
                
            if M0_m2_t[i]* NT < 1:
                dp_m2_t[i] = 0
            else:
                dp_m2_t[i] = (6.0/ma.pi * (M1_m2_t[i]**(2.0)/(M0_m2_t[i]**(3.0/2) * M2_m2_t[i]**(1.0/2)))  * v_zero_0)**(1.0/3)
            
            if M0_m1_t[i] * NT < 1:
                dp_m1_t[i] = 0
            else:
                dp_m1_t[i] = (6.0/ma.pi * (M1_m1_t[i]**(2.0)/(M0_m1_t[i]**(3.0/2) * M2_m1_t[i]**(1.0/2)))  * v_zero_0)**(1.0/3)
            
            if M0_0_t[i] * NT< 1:
                dp_0_t[i] = 0
            else:
                dp_0_t[i] = (6.0/ma.pi *  (M1_0_t[i]**(2.0)/(M0_0_t[i]**(3.0/2) * M2_0_t[i]**(1.0/2)))  * v_zero_0)**(1.0/3)
            
            if M0_p1_t[i] * NT< 1:
                dp_p1_t[i] = 0
            else:
                dp_p1_t[i] = (6.0/ma.pi * (M1_p1_t[i]**(2.0)/(M0_p1_t[i]**(3.0/2) * M2_p1_t[i]**(1.0/2)))  * v_zero_0)**(1.0/3)
                
            if M0_p2_t[i] * NT< 1:
                dp_p2_t[i] = 0
            else:
                dp_p2_t[i] = (6.0/ma.pi * (M1_p2_t[i]**(2.0)/(M0_p2_t[i]**(3.0/2) * M2_p2_t[i]**(1.0/2)))  * v_zero_0)**(1.0/3)
                
            if M0_p3_t[i] * NT< 1:
                dp_p3_t[i] = 0
            else:
                dp_p3_t[i] = (6.0/ma.pi * (M1_p3_t[i]**(2.0)/(M0_p3_t[i]**(3.0/2) * M2_p3_t[i]**(1.0/2)))  * v_zero_0)**(1.0/3)
                
            if M0_p4_t[i] * NT< 1:
                dp_p4_t[i] = 0
            else:
                dp_p4_t[i] = (6.0/ma.pi * (M1_p4_t[i]**(2.0)/(M0_p4_t[i]**(3.0/2) * M2_p4_t[i]**(1.0/2)))  * v_zero_0)**(1.0/3)


                
                
                
                
            if M0_m4_t[i] * NT < 1 or M0_m4_t[i]*M2_m4_t[i]/M1_m4_t[i]**2.0 < 1:
                sigmag_m4_t[i] = 1.0
            else:
                sigmag_m4_t[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0_m4_t[i]*M2_m4_t[i]/M1_m4_t[i]**2.0)))
            
            if M0_m3_t[i]* NT < 1 or M0_m3_t[i]*M2_m3_t[i]/M1_m3_t[i]**2.0 < 1:
                sigmag_m3_t[i]  = 1.0      
            else:
                sigmag_m3_t[i]  = ma.exp(ma.sqrt(1.0/9*ma.log(M0_m3_t[i]*M2_m3_t[i]/M1_m3_t[i]**2.0)))             
            
            if M0_m2_t[i]* NT < 1 or M0_m2_t[i]*M2_m2_t[i]/M1_m2_t[i]**2.0 < 1:
                sigmag_m2_t[i]  = 1.0
            else:
                sigmag_m2_t[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0_m2_t[i]*M2_m2_t[i]/M1_m2_t[i]**2.0)))
            
            if M0_m1_t[i] * NT < 1 or M0_m1_t[i]*M2_m1_t[i]/M1_m1_t[i]**2.0 < 1:
                sigmag_m1_t[i]  = 1.0
            else:
                sigmag_m1_t[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0_m1_t[i]*M2_m1_t[i]/M1_m1_t[i]**2.0)))
            
            if M0_0_t[i] * NT< 1 or M0_0_t[i]*M2_0_t[i]/M1_0_t[i]**2.0 < 1:
                sigmag_0_t[i]  = 1.0
            else:
                sigmag_0_t[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0_0_t[i]*M2_0_t[i]/M1_0_t[i]**2.0)))
            
            if M0_p1_t[i] * NT< 1 or M0_p1_t[i]*M2_p1_t[i]/M1_p1_t[i]**2.0 < 1:
                sigmag_p1_t[i]  = 1.0
            else:
                sigmag_p1_t[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0_p1_t[i]*M2_p1_t[i]/M1_p1_t[i]**2.0)))
                
            if M0_p2_t[i] * NT< 1 or M0_p2_t[i]*M2_p2_t[i]/M1_p2_t[i]**2.0 < 1:
                sigmag_p2_t[i]  = 1.0
            else:
                sigmag_p2_t[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0_p2_t[i]*M2_p2_t[i]/M1_p2_t[i]**2.0)))
                
            if M0_p3_t[i] * NT< 1 or M0_p3_t[i]*M2_p3_t[i]/M1_p3_t[i]**2.0 < 1:
                sigmag_p3_t[i]  = 1.0
            else:
                sigmag_p3_t[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0_p3_t[i]*M2_p3_t[i]/M1_p3_t[i]**2.0)))
                
            if M0_p4_t[i] * NT< 1 or M0_p4_t[i]*M2_p4_t[i]/M1_p4_t[i]**2.0 < 1:
                sigmag_p4_t[i]  = 1.0
            else:
                sigmag_p4_t[i] = ma.exp(ma.sqrt(1.0/9*ma.log(M0_p4_t[i]*M2_p4_t[i]/M1_p4_t[i]**2.0)))
            
            
            temp_dp =  NT*(dp_m4_t[i]*M0_m4_t[i] + dp_m3_t[i]*M0_m3_t[i]+ dp_m2_t[i]*M0_m2_t[i]\
                    + dp_m1_t[i]*M0_m1_t[i] + dp_0_t[i]*M0_0_t[i] + dp_p1_t[i]*M0_p1_t[i] + dp_p2_t[i]*M0_p2_t[i]\
                     + dp_p3_t[i]*M0_p3_t[i]  + dp_p4_t[i]*M0_p4_t[i]  )

            temp_sigmag =  NT*(sigmag_m4_t[i]*M0_m4_t[i] + sigmag_m3_t[i]*M0_m3_t[i]+ sigmag_m2_t[i]*M0_m2_t[i]\
                    + sigmag_m1_t[i]*M0_m1_t[i] + sigmag_0_t[i]*M0_0_t[i] + sigmag_p1_t[i]*M0_p1_t[i] + sigmag_p2_t[i]*M0_p2_t[i]\
                     + sigmag_p3_t[i]*M0_p3_t[i]  + sigmag_p4_t[i]*M0_p4_t[i]  )
            
            M0_tot_t[i] =  NT* (M0_m4_t[i] + M0_m3_t[i] + M0_m2_t[i] + M0_m1_t[i] + M0_0_t[i] + M0_p1_t[i] \
                          + M0_p2_t[i] + M0_p3_t[i]  + M0_p4_t[i]) 
            dp_avg_t[i] = temp_dp/M0_tot_t[i]
            sigmag_avg_t[i] = temp_sigmag/M0_tot_t[i]

            average_charge[i] = NT*(-4.0*M0_m4_t[i] + -3.0*M0_m3_t[i] + -2.0*M0_m2_t[i] + -1.0*M0_m1_t[i] + 0.0*M0_0_t[i] + 1.0*M0_p1_t[i] \
                          + 2.0*M0_p2_t[i] + 3.0*M0_p3_t[i]  + 4.0*M0_p4_t[i] )/M0_tot_t[i]
#             #std_dp = np.std([dp_m4_t[i],dp_m3_t[i],dp_m2_t[i],dp_m1_t[i],dp_0_t[i],dp_p1_t[i],dp_p2_t[i],dp_p3_t[i],dp_p4_t[i] ])
#             #std_qq = np.std([])
#             #print N_tot_t[i]
            print >> f, t[i]*tau, M0_tot_t[i], dp_avg_t[i], sigmag_avg_t[i], average_charge[i]
    
    # Remove the brackets from temp file [ and ]
    with open('temp', 'r') as my_file:
        text = my_file.read()
        text = text.replace("[", "")
        text = text.replace("]", "")     
   
    # Use proper name for proper sub-case for CASE A
    file2 = 'bi_ion_' + str(int(ion)) +   '_par_' + str(int(par)) + '_' + str(int(T)) + '_' + file_add +'.txt'
    with open(file2, 'w') as my_file:
        my_file.write(text)
    
    # Remove the temp file from the folder
    os.remove('temp')



if __name__ == '__main__':
    main()
#     t = 0.01
#     y = Initial

# Checking with the monodisperse simplest model from SDH
#     print lam
#     N_simple = np.zeros((10))
#     t_simple = [0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.1, 0.5, 1.0]
#     for i in range(0,9):
#         N_simple[i] = N_0_0/(1.0 + KK*N_0_0*t_simple[i]/2)
#     print t_simple, N_simple
#     ans = rhs(t,Initial)
#     print ans


#     M0 = y[0:2*qm + 1]
#     M1 = y[2*qm + 1:2*(2*qm +1)]
#     M2 = y[2*(2*qm +1):3*(2*qm +1)]
#     print M1[qm+0]
#     print mom(M0,M1,M2,0,1.0)
#     print tau, NT


############################################ TRY_MOMENT_1 SUMMARY ###################################################
# Included:
#     1. Raw data printing in file
#     2. rhs is working fine (coagulation part, charging part and ion balance part)
#     3. A and B matrices function running fine
#     4. moment function working fine
#     5. Dimensional form
    
# Still needs to be done:
#     1. run the code, ode solver
#     2. include the post processing section

########################################## TRY_MOMENT_2 SUMMARY ###################################################
# Included:
# Non-dimensional
# fixed a lot of errors in the non-dimensional form
# post processing not doen in safe.

######################################### TRY_MOMENT_#5 SUMMARY ####################################################
# Checked the charging part and it is working fine.
# Correct non-dimensionalisation of coagulation part is done. This is done by doing three different non-dimensionalisation of the coefficients.
# The coagulation part is checked for neutral particles.
# Sum of dm1 for all the charges is found to be 0 for coagulation part as well as charging part.
# Code is working fine. Compares well with the unimodal code especially for ion concentration changing, initial ion concentration
# = 10^16 and initial particle concentration = 10^19.

######################################### TRY_MOMENT_#5 SUMMARY ####################################################
# 1. included sigmag part for b2 (coag_coeff_2) 
# 2. Works for unipolar as well as bipolar case


  
  

    

    
    
    
    
    
    
    


# In[ ]:




