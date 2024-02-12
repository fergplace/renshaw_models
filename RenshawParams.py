import scipy.integrate as si
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from dataclasses import dataclass

@dataclass 
class RenshawParams :
    ' class for storing default data '
    
########## Single Compartment Renshaw (Randy Power Model) ##########
    soma_L_um       : float         = 12. #Jankowska 10-15 
    soma_diam_um    : float         = 12. 
    soma_area_um2   : float         = np.pi * soma_diam_um * soma_L_um
    soma_area_cm2   : float         = soma_area_um2 * 1e-8  # area is 4.524e-4 cm2
# the goal would be to do mV * uS = nA

    E_L             : float         =   -70.0 # mV
    g_L             : float         = 0.00033  # S/cm2 
# Memmbrane Potential Voltage(V) 
    V0              : float         =  E_L     # mV was at -115, but changed to -65 (unsure where -115 came from)
    # Neuron Membrane Capacitance 
    Cm              : float         = 1. * soma_area_cm2 * 1e3  # uF/cm2 * cm2 * 1e3 = nF

    ## Sodium Variable Table: ##      
    # Global Values
    g_Na            : float         =     0.03   # S/cm^2, g_Na is for Sodium conductance, note this is per-unit area so if we want a bigger cell add here 
    E_Na            : float         =    55 #50.0    # mV, E_Na is for Reversal potential of sodium
    sh_Na           : float         =   -2.0

    # Gate m
    Na_m_alpha_rate : float         =     0.8    # Activation rate Na_m_α_rate ms^(-1), aka am_rate
    Na_m_alpha_half  : float        =   -47.5    # V for ½ activation Na_(m,α,half) mV, aka am_half
    Na_m_alpha_slope : float        =     4.0    # Activation slope Na_(m,α,slope) mV, aka am_slp
    Na_m_beta_rate  : float         =     0.7    # Deactivation rate Na_(m,β,rate) ms^(-1), aka bm_rate
    Na_m_beta_half  : float         =   -20.0    # V for ½ deactivation Na_(m,β,half) mV, aka bm_half
    Na_m_beta_slope : float         =     5.0    # Deactivation slope Na_(m,β,slope) mV, aka bm_slp
    # Gate h
    Na_h_alpha_rate : float         =     0.32   # Activation rate Na_(h,α,rate) ms^(-1), aka ah_rate
    Na_h_alpha_half : float         =   -35.0    # V for ½ activation Na_(h,α,half) mV, aka ah_half
    Na_h_alpha_slope: float         =    18.0    # Activation slope Na_(h,α,slope) mV, aka ah_slp
    Na_h_beta_rate  : float         =    10.0    # Deactivation rate ms^(-1), aka bh_rate
    Na_h_beta_half  : float         =   -20.0    # V for ½ deactivation mV, aka bh_half
    Na_h_beta_slope : float         =     5.0    # Deactivation slope mV, aka bh_slp

    ## Potassium Variable Table: ##
            # Global Parameter
    g_K             : float         =     0.01   # S/cm^2, Potassium Conductance
    E_K             : float         =   -80.0    # mV, Reversal potential of potassium
    sh_K            : float         =  -2.0  # mV, shift of potassium, K_traub.sh
            # Gate m
    K_m_alpha_rate  : float         =     0.03   # Activation rate ms^(-1), aka an_rate
    K_m_alpha_half  : float         =   -45.0    # V for half-activation mV, aka an_half
    K_m_alpha_slope : float         =     5.0    # Activation slope mV, aka an_slp
    K_m_beta_rate   : float         =     0.5    # Deactivation rate ms^(-1), aka bn_rate
    K_m_beta_half   : float         =   -50.0    # V for half-activation mV, aka bn_half
    K_m_beta_slope  : float         =    40.0    # Deactivation slope mV, aka bn_slp

    ## AHP (Afterhyperpolarization) gCa Channel Variable Table: ##
    g_Ca            : float         =   1e-5   # S/cm^2, Calcium Conductance
    E_Ca            : float         =   132.458    # mV, Reversal potential of Calcium
    # Gate m
    Ca_m_alpha_half : float         =   -20.0    # mV, aka mvhalfca
    Ca_m_alpha_slope: float         =     4.0    # mV/s, aka mslpca
    Ca_m_alpha_rate : float         =     1.0    # ms, aka mtauca
    Ca_m_beta_half  : float         =   -20.0    # mV
    Ca_m_beta_slope : float         =     4.0    # mV/s
    Ca_m_beta_rate  : float         =     1.0    # ms
    tau_m_Ca        : float         =   Ca_m_alpha_rate    # given by eq., aka mtauca

    ## AHP (Afterhyperpolarization) gKCa Channel Variable Table: ##
    g_KCa           : float         =     0.002  # S/cm^2 Calcium mediated potassium conductance
    b_K_Ca          : float         =     0.1    # Maximum activation rate, aka bKCa
    f_K_Ca          : float         =     0.1    # Maximum deactivation rate (something something uM), aka fKCa
    Ca2_depth       : float         =     0.1    # Calcium concentration diffusion depth, aka depth
    Ca2_inf         : float         =     1e-4   # Calcium maximum concentration limit mM, aka cainf
    F               : float         = 96485.3415 # Faraday Constant, aka FARADAY
    beta_m_AHP      : float         = b_K_Ca           # given by eq., aka bKCa

    Ca2_tau         : float         = 20.                 # rate of calcium removal, aka taur

    m_Na0           : float         = 0.0  # calc_gate_inf( calc_alpha_m_Na(V0), calc_beta_m_Na(V0))
    h_Na0           : float         = 0.0  # calc_gate_inf( calc_alpha_h_Na(V0), calc_beta_h_Na(V0))
    m_K0            : float         = 0.0
    Ca20            : float         = Ca2_inf 
    m_AHP0          : float         = 0.0  # calc_gate_inf( calc_alpha_m_AHP(Ca20), beta_m_AHP)
    m_Ca0           : float         = 0.0  # calc_gate_inf( calc_alpha_m_Ca(V0), calc_beta_m_Ca(V0))

    
    #K_nM
    ca              : float         = 0.036 #ms^-1
    RT              : float         = 8.31 * 308 #R gas constatn and T temp in K (35C)
    cb              : float         = 0.002 #ms^-1
    zeta_a          : float         = 0.909#mV
    zeta_b          : float         = 1.102 #mV 
    g_bar_K_nM      : float         =.0001  #mho/cm^2 
    #for exponents in current eq. in order of Na, I_K, I_Ca, I_AHP, I_L, I_K_nM
    channel_powers = np.array([ [3,1], [4,1] , [1,1] , [1,1] , [1,1] , [1,1] ])
    # Y format is: V, m_Na, h_Na, m_k, m_Ca, Ca2, m_AHP, I_app, I_Na, I_K, I_L , K_nM, I_K_nM
    Y0              : np.ndarray    = np.array((V0, m_Na0, h_Na0, m_K0, m_Ca0, Ca20, m_AHP0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))