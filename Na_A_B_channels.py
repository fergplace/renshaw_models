import numpy as np 
from expm1 import expm1
##could do a class for each channel, but might do a dict 
#

#dict:


class channel:
    def __init__(self, 
                alpha_rate :float, #ms^(-1)
                alpha_slope : float, #mV
                alpha_half : float , #mV
                
                beta_rate :float,
                beta_slope : float,
                beta_half : float ,
                
                reversal_potential : float, #mV
                conductance         : float, #S/cm^2
                
                channel_powers : np.array , 
                
                tau : float ,
                shift:float  =0.0 #mV
                
                ):
        
        self.alpha_rate = alpha_rate
        self.alpha_slope = alpha_slope
        self.alpha_half = alpha_half
        
    
        self.beta_rate = beta_rate
        self.beta_slope = beta_slope
        self.beta_half = beta_half
        
        self.reversal_potential = reversal_potential
        
        self.conductance = conductance
        
        self.channel_powers = channel_powers
        self.tau = tau 
        self.shift = shift 
        
####################THINK ABOUT ADDING THESE INTO CHANNEL CLASS
# Forward Dynamics alpha_m_K:
#3.1
def calc_alpha_m_K(V: float, sh: float, K_m_alpha_half:float ,K_m_alpha_rate :float , K_m_alpha_slope:float ) -> float:
    temp1 = K_m_alpha_half + sh - V 
    alpha_m_K = K_m_alpha_rate * expm1(temp1, K_m_alpha_slope)
    # alpha_m_K = K_m_alpha_rate* (temp1  / (np.exp( temp1 / K_m_alpha_slope) -1) )
    return alpha_m_K

# Reverse Dynamics Beta_m_K:
#3.2
def calc_beta_m_K(V: float, sh: float, K_m_beta_half:float ,K_m_beta_rate :float , K_m_beta_slope:float) -> float:
    temp1 = (K_m_beta_half + sh - V) / K_m_beta_slope 
    beta_m_K = K_m_beta_rate * np.exp(temp1)
    return beta_m_K

def calc_alpha_m_Na(V: float, sh: float, Na_m_alpha_half : float , Na_m_alpha_slope: float, Na_m_alpha_rate :float ) -> float:
    temp1 = Na_m_alpha_half + sh - V 
    alpha_m_Na = Na_m_alpha_rate * expm1(temp1, Na_m_alpha_slope)
    # alpha_m_Na = Na_m_alpha_rate * (temp1/  ( np.exp(temp1/Na_m_alpha_slope)- 1  ) )
    return alpha_m_Na
    
# ---------Reverse Dynamics beta_m_Na: #
#eq. 2.1.1
def calc_beta_m_Na(V: float, sh: float,  Na_m_beta_half : float , Na_m_beta_slope: float, Na_m_beta_rate :float) -> float:
    temp1 = V - Na_m_beta_half - sh
    beta_m_Na = Na_m_beta_rate * expm1(temp1, Na_m_beta_slope)
    # beta_m_Na = Na_m_beta_rate* (  temp1 / ( np.exp( temp1/Na_m_beta_slope)-1) ) 
    return beta_m_Na


# Forward Dynamics alpha_h_Na: #
#eq. 2.2.1
def calc_alpha_h_Na(V: float, sh: float,  Na_h_alpha_half : float , Na_h_alpha_slope: float, Na_h_alpha_rate :float) -> float:
    temp1 = (Na_h_alpha_half + sh - V) / Na_h_alpha_slope
    alpha_h_Na = Na_h_alpha_rate * np.exp(temp1)
    return alpha_h_Na

# Reverse Dynamics beta_h_Na:
#eq. 2.2.1
def calc_beta_h_Na(V: float, sh: float, Na_h_beta_half : float , Na_h_beta_slope: float, Na_h_beta_rate :float) -> float:
    temp1 = (Na_h_beta_half + sh - V) / Na_h_beta_slope
    beta_h_Na = Na_h_beta_rate / (np.exp(temp1) + 1.)
    return beta_h_Na
