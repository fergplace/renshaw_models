import scipy.integrate as si
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from dataclasses import dataclass
"""
https://github.com/spines-center/gpu_ode_playground    
"""
import new_channels
import RenshawParams as rp
#from general_current_eq import general_current as gen_I
#import general_gate_eq  as gen_gate
#import default_Na_channel as default_Na
#import Potassium_AHP_channel as Ca_ch
#import K_nM_channel as K_nM_ch
#import delayed_rectifier_I_K as default_K

##eq. 1
def current_membrane(I_app :float, Cm :float, I_all : np.array)-> float:    
    """
    dV_dt is the rate of change of potential membrane with respect to time
    Args:
        I_Na (float):  Current sodium depolarization of membrance
        I_K (float):   Current Potassium
        I_AHP (float): Current After-Hyperpolarization membrane after action potential and 
                       contributation of refractory period.
        I_L (float):   Current Leak

    Returns:
        float: _description_
    """
    # mV/ms = 1/nF * (nA)
    dV_dt = (1. / Cm) * (I_app - (np.sum(I_all)  )) 
    return dV_dt


def renshaw_model(t: float, Y: np.ndarray, I_app_fn: float, delay_ms: float, 
                max_I_mA: float, duration_ms: float, Cm :float 
                ,channels ,desired_channels_name 
                ,channel_conduct, num_gates, E_L, soma_area_cm2
                ) -> np.ndarray:

    
    now_V = potential = Y[0]
    #potential = np.array([potential])
    # now_m_Na    = Y[1]
    # now_h_Na    = Y[2]
    # now_m_K     = Y[3]
    # now_m_AHP   = Y[4]
    # now_m_K_nM  = Y[5]
    # now_Ca2     = Y[6]
    # now_m_Ca    = Y[7]
    
    # now_I_Na    = Y[8]
    # now_I_K     = Y[9]
    # now_I_L     = Y[10]
    # now_I_K_nM  = Y[11]
    # now_I_AHP   = Y[12]

    num_crrents = len(desired_channels_name)+1 #add one for leak 
    #take input num channels and num gates : 
    all_I = np.zeros(num_crrents) 
    new_Y= np.zeros(num_crrents + num_gates +1)  #add one for V_new
    
    idx=1 
    for ch_num , channel in enumerate(desired_channels_name) :
        channel_in_use =  channels[channel]
        
        num_gates = channel_in_use.num_gates
        new_Y[idx:( idx+ num_gates) ] = \
            channel_in_use.calc_dgate_dt(potential, Y[idx: idx+ num_gates ])
            
        all_I[ch_num] = channel_in_use.current_channel_calc(
                                Y[idx: idx+ num_gates]
                                ,potential, channel_conduct[ch_num]
                                ) 
           
        idx = idx+ num_gates #update 
    #deal with the leak.      
    all_I[-1] = channel_conduct[-1]*(potential - E_L)*1e6 *soma_area_cm2
    #TODO fix the Y
    #dgate_dt_all = gen_gate.calc_dgate_dt(inf_all, tau_all, Y[1:6])
    
    I_app = I_app_fn(t, delay_ms, max_I_mA, duration_ms)
    
    #sum all currents 
    #currents_all = np.sum(all_I)

    dV_dt = current_membrane( I_app, Cm, all_I )
    

    #want to return new_Y
    new_Y[0] = dV_dt
    new_Y[idx:] = Y[idx:] - all_I
    return (new_Y)



def square_I_pulse(t: np.ndarray, delay_ms: float, max_I_nA: float, duration_ms: float) -> float:
    """
    Returns I_app value depending on the time t.
    :param t: current time in ms
    :param delay_ms: time delay in ms
    :param max_I_nA: current applied (nA)
    :param duration_ms: time duration of applied current in ms
    """
    return (max_I_nA *
            (np.heaviside(t - delay_ms, 1.) - 
             np.heaviside(t - (delay_ms + duration_ms), 1.) ) )


#only have I_app , duration, G

#move all G and make none default. 
##TODO add new channel in inputs here
def r_m_solver(
    square_I_pulse : float
    ,total_time_ms : float 
    ,delay_ms : float 
    ,I_app : float 
    ,duration_ms : float 
    ,channels : dict
    ,desired_channels_name 
    ,channel_conduct
    ,num_gates
    ,V0 : float  = rp.RenshawParams.Y0[0]
    ,Y0 : np.ndarray = rp.RenshawParams.Y0 
    ,soma_area_cm2 : float = rp.RenshawParams.soma_area_cm2
    ,Cm : float = rp.RenshawParams.Cm 
    ,E_L : float = rp.RenshawParams.E_L
   
    
            
    ) : 

    'positional args square_I_pulse, total_time,  delay_ms, I_app , duration_ms'
    t = np.linspace(0, total_time_ms, 10000)
    Y0[0] = V0
    sol = si.solve_ivp(
        renshaw_model, (0, t[-1]), Y0, t_eval = t, 
                       args=(
                           square_I_pulse, delay_ms, I_app, duration_ms
                            , Cm
                           ,channels ,desired_channels_name 
                           ,channel_conduct, num_gates,
                           E_L, soma_area_cm2
                           ), 
        method="BDF"
        )
    return sol 



    
def main()-> None:
    time_start = time.time()
    
    ### Applied Neuron Current "Injected electrode/generator" synaptic inputs
    I_app1 = 300e-3 # 350e-3  # nA 
    delay_ms1 = 30.
    duration_ms1 = 450.
    total_time_ms1 = 600
     # first(initial time), end time, number of step to reach the end time
 
    
    channels = new_channels.main()
    #list the channels we want, can see all options in new_channels.py
    desired_channels_name = ["Kv_1_2", "Kv_3_1", "Nav_1_6"]
    channel_conduct = np.array([200e-3,300e-3,10e-3 , 3.3e-4])
    
    num_gates  = 0
    for ch_name in desired_channels_name:
        num_gates = channels[ch_name].num_gates +num_gates 
    Y0 = np.zeros( len(desired_channels_name) + num_gates +2  )
    
    #TODO generate a Y0 
    #TODO fix new_channels.main to only create dict of channels we want. 
    sol = r_m_solver(
        square_I_pulse, total_time_ms = total_time_ms1, Y0=Y0,
        delay_ms = delay_ms1  , I_app = I_app1  , duration_ms = duration_ms1 
        ,channels = channels, desired_channels_name = desired_channels_name,
        channel_conduct =channel_conduct, num_gates = num_gates) 

    
    
    #print( get_mean_data_sweep_temporal_freq_adj( sol.t  , sol.y[0]  ,  -10 , duration_ms1   ) )


    time_end = time.time()
    print(f"Completed solution in {time_end - time_start:.2} seconds")
    plt.figure()
    # now_m_Na    = Y[1]
    # now_h_Na    = Y[2]
    # now_m_K     = Y[3]
    # now_m_AHP   = Y[4]
    # now_m_K_nM  = Y[5]
    # now_Ca2     = Y[6]
    # now_m_Ca    = Y[7] 
    #plt.plot(sol.t, sol.y[3, :], label='m_K')
    plt.plot(sol.t, sol.y[1, :], label='m_Kv_1_2') 
    plt.plot(sol.t, sol.y[2, :], label='h_Kv_1_2')
    plt.plot(sol.t, sol.y[3, :], label='m_Kv_3_1') 
    plt.plot(sol.t, sol.y[4, :], label='h_Kv_3_1')
    plt.plot(sol.t, sol.y[5, :], label='m_Nav_1_6') 
    plt.xlim(0, 200)
    plt.ylim(-1,1)
    plt.legend() 
    plt.title("newer")


    
    plt.figure()                             # creates a new graph with plot
    plt.plot(sol.t, sol.y[0, :])
    plt.xlabel('Time step (ms)')
    plt.ylabel('Membrane Voltage (mV)')
    #plt.xlim(0, total_time_ms1)
    #plt.xlim(25,45)
    plt.ylim(-80., 100.)
    plt.legend()                             # to label the slopes
    plt.show()                               # to show second graph
    
    
    return

if __name__ == "__main__":
    main()