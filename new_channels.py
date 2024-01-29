import numpy as np 
class Kv_x :
    def __init__(self,
                 potential : float , 
                 reversal_potential :float,
                 conductance : float ,
                 T : float ,
                 animal : str
                 ) :
        self.potential = potential
        self.reversal_potential = reversal_potential
        self.gates = {} 
        self.conductance = conductance
        self.T = T #temp in C
        self.animal = animal 
    def add_gate(self, gate_name : str, 
                 gate_constants_inf : np.array,
                 gate_constants_tau : np.array,
                 gate_power : int,
                 inf_function ,
                 tau_function 
                  ) :
        self.gates[gate_name] = {
            "constants_inf" :gate_constants_inf,
            "constants_tau" : gate_constants_tau,
            "_inf" : inf_function(self.potential, gate_constants_inf  ), 
            "_tau" : tau_function(self.potential , gate_constants_tau ),
            "_pow" : gate_power,
            "inf_func" : inf_function ,
            "tau_func" : tau_function
            }
        
        
        
        
def inf_Kv1__12( potential ,gate_constants  ) :
    gate_inf = gate_constants[0] / ( 1 + np.exp( (potential - (-gate_constants[1] ) ) / gate_constants[2])    )
#1.0000/(1+ exp((v - -30.5000)/-11.3943))
    return gate_inf
def tau_Kv1__12(potential, gate_constants):
#1.0000/(1+ exp((v - -30.0000)/27.3943)) 
    gate_tau = gate_constants[0] / ( 1 + np.exp( (potential - (-gate_constants[1] ) ) / gate_constants[2])    )
    return gate_tau


V =30 #random voltage : 
## Kv_1_1
'''
Animal	rat
CellType	Oocyte
Age	0 Days
Temperature	24.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[271] J P Adelman et. al; Science 1989 Apr 14
mpower	1.0
m Inf	1.0000/(1+ exp((v - -30.5000)/-11.3943))
m Tau	30.0000/(1+ exp((v - -76.5600)/26.1479))
hpower	2.0
h Inf	1.0000/(1+ exp((v - -30.0000)/27.3943))
h Tau	15000.0000/(1+ exp((v - -160.5600)/-100.0000))
'''
#https://channelpedia.epfl.ch/wikipages/1
Kv_1_1 = Kv_x(V, -65, 0.00001, 24, "rat")
Kv_1_1.add_gate("m", 
                [1.0, 30.5, -11.39343], 
                [30.0, 76.56, 26.1479], 
                1,
                inf_Kv1__12, 
                tau_Kv1__12)
Kv_1_1.add_gate("h", 
                [1.0, 30.0, 27.393],
                [15000.0000, 160.5600, -100.0000], 
                2, 
                inf_Kv1__12, 
                tau_Kv1__12)


## Kv_1_2
'''
Animal	rat
CellType	Oocyte
Age	0 Days
Temperature	20.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[272] L K Sprunger et. al; Eur. J. Pharmacol. 1996 Oct 31
mpower	1.0
m Inf	1.0000/(1+ exp((v - -21.0000)/-11.3943))
m Tau	150.0000/(1+ exp((v - -67.5600)/34.1479))
hpower	1.0
h Inf	1.0000/(1+ exp((v - -22.0000)/11.3943))
h Tau	15000.0000/(1+ exp((v - -46.5600)/-44.1479))
'''
#https://channelpedia.epfl.ch/wikipages/2
Kv_1_2 = Kv_x(V, -65, 0.00001, 20, "rat")
Kv_1_2.add_gate("m", 
                [1.0000, 21.0000, -11.3943],
                [150.0000, 67.5600, 34.1479],
                1,
                inf_Kv1__12, 
                tau_Kv1__12
                )
Kv_1_2.add_gate("h", 
                [1.0000, 22.0000, 11.3943], 
                [15000.0000, 46.5600, -44.1479],
                1,
                inf_Kv1__12, 
                tau_Kv1__12
                )


#Kv_1_3 need new fuinctions 
'''
Animal	rat
CellType	Oocyte
Age	0 Days
Temperature	0.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[273] J P Adelman et. al; J. Immunol. 1990 Jun 15
mpower	1.0
m Inf	1.0000/(1+ exp((v - -14.1000)/-10.3000))
m Tau	(-0.2840 * v) + 19.1600 If v lt 50
m Tau	5 If v gteq 50
hpower	1.0
h Inf	1.0000/(1+ exp((v - -33.0000)/3.7000))
h Tau	(-13.7600 * v) + 1162.4000 If v lt 80
h Tau	60 If v gteq 80
'''
#TODO update this inf 
def inf_Kv1__123( potential ,gate_constants  ) :
    gate_inf = gate_constants[0] / ( 1 + np.exp( (potential - (-gate_constants[1] ) ) / gate_constants[2])    )
    return gate_inf
def tau_Kv1__3(potential, gate_constants):
#(-13.7600 * v) + 1162.4000 If v lt 80
    if potential >= gate_constants[2] :
        return gate_constants[3]
    gate_tau = (potential * gate_constants[0]) + gate_constants[1]
    return gate_tau

Kv_1_3 = Kv_x(V, -65, 0.00001, 0, "rat")
Kv_1_3.add_gate("m", 
                [1.0, 14.1, -10.3],
                [-0.2840, 19.1600, 50., 5],
                1,
                inf_Kv1__123, 
                tau_Kv1__3
                )
Kv_1_3.add_gate("h", 
                [1.0, 33.0, 3.7], 
                [-13.7600, 1162.4000, 80,60],
                1,
                inf_Kv1__123, 
                tau_Kv1__3
                )


##Kv1_4
'''
Animal	rat
CellType	Oocyte
Age	0 Days
Temperature	0.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[274] O Pongs et. al; EMBO J. 1989 Nov
mpower	1.0
m Inf	1.0000/(1+ exp((v - -21.7000)/-16.9000))
m Tau	3.0000
hpower	1.0
h Inf	1.0000/(1+ exp((v - -73.6000)/12.8000))
h Tau	119.0000
'''
def inf_Kv1__1234( potential ,gate_constants  ) :
    gate_inf = gate_constants[0] / ( 1 + np.exp( (potential - (-gate_constants[1] ) ) / gate_constants[2])    )
    return gate_inf
def tau_Kv1__4(potential, gate_constants):
    gate_tau = gate_constants[0]
    return gate_tau
Kv_1_4 = Kv_x(V, -65, 0.00001, 0, "rat")
Kv_1_4.add_gate("m", 
                [1.0, 21.7, -16.9],
                [3.],
                1,
                inf_Kv1__1234, 
                tau_Kv1__3
                )
Kv_1_4.add_gate("h", 
                [1.0, 73.6, 12.8], 
                [119.],
                1,
                inf_Kv1__1234, 
                tau_Kv1__4
                )

