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
    #change to have alpha beta per gate... 
    def add_gate(self, gate_name : str, 
                 gate_constants_inf : np.array =None,
                 gate_constants_tau : np.array =None,
                 gate_power : int =1 ,
                 inf_function =None,
                 tau_function =None,
                 alpha_parms : np.array =None,
                 beta_parms : np.array=None,
                 alpha_function=None,
                 beta_function =None
                  ) :
        if gate_name in  self.gates.keys() :
             self.gates[gate_name] = {} #clear gate if trying to add to it again 

        tmp_dict= {
            "constants_inf" :gate_constants_inf,
            "constants_tau" : gate_constants_tau,
            "_pow" : gate_power,
            "inf_func" : inf_function ,
            "tau_func" : tau_function,
            "alpha_func" : alpha_function,
            "beta_func" : beta_function,
            "alpha_parms" :alpha_parms,
            "beta_parms" : beta_parms 
            }
        self.gates[gate_name]  = {k : v for k,v  in tmp_dict.items() if v is not None }
        
        #if inf_constants + tau constants 
        if "constants_inf" in self.gates[gate_name] :
            self.gates[gate_name]["_inf" ] = inf_function(self.potential, gate_constants_inf)
        if "constants_tau" in self.gates[gate_name] :
            self.gates[gate_name]["_tau" ] = tau_function(self.potential, gate_constants_inf)
        # if alpha and beta 
        if "alpha_parms" and "beta_parms" in self.gates[gate_name] :
            alpha = self.gates[gate_name]["alpha_func"](self.potential, self.gates[gate_name]["alpha_parms"])
            beta =  self.gates[gate_name]["beta_func"](self.potential, self.gates[gate_name]["beta_parms"])
            self.gates[gate_name]["_tau" ] = calc_tau_gate(alpha, beta)
            self.gates[gate_name]["_inf" ] = calc_gate_inf(alpha, beta)
        

def calc_gate_inf(alpha_gate: np.array, beta_gate: np.array) -> np.array:
    gate_inf = alpha_gate / (alpha_gate + beta_gate)
    return gate_inf

#eq. 2.1.4, 2.2.4 , 3.4 , 4.2.6
def calc_tau_gate(alpha_gate: np.array, beta_gate: np.array) -> np.array: 
    tau = 1. / (alpha_gate + beta_gate) 
    return tau  
        
def inf_Kv__default( potential ,gate_constants  ) :
    gate_inf = gate_constants[0] / ( 1 + np.exp( (potential - (gate_constants[1] ) ) / gate_constants[2])    )
#1.0000/(1+ exp((v - -30.5000)/-11.3943))
    return gate_inf
def tau_Kv__default(potential, gate_constants):
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
                [1.0, -30.5, -11.39343], 
                [30.0, -76.56, 26.1479], 
                1,
                inf_Kv__default, 
                tau_Kv__default)
Kv_1_1.add_gate("h", 
                [1.0, -30.0, 27.393],
                [15000.0000, -160.5600, -100.0000], 
                2, 
                inf_Kv__default, 
                tau_Kv__default)


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
                [1.0000, -21.0000, -11.3943],
                [150.0000, -67.5600, 34.1479],
                1,
                inf_Kv__default, 
                tau_Kv__default
                )
Kv_1_2.add_gate("h", 
                [1.0000, -22.0000, 11.3943], 
                [15000.0000, -46.5600, -44.1479],
                1,
                inf_Kv__default, 
                tau_Kv__default
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


def tau_Kv__clipping(potential, gate_constants):
#(-13.7600 * v) + 1162.4000 If v lt 80
    if potential >= gate_constants[2] :
        return gate_constants[3]
    gate_tau = (potential * gate_constants[0]) + gate_constants[1]
    return gate_tau

Kv_1_3 = Kv_x(V, -65, 0.00001, 0, "rat")
Kv_1_3.add_gate("m", 
                [1.0, -14.1, -10.3],
                [-0.2840, 19.1600, 50., 5],
                1,
                inf_Kv__default, 
                tau_Kv__clipping
                )
Kv_1_3.add_gate("h", 
                [1.0, -33.0, 3.7], 
                [-13.7600, 1162.4000, 80,60],
                1,
                inf_Kv__default, 
                tau_Kv__clipping
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

def tau_Kv__const(potential, gate_constants):
    gate_tau = gate_constants[0]
    return gate_tau
Kv_1_4 = Kv_x(V, -65, 0.00001, 0, "rat")
Kv_1_4.add_gate("m", 
                [1.0, -21.7, -16.9],
                [3.],
                1,
                inf_Kv__default, 
                tau_Kv__const
                )
Kv_1_4.add_gate("h", 
                [1.0, -73.6, 12.8], 
                [119.],
                1,
                inf_Kv__default, 
                tau_Kv__const
                )

##Kv1_5
'''
Animal	Human
CellType	Oocyte
Age	0 Days
Temperature	9.0°C
Reversal	-90.0 mV
Ion	K +
Ligand ion	
Reference	[275] L H Philipson et. al; Proc. Natl. Acad. Sci. U.S.A. 1991 Jan 1
mpower	1.0
m Inf	1.0000/(1+ exp((v - -6.0000)/-6.4000))
m Tau	(-0.1163 * v) + 8.3300 If v lt 50
m Tau	2 If v gteq 50
hpower	1.0
h Inf	1.0000/(1+ exp((v - -25.3000)/3.5000))
h Tau	(-15.5000 * v) + 1620.0000 If v lt 100
h Tau	50 If v gteq 100
'''


Kv_1_5 = Kv_x(V, -90, 0.00001, 9, "human")
Kv_1_5.add_gate("m", 
                [1.0, -6.0, -6.4],
                [-0.1163, 8.3300, 50., 2.],
                1,
                inf_Kv__default, 
                tau_Kv__clipping
                )
Kv_1_5.add_gate("h", 
                [1.0, -25.3, 3.5], 
                [-15.5000, 1620.0000, 100. , 50.],
                1,
                inf_Kv__default, 
                tau_Kv__clipping)

#Kv1_6
'''
Animal	Xenopus
CellType	oocyte
Age	0 Days
Temperature	23.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[276] O Pongs et. al; EMBO J. 1990 Jun
mpower	1.0
m Inf	1/(1+exp(((v -(-20.800))/(-8.100))))
m Tau	30.000/(1+exp(((v -(-46.560))/(44.140))))
hpower	1.0
h Inf	1/(1+exp(((v -(-22.000))/(11.390))))
h Tau	5000.000/(1+exp(((v -(-46.560))/(-44.140))))
'''
Kv_1_6 = Kv_x(V, -90, 0.00001, 9, "human")
Kv_1_6.add_gate("m", 
                [1.0, -20.800, -8.100],
                [30.000, -46.560, 44.140],
                1,
                inf_Kv__default, 
                tau_Kv__default
                )
Kv_1_6.add_gate("h", 
                [1.0, -22.000, 11.390], 
                [5000.000, -46.560, -44.140],
                1,
                inf_Kv__default, 
                tau_Kv__default)

#Kv1_7 , 1_8 NA


#Kv2.1
'''
Animal	Xenopus
CellType	oocyte
Age	0 Days
Temperature	23.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[277] A M VanDongen et. al; Neuron 1990 Oct
mpower	1.0
m Inf	1/(1+exp(((v -(-9.200))/(-6.600))))
m Tau	100.000/(1+exp(((v -(-46.560))/(44.140))))
hpower	1.0
h Inf	1/(1+exp(((v -(-19.000))/(5.000))))
h Tau	10000.000/(1+exp(((v -(-46.560))/(-44.140))))
'''
Kv_2_1 = Kv_x(V, -65., 0.00001, 23, "Xenopus")
Kv_2_1.add_gate("m", 
                [1.0,- 9.200, -6.600],
                [100.000, -46.560, 44.140],
                1,
                inf_Kv__default, 
                tau_Kv__default
                )
Kv_2_1.add_gate("h", 
                [1.0, -19.000, 5.0], 
                [10000.000, -46.560, -44.140],
                1,
                inf_Kv__default, 
                tau_Kv__default)



#Kv2.2
'''
Animal	Xenopus
CellType	oocyte
Age	0 Days
Temperature	23.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[278] S D Koh et. al; Am. J. Physiol. 1998 May
mpower	1.0
m Inf	1/(1+exp(((v -(5.000))/(-12.000))))
m Tau	130.000/(1+exp(((v -(-46.560))/(-44.140))))
hpower	1.0
h Inf	1/(1+exp(((v -(-16.300))/(4.800))))
h Tau	10000.000/(1+exp(((v -(-46.560))/(-44.140))))
'''

Kv_2_2 = Kv_x(V, -65., 0.00001, 23, "Xenopus")
Kv_2_2.add_gate("m", 
                [1.0, 5.00, -12.00],
                [130.000, -46.560, -44.140],
                1,
                inf_Kv__default, 
                tau_Kv__default
                )
Kv_2_2.add_gate("h", 
                [1.0, 16.300, 4.8], 
                [10000.000, 46.560, -44.140],
                1,
                inf_Kv__default, 
                tau_Kv__default)


#Kv3.1
'''
Animal	Xenopus
CellType	oocyte
Age	0 Days
Temperature	23.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[279] M Stocker et. al; EMBO J. 1992 Jul
mpower	1.0
m Inf	1/(1+exp(((v -(18.700))/(-9.700))))
m Tau	20.000/(1+exp(((v -(-46.560))/(-44.140))))
'''
Kv_3_1 = Kv_x(V, -65., 0.00001, 23, "Xenopus")
Kv_3_1.add_gate("m", 
                [1.0, 18.700, -9.700],
                [20., -46.560, -44.140],
                1,
                inf_Kv__default, 
                tau_Kv__default
                )

#kv3_2
'''
Animal	CH
CellType	CHO
Age	0 Days
Temperature	22.0°C
Reversal	-80.0 mV
Ion	K +
Ligand ion	
Reference	[18] R Hernandez-Pineda et. al; J. Neurophysiol. 1999 Sep
mpower	2.0
m Inf	1/(1+exp((v- -0.373267)/-8.568187))
m Tau	3.241643 +( 19.106496 / (1 + exp((v - 19.220623)/4.451533)))
'''

def tau_Kv_shift(potential , gate_constants):
    
    val = tau_Kv__default(potential, gate_constants[1:]) 
    val_shifted = val + gate_constants[0]
    return val_shifted

Kv_3_2 = Kv_x(V, -80.0, 0.00001, 22, "CH")
Kv_3_2.add_gate("m", 
                [1.0, -0.373267, -8.568187],
                [3.241643, 19.106496, 19.220623, 4.451533],
                2,
                inf_Kv__default, 
                tau_Kv_shift
                )

#Kv7_1:5 NA have to check in the papers


#Kv4.1 NA, 
#KV4.2
'''
Animal	rat
CellType	Neocortical L5PC
Age	14 Days
Temperature	23.0°C
Reversal	-68.7 mV
Ion	K +
Ligand ion	
Reference	[296] J M Bekkers et. al; J. Physiol. (Lond.) 2000 Jun 15
mpower	3.0
m Inf	(1/(1 + exp((v- -18.8)/-16.6)))
m Tau	1.0/((0.026* exp(-0.026*v)) + (35* exp(0.136*v))) If v lt -50
m Tau	1.7/(1+ exp((-42 - v)/-26)) + 0.34 If v gteq -50
hpower	1.0
h Inf	1/(1 + exp((v - -81.6)/6.7))
h Tau	0.01*v + 6.7
'''
def tau_Kv__4m( potential ,gate_constants  ) :
    if potential < gate_constants[0] :
        tmp1 = gate_constants[1]
        tmp2 = gate_constants[2] * np.exp( (gate_constants[3] - potential)  )
        tmp3 = gate_constants[4] *  np.exp(gate_constants[5] * potential)
        return tmp1/(tmp2 + tmp3)
    else : 
        tmp1 = gate_constants[6]
        tmp2 = 1 + np.exp( ( gate_constants[7] - potential)/gate_constants[8] ) + gate_constants[9]
        return tmp1/tmp2
def tau_Kv__4h(potential, gate_constants):
    return (gate_constants[0] * potential ) + gate_constants[1]
Kv_4_2 = Kv_x(V, -68.7, 0.00001, 23, "rat")
Kv_4_2.add_gate("m", 
                [1.0, -18.8, -16.6],
                [-50., 1., 0.026, -0.026, 35., 0.136, 1.7, -42., -26, 0.34],
                3,
                inf_Kv__default, 
                tau_Kv__4m
                )
Kv_4_2.add_gate("h", 
                [1.0, -81.6, 6.7],
                [0.01, 6.7],
                1,
                inf_Kv__default, 
                tau_Kv__4h
                )

#Kv_4_3 NA

#NAV: 1.1, 2 in neuron... 

#NAV 1.3:
'''
Animal	rat
CellType	Neocortical
Age	0 Days
Temperature	23.0°C
Reversal	50.0 mV
Ion	Na +
Ligand ion	
Reference	[43] T R Cummins et. al; J. Neurosci. 2001 Aug 15
mpower	3.0
m Alpha	(0.182 * ((v)- -26))/(1-(exp(-((v)- -26)/9))) If v neq -26
m Beta	(0.124 * (-(v) -26))/(1-(exp(-(-(v) -26)/9))) If v neq -26
hpower	1.0
h Inf	1 /(1+exp((v-(-65.0))/8.1))
h Tau	0.40 + (0.265 * exp(-v/9.47))
'''
#TODO make flag for alpha beta vs. minf and mtau , then change class to be for optinal calling..
def inf_Nav__default( potential ,gate_constants  ) :
    if potential == gate_constants[1] : 
        potential = potential + 0.000001 #from mod file 

    alpha, beta = alpha_beta_Nav(potential ,gate_constants )
    
    return 

def alpha_beta_Nav(potential ,gate_constants ) :
    tmp1 = (potential)- gate_constants[1]
    tmp2 = gate_constants[0] * tmp1
    tmp3 = 1- np.exp(-(tmp1/gate_constants[2] )) 
    alpha = tmp2/ tmp3
     
    tmp4 = (-potential)- gate_constants[1]
    tmp5 = gate_constants[3] * tmp4
    tmp6 = 1- np.exp(-(tmp4/gate_constants[2] )) 
    beta= tmp5/tmp6
    return alpha, beta 
 

[ 0.182, -26, 9., 0.124]
def tau_Nav__default(potential, gate_constants):
    
    return 

