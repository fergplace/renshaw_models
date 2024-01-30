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
        #used to make a dict of only keys we want to use      
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
                gate_constants_inf= [1.0, -30.5, -11.39343], 
                gate_constants_tau= [30.0, -76.56, 26.1479], 
                gate_power= 1,
                inf_function= inf_Kv__default, 
                tau_function= tau_Kv__default)
Kv_1_1.add_gate("h", 
                gate_constants_inf= [1.0, -30.0, 27.393],
                gate_constants_tau= [15000.0000, -160.5600, -100.0000], 
                gate_power= 2, 
                inf_function= inf_Kv__default, 
                tau_function= tau_Kv__default)


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
                gate_constants_inf= [1.0000, -21.0000, -11.3943],
                gate_constants_tau= [150.0000, -67.5600, 34.1479],
                gate_power= 1,
                inf_function= inf_Kv__default, 
                tau_function= tau_Kv__default
                )
Kv_1_2.add_gate("h", 
                gate_constants_inf= [1.0000, -22.0000, 11.3943], 
                gate_constants_tau = [15000.0000, -46.5600, -44.1479],
                gate_power = 1,
                inf_function= inf_Kv__default, 
                tau_function =tau_Kv__default
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
                gate_constants_inf= [1.0, -14.1, -10.3],
                gate_constants_tau= [-0.2840, 19.1600, 50., 5],
                gate_power=1,
                inf_function=inf_Kv__default, 
                tau_function= tau_Kv__clipping
                )
Kv_1_3.add_gate("h", 
                gate_constants_inf= [1.0, -33.0, 3.7], 
                gate_constants_tau= [-13.7600, 1162.4000, 80,60],
                gate_power=1,
                inf_function=inf_Kv__default, 
                tau_function=  tau_Kv__clipping
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
                gate_constants_inf= [1.0, -21.7, -16.9],
                gate_constants_tau=[3.],
                gate_power=1,
                inf_function= inf_Kv__default, 
                tau_function= tau_Kv__const
                )
Kv_1_4.add_gate("h", 
                gate_constants_inf=[1.0, -73.6, 12.8], 
                gate_constants_tau=[119.],
                gate_power=1,
                inf_function= inf_Kv__default, 
                tau_function= tau_Kv__const
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
                gate_constants_inf= [1.0, -6.0, -6.4],
                gate_constants_tau=[-0.1163, 8.3300, 50., 2.],
                 gate_power=1,
                inf_function=inf_Kv__default, 
                tau_function=tau_Kv__clipping
                )
Kv_1_5.add_gate("h", 
                gate_constants_inf=[1.0, -25.3, 3.5], 
                gate_constants_tau=[-15.5000, 1620.0000, 100. , 50.],
                 gate_power=1,
                inf_function=inf_Kv__default, 
                tau_function=tau_Kv__clipping)

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
                gate_constants_inf=[1.0, -20.800, -8.100],
                gate_constants_tau= [30.000, -46.560, 44.140],
                gate_power=1,
                inf_function=inf_Kv__default, 
                tau_function=tau_Kv__default
                )
Kv_1_6.add_gate("h", 
                gate_constants_inf=[1.0, -22.000, 11.390], 
                gate_constants_tau= [5000.000, -46.560, -44.140],
                gate_power=1,
                inf_function=inf_Kv__default, 
                tau_function=tau_Kv__default)

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
                gate_constants_inf= [1.0,- 9.200, -6.600],
                gate_constants_tau=[100.000, -46.560, 44.140],
                gate_power=1,
                inf_function= inf_Kv__default, 
                tau_function= tau_Kv__default
                )
Kv_2_1.add_gate("h", 
                gate_constants_inf= [1.0, -19.000, 5.0], 
                gate_constants_tau=[10000.000, -46.560, -44.140],
                gate_power=1,
                inf_function= inf_Kv__default, 
                tau_function= tau_Kv__default)



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
                gate_constants_inf=[1.0, 5.00, -12.00],
                gate_constants_tau=[130.000, -46.560, -44.140],
                gate_power=1,
                inf_function=inf_Kv__default, 
                tau_function=tau_Kv__default
                )
Kv_2_2.add_gate("h", 
                gate_constants_inf=[1.0, 16.300, 4.8], 
                gate_constants_tau=[10000.000, 46.560, -44.140],
                gate_power=1,
                inf_function=inf_Kv__default, 
                tau_function=tau_Kv__default)


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
                gate_constants_inf=[1.0, 18.700, -9.700],
                gate_constants_tau= [20., -46.560, -44.140],
                gate_power= 1,
                inf_function=inf_Kv__default, 
                tau_function=tau_Kv__default
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
                gate_constants_inf=[1.0, -0.373267, -8.568187],
                gate_constants_tau=[3.241643, 19.106496, 19.220623, 4.451533],
                gate_power=2,
                inf_function=inf_Kv__default, 
                tau_function=tau_Kv_shift
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
                gate_constants_inf=[1.0, -18.8, -16.6],
                gate_constants_tau= [-50., 1., 0.026, -0.026, 35., 0.136, 1.7, -42., -26, 0.34],
                gate_power=3,
                inf_function= inf_Kv__default, 
                tau_function=tau_Kv__4m
                )
Kv_4_2.add_gate("h", 
                gate_constants_inf=[1.0, -81.6, 6.7],
                gate_constants_tau= [0.01, 6.7],
                gate_power=1,
                inf_function= inf_Kv__default, 
                tau_function=tau_Kv__4h
                )

#Kv_4_3 NA

#NAV: 1.1, 2 in neuron... 

#NAV 1.3: ######NOTE:reversal is >0 !!!
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
def alpha_Nav(potential ,gate_constants ) :
    if potential == gate_constants[1] : 
        potential = potential + 0.000001 #from mod file 
    tmp1 = (potential)- gate_constants[1]
    tmp2 = gate_constants[0] * tmp1
    tmp3 = 1- np.exp(-(tmp1/gate_constants[2] )) 
    alpha = tmp2/ tmp3
    return alpha

def beta_Nav(potential ,gate_constants ) :
    if potential == gate_constants[1] : 
        potential = potential + 0.000001 #from mod file 
    tmp4 = (-potential)- gate_constants[1]
    tmp5 = gate_constants[3] * tmp4
    tmp6 = 1- np.exp(-(tmp4/gate_constants[2] )) 
    beta= tmp5/tmp6
    return beta 
 

[ 0.182, -26, 9., 0.124]
def tau_Nav__default(potential, gate_constants):
    tau = gate_constants[0] + ( gate_constants[1] * np.exp((-potential )/ gate_constants[2]))
    return tau 

Nav_1_3 = Kv_x(V, 50, 0.00001, 23, "rat")
Nav_1_3.add_gate("m", 
                alpha_parms=[ 0.182, -26, 9., 0.124],
                beta_parms= [ 0.182, -26, 9., 0.124],
                gate_power=3,
                alpha_function= alpha_Nav, 
                beta_function=beta_Nav
                )
Nav_1_3.add_gate("h", 
                gate_constants_inf=[1. , -65.0, 8.1],
                gate_constants_tau= [0.40, 0.265, 9.47],
                gate_power=1,
                inf_function= inf_Kv__default, 
                tau_function=tau_Nav__default
                )

#Nav1.6
'''
Animal	rat
CellType	L5PC
Age	21 Days
Temperature	23.0°C
Reversal	50.0 mV
Ion	Na +
Ligand ion	
Reference	[288] A L Goldin et. al; J. Neurosci. 1998 Aug 15
mpower	1.0
m Inf	1.0000/(1+ exp(-0.03937*4.2*(v - -17.000)))
m Tau	1
'''
def Nav_1_6_inf( potential, gate_constants):
    inf = gate_constants[0] / (1 + np.exp(gate_constants[1] * gate_constants[2] * (potential - gate_constants[3])  ))
    return inf 

Nav_1_6 = Kv_x(V, 50, 0.00001, 23, "rat")
Nav_1_6.add_gate("m", 
                gate_constants_inf=[1.0,-0.03937, 4.2, -17. ],
                gate_constants_tau= [1.],
                gate_power=1,
                inf_function= Nav_1_6_inf, 
                tau_function=tau_Kv__const
                )

#Cav2_2
'''
Animal	Cat
CellType	RGC
Age	34 Days
Temperature	36.0°C
Reversal	135.0 mV
Ion	Ca +
Ligand ion	
Reference	[259] S J Huang et. al; Neuroscience 1998 Jul
mpower	2.0
m Alpha	(0.1*(v-20)/(1-exp(-(v-20)/10))) If v neq 20
m Beta	0.4*exp(-(v+25)/18)
hpower	1.0
h Alpha	0.01*exp(-(v+50)/10)
h Beta	0.1/(1+exp(-(v+17)/17))
'''
def cav_2_2_alpha_clip(potential, gate_constants) :
    if potential == gate_constants[0] : 
        potential = potential + 0.000001 #from mod 
    tmp1 = gate_constants[1] * (potential -gate_constants[0] )
    tmp2 = 1 - np.exp( - ( potential - gate_constants[0]) /gate_constants[2]  )
    return tmp1/tmp2

def cav_2_2_alpha_no_clip(potential, gate_constants) :
    alpha = gate_constants[0] * np.exp( - (potential + gate_constants[1] )/gate_constants[2])
    return alpha 

def cav_2_2_beta_m(potential, gate_constants) :
    beta = gate_constants[0] * np.exp( - (potential + gate_constants[1] )/gate_constants[2])
    return beta 

def cav_2_2_beta_h(potential, gate_constants) :
    beta = gate_constants[0] / (1 + np.exp( - (potential + gate_constants[1] )/gate_constants[2]))
    return beta 

Cav_2_2 = Kv_x(V, 135., 0.00001, 36, "Cat")
Cav_2_2.add_gate("m", 
                alpha_parms=[20., 0.1, 10. ],
                beta_parms= [ 0.4, 25., 18.],
                gate_power=2,
                alpha_function= cav_2_2_alpha_clip, 
                beta_function=cav_2_2_beta_m
                )
Cav_2_2.add_gate("h", 
                gate_constants_inf=[0.01 , 50., 10.],
                gate_constants_tau= [0.1, 17, 17],
                gate_power=1,
                inf_function= cav_2_2_alpha_no_clip, 
                tau_function=cav_2_2_beta_h
                )
#Cav3.1 
'''
Animal	CH
CellType	CHO
Age	0 Days
Temperature	0.0°C
Reversal	30.0 mV
Ion	Ca +
Ligand ion	
Reference	[103] Achraf Traboulsie et. al; J. Physiol. (Lond.) 2007 Jan 1
mpower	1.0
m Inf	1 /(1+exp((v-(-42.921064))/-5.163208))
m Tau	-0.855809 + (1.493527 * exp(-v/27.414182)) If v lt -10
m Tau	1.0 If v gteq -10
hpower	1.0
h Inf	1 /(1+exp((v-(-72.907420))/4.575763))
h Tau	9.987873 + (0.002883 * exp(-v/5.598574))
'''
def tau_Cav_3_1_clip(potential, gate_constants) :
    if potential < -10 :
        tau = tau_Cav_3_1_no_clip(potential, gate_constants[1:-1]) 
    else :
        tau =gate_constants[4]
    return tau 

def tau_Cav_3_1_no_clip(potential, gate_constants) :
    tau = gate_constants[0] + (gate_constants[1] * np.exp((-potential)/gate_constants[2]) )
    return tau 
    
Cav_3_1 = Kv_x(V, 30, 0.00001, 0, "CH")
Cav_3_1.add_gate("m", 
                gate_constants_inf=[1.0,-42.921064, -5.163208],
                gate_constants_tau= [-10, -0.855809, 1.493527, 27.414182, 1],
                gate_power=1,
                inf_function= inf_Kv__default, 
                tau_function=tau_Cav_3_1_clip
                )
Cav_3_1.add_gate("h", 
                gate_constants_inf=[1.0, -72.907420, 4.575763],
                gate_constants_tau= [9.987873, 0.002883, 5.598574],
                gate_power=1,
                inf_function= inf_Kv__default, 
                tau_function=tau_Cav_3_1_no_clip
                )
