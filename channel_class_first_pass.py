import numpy as np 

##Kv_x class 
class Kv_x :
    def __init__(self,
                 potential : float , 
                 reversal_potential :float) :
        self.potential = potential
        self.reversal_potential = reversal_potential
        self.gates = {} 
    
    def add_gate(self, gate_name : str, 
                 gate_constants_inf : np.array,
                 gate_constants_tau : np.array,
                 gate_power : int ) :
        self.gates[gate_name] = {
            "constants_inf" :gate_constants_inf,
            "constants_tau" : gate_constants_tau,
            "_inf" : inf(gate_constants_inf,  self.potential), 
            "_tau" : tau(gate_constants_tau,  self.potential),
            "_pow" : gate_power
            }


def inf( gate_constants, potential ) :
#1.0000/(1+ exp((v - -30.5000)/-11.3943))
    return

def tau(gate, potential):
#1.0000/(1+ exp((v - -30.0000)/27.3943)) 
    return 
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
#T=24 C
Kv_1_1 = Kv_x(V, -65)
#1.0000/(1+ exp((v - -30.5000)/-11.3943))
Kv_1_1.add_gate("m", [1.0, 30.5, 11.39343], [30.0, 76.56, 26.1479], 1)
#1.0000/(1+ exp((v - -30.0000)/27.3943)) 
Kv_1_1.add_gate("h", [1.0, 30.0, 27.393], [15000.0000, 160.5600, 100.0000], 2)





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
Kv_1_2 = Kv_x(V, -65)
Kv_1_2.add_gate("m", [1.0000, 21.0000, -11.3943], [150.0000, 67.5600, 34.1479], 1)
Kv_1_2.add_gate("h", [1.0000, 22.0000, 11.3943], [15000.0000, 46.5600, -44.1479], 1)
































### TODO add class for this type : 
## Kv_1_3 same type as 1.5
'''
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
Kv_1_2 = Kv_x(V, -65)
Kv_1_2.add_gate("m", [1.0, 14.1, -10.3], [150.0, 67.56, 34.1479], 1)
Kv_1_2.add_gate("h", [1.0000, 22.0000, 11.3943], [15000.0000, 46.5600, -44.1479], 1)


#Kv1.4 TODO add class for this 
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

#1.6 is frog 

#2.1 and 2.2 frog 

#3.1 unknown animal eq. similar to 


##3.2 DC added to tau gate:  TODO maybe class for this. 
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


##3.3
#frog but same format as our class  
'''
Model Kv3.3 (ID=28)      
Animal	Xenopus
CellType	oocyte
Age	0 Days
Temperature	23.0°C
Reversal	-65.0 mV
Ion	K +
Ligand ion	
Reference	[280] A J Rashid et. al; J. Neurosci. 2001 Jan 1
mpower	1.0
m Inf	1/(1+exp(((v -(18.700))/(-9.700))))
m Tau	20.000/(1+exp(((v -(-46.560))/(-44.140))))
'''


##4
#4.2 , 4.1 no HH model on pedia, looking at paper might be some info
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
#clipping here as well



##7

##7.x NA eq on first pass 



##NAV: mod files look like this :https://modeldb.science/230137?tab=2&file=SinKin_ModelDB/Nav11_a.mod

#1.3
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
#clipping class and v dependent numerator 

#1.6
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
#different eq. 
##TODO add equtions as an input into the channels (could work)


#CAV2.2 N
'''
https://channelpedia.epfl.ch/icdata/MOD/6.mod
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
#alpha beta def channel here 

#CAV3.1 T 
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
#clipping, DC, 

'''
Model Cav3.3 (ID=42)      
Animal	CH
CellType	CHO
Age	0 Days
Temperature	0.0°C
Reversal	30.0 mV
Ion	Ca +
Ligand ion	
Reference	[103] Achraf Traboulsie et. al; J. Physiol. (Lond.) 2007 Jan 1
mpower	1.0
m Inf	1/(1+exp((v- -45.454426)/-5.073015))
m Tau	3.394938 +( 54.187616 / (1 + exp((v - -40.040397)/4.110392)))
hpower	1.0
h Inf	1 /(1+exp((v-(-74.031965))/8.416382))
h Tau	109.701136 + (0.003816 * exp(-v/4.781719))
'''


#KCA: SK
#no ez HH model


#SLo1 : (BK) no HH model
