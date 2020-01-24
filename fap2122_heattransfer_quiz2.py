#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 17:36:43 2018

@author: nandopas
"""

import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
import numpy as np
import scipy as sci
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#from pylab import figure, plot, xlabel, grid, hold, legend, title, savefig
from pylab import *
from scipy.optimize import curve_fit
import math

material = 'water'
Pressure_condensor = 8000
Pressure_condensor_out = Pressure_condensor
Inlet_quality = 0.9
Outlet_quality = 0.0
T_2_in = 293 #Temperature cooling water inlet
T_2_out = 306 #Temperature cooling water outlet
efficiency = .4
W_net = 900e6
Q_in = W_net/efficiency #Heat input to cycle as whole (Boiler)
Q_out = Q_in - W_net #Heat output from condenser
F = 1 #initial correction factor
fouling = .0003

#dimensions of tubes
L = 16.5 #meters
d_o = .022 
d_i = .01
n_passes = 2 #number of tube passes

'''
---notes to self---
W_net/Q_in = efficiency
so Q_in = W_net/efficiency
W_net = Q_in-Q_out
What i need:
    steam inlet temperature
        quality 0.9 and pressure 8000 Pa
    condensed steam outlet temperature
        quality 0.0 and pressure 8000 Pa
    Enthalpy values for all points
    
    -steady state energy equation
    
    0 = sum(mh)_inlets - sum(mh)_outlets
    0 = (m_1*h_1_in + m_2*h_2_in) - (m_1*h_1_out + m_2*h_2_out)
    
    -Heat exchanger heat transfer with LMTD
    Q_out = U(pi*D*L)*LMTD
    where LMTD = ((T_1_0-T_2_0)-(T_1_L-T_2_L))/(ln((T_1_0-T_2_0)/(T_1_L-T_2_L)))
    
'''

#Many initial dimensions based off of GE 1000 MW condensor
#https://www.gepower.com/steam/heat-exchange/condenser
#Calculation design based off of Kern method
#https://nptel.ac.in/courses/103103027/pdf/mod1.pdf
#other important sources of information:
#http://www.cit-wulkow.chemstations.com/content/documents/Technical_Articles/shell.pdf
#Chapter 3 of A Heat Transfer Textbook (Lienhardt)
#http://web.mit.edu/lienhard/www/ahttv211.pdf




#Inlet and outlet Temperatures
T_1_in = PropsSI('T','P', Pressure_condensor, 'Q', Inlet_quality, material)
T_1_out = PropsSI('T','P', Pressure_condensor_out, 'Q', Outlet_quality, material)

#mean density of flow
density_in = PropsSI('D','P', Pressure_condensor, 'Q', Inlet_quality, material)
density_out = PropsSI('D','P', Pressure_condensor_out, 'Q', Outlet_quality, material)
density_1 = (density_in + density_out)/2

#mean viscosity of flow
viscosity_in = PropsSI('V','P', Pressure_condensor, 'Q', Inlet_quality, material)
viscosity_out = PropsSI('V','P', Pressure_condensor_out, 'Q', Outlet_quality, material)
viscosity_1 = (viscosity_in + viscosity_out)/2

#LMTD
LMTD = ((T_1_in-T_2_out)-(T_1_out - T_2_in))/(math.log((T_1_in-T_2_out)/(T_1_out - T_2_in)))


#UA based on values we know
UA = Q_out/(LMTD*F)

#based off http://www.cit-wulkow.chemstations.com/content/documents/Technical_Articles/shell.pdf
#Typical values of U range from 141 - 264 (page 24)
#so use the average at first
U_mean = (800+1500)/2
A = UA/U_mean


#Enthalpies of steam and condensed water
h_1_in = PropsSI('H','P', Pressure_condensor, 'Q', Inlet_quality, material)
h_1_out = PropsSI('H','P', Pressure_condensor, 'Q', Outlet_quality, material)
C_2 = PropsSI('C', 'T', (T_2_in+T_2_out)/2, 'Q', Inlet_quality , material)

#Q = m(delta H)

#mass flow rate of flow 1 (tube)
m_1 = Q_out/(h_1_in-h_1_out)

#number of tubes
n_tubes = A/(math.pi*d_o*L)

#tube fluid velocity
u = (4*m_1*(n_passes/n_tubes))/(math.pi*density_1*(d_i**2))

#Reynolds number tube
Re_tube = u*d_i*density_1/viscosity_1

while Re_tube < 1e4:
    if Re_tube < 1e4:
        print("The Reynolds number within the tubes is less than 100000")
        print("fix the number of passes")
    elif Re_tube >= 1e4:
        print("The Reynolds number within the tubes is greater than 100000")
        print("calculate pressure")
    n_passes = n_passes + 2
    print("Number of passes", n_passes)
    u = (4*m_1*(n_passes/n_tubes))/(math.pi*density_1*(d_i**2))
    Re_tube = u*d_i*density_1/viscosity_1
    print("Reynolds Number", Re_tube)
    print("")

print("The Reynolds number within the tubes is greater than 100000")

#tube pressure drop
f = 64/Re_tube #Darcy friction
delta_p_tube = f*((density_1*(u**2))/2)


    
print("")
print("Reynolds Number", Re_tube)
print("velocity", u)
print("")

print("calculate pressure")
print("pressure drop", delta_p_tube)
print("now reiterate")


#------------------------------REITERATION-------------------------------------
Pressure_condensor_out = Pressure_condensor_out - delta_p_tube

#Inlet and outlet Temperatures
T_1_in = PropsSI('T','P', Pressure_condensor, 'Q', Inlet_quality, material)
T_1_out = PropsSI('T','P', Pressure_condensor_out, 'Q', Outlet_quality, material)

#mean density of flow
density_in = PropsSI('D','P', Pressure_condensor, 'Q', Inlet_quality, material)
density_out = PropsSI('D','P', Pressure_condensor_out, 'Q', Outlet_quality, material)
density_1 = (density_in + density_out)/2

#mean viscosity of flow
viscosity_in = PropsSI('V','P', Pressure_condensor, 'Q', Inlet_quality, material)
viscosity_out = PropsSI('V','P', Pressure_condensor_out, 'Q', Outlet_quality, material)
viscosity_1 = (viscosity_in + viscosity_out)/2

LMTD = ((T_1_in-T_2_out)-(T_1_out - T_2_in))/(math.log((T_1_in-T_2_out)/(T_1_out - T_2_in)))

#correction factor calculation:
P = (T_1_out-T_1_in)/(T_2_in-T_1_in)
R = (T_2_in-T_2_out)/(T_1_out-T_1_in)
F = .99
U = Q_out/(A*F*LMTD)

#Tube side heat transfer coefficient
k = 10 #conductivity Stainless Steel
C_in = PropsSI('C', 'P', Pressure_condensor, 'Q', Inlet_quality , material)
C_out = PropsSI('C', 'P', Pressure_condensor, 'Q', Outlet_quality , material)
C_1 = (C_in+C_out)/2

h_tube = (k*42.0)/((((viscosity_1 * C_1)/k)**(-.33))*d_o)


#stuff for shell
n_tubes = int(n_tubes)
pitch = .001
r_shell = 5
d_shell = 10
n = (r_shell*math.sqrt(math.pi))/(.022+pitch) #number of tubes per square side


viscosity_2 = PropsSI('V', 'T', (T_2_in+T_2_out)/2, 'Q', Inlet_quality , material)



#Mass flow rate cooling flow 2
m_2 = Q_out/(C_2*(T_2_out-T_2_in))

#h_shell = 2400.0
k = 10 #for stainless
#equivalent diameter for the shell
d_equivalent = abs(4*((pitch**2)-((math.pi *(d_o**2))/4))/(math.pi*d_o))

h_shell = (k*20)/((((viscosity_2 * C_2)/k)**(-.33))*d_equivalent)
#U_overall = 1/(((1/1) + fouling + (d_o**2/d_i**2)*(d_o-d_i/2*k) + (d_o**2/d_i**2)*(1/h_tube) + (d_o**2/d_i**2)*fouling))
U_overall = 1/((1/h_tube) + (1/h_shell))

percent = (U-U_overall)/U_overall




print("")
print("")
print("")
print("")
print("list of important values")
print("Efficiency assumed", efficiency*100, "%")
print("U assumed", U_mean, "W/m^2K")
print("Heat transfered", Q_out, "W")
print("Area of transfer", A, "m^2")
print("mass flow rate steam/ rate of condensation", m_1, "kg/s")
print("mass flow rate shell", m_2, "kg/s")
print("pressure drop tube", delta_p_tube,"Pa")
print("number of tubes", n_tubes)
print("tube fluid velocity", u, "m/s")
print("Reynolds number of tube", Re_tube)
print("Tube passes", n_passes)
print("Length per pass", L, "m")
print("Tube inner diameter",d_i,"m")
print("Tube outer diameter",d_o,"m")
print("Reiterated U", U,"W/m^2K")
