#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 13:51:28 2018

@author: vanbeers
"""

import numpy as np

###################### KNOWN PARAMETERS ################################

mu_e = 398600.441 #[km^3/s^2]
au = 149597870.88 #[km]


###################### CARTESIAN TO KEPLERIAN #########################

#assign array for cartesian coordinates
keplerian = []

#input variables

x = 10157768.1264                                   #[m]                 
y = -6475997.0091                                  #[m]
z = 2421205.9518                                    #[m]
v_x = 1099.2953996                                 #[m/s]
v_y = 3455.1059240                                #[m/s]
v_z = 4355.0978095                                #[m/s]

#input variables in vector and scalar forms in m and m/s

r_v = np.array([x,                                #vector form [m]
                y,
                z])
r_v = r_v/1000                                    #vector form [km]
r_s = (np.sqrt(x**2 + y**2 + z**2))/1000          #scalar form [km]
r_hat = r_v/r_s                                   #unit vector form [km]

V_v = np.array([v_x,                              #vector form [m/s]
                v_y,
                v_z])
V_s = np.sqrt(v_x**2 + v_y**2 + v_z**2)           #scalar form [m/s]

#calculate angular momentum 
h_v = np.cross(r_v,V_v/1000)                      #vector form [km^2/s]
h_s = np.sqrt(h_v[0]**2 + h_v[1]**2 + h_v[2]**2)  #scalar form [km^2/s]

#calculate a in km
a = (1/((2/r_s)-((V_s/1000)**2/mu_e)))*1000
keplerian.append(a)
print "a =", a

#calculate e
e_v = (np.cross(V_v/1000,h_v)/mu_e) - (r_v/r_s)   #vector form
e_s = np.sqrt(e_v[0]**2 + e_v[1]**2 + e_v[2]**2)  #scalar form
e_hat = e_v/e_s                                   #unit vector form
keplerian.append(e_s)                                  
print "e =", e_s

#calculate i in radians
i = np.arccos(h_v[2]/h_s)
keplerian.append(i)
print "i =", i* (180/np.pi)

#calculate quardrant parameters to check the sign
nvector = np.array([0,
                    0,
                    1])
N_v = np.cross(nvector,h_v)                            #vector form
N_xy = np.sqrt(N_v[0]**2 + N_v[1]**2)                  #vector in xy plane
N_hat = N_v/N_xy                                       #unit vector form

#calculate RAAN in radians
omega = np.arctan2((N_v[1]/N_xy),(N_v[0]/N_xy))
if omega < 0:
    omega = (2*np.pi)+omega
keplerian.append(omega)
print "omega = ", omega*(180/np.pi)

#calculate argument of perigee in radians
if np.dot((np.cross(N_hat,e_hat)),h_v) > 0:
    w = np.arccos(np.dot(e_hat,N_hat))
    keplerian.append(w)
else:
    w = -np.arccos(np.dot(e_hat,N_hat))
    keplerian.append(w)
print "w = ", w*(180/np.pi)

#calculate true anomaly in radians
if np.dot((np.cross(e_v,r_v)),h_v) > 0:
    theta = np.arccos(np.dot(r_hat,e_hat))
    keplerian.append(theta)
else:
    theta = -np.arccos(np.dot(r_hat,e_hat))
    keplerian.append(theta)
if theta < 0:
    theta = (2*np.pi)+theta
print "theta = ", theta*(180/np.pi)

#calculate eccentric anomaly in radians
E = 2*np.arctan(np.tan(theta/2)/(np.sqrt((1+e_s)/(1-e_s))))
keplerian.append(E)
if E < 0:
    E = (2*np.pi)+E
print "E = ", E * (180/np.pi)

#calculate mean anomaly
M = E - e_s*np.sin(E)
if M < 0:
    M = (2*np.pi)+M
keplerian.append(M)
print "M = ", M * (180/np.pi)








