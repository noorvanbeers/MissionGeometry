#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 22:57:47 2018

@author: vanbeers
"""

import numpy as np
import matplotlib.pyplot as plt

##########  constants   ##################################

pi = np.pi

G = 6.67259*10**-20             # gravitational constant [km^3 kg^-1 s^-2]
mu_e = 3.98600441*10**5         #mu earth [km^3 s^-2]
R_e = 6378.136              #radius of the earth [km]
J2 = 1082.62*10**-6             #J2 coefficient
T_e = 86164.1004            #length of a day [s]
k = 3                       #number of days before repeat

###########  define functions from Wertz #################################

k1 = mu_e**(1./3)*(2*pi/T_e)**(-2./3)
k2 = 1.029549648*10**14  #0.75*J2*(mu_e**0.5)*(R_e**2)
L_dot = 360.
e = 0.


def H0(j):
    H_0 = k1*(j/k)**(-2./3) #km
    return H_0

def getomegadot(a,i):  #deg per day
    omega_dot = -2.*k2*(a**(-7./2))*np.cos(i)*(1-e**2)**-2
    return omega_dot

def getwdot(a,i):  #deg per day
    w_dot = k2*(a**(-7./2))*((5*np.cos(i)**2)-1)*(1-e**2)**-2
    return w_dot

def getmdot(a,i):  #deg per day
    m_dot = k2*(a**(-7./2))*((3*np.cos(i)**2)-1)*(1-e**2)**-2
    return m_dot

def getn(j,omega_dot,w_dot,m_dot):   #deg per day
    n = (j/k)*(L_dot-omega_dot) - (w_dot + m_dot)
    return n
    
def Hrev(n):   #H revised in km
    nrads = (n/T_e)*(pi/180.)
    H = (mu_e/(nrads**2))**(1./3)
    return H

def iterate(pair):
    jj = pair[0]
    ii = pair[1]
    H = H0(jj)
    calculating = True
    while calculating:
        omega_dot = getomegadot(H,ii)
        w_dot = getwdot(H,ii)
        m_dot = getmdot(H,ii)
        n = getn(jj,omega_dot,w_dot,m_dot)   #deg per sidereal day
        Hit = Hrev(n)
        
        #accuracy
        if abs(Hit-H) > 0.00001:
            H = Hit
        else:
            return [Hit,ii,jj]
            calculating = False

####################### verification with slides
#
#
#k = 1
#e = 0
#pair = [14,(28.*pi/180)]
#jj = pair[0]
#ii = pair[1]
#H = H0(jj) 
#print H
#calculating = True
#while calculating:        
#    omega_dot = getomegadot(H,ii)
#    print omega_dot
#    w_dot = getwdot(H,ii)
#    print w_dot
#    m_dot = getmdot(H,ii)
#    print m_dot
#    n = getn(jj,omega_dot,w_dot,m_dot)
#    print n
#    nrads = (n/T_e)*(pi/180.)
#    print nrads
#    Hit = (mu_e/(nrads**2))**(1./3)
#    print Hit
#        
#    #need accuracy
#    if abs(Hit - H) > 0.0001:    
#        H = Hit
#    else:
#        print [Hit,ii]
#        calculating = False


###################### start

j = np.arange(39.,49.)     #number of orbits before repeat
alli = np.arange(0,pi,0.0174533)    #inclination range, rad


#get initial j + inclination pairs for iteration
pairs = []

for js in j:
    for i in alli:
        pairs.append([js,i])


#calculate the combinations for repeat orbits, append them if  between 200 and 1200km

repeats = []   #list of repeat orbit combinations

for pair in pairs:
    repeat = iterate(pair)
    if repeat[0]<(1200+R_e) and repeat[0]>(200+R_e):
        repeats.append(repeat)
    
#split data for plotting    
repeats = np.array(repeats)
   
splt = np.split(repeats,np.where(np.diff(repeats[:,2]))[0]+1)

    
#plot the data, yes its bad programming this way
       
plt.figure()
plt.plot((splt[0][:,0]-R_e),splt[0][:,1]*(180/pi),label='39,3')
plt.plot((splt[1][:,0]-R_e),splt[1][:,1]*(180/pi),label='40,3')
plt.plot((splt[2][:,0]-R_e),splt[2][:,1]*(180/pi),label='41,3')
plt.plot((splt[3][:,0]-R_e),splt[3][:,1]*(180/pi),label='42,3')
plt.plot((splt[4][:,0]-R_e),splt[4][:,1]*(180/pi),label='43,3')
plt.plot((splt[5][:,0]-R_e),splt[5][:,1]*(180/pi),label='44,3')
plt.plot((splt[6][:,0]-R_e),splt[6][:,1]*(180/pi),label='45,3')
plt.plot((splt[7][:,0]-R_e),splt[7][:,1]*(180/pi),label='46,3')
plt.plot((splt[8][:,0]-R_e),splt[8][:,1]*(180/pi),label='47,3')
plt.plot((splt[9][:,0]-R_e),splt[9][:,1]*(180/pi),label='48,3')
plt.xlabel('altitude [km]')
plt.ylabel('inclination [deg]')
plt.grid(which='major',axis='both')
plt.xlim(180,1500)
plt.legend(loc="upper right")
plt.show()

app2 = []
y = np.arange(0,len(repeats),50)

for i in y:
    app2.append(repeats[i])
    

np.savetxt('app2.txt', app2, fmt='%10f' ,header= "       a          i          j")

#
##### THIS ONE WORKS ######
#k = 1
#e = 0
#pair = [14,(28.*pi/180)]
#jj = pair[0]
#ii = pair[1]
#H = H0(jj) 
#calculating = True
#while calculating:        
#    omega_dot = getomegadot(H,ii)
#    print omega_dot
#    w_dot = getwdot(H,ii)
#    print w_dot
#    m_dot = getmdot(H,ii)
#    print m_dot
#    n = getn(jj,omega_dot,w_dot,m_dot)
#    print n
#    nrads = (n/T_e)*(pi/180.)
#    print nrads
#    Hit = (mu_e/(nrads**2))**(1./3)
#    print Hit
#        
#    #need accuracy
#    if abs(Hit - H) > 0.0001:    
#        H = Hit
#    else:
#        print [Hit,ii]
#        calculating = False


    
    
    
    
    

