#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 17:51:11 2018

@author: vanbeers
"""

import matplotlib.pyplot as plt
import numpy as np
import time
from car2kep import cart2kep


######## START TIME
start_time = time.time()


####### DEFINE PARAMETERS AND CONSTANTS FROM WERTZ

mu_e = 3.98600441*10**14                #gravitational parameters, [m^3/s^2]
a = 7500000.                           #semi major axis [m]
e = 0.1                                 #eccentricity
i = 0                                   #inclination [rad]


####### DEFINE INITIAL CONDITIONS (at periapsis intersection with x axis)

r_p = a*(1-e)                          #radius of periapsis [m]
v_y = (mu_e*((2./r_p)-(1./a)))**0.5       #velocity [m/s]

y0_o = [r_p,0.,0.,0.,v_y,0.]                  #x, y, z, xdot, ydot, zdot (orbit)
y0_c = [0., 0.]                           #y, ydot (car)


###### FUNCTIONS TO INTEGATE

def orbit(t,y0):      
   r = np.sqrt(y0[0]**2. + y0[1]**2. + y0[2]**2.)       #magnitude of r
   dVdtx = -(mu_e/r**3.)*y0[0]                           #acc in x/y/z
   dVdty = -(mu_e/r**3.)*y0[1]
   dVdtz = -(mu_e/r**3.)*y0[2]
   return np.array([y0[3], y0[4], y0[5], dVdtx, dVdty, dVdtz])

                         
def car(t,y0):          
    dVdt = 2.                       #acceleration is 2 m/s^2,
    return np.array([y0[1], dVdt])  #velocity, acceleration


####### RUNGE-KUTTA 4

def rk4(f,initial,t_i,t_e,h):
    y0 = initial                #initial conditions
    t = t_i                     #initial time
    tl = [t]                    #plot times
    plot = np.array([y0])       #plot results
    while t < t_e:
        k1 = f(t,y0)
        k2 = f(t+(h*0.5),(y0+h*k1*0.5))
        k3 = f(t+(h*0.5),(y0+h*k2*0.5))
        k4 = f(t+h,y0+h*k3)
        sigma = (1./6.)*(k1 + 2.*k2 + 2.*k3 + k4)
        y0 = y0 + sigma * h
        t = t + h
        tl.append(t)
        plot = np.append(plot, np.array([y0]), axis=0)
    plt.plot(plot[:,0],plot[:,1])   #change plot according to integrated function
    plt.grid(which='major', axis='both')
    plt.xlabel('x [m]')             #change axis labels accordingly
    plt.ylabel('y [m]')
    plt.show()
    return y0


######## GET RESULTS change return plot to return y0

results = rk4(orbit,y0_o,0.,604800.,23.)


results = results.tolist()
               
transformed = cart2kep(results)                 #get kepler elements 

print "a =",transformed[0],"e =", transformed[1]  #print a and e


print transformed[0] - 7500000
######## GET RUN TIME

print "run time = ", (time.time()-start_time)



####### ERRORS change return plot to return y0


#
#errors = []
#evaluations = []
#for i in [150.,100.,50.,10.,5.]:  #for step size x
#    errors.append(abs((rk4(orbit,y0_o,0,604800.,i))[0]-7500000.))  #get error
#    evaluations.append(604800./i)                               #get # evals
#
#plt.figure()
#plt.scatter(errors,evaluations)
#plt.xlabel('error [m]')
#plt.ylabel('step size [s]')
#plt.yscale('log')
#plt.grid(which='major',axis='both')
#plt.show()
#
#print [errors, evaluations]



#### PLOTTING change return y0 to return plot
#
#
#plt.figure()
#plt.plot(np.arange(0,70,10.),(rk4(carvel,0,0,0,60,10.)[:,0]),label='h=10')
#plt.plot((rk4(carvel,0,0,0,60,1.)[:,0]),label='h=1')
#plt.plot(np.arange(0,60.1,0.1),(rk4(carvel,0,0,0,60,0.1)[:,0]),label='h=0.1')
#plt.plot(np.arange(0,60.02,0.01),(rk4(carvel,0,0,0,60,0.01)[:,0]),label='h=0.01')
#plt.plot(np.arange(0,60.002,0.001),(rk4(carvel,0,0,0,60,0.001)[:,0]),label='h=0.001')
#plt.xlabel('time [s]')
#plt.ylabel('distance [m]')
#plt.grid(which='major',axis='both')
#plt.legend()
#plt.xlim((0,85))
#plt.show()





