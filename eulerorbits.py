#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 16:27:18 2018

@author: vanbeers
"""

import numpy as np
import time
import matplotlib.pyplot as plt
from car2kep import cart2kep

######## START TIME

start_time = time.time()


####### DEFINE PARAMETERS AND CONSTANTS FROM WERTZ

mu_e = 3.98600441*10**14.               #gravitational parameters, [m^3/s^2]
a = 7500000.                             #semi major axis [m]
e = 0.1                                 #eccentricity
i = 0                                   #inclination [rad]

####### DEFINE INITIAL CONDITIONS (at periapsis intersection with x axis)

r_p = a*(1.-e)                            #radius of periapsis [m]
v_y = (mu_e*((2./r_p)-(1./a)))**0.5         #velocity [m/s]

y0_o = [r_p,0.,0.,0.,v_y,0.]                      #initial conditions, x, y, z, xdot, ydot, zdot
y0_c = [0., 0.]

######## DEFINE FUNCTIONS TO INTEGRATE

def orbit(t,y0):      
   r = np.sqrt(y0[0]**2. + y0[1]**2. + y0[2]**2.)
   dVdtx = -(mu_e/r**3.)*y0[0]
   dVdty = -(mu_e/r**3.)*y0[1]
   dVdtz = -(mu_e/r**3.)*y0[2]
   return np.array([y0[3], y0[4], y0[5], dVdtx, dVdty, dVdtz])

                 
def carvel(t,y0):      
    dVdt = 2.   #acceleration is 2 m/s^2, does not depend on t
    return np.array([y0[1], dVdt])
 


####### EULER INTEGRATOR


def euler(f,initial,t_i,t_e,h):
    y0 = initial
    t = t_i
    tl = [t]
    plot = np.array([y0])
    while t <= t_e:  
        y0 = y0 + h * f(t,y0)  #integrate acceleration and velocity
        t = t + h 
        tl.append(t)
        plot = np.append(plot, np.array([y0]), axis=0)
        
    plt.figure()
    plt.plot(plot[:,0],plot[:,1])
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.grid(which='major',axis='both')
    plt.show()
    return y0


# f = function to integrate in ODE form (derivative of y)
# initial = initial conditions, predefined
# t_i = initial time [s]
# t_e = end time  [s]
# h = time step [s]

######## GET RESULTS

results = euler(orbit,y0_o,0.,604800.,20.)
    
results = results.tolist()
               
transformed = cart2kep(results)                  #get kepler elements 

print "a =",transformed[0],"e =", transformed[1]   #print a and e
               
########## GET RUNTIME

print "runtime =", (time.time()-start_time)  #get runtime

      
########### ERRORS change return plot to return y0
#errors = []
#evaluations = []
#for i in [10.,5.,1.,0.5,0.1,0.05,0.01,0.005,0.001,0.0005]:  #for step size x
#    errors.append(abs((euler(carvel,y0_c,0,60,i))[0]-3600.))  #get error
#    evaluations.append(60./i)                               #get # evals
#
#plt.figure()
#plt.scatter(errors,evaluations)
#plt.xlabel('error [m]')
#plt.ylabel('step size [s]')
#plt.yscale('log')
#plt.grid(which='major',axis='both')
#plt.show()


    
#            
### PLOTTING, change return y0 to return plot
#
#
#
#plt.figure()
#plt.plot((euler(orbit,y0,0,6000,1.)[:,0]),(euler(orbit,y0,0,6000,1.)[:,1]),label='h=10')
#plt.plot((euler((orbit,y0,0,60,1.)[:,1]),label='h=1'))
#plt.plot(np.arange(0,60.1,0.1),(euler(orbit,y0,0,60,1.)[:,0]),label='h=0.1')
#plt.plot(np.arange(0,60.02,0.01),(euler(orbit,y0,0,60,1.)[:,0]),label='h=0.01')
#plt.plot(np.arange(0,60.002,0.001),(euler(orbit,y0,0,60,1.)[:,0]),label='h=0.001')
#plt.xlabel('x [m]')
#plt.ylabel('y [m]')
#plt.grid(which='major',axis='both')
#plt.legend()
#plt.xlim((0,85))
#plt.show()



    
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           
           

#def carvel(t):      
#    dVdt = 2.   #acceleration is 2 m/s^2, does not depend on t
#    return dVdt
#            
#
######## EULER INTEGRATOR, FOR TWO ODE
#
#
#def euler(f,y10,y20,t_i,t_e,h):
#    (t,y1) = (t_i,y10)
#    (t,y2) = (t_i,y20)
#    while t < t_e:  
#        y1 = y1 + h * y2  #get distance from velocity
#        y2 = y2 + h * f(t)
#        t = t + h
#
#    
#    return y1, y2
#
#print euler(carvel,0,0,0,60,0.01)  #0 initial velocity, 0-60 seconds, 1 s time step
#
#           
#
#           
           
           
           
           
           
           
           
           
           