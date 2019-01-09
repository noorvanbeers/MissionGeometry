#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 19:29:20 2018

@author: vanbeers
"""


import matplotlib.pyplot as plt
import numpy as np
import time

##################### NEWTON-RAPHSON SIMULATION #############################


########## DEFINE FUNCTIONS

def himmelblau(x):
    f = (x**4 - 22.*x**2 + x + 114)
    return f

def poly(x):
    f = x**2 - 3*x - 20
    return f

def newtonraphson(x0,h,margin):
    f0 = himmelblau(x0)
    iterations = 0
    while abs(f0) >= margin:
        xnew = x0 - (himmelblau(x0)/((himmelblau(x0+h)-himmelblau(h))/h))
        f0 = himmelblau(xnew)
        x0 = xnew
        iterations = iterations + 1
    return (f0, x0, iterations)
        
    
#h = step size
#margin = tolerance margin
#x0 = initial guess 

######### GET START TIME

start_time = time.time()


######## RUN NEWTON RAPHSON

guess = 2.                #initial guess
step = 0.1                  #step size
margin = 0.01             #tolerance allowed


a = newtonraphson(guess,step,margin)
print "optimum =", a[0]
print "x_1 =", a[1]
print "x_2 =", (7.-a[1])**0.5 
print "number of iterations =", a[2]


######### GET RUN TIME

run_time = time.time()-start_time
print "run-time is", run_time, "s"



