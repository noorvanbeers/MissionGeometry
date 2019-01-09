#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:47:09 2018

@author: vanbeers
"""

import random
import matplotlib.pyplot as plt
import numpy as np
import time

##################### MONTE CARLO SIMULATION #############################


########## DEFINE HIMMELBLAU FUNCTION

def himmelblau(x,y):
    f = (x**2 + y - 11)**2 + (x + y**2 - 7)**2
    return f

def surface(x,y):
    f = 2.*x**2 + 2.*y**2
    return f
    
######### GET START TIME

start_time = time.time()

######### GENERATE RANDOM NUMBERS

ss = 500. #sample size
xy = []


for i in np.arange(ss):
    x = random.uniform(0.,5.)
    y = random.uniform(0.,5.)
    xy.append([x ,y])
    
    
########## RUN MONTE CARLO

answers = []

for i in np.arange(int(ss)):
    x = xy[i][0]
    y = xy[i][1]
    f = himmelblau(x,y)
    answers.append(f)


########## GET OPTIMUM
xy = np.array(xy)

optimum = min(answers)
location = answers.index(optimum)

print "optimum =", optimum
print "x =", xy[location,0], "y =", xy[location,1]

######### GET RUN TIME

run_time = time.time()-start_time
print "run-time is", run_time, "s"



######### PLOTTING
    

plt.figure()
plt.scatter(xy[:,0],xy[:,1],c=answers,s=70)
plt.grid(which='major',axis='both')
plt.xlabel('x1')
plt.ylabel('x2')
plt.colorbar()
plt.rcParams.update({'font.size': 22})
plt.show()















