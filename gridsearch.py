#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 11:12:39 2018

@author: vanbeers
"""


import matplotlib.pyplot as plt
import numpy as np
import time 

##################### GRID SEARCH SIMULATION #############################


########## DEFINE FUNCTIONS

def himmelblau(x,y):
    f = (x**2 + y - 11.)**2 + (x + y**2 - 7.)**2
    return f

def surface(x,y):
    f = 2.*x**2 + 2.*y**2
    return f

def grid(ss):                           #get ss number of locations for function 
    xy = []
    side = int(np.sqrt(ss))-1           #number of gridlines per side
    jump = 5./side                      #distance between gridlines
    gridlines = np.arange(0.,5.+jump,jump) 
    for i in gridlines:
        for j in gridlines:
            xy.append([i,j])
    return xy
    


######### GET START TIME

start_time = time.time()



########## RUN GRID SEARCH

ss = 50.       #sample size

answers = []

xy = grid(ss)  #get pairs

for i in range(len(xy)): 
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

xy = np.array(xy)

plt.figure()
plt.scatter(xy[:,0],xy[:,1],c=answers,s=70)
plt.grid(which='major',axis='both')
plt.xlabel('x1')
plt.ylabel('x2')
plt.colorbar()
plt.rcParams.update({'font.size': 22})
plt.show()








