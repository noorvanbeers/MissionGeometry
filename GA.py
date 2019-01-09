#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 10:43:33 2019

@author: vanbeers
"""

import random
import numpy as np
import time


#################### GENETIC ALGORITHM ####################


##### DEFINE NECESSARY FUNCTIONS

def himmelblau(x,y):              #execute himmelblau function
    f = (x**2. + y - 11.)**2. + (x + y**2. - 7.)**2.
    return f

def getbinary(denary):            #convert denary number to binary
    res = denary - 1
    binary= ""  
    while res > 0:   
        binary = str(res % 2) + binary  
        res = res//2  
    return int(binary)

def getdenary(binary):           #convert binary to denary number
    binary = str(binary)
    denary = 0  
    for digit in binary:  
        denary = denary*2 + int(digit)  
    return int(denary) 

def bitcount(res):               #count number of bits in binary number
    bitnum = len(str(res))
    return bitnum

def getdelta(UB1,UB2,LB1,LB2):
    d1 = (UB1-LB1)/((2**split)-1)
    d2 = (UB2-LB2)/((2**(bitnum-split)-1))
    return d1, d2

def getparent(bits):             #generate random parent from bitnumber
    parent = ""
    for i in range(bits):
        parent = parent + str(random.randint(0,1))
    return parent
    
def populate(resolution):        #get initial population
    population = []
    bits = bitcount(getbinary(resolution))
    for i in range(resolution):
        population.append(getparent(bits))
    return np.array(population)

def mutate(parent):              #mutate the digits in the parents 
    chance = 0.0001
    for digit in parent:
        selector = random.randint(0,1./chance)
        if selector == 1:
            if digit == 0:
                digit = 1
            else:
                digit = 0
    return parent
               
def getpairs(population):         #split population into parent pairs
    pairs = []
    i = 0
    for i in range(0,len(population),2):
        pairs.append([population[i],population[i+1]])
    return pairs
        

def getchildren(pairs):            #make children from parents by splitting
    children = []
    for pair in pairs:
        mother = pair[0]
        father = pair[1]
        child1 = str(mother[:split]) + str(father[split:])
        child2 = str(father[:split]) + str(mother[split:])
        children.append(child1)
        children.append(child2)
    return children


def getfittest(evaluations,population):           #get fittest (minimum)
    while len(evaluations) > lenevals/2:
           indx = evaluations.index(max(evaluations))   #find index of max
           population.remove(population[indx])          #remove from list
           evaluations.remove(max(evaluations))         #remove from list
    return evaluations, population



############ INITIATE RUN TIME AND GENERATION COUNT

start_time = time.time()
generation = 1



############# INITIATE POPULATION 


LB1 = 0.                           #upper and lower boundaries of function
UB1 = 5.
LB2 = 0.
UB2 = 5.


res = 32                           #state resolution needed, input a integer!


population = populate(res)         #initiate population

                     
bitnum = bitcount(population[0])   #get number of bits for future functions
split = bitnum/2                   #get where to split for future functions
                  
         
for parent in population:          #random mutations
    mutate(parent)

evaluations = []

for parent in population:           #evaluate himmelblau for parents
    d1,d2 = getdelta(UB1,UB2,LB1,LB2)
    x1 = LB1 + d1*(getdenary(parent[:split]))
    x2 = LB2 + d2*(getdenary(parent[split:]))
    f = himmelblau(x1,x2)    
    evaluations.append(f)

best1 = min(evaluations)


running = True

while running:
    
################# START GENERATING CHILDREN 

    pairs = getpairs(population)         #sort population into parent pairs
    
    children = getchildren(pairs)        #generate children
    
    population = list(population) + children   #add children to population
    
    for child in children:               #evaluate himmelblau for children
        child = mutate(child)            #mutations
        d1,d2 = getdelta(UB1,UB2,LB1,LB2)
        x1 = LB1 + d1*(getdenary(child[:split]))    #apply constraints
        x2 = LB2 + d2*(getdenary(child[split:]))
        f = himmelblau(x1,x2)
        evaluations.append(f)
        
    lenevals = len(evaluations)         #get number of evaluations for future
    
                  
                  
################## SORT THE FITTEST
    
    evaluations, population = getfittest(evaluations,population)
    
    best2 = min(evaluations)
    
    
################ EVALUATE STOP OR GO
    
    generation = generation + 1
    
    if (best1-best2)/best1 <= 0.0001:
    
        running = False
        
        varindex = evaluations.index(best2)
        solution = population[varindex]
        x1 = getdenary(solution[:split])
        x2 = getdenary(solution[split:])
    
    else:
        best1 = best2

    
    
################# RESULTS 

print "runtime is", time.time()- start_time
print best2, "at x1 =", x1*d1, "and x2 =", x2*d2
print "generations = ", generation



    







