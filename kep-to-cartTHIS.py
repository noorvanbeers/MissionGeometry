#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 13:50:19 2018

@author: vanbeers
"""


import numpy as np

###################### KNOWN PARAMETERS ################################

mu_e = 398600.441 #[km^3/s^2]
au = 149597870.88 #[km]



###################### KEPLERIAN TO CARTESIAN #########################


def kep2cart(a,e,i,omega,w,theta,E,M):  #all [m] or [degrees]
    mu_e = 398600.441 #[km^3/s^2]
    
    #degrees to radians
    i = i*(np.pi/180)
    omega = omega*(np.pi/180)
    w = w*(np.pi/180)
    theta = theta*(np.pi/180)
    E = E*(np.pi/180)
    M = M*(np.pi/180)
    
    
    #if only M is known, calculate E then theta
    if E == 0:
        E=M
        diff=1
        while abs(diff) > 1e-10:
            diff=(M-E+e*np.sin(E))/(1-e*np.cos(E))
            E=E+diff
        
        theta = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
        if theta < 0:
            theta = (2*np.pi) + theta
    
    r = a*(1-(e*np.cos(E)))

    l1 = np.cos(omega)*np.cos(w) - np.sin(omega)*np.sin(w)*np.cos(i)
    l2 = -np.cos(omega)*np.sin(w) - np.sin(omega)*np.cos(w)*np.cos(i)
    m1 = np.sin(omega)*np.cos(w) + np.cos(omega)*np.sin(w)*np.cos(i)
    m2 = -np.sin(omega)*np.sin(w) + np.cos(omega)*np.cos(w)*np.cos(i)
    n1 = np.sin(w)*np.sin(i)
    n2 = np.cos(w)*np.sin(i)
    mat1 = np.array([[l1, l2], 
                     [m1, m2],
                     [n1, n2]])

    #calculate the second matrix for position

    curlye = r*np.cos(theta)  #curly e to perigee
    longn = r*np.sin(theta)   #long n perpendicular to curly e
    mat2 = np.array([[curlye],
                 [longn]])

    #calculate the position coordinates in km
    coordinates_c = np.dot(mat1,mat2)
    x = coordinates_c[0]
    y = coordinates_c[1]
    z = coordinates_c[2]
    
    coordinates = []
    coordinates.append(x)
    coordinates.append(y)
    coordinates.append(z)
    
    #calculate the angular momentum
    H = np.sqrt(mu_e*(a/1000)*(1-e**2))
 

    #calculate the velocities in m/s
    coordinates.append((mu_e*1000/H)*(-l1*np.sin(theta) + l2*(e+np.cos(theta))))   
    coordinates.append((mu_e*1000/H)*(-m1*np.sin(theta) + m2*(e+np.cos(theta))))  
    coordinates.append((mu_e*1000/H)*(-n1*np.sin(theta) + n2*(e+np.cos(theta))))  

    return coordinates

this = kep2cart(12269687.5912,0.004932091570,109.823277603,134.625563565,106.380426142,0,0,301.149932402)
print this
#
#
##assign array for cartesian coordinates
#cartesian = []
#
##input variables
#
#a = 12269687.5912                         #[m], semi-major axis
#e = 0.004932091570                        #[-], eccentricity
#i = 109.823277603                      #[deg], inclination
#omega = 134.625563565                     #[deg], RAAN
#w = 106.380426142                        #[deg], argument of peri
##theta = 239.5437                     #[deg], true anomaly
##E = 239.5991                         #[deg], eccentric anomaly
#M = 301.149932402                         #[deg], mean anomaly
#
##transform degrees to radians
#i = i*(np.pi/180)
#omega = omega*(np.pi/180)
#w = w*(np.pi/180)
##theta = theta*(np.pi/180)
##E = E*(np.pi/180)
#M = M*(np.pi/180)
#
##calculate E from M through iteration
#E=M
#diff=1
#while abs(diff) > 1e-10:
#            diff=(M-E+e*np.sin(E))/(1-e*np.cos(E))
#            E=E+diff
#
##calculate theta from E
#theta = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
#if theta < 0:
#    theta = (2*np.pi) + theta
#
##calculate r from input variables
#r = a*(1-(e*np.cos(E)))                 #[m]
#
#
########## calculate the cartesian coordinates #######
#
#
##calculate the first matrix for position
#
#l1 = np.cos(omega)*np.cos(w) - np.sin(omega)*np.sin(w)*np.cos(i)
#l2 = -np.cos(omega)*np.sin(w) - np.sin(omega)*np.cos(w)*np.cos(i)
#m1 = np.sin(omega)*np.cos(w) + np.cos(omega)*np.sin(w)*np.cos(i)
#m2 = -np.sin(omega)*np.sin(w) + np.cos(omega)*np.cos(w)*np.cos(i)
#n1 = np.sin(w)*np.sin(i)
#n2 = np.cos(w)*np.sin(i)
#mat1 = np.array([[l1, l2], 
#                 [m1, m2],
#                 [n1, n2]])
#
#
##calculate the second matrix for position
#
#curlye = r*np.cos(theta)  #curly e to perigee
#longn = r*np.sin(theta)   #long n perpendicular to curly e
#mat2 = np.array([[curlye],
#                 [longn]])
#
##calculate the position coordinates in km
#coordinates_c = np.dot(mat1,mat2)
#x = coordinates_c[0]
#y = coordinates_c[1]
#z = coordinates_c[2]
#
#print "x =", x
#print "y =", y
#print "z =",  z #m
#
#cartesian.append(x)
#cartesian.append(y)
#cartesian.append(z)
#
#
##calculate the angular momentum
#H = np.sqrt(mu_e*(a/1000)*(1-e**2))
# 
#
##calculate the velocities in km/s
#v_x = (mu_e/H)*(-l1*np.sin(theta) + l2*(e+np.cos(theta)))   
#v_y = (mu_e/H)*(-m1*np.sin(theta) + m2*(e+np.cos(theta)))   
#v_z = (mu_e/H)*(-n1*np.sin(theta) + n2*(e+np.cos(theta)))  
#
#cartesian.append(v_x*1000)
#cartesian.append(v_y*1000)
#cartesian.append(v_z*1000) 
#
#print "v_x =", v_x * 1000 
#print "v_y =", v_y*1000 
#print"v_z =", v_z*1000 #m/s
#
#
#
#
#
#
