# -*- coding: utf-8 -*-
"""
Created on Mon May 20 12:18:38 2019

@author: woute
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math

theta = math.radians(75)
omega = math.radians(100)

L = 0.6
g= 9.81
m =11
t = np.linspace(0,2,1000)

Msmax = 20
f =2

def motor(t):
    return Msmax * np.cos(2*np.pi*f*t)/((1/3)*m*L**2)

def fq(q,t):
    dydt = [q[1],  -3*g/(2*L)* np.sin(q[0]) + motor(t)]
    return dydt




sol = odeint(fq, [omega, theta], t)



plt.plot(t, sol[:,0] , 'b', label = 'theta(t)')
plt.plot(t, sol[:,1] , 'g', label = 'omega(t)')
plt.xlabel('t')
plt.ylabel('rad')
plt.legend(loc = 'best')
plt.grid()
plt.show()
