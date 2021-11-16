# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:23:49 2019

@author: woute
"""

from sympy.physics.mechanics import dynamicsymbols, ReferenceFrame, Point, RigidBody, inertia 
from sympy import symbols 
from sympy.physics.mechanics.kane import KanesMethod 
from numpy import linspace, deg2rad, rad2deg, sin
import matplotlib.pyplot as plt 
from pydy.system import System


# Symbolen
q1, q2 = dynamicsymbols('q1 q2')
u1, u2 = dynamicsymbols('u1 u2')
F = dynamicsymbols('F')
l1, l2, m1, m2, g = symbols('l1 l2 m1 m2 g')
t = symbols('time')

# Assenstelsels
N = ReferenceFrame('N')
A = N.orientnew('A', 'Axis', [q1, N.z])
B = A.orientnew('B', 'Axis', [-q2, A.z])

# Joints locatie en snelheid
J0 = Point('J0')
J1 = J0.locatenew('J1', -l1 * A.x)#J1 op - lengte A.x
J0.set_vel(N, 0)# J0 op 0 in N assenstelsel
J1.v2pt_theory(J0,N,A)

# Massamiddelpunten locatie en snelheid
P = J0.locatenew('P', 0.5 * l1 * A.x)
R = J1.locatenew('R', 0.5 * l2 * B.y)
P.v2pt_theory(J0, N, A)
R.v2pt_theory(J1, N, B)

# Starre lichamen
IP = inertia(A,1/12 * m1 * l1**2 ,0 ,1/12 * m2 * l2**2) 
IP_tuple = (IP,P) 
IR = inertia(B,1.068, 0.496, 1.2343)
IR_tuple = (IR,R) 
BodyP = RigidBody('BodyP', P, A, m1, IP_tuple) 
BodyR = RigidBody('BodyR', R, B, m1, IR_tuple) 
BL = [BodyP, BodyR]

# Model krachten
FL = [(P, - m1 * g * N.y), (R, -m2 * g * N.y), (R, F* N.y)]

# Opstellen bewegingsvergelijking
kd = [A.ang_vel_in(N).dot(N.z) + u1, B.ang_vel_in(N).dot(N.z) + u2] 
KM = KanesMethod(N, q_ind = [q1, q2], u_ind = [u1, u2], kd_eqs = kd) 
KM.kanes_equations(FL, BL)

# Integreren bewegingsvergelijking
sys = System(KM, constants = {m1: 4.2, m2: 42, g: 9.81, l1:0.55, l2: 0.3}, 
             specifieds={F: lambda x, t: sin(t)},
             initial_conditions = {q1: deg2rad(43), q2: deg2rad(43), u1: 0, u2:0}, 
             times = linspace(0,0.3,1000)) 
result = sys.integrate()

# Plotten van grafieken
plt.grid(True)
plt.plot(sys.times, rad2deg(result[:,0])) 
#plt.plot(sys.times, rad2deg(result[:,1])) 
plt.plot(sys.times, rad2deg(result[:,2])) 
#plt.plot(sys.times, rad2deg(result[:,3])) 
#plt.plot(sys.times, result)
plt.legend((str(q1), str(q2))) 
plt.title('Hoeken van een dubbele slinger') 
plt.xlabel('Tijd (sec)')
plt.ylabel('Hoek (graden)')
plt.show()

