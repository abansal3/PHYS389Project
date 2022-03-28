import math as m
import matplotlib.pyplot as plt
import numpy as np


class PendulumSpring:
    def __init__(self,
                 name = 'Pendulum',
                 mass = 1.0,
                 theta = m.pi/9,
                 g = -9.81,
                 nl = 1.0,
                 hc = 1,
                 pivot=np.array([0,0], dtype=np.float64),
                 pos=np.array([0,0], dtype=np.float64),
                 angularVel= 0,
                 angularAcc = 0,
                 spring=False):

        pos=[nl*m.sin(theta)+pivot[0],-nl*m.cos(theta)+pivot[1]]
        self.name = name
        self.mass = 1
        self.pos = np.array(pos, dtype=np.float64)
        self.angVel = np.array(angularVel, dtype=np.float64)
        self.angAcc = np.array(angularAcc, dtype=np.float64)
        self.g = g
        self.nl = nl
        self.hc = hc
        self.theta = theta
        self.spring = spring

    def UpdateAcceleration(self):
        self.theta = m.atan(self.pos[0]/self.pos[1])
        self.angAcc = (self.g / self.nl) * np.sin(self.theta)

    def ER_update(self,deltaT):
        self.theta += 0.5 * self.angVel * deltaT
        self.angVel += 0.5 * self.angAcc * deltaT
        self.pos = [self.nl * np.sin(self.theta), self.nl * np.cos(self.theta)]
            
            
        


Pendulum1 = PendulumSpring(name='Pendulum1')
deltaT = 0.01
T = []
posx = []
posy = []
length = []
for i in range(1,300):
    Pendulum1.UpdateAcceleration()
    Pendulum1.ER_update(deltaT)
    T.append(i)
    posx.append(Pendulum1.pos[0])
    posy.append(Pendulum1.pos[1])
    length.append(m.sqrt(sum(x**2 for x in Pendulum1.pos)))
    

fig = plt.figure()
ax = plt.axes()
ax.set_ylabel('y')
ax.set_xlabel('x')
plt.scatter(posx, posy, marker = '.')
plt.show()
