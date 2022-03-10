import math as m
import matplotlib.pyplot as plt
import numpy as np


class PendulumSpring:
    def __init__(self,
                 name = 'Pendulum',
                 mass = 1.0,
                 theta = m.pi/9,
                 g = 9.81,
                 nl = 1.0,
                 hc = 1,
                 pivot=np.array([0,0,0], dtype=np.float64),
                 pos=np.array([0,0,0], dtype=np.float64),
                 vel=np.array([0,0,0], dtype=np.float64),
                 acln = np.array([0,0,0], dtype=np.float64),
                 spring=False):

        pos=[nl*m.sin(theta)+pivot[0],-nl*m.cos(theta)+pivot[1],0]
        self.name = name
        self.mass = 1
        self.pos = np.array(pos, dtype=np.float64)
        self.vel = np.array(vel, dtype=np.float64)
        self.acln = np.array(acln, dtype=np.float64)
        self.g = g
        self.nl = nl
        self.hc = hc
        self.theta = theta
        self.spring = spring

    def UpdateAcceleration(self):
        self.theta = m.atan(self.pos[0]/self.pos[1])
        velocity = m.sqrt(sum(x**2 for x in self.vel))
        length = m.sqrt(sum(x**2 for x in self.pos))
        self.acln = np.array([((velocity**2)/self.nl + self.g*np.cos(self.theta))*np.sin(self.theta), velocity**2*np.cos(self.theta)/self.nl + self.g*(np.cos(self.theta)**2-1) , 0], dtype=np.float64)
        if self.spring == True:
            ext = length - self.nl
            acc = self.hc * ext / self.mass
            self.acln += [acc*m.sin(self.theta), acc*m.cos(self.theta), 0]

    def ER_update(self,deltaT):
        self.pos += 0.5 * self.vel * deltaT
        self.vel += 0.5 * self.acln * deltaT
            
            
        


Pendulum1 = PendulumSpring(name='Pendulum1')
deltaT = 0.01
T = []
posx = []
posy = []
length = []
for i in range(1,1000):
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
