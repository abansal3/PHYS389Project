import math as m
import matplotlib.pyplot as plt


class PendulumSpring:
    def __init__(self, name = 'Pendulum', mass = 1.0, theta = m.pi/9, g = -9.81, nl = 1.0, hc = 1, pos=[0,0,0], vel=[0,0,0], acln = [0,0,0], spring=False):

        pos=[nl*m.sin(theta),nl*m.cos(theta),0]
        self.name = name
        self.mass = 1
        self.pos = pos
        self.vel = vel
        self.acln = acln
        self.g = g
        self.nl = nl
        self.hc = hc
        self.thetta = theta
        self.spring = spring

    def UpdateAcceleration():
        self.theta = abs(m.atan(self.pos[0]/self.pos[1]))
        self.acln = [self.g*m.tan(self.theta), self.g*1-m.cos(self.theta), 0]
        if self.spring == True:
            ext = sqrt(sum(x**2 for x in self.pos)) - self.nl
            acln = self.hc * ext / self.mass
            self.acln += [acln*m.sin(self.theta), acln*m.cos(self.theta), 0]

    def ER_update(self,deltaT):
        self.position += 0.5 * self.velocity * deltaT
        self.velocity += 0.5 * self.acceleration * deltaT
            
            
        


Pendulum1 = PendulumSpring(name='Pendulum1')
deltaT = 1
Data = []
for i in range(1,100):
    Pendulum1.UpdateAcceleration()
    Pendulum1.ER_update(deltaT)
    T.append(i)
    posx.append(Pendulum1.pos[0])
    posy.append(Pendulum1.pos[1])
    length.append(sqrt(sum(x**2 for x in self.pos)))
    

fig = plt.figure()
ax = plt.axes()
ax.set_ylabel('y')
ax.set_xlabel('x')
plt.scatter(posx, posy)
plt.shwow()
