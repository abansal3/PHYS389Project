import sympy as smp
import numpy as np
from scipy.integrate import odeint
from Particle import System

#creates particles
bob1 = System(name = 'Spring', IV = [1, 0, (np.pi/180)*0, 0])
bob2 = System(name = 'Pendulum', IV = [0.3*9.81, 0, (np.pi/180)*30, 0], hooke = 1000000)
bob3 = System(name = 'Stationary Elastic Pendulum', IV = [1, 0, (np.pi/180)*20, 0])
bob4 = System(name = 'Moving Elastic Pendulum', IV = [1, 0, (np.pi/180)*20, 0], pvtV = 0.2)
bob5 = System(name = 'Damped Stationary Elastic Pendulum', IV = [1, 0, (np.pi/180)*20, 0], dampingThe = 0.5, dampingR = 0.2)
bob6 = System(name = 'Damped Moving Elastic Pendulum', IV = [1, 0, (np.pi/180)*20, 0], pvtV = 0.2, dampingThe = 0.5, dampingR = 0.2)

bob7 = System(name = '1', IV = [3.04, 0, (np.pi/180)*20.2, 0], naturalLength = 0.4)
bob8 = System(name = '2', IV = [5.8, 0, (np.pi/180)*19.5, 0], naturalLength = 2)
bob9 = System(name = '3', IV = [8.77, 0, (np.pi/180)*20.5, 0], naturalLength = 4)
bob10 = System(name = '4', IV = [11.2, 0, (np.pi/180)*45.5, 0], naturalLength = 4)
bob11 = System(name = '5', IV = [19.1, 0, (np.pi/180)*46.8, 0], naturalLength = 9)
bob12 = System(name = '6', IV = [26.3, 0, (np.pi/180)*46.9, 0], naturalLength = 14)

#add particles to array
Particles = [bob7, bob8, bob9, bob10, bob11, bob12]
for p in Particles:
    #equations of motion for the particle
    p.x = p.pvtV * p.t + p.r * smp.sin(p.the)
    p.y = - p.r * smp.cos(p.the)
    #create Langragian and ODEs for the system
    p.lagrangian()
    p.F = - p.dampingThe * p.the_d ** 2 - p.dampingR * p.r_d **2
    p.eulerLagrange()
    #assign initial values to the particle such that they can be entered into the ODEs
    p.lambdifyFuncs()
    p.assignValues()
#solve ODE
for p in Particles:
    p.ans = odeint(p.dZdt, y0=p.IV, t=p.t)
    #gets data for position, velocities and energies for all the particles
    def getXY(t, r, r_d, the, the_d):
        return (p.x_f(t, r, the),
                p.y_f(r, the),
                p.x_df(r, r_d, the, the_d),
                p.y_df(r, r_d, the, the_d),
                p.T_f(p.m, r, r_d, the, the_d),
                p.V_f(p.m, p.g, p.k, p.l, r, the),
                p.E_f(p.m, p.g, p.k, p.l, r, r_d, the, the_d))
    p.x, p.y, p.x_d, p.y_d, p.T, p.V, p.E = getXY(p.t, p.ans.T[0], p.ans.T[1], p.ans.T[2], p.ans.T[3])
    #calculate R and mu
    E_start = p.E_f(p.m, p.g, p.k, p.l, p.IV[0], p.IV[1], p.IV[2], p.IV[3])
    E_min = p.E_f(p.m, p.g, p.k, p.l, p.l + p.m * p.g / p.k, 0, 0, 0)
    p.R = E_start / E_min
    p.mu = 1 + (p.k * p.l) / (p.m * p.g)
    print('bob number {} has R={} and mu={}'.format(p.name, p.R, p.mu))

#stores the information found in an external file
    p.Data = np.array([p.t, p.x, p.y, p.x_d, p.y_d, p.T, p.V, p.E, p.ans.T[0], p.ans.T[1],p.ans.T[2], p.ans.T[3], p.R, p.mu], dtype=object)

Data = ['Data1', 'Data2', 'Data3', 'Data4', 'Data5', 'Data6']
for i in range(0, len(Particles)):
    np.save(Data[i], Particles[i].Data, allow_pickle=True)

