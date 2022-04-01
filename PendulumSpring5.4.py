import sympy as smp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

class System:
    def __init__(self,
                 name = 'Particle', #name of particle
                 time = np.linspace(0, 8, 400), # time over which simulation takes place
                 mass = 1, #mass in terms of kilograms
                 gravity = 9.81, #force of gravity
                 hooke = 10, #hooke's constant for the spring
                 naturalLength = 0.5, #natural length of the spring
                 IV = [0.5, 0, 0, 0], #initial r, v, theta, omega
                 pvtV = 0.2, # pivot velocity
                 dampingThe = 0.1, #pendulum damping factor
                 dampingR = 0.1, #spring damping factor
                 cSymbols = 't m k l g', #constant symbols
                 fSymbols = 'theta r x y F'): #function symbols
        #set inital values to class attributes
        self.name, self.time, self.mass, self.gravity, self.hooke, self.naturalLength = name, time, mass, gravity, hooke, naturalLength
        #define input constants as sympy constants
        self.t, self.m, self.k, self.l, self.g = smp.symbols(cSymbols)
        #define input variables as functions which can be acted on
        self.the, self.r, self.x, self.y, self.F = smp.symbols(fSymbols, cls = smp.Function)
        self.the = self.the(self.t)
        self.r = self.r(self.t)
        self.F = self.F(self.t)
        self.IV = IV
        self.ans = 0
        self.pvtV = pvtV
        self.dampingThe = dampingThe
        self.dampingR = dampingR
    def diffrentiate(self):
    #define functions as what they are a function of, introduce derivatves and second derivatives of functions
        self.the_d = smp.diff(self.the, self.t)
        self.r_d = smp.diff(self.r, self.t)
        self.x_d = smp.diff(self.x, self.t)
        self.y_d = smp.diff(self.y, self.t)
        self.the_dd = smp.diff(self.the, (self.t, 2))
        self.r_dd = smp.diff(self.r, (self.t, 2))
    #produces Euler-Lagrange equations from the Lagrangian which is = 0 because dL/dt = 0 and this is chain ruled to produce the E-L eqn.
    def eulerLagrange(self, L):
        self.EL_the = smp.diff(L, self.the) - smp.diff(smp.diff(L, self.the_d), self.t) + smp.diff(self.F, self.the_d)
        self.EL_the = self.EL_the.simplify()
        self.EL_r = smp.diff(L, self.r) - smp.diff(smp.diff(L, self.r_d), self.t) + smp.diff(self.F, self.r_d)
        self.EL_r = self.EL_r.simplify()
    #assigns the values given to the variables initialisation to the Sympy variables
    def assignValues(self):
        self.t = self.time
        self.m = self.mass
        self.k = self.hooke
        self.l = self.naturalLength
        self.g = self.gravity
    def lambdify_xy(self):
        self.x_f = smp.lambdify([self.r, self.the], self.x)
        self.y_f = smp.lambdify([self.r, self.the], self.y)
    #writes derivates as an array that can be solved using the odeint function from scipy.integrate which was imported at the beginning
    def dZdt(self,Z,t):
        return np.array([self.v(Z[1]),
               self.a(self.k, self.m, self.l, Z[0], self.g, Z[1], Z[2], Z[3]),
               self.w(Z[3]),
               self.dwdt(self.k, self.m, self.l, self.g, Z[0], Z[2], Z[1], Z[3])], dtype=np.float64)

#create particles: all subsequent particles need to be functions of the angles and radiuses of the particles before them
bob1 = System(name = 'p1', IV = [1, 0, np.pi/4, 0], cSymbols = 't m1 k1 l1 g1', fSymbols = 'theta1 r1 x1 y1 F1')
bob2 = System(name = 'p2', hooke = 30, IV = [0.5, 0, np.pi/6, 0], pvtV=0, cSymbols = 't m2 k2 l2 g2', fSymbols = 'theta2 r2 x2 y2 F2')
#add particles to array
Particles = [bob1, bob2]

for i in range(0, len(Particles)):
    dependant = []
    for k in range(0, len(Particles)):
        if k <= i:
            dependant.append(Particles[k].r)
            dependant.append(Particles[k].the)
    dependant = tuple(dependant)
    Particles[i].x = Particles[i].x(dependant)
    Particles[i].y = Particles[i].y(dependant)
#equations of motion for the particle
bob1.x = bob1.r * smp.sin(bob1.the)
bob1.y = - bob1.r * smp.cos(bob1.the)
bob2.x = bob1.x + bob2.r * smp.sin(bob2.the)
bob2.y = - bob1.y + bob2.r * smp.cos(bob2.the)
#diffrentiate
for p in Particles:
    p.diffrentiate()
#create lagrangian

T = 0
V = 0
for p in Particles:
    T += 1/2 * p.m * (p.x_d**2 + p.y_d**2)
    V += p.m * p.g * p.y + 1/2 * p.k * (p.r-p.l)**2
L = T - V
#creates euler-lagrange equations for all the particles in Particles wrt theta are r and adds them to an array
elEquations = []
rtheta_dd = []
constants = []
for p in Particles:
    p.eulerLagrange(L)
    elEquations.append(p.EL_the)
    elEquations.append(p.EL_r)
    rtheta_dd.append(p.the_dd)
    rtheta_dd.append(p.r_dd)
    constants.append(p.m)
    constants.append(p.k)
    constants.append(p.l)
    constants.append(p.g)
    constants.append(p.r)
    constants.append(p.r_d)
    constants.append(p.the)
    constants.append(p.the_d)
rtheta_dd = tuple(rtheta_dd)
constants = tuple(constants)
#rearranges these to form ODEs with second derivative as the subject of the formula.
solution = smp.solve(elEquations, rtheta_dd, simplify=False, rational=False)
#turn these ode to functions that can take numerical and numpy values
solutionsLambdify = []
for key in solution:
    solutionsLambdify.append(smp.lambdify(constants, solution[key]))
#writes derivates as an array that can be solved using the odeint function from scipy.integrate which was imported at the beginning
def dZdt(Z, t):
    r1, v1, the1, w1, r2, v2, the2, w2 = Z
    c = []
    for i in range(0, len(Particles)):
        c.append(Particles[i].m)
        c.append(Particles[i].k)
        c.append(Particles[i].l)
        c.append(Particles[i].g)
        c.append(Particles[i].r)
        c.append(Z[4*i])
        c.append(Z[1+4*i])
        c.append(Z[2+4*i])
        c.append(Z[3+4*1])
    #print(c)
    return [ solutionsLambdify[0](c[1:-1]),
             solutionsLambdify[1](c),
             solutionsLambdify[2](c),
             solutionsLambdify[3](c)]
IV = []
for p in Particles:
    p.assignValues()
    IV += p.IV
    p.lambdify_xy()

ans = odeint(dZdt, y0=IV, t=bob1.t)

def get_x1y1x2y2(r1, r2, the1, the2):
    return (bob1.x_f(r1, the1),
            bob1.y_f(r1, the1),
            bob2.x_f(r1, the1, r2, the2),
            bob2.y_f(r1, the1, r2, the2))
bob1.x, bob1.y, bob2.x, bob2.y = get_x1y1x2y2(ans.T[0], ans.T[4], ans.T[2], ans.T[6])

def animate(i):
            ln1.set_data([0, bob1.x[i]], [0, bob1.y[i]])
            ln2.set_data([bob1.x[i], bob2.x[i]], [bob1.y[i], bob2.y[i]])
fig, ax = plt.subplots(1,1, figsize = (8,8))
ax.grid()
ln1, = plt.plot([], [], 'ro--', lw=3, markersize=8)
ln2, = plt.plot([], [], 'bo--', lw=3, markersize=8)
ax.set_ylim(-3, 1)
ax.set_xlim(-1, 5)
ani = animation.FuncAnimation(fig, animate, frames=400, interval=50)
ani.save('DoubleDampedPendulumSpring.gif', writer='pillow', fps=50)
#plt.plot(p2.x, p2.y)
print('done')
#plt.show()
