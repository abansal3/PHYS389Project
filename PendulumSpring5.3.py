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
                 pvtV = 0.2,
                 dampingThe = 0.1,
                 dampingR = 0.1): 
        #set inital values to class attributes
        self.name, self.time, self.mass, self.gravity, self.hooke, self.naturalLength = name, time, mass, gravity, hooke, naturalLength
        #define input constants as sympy constants
        self.t, self.m, self.k, self.l, self.g = smp.symbols('t m k l g')
        #define input variables as functions which can be acted on
        self.the, self.r, self.x, self.y, self.F = smp.symbols(r'theta r x y F', cls = smp.Function)
        self.the = self.the(self.t)
        self.r = self.r(self.t)
        self.x = self.x(self.t)
        self.y = self.y(self.t)
        self.F = self.F(self.t)
        self.IV = IV
        self.ans = 0
        self.pvtV = pvtV
        self.dampingThe = dampingThe
        self.dampingR = dampingR
    #form lagrange equation using the summation of all the kinetic and potential energies
    def lagrangian(self):
        #define functions as what they are a function of, introduce derivatves and second derivatives of functions
        self.the_d = smp.diff(self.the, self.t)
        self.r_d = smp.diff(self.r, self.t)
        self.x_d = smp.diff(self.x, self.t)
        self.y_d = smp.diff(self.y, self.t)
        self.the_dd = smp.diff(self.the, (self.t, 2))
        self.r_dd = smp.diff(self.r, (self.t, 2))
        T = 1/2 * self.m * (self.x_d ** 2 + self.y_d ** 2)
        V =  self.m * self.g * self.y + 1/2 * self.k * (self.r - self.l)**2
        L = T - V
        return L
    #produces Euler-Lagrange equations from the Lagrangian which is = 0 because dL/dt = 0 and this is chain ruled to produce the E-L eqn.
    def eulerLagrange(self, L):
        LE_the = smp.diff(L, self.the) - smp.diff(smp.diff(L, self.the_d), self.t) + smp.diff(self.F, self.the_d)
        LE_r = smp.diff(L, self.r) - smp.diff(smp.diff(L, self.r_d), self.t) + smp.diff(self.F, self.r_d)
        #rearranges these to form ODEs with second derivative as the subject of the formula.
        derivs2 = smp.solve([LE_the, LE_r], (self.the_dd, self.r_dd))
        print(self.name, LE_the, LE_r, derivs2, sep='\n')
        r_deriv1 = self.r_d
        the_deriv1 = self.the_d
        #turn these ode to functions that can take numerical and numpy values
        self.dwdt = smp.lambdify([self.k, self.m, self.l, self.g,self.r, self.the, self.r_d, self.the_d],derivs2[self.the_dd])
        self.w = smp.lambdify(self.the_d, self.the_d)
        self.a = smp.lambdify([self.k, self.m, self.l, self.r, self.g, self.r_d, self.the, self.the_d],derivs2[self.r_dd])
        self.v = smp.lambdify(self.r_d,self.r_d)
    #assigns the values given to the variables initialisation to the Sympy variables
    def assignValues(self):
        self.t = self.time
        self.m = self.mass
        self.k = self.hooke
        self.l = self.naturalLength
        self.g = self.gravity
    def lambdify_xy(self):
        self.x_f = smp.lambdify([self.t, self.r, self.the], self.x)
        self.y_f = smp.lambdify([self.r, self.the], self.y)
    #writes derivates as an array that can be solved using the odeint function from scipy.integrate which was imported at the beginning
    def dZdt(self,Z,t):
        return np.array([self.v(Z[1]),
               self.a(self.k, self.m, self.l, Z[0], self.g, Z[1], Z[2], Z[3]),
               self.w(Z[3]),
               self.dwdt(self.k, self.m, self.l, self.g, Z[0], Z[2], Z[1], Z[3])], dtype=np.float64)
#creates particles
Spring = System()
p = System(name = 'p1', IV = [1, 0, np.pi/4, 0])
p2 = System(name = 'p2', hooke = 30, IV = [0.5, 0, np.pi/6, 0], pvtV=0, dampingThe = 0, dampingR  = 0)
#sets any other inital condition values of the system (length, velocity, angle, omega)
pvtV, dampingThe, dampingR = 0.2, 0.1, 0.1
#equations of motion for the particle
p.x = pvtV * p.t + p.r * smp.sin(p.the)
p.y = - p.r * smp.cos(p.the)
p2.x = p.x + p2.r * smp.sin(p2.the)
p2.y = p.y + p2.r * smp.cos(p2.the)
#create Langragian and ODEs for the system
P_Lagrange = p.lagrangian()
P2_Lagrange = p2.lagrangian()
p.F = - dampingThe * p.the_d ** 2 - dampingR * p.r_d **2
p2.F = - dampingThe * p2.the_d ** 2 - dampingR * p2.r_d **2
p.eulerLagrange(P_Lagrange)
p2.eulerLagrange(P2_Lagrange)
#assign initial values to the particle such that they can be entered into the ODEs
p.lambdify_xy()
p.assignValues()
p2.lambdify_xy()
p2.assignValues()
#solve ODE
for i in [p, p2]:
    i.ans = odeint(i.dZdt, y0=i.IV, t=i.t)
def getXY(p, t, r, the):
    return (p.x_f(t, r, the),
            p.y_f(r, the))
p.x, p.y = getXY(p, p.t, p.ans.T[0], p.ans.T[2])
p2.x, p2.y = getXY(p2, p2.t, p2.ans.T[0], p2.ans.T[2])
print('calculated')
def animate(i):
            ln1.set_data([pvtV*p.t[i], p.x[i]], [0, p.y[i]])
            ln2.set_data([p.x[i], p2.x[i]], [p.y[i], p2.y[i]])
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
