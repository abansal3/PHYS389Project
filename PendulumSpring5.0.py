import sympy as smp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

class System:
    def __init__(self,
                 name = 'Particle', #name of particle
                 time = np.linspace(0, 20, 1000), # time over which simulation takes place
                 mass = 1, #mass in terms of kilograms
                 gravity = 9.81, #force of gravity
                 hooke = 10, #hooke's constant for the spring
                 naturalLength = 0.5): #natural length of the spring
        #set inital values to class attributes
        self.name, self.time, self.mass, self.gravity, self.hooke, self.naturalLength = name, time, mass, gravity, hooke, naturalLength
        #define input constants as sympy constants
        self.t, self.m, self.k, self.l, self.g = smp.symbols('t m k l g')
        #define input variables as functions which can be acted on
        self.the, self.r, self.x, self.y = smp.symbols(r'theta r x y', cls = smp.Function)
        self.the = self.the(self.t)
        self.r = self.r(self.t)
        self.x = self.x(self.t)
        self.y = self.y(self.t)
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
        LE_the = smp.diff(L, self.the) - smp.diff(smp.diff(L, self.the_d), self.t)
        LE_r = smp.diff(L, self.r) - smp.diff(smp.diff(L, self.r_d), self.t)
        #rearranges these to form ODEs with second derivative as the subject of the formula.
        the_deriv2 = smp.solve(LE_the, self.the_dd)
        the_deriv1 = self.the_d
        r_deriv2 = smp.solve(LE_r, self.r_dd)
        if r_deriv2 == []: r_deriv2 = [0]
        r_deriv1 = self.r_d
        #turn these ode to functions that can take numerical and numpy values
        the_deriv2_f = smp.lambdify([self.l, self.g,self.r, self.the, self.r_d, self.the_d],the_deriv2)
        the_deriv1_f = smp.lambdify(self.the_d, self.the_d)
        r_deriv2_f = smp.lambdify([self.k, self.m, self.l, self.r, self.g, self.the, self.the_d],r_deriv2)
        r_deriv1_f = smp.lambdify(self.r_d,self.r_d)
        return the_deriv1_f, the_deriv2_f, r_deriv1_f, r_deriv2_f
    #assigns the values given to the variables initialisation to the Sympy variables
    def assignValues(self):
        self.t = self.time
        self.m = self.mass
        self.k = self.hooke
        self.l = self.naturalLength
        self.g = self.gravity
    def lambdify_xy(self):
        self.x_f = smp.lambdify(self.the, self.x)
        self.y_f = smp.lambdify(self.the, self.y)
        
#creates particles
Spring = System()
p = System(hooke = 100000000000000000000000000000000000000000)
#sets intial condition values of the system (length, velocity, angle, omega)
p_IV = [0.5, 0, np.pi/3, 0]
#equations of motion for the particle
p.x = p.naturalLength * smp.cos(p.the)
p.y = p.naturalLength * smp.sin(p.the)
#create Langragian and ODEs for the system
P_Lagrange = p.lagrangian()
p_omega, p_domegadt, p_v, p_a = p.eulerLagrange(P_Lagrange)
#assign initial values to the particle such that they can be entered into the ODEs
p.assignValues()
p.lambdify_xy()
#writes derivates as an array that can be solved using the odeint function from scipy.integrate which was imported at the beginning
def dZdt(Z,t):
    return np.array([p_v(Z[1]),
           p_a(p.k, p.m, p.l, Z[0], p.g, Z[2], Z[3])[0],
           p_omega(Z[3]),
           p_domegadt(p.l, p.g, Z[0], Z[2], Z[1], Z[3])[0]], dtype=np.float64)
#solve ODE
ans = odeint(dZdt, y0=p_IV, t=p.t)
def getXY(the):
    return (p.x_f(the),
            p.y_f(the))
p.x, p.y = getXY(ans.T[2])
fig, ax = plt.subplots(1,1, figsize = (8,8))
ax.grid()
print(p.x, p.y, sep='\n')
plt.plot(p.x, p.y)
plt.show()
