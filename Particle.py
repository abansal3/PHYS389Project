import sympy as smp
import numpy as np

class System:
    def __init__(self,
                 name = 'Particle', #name of particle
                 time = np.linspace(0, 50, 5000), # time over which simulation takes place
                 mass = 1, #mass in terms of kilograms
                 gravity = 9.81, #force of gravity
                 hooke = 9.81, #hooke's constant for the spring
                 naturalLength = 1, #natural length of the spring
                 IV = [0.5, 0, 0, 0], #initial r, v, theta, omega
                 pvtV = 0., # pivot velocity
                 dampingThe = 0, #pendulum damping factor
                 dampingR = 0): #spring damping factor
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
        self.R = 0
        self.mu = 0
        self.Data = []  
    #form lagrange equation using the summation of all the kinetic and potential energies
    def lagrangian(self):
        #define functions as what they are a function of, introduce derivatves and second derivatives of functions
        self.the_d = smp.diff(self.the, self.t)
        self.r_d = smp.diff(self.r, self.t)
        self.x_d = smp.diff(self.x, self.t)
        self.y_d = smp.diff(self.y, self.t)
        self.the_dd = smp.diff(self.the, (self.t, 2))
        self.r_dd = smp.diff(self.r, (self.t, 2))
        self.T = 1/2 * self.m * (self.x_d ** 2 + self.y_d ** 2)
        self.V =  self.m * self.g * self.y + 1/2 * self.k * (self.r - self.l)**2
        self.L = self.T - self.V
        self.E = self.T + self.V
    #produces Euler-Lagrange equations from the Lagrangian which is = 0 because dL/dt = 0 and this is chain ruled to produce the E-L eqn.
    def eulerLagrange(self):
        LE_the = smp.diff(self.L, self.the) - smp.diff(smp.diff(self.L, self.the_d), self.t) + smp.diff(self.F, self.the_d)
        LE_r = smp.diff(self.L, self.r) - smp.diff(smp.diff(self.L, self.r_d), self.t) + smp.diff(self.F, self.r_d)
        #rearranges these to form ODEs with second derivative as the subject of the formula.
        the_deriv2 = smp.solve(LE_the, self.the_dd)
        the_deriv1 = self.the_d
        r_deriv2 = smp.solve(LE_r, self.r_dd)
        if r_deriv2 == []: r_deriv2 = [0]
        r_deriv1 = self.r_d
        #turn these ode to functions that can take numerical and numpy values
        self.dwdt = smp.lambdify([self.m, self.l, self.g,self.r, self.the, self.r_d, self.the_d],the_deriv2)
        self.w = smp.lambdify(self.the_d, self.the_d)
        self.a = smp.lambdify([self.k, self.m, self.l, self.r, self.g, self.r_d, self.the, self.the_d],r_deriv2)
        self.v = smp.lambdify(self.r_d,self.r_d)
    #assigns the values given to the variables initialisation to the Sympy variables
    def assignValues(self):
        self.t = self.time
        self.m = self.mass
        self.k = self.hooke
        self.l = self.naturalLength
        self.g = self.gravity
    def lambdifyFuncs(self):
        self.x_f = smp.lambdify([self.t, self.r, self.the], self.x)
        self.y_f = smp.lambdify([self.r, self.the], self.y)
        self.x_df = smp.lambdify([self.r, self.r_d, self.the, self.the_d], self.x_d)
        self.y_df = smp.lambdify([self.r, self.r_d, self.the, self.the_d], self.y_d)
        self.T_f = smp.lambdify([self.m, self.r, self.r_d, self.the, self.the_d], self.T)
        self.V_f = smp.lambdify([self.m, self.g, self.k, self.l, self.r, self.the], self.V)
        self.E_f = smp.lambdify([self.m, self.g, self.k, self.l, self.r, self.r_d, self.the, self.the_d], self.E)
    #writes derivates as an array that can be solved using the odeint function from scipy.integrate which was imported at the beginning
    def dZdt(self, Z,t):
        return np.array([self.v(Z[1]),
               self.a(self.k, self.m, self.l, Z[0], self.g, Z[1], Z[2], Z[3])[0],
               self.w(Z[3]),
               self.dwdt(self.m, self.l, self.g, Z[0], Z[2], Z[1], Z[3])[0]], dtype=np.float64)
