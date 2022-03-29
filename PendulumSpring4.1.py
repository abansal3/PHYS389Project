import sympy as smp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter


global self
self.t, self.m1, self.m2, self.k, self.l, self.g, self.x1_d = smp.symbols('t m1 m2 k l g x1_d')
#define input variables as functions which can be acted on
self.the, self.x1, self.x2, self.y2, self.r = smp.symbols(r'theta x_{1} x_{2} y_{2} r', cls = smp.Function)
#define functions as what they are a function of, introduce derivatves and second derivatives of functions
self.the = self.the(self.t)
self.x1 = self.x1(self.x1_d, self.t)
self.r = self.r(self.t)
self.x2 = self.x2(self.r, self.the)
self.y2 = self.y2(self.r, self.the, self.x1)
self.the_d = smp.diff(self.the, self.t)
self.r_d = smp.diff(self.r, self.t)
self.the_dd = smp.diff(self.the, (self.t, 2))
self.r_dd = smp.diff(self.r, (self.t, 2))
#define equations describing these variables and their derivatives
self.x1 = self.x1_d * self.t
self.x2 = self.x1 + self.r * smp.sin(self.the)
self.y2 = - self.r * smp.cos(self.the)
self.x2_d = smp.diff(self.x2, self.t)
self.y2_d = smp.diff(self.y2, self.t)
#turn these symbol equations into ones that can take values
self.x1_f = smp.lambdify([self.x1_d, self.t], self.x1)
self.x2_f = smp.lambdify([self.x1_d, self.t, self.r, self.the], self.x2)
self.y2_f = smp.lambdify([self.r, self.the], self.y2)

#form Langrangian
self.T = 1/2 * (self.m1 * self.x1_d ** 2 + self.m2 * (self.x2_d ** 2 + self.y2_d ** 2))
self.V = - self.m2 * self.g * self.r * smp.cos(self.the) + 1/2 * self.k * (self.r - self.l)**2
self.L = self.T - self.V
self.L = self.L.simplify()
self.LE_the = smp.diff(self.L, self.the) - smp.diff(smp.diff(self.L, self.the_d), self.t)
self.LE_r = smp.diff(self.L, self.r) - smp.diff(smp.diff(self.L, self.r_d), self.t)
self.LE_the = self.LE_the.simplify()
self.LE_r = self.LE_r.simplify()
#represent langrange-euler equations as a ODE
self.the_deriv2 = smp.solve(self.LE_the, self.the_dd)
self.the_deriv1 = self.the_d
self.r_deriv2 = smp.solve(self.LE_r, self.r_dd)
self.r_deriv1 = self.r_d
#turn these ode to functions that can take values
self.the_deriv2_f = smp.lambdify([self.g,self.r, self.the, self.r_d, self.the_d],self.the_deriv2)
self.the_deriv1_f = smp.lambdify(self.the_d, self.the_d)
self.r_deriv2_f = smp.lambdify([self.k, self.m2, self.l, self.r, self.g, self.the, self.the_d],self.r_deriv2)
self.r_deriv1_f = smp.lambdify(self.r_d,self.r_d)
#give values to constants
t = np.linspace(0,20, 1000)
self.m1 = 1
self.m2 = 1
self.k = 3
self.l = 0.5
self.g = 9.81
self_x1_d = 0
 
#create Z as a function of r, dr/dt, theta, dtheta/dt
def dZdt(Z,t):
    return np.array([self.r_deriv1_f(Z[1]),
           self.r_deriv2_f(self.k, self.m2, self.l, Z[0], self.g, Z[2], Z[3])[0],
           self.the_deriv1_f(Z[3]),
           self.the_deriv2_f(self.g, Z[0], Z[2], Z[1], Z[3])[0]], dtype=np.float64)
#solve ODE
self.ans = odeint(dZdt, y0=[0.5, 0, np.pi/9, 0], t=t)

#convert this to x and y coordinates for particles
def getx1x2y2(t, the, r):
    return (self.x1_f(self.x1_d, self.t),
            self.x2_f(self.x1_d, self.t, r, the),
            self.y2_f(r, the))
self.x1, self.x2, self.y2 = getx1x2y2(self.t, self.ans.T[2], self.ans.T[0])
def animate(i):
    ln1.set_data([0, x2[i]], [0, y2[i]])
fig, ax = plt.subplots(1,1, figsize = (8,8))
ax.grid()
ln1, = plt.plot([], [], 'ro--', lw=3, markersize=8)
ax.set_ylim(-10, -10)
ax.set_xlim(-10, -10)
self.ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
