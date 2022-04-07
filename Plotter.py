import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter
import numpy as np

#load data plots into plotter
DataNames = ['Data1.npy', 'Data2.npy', 'Data3.npy', 'Data4.npy', 'Data5.npy', 'Data6.npy']
Data1, Data2, Data3, Data4, Data5, Data6 = [], [], [], [], [], []
Data = [Data1, Data2, Data3, Data4, Data5, Data6]
for i in range(0, len(Data)):
    Data[i].append(list(np.load(DataNames[i], allow_pickle=True)))

def createPlots(Data1):
    #creat plot of motion of the particle (x against y)
    fig1, ax1 = plt.subplots(1, 1)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.plot(Data1[1], Data1[2])

    #create plot of energies of the system against time
    fig2, ax2 = plt.subplots(1, 1)
    ax2.set_xlabel('t')
    ax2.set_ylabel('Energy')
    ax2.plot(Data1[0], Data1[5])
    ax2.plot(Data1[0], Data1[6])
    ax2.plot(Data1[0], Data1[7])
    
    #create Poincare sections for x, y, r, and theta respectively.
    fig3, ((ax3, ax4), (ax5, ax6)) = plt.subplots(2, 2)
    #fig3.suptitle('Poincare sections for R = {}, and mu = {}'.format(Data1[12], Data1[13]))
    ax3.set_xlabel('x')
    ax3.set_ylabel('dx / dt')
    ax3.scatter(Data1[1], Data1[3], marker='.')

    ax4.set_xlabel('y')
    ax4.set_ylabel('dy / dt')
    ax4.scatter(Data1[2], Data1[4], marker='.')

    ax5.set_xlabel('r')
    ax5.set_ylabel('dr / dt')
    ax5.scatter(Data1[8], Data1[9], marker='.')

    ax6.set_xlabel('theta')
    ax6.set_ylabel('d theta / dt')
    ax6.scatter(Data1[10], Data1[11], marker='.')

    plt.show()

#create gif of the motion of the particle
def createGif(Data):
    def animate(i):
            ln1.set_data([0, Data[1][i]], [0, Data[2][i]])
    fig, (ax1, ax2) = plt.subplots(2, 1)
    plt.grid()
    ln1, = plt.plot([], [], 'ro--', lw=3, markersize=8)
    ax.set_ylim(-7, 1)
    ax.set_xlim(-3, 3)
    ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50) #frames is amount of time slices
    ani.save('muis4.gif', writer='pillow', fps=50) # interval = fames = slices per second

def poincares():
    fig4, ((ax7, ax8, ax9), (ax10, ax11, ax12)) = plt.subplots(2,3)
    axes = [ax7, ax8, ax9, ax10, ax11, ax12]

    for i in range(0, len(axes)):
        axes[i].set_title('R={}, mu = {}'.format(Data[i][0][12], Data[i][0][13]))
        axes[i].set_xlabel('r')
        axes[i].set_ylabel('d r / dt')
        axes[i].scatter(Data[i][0][10], Data[i][0][11], marker='.')

##    ax7.set_title('R={}, mu = {}'.format(Data1[0][12], Data1[0][13]))
##    ax7.set_xlabel('theta')
##    ax7.set_ylabel('d theta / dt')
##    ax7.scatter(Data1[0][10], Data1[0][11], marker='.')
##
##    ax8.set_title('R={}, mu = {}'.format(Data2[0][12], Data2[0][13]))
##    ax8.set_xlabel('theta')
##    ax8.set_ylabel('d theta / dt')
##    ax8.scatter(Data2[0][10], Data2[0][11], marker='.')
##
##    ax9.set_title('R={}, mu = {}'.format(Data3[0][12], Data3[0][13]))
##    ax9.set_xlabel('theta')
##    ax9.set_ylabel('d theta / dt')
##    ax9.scatter(Data3[0][10], Data3[0][11], marker='.')
##
##    ax10.set_title('R={}, mu = {}'.format(Data4[0][12], Data4[0][13]))
##    ax10.set_xlabel('theta')
##    ax10.set_ylabel('d theta / dt')
##    ax10.scatter(Data4[0][10], Data4[0][11], marker='.')
##
##    ax11.set_title('R={}, mu = {}'.format(Data5[0][12], Data5[0][13]))
##    ax11.set_xlabel('theta')
##    ax11.set_ylabel('d theta / dt')
##    ax11.scatter(Data5[0][10], Data5[0][11], marker='.')
##
##    ax12.set_title('R={}, mu = {}'.format(Data6[0][12], Data6[0][13]))
##    ax12.set_xlabel('theta')
##    ax12.set_ylabel('d theta / dt')
##    ax12.scatter(Data6[0][10], Data6[0][11], marker='.')

    plt.show()

poincares()
