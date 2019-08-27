from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import integrate

def solve_circuit(A=1,phi=0):
    fig = plt.figure()
    ax = fig.add_axes([0,0,1,1])
    t = np.linspace(0, 1, 250)
    y = A*np.sin(2*np.pi*t+phi)
    lines = plt.semilogx(t, y)
    plt.setp(lines, linewidth=2, color='r')
    plt.show()
    return t, y

def solve_linSystem(R1 = 100e-3, uL1 = 100, R2 = 100e-3, uL2 = 100, k = 0.8, Rs = 1, U1 = 100, f0 = 20e3, f1 = 100e3, f2 = 10e3):
    w0 = 2*np.pi*f0;
    L1 = 1e-6*uL1
    L2 = 1e-6*uL2
    C1 = 1/(np.square(w0)*L1)
    C2 = 1/(np.square(w0)*L2)
    M = k*np.sqrt(L1*L2)
    
    w = np.linspace(2*np.pi*f1, 2*np.pi*f2, 500)
    Z11 = R1 + 1j*w*L1 + 1/(1j*w*C1)
    Z22 = R2 + Rs + 1j*w*L2 + 1/(1j*w*C2)
    Z12 = -1j*w*M
    Z21 = Z12
    A = 1/(Z11*Z22-Z12*Z21); #A = 1/(Z11*Z22-Z12*Z21);
    
    I1 = A*Z22*U1      #I1 = A*(Z22*V1  - Z12*V2);
    I2 = -A*Z12*U1     #I2 = A*(-Z12*V1 + Z11*V2);
    
    Pu = Rs*I2*np.conj(I2)
    Pu = Pu.real #puissance utile
    theta = np.angle(I1) 
    Pf = U1*np.absolute(I1)*np.cos(theta) #puissance fournie
    eta = Pu.real/Pf #rendement

    plt.close('all')
    fig = plt.figure()
    plt.rcParams["figure.figsize"] = (20,20)
   
    lines = plt.semilogx(w/2*np.pi, eta)
    plt.title('Rendement')
    plt.setp(lines, linewidth=2, color='r')
    plt.xlabel('Frequence [Hz]')
    
    fig = plt.figure()
    line_down, = plt.semilogx(w/2*np.pi, Pu)
    line_up, = plt.semilogx(w/2*np.pi, Pf, 'r')
    plt.legend([line_up, line_down], ['Puissance fournie', 'Puissance utile'])
    plt.xlabel('Frequence [Hz]')
    plt.ylabel('Puissance [W]')
    plt.title('Puissance fournie & utile')

    
    fig = plt.figure()
    lines = plt.semilogx(w/2*np.pi, theta)
    plt.xlabel('Frequence [Hz]')
    plt.ylabel('Dephasage [rad]')
    plt.title('Dephasage courant $I_1$ (reference Us)')
    plt.show()
    return w, eta