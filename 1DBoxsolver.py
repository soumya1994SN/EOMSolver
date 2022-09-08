import numpy as np
from scipy.integrate import ode
from scipy.integrate import complex_ode
import scipy.integrate as spi
from scipy import fftpack as sp
from numpy import ma
import time
import math
Configurations = np.genfromtxt('Configurations1Dbox.txt')[1, :]
x_ini = Configurations[0]
x_fin = Configurations[1]
N_sp = Configurations[2]
v_min = Configurations[3]
v_max = Configurations[4]
N_vel = Configurations[5]
t_min = Configurations[6]
t_max = Configurations[7]
N_t = Configurations[8]
N = N_sp*N_vel*3
def func(t, tot_function):   
    Interval = x_fin-x_ini
    rhs = np.reshape(tot_function, N).astype(np.complex64)
    rho = np.reshape(rhs, (N_sp, N_vel, 3)).astype(np.complex64)
    gradient0 = np.zeros(N).astype(np.complex64)
    gradient1 = np.reshape(gradient0, (N_sp, N_vel, 3)).astype(np.complex64)
    vel_arr = np.reshape(np.zeros(N_vel*2).astype(np.complex64), (N_vel, 2 ))
    vel_arr1 = np.linspace(-1, 1, 128).astype(np.complex64)
    for i in range(128):
        vel_arr[i, :] = vel_arr1[i]
    vel_arr[:, 0] = 1j
    Selfpart0 = np.einsum('ij,ki...->kj...', vel_arr, rho)
    Selfpart1 = -np.einsum('ij,kj...->ki...', vel_arr, Selfpart0)*(vel_arr1[1]-vel_arr1[0])
    tot_Ham = Selfpart1
    Commutator = np.cross(tot_Ham, rho)
    for m in range(N_vel):
        for n in range(3):
            gradient1[:,m,n] = vel_arr[m,1]*(-1.0)*sp.diff(rho[:,m,n], period = Interval)
    gradient2 = np.reshape(gradient1, N).astype(np.complex64)
    Commutator1 = np.reshape(Commutator, N).astype(np.complex64)
    rhs =  gradient2 + Commutator1
    return rhs
x_arr = np.linspace(x_ini, x_fin, N_sp)
initial = np.reshape(np.zeros(N).astype(np.complex64), (Nxdiv, Nof, Ndim))
vel_arr = np.linspace(-1, 1, 128).astype(np.complex64)
for i in range(Nxdiv):
    initial[i, 0:64, 2] = (2-2*0.8)*vel_arr[0:64]
    initial[i, 64:128, 2] = 2*vel_arr[64:128]
t = np.linspace(0.0, 50.0, Ntimediv)
x_arr = np.linspace(0.0, 115, Nxdiv)
initial[L, :, 0] = 10**(-6)
initial[L, :, 1] = -10**(-6)
initial2 = np.reshape(initial, N).astype(np.complex64)
Ntimediv = 500
delta_t = t[1]-t[0]
solver = spi.ode(func)
solver.set_integrator(name = 'zvode', method = 'BDF', nsteps = 5000000000, atol = 1e-6, rtol = 1e-6)
#tim2 = time.time() - start
arr1 = []
solver.set_initial_value(initial2, 0.0)
for i in range(300):
    solver.integrate(solver.t+delta_t)
    arr1.append(np.reshape(solver.y, (Nxdiv, 128, 3))[:, :, :])
    np.save("LinearApt8.npy", np.array(arr1))

 


