import numpy as np
from scipy.integrate import ode
from scipy.integrate import complex_ode
import scipy.integrate as spi
from scipy import fftpack as sp
from numpy import ma
import time
import math
Configurations = np.genfromtxt('Configurations1Dbox.txt')[1, :]
x_min = Configurations[0]
x_max = Configurations[1]
N_sp = Configurations[2]
v_min = Configurations[3]
v_max = Configurations[4]
N_vel = Configurations[5]
t_min = Configurations[6]
t_max = Configurations[7]
N_t = Configurations[8]
N = N_sp*N_vel*3
x_arr = np.linspace(x_min, x_max, N_sp)
v_x = np.linspace(v_min, v_max, N_vel)
t = np.linspace(t_min, t_max, N_t)
delta_t = t[1]-t[0]
def func(t, Polarization_vector):   
    Interval = x_max-x_min
    S = np.reshape(Polarization_vector, N).astype(np.complex64)
    Stensor = np.reshape(S, (N_sp, N_vel, 3)).astype(np.complex64)
    Advection = np.zeros(N).astype(np.complex64)
    Advectiontensor = np.reshape(Advection, (N_sp, N_vel, 3)).astype(np.complex64)
    veltensor = np.reshape(np.zeros(N_vel*2).astype(np.complex64), (N_vel, 2 ))
    for i in range(N_vel):
        veltensor[i, :] = v_x[i]
    veltensor[:, 0] = 1j
    Hamiltonian = np.einsum('ij,ki...->kj...', veltensor, Stensor)
    Hamiltoniantensor = -np.einsum('ij,kj...->ki...', veltensor, Stensor)*(v_x[1]-v_x[0])
    Crossproduct = np.cross(Hamiltoniantensor, Stensor)
    for m in range(N_vel):
        for n in range(3):
            Advectiontensor[:,m,n] = veltensor[m,1]*(-1.0)*sp.diff(Stensor[:,m,n], period = Interval)
        RHS = np.reshape(Advectiontensor+Crossproduct, N).astype(np.complex64)
    return RHS
InitialPolarization_vector = np.load('.npy', mmap = 'r')
solver = spi.ode(func)
solver.set_integrator(name = 'zvode', method = 'BDF', nsteps = 5000000000, atol = 1e-9, rtol = 1e-9)
arr1 = []
solver.set_initial_value(InitialPolarization_vector, t_min)
for i in range(N_t):
    solver.integrate(solver.t+delta_t)
    arr1.append(np.reshape(solver.y, (N_sp, N_vel, 3))[:, :, :])
    np.save("LinearApt8.npy", np.array(arr1))

 


