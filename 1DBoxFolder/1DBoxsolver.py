#Line No. 1-9 : Importing various Python libraries and solvers 
import numpy as np
from scipy.integrate import ode
from scipy.integrate import complex_ode
import scipy.integrate as spi
from scipy import fftpack as sp
from numpy import ma
import time
import math
#Line No. 11-25 : Defining the configuration of the system using input parameters from Configurations1Dbox.txt 
Configurations = np.genfromtxt('Configurations1Dbox.txt')[1, :]
xmin = Configurations[0]
xmax = Configurations[1]
Nsp = Configurations[2]
vmin = Configurations[3]
vmax = Configurations[4]
Nvel = Configurations[5]
tmin = Configurations[6]
tmax = Configurations[7]
Nt = Configurations[8]
N = Nsp*Nvel*3
x_arr = np.linspace(xmin, xmax, Nsp)
v = np.linspace(vmin, vmax, Nvel)
t = np.linspace(tmin, tmax, Nt)
delta_t = t[1]-t[0]
Interval = xmax-xmin
#Line No. 28-44 : Defining the right hand side for the discretized version of Eq.(3.2.1) (after shifting the advective term to right hand side) in terms of time and the components of the polarization vector  
def func(t, Polarization_vector):   
    S = np.reshape(Polarization_vector, N).astype(np.complex64) #
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
InitialPolarization_vector = np.load('Initialval.npy', mmap = 'r')
solver = spi.ode(func)
solver.set_integrator(name = 'zvode', method = 'BDF', nsteps = 5000000000, atol = 1e-9, rtol = 1e-9)
arr1 = []
solver.set_initial_value(InitialPolarization_vector, t_min)
for i in range(N_t):
    solver.integrate(solver.t+delta_t)
    arr1.append(np.reshape(solver.y, (N_sp, N_vel, 3))[:, :, :])
    np.save("LinearApt8.npy", np.array(arr1))

 



