import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mline
from scipy.integrate import ode
from scipy.integrate import complex_ode
import scipy.integrate as spi
from scipy import fftpack as sp
from numpy import ma
from matplotlib import rc
import matplotlib.lines as mline
import time
import math
def func(t, tot_function):
    Nof = 128
    Ndim = 3
    Nxdiv = 2**(12)
    N = Ndim*Nxdiv*Nof
    rhs = np.reshape(tot_function, N).astype(np.complex64)
    rho = np.reshape(rhs, (Nxdiv, Nof, Ndim)).astype(np.complex64)
    gradient0 = np.zeros(N).astype(np.complex64)
    gradient1 = np.reshape(gradient0, (Nxdiv, Nof, Ndim)).astype(np.complex64)
    vel_arr = np.reshape(np.zeros(Nof*2).astype(np.complex64), (Nof, 2 ))
    vel_arr1 = np.linspace(-1, 1, 128).astype(np.complex64)
    for i in range(128):
        vel_arr[i, :] = vel_arr1[i]
    vel_arr[:, 0] = 1j
    Selfpart0 = np.einsum('ij,ki...->kj...', vel_arr, rho)
    Selfpart1 = -np.einsum('ij,kj...->ki...', vel_arr, Selfpart0)*(vel_arr1[1]-vel_arr1[0])
    tot_Ham = Selfpart1*5
    Commutator = np.cross(tot_Ham, rho)
    for m in range(Nof):
        for n in range(Ndim):
            gradient1[:,m,n] = vel_arr[m,1]*(-1.0)*sp.diff(rho[:,m,n], period = 115)
    gradient2 = np.reshape(gradient1, N).astype(np.complex64)
    Commutator1 = np.reshape(Commutator, N).astype(np.complex64)
    rhs =  gradient2 + Commutator1
    return rhs
Nof = 128
Ndim = 3
a = 12
Nxdiv = 2**a
Ntimediv = 500
N = Ndim*Nof*Nxdiv
L = 2**(a-1)
x_arr = np.linspace(0.0, 115, Nxdiv)
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
#    print(arr1)
#tim3 = time.time() - start
#arr9 = np.array(["--- %s seconds ---" %(start), "--- %s seconds ---" %(tim1), "--- %s seconds ---" %(tim2), "--- %s seconds ---" %(tim3)])
#np.save("10000Continuousmodesabsnew7time100.npy", arr9)
 


