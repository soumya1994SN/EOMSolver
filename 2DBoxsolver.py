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
def func(t, tot_function):
    Nof = 64
    Ndim = 3
    Nxdiv = 480
    Nydiv = 480
    N = Ndim*Nxdiv*Nydiv*Nof
    rhs = np.reshape(tot_function, N).astype(np.complex64)
    rho = np.reshape(rhs, (Nxdiv, Nydiv, Nof, Ndim)).astype(np.complex64)
    gradient0 = np.zeros(N).astype(np.complex64)
    gradientx = np.reshape(gradient0, (Nxdiv, Nydiv, Nof, Ndim)).astype(np.complex64)
    gradienty = np.reshape(gradient0, (Nxdiv, Nydiv, Nof, Ndim)).astype(np.complex64)
    vel_arr = np.reshape(np.zeros(Nof*3).astype(np.complex64), (Nof, 3))
    vel_arr1 = np.linspace(-1, 1, 32).astype(np.complex64)
    vel_arr2 = np.concatenate([vel_arr1[0:16], vel_arr1[16:32], vel_arr1[16:32], vel_arr1[0:16]])
    vel_arr3 = np.concatenate([np.sqrt(1-vel_arr1[0:16]**2), np.sqrt(1-vel_arr1[16:32]**2), -np.sqrt(1-vel_arr1[16:32]**2), -np.sqrt(1-vel_arr1[0:16]**2)])
    vel_arr[:, 0] = 1j
    vel_arr[0:64, 1] = vel_arr2
    vel_arr[0:64, 2] = vel_arr3
    Selfpart0 = np.einsum('ij,kli...->klj...', vel_arr, rho)
    Selfpart1 = -np.einsum('ij,klj...->kli...', vel_arr, Selfpart0)*(vel_arr1[1]-vel_arr1[0])
    tot_Ham = Selfpart1
    Commutator = np.cross(tot_Ham, rho)
    for p in range(480):
        for m in range(Nof):
            for n in range(Ndim):
                gradientx[:,p,m,n] = vel_arr[m,1]*(-1.0)*sp.diff(rho[:,p,m,n], period = 18)
                gradienty[p,:,m,n] = vel_arr[m,2]*(-1.0)*sp.diff(rho[p,:,m,n], period = 18)
    gradientx = np.reshape(gradientx, N).astype(np.complex64)
    gradienty = np.reshape(gradienty, N).astype(np.complex64)
    Commutator1 = np.reshape(Commutator, N).astype(np.complex64)
    rhs =  gradientx+Commutator1+gradienty
    return rhs
start = time.time()
Nof = 64
Ndim = 3
Nxdiv = 480
Nydiv = 480
Ntimediv = 80
N = Ndim*Nof*Nxdiv*Nydiv
x_arr = np.linspace(0.0, 18.0, Nxdiv)
y_arr = np.linspace(0.0, 18.0, Nydiv)
initial = np.reshape(np.zeros(N).astype(np.complex64), (Nxdiv, Nydiv, Nof, Ndim))
initial[:, :, 0:32, 2] =  3
initial[:, :, 32:64, 2] = -3
t = np.linspace(0.0, 8, Ntimediv)
x_arr = np.linspace(0.0, 18.0, Nxdiv)
initial[240, 240, 0:64, 0] = 1*10**(-7)
initial[240, 240, 0:64, 1] = -1*10**(-7)
initial2 = np.reshape(initial, N).astype(np.complex64)
Ntimediv = 80
delta_t = t[1]-t[0]
solver = spi.ode(func)
solver.set_integrator(name = 'zvode', method = 'BDF', nsteps = 500000000, atol = 1e-9, rtol = 1e-9)
arr1 = []
solver.set_initial_value(initial2, 0.0)
for i in range(20):
    solver.integrate(solver.t+delta_t)
    arr1.append(np.reshape(solver.y, (Nxdiv, Nydiv, 64, 3))[:, :, :, :])
    np.save("2+1Dsimplified5NEW7752NEW100excel.npy", np.array(arr1))
#tim3 = time.time() - start
#arr9 = np.array(["--- %s seconds ---" %(tim3)])
#np.save("2+1Dtimesimplified5NEW7752NEW100.npy", arr9)
 

