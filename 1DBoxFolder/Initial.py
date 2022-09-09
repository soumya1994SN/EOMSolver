import numpy as np
#Here we define the system configuration calling values from "Configurations1Dbox.txt" file
Configurations = np.genfromtxt('Configurations1Dbox.txt')[1, :]
zmin = int(Configurations[0])
zmax = int(Configurations[1])
Nsp = int(Configurations[2])
vmin = int(Configurations[3])
vmax = int(Configurations[4])
Nvel = int(Configurations[5])
interval = (zmax-zmin)*0.5
def ELN(v): #Definition of the ELN distribution
    return 6*v
def Strans0ini(z, v): # e^1 component of S_v[z, t] as a function of z and v at t=0 
    return 10**(-6)*np.exp(-(z-interval)**2*(0.5)*(1/(10**(-6))))
def Strans1ini(z, v): # e^2 component of S_v[z, t] as a function of z and v at t=0 
    return (-1)*10**(-6)*np.exp(-(z-interval)**2*(0.5)*(1/(10**(-6))))
def Sparaini(z, v): # e^3 component of S_v[z, t] as a function of z and v at t=0 
    return 1
arr = np.reshape(np.zeros(Nsp*Nvel*3).astype(np.complex64), (Nsp, Nvel, 3)) #This tensor array stores the initial value of three components of the polarization vector at each z and v point
z_arr = np.linspace(zmin, zmax, Nsp)
v_arr = np.linspace(vmin, vmax, Nvel)
#The loops below insert the initial values into the tensor array defined as "arr"
for i in range(Nsp):
    for j in range(Nvel):
        arr[i, j, 0] = Strans0ini(z_arr[i], v_arr[j])*ELN(v_arr[j])
        arr[i, j, 1] = Strans1ini(z_arr[i], v_arr[j])*ELN(v_arr[j])
        arr[i, j, 2] = Sparaini(z_arr[i], v_arr[j])*ELN(v_arr[j])
np.save("Initialval.npy", np.reshape(arr, Nsp*Nvel*3)) 


