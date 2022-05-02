import os  
import numpy as np
import matplotlib.pyplot as plt
import time
import glob
import dpdata

from monty.serialization import loadfn,dumpfn
# ==================================== #
#  Constants
# ==================================== #

e = 1.602e-19                           #charge of a electron, C / electron volt, J
epsilon0 = 8.854e-12                    #vaccum dielectric constant, F/m
me = 9.109e-31                          #mass of a electron , kg
h = 6.626e-34                           #Planck constant, J*s 
hbar = h/(2*np.pi)
NA = 6.022e23
kb = 1.38e-23                            # Boltzmann constant, J/K
aB = 4*np.pi*epsilon0*hbar**2 / (me*e**2)# Bohr radius, m
Ry = 13.6                                # Ryberg, eV

# ==================================== #
#  Metal Units <---> S.I. Units
# ==================================== #

J2eV = 6.2442e18
m2A = 1e10
gmol2g = 1/NA
kg2gmol = 1e3/gmol2g
s2ps = 1e12
s2fs = 1e15

cm2m = 1e-2
gcc2kgm3 = 1e-3/(cm2m**3)
bar2Pa = 1e5
kbar2GPa = 1e3*bar2Pa/1e9

# ==================================== #
#  Atomic Units <---> Metal/S.I. Units
# ==================================== #

m2bohr = 1/aB
A2bohr = m2bohr/m2A
eV2Ry = 1/Ry

# ==================================== #
#  PHONON
# ==================================== #

c = 3*10**8                               # speed of light, m/s
omega2k = s2ps*1/(c) *cm2m                # 1/ps -> 1/cm

# ==================================== #
#  Basic Function 
# ==================================== #

def _estimate_sel(rcut, num_rho):
    return 4/3 * np.pi * rcut**3 * num_rho

def _plot_phonon(FILE, color, ll, lw, mk, INTERVAL, label):
    
    data = np.loadtxt(FILE)
    nbranch = data.shape[1] - 1
    
    k_max = data[-1,0]

    for i in range(1,nbranch+1):

        inver_cm2THz = 1/omega2k
        inver_cm2meV = 1/omega2k*s2ps*hbar*J2eV*1e3*2*np.pi
        #INTERVAL = 5
        if i == 1:
            ax.plot(data[::INTERVAL,0]/k_max,data[::INTERVAL,i]*inver_cm2THz,ls=ll,marker=mk, 
                    markersize=1.5, mew=0.5, color=color, linewidth=lw,  label=label)

        else:
            ax.plot(data[::INTERVAL,0]/k_max,data[::INTERVAL,i]*inver_cm2THz,ls=ll,marker=mk, 
                    markersize=1.5, mew=0.5, color=color, linewidth=lw, )        
              
    k_path = data[::100,0]/k_max
    
    #ax.set_xlim(data[0,0],data[-1,0])
        
    return k_path      

# ==================================== #
#  Basic For Plot
# ==================================== #

sci_color = np.array(['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e'])

def _gen_colormap(N, colormap=plt.cm.rainbow):

    return [colormap(int(x*colormap.N/N)) for x in range(N)]   
