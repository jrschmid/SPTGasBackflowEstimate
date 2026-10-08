##
## jschmidt Wed Sep  2 06:55:34 EDT 2026
##

import numpy as np
import matplotlib.pyplot as plt

###################################################
############### AUX ###############################
###################################################
def MeanFreePath(r0,n0,sigma,dd,rr):
# input:
# r0 is the linear dimension of the source in meters
# (radius for circular source, width for slit)
# d=2: expansion from circular source
# d=1: expansion from linear slit
# rr is the vertical distance from the source
# n0: gas number density at outlet
# sigma: collisional cross section of a water molecule

    # mean free path directly at the source (outlet)
    lambda0=1e0/np.sqrt(2e0)/n0/sigma
    # use geometric expansion:
    return lambda0*(rr/r0)**dd

###################################################
def FreezeLength(r0, n0, sigma, dd, MM, gamma):
# input:
# r0 is the linear dimension of the source in meters
# (radius for circular source, width for slit)
# d=2: expansion from circular source
# d=1: expansion from linear slit
# n0: gas number density at outlet
# sigma: collisional cross section of a water molecule
# MM: Mach number
# gamma: adiabatic index
    apu=1e0/(dd/2e0*(gamma+1e0)-1e0)
    aux=r0*sigma*n0/MM/dd
    return r0*aux**apu

###################################################
############### MAIN ##############################
###################################################

# number density: 1/10 of the triple
n0=1.6e22
# collisional cross section of a water molecule
sigma=4e-19
# linear dimension of outlet (radius for
# localized source, width for a slit
r0=1e0
# Mach number
MM=1e0
# adiabatic index for water
gamma=4e0/3e0

#
# plot settings
#
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["figure.dpi"] = 100
plt.rcParams["savefig.dpi"] = 300
# thicker axes and ticks
plt.rcParams["axes.linewidth"] = 2
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["ytick.major.width"] = 2
plt.rcParams["xtick.minor.width"] = 1.5
plt.rcParams["ytick.minor.width"] = 1.5
plt.rcParams["xtick.major.size"] = 7
plt.rcParams["ytick.major.size"] = 7
plt.rcParams["xtick.minor.size"] = 4
plt.rcParams["ytick.minor.size"] = 4
plt.rcParams["lines.linewidth"] = 2.5
# settings for text
plt.rcParams["font.size"] = 12
#plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["axes.labelsize"] = 14
plt.rcParams["axes.titlesize"] = 14

mode = 2

if mode==1:
    fig, p1 = plt.subplots(1)
    hh=np.logspace(3,6,100)
    p1.plot(hh/1e3,MeanFreePath(1e0,n0,sigma,2,hh),color='red',label='d=2 (circular gas source)')
    p1.set_xscale('log')
    p1.set_yscale('log')
    p1.plot(hh/1e3,MeanFreePath(1e0,n0,sigma,1,hh),color='blue',label='d=1 (linear gas source)')
    p1.set_xscale('log')
    p1.set(xlabel='ALTITUDE ABOVE SOURCE [km]')
    p1.set(ylabel='MEAN FREE PATH [m]')
    p1.set(title='GEOMETRIC EXPANSION (NO PRESSURE, NO GRAVITY)')
    hill_radius_altitude=960-252
    p1.axvline(hill_radius_altitude, color='k', linestyle='--',linewidth=1.)
    p1.text(hill_radius_altitude*1.1,.1,'HILL RADIUS',rotation=90, fontweight='bold')
    spt_radius_altitude=130
    p1.axvline(spt_radius_altitude, color='k', linestyle='--',linewidth=1.)
    p1.text(spt_radius_altitude*1.1,.1,'SPT DIAM.',rotation=90, fontweight='bold')
    p1.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98),prop={"weight": "bold"})
    p1.text(0.04, 0.78, r"$\mathbf{\lambda=\frac{1}{\sqrt{2}\,n_0\,\sigma}\left(\frac{h}{r_0}\right)^d}$", transform=p1.transAxes, ha="left", va="top")

elif mode==2:

    fig, p1 = plt.subplots(1)
    r0= np.logspace(-2, 1, 100)
    p1.plot(r0, FreezeLength(r0, n0, sigma, 2, MM, gamma),
            color='red', label='d=2 (circular gas source), M=1')
    p1.plot(r0, FreezeLength(r0, n0, sigma, 2, 4e0, gamma),
            color='red', label='d=2 (circular gas source), M=4', linestyle='--')
    p1.plot(r0, FreezeLength(r0, n0, sigma, 1, MM, gamma),
            color='blue', label='d=1 (linear gas source), M=1')
    p1.set_xscale('log')
    p1.set_yscale('log')
    p1.set(xlabel='OUTLET SIZE '+r"$\mathbf{r_0}$"+ ' [m]')
    p1.set(ylabel='COLLISIONAL DECOUPLING LENGTH [m]')
    p1.set(title='GEOMETRIC, ADIABATIC EXPANSION (NO PRESSURE, NO GRAVITY)')
    p1.legend(loc="upper left", bbox_to_anchor=(0.02, 0.98),prop={"weight": "bold"})
    p1.text(0.04, 0.68, r"$\mathbf{L=r_0\left[\frac{r_0\,\sigma\,n_0}{d\,M}\right]^{\left(\frac {d} {2} (\gamma +1)-1\right)^{-1}}}$", transform=p1.transAxes,
        ha="left", va="top")
#p1.xaxis.set_major_locator(MultipleLocator(2))
#p1.xaxis.set_minor_locator(MultipleLocator(1))
#p1.axhline(80e3, color='k', linestyle='--',linewidth=0.5)
#idx = np.where(np.diff(np.sign(phitab)) != 0)[0][0]
#equinox=yeartab[idx]

# bold tick labels (there is no rcParam for this)
plt.setp(p1.get_xticklabels(which='both') + p1.get_yticklabels(which='both'), fontweight='bold')

fig.subplots_adjust(bottom=0.15)
plt.figtext(0.04, 0.015,
            str(__file__) + " (mode:" + str(mode) + ")",
#        str(__file__) ,
             ha="left", fontsize=6, color="gray")

plt.show()




