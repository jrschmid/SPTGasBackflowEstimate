from scipy.special import erfc

import numpy as np
import matplotlib.pyplot as plt
import math

def Down_Flow_Fraction(uu):
    return (np.exp(-uu*uu)*uu/np.sqrt(np.pi)+erfc(uu)/2e0)

def Vtherm_func(temp):
    kb = 1.38064852e-23
    mh2o=2.99e-26
    return np.sqrt(2e0*kb*temp/mh2o)

logplot=False
# logplot=True

if logplot:
    ugas = np.logspace(2, 3, 100)
else:
    ugas=np.linspace(1,1000,100)

fig, p1 = plt.subplots(1,1,layout="constrained")
temp=200e0
p1.plot(ugas,100*Down_Flow_Fraction(ugas/Vtherm_func(temp)),color="red",label='T=200K')
temp=240e0
p1.plot(ugas,100*Down_Flow_Fraction(ugas/Vtherm_func(temp)),color="blue",label='T=240K')
if logplot: plt.xscale("log")
plt.xlabel("GAS SPEED [m/s]")
plt.ylabel("FRACTION OF GAS BACKFLOW TO SPT [%]")
plt.title("ESTIMATE WITH MAXWELL-BOLTZMAN, NO GRAVITY")
p1.axvline(x=240,color="black",linestyle="--")
p1.text(240,0,'escape speed',rotation=90)
plt.legend()
plt.show()