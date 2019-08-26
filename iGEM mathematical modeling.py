# -*- coding: utf-8 -*-
" iGEM 2019 Mathematical Modelling"

#Calling in modules
import numpy as np
from scipy.integrate import odeint
from math import e, sin, pi, cos
import matplotlib.pyplot as plt

#Defining Fuctions
def model(x,t):
    alpha_d = 0.15
    alpha_g = 0.1
    alpha_r = 0.06
    beta_d = 0.001
    P_d = 1e-12
    P_g = 1e-12
    P_r = 1e-12
    gamma = 0.1
    gamma_m = 0.001
    gamma_g = gamma_m
    kf_gc = 0.6
    kr_gc = 0.1
    kr_cp = 0.1
    kf_cp = 0.5
    md = x[0]
    d = x[1]
    Sg = x[2]
    Sr = x[3]
    Cg = x[4]
    Cr = x[5]
    dmddt = alpha_d*P_d - (gamma + gamma_m)*md
    dddt = beta_d*md - kf_gc*(Sg+Sr) + kr_gc*(Cg+Cr) - gamma*d
    dSgdt = alpha_g*P_g + kr_gc*Cg - kf_gc*d - (gamma+gamma_g)*Sg    
    dSrdt = alpha_r*P_r + kr_gc*Cr - kf_gc*d - (gamma+gamma_g)*Sr
    dCgdt = kf_gc*Sg*d - kr_gc*Cg- kf_cp*P_g - gamma*Cg    
    dCrdt = kf_gc*Sr*d - kr_gc*Cr + kf_cp*P_r - gamma*Cr
    return [dmddt, dddt, dSgdt, dSrdt, dCgdt, dCrdt]


#Calling in the functions and initial variable
x0 = [0, 0, 0, 0, 0, 0]
t = np.linspace(0, 100, 1000000)
x = odeint(model, x0, t)

md = np.array(x)[:,0]
d = np.array(x)[:,1]
Sg = np.array(x)[:,2]
Sr = np.array(x)[:,3]
Cg = np.array(x)[:,4]
Cr = np.array(x)[:,5]

plt.plot(t, md)
plt.plot(t, d)
plt.plot(t, Sg)
plt.plot(t, Sr)
plt.plot(t, Cg)
plt.plot(t, Cr)
plt.xlabel("Time/unit")
plt.ylabel("Variables/unit")
plt.title("Combined Plot")
plt.gca().legend(("md", "d", "Sg", "Sr", "Cg", "Cr"))
plt.show()

