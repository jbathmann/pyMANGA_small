#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 13:34:53 2018

@author: bathmann
"""

import numpy as np
from matplotlib import pyplot as plt

Diff= 0#2e-9
poro = 1e-2
w_0=0.035
kappa=1e-11
mu=1e-3
dp = -1e10
L = 10.


def v_darcy(kappa,mu,dp):
    return -kappa/mu*1000.*9.81*dp

print(v_darcy(kappa,mu,dp))

for dp in np.logspace(-5,5,11):
    print("dp: ", dp)
    Dh= poro*Diff+.5*np.sqrt(v_darcy(kappa,mu,-dp)**2)
    print(Dh)
    Kstar = - v_darcy(kappa,mu,-dp)/Dh
    print("kstar = ", Kstar)
    prefac = w_0*L/(L-1/Kstar*(1-np.exp(-Kstar*L)))
    print(prefac)
coordsX = np.linspace(0,L,1001)
w_x = prefac*(1-np.exp(-Kstar*coordsX))

plt.plot(coordsX,w_x)