#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 16:01:07 2019

@author: prateek
calculate logarithmic derivatives of Lambda w.r.t. Z and T
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt

def Lambda(temp,ZbyZsun): #interpolated in T but appropriately scaled in Z
    klo=0; khi=tab_sz-1
    while (klo != (khi-1)):
        kmid = int((khi+klo)/2)
        Tmid = Ttab[kmid]
        if (temp<Tmid):
            khi = kmid
        if (temp>Tmid):
            klo = kmid
    dT = Ttab[khi] - Ttab[klo]
    Lambda_lo = L_HHe_tab[klo]+ZbyZsun*L_Z_tab[klo]
    Lambda_hi = L_HHe_tab[khi]+ZbyZsun*L_Z_tab[khi]
    scrh = Lambda_lo*(Ttab[khi]-temp)/dT + Lambda_hi*(temp-Ttab[klo])/dT; #linear interpolation
    return scrh

global Ttab, L_HHe_tab, L_Z_tab, tab_sz
D = np.loadtxt('CT_WSS09.dat')
Ttab = D[:,0]; L_HHe_tab = D[:,1]; L_Z_tab = D[:,2]; tab_sz=np.size(Ttab)
#now start with derivative calculation

nT = 60; nZ = 50
T = np.logspace(6,8.1,nT); Z = np.logspace(-1,0.6,nZ)
dlnLdlnT = np.zeros((nT,nZ)); dlnLdlnZ = np.zeros((nT,nZ))
for i in range(nT-1):
    for j in range(nZ-1):
        dlnLdlnT[i,j] = m.log(Lambda(T[i+1],Z[j])/Lambda(T[i],Z[j]))/m.log(T[i+1]/T[i])
        dlnLdlnZ[i,j] = m.log(Lambda(T[i],Z[j+1])/Lambda(T[i],Z[j]))/m.log(Z[j+1]/Z[j])
'''
plt.contourf(np.log10(T),np.log10(Z),np.transpose(dlnLdlnT),80)
plt.title(r'$d\ln \Lambda/d\ln T$')
plt.xlim(6,8); plt.ylim(-1,0.5)
plt.colorbar()
CS = plt.contour(np.log10(T),np.log10(Z),np.transpose(dlnLdlnT), 8,colors='white')
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
plt.xlabel(r'$\log_{10}T(K)$'); plt.ylabel(r'$\log_{10} Z/Z_\odot$')
plt.show()
'''
plt.contourf(np.log10(T),np.log10(Z),np.transpose(dlnLdlnZ),80)
plt.title(r'$d\ln \Lambda/d\ln Z$')
plt.xlim(6,8); plt.ylim(-1,0.5)
plt.colorbar()
CS = plt.contour(np.log10(T),np.log10(Z),np.transpose(dlnLdlnZ), 8,colors='white')
plt.clabel(CS, fontsize=12, inline=1, fmt='%1.1f')
plt.xlabel(r'$\log_{10}T(K)$'); plt.ylabel(r'$\log_{10} Z/Z_\odot$')
plt.show()