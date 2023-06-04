#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 11:22:54 2022

@author: ravikularaman
"""

import numpy as np
from math import pi
from ionization_tools import integrate, velocity


m_p = 938.272 # in MeV, proton mass
m_e = 0.511  # MeV
M_sol = 2e30 # kg
M_gcr = 4.4e7*M_sol
m_H = 1.67e-27 # kg
m_avg = 1.4*m_H
D_gcr = 8.5*(3e21) # cm


#%% Extract cross-section data

p64_data = np.genfromtxt(open("sigma_6p4_p.txt", "r"), dtype=None)
E_p_list = p64_data[:, 0]*(1e-6)
sigp64_list = p64_data[:, 1]

e64_data = np.genfromtxt(open("sigma_6p4_e.txt", "r"), dtype=None)
E_e_list = e64_data[:, 0]*(1e-6)
sige64_list = e64_data[:, 1]


def sigma_64_p(E):
    # take E in MeV
    i=0
    while E >= E_p_list[i] and i<len(p64_data)-1:
        i += 1
    m = (sigp64_list[i]-sigp64_list[i-1])/(E_p_list[i]-E_p_list[i-1])
    sig = sigp64_list[i-1]+(m*(E-E_p_list[i-1]))
    return sig


def sigma_64_e(E):
    # take E in MeV
    i=0
    while E >= E_e_list[i] and i<len(e64_data)-1:
        i += 1
    m = (sige64_list[i]-sige64_list[i-1])/(E_e_list[i]-E_e_list[i-1])
    sig = sige64_list[i-1]+(m*(E-E_e_list[i-1]))
    return sig


#%% Fe K Alpha emissions

def Flux_Fe_64_p(T_c, f_p, norm):
    def integrand(T):
        v,_,_ = velocity(T, m_p)
        return sigma_64_p(T)*v*f_p(T)
    integ = integrate(integrand, T_c, 1e9)
    return norm*integ


def Flux_Fe_64_e(T_c, f_e, norm):
    def integrand(T):
        v,_,_ = velocity(T, m_e)
        return sigma_64_e(T)*v*f_e(T)
    integ = integrate(integrand, T_c, 1e9)
    return norm*integ
