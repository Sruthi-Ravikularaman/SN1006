#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 11:22:54 2022

@author: ravikularaman
"""

import numpy as np
from math import pi



m_p = 938.272 # in MeV, proton mass
m_e = 0.511  # MeV
M_sol = 2e30 # kg
M_gcr = 4.4e7*M_sol
m_H = 1.67e-27 # kg
m_avg = 1.4*m_H
D_gcr = 8.5*(3e21) # cm
c = 3e10


def velocity(E_kin, rest_mass):
    '''
    Gives the velocity, it's ratio with c and Lorentz factor.

    Parameters
    ----------
    E_kin : float
        Kinetic energy of particle in MeV.
    rest_mass : float
        Rest mass energy of particle in MeV.

    Returns
    -------
    v : float
        Particle velocity in cm s^-1.
    beta : float
        v/c.
    gamma : float
        Lorentz factor.

    '''
    E_kin = E_kin
    rest_mass = rest_mass
    gamma = 1+(E_kin/rest_mass)
    beta = np.sqrt(1-(1/(gamma**2)))
    v = c*beta
    return v, beta, gamma

def integrate(f, x_min, x_max, N_pts = 500):
    g = np.vectorize(f)
    if x_min==x_max:
        int_value=0
    else:
        if x_min == 0:
            z_min = np.log10(1e-50)
        else:
            z_min = np.log10(x_min)
        z_max = np.log10(x_max)
        int_range = np.log(np.logspace(z_min, z_max, N_pts))
        mid_points = (int_range[:-1]+int_range[1:])/2
        y_mid_points = np.exp(mid_points)*g(np.exp(mid_points))
        int_value = np.sum(y_mid_points*np.diff(int_range))
    return int_value


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


def Phi_Fe_64_p(T_p_i, J_CRp, del_p, Ap, Tp_c, M_cl, eta_Fe):
    def integrand(T):
        v,_,_ = velocity(T, m_p)
        return sigma_64_p(T) * v * J_CRp(T, del_p, Ap, Tp_c)
    integ = integrate(integrand, T_p_i, 1e12)
    norm = M_cl * eta_Fe / m_avg
    return norm*integ


def Phi_Fe_64_e(T_e_i, J_CRe, del_e, Ae, Te_c, M_cl, eta_Fe):
    def integrand(T):
        v,_,_ = velocity(T, m_e)
        return sigma_64_e(T) * v * J_CRe(T, del_e, Ae, Te_c)
    integ = integrate(integrand, T_e_i, 1e12)
    norm = M_cl * eta_Fe / m_avg
    return norm*integ
