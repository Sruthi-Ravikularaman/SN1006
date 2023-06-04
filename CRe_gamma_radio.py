import numpy as np
from math import pi, sqrt, exp
from scipy.integrate import quad

#%% Constants

c = 3e10 # in cm/s, light velocity
m_p = 938.272 # in MeV, proton mass
m_e = 0.511  # in MeV
n_H = 1 # in cm-3, proton density target
B_0 = 1e-3 #G so cm⁻1/2 g1/2 s-1
e = 4.797e-10 # in cm^3/2 g^1/2 s-1
h_bar = 6.582e-22 # MeV s
e3B_0_MeV = 4.3e-20 # MeV2
M_sol = 2e30 # kg
m_H = 1.67e-27 #u.kg
m_avg = 1.4*m_H
alpha = 1/137
h = 4.1357e-21 # MeV s
h_cgs = 6.6261e-27 # cm2 g s-1
h_bar_cgs = 1.0546e-27 #cm2 g s-1
k_B = 8.6173e-11 # MeV K-1
k_B_cgs = 1.3807e-16 # cm2 g s-2 K-1
m_el = 9.1094e-28 #g
r_0 = e**2/(m_el*(c**2)) # cm electron radius


#%% Synchroton radiation

def syn_photon_E(T_e):
    T_e_TeV = T_e*(1e-6)
    return 0.02*(T_e_TeV**2)*1e-3 #in MeV


def E_c(T_e):
    return syn_photon_E(T_e)/0.29 #MeV


def R(x):
    num = 1.81*exp(-x)
    den_2 = pow(x, -2/3)+((3.62/pi)**2)
    return num/sqrt(den_2)


def Phi_syn(E, J_CRe, n):
    def integrand(T_e):
        x = E/E_c(T_e)
        J_val = J_CRe(T_e) # MeV-1 cm-3
        R_val = R(x)
        return J_val*R_val  #MeV-1 cm-3
    integ = quad(integrand, 1e2, 1e20) #cm-3
    pre = (sqrt(3)/(2*pi))*(e3B_0_MeV/m_e)*(1/(h_bar*E)) #MeV-1 s-1
    emi = pre * integ
    return emi #MeV-1 s-1 cm-3


#%% Non-thermal Bremsstrahlung emission

def sigma_scat_simple(E_gamma):
    m = 2e-24 # g
    tau_r = 66 # g/cm2
    return m/(tau_r*E_gamma)


def sigma_scat_diff(E_gamma, E_e):
    phi_1 = np.log((2*E_e/m_e)*(E_e-E_gamma)/E_gamma)-(1/2)
    phi_2 = -2*phi_1/3
    term_1 = (1+(1-(E_gamma/E_e))**2)*phi_1
    term_2 = (1-(E_gamma/E_e))*phi_2
    pre = 4*alpha*(r_0**2)
    return pre*(term_1+term_2)/E_gamma


def Phi_e_rel_brem(E, J_CRe, n):
    def integrand(T_e):
        sig = sigma_scat_diff(E, T_e) # in cm2 MeV-1
        F = J_CRe(T_e) # in MeV-1 cm-3
        return sig * F # in cm-1 MeV-2
    integ = quad(integrand, E, 1e9) # in MeV-1 cm-1
    emi = n * c * integ  # MeV-1 cm-3 s-1
    return emi


#%% Inverse Compton gamma emissivity

def G_3_0(x):
    c_3 = 0.319
    return ((pi**2)/6)*(1+c_3*x)*np.exp(-x)/(1+((pi**2)/6)*c_3*x)

def G_4_0(x):
    c_4 = 6.62
    return (pi**2/6)*(1+c_4*x)*np.exp(-x)/(1+(pi**2/6)*c_4*x)

def g_3(x):
    a_3 = 0.443
    b_3 = 0.54
    alpha_3 = 0.606
    beta_3 = 1.481
    g = a_3*np.power(x, alpha_3)/(1 + b_3*np.power(x, beta_3))
    return 1/(1+g)

def g_4(x):
    a_4 = 0.726
    b_4 = 0.382
    alpha_4 = 0.461
    beta_4 = 1.457
    g = a_4*np.power(x, alpha_4)/(1 + b_4*np.power(x, beta_4))
    return 1/(1+g)

def G_3(x):
    return G_3_0(x)*g_3(x)

def G_4(x):
    return G_4_0(x)*g_4(x)

def N_iso(E_ph, E_e, T, k_dil):
    "E_ph and E_e in units MeV, T in K"
    # Convert E_ph and E_e to units of m_e*c^2
    E_ph = E_ph/m_e
    E_e = E_e/m_e
    # Convert T to units of m_e*c^2
    T = k_B*T/m_e

    z = E_ph/E_e

    t = 4*E_e*T
    x_0 = z/((1-z)*t)

    tmp=(T/E_e)**2
    pre_1 = 2*(r_0**2)*(m_el**3)*(c**4)*k_dil
    pre_2 = pi*(h_bar_cgs**3)
    pre = pre_1/pre_2
    pre_tmp = pre*tmp

    term_1 = (z**2/(2*(1-z)))*G_3(x_0)
    term_2 = G_4(x_0)

    return pre_tmp*(term_1+term_2)/m_e # MeV-1 s-1

w_star = h_cgs*c/(1e-4) #g cm2 s-2
T_star = 0.255*w_star/k_B_cgs


def N_iso_ph_N(E_ph, J_CRe, T, k_dil):
    def integrand(E_e):
        return N_iso(E_ph, E_e, T, k_dil) * J_CRe(E_e) #  MeV-2 cm-3 s-1
    E_e_min = E_ph
    E_e_max = 1e12 #MeV
    return quad(integrand, E_e_min, E_e_max) # MeV-1 s-1


def Phi_e_IC(E, J_CRe, T, k_dil):
    N_iso_term = N_iso_ph_N(E, J_CRe, T, k_dil) # MeV-1 cm-3 s-1
    return N_iso_term # MeV-1 cm-3 s-1

