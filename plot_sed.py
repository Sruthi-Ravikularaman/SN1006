import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc("text",usetex=True)

if __name__ == "__main__":

    hP=4.1357e-15;# eV s

    fs=24
    
    fig=plt.figure(figsize=(10, 8))
    ax=plt.subplot(111)

    alpha_p=2.4
    alpha_e=2.4
    ep_p=0.05
    ep_e=0.01
    Ec_e=18.0e12
    B0=18.0

    tag_ep_p="_epp="+str(int(ep_p*1.0e3))
    tag_ep_e="_epe="+str(int(ep_e*1.0e3))
    tag_alpha_p="_alp="+str(int(alpha_p*100.0))
    tag_alpha_e="_ale="+str(int(alpha_e*100.0))
    tag_Ec_e="_Ec="+str(int(Ec_e*1.0e-12))+"TeV"
    tag_B0="_B0="+str(int(B0))

    filename="Result/gamma_SN1006"+tag_ep_p+tag_alpha_p+tag_ep_e+tag_alpha_e+tag_Ec_e+tag_B0+"_SW.dat"

    Eg, phi_SYN, phi_ICS, phi_PPI, phi_cl=np.loadtxt(filename,unpack=True,usecols=[0, 1, 2, 3, 4])
    phi_PPI=np.nan_to_num(phi_PPI)
    phi_cl=np.nan_to_num(phi_cl)
    ax.plot(Eg,phi_SYN/6.242e+11,'y-.',label=r'$\textrm{SYN}$')
    ax.plot(Eg,phi_ICS/6.242e+11,'b--',label=r'$\textrm{ICS}$')
    ax.plot(Eg,phi_PPI/6.242e+11,'r-.',label=r'$\textrm{Shell}$')
    ax.plot(Eg,phi_cl/6.242e+11,'g:',label=r'$\textrm{Cloud}$')
    ax.plot(Eg,(phi_ICS+phi_PPI+phi_SYN+phi_cl)/6.242e+11,'k-')

    A=6.0e-6
    ax.plot(Eg,Eg*Eg*A*(Eg/1.0e3)**(-3)/6.242e+11,'--',color="purple")

    Eg, phi=np.loadtxt("Data/SN1006_data_SW_upper_FERMI_Acero2015.dat",unpack=True,usecols=[0, 1])
    ax.errorbar(Eg,phi,yerr=0.5*phi,linestyle='none',marker='_',color="green",uplims=True)

    Eg, phi_upper, phi_lower=np.loadtxt("Data/SN1006_data_SW_HESS.dat",unpack=True,usecols=[0, 1, 2])
    ax.fill_between(Eg,phi_lower,phi_upper,color=[(236.0/255.0,92.0/255.0,92.0/255.0)])

    Eg, Eg_lower, Eg_upper, phi, phi_upper, phi_lower=np.loadtxt("Data/SN1006_data_SW_FERMI.dat",unpack=True,usecols=[0, 1, 2, 3, 4, 5])
    Eg*=1.0e9
    Eg_upper*=1.0e9
    Eg_lower*=1.0e9
    ax.errorbar(Eg, phi, yerr=[phi-phi_lower,phi_upper-phi], xerr=[Eg-Eg_lower,Eg_upper-Eg], fmt='ko')

    Eg, phi=np.loadtxt("Data/SN1006_data_radio.dat",unpack=True,usecols=[0, 1])
    phi*=1.0/Eg
    Eg*=hP
    phi*=Eg*6.242e-12/hP
    ax.plot(Eg,phi/6.242e+11,'y^') 

    Eg, phi=np.loadtxt("Data/SN1006_data_Suzaku.dat",unpack=True,usecols=[0, 1])
    ax.plot(Eg,phi/6.242e+11,'bs') 

    Eg_test=np.array([3.0e9,9.0e9])
    phi_test=np.array([5.0e-13,3.5e-13])
    ax.plot(Eg_test,phi_test,'go') 

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel(r'$E_{\gamma} \textrm{ (eV)}$',fontsize=fs)
    ax.set_ylabel(r'$E^2_\gamma\phi(E_\gamma) \textrm{ (erg cm$^{-2}$ s$^{-1}$)}$',fontsize=fs)
    for label_axi in (ax.get_xticklabels() + ax.get_yticklabels()):
        label_axi.set_fontsize(fs)
    ax.set_xlim(1.0e-7,1.0e14)
    ax.set_ylim(1.0e-15,1.0e-10)
    # ax.set_xlim(1.0e2,1.0e5)
    # ax.set_ylim(1.0e-12,1.0e-11)
    ax.legend(loc='lower left', prop={"size":22})
    ax.grid(linestyle='--')

    plt.savefig("Result/SN1006"+tag_ep_p+tag_alpha_p+tag_ep_e+tag_alpha_e+tag_Ec_e+tag_B0+"_gamma.png")

    fig=plt.figure(figsize=(10, 8))
    ax=plt.subplot(111)


    fig=plt.figure(figsize=(10, 8))
    ax=plt.subplot(111)

    E, Jp_p, Jp_e, Jp_Vp, Jp_Ve, Jp_W28=np.loadtxt("Result/spectrum.dat",unpack=True,usecols=[0, 1, 2, 3, 4, 5])

    ax.plot(E,E**alpha_p*Jp_p,'g--',label=r'$\textrm{Proton SN 1006}$') 
    ax.plot(E,E**alpha_e*Jp_e,'r:',label=r'$\textrm{Electron SN 1006}$') 
    ax.plot(E,E**alpha_p*Jp_W28,'k-',label=r'$\textrm{SNR W28}$') 

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel(r'$T_i \textrm{ (eV)}$',fontsize=fs)
    ax.set_ylabel(r'$T_{i}^{\alpha}J_{d}(T_{i}) \textrm{ (eV$^{\alpha-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$)}$',fontsize=fs)
    for label_axi in (ax.get_xticklabels() + ax.get_yticklabels()):
        label_axi.set_fontsize(fs)
    ax.set_xlim(1.0e5,1.0e11)
    ax.set_ylim(1.0e10,2.0e16)
    ax.legend(loc='lower left', prop={"size":22})
    ax.grid(linestyle='--')

    plt.savefig("Result/SN1006_spectrum.png")


    fig=plt.figure(figsize=(10, 8))
    ax=plt.subplot(111)

    Ep, sigma_p=np.loadtxt("Result/diff_6p4"+tag_alpha_e+"_p.dat",unpack=True,usecols=[0, 1])
    Ee, sigma_e=np.loadtxt("Result/diff_6p4"+tag_alpha_e+"_e.dat",unpack=True,usecols=[0, 1])

    ax.plot(Ep,sigma_p*1.0e24,'g--',label=r'$\textrm{Proton}$',linewidth=2) 
    ax.plot(Ee,sigma_e*1.0e24,'r:',label=r'$\textrm{Electron}$',linewidth=2) 

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel(r'$T_i \textrm{ (eV)}$',fontsize=fs)
    ax.set_ylabel(r'$\sigma_{K\alpha}(T_i)\textrm{ (b H atom$^{-1}$)}$',fontsize=fs)
    for label_axi in (ax.get_xticklabels() + ax.get_yticklabels()):
        label_axi.set_fontsize(fs)
    ax.set_xlim(1.0e4,1.0e10)
    ax.set_ylim(5.0e-4,2.0e-2)
    ax.legend(loc='upper right', prop={"size":22})
    ax.grid(linestyle='--')

    plt.savefig("Result/SN1006_sigma_6p4.png")

    fig=plt.figure(figsize=(10, 8))
    ax=plt.subplot(111)

    ax.fill_between(E,2.13e-7*E**0,4.33e-7*E**0,color=[(255.0/255.0,102.0/255.0,102.0/255.0)],label=r'$\textrm{Observed Excess}$')

    Ep, diff_p=np.loadtxt("Result/diff_6p4_ale=220_p.dat",unpack=True,usecols=[0, 2])
    Ee, diff_e=np.loadtxt("Result/diff_6p4_ale=220_e.dat",unpack=True,usecols=[0, 2])

    ax.plot(Ep,diff_p,'g--',label=r'$\alpha=2.2$',linewidth=2) 
    ax.plot(Ee,diff_e,'g:',linewidth=2) 

    Ep, diff_p=np.loadtxt("Result/diff_6p4_ale=240_p.dat",unpack=True,usecols=[0, 2])
    Ee, diff_e=np.loadtxt("Result/diff_6p4_ale=240_e.dat",unpack=True,usecols=[0, 2])

    ax.plot(Ep,diff_p,'r--',label=r'$\alpha=2.4$',linewidth=2) 
    ax.plot(Ee,diff_e,'r:',linewidth=2) 

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend()
    ax.set_xlabel(r'$T_i \textrm{ (eV)}$',fontsize=fs)
    ax.set_ylabel(r'$T_i\sigma_{K\alpha}(T_i)J_d(T_i)M_{cl}/d_s^2  \textrm{ (photon cm$^{-2}$ s$^{-1}$)}$',fontsize=fs)
    for label_axi in (ax.get_xticklabels() + ax.get_yticklabels()):
        label_axi.set_fontsize(fs)
    ax.set_xlim(1.0e4,1.0e10)
    ax.set_ylim(1.0e-12,1.0e-6)
    ax.legend(loc='upper right', prop={"size":22})
    ax.grid(linestyle='--')

    plt.savefig("Result/SN1006_diff_sigma_6p4.png")

