#include <iostream>
#include <math.h>
#include "gato.h"

using namespace std;

const double pi=3.14159265359;
const double me=0.510998e6;// eV
const double mp=938.272e6;// eV

const double meCGS=9.10938356e-28;// g
const double sigma_sb=3.5394474508e7;// erg cm^-2 s^-1 K^-4
const double qe=4.80320451e-10;// cgs unit
const double hP=4.1357e-15;// eV s
const double Ktomec2=1.6863699549e-10;// -> To convert temperture T from K/kB to eV/me 
const double mpi=134.9766e6;// eV
const double Tpth=2.0*mpi+(pow(mpi,2)/(2.0*mp));

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cosmic-ray spectra
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Normalization factor for the source function.
double func_Gam(double mass, double alpha, double Emin, double Emax){
// mass (eV) -> mass of CR species,
// alpha -> index of CR spectrum,
// Emin (eV) -> minimum energy of CRs,
// Emax (eV) -> maximum energy of CRs.

    double Gam, p0=1.0e9, x, dx, xmin=sqrt(pow(Emin+mass,2)-mass*mass)/p0, xmax=sqrt(pow(Emax+mass,2)-mass*mass)/p0;
    Gam=0.0;
    x=xmin;
    while(x<xmax){
        dx=min(0.001*x,xmax-x);
        Gam+=pow(x,-alpha)*(sqrt(x*x+pow(mass/p0,2))-(mass/p0))*dx;

        x+=dx;
    }

    return Gam;
}

// Volume integrated spectrum of CR protons from the SNR.
double func_Qp(double xiSNR, double Gam, double alpha, double E){
// xiSNR -> acceleration efficiency, 
// Gam -> normalization factor,
// alpha -> index of CR spectrum,
// E (eV) -> CR energy.

    double ESNR, p0, p, vp;
    ESNR=1.0e51*6.242e+11;// eV
    p0=1.0e9;// eV
    p=sqrt(pow(E+mp,2)-pow(mp,2));// eV
    vp=p/(E+mp);
    
    double Q0, Q;
    Q0=xiSNR*ESNR/(pow(p0,2)*vp*Gam);
    Q=Q0*pow(p/p0,-alpha);
    
    return Q;// eV^-1 
}

// Volume integrated spectrum of CR electrons from the SNR.
double func_Qe(double xiSNR, double Gam, double alpha, double E){
// xiSNR -> acceleration efficiency, 
// Gam -> normalization factor,
// alpha -> index of CR spectrum,
// E (eV) -> CR energy.
    
    double ENSR, p0, p , vp;
    ENSR=1.0e51*6.242e+11;// eV
    p0=1.0e9;// eV
    p=sqrt(pow(E+me,2)-pow(me,2));// eV
    vp=p/(E+me);
    
    double Q0, Q;
    Q0=xiSNR*ENSR/(pow(p0,2)*vp*Gam);
    Q=Q0*pow(p/p0,-alpha);
    
    // Q(E) is the source term for the transport equation of f(E) -> unit of Q(E) should be eV^-1
    // Q(E)=xiSNR*ESNR*pow(p(E)/p0,2-delta)/((p0*c)^2*Gamma*beta)
    // Gamma=Integrate[x^(2-delta)*(sqrt(x^2+(m^2c^2/p0^2))-mc/p0),{x,p_min/p0,p_max/p0}]
    // The mass m should be adjusted appropriately for protons and electrons
    // Gamma_p=4.74166 -> delta=4.2, Emin=1.0e6, and Emax=5.0e15 eV
    // Gamma_e=17.548 -> delta=4.2, Emin=1.0e6 and Emax=1.0e14 eV
    
    return Q;// eV^-1
}

// Proton spectrum in SNR W28
double func_Jp_W28(double E){
// E (eV) -> CR energy.

    double Ap=3.15e-17, vp=sqrt(pow(E+mp,2)-mp*mp)*3.0e10/(E+mp), Jp;
    Jp=Ap*(E+mp)*1.0e-9*pow((E*E+2.0*E*mp)*1.0e-18,-0.5*(2.76+1.0));
    Jp*=vp/(4.0*pi);

    return Jp;// eV^-1 cm^-2 s^-1 sr^-1
}

// Proton spectrum form Voyager  
double func_Jp_Vp(double E){
// E (eV) -> CR energy.

    double Jp=1.882e-9*pow(E*1.0e-6,0.129)*pow(1.0+(E/624.5e6),-2.829);

    return Jp;// eV^-1 cm^-2 s^-1 sr^-1
}

// Electron spectrum form Voyager  
double func_Jp_Ve(double E){
// E (eV) -> CR energy.

    double Je=4.658e-7*pow(E*1.0e-6,-1.123)*pow(1.0+(E/736.2e6),-2.033);

    return Je;// eV^-1 cm^-2 s^-1 sr^-1
}

// CR spectrum for SN 1006
double func_Jp(double Anorm, double alpha, double mass, double Ec, double E){
// Anorm (eV^-1 cm^-3) -> normalization factor,
// alpha -> index of CR spectrum,
// mass (eV) -> mass of CR species,
// Ec (eV) -> cut-off energy of CR spectrum,
// E (eV) -> CR energy.

    double p0=1.0e11, p=sqrt(pow(E+mass,2)-mass*mass), Jp;
    Jp=Anorm*pow(p/p0,-alpha)*exp(-E/Ec);

    // Uncomment this if you want to try a broken power law.
    // double pbr=1.0e9, beta=-0.6;
    // Jp*=pow(pbr/p0,-beta)*pow(p/p0,beta)*pow(1.0+(p/pbr),-beta)*exp(-E/Ec);

    return Jp;// eV^-1 cm^-2 s^-1 sr^-1
}

// Function to calculate the total energy of CRs 
double func_ECR(double Anorm, double alpha, double mass, double Ec, double Emin, double Emax){
// Anorm (eV^-1 cm^-3) -> normalization factor,
// alpha -> index of CR spectrum,
// mass (eV) -> mass of CR species,
// Ec (eV) -> cut-off energy of CR spectrum,
// Emin (eV) -> minimum energy of CRs,
// Emax (eV) -> maximum energy of CRs.

    double E, dE, vp, ECR;
    E=Emin, ECR=0.0;
    while(E<Emax){
        dE=min(0.01*E,Emax-E);
        vp=sqrt(pow(E+mass,2)-mass*mass)*3.0e10/(E+mass);
        ECR+=dE*E*4.0*pi*func_Jp(Anorm,alpha,mass,Ec,E)/vp;

        E+=dE;
    }
    
    return ECR*1.60218e-12;// erg cm^-3
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gamma rays from proton-proton interaction -> Kafexhiu et al. 2014
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Dimensionless Breit-Wigner distribution.
double func_fBW(double sqrt_s){
// Eq. 4.

    double K, M_res, Gamma_res, gamma;
    M_res=1.1883*pow(10.0,9);// eV
    Gamma_res=0.2264*pow(10.0,9);// eV
    gamma=sqrt(pow(M_res,2)*(pow(M_res,2)+pow(Gamma_res,2)));
    K=sqrt(8.0)*M_res*Gamma_res*gamma/(pi*sqrt(pow(M_res,2)+gamma));
    
    double fBW=mp*K/(pow(pow(sqrt_s-mp,2)-pow(M_res,2),2)+pow(M_res*Gamma_res,2));
    
    return fBW;
}

// Cross-section for p+p -> p+p+pi0.
double func_sigma_1pi(double Tp){
// Eq. 2,
// Tp (eV) -> kinetic energy of CR proton.

    double sigma_0, sqrt_s, eta;
    sigma_0=7.66e-3;
    sqrt_s=sqrt(2.0*mp*(Tp+2.0*mp));
    eta=sqrt(pow(pow(sqrt_s,2)-pow(mpi,2)-4.0*pow(mp,2),2)-16.0*pow(mpi*mp,2))/(2.0*mpi*sqrt_s);
    
    double sigma_1pi=sigma_0*pow(eta,1.95)*(1.0+eta+pow(eta,5))*pow(func_fBW(sqrt_s),1.86);
    
    return sigma_1pi;// mb
}

// Cross-section for p+p -> p+p+2pi0.
double func_sigma_2pi(double Tp){
// Eq. 5,
// Tp (eV) -> kinetic energy of CR proton.
 
    double sigma_2pi;
    if(Tp<0.56*pow(10.0,9)){
        sigma_2pi=0.0;
    }
    if(Tp>=0.56*pow(10.0,9)){
        sigma_2pi=5.7/(1.0+exp(-9.3*(Tp*pow(10.0,-9)-1.4)));
    }
    
    return sigma_2pi;// mb
}

// Proton-proton total inelastic cross-section. 
double func_sigma_inel(double Tp){
// Eq. 1,
// Tp (eV) -> kinetic energy of CR proton.

    double sigma_inel;
    sigma_inel=30.7-0.96*log(Tp/Tpth)+0.18*pow(log(Tp/Tpth),2);
    sigma_inel*=pow(1.0-pow(Tpth/Tp,1.9),3);
    
    return sigma_inel;// mb
}

// Average pi0 multiplicity.
double func_npi(double Tp){
// Eq. 7 and Table 4 (GEANT 4 model),
// Tp (eV) -> kinetic energy of CR proton.

    double Qp, xip, npi, a1, a2, a3, a4, a5;
    a1=0.728;
    a2=0.596;
    a3=0.491;
    a4=0.2503;
    a5=0.117;
    
    Qp=(Tp-Tpth)/mp;
    xip=(Tp-3.0e9)/mp;
    
    if((Tp>=1.0e9) && (Tp<5.0e9)){
        npi=-6.0*pow(10.0,-3)+0.237*Qp-0.023*pow(Qp,2);
    }
    if(Tp>=5.0e9){
        npi=a1*pow(xip,a4)*(1.0+exp(-a2*pow(xip,a5)))*(1.0-exp(-a3*pow(xip,0.25)));
    }
    
    return npi;
}

// Pi0 production cross-section.
double func_sigma_pi(double Tp){
// See paragraph above Table 4,
// Tp (eV) -> kinetic energy of CR proton.

    double sigma_pi;
    
    if((Tp>=Tpth) && (Tp<2.0e9)){
        sigma_pi=func_sigma_1pi(Tp)+func_sigma_2pi(Tp);
    }
    if(Tp>=2.0e9){
        sigma_pi=func_sigma_inel(Tp)*func_npi(Tp);
    }
    
    return sigma_pi;// mb
}

// Complementary function for the differential cross-section.
double func_Amax(double Tp){
// Eq. 12 and Table 7 (GEANT 4 model),
// Tp (eV) -> kinetic energy of CR proton.

    double b0, b1, b2, b3, Amax, Epi_max, Epi_min, sqrt_s, gamma_CM, Epi_CM, Ppi_CM, beta_CM, theta_p;
    Amax=0.0;
    if(Tp>Tpth){
        b0=5.9;
        
        sqrt_s=sqrt(2.0*mp*(Tp+2.0*mp));
        gamma_CM=(Tp+2.0*mp)/sqrt_s;
        beta_CM=sqrt(1.0-pow(gamma_CM,-2));
        Epi_CM=(pow(sqrt_s,2)-4.0*pow(mp,2)+pow(mpi,2))/(2.0*sqrt_s);
        Ppi_CM=sqrt(pow(Epi_CM,2)-pow(mpi,2));
        Epi_max=gamma_CM*(Epi_CM+Ppi_CM*beta_CM);
        Epi_min=gamma_CM*(Epi_CM-Ppi_CM*beta_CM);
        theta_p=Tp/mp;
        
        if((Tp>=Tpth) && (Tp<1.0e9)){
            Amax+=b0*func_sigma_pi(Tp)/Epi_max;
        }
        if(Tp>=1.0e9){
            
            if((Tp>=1.0e9) && (Tp<5.0e9)){
                b1=9.53;
                b2=0.52;
                b3=0.054;
            }
            if(Tp>=5.0e9){
                b1=9.13;
                b2=0.35;
                b3=9.7e-3;
            }
            
            Amax+=b1*pow(theta_p,-b2)*exp(b3*pow(log(theta_p),2))*func_sigma_pi(Tp)/mp;
        }
    }
    
    return Amax;// mb/eV
}

// Complementary function for the differential cross-section.
double func_alpha(double Tp){
// Table 5, Eq. 14, and Eq. 15,
// Tp (eV) -> kinetic energy of CR proton.

    double alpha;
    if((Tp>=Tpth) && (Tp<=20.0e9)){
        alpha=1.0;
    }
    if(Tp>20.0e9){
        alpha=0.5;
    }
    
    return alpha;
}

// Complementary function for the differential cross-section.
double func_beta(double Tp){
// Table 5, Eq. 14, and Eq. 15,
// Tp (eV) -> kinetic energy of CR proton.

    double beta, mu, q, theta_p;
    q=(Tp-1.0e9)/mp;
    mu=1.25*pow(q,1.25)*exp(-1.25*q);
    theta_p=Tp/mp;
    
    if((Tp>=Tpth) && (Tp<=1.0e9)){
        beta=3.29-0.2*pow(theta_p,-1.5);
    }
    if((Tp>1.0e9) && (Tp<=4.0e9)){
        beta=mu+2.45;
    }
    if((Tp>4.0e9) && (Tp<=20.0e9)){
        beta=1.5*mu+4.95;
    }
    if((Tp>20.0e9) && (Tp<=100.0e9)){
        beta=4.2;
    }
    if(Tp>100.0e9){
        beta=4.9;
    }
    
    return beta;
}

// Complementary function for the differential cross-section.
double func_gamma(double Tp){
// Table 5, Eq. 14, and Eq. 15,
// Tp (eV) -> kinetic energy of CR proton.
    
    double gamma, mu, q;
    q=(Tp-1.0e9)/mp;
    mu=1.25*pow(q,1.25)*exp(-1.25*q);
    
    if((Tp>=Tpth) && (Tp<1.0e9)){
        gamma=0.0;
    }
    if((Tp>=1.0e9) && (Tp<=4.0e9)){
        gamma=mu+1.45;
    }
    if((Tp>4.0e9) && (Tp<=20.0e9)){
        gamma=mu+1.5;
    }
    if(Tp>20.0e9){
        gamma=1.0;
    }
    
    return gamma;
}

// Complementary function for the differential cross-section.
double func_F(double Tp, double Eg){
// Eq. 11 and Table 5,    
// Tp (eV) -> kinetic energy of CR proton,
// Eg (eV) -> energy of gamma ray.

    double Yg, Yg_max, Xg, sqrt_s, gamma_CM, gammapi_LAB, beta_CM, betapi_LAB, Epi_CM, Ppi_CM, Epi_max_LAB, F, C, Eg_max, Eg_min;
    sqrt_s=sqrt(2.0*mp*(Tp+2.0*mp));
    gamma_CM=(Tp+2.0*mp)/sqrt_s;
    beta_CM=sqrt(1.0-pow(gamma_CM,-2));
    Epi_CM=(pow(sqrt_s,2)-4.0*pow(mp,2)+pow(mpi,2))/(2.0*sqrt_s);
    Ppi_CM=sqrt(pow(Epi_CM,2)-pow(mpi,2));
    
    Epi_max_LAB=gamma_CM*(Epi_CM+Ppi_CM*beta_CM);
    gammapi_LAB=Epi_max_LAB/mpi;
    betapi_LAB=sqrt(1.0-pow(gammapi_LAB,-2));
    Eg_max=mpi*gammapi_LAB*(1.0+betapi_LAB)/2.0;
    Eg_min=mpi*gammapi_LAB*(1.0-betapi_LAB)/2.0;
    
    F=0.0;
    if((Eg>=Eg_min) && (Eg<=Eg_max)){
        Yg=Eg+pow(mpi,2)/(4.0*Eg);
        
        Yg_max=Eg_max+pow(mpi,2)/(4.0*Eg_max);// Yg_max=mpi*gammapi_LAB
        Xg=(Yg-mpi)/(Yg_max-mpi);
        C=3.0*mpi/Yg_max;

        F+=pow(1.0-pow(Xg,func_alpha(Tp)),func_beta(Tp))/pow(1.0+Xg/C,func_gamma(Tp));
    }
    
    return F;
}

// Complementary function for the nuclear enhancement factor.
double func_GTp(double Tp){
// Eq. 19,
// Tp (eV) -> kinetic energy of CR proton.

    double GTp;
    GTp=1.0+log(max(1.0,func_sigma_inel(Tp)/func_sigma_inel(pow(10.0,12))));
    
    return GTp;
}

// Nuclear enhancement factor. 
double func_enhancement(double Tp){
// Eq. 24,
// Tp (eV) -> kinetic energy of CR proton.

    double eps_nucl;
    eps_nucl=0.0;
    if(Tp>=Tpth){
        eps_nucl=1.7;
    }
    if(Tp>=pow(10.0,9)){
        eps_nucl=1.37+0.39*10.0*pi*func_GTp(Tp)/func_sigma_inel(Tp);
    }
    
    return eps_nucl;
}

// Differential cross-section of gamma-ray from pi0 decay.
double func_d_sigma_g(double Tp, double Eg){
// Eq. 8,  
// Tp (eV) -> kinetic energy of CR proton,
// Eg (eV) -> energy of gamma ray.

    double d_sigma_g;
    
    if(Tp>=Tpth){
        d_sigma_g=func_Amax(Tp)*func_F(Tp,Eg)*pow(10.0,-27);
    }
    if(Tp<Tpth){
        d_sigma_g=0.0;
    }
    
    return d_sigma_g;// cm^2/eV
}

// Volume emissitivity of hadronic gamma rays.
double func_eps_PPI(double nH, double Tp, double Eg){
// Notice that we need the hydrogen denisty but for clouds this wil cancel
// out and the spectrum scales only with M_cloud/d^2.
// nH (cm^-3) -> hydrogen density,
// Tp (eV) -> kinetic energy of CR proton,
// Eg (eV) -> energy of gamma ray.

    double vp, eps_PPI;
    vp=sqrt(pow(Tp+mp,2)-mp*mp)*3.0e10/(Tp+mp);
    eps_PPI=nH*vp*func_d_sigma_g(Tp,Eg);
    eps_PPI*=func_enhancement(Tp);

    return eps_PPI;// eV^-1 s^-1
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gamma rays from inverse Compton scattering -> Khangulyan et al. 2014
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Complementary functions for gamma rays from IC.
double func_G34(double x0, double* a){
// Eq. 20, Eq. 24, and Eq.25.

    double alphai=a[0], ai=a[1], betai=a[2], bi=a[3], ci=a[4], tmp, G, g;
    tmp=(1.0+ci*x0)/(1.0+(pi*pi*ci*x0/6.0));
    G=pi*pi*tmp*exp(-x0)/6.0;
    tmp=1.0+bi*pow(x0,betai);
    g=1.0/((ai*pow(x0,alphai)/tmp)+1.0);

    return G*g;
}

// Volume emissitivity of upscattered thermal photons by one electron (IC of a black-body or a gray-body photon distribution). 
double func_eps_ICS(double Tph, double kappa, double E, double Eg){
// Eq. 14, 
// Tph (K) -> soft photon temperature,
// kappa (eV/cm^3) -> photon energy density,
// E (eV) -> cosmic ray electron's energy,
// Eg (eV) -> gamma energy.

    // The dilution factor for a grey-body photon spectrum (ratio between energy density in the grey-body photons to the energy density
    // of black-body phtons of the same temperature). See also the definition in Eq. 29 of Khangulyan et al. 2014.
    kappa*=1.0/(4.0*sigma_sb*pow(Tph,4)/3.0e10);
    Tph*=Ktomec2, Eg*=1.0/me, E*=1.0/me;

    // Parameters from Eqs 26, 27
    double a3[5]={0.606, 0.443, 1.481, 0.540, 0.319};
    double a4[5]={0.461, 0.726, 1.457, 0.382, 6.620};
    
    double z, x0, eps_IC;
    z=Eg/E;
    x0=z/((1.0-z)*(4.0*E*Tph));

    // Eq. 14
    eps_IC=(pow(z,2)*func_G34(x0,a3)/(2.0*(1.0-z)))+func_G34(x0,a4);
    eps_IC*=2.6318735743809104e16*kappa*pow((Tph/E),2)/me;
    
    if((Eg>E) || (E<1.0)){
        eps_IC=0.0;
    }

    return eps_IC;// eV^-1 s^-1
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Radio and X-ray from synchrotron radiation -> Zirakashvili & Aharonian 2007
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Volume emissitivity of synchrotron photons. 
double func_eps_SYN(double B, double E, double Eg){
// Eq. 30,
// B (uG) -> magnetic field of the synchrotron emitter, 
// E (eV) -> kinetic energy of the electrons,
// Eg (eV) -> energy of synchrotron photon.

    double nu=Eg/hP;// Hz
    double p, gamma, vp, nu_c, nu_g, eps_SYN;
    nu_g=qe*B*1.0e-6/(2.0*pi*meCGS*3.0e10);// Cyclotron frequency -> Hz
    p=sqrt(pow(E+me,2)-pow(me,2));
    vp=p/(E+me);
    gamma=(E+me)/me;
    nu_c=1.5*pow(gamma,2)*nu_g/vp;// Critical frequency of synchrotron emission -> Hz
    eps_SYN=1.81*exp(-nu/nu_c)/(sqrt(pow(nu_c/nu,2.0/3.0)+pow(3.62/pi,2)));

    B*=1.0e-6;
    eps_SYN*=6.24150913e11*sqrt(3.0)*pow(qe,3)*B/(hP*Eg*meCGS*pow(3.0e10,2)); 
    
    return eps_SYN;// eV^-1 s^-1
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Photoionization cross-ection of Iron Kalpha
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double func_sigma_Kalpha_Xray(double Eg){

    double sigma;
    sigma=6.0e-18*pow(Eg/1.0e3,-2.6);

    return sigma;// cm^2
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Photoelectric absorption cross-section per H atom -> Tatischeff et al. 2012
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

double func_sigma_abs(double E){
    
    E*=1.0e-3;// keV
    
    double Eedge[]={0.03, 0.1, 0.284, 0.4, 0.532, 0.707, 0.867,1.303,1.840,2.471, 3.210, 4.038, 7.111, 8.331, 10.0};
    double c0[]={17.3, 34.6, 78.1, 71.4, 95.5, 308.9, 120.6, 141.3, 202.7, 342.7, 352.2, 433.9, 629.0, 701.2};
    double c1[]={608.1, 267.9, 18.8, 66.8, 145.8, -380.6, 169.3, 146.8, 104.7, 18.7, 18.7, -2.4, 30.9, 25.2};
    double c2[]={-2150.,-476.1,4.3,-51.4,-61.1, 294.0,-47.7, -31.5, -17.0, 0.0, 0.0, 0.75, 0.0, 0.0};
    
    int Nedge=15;
    
    double sigma_abs, C0, C1, C2;
    int j;
    
    if(E<0.03){
        sigma_abs=3.2*pow(10.0,6);
    }
    if(E>=10.0){
        C0=c0[Nedge-2];
        C1=c1[Nedge-2];
        C2=c2[Nedge-2];
        sigma_abs=C0+C1*E+C2*pow(E,2);
    }
    if((E>=0.03) && (E<10.0)){
        for(int i=0;i<Nedge-1;i++){
            if((E>=Eedge[i]) && (E<Eedge[i+1])){
                j=i;
            }
        }
        C0=c0[j];
        C1=c1[j];
        C2=c2[j];
        sigma_abs=C0+C1*E+C2*pow(E,2);
    }
    
    return sigma_abs*pow(10.0,-24)/pow(E,3);// cm^2
}