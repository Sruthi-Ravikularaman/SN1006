#include <iostream>
#include <complex>
#include <math.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include "gato.h"

using namespace std;

const double pi=3.14159265359;
const double me=0.510998e6;// eV
const double sigma_sb=3.5394474508e7;// erg cm^-2 s^-1 K^-4
const double mp=938.272e6;// eV
const double qe=1.602176634e-19;// SI unit

int main(){

    // Input files of the cross-section for Kalpha.    
    ifstream input1, input2;
    input1.open("sigma_6p4_p.dat");
    input2.open("sigma_6p4_e.dat");

    // Arrays to store input data.
    int N_data=13823;
    double *Xp_data=new double[N_data];
    double *Xe_data=new double[N_data];
    double *sigma_p_data=new double[N_data];
    double *sigma_e_data=new double[N_data];

    double E, sigma;

    int count=0;
    while(input1 >> E >> sigma){
        Xp_data[count]=E;// eV
        sigma_p_data[count]=sigma;// cm^2 per H (the cross section was multiplied by nFe/nH)
        
        count+=1;
    }

    count=0;
    while(input2 >> E >> sigma){
        Xe_data[count]=E;// eV
        sigma_e_data[count]=sigma;// cm^2 per H (the cross section was multiplied by nFe/nH)
        
        count+=1;
    }

    // Temperature and energy density of the radiation field including CMB, NIR, and FIR.
    double T_CMB=2.72, kappa_CMB=4.0*sigma_sb*pow(2.72,4)/3.0e10;
    cout << "CMB energy density: " << kappa_CMB << " eV/cm^3" << endl;
    cout << " " << endl;

    // Change the value for kappa to include the NIR and FIR for the target radiation field.
    double T_NIR=30.0, kappa_NIR=0.0;//0.5; 
    double T_FIR=3000.0, kappa_FIR=0.0;//1.0;

    // double u1=4.3e8, u2=u1/4.0;// u1 (cm/s) is shock speed and Rcl (cm) is cloud radius  
    // double tage=1000.0*365.0*86400.0;
    // double Vsh=filling_factor*4.0*pi*Rsh*Rsh*u2*tage;// cm^3 -> Volume of the shell 

    double ncl=0.5, Rcl=2.0*3.086e18;// cm^-3 and cm -> Cloud density and radius
    double nsh=0.03*4.0, Rsh=1.1*7.7*3.086e18;// cm^-3 and cm -> Shock density and radius 

    // Assuming only half of the shell has CRs.
    double filling_factor=0.5;

    double Vsh=filling_factor*pi*pow(Rsh,3)/3.0;// cm^3 -> Volume of the shell 
    double Vcl=4.0*pi*pow(Rcl,3)/3.0;// cm^3 -> Volume of the cloud
    double ds=2.2e3*3.086e18;// cm -> Distance from SN 1006 to Earth

    double Eg=1.0e-7, Egmax=1.0e14, dEg;// eV -> Energy of photons
    double Emin=1.0e4, Emax_p=5.0e15, Emax_e=1.0e14, dE;// eV -> Energy of CRs
    double phi_SYN, phi_ICS, phi_PPI, phi_cl, vp_p, vp_e, Db, Dk;

    // Main parameters 
    double B0=18.0;// uG -> Magnetic field in the emission region.
    double Ec_p=1.0e15, Ec_e=18.0e12;// eV -> Cut-off energies for protons and electrons
    double alpha_p=2.4, ep_p=0.05, Gam_p=func_Gam(mp,alpha_p,Emin,Emax_p);
    double alpha_e=2.4, ep_e=0.01, Gam_e=func_Gam(me,alpha_e,Emin,Emax_e); 

    // Normalization of the CR spectrum
    double Anorm_p=func_Qp(ep_p,Gam_p,alpha_p,1.0e11)*3.0e10/(4.0*pi*Vsh);// eV^-1 cm^-2 s^-1 sr^-1
    double Anorm_e=func_Qe(ep_e,Gam_e,alpha_e,1.0e11)*3.0e10/(4.0*pi*Vsh);// eV^-1 cm^-2 s^-1 sr^-1

    // Considering the SW emission region only
    Vsh*=0.5;

    // Name tags for output files
    string tag_ep_p, tag_ep_e, tag_alpha_p, tag_alpha_e, tag_Ec_e, tag_B0;
    tag_ep_p="_epp="+to_string(int(ep_p*1.0e3));
    tag_ep_e="_epe="+to_string(int(ep_e*1.0e3));
    tag_alpha_p="_alp="+to_string(int(alpha_p*100.0));
    tag_alpha_e="_ale="+to_string(int(alpha_e*100.0));
    tag_Ec_e="_Ec="+to_string(int(Ec_e*1.0e-12))+"TeV";
    tag_B0="_B0="+to_string(int(B0));

    ofstream output1, output2, output3, output4;
    output1.open("Result/gamma_SN1006"+tag_ep_p+tag_alpha_p+tag_ep_e+tag_alpha_e+tag_Ec_e+tag_B0+"_SW.dat");
    output2.open("Result/spectrum.dat");
    output3.open("Result/diff_6p4"+tag_alpha_e+"_p.dat");
    output4.open("Result/diff_6p4"+tag_alpha_e+"_e.dat");

    cout << "alpha_p = " << alpha_p << endl;
    cout << "alpha_e = " << alpha_e << endl;
    cout << " " << endl;

    cout << "A_p = " << func_Jp(Anorm_p,alpha_p,mp,Ec_p,1.0e9)/pow(1.0+2.0*1.0*mp/1.0e9,-alpha_p/2.0) << " eV^-1 cm^-2 s^-1 sr^-1"  << endl;
    cout << "A_e = " << func_Jp(Anorm_e,alpha_e,me,Ec_e,1.0e9)/pow(1.0+2.0*1.0*me/1.0e9,-alpha_e/2.0) << " eV^-1 cm^-2 s^-1 sr^-1"  << endl;
    cout << " " << endl;

    double ECR_p=func_ECR(Anorm_p,alpha_p,mp,Ec_p,Emin,Emax_p);
    double ECR_e=func_ECR(Anorm_e,alpha_e,me,Ec_e,Emin,Emax_e);
    cout << "Kinetic energy into CR protons: " << ECR_p*Vsh << " erg (" << ECR_p*Vsh*100.0/1.0e51 << "%)" << endl;
    cout << "Kinetic energy into CR electrons: " << ECR_e*Vsh << " erg (" << ECR_e*Vsh*100.0/1.0e51 << "%)" << endl;
    cout << " " << endl;

    double L=Rcl/2.0;

    cout << "Normalization of the proton spectrum at 100 GeV: " << func_Jp(Anorm_p,alpha_p,mp,Ec_p,1.0e11) << " eV^-1 cm^-2 s^-1 sr^-1" << endl;
    cout << "Normalization of the electron spectrum at 100 GeV: " << func_Jp(Anorm_e,alpha_e,me,Ec_p,1.0e11) << " eV^-1 cm^-2 s^-1 sr^-1" << endl;
    cout << " " << endl;

    string yn;
    cout << "Do gamma rays?" << endl;
    cin >> yn;

    if(yn=="y"){
        while(Eg<Egmax){
            dEg=min(0.1*Eg,Egmax-Eg);
            E=Emin;
            phi_SYN=0.0, phi_ICS=0.0, phi_PPI=0.0, phi_cl=0.0;
            while(E<Emax_p){
                dE=min(0.001*E,Emax_p-E);

                vp_p=sqrt(pow(E+mp,2)-mp*mp)*3.0e10/(E+mp);// cm/s
                vp_e=sqrt(pow(E+me,2)-me*me)*3.0e10/(E+me);// cm/s

                phi_SYN+=dE*(4.0*pi*func_Jp(Anorm_e,alpha_e,me,Ec_e,E)/vp_e)*Vsh*func_eps_SYN(B0,E,Eg);
                phi_ICS+=dE*(4.0*pi*func_Jp(Anorm_e,alpha_e,me,Ec_e,E)/vp_e)*Vsh*(func_eps_ICS(T_CMB,kappa_CMB,E,Eg)+func_eps_ICS(T_NIR,kappa_NIR,E,Eg)+func_eps_ICS(T_FIR,kappa_FIR,E,Eg));
                phi_PPI+=dE*(4.0*pi*func_Jp(Anorm_p,alpha_p,mp,Ec_p,E)/vp_p)*Vsh*func_eps_PPI(nsh,E,Eg);
                phi_cl+=dE*(4.0*pi*func_Jp(Anorm_p,alpha_p,mp,Ec_p,E)/vp_p)*Vcl*func_eps_PPI(ncl,E,Eg);

                if(Eg==1.0e-7){
                    output2 << E;
                    output2 << " " << func_Jp(Anorm_p,alpha_p,mp,Ec_p,E); 
                    output2 << " " << func_Jp(Anorm_e,alpha_p,me,Ec_e,E); 
                    output2 << " " << func_Jp_Vp(E);
                    output2 << " " << func_Jp_Ve(E);
                    output2 << " " << func_Jp_W28(E);
                    output2 << endl;
                }

                E+=dE;
            }

            phi_SYN*=1.0/(4.0*pi*ds*ds);// eV^-1 cm^-2 s^-1
            phi_ICS*=1.0/(4.0*pi*ds*ds);// eV^-1 cm^-2 s^-1
            phi_PPI*=1.0/(4.0*pi*ds*ds);// eV^-1 cm^-2 s^-1
            phi_cl*=1.0/(4.0*pi*ds*ds);// eV^-1 cm^-2 s^-1
            
            output1 << Eg << " " << Eg*Eg*phi_SYN << " " << Eg*Eg*phi_ICS << " " << Eg*Eg*phi_PPI << " " << Eg*Eg*phi_cl << endl;

            Eg+=dEg;
        }
    }

    int beg_p, beg_e;
    beg_p=0, beg_e=0;
    for(int i=0;i<N_data;i++){
        if((Xp_data[i]>=Emin) && (beg_p==0)){
            beg_p=i;
        }
    }
    for(int i=0;i<N_data;i++){
        if((Xe_data[i]>=Emin) && (beg_e==0)){
            beg_e=i;
        }
    }

    cout << "Minimum energy of cosmic rays: " << Emin << " eV" << endl;
    cout << " " << endl;

    cout << "Minimum proton energy for cross-section: " << Xp_data[beg_p] << endl;
    cout << "Minimum electron energy for cross-section: " << Xe_data[beg_e] << endl;
    cout << " " << endl;

    cout << "Distance from the enhanced emission to the shock: " << L/3.086e18 << " pc" << endl;
    cout << "Mass of the SN shell: " << nsh*Vsh*1.6726219e-27/2.0e30 << " Msol" << endl;
    cout << "Mass of the atomic cloud: " << ncl*Vcl*1.6726219e-27/2.0e30 << " Msol" << endl;
    cout << " " << endl;

    double Mcl=ncl*Vcl*1.6726219e-27/2.0e30;// Msol

    double I64_X, EX_max=8.0e3;
    I64_X=0.0;
    E=7.1e3;
    while(E<EX_max){
        dE=min(0.01*E,EX_max-E);
        I64_X+=3.0e-5*dE*ncl*func_sigma_Kalpha_Xray(E)*Vcl*pow(2.0*ds/Rcl,2)*6.0e-6*pow(E/1.0e3,-3);//*exp(-ncl*func_sigma_abs(E)*Rcl);
        //cout << E << " " << 6.0e-6*pow(E/1.0e3,-3) << endl;

        E+=dE;
    }
    I64_X*=1.0/(4.0*pi*ds*ds);

    double I64_p, I64_e, diff_p, diff_e, I64_Vp, I64_Ve, IW28;

    I64_p=0.0, I64_Vp=0.0;
    for(int i=beg_p;i<N_data-1;i++){
        dE=Xp_data[i+1]-Xp_data[i];// eV
        E=Xp_data[i];// eV
        vp_p=sqrt(pow(E+mp,2)-mp*mp)*3.0e10/(E+mp);// cm/s

        output3 << Xp_data[i];

        I64_p+=dE*ncl*sigma_p_data[i]*Vcl*4.0*pi*func_Jp(Anorm_p,alpha_p,mp,Ec_p,E);// s^-1
        if(E>=30.0e6){
            IW28+=dE*sigma_p_data[i]*4.0*pi*func_Jp_W28(E);// s^-1
        }

        output3 << " " << sigma_p_data[i];
        output3 << " " << E*sigma_p_data[i]*4.0*pi*func_Jp(Anorm_p,alpha_p,mp,Ec_p,E)*Mcl*2.0e30/(1.4*1.6726219e-27*4.0*pi*ds*ds);
        output3 << " " << E*sigma_p_data[i]*4.0*pi*func_Jp_W28(E)*5.0e4*2.0e30/(1.4*1.6726219e-27*4.0*pi*pow(2.0e3*3.086e18,2));

        I64_Vp+=dE*ncl*sigma_p_data[i]*Vcl*4.0*pi*func_Jp_Vp(E);// eV^-1 s^-1
        
        output3 << endl;
    }

    I64_e=0.0, I64_Ve=0.0;
    for(int i=beg_e;i<N_data-1;i++){
        dE=Xe_data[i+1]-Xe_data[i];// eV
        E=Xe_data[i];// eV
        vp_e=sqrt(pow(E+me,2)-me*me)*3.0e10/(E+me);// cm/s

        output4 << Xe_data[i];

        I64_e+=dE*ncl*sigma_e_data[i]*Vcl*4.0*pi*func_Jp(Anorm_e,alpha_e,me,Ec_e,E);// s^-1

        output4 << " " << sigma_e_data[i];
        output4 << " " << E*sigma_e_data[i]*4.0*pi*func_Jp(Anorm_e,alpha_e,me,Ec_e,E)*Mcl*2.0e30/(1.4*1.6726219e-27*4.0*pi*ds*ds);

        I64_Ve+=dE*ncl*sigma_e_data[i]*Vcl*4.0*pi*func_Jp_Ve(E);// eV^-1 s^-1
        
        output4 << endl;
    }

    IW28*=5.0e4*2.0e30/(1.4*1.6726219e-27*4.0*pi*pow(2.0e3*3.086e18,2));// photon cm^-2 s^-1
    I64_p*=1.0/(4.0*pi*ds*ds);// photon cm^-2 s^-1
    I64_e*=1.0/(4.0*pi*ds*ds);// photon cm^-2 s^-1
    I64_Vp*=1.0/(4.0*pi*ds*ds);// photon cm^-2 s^-1
    I64_Ve*=1.0/(4.0*pi*ds*ds);// photon cm^-2 s^-1

    cout << "Testing the SNR W28 case: " << endl;
    cout << " Observed -> " << (3.14-2.33)*1.0e-8*15.0*15.0 << endl;
    cout << " Predicted -> " << IW28 << endl;
    cout << " " << endl;

    cout << "Intensity of 6.4 keV line by X-ray: " << I64_X << " photon cm^-2 s^-1" << endl;
    cout << "Intensity of 6.4 keV line by Voyager protons: " << I64_Vp << " photon cm^-2 s^-1" << endl;
    cout << "Intensity of 6.4 keV line by Voyager electrons: " << I64_Ve << " photon cm^-2 s^-1" << endl;
    cout << " " << endl;

    cout << "Intensity of 6.4 keV line by proton: " << I64_p << " photon cm^-2 s^-1" << endl;
    cout << "Intensity of 6.4 keV line by electron: " << I64_e << " photon cm^-2 s^-1" << endl;
    cout << " " << endl;

    output1.close();
    output2.close();
    output3.close();
    input1.close();
    input1.close();

    return 0;
}
