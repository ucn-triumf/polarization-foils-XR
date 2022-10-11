//積分する時に面倒なため、ビン幅の中央値にプロットするのをやめる

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TMath.h"
#include "TFile.h"

const double c=299792458; //[m/s]
const double m1=939.5654133; //[MeV/c^2]
const double m2=939.5654133e9; //[meV/c^2]
//Double_t m=939.5654133e9/(pow(c,2)); //[meV/(m/s)^2]
Double_t m=939.5654133e6/(pow(c,2)); //[eV/(m/s)^2]

const double h1=4.135667696e-15; //[eV.s]
const double h=4.135667696e-12; //[meV.s]
const double hbar=4.135667696e-12/(2*TMath::Pi());

//
const double length=1.;//[m] 距離
//const double length=2.794;//[m] 距離 検出器から発生源？
//const double length=1.33;//[m] 距離

const double str=1.e-6/length;
const double area=1.;//cm^2

//積分範囲の設定
const double Emax=5e-4;//[eV]
const double vmax=pow(2*Emax/m,0.5);
const double TOFmin=length/vmax;

const double Emax_=1e-2;//[eV]
const double vmax_=pow(2*Emax_/m,0.5);
const double TOFmin_=length/vmax_;

const double Emax_1=1.;//[eV]
const double vmax_1=pow(2*Emax_1/m,0.5);
const double TOFmin_1=length/vmax_1;

Double_t NA=6.022e23; //アボガドロ数
Double_t rho=2.648;//[g/cm^3]
Double_t rho_Fe=7.874;//[g/cm^3]
Double_t M=60.90;//[g/mol]
Double_t M_Fe=55.85;//[g/mol]
Double_t sigma_300neV=743e-24;//[cm^2]
Double_t sigma_Fe_2200m_s=2.56;//barn

Double_t sigma_Al_2200m_s=0.231;
Double_t rho_Al=2.70;//[g/cm^3]
Double_t M_Al=26.9815386;

Double_t sigma_Si_2200m_s=0.171;
Double_t rho_Si=2.3290;
Double_t M_Si=28.0855;

Double_t sigma_Mg_2200m_s=0.063;
Double_t rho_Mg=1.738;
Double_t M_Mg=24.305;

Double_t sigma_Cu_2200m_s=3.78;
Double_t rho_Cu=8.960;
Double_t M_Cu=63.546;
Double_t sigma_Mn_2200m_s=13.3;
Double_t rho_Mn=7.440;
Double_t M_Mn=54.938044;
Double_t sigma_Cr_2200m_s=3.05;
Double_t rho_Cr=7.180    ;
Double_t M_Cr=51.9961;
Double_t sigma_Zn_2200m_s=1.11;
Double_t rho_Zn=7.133;
Double_t M_Zn=65.38;
Double_t sigma_Ti_2200m_s=6.09;
Double_t rho_Ti=4.540;
Double_t M_Ti=47.867;

Double_t sigma_Cd_2200m_s=2520.;
Double_t rho_Cd=112.414;
Double_t M_Cd=8.650;

Double_t sigma_Gd_2200m_s=49700.;
Double_t rho_Gd=157.25;
Double_t M_Gd=7.901;

Double_t rho_He=51.9;
Double_t M_He=4.002602;



TF1 *fa = new TF1("fa","[0]*x*sin([1]*x)",0,6);


//追加で書き加えた部分 end//////////////////

void X_ray_reflect(){
    
    TCanvas* c = new TCanvas("c","c",800,600);
    string filepath="4_20210607.ras";
    ifstream ifs(filepath);
    double aa,bb,cc;
    
    int nnn=1500;
    
    double theta[nnn],theta_error[nnn],X_ray_R[nnn],X_ray_R_error[nnn];
    
    int i1=0;
   
    
    
    while(ifs >> aa >> bb >> cc )
    {
        theta[i1]=aa;
        X_ray_R[i1]=bb/85000;
        X_ray_R_error[i1]=cc;//?? 85000 全強度の根拠がない
    
        cout<<theta[i1]<< "_"<<X_ray_R[i1] <<endl;
        i1++;
        if(i1>nnn){
            break;
        }

    }
    
        
    
    //TGraphErrors *grR = new TGraphErrors(nnn,theta,X_ray_R,theta_error,X_ray_R_error);
    
    TGraph *grR = new TGraph(nnn,theta,X_ray_R);
    gPad->SetLogy();
    grR->Draw("AL");
    
        
    
    

  return;
}
