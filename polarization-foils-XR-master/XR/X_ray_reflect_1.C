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

double M_O=15.999;

//パラメータは、薄膜密度、ラフネス
//組成
double w1_Fe2O3=2;
double w2_Fe2O3=3;//Fe2O3
//double rho1;//=4.95;
double t_1=1.4;//nm 厚さ

double w1_FeFe2O3=7;//Fe
double w2_FeFe2O3=3;//Fe2O3     Fe+Fe2O3
//double rho2;//=4.55;
double t_2=0.694;

double w1_Fe=1;//Fe
//double rho3;//=7;
double t_3=86.04;

//double w1_FeSi=1;//Fe
double w1_Si_thick=1;//Si
//double rho4;//=0.23;
double t_4=4.1;

double w1_Si=1;//Si
//double rho5;//=2.3;

//原子番号
double Z_Si=14.;
double Z_Fe=56.;
double Z_O=8.;

//X線波長
double lambda_X=1.540593e-8;//cm

//古典電子半径
double r_e=2.8e-13;//cm

//原子散乱因子の異常項　実部　虚部
double f_Si=0.244;
double f__Si=0.330;

double f_Fe=-1.179;
double f__Fe=3.204;

double f_O=0.047;
double f__O=0.032;


//sigma 粗さ
double delta_4_1;
//double rho1,rho2,rho3,rho4,rho5,d1,d2,d3,d4,d5,sigma1,sigma2,sigma3,sigma4,sigma5;

/*
double RR_4(double *two_theta1,double *par){
//,double *rho1,double *rho2,double *rho3,double *rho4,double *rho5,double *d1,double *d2,double *d3,double *d4,double *d5,double *sigma_1,double *sigma_2,double *sigma_3,double *sigma_4,double *sigma_5,

    //delta_4_1=8*r_e*(rho1*d1*rho2*d2*rho3*d3*rho4*d4*rho5*d5)/(d1+d2+d3+d4+d5)*NA*pow(lambdaTMath::Pi(),2)*(2*);
    
    double two_theta=two_theta1[0];
    double rho1=par[0];
    double rho2=par[1];
    double rho3=par[2];
    double rho4=par[3];
    double rho5=par[4];

    double d1=par[5]*1.e-7;//cm
    double d2=par[6]*1.e-7;
    double d3=par[7]*1.e-7;
    double d4=par[8]*1.e-7;
    double d5=par[9]*1.e-7;

    double sigma_1=par[10]*1.e-7;
    double sigma_2=par[11]*1.e-7;
    double sigma_3=par[12]*1.e-7;
    double sigma_4=par[13]*1.e-7;
    double sigma_5=par[14]*1.e-7;
   
    double delta_4_1=8*r_e*(rho1)*NA*pow(lambda_X*TMath::Pi(),2)*(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O))/(2/5*M_Fe+3/5*M_O);
    double delta_4_2=8*r_e*(rho2)*NA*pow(lambda_X*TMath::Pi(),2)*(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O))/(2/5*M_Fe+3/5*M_O);
    double delta_4_3=8*r_e*(rho3)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double delta_4_4=8*r_e*(rho4)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double delta_4_5=8*r_e*(rho5)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double beta_4_1=8*r_e*(rho1)*NA*pow(lambda_X*TMath::Pi(),2)*(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O))/(2/5*M_Fe+3/5*M_O);
    double beta_4_2=8*r_e*(rho2)*NA*pow(lambda_X*TMath::Pi(),2)*(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O))/(2/5*M_Fe+3/5*M_O);
    double beta_4_3=8*r_e*(rho3)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double beta_4_4=8*r_e*(rho4)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double beta_4_5=8*r_e*(rho5)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double D_delta_4_2=delta_4_1-delta_4_2;
    double D_delta_4_3=delta_4_2-delta_4_3;
    double D_delta_4_4=delta_4_3-delta_4_4;
    double D_delta_4_5=delta_4_4-delta_4_5;
    

    double D_beta_4_2=beta_4_1-beta_4_2;
    double D_beta_4_3=beta_4_2-beta_4_3;
    double D_beta_4_4=beta_4_3-beta_4_4;
    double D_beta_4_5=beta_4_4-beta_4_5;

    double qq=lambda_X/2*sin(two_theta/2);

    double tan_3_2=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_4_2=(D_delta_4_4*D_delta_4_2-D_beta_4_4*D_beta_4_2)/(D_delta_4_4*D_delta_4_2+D_beta_4_4*D_beta_4_2);
    double tan_4_3=(D_delta_4_4*D_delta_4_3-D_beta_4_4*D_beta_4_3)/(D_delta_4_4*D_delta_4_3+D_beta_4_4*D_beta_4_3);
    double tan_5_2=(D_delta_4_5*D_delta_4_2-D_beta_4_5*D_beta_4_2)/(D_delta_4_5*D_delta_4_2+D_beta_4_5*D_beta_4_2);
    double tan_5_3=(D_delta_4_5*D_delta_4_3-D_beta_4_5*D_beta_4_3)/(D_delta_4_5*D_delta_4_3+D_beta_4_5*D_beta_4_3);
    //double tan_5_3=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    //double tan_5_4=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_5_4=(D_delta_4_5*D_delta_4_4-D_beta_4_5*D_beta_4_4)/(D_delta_4_5*D_delta_4_4+D_beta_4_5*D_beta_4_4);



    double cosphi_3_2=sqrt(1/(1+pow(tan_3_2,2)));
    double cosphi_4_2=sqrt(1/(1+pow(tan_4_2,2)));
    double cosphi_4_3=sqrt(1/(1+pow(tan_4_3,2)));
    double cosphi_5_2=sqrt(1/(1+pow(tan_5_2,2)));
    double cosphi_5_3=sqrt(1/(1+pow(tan_5_3,2)));
    double cosphi_5_4=sqrt(1/(1+pow(tan_5_4,2)));

    double phi3_2=qq*(d2);
    double phi4_2=qq*(d2+d3);
    double phi4_3=qq*(d2+d3);
    double phi5_2=qq*(d2+d3+d4);
    double phi5_3=qq*(d2+d3+d4);
    double phi5_4=qq*(d2+d3+d4);

    double ccss_3_2=cosphi_3_2*(cos(phi3_2)-sin(phi3_2)*tan_3_2);
    double ccss_4_2=cosphi_4_2*(cos(phi4_2)-sin(phi4_2)*tan_4_2);
    double ccss_4_3=cosphi_4_3*(cos(phi4_3)-sin(phi4_3)*tan_4_3);
    double ccss_5_2=cosphi_5_2*(cos(phi5_2)-sin(phi5_2)*tan_5_2);
    double ccss_5_3=cosphi_5_3*(cos(phi5_3)-sin(phi5_3)*tan_5_3);
    double ccss_5_4=cosphi_5_4*(cos(phi5_4)-sin(phi5_4)*tan_5_4);

    double H_3_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_3,2)+sigma_2,2))*ccss_3_2;
    double H_4_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_2,2))*ccss_4_2;
    double H_4_3=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_3,2))*ccss_4_3;
    double H_5_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_2;
    double H_5_3=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_3;
    double H_5_4=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_4;

    double RR_Sigma_2=(pow(D_delta_4_2,2)+pow(D_beta_4_2,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_2,2));
    double RR_Sigma_3=(pow(D_delta_4_3,2)+pow(D_beta_4_3,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_3,2));
    double RR_Sigma_4=(pow(D_delta_4_4,2)+pow(D_beta_4_4,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_4,2));
    double RR_Sigma_5=(pow(D_delta_4_5,2)+pow(D_beta_4_5,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_5,2));

    //double RR1=RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma5+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;

    double RR=pow(4*TMath::Pi()/lambda_X,4)*(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    return RR;

}

*/

double RR_4(double *two_theta1,double *par){
//,double *rho1,double *rho2,double *rho3,double *rho4,double *rho5,double *d1,double *d2,double *d3,double *d4,double *d5,double *sigma_1,double *sigma_2,double *sigma_3,double *sigma_4,double *sigma_5,

    //delta_4_1=8*r_e*(rho1*d1*rho2*d2*rho3*d3*rho4*d4*rho5*d5)/(d1+d2+d3+d4+d5)*NA*pow(lambdaTMath::Pi(),2)*(2*);
    
    double two_theta=two_theta1[0];
    double rho1=par[0];
    double rho2=par[1];
    double rho3=par[2];
    double rho4=par[3];
    double rho5=par[4];

    double d1=par[5]*1.e-7;//cm
    double d2=par[6]*1.e-7;
    double d3=par[7]*1.e-7;
    double d4=par[8]*1.e-7;
    double d5=par[9]*1.e-7;

    double sigma_1=par[10]*1.e-7;
    double sigma_2=par[11]*1.e-7;
    double sigma_3=par[12]*1.e-7;
    double sigma_4=par[13]*1.e-7;
    double sigma_5=par[14]*1.e-7;

    double aaaa=par[15];

    
    double delta_4_1=8*r_e*(rho1)*NA*pow(lambda_X/4/TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_2=8*r_e*(rho2)*NA*pow(lambda_X/4/TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_3=8*r_e*(rho3)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double delta_4_4=8*r_e*(rho4)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double delta_4_5=8*r_e*(rho5)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double beta_4_1=8*r_e*(rho1)*NA*pow(lambda_X/4/TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_2=8*r_e*(rho2)*NA*pow(lambda_X/4/TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_3=8*r_e*(rho3)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double beta_4_4=8*r_e*(rho4)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double beta_4_5=8*r_e*(rho5)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double D_delta_4_2=delta_4_1-delta_4_2;
    double D_delta_4_3=delta_4_2-delta_4_3;
    double D_delta_4_4=delta_4_3-delta_4_4;
    double D_delta_4_5=delta_4_4-delta_4_5;
    
    double D_beta_4_2=beta_4_1-beta_4_2;
    double D_beta_4_3=beta_4_2-beta_4_3;
    double D_beta_4_4=beta_4_3-beta_4_4;
    double D_beta_4_5=beta_4_4-beta_4_5;

    double qq=lambda_X/2.*sin(TMath::Pi()*two_theta/2./180.);

    double tan_3_2=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_4_2=(D_delta_4_4*D_delta_4_2-D_beta_4_4*D_beta_4_2)/(D_delta_4_4*D_delta_4_2+D_beta_4_4*D_beta_4_2);
    double tan_4_3=(D_delta_4_4*D_delta_4_3-D_beta_4_4*D_beta_4_3)/(D_delta_4_4*D_delta_4_3+D_beta_4_4*D_beta_4_3);
    double tan_5_2=(D_delta_4_5*D_delta_4_2-D_beta_4_5*D_beta_4_2)/(D_delta_4_5*D_delta_4_2+D_beta_4_5*D_beta_4_2);
    double tan_5_3=(D_delta_4_5*D_delta_4_3-D_beta_4_5*D_beta_4_3)/(D_delta_4_5*D_delta_4_3+D_beta_4_5*D_beta_4_3);
    //double tan_5_3=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    //double tan_5_4=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_5_4=(D_delta_4_5*D_delta_4_4-D_beta_4_5*D_beta_4_4)/(D_delta_4_5*D_delta_4_4+D_beta_4_5*D_beta_4_4);



    double cosphi_3_2=sqrt(1/(1+pow(tan_3_2,2)));
    double cosphi_4_2=sqrt(1/(1+pow(tan_4_2,2)));
    double cosphi_4_3=sqrt(1/(1+pow(tan_4_3,2)));
    double cosphi_5_2=sqrt(1/(1+pow(tan_5_2,2)));
    double cosphi_5_3=sqrt(1/(1+pow(tan_5_3,2)));
    double cosphi_5_4=sqrt(1/(1+pow(tan_5_4,2)));

    double phi3_2=qq*(d2);
    double phi4_2=qq*(d2+d3);
    double phi4_3=qq*(d2+d3);
    double phi5_2=qq*(d2+d3+d4);
    double phi5_3=qq*(d2+d3+d4);
    double phi5_4=qq*(d2+d3+d4);

    double ccss_3_2=cosphi_3_2*(cos(phi3_2)-sin(phi3_2)*tan_3_2);
    double ccss_4_2=cosphi_4_2*(cos(phi4_2)-sin(phi4_2)*tan_4_2);
    double ccss_4_3=cosphi_4_3*(cos(phi4_3)-sin(phi4_3)*tan_4_3);
    double ccss_5_2=cosphi_5_2*(cos(phi5_2)-sin(phi5_2)*tan_5_2);
    double ccss_5_3=cosphi_5_3*(cos(phi5_3)-sin(phi5_3)*tan_5_3);
    double ccss_5_4=cosphi_5_4*(cos(phi5_4)-sin(phi5_4)*tan_5_4);



    double H_3_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_3,2)+sigma_2,2))*ccss_3_2;
    double H_4_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_2,2))*ccss_4_2;
    double H_4_3=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_3,2))*ccss_4_3;
    double H_5_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_2;
    double H_5_3=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_3;
    double H_5_4=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_4;


    double RR_Sigma_2=pow(4*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_2,2)+pow(D_beta_4_2,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_2,2));
    double RR_Sigma_3=pow(4*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_3,2)+pow(D_beta_4_3,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_3,2));
    double RR_Sigma_4=pow(4*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_4,2)+pow(D_beta_4_4,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_4,2));
    double RR_Sigma_5=pow(4*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_5,2)+pow(D_beta_4_5,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_5,2));

    //double RR1=RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma5+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;

    //double aaa=(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O));
    double aaa1=((2/5)+(3/5));
    //cout<<D_delta_4_2<<"_tan_3_2_"<<tan_3_2<<endl;//delta_4_1,2がnan
    //cout<<D_delta_4_2<<"_"<<D_delta_4_3<<endl;  

    //cout<<delta_4_1<<"_"<<delta_4_2<<endl;  

    //double RR=pow(4*TMath::Pi()/lambda_X,4)*two_theta;//(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    double RR=aaaa*(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    return RR;
    

}

double RR_4_1(double *two_theta1,double *par){
//,double *rho1,double *rho2,double *rho3,double *rho4,double *rho5,double *d1,double *d2,double *d3,double *d4,double *d5,double *sigma_1,double *sigma_2,double *sigma_3,double *sigma_4,double *sigma_5,

    //delta_4_1=8*r_e*(rho1*d1*rho2*d2*rho3*d3*rho4*d4*rho5*d5)/(d1+d2+d3+d4+d5)*NA*pow(lambdaTMath::Pi(),2)*(2*);
    
    double two_theta=two_theta1[0]*TMath::Pi()/180.;

    double rho1=par[0];
    double rho2=par[1];
    double rho3=par[2];
    double rho4=par[3];
    double rho5=par[4];

    double d1=par[5]*1.e-7;//cm
    double d2=par[6]*1.e-7;
    double d3=par[7]*1.e-7;
    double d4=par[8]*1.e-7;
    double d5=par[9]*1.e-7;

    double sigma_1=par[10]*1.e-7;
    double sigma_2=par[11]*1.e-7;
    double sigma_3=par[12]*1.e-7;
    double sigma_4=par[13]*1.e-7;
    double sigma_5=par[14]*1.e-7;

    double aaaa=par[15];

    /*
    double rho1=5.;
    double rho2=5.1;
    double rho3=7.;
    double rho4=0.23;
    double rho5=2.3;

    double d1=1.44e-7;//cm
    double d2=0.694e-7;
    double d3=86.e-7;
    double d4=4.1e-7;
    double d5=1000000e-7;

    */

/*
    double d1=1.44e-7;//cm
    double d2=0.694e-7;
    double d3=86.e-7;
    double d4=4.1e-7;
    double d5=1000000e-7;


    double sigma_1=0.5e-7;
    double sigma_2=0.;
    double sigma_3=0.3e-7;
    double sigma_4=2.8e-7;
    double sigma_5=0.3e-7;

    
*/
    
    double delta_4_1=8.*r_e*(rho1)*NA*pow(lambda_X/4./TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_2=8.*r_e*(rho2)*NA*pow(lambda_X/4./TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_3=8.*r_e*(rho3)*NA*pow(lambda_X/4./TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double delta_4_4=8.*r_e*(rho4)*NA*pow(lambda_X/4./TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double delta_4_5=8.*r_e*(rho5)*NA*pow(lambda_X/4./TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double beta_4_1=8.*r_e*(rho1)*NA*pow(lambda_X/4./TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_2=8.*r_e*(rho2)*NA*pow(lambda_X/4./TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_3=8.*r_e*(rho3)*NA*pow(lambda_X/4./TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double beta_4_4=8.*r_e*(rho4)*NA*pow(lambda_X/4./TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double beta_4_5=8.*r_e*(rho5)*NA*pow(lambda_X/4./TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double D_delta_4_2=delta_4_1-delta_4_2;
    double D_delta_4_3=delta_4_2-delta_4_3;
    double D_delta_4_4=delta_4_3-delta_4_4;
    double D_delta_4_5=delta_4_4-delta_4_5;
    

    double D_beta_4_2=beta_4_1-beta_4_2;
    double D_beta_4_3=beta_4_2-beta_4_3;
    double D_beta_4_4=beta_4_3-beta_4_4;
    double D_beta_4_5=beta_4_4-beta_4_5;
/*
    complex<double> n_4_1(1-delta_4_1,-beta_4_1);
    complex<double> n_4_2(1-delta_4_2,-beta_4_2);
    complex<double> n_4_3(1-delta_4_3,-beta_4_3);
    complex<double> n_4_4(1-delta_4_4,-beta_4_4);
    complex<double> n_4_5(1-delta_4_5,-beta_4_5);
    //complex<double> n_4_1(1-delta_4_1,-beta_4_1);

    complex<double>g_4_1=sqrt(pow(n_4_1,2)-pow(cos(two_theta/2.),2));
    complex<double>g_4_2=sqrt(pow(n_4_2,2)-pow(cos(two_theta/2.),2));
    complex<double>g_4_3=sqrt(pow(n_4_3,2)-pow(cos(two_theta/2.),2));
    complex<double>g_4_4=sqrt(pow(n_4_4,2)-pow(cos(two_theta/2.),2));
    complex<double>g_4_5=sqrt(pow(n_4_5,2)-pow(cos(two_theta/2.),2));

    complex<double>gamma_4_1=pow(4*TMath::Pi()*g_4_1,2);

    //complex<double> 
*/
    //double qq=lambda_X/2*sin(two_theta/2);
    //double qq=4.*TMath::Pi()/lambda_X*(pow()-pow(cos(two_theta/2.),2));//sin(two_theta/2.);
    double qq=4.*TMath::Pi()/lambda_X*sin(two_theta/2.);

    double tan_3_2=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_4_2=(D_delta_4_4*D_delta_4_2-D_beta_4_4*D_beta_4_2)/(D_delta_4_4*D_delta_4_2+D_beta_4_4*D_beta_4_2);
    double tan_4_3=(D_delta_4_4*D_delta_4_3-D_beta_4_4*D_beta_4_3)/(D_delta_4_4*D_delta_4_3+D_beta_4_4*D_beta_4_3);
    double tan_5_2=(D_delta_4_5*D_delta_4_2-D_beta_4_5*D_beta_4_2)/(D_delta_4_5*D_delta_4_2+D_beta_4_5*D_beta_4_2);
    double tan_5_3=(D_delta_4_5*D_delta_4_3-D_beta_4_5*D_beta_4_3)/(D_delta_4_5*D_delta_4_3+D_beta_4_5*D_beta_4_3);
    //double tan_5_3=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    //double tan_5_4=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_5_4=(D_delta_4_5*D_delta_4_4-D_beta_4_5*D_beta_4_4)/(D_delta_4_5*D_delta_4_4+D_beta_4_5*D_beta_4_4);



    double cosphi_3_2=sqrt(1./(1.+pow(tan_3_2,2)));
    double cosphi_4_2=sqrt(1./(1.+pow(tan_4_2,2)));
    double cosphi_4_3=sqrt(1./(1.+pow(tan_4_3,2)));
    double cosphi_5_2=sqrt(1./(1.+pow(tan_5_2,2)));
    double cosphi_5_3=sqrt(1./(1.+pow(tan_5_3,2)));
    double cosphi_5_4=sqrt(1./(1.+pow(tan_5_4,2)));

    double phi3_2=qq*(d2);
    double phi4_2=qq*(d2+d3);
    double phi4_3=qq*(d2+d3);
    double phi5_2=qq*(d2+d3+d4);
    double phi5_3=qq*(d2+d3+d4);
    double phi5_4=qq*(d2+d3+d4);

    double ccss_3_2=abs(cosphi_3_2*(cos(phi3_2)-sin(phi3_2)*tan_3_2));
    double ccss_4_2=abs(cosphi_4_2*(cos(phi4_2)-sin(phi4_2)*tan_4_2));
    double ccss_4_3=abs(cosphi_4_3*(cos(phi4_3)-sin(phi4_3)*tan_4_3));
    double ccss_5_2=abs(cosphi_5_2*(cos(phi5_2)-sin(phi5_2)*tan_5_2));
    double ccss_5_3=abs(cosphi_5_3*(cos(phi5_3)-sin(phi5_3)*tan_5_3));
    double ccss_5_4=abs(cosphi_5_4*(cos(phi5_4)-sin(phi5_4)*tan_5_4));



    double H_3_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1./2.)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_3,2)+pow(sigma_2,2)))*ccss_3_2;
    double H_4_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1./2.)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_4,2)+pow(sigma_2,2)))*ccss_4_2;
    double H_4_3=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1./2.)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_4,2)+pow(sigma_3,2)))*ccss_4_3;
    double H_5_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1./2.)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_5,2)+pow(sigma_2,2)))*ccss_5_2;
    double H_5_3=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1./2.)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_5,2)+pow(sigma_2,2)))*ccss_5_3;
    double H_5_4=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1./2.)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_5,2)+pow(sigma_2,2)))*ccss_5_4;


    double RR_Sigma_2=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_2,2)+pow(D_beta_4_2,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_2,2));
    double RR_Sigma_3=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_3,2)+pow(D_beta_4_3,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_3,2));
    double RR_Sigma_4=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_4,2)+pow(D_beta_4_4,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_4,2));
    double RR_Sigma_5=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_5,2)+pow(D_beta_4_5,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_5,2));

    //double RR1=RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma5+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;

    //double aaa=(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O));
    double aaa1=((2./5.)+(3./5.));
    //cout<<D_delta_4_2<<"_tan_3_2_"<<tan_3_2<<endl;//delta_4_1,2がnan
    //cout<<D_delta_4_2<<"_"<<D_delta_4_3<<endl;  

    //cout<<delta_4_1<<"_"<<delta_4_2<<endl;  

    //complex<double> II(0.,1.);
    //double RR=pow(4*TMath::Pi()/lambda_X,4)*two_theta;//(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    //complex<double>RR1=pow(II,2);

    double RR=aaaa*((RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4);
    //double RR=real(RR1);
    return RR;

}

///double two_theta1;

/*20210623 近似式であることに気づき訂正
double RR_4a(double two_theta1){
//,double *rho1,double *rho2,double *rho3,double *rho4,double *rho5,double *d1,double *d2,double *d3,double *d4,double *d5,double *sigma_1,double *sigma_2,double *sigma_3,double *sigma_4,double *sigma_5,

    //delta_4_1=8*r_e*(rho1*d1*rho2*d2*rho3*d3*rho4*d4*rho5*d5)/(d1+d2+d3+d4+d5)*NA*pow(lambdaTMath::Pi(),2)*(2*);
    
    double two_theta=two_theta1*TMath::Pi()/180.;
    double rho1=5.;
    double rho2=5.1;
    double rho3=7.;
    double rho4=0.23;
    double rho5=2.3;

    double d1=1.44e-7;//cm
    double d2=0.694e-7;
    double d3=86.e-7;
    double d4=4.1e-7;
    double d5=1000000e-7;

    double sigma_1=0.5e-7;
    double sigma_2=0.;
    double sigma_3=0.3e-7;
    double sigma_4=2.8e-7;
    double sigma_5=0.3e-7;

    
    double delta_4_1=8.*r_e*(rho1)*NA*pow(lambda_X*TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_2=8.*r_e*(rho2)*NA*pow(lambda_X*TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_3=8.*r_e*(rho3)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double delta_4_4=8.*r_e*(rho4)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double delta_4_5=8.*r_e*(rho5)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double beta_4_1=8.*r_e*(rho1)*NA*pow(lambda_X*TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_2=8.*r_e*(rho2)*NA*pow(lambda_X*TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_3=8.*r_e*(rho3)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double beta_4_4=8.*r_e*(rho4)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double beta_4_5=8.*r_e*(rho5)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double D_delta_4_2=delta_4_1-delta_4_2;
    double D_delta_4_3=delta_4_2-delta_4_3;
    double D_delta_4_4=delta_4_3-delta_4_4;
    double D_delta_4_5=delta_4_4-delta_4_5;
    

    double D_beta_4_2=beta_4_1-beta_4_2;
    double D_beta_4_3=beta_4_2-beta_4_3;
    double D_beta_4_4=beta_4_3-beta_4_4;
    double D_beta_4_5=beta_4_4-beta_4_5;

    //double qq=lambda_X/2*sin(two_theta/2);
    double qq=4*TMath::Pi()/lambda_X*sin(two_theta/2);

    double tan_3_2=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_4_2=(D_delta_4_4*D_delta_4_2-D_beta_4_4*D_beta_4_2)/(D_delta_4_4*D_delta_4_2+D_beta_4_4*D_beta_4_2);
    double tan_4_3=(D_delta_4_4*D_delta_4_3-D_beta_4_4*D_beta_4_3)/(D_delta_4_4*D_delta_4_3+D_beta_4_4*D_beta_4_3);
    double tan_5_2=(D_delta_4_5*D_delta_4_2-D_beta_4_5*D_beta_4_2)/(D_delta_4_5*D_delta_4_2+D_beta_4_5*D_beta_4_2);
    double tan_5_3=(D_delta_4_5*D_delta_4_3-D_beta_4_5*D_beta_4_3)/(D_delta_4_5*D_delta_4_3+D_beta_4_5*D_beta_4_3);
    //double tan_5_3=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    //double tan_5_4=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_5_4=(D_delta_4_5*D_delta_4_4-D_beta_4_5*D_beta_4_4)/(D_delta_4_5*D_delta_4_4+D_beta_4_5*D_beta_4_4);



    double cosphi_3_2=sqrt(1./(1+pow(tan_3_2,2)));
    double cosphi_4_2=sqrt(1./(1+pow(tan_4_2,2)));
    double cosphi_4_3=sqrt(1./(1+pow(tan_4_3,2)));
    double cosphi_5_2=sqrt(1./(1+pow(tan_5_2,2)));
    double cosphi_5_3=sqrt(1./(1+pow(tan_5_3,2)));
    double cosphi_5_4=sqrt(1./(1+pow(tan_5_4,2)));

    double phi3_2=qq*(d2);
    double phi4_2=qq*(d2+d3);
    double phi4_3=qq*(d2+d3);
    double phi5_2=qq*(d2+d3+d4);
    double phi5_3=qq*(d2+d3+d4);
    double phi5_4=qq*(d2+d3+d4);

    double ccss_3_2=cosphi_3_2*(cos(phi3_2)-sin(phi3_2)*tan_3_2);
    double ccss_4_2=cosphi_4_2*(cos(phi4_2)-sin(phi4_2)*tan_4_2);
    double ccss_4_3=cosphi_4_3*(cos(phi4_3)-sin(phi4_3)*tan_4_3);
    double ccss_5_2=cosphi_5_2*(cos(phi5_2)-sin(phi5_2)*tan_5_2);
    double ccss_5_3=cosphi_5_3*(cos(phi5_3)-sin(phi5_3)*tan_5_3);
    double ccss_5_4=cosphi_5_4*(cos(phi5_4)-sin(phi5_4)*tan_5_4);



    double H_3_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_3,2)+sigma_2,2))*ccss_3_2;
    double H_4_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_2,2))*ccss_4_2;
    double H_4_3=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_3,2))*ccss_4_3;
    double H_5_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_2;
    double H_5_3=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_3;
    double H_5_4=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_4;


    double RR_Sigma_2=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_2,2)+pow(D_beta_4_2,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_2,2));
    double RR_Sigma_3=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_3,2)+pow(D_beta_4_3,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_3,2));
    double RR_Sigma_4=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_4,2)+pow(D_beta_4_4,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_4,2));
    double RR_Sigma_5=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_5,2)+pow(D_beta_4_5,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_5,2));

    //double RR1=RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma5+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;

    //double aaa=(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O));
    double aaa1=((2/5)+(3/5));
    //cout<<D_delta_4_2<<"_tan_3_2_"<<tan_3_2<<endl;//delta_4_1,2がnan
    //cout<<D_delta_4_2<<"_"<<D_delta_4_3<<endl;  

    //cout<<delta_4_1<<"_"<<delta_4_2<<endl;  

    //double RR=pow(4*TMath::Pi()/lambda_X,4)*two_theta;//(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    double RR=(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    return RR;

}
*/

double RR_4a(double two_theta1){
//,double *rho1,double *rho2,double *rho3,double *rho4,double *rho5,double *d1,double *d2,double *d3,double *d4,double *d5,double *sigma_1,double *sigma_2,double *sigma_3,double *sigma_4,double *sigma_5,

    //delta_4_1=8*r_e*(rho1*d1*rho2*d2*rho3*d3*rho4*d4*rho5*d5)/(d1+d2+d3+d4+d5)*NA*pow(lambdaTMath::Pi(),2)*(2*);
    
    double two_theta=two_theta1*TMath::Pi()/180.;
    double rho1=5.;
    double rho2=5.1;
    double rho3=7.;
    double rho4=0.23;
    double rho5=2.3;

    double d1=1.44e-7;//cm
    double d2=0.694e-7;
    double d3=86.e-7;
    double d4=4.1e-7;
    double d5=1000000e-7;
/*
    double d1=1.44e-7;//cm
    double d2=0.694e-7;
    double d3=86.e-7;
    double d4=4.1e-7;
    double d5=1000000e-7;
*/

    double sigma_1=0.5e-7;
    double sigma_2=0.;
    double sigma_3=0.3e-7;
    double sigma_4=2.8e-7;
    double sigma_5=0.3e-7;

    
    double delta_4_1=8.*r_e*(rho1)*NA*pow(lambda_X/4/TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_2=8.*r_e*(rho2)*NA*pow(lambda_X/4/TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_3=8.*r_e*(rho3)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double delta_4_4=8.*r_e*(rho4)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double delta_4_5=8.*r_e*(rho5)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double beta_4_1=8.*r_e*(rho1)*NA*pow(lambda_X/4/TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_2=8.*r_e*(rho2)*NA*pow(lambda_X/4/TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_3=8.*r_e*(rho3)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double beta_4_4=8.*r_e*(rho4)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double beta_4_5=8.*r_e*(rho5)*NA*pow(lambda_X/4/TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double D_delta_4_2=delta_4_1-delta_4_2;
    double D_delta_4_3=delta_4_2-delta_4_3;
    double D_delta_4_4=delta_4_3-delta_4_4;
    double D_delta_4_5=delta_4_4-delta_4_5;
    

    double D_beta_4_2=beta_4_1-beta_4_2;
    double D_beta_4_3=beta_4_2-beta_4_3;
    double D_beta_4_4=beta_4_3-beta_4_4;
    double D_beta_4_5=beta_4_4-beta_4_5;
/*
    complex<double> n_4_1(1-delta_4_1,-beta_4_1);
    complex<double> n_4_2(1-delta_4_2,-beta_4_2);
    complex<double> n_4_3(1-delta_4_3,-beta_4_3);
    complex<double> n_4_4(1-delta_4_4,-beta_4_4);
    complex<double> n_4_5(1-delta_4_5,-beta_4_5);
    //complex<double> n_4_1(1-delta_4_1,-beta_4_1);

    complex<double>g_4_1=sqrt(pow(n_4_1,2)-pow(cos(two_theta/2.),2));
    complex<double>g_4_2=sqrt(pow(n_4_2,2)-pow(cos(two_theta/2.),2));
    complex<double>g_4_3=sqrt(pow(n_4_3,2)-pow(cos(two_theta/2.),2));
    complex<double>g_4_4=sqrt(pow(n_4_4,2)-pow(cos(two_theta/2.),2));
    complex<double>g_4_5=sqrt(pow(n_4_5,2)-pow(cos(two_theta/2.),2));

    complex<double>gamma_4_1=pow(4*TMath::Pi()*g_4_1,2);

    //complex<double> 
*/
    //double qq=lambda_X/2*sin(two_theta/2);
    double qq=4*TMath::Pi()/lambda_X*sin(two_theta/2.);

    double tan_3_2=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_4_2=(D_delta_4_4*D_delta_4_2-D_beta_4_4*D_beta_4_2)/(D_delta_4_4*D_delta_4_2+D_beta_4_4*D_beta_4_2);
    double tan_4_3=(D_delta_4_4*D_delta_4_3-D_beta_4_4*D_beta_4_3)/(D_delta_4_4*D_delta_4_3+D_beta_4_4*D_beta_4_3);
    double tan_5_2=(D_delta_4_5*D_delta_4_2-D_beta_4_5*D_beta_4_2)/(D_delta_4_5*D_delta_4_2+D_beta_4_5*D_beta_4_2);
    double tan_5_3=(D_delta_4_5*D_delta_4_3-D_beta_4_5*D_beta_4_3)/(D_delta_4_5*D_delta_4_3+D_beta_4_5*D_beta_4_3);
    //double tan_5_3=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    //double tan_5_4=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_5_4=(D_delta_4_5*D_delta_4_4-D_beta_4_5*D_beta_4_4)/(D_delta_4_5*D_delta_4_4+D_beta_4_5*D_beta_4_4);



    double cosphi_3_2=sqrt(1./(1+pow(tan_3_2,2)));
    double cosphi_4_2=sqrt(1./(1+pow(tan_4_2,2)));
    double cosphi_4_3=sqrt(1./(1+pow(tan_4_3,2)));
    double cosphi_5_2=sqrt(1./(1+pow(tan_5_2,2)));
    double cosphi_5_3=sqrt(1./(1+pow(tan_5_3,2)));
    double cosphi_5_4=sqrt(1./(1+pow(tan_5_4,2)));

    double phi3_2=qq*(d2);
    double phi4_2=qq*(d2+d3);
    double phi4_3=qq*(d2+d3);
    double phi5_2=qq*(d2+d3+d4);
    double phi5_3=qq*(d2+d3+d4);
    double phi5_4=qq*(d2+d3+d4);

    double ccss_3_2=abs(cosphi_3_2*(cos(phi3_2)-sin(phi3_2)*tan_3_2));
    double ccss_4_2=abs(cosphi_4_2*(cos(phi4_2)-sin(phi4_2)*tan_4_2));
    double ccss_4_3=abs(cosphi_4_3*(cos(phi4_3)-sin(phi4_3)*tan_4_3));
    double ccss_5_2=abs(cosphi_5_2*(cos(phi5_2)-sin(phi5_2)*tan_5_2));
    double ccss_5_3=abs(cosphi_5_3*(cos(phi5_3)-sin(phi5_3)*tan_5_3));
    double ccss_5_4=abs(cosphi_5_4*(cos(phi5_4)-sin(phi5_4)*tan_5_4));



    double H_3_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1./2.)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_3,2)+pow(sigma_2,2)))*ccss_3_2;
    double H_4_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1./2.)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_4,2)+pow(sigma_2,2)))*ccss_4_2;
    double H_4_3=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1./2.)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_4,2)+pow(sigma_3,2)))*ccss_4_3;
    double H_5_2=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1./2.)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_5,2)+pow(sigma_2,2)))*ccss_5_2;
    double H_5_3=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1./2.)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_5,2)+pow(sigma_2,2)))*ccss_5_3;
    double H_5_4=1./2.*pow(qq,-4)*pow(4.*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1./2.)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1./2.)*exp(-pow(qq/2.,2)*(pow(sigma_5,2)+pow(sigma_2,2)))*ccss_5_4;


    double RR_Sigma_2=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_2,2)+pow(D_beta_4_2,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_2,2));
    double RR_Sigma_3=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_3,2)+pow(D_beta_4_3,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_3,2));
    double RR_Sigma_4=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_4,2)+pow(D_beta_4_4,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_4,2));
    double RR_Sigma_5=pow(4.*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_5,2)+pow(D_beta_4_5,2))/(4.*pow(qq,4))*exp(-pow(qq*sigma_5,2));

    //double RR1=RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma5+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;

    //double aaa=(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O));
    double aaa1=((2/5)+(3/5));
    //cout<<D_delta_4_2<<"_tan_3_2_"<<tan_3_2<<endl;//delta_4_1,2がnan
    //cout<<D_delta_4_2<<"_"<<D_delta_4_3<<endl;  

    //cout<<delta_4_1<<"_"<<delta_4_2<<endl;  

    //complex<double> II(0.,1.);
    //double RR=pow(4*TMath::Pi()/lambda_X,4)*two_theta;//(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    //complex<double>RR1=pow(II,2);

    double RR=(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    //double RR=real(RR1);
    return RR;

}

/* Tgraphで描いてみた
double RR_4a(double two_theta1){
//,double *rho1,double *rho2,double *rho3,double *rho4,double *rho5,double *d1,double *d2,double *d3,double *d4,double *d5,double *sigma_1,double *sigma_2,double *sigma_3,double *sigma_4,double *sigma_5,

    //delta_4_1=8*r_e*(rho1*d1*rho2*d2*rho3*d3*rho4*d4*rho5*d5)/(d1+d2+d3+d4+d5)*NA*pow(lambdaTMath::Pi(),2)*(2*);
    
    double two_theta=two_theta1;
    double rho1=5.;
    double rho2=5.1;
    double rho3=7.;
    double rho4=0.23;
    double rho5=2.3;

    double d1=1.44e-7;//cm
    double d2=0.694e-7;
    double d3=86.e-7;
    double d4=4.1e-7;
    double d5=1000000e-7;

    double sigma_1=0.5e-7;
    double sigma_2=0.;
    double sigma_3=0.3e-7;
    double sigma_4=2.8e-7;
    double sigma_5=0.3e-7;

    
    double delta_4_1=8*r_e*(rho1)*NA*pow(lambda_X*TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_2=8*r_e*(rho2)*NA*pow(lambda_X*TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double delta_4_3=8*r_e*(rho3)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double delta_4_4=8*r_e*(rho4)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double delta_4_5=8*r_e*(rho5)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double beta_4_1=8*r_e*(rho1)*NA*pow(lambda_X*TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_2=8*r_e*(rho2)*NA*pow(lambda_X*TMath::Pi(),2)*(2./5.*(Z_Fe+f_Fe)+3./5.*(Z_O+f_O))/(2./5.*M_Fe+3./5.*M_O);
    double beta_4_3=8*r_e*(rho3)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double beta_4_4=8*r_e*(rho4)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double beta_4_5=8*r_e*(rho5)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double D_delta_4_2=delta_4_1-delta_4_2;
    double D_delta_4_3=delta_4_2-delta_4_3;
    double D_delta_4_4=delta_4_3-delta_4_4;
    double D_delta_4_5=delta_4_4-delta_4_5;
    

    double D_beta_4_2=beta_4_1-beta_4_2;
    double D_beta_4_3=beta_4_2-beta_4_3;
    double D_beta_4_4=beta_4_3-beta_4_4;
    double D_beta_4_5=beta_4_4-beta_4_5;

    double qq=lambda_X/2*sin(two_theta/2);

    double tan_3_2=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_4_2=(D_delta_4_4*D_delta_4_2-D_beta_4_4*D_beta_4_2)/(D_delta_4_4*D_delta_4_2+D_beta_4_4*D_beta_4_2);
    double tan_4_3=(D_delta_4_4*D_delta_4_3-D_beta_4_4*D_beta_4_3)/(D_delta_4_4*D_delta_4_3+D_beta_4_4*D_beta_4_3);
    double tan_5_2=(D_delta_4_5*D_delta_4_2-D_beta_4_5*D_beta_4_2)/(D_delta_4_5*D_delta_4_2+D_beta_4_5*D_beta_4_2);
    double tan_5_3=(D_delta_4_5*D_delta_4_3-D_beta_4_5*D_beta_4_3)/(D_delta_4_5*D_delta_4_3+D_beta_4_5*D_beta_4_3);
    //double tan_5_3=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    //double tan_5_4=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_5_4=(D_delta_4_5*D_delta_4_4-D_beta_4_5*D_beta_4_4)/(D_delta_4_5*D_delta_4_4+D_beta_4_5*D_beta_4_4);



    double cosphi_3_2=sqrt(1/(1+pow(tan_3_2,2)));
    double cosphi_4_2=sqrt(1/(1+pow(tan_4_2,2)));
    double cosphi_4_3=sqrt(1/(1+pow(tan_4_3,2)));
    double cosphi_5_2=sqrt(1/(1+pow(tan_5_2,2)));
    double cosphi_5_3=sqrt(1/(1+pow(tan_5_3,2)));
    double cosphi_5_4=sqrt(1/(1+pow(tan_5_4,2)));

    double phi3_2=qq*(d2);
    double phi4_2=qq*(d2+d3);
    double phi4_3=qq*(d2+d3);
    double phi5_2=qq*(d2+d3+d4);
    double phi5_3=qq*(d2+d3+d4);
    double phi5_4=qq*(d2+d3+d4);

    double ccss_3_2=cosphi_3_2*(cos(phi3_2)-sin(phi3_2)*tan_3_2);
    double ccss_4_2=cosphi_4_2*(cos(phi4_2)-sin(phi4_2)*tan_4_2);
    double ccss_4_3=cosphi_4_3*(cos(phi4_3)-sin(phi4_3)*tan_4_3);
    double ccss_5_2=cosphi_5_2*(cos(phi5_2)-sin(phi5_2)*tan_5_2);
    double ccss_5_3=cosphi_5_3*(cos(phi5_3)-sin(phi5_3)*tan_5_3);
    double ccss_5_4=cosphi_5_4*(cos(phi5_4)-sin(phi5_4)*tan_5_4);



    double H_3_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_3,2)+sigma_2,2))*ccss_3_2;
    double H_4_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_2,2))*ccss_4_2;
    double H_4_3=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_3,2))*ccss_4_3;
    double H_5_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_2;
    double H_5_3=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_3;
    double H_5_4=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_4;


    double RR_Sigma_2=pow(4*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_2,2)+pow(D_beta_4_2,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_2,2));
    double RR_Sigma_3=pow(4*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_3,2)+pow(D_beta_4_3,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_3,2));
    double RR_Sigma_4=pow(4*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_4,2)+pow(D_beta_4_4,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_4,2));
    double RR_Sigma_5=pow(4*TMath::Pi()/lambda_X,4)*(pow(D_delta_4_5,2)+pow(D_beta_4_5,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_5,2));

    //double RR1=RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma5+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;

    //double aaa=(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O));
    double aaa1=((2/5)+(3/5));
    //cout<<D_delta_4_2<<"_tan_3_2_"<<tan_3_2<<endl;//delta_4_1,2がnan
    //cout<<D_delta_4_2<<"_"<<D_delta_4_3<<endl;  

    //cout<<delta_4_1<<"_"<<delta_4_2<<endl;  

    //double RR=pow(4*TMath::Pi()/lambda_X,4)*two_theta;//(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    double RR=(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    return RR;

}


*/


/*
double RR_4(double *two_theta,double *rho1,double *rho2,double *rho3,double *rho4,double *rho5,double *d1,double *d2,double *d3,double *d4,double *d5,double *sigma_1,double *sigma_2,double *sigma_3,double *sigma_4,double *sigma_5,double *par){
    //delta_4_1=8*r_e*(rho1*d1*rho2*d2*rho3*d3*rho4*d4*rho5*d5)/(d1+d2+d3+d4+d5)*NA*pow(lambdaTMath::Pi(),2)*(2*);
    
    

    
    double delta_4_1=8*r_e*(rho1)*NA*pow(lambda_X*TMath::Pi(),2)*(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O))/(2/5*M_Fe+3/5*M_O);
    double delta_4_2=8*r_e*(rho2)*NA*pow(lambda_X*TMath::Pi(),2)*(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O))/(2/5*M_Fe+3/5*M_O);
    double delta_4_3=8*r_e*(rho3)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double delta_4_4=8*r_e*(rho4)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double delta_4_5=8*r_e*(rho5)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double beta_4_1=8*r_e*(rho1)*NA*pow(lambda_X*TMath::Pi(),2)*(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O))/(2/5*M_Fe+3/5*M_O);
    double beta_4_2=8*r_e*(rho2)*NA*pow(lambda_X*TMath::Pi(),2)*(2/5*(Z_Fe+f_Fe)+3/5*(Z_O+f_O))/(2/5*M_Fe+3/5*M_O);
    double beta_4_3=8*r_e*(rho3)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    //delta_4_4=8*r_e*(rho4)*NA*pow(lambdaTMath::Pi(),2)*((Z_Fe+f_Fe))/(M_Fe);
    double beta_4_4=8*r_e*(rho4)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);
    double beta_4_5=8*r_e*(rho5)*NA*pow(lambda_X*TMath::Pi(),2)*((Z_Si+f_Si))/(M_Si);

    double D_delta_4_2=delta_4_1-delta_4_2;
    double D_delta_4_3=delta_4_2-delta_4_3;
    double D_delta_4_4=delta_4_3-delta_4_4;
    double D_delta_4_5=delta_4_4-delta_4_5;
    

    double D_beta_4_2=beta_4_1-beta_4_2;
    double D_beta_4_3=beta_4_2-beta_4_3;
    double D_beta_4_4=beta_4_3-beta_4_4;
    double D_beta_4_5=beta_4_4-beta_4_5;

    double qq=lambda_X/2*sin(two_theta/2);

    double tan_3_2=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_4_2=(D_delta_4_4*D_delta_4_2-D_beta_4_4*D_beta_4_2)/(D_delta_4_4*D_delta_4_2+D_beta_4_4*D_beta_4_2);
    double tan_4_3=(D_delta_4_4*D_delta_4_3-D_beta_4_4*D_beta_4_3)/(D_delta_4_4*D_delta_4_3+D_beta_4_4*D_beta_4_3);
    double tan_5_2=(D_delta_4_5*D_delta_4_2-D_beta_4_5*D_beta_4_2)/(D_delta_4_5*D_delta_4_2+D_beta_4_5*D_beta_4_2);
    double tan_5_3=(D_delta_4_5*D_delta_4_3-D_beta_4_5*D_beta_4_3)/(D_delta_4_5*D_delta_4_3+D_beta_4_5*D_beta_4_3);
    //double tan_5_3=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    //double tan_5_4=(D_delta_4_3*D_delta_4_2-D_beta_4_3*D_beta_4_2)/(D_delta_4_3*D_delta_4_2+D_beta_4_3*D_beta_4_2);
    double tan_5_4=(D_delta_4_5*D_delta_4_4-D_beta_4_5*D_beta_4_4)/(D_delta_4_5*D_delta_4_4+D_beta_4_5*D_beta_4_4);



    double cosphi_3_2=sqrt(1/(1+pow(tan_3_2,2)));
    double cosphi_4_2=sqrt(1/(1+pow(tan_4_2,2)));
    double cosphi_4_3=sqrt(1/(1+pow(tan_4_3,2)));
    double cosphi_5_2=sqrt(1/(1+pow(tan_5_2,2)));
    double cosphi_5_3=sqrt(1/(1+pow(tan_5_3,2)));
    double cosphi_5_4=sqrt(1/(1+pow(tan_5_4,2)));

    double phi3_2=qq*(d2);
    double phi4_2=qq*(d2+d3);
    double phi4_3=qq*(d2+d3);
    double phi5_2=qq*(d2+d3+d4);
    double phi5_3=qq*(d2+d3+d4);
    double phi5_4=qq*(d2+d3+d4);

    double ccss_3_2=cosphi_3_2*(cos(phi3_2)-sin(phi3_2)*tan_3_2);
    double ccss_4_2=cosphi_4_2*(cos(phi4_2)-sin(phi4_2)*tan_4_2);
    double ccss_4_3=cosphi_4_3*(cos(phi4_3)-sin(phi4_3)*tan_4_3);
    double ccss_5_2=cosphi_5_2*(cos(phi5_2)-sin(phi5_2)*tan_5_2);
    double ccss_5_3=cosphi_5_3*(cos(phi5_3)-sin(phi5_3)*tan_5_3);
    double ccss_5_4=cosphi_5_4*(cos(phi5_4)-sin(phi5_4)*tan_5_4);



    double H_3_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_3,2)+sigma_2,2))*ccss_3_2;
    double H_4_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_2,2))*ccss_4_2;
    double H_4_3=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_4,2)+sigma_3,2))*ccss_4_3;
    double H_5_2=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_2,2)+pow(D_beta_4_2,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_2;
    double H_5_3=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_3,2)+pow(D_beta_4_3,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_3;
    double H_5_4=1/2*pow(qq,4)*pow(4*TMath::Pi()/lambda_X,4)*pow((pow(D_delta_4_5,2)+pow(D_beta_4_5,2)),1/2)*pow((pow(D_delta_4_4,2)+pow(D_beta_4_4,2)),1/2)*exp(-pow(qq,2)*(pow(sigma_5,2)+sigma_2,2))*ccss_5_4;


    double RR_Sigma_2=(pow(D_delta_4_2,2)+pow(D_beta_4_2,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_2,2));
    double RR_Sigma_3=(pow(D_delta_4_3,2)+pow(D_beta_4_3,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_3,2));
    double RR_Sigma_4=(pow(D_delta_4_4,2)+pow(D_beta_4_4,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_4,2));
    double RR_Sigma_5=(pow(D_delta_4_5,2)+pow(D_beta_4_5,2))/(4*pow(qq,4))*exp(-pow(qq*sigma_5,2));

    //double RR1=RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma5+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;

    double RR=pow(4*TMath::Pi()/lambda_X,4)*(RR_Sigma_2+RR_Sigma_3+RR_Sigma_4+RR_Sigma_5)+H_3_2+H_4_2+H_4_3+H_5_2+H_5_3+H_5_4;
    return RR;

}
*/

TF1 *fa = new TF1("fa",RR_4,0.5,6.0,15);
TF1 *fa1 = new TF1("fa1",RR_4_1,0.5,6.0,16);

//TF1 *faa = new TF1("fa",RR_41,0,6,15);


//追加で書き加えた部分 end//////////////////

void X_ray_reflect_1(){
    
    TCanvas* c = new TCanvas("c","c",800,600);
    string filepath="center1.ras";
    ifstream ifs(filepath);
    double aa,bb,cc;
    
    int nnn=1500;
    
    double theta[nnn],theta_error[nnn],X_ray_R[nnn],X_ray_R_attenuator[nnn],X_ray_R1[nnn];
    
    int i1=0;
   
    
    
    while(ifs >> aa >> bb >> cc )
    {
        theta[i1]=aa;
        X_ray_R[i1]=bb;
        X_ray_R_attenuator[i1]=cc;
        X_ray_R1[i1]=bb*cc;//(2.e7);
    
        //cout<<theta[i1]<< "_"<<X_ray_R[i1] <<endl;
        i1++;
        if(i1>nnn){
            break;
        }

    }
    
        
    
    //TGraphErrors *grR = new TGraphErrors(nnn,theta,X_ray_R,theta_error,X_ray_R_error);
    
    //double X_ray_R2=X_ray_R1/(2.*e7);
    TGraph *grR = new TGraph(nnn,theta,X_ray_R1);
    gPad->SetLogy();
    //

    grR->Draw("AL");


    //fa->SetParLimits(1,0.0,1.0);
    //fa->SetParLimits(1,0.0,1.0);
/*
    fa->SetParameters(0,5.);
    fa->SetParameters(1,5.);
    fa->SetParameters(2,7.);
    fa->SetParameters(3,0.2);
    fa->SetParameters(4,2);

    fa->SetParameters(5,1.4);
    fa->SetParameters(6,0.7);
    fa->SetParameters(7,90);
    fa->SetParameters(8,4.);
    fa->SetParameters(9,0.);

    fa->SetParameters(10,0.5);
    fa->SetParameters(11,0.);
    fa->SetParameters(12,0.3);
    fa->SetParameters(13,3.);
    fa->SetParameters(14,0.3);
*/
///*固定
    fa1->FixParameter(0,4.95);
    fa1->FixParameter(1,4.55);
    fa1->FixParameter(2,7.);
    fa1->FixParameter(3,0.23);
    fa1->FixParameter(4,2.33);

    fa1->FixParameter(5,1.4);
    fa1->FixParameter(6,0.694);
    fa1->FixParameter(7,86.04);
    fa1->FixParameter(8,4.1);
    fa1->FixParameter(9,380.);

    fa1->FixParameter(10,0.440);
    fa1->FixParameter(11,0.);
    fa1->FixParameter(12,0.3);
    //fa1->FixParameter(13,2.8);
    fa1->SetParLimits(13,0.,100.);
    fa1->FixParameter(14,0.3);
    fa1->SetParLimits(15,1.e5,1.e9);

    fa1->SetNpx(100000);
    //fa->FixParLimits(1,0.0);
//*/

/*
    fa1->SetParLimits(0,3.0,6.0);
    fa1->SetParLimits(1,3.0,6.0);
    fa1->SetParLimits(2,7.0,8.0);
    fa1->SetParLimits(3,0.1,0.3);
    fa1->SetParLimits(4,2.0,3.0);
    fa1->SetParLimits(5,1.0,2.0);
    fa1->SetParLimits(6,0.5,1.0);
    fa1->SetParLimits(7,80.0,95.0);
    fa1->SetParLimits(8,4.0,5.0);
    fa1->SetParLimits(9,0.0,1.0);
    fa1->SetParLimits(10,0.3,1.0);
    fa1->SetParLimits(11,0.0,1.0);
    fa1->SetParLimits(12,0.0,1.0);
    fa1->SetParLimits(13,2.0,3.0);
    fa1->SetParLimits(14,0.0,1.0);
    
*/



    
    grR->Fit("fa1","R","10000",1.,4.0);



    //fa->GetXaxis()->SetLimits(0.5,3.);
    //fa->Draw("AL");
    
   
TCanvas* c1 = new TCanvas("c1","c1",800,600);
const int num=1000;
double xx[num],yy[num];
for(int i=0; i<num; i++){
    xx[i]=0.5+5.5*i/num;
    yy[i]=RR_4a(xx[i]);
    TGraph *gr = new TGraph(num,xx,yy);
    //cout<<yy[i]<<endl;
    gr->Draw("AL");

    //cout<<xxx1<<"_"<<yyy1<<"_"<<yyy2<<endl;
}
  return;
}
