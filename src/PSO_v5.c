// #include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "pso_bi.h"
#include "pso_ci.h"
#include "pso_ni.h"
#include "pso_mi.h"
#include "pso_di.h"
#include "pso_b.h"
#include "pso_c.h"
#include "pso_n.h"
#include "pso_m.h"
#include "pso_d.h"


static int swarm_size = 10000; // 50000
// PSO setup for the different pools:
void pso_set_BM_settings_bi(pso_settings_t_bi *settings_bi) {
    //Bound pool initial
    /*  M0 */     settings_bi->x_lo[0] =  0;    settings_bi->x_hi[0] =  0.15; 
    /*  k  */     settings_bi->x_lo[1] = 1;     settings_bi->x_hi[1] = 100;
    /*  T2 */     settings_bi->x_lo[2] = 0.5e-6;settings_bi->x_hi[2] =  200e-6; 
    /*  dw */     settings_bi->x_lo[3] =   -720;settings_bi->x_hi[3] =   -690; 
    /*T1obs*/     settings_bi->x_lo[4] = 0.5;   settings_bi->x_hi[4] = 5.0;                          
    /* T2f */     settings_bi->x_lo[5] = 0.020; settings_bi->x_hi[5] = 2.0;
    /*Options*/   settings_bi->size = swarm_size; settings_bi->goal = 1e-50; settings_bi->steps = 250; settings_bi->print_every = 50;
}
void pso_set_BM_settings_b(pso_settings_t_b *settings_b) {
    //Bound pool 
    /*  M0 */     settings_b->x_lo[0] =  0;    settings_b->x_hi[0] =  0.15; 
    /*  k  */     settings_b->x_lo[1] = 1;     settings_b->x_hi[1] = 100;
    /*  T2 */     settings_b->x_lo[2] = 0.5e-6;settings_b->x_hi[2] =  200e-6; 
    /*  dw */     settings_b->x_lo[3] =   -720;settings_b->x_hi[3] =   -690; 
    /*Options*/   settings_b->size = swarm_size; settings_b->goal = 1e-50; settings_b->steps = 250; settings_b->print_every = 50;
}
void pso_set_BM_settings_ci(pso_settings_t_ci *settings_ci) {
    //Amide initial
    /*  M0 */     settings_ci->x_lo[0] =  0;    settings_ci->x_hi[0] =  0.05; 
    /*  k  */     settings_ci->x_lo[1] = 1;     settings_ci->x_hi[1] = 5000;
    /*  T2 */     settings_ci->x_lo[2] = 0.5e-3;settings_ci->x_hi[2] =  100e-3; 
    /*  dw */     settings_ci->x_lo[3] =   1020;settings_ci->x_hi[3] =   1080; 
    /*Options*/   settings_ci->size = swarm_size; settings_ci->goal = 1e-50; settings_ci->steps = 250; settings_ci->print_every = 50;
}
void pso_set_BM_settings_c(pso_settings_t_c *settings_c) {
    //Amide 
    /*  M0 */     settings_c->x_lo[0] =  0;    settings_c->x_hi[0] =  0.05; 
    /*  k  */     settings_c->x_lo[1] = 1;     settings_c->x_hi[1] = 5000;
    /*  T2 */     settings_c->x_lo[2] = 0.5e-3;settings_c->x_hi[2] =  100e-3; 
    /*  dw */     settings_c->x_lo[3] =   1020;settings_c->x_hi[3] =   1080; 
    /*Options*/   settings_c->size = swarm_size; settings_c->goal = 1e-50; settings_c->steps = 250; settings_c->print_every = 50;
}
void pso_set_BM_settings_ni(pso_settings_t_ni *settings_ni) {
    //NOE -3.5ppm initial
    /*  M0 */     settings_ni->x_lo[0] =  0;    settings_ni->x_hi[0] =  0.05; 
    /*  k  */     settings_ni->x_lo[1] = 1;     settings_ni->x_hi[1] = 5000;
    /*  T2 */     settings_ni->x_lo[2] = 0.5e-3;settings_ni->x_hi[2] =  100e-3; 
    /*  dw */     settings_ni->x_lo[3] =   -1080;settings_ni->x_hi[3] =  -1020; 
    /*Options*/   settings_ni->size = swarm_size; settings_ni->goal = 1e-50; settings_ni->steps = 250; settings_ni->print_every = 50;
}
void pso_set_BM_settings_n(pso_settings_t_n *settings_n) {
    //NOE -3.5ppm
    /*  M0 */     settings_n->x_lo[0] =  0;    settings_n->x_hi[0] =  0.05; 
    /*  k  */     settings_n->x_lo[1] = 1;     settings_n->x_hi[1] = 5000;
    /*  T2 */     settings_n->x_lo[2] = 0.5e-3;settings_n->x_hi[2] =  100e-3; 
    /*  dw */     settings_n->x_lo[3] =   -1080;settings_n->x_hi[3] =  -1020; 
    /*Options*/   settings_n->size = swarm_size; settings_n->goal = 1e-50; settings_n->steps = 250; settings_n->print_every = 50;
}
void pso_set_BM_settings_mi(pso_settings_t_mi *settings_mi) {
    //NOE -1.7ppm initial
    /*  M0 */     settings_mi->x_lo[0] =  0;    settings_mi->x_hi[0] =  0.05; 
    /*  k  */     settings_mi->x_lo[1] = 1;     settings_mi->x_hi[1] = 5000;
    /*  T2 */     settings_mi->x_lo[2] = 0.5e-3;settings_mi->x_hi[2] =  100e-3; 
    /*  dw */     settings_mi->x_lo[3] =   -540;settings_mi->x_hi[3] =  -450; 
    /*Options*/   settings_mi->size = swarm_size; settings_mi->goal = 1e-50; settings_mi->steps = 250; settings_mi->print_every = 50;
}
void pso_set_BM_settings_m(pso_settings_t_m *settings_m) {
    //NOE -1.7ppm
    /*  M0 */     settings_m->x_lo[0] =  0;    settings_m->x_hi[0] =  0.05; 
    /*  k  */     settings_m->x_lo[1] = 1;     settings_m->x_hi[1] = 5000;
    /*  T2 */     settings_m->x_lo[2] = 0.5e-3;settings_m->x_hi[2] =  100e-3; 
    /*  dw */     settings_m->x_lo[3] =   -540;settings_m->x_hi[3] =  -450; 
    /*Options*/   settings_m->size = swarm_size; settings_m->goal = 1e-50; settings_m->steps = 250; settings_m->print_every = 50;
}
void pso_set_BM_settings_di(pso_settings_t_di *settings_di) {
    //NOE creatine initial
    /*  M0 */     settings_di->x_lo[0] =  0;    settings_di->x_hi[0] =  0.05; 
    /*  k  */     settings_di->x_lo[1] = 1;     settings_di->x_hi[1] = 5000;
    /*  T2 */     settings_di->x_lo[2] = 0.5e-3;settings_di->x_hi[2] =  100e-3; 
    /*  dw */     settings_di->x_lo[3] =   590;settings_di->x_hi[3] =  690; 
    /*Options*/   settings_di->size = swarm_size; settings_di->goal = 1e-50; settings_di->steps = 250; settings_di->print_every = 50;
}
void pso_set_BM_settings_d(pso_settings_t_d *settings_d) {
    //NOE creatine
    /*  M0 */     settings_d->x_lo[0] =  0;    settings_d->x_hi[0] =  0.05; 
    /*  k  */     settings_d->x_lo[1] = 1;     settings_d->x_hi[1] = 5000;
    /*  T2 */     settings_d->x_lo[2] = 0.5e-3;settings_d->x_hi[2] =  100e-3; 
    /*  dw */     settings_d->x_lo[3] =   590;settings_d->x_hi[3] =  690; 
    /*Options*/   settings_d->size = swarm_size; settings_d->goal = 1e-50; settings_d->steps = 250; settings_d->print_every = 50;
}
double BMsim_bi(double *x, int dim, void *params, double *specs, double *B1inh) {

    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = x[0];  double kb = x[1];  double T2b = x[2];  double wb = x[3];
    double T1obs = x[4]; double T2f = x[5];

    double R2f = 1.0/T2f;    double R1obs = 1/T1obs;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R1f = 0.5;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);

            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));
        
            // Calculate spectrum

            double R1p = Reff + Rexb;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }

    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }

    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};


    double ffweighting[5*63];
      int wx; int wy;
      for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
        ffweighting[wx*63+wy] = ffweights[wy];
        }
      }


        int j;
        for(j = 0; j < 5*63; j++){
            fftemp[j] = Zarray[j]-specs[j];
            if(fftemp[j] < 0){
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }else{
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }
        }


    double ff = 0;
       int k;
       for (k=0; k<5*63; k++)
       {
             ff = ff + ffvec[k];
       }

    return ff;
}
double BMsim_ci(double *x, int dim, void *params, double *specs, double *B1inh, double *MT) {

    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = MT[0]; double kb = MT[1]; double T2b = MT[2];  double wb = MT[3];
    double M0c = x[0];  double kc = x[1];  double T2c = x[2];  double wc = x[3];
    double T1f = 2.0; double T2f = MT[5];
    double R1obs = 1/MT[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);

            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));


            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
      int wx; int wy;
      for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
        ffweighting[wx*63+wy] = ffweights[wy];
      }}


        int j;
        for(j = 0; j < 5*63; j++){
            fftemp[j] = Zarray[j]-specs[j];
            if(fftemp[j] < 0){
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }else{
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }   
        }


    double ff = 0;
       int k;
       for (k=0; k<5*63; k++)
       {
             ff = ff + ffvec[k];
       }


    return ff;
}
double BMsim_ni(double *x, int dim, void *params, double *specs, double *B1inh, double *MT, double *APT) {

    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = MT[0];  double kb = MT[1];  double T2b = MT[2];  double wb = MT[3];
    double M0c = APT[0]; double kc = APT[1]; double T2c = APT[2]; double wc = APT[3];
    double M0n = x[0];   double kn = x[1];   double T2n = x[2];   double wn = x[3];
    double T1f = 2.0; double T2f = MT[5];
    double R1obs = 1/MT[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;
    double R2n = 1.0/T2n;    double dwn=wn*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);

            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));

            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // NOE -3.5ppm

            double Ln = 2.0 * pow( ((kn+R2n)/kn) * pow(w1,2) + pow(kn+R2n,2) ,0.5);
            double AD = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwn,2) + (pow(dw,2)+pow(w1,2)) * R2n/kn + R2n * (kn+R2n));
            double Rexn = M0n * kn * AD / (pow(Ln/2,2)+pow(dw-dwn,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc + Rexn;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
      int wx; int wy;
      for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
        ffweighting[wx*63+wy] = ffweights[wy];
        }
      }


        int j;
        for(j = 0; j < 5*63; j++){
            fftemp[j] = Zarray[j]-specs[j];
            if(fftemp[j] < 0){
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }else{
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }
        }


    double ff = 0;
       int k;
       for (k=0; k<5*63; k++)
       {
             ff = ff + ffvec[k];
       }


    return ff;
}
double BMsim_mi(double *x, int dim, void *params, double *specs, double *B1inh, double *MT, double *APT, double *NOE) {
    
    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = MT[0]; double kb = MT[1]; double T2b = MT[2]; double wb = MT[3];
    double M0c = APT[0];double kc = APT[1];double T2c = APT[2];double wc = APT[3];
    double M0n = NOE[0];double kn = NOE[1];double T2n = NOE[2];double wn = NOE[3];
    double M0m = x[0];  double km = x[1];  double T2m = x[2];  double wm = x[3];
    double T1f = 2.0; double T2f = MT[5];
    double R1obs = 1/MT[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;
    double R2n = 1.0/T2n;    double dwn=wn*2.0*M_PI;
    double R2m = 1.0/T2m;    double dwm=wm*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);

            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));

            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // NOE -3.5ppm

            double Ln = 2.0 * pow( ((kn+R2n)/kn) * pow(w1,2) + pow(kn+R2n,2) ,0.5);
            double AD = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwn,2) + (pow(dw,2)+pow(w1,2)) * R2n/kn + R2n * (kn+R2n));
            double Rexn = M0n * kn * AD / (pow(Ln/2,2)+pow(dw-dwn,2));

            // NOE -1.7ppm

            double Lm = 2.0 * pow( ((km+R2m)/km) * pow(w1,2) + pow(km+R2m,2) ,0.5);
            double AE = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwm,2) + (pow(dw,2)+pow(w1,2)) * R2m/km + R2m * (km+R2m));
            double Rexm = M0m * km * AE / (pow(Lm/2,2)+pow(dw-dwm,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc + Rexn + Rexm;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
      int wx; int wy;
      for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
            ffweighting[wx*63+wy] = ffweights[wy];
        }
      }


        int j;
        for(j = 0; j < 5*63; j++){
            fftemp[j] = Zarray[j]-specs[j];
            if(fftemp[j] < 0){
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }else{
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }
        }


    double ff = 0;
       int k;
       for (k=0; k<5*63; k++)
       {
             ff = ff + ffvec[k];
       }


    return ff;
}
double BMsim_di(double *x, int dim, void *params, double *specs, double *B1inh, double *MT, double *APT, double *NOE, double *NOE2) {
        
    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = MT[0];  double kb = MT[1];  double T2b = MT[2];  double wb = MT[3];
    double M0c = APT[0]; double kc = APT[1]; double T2c = APT[2]; double wc = APT[3];
    double M0d = x[0];   double kd = x[1];   double T2d = x[2];   double wd = x[3];
    double M0n = NOE[0]; double kn = NOE[1]; double T2n = NOE[2]; double wn = NOE[3];
    double M0m = NOE2[0];double km = NOE2[1];double T2m = NOE2[2];double wm = NOE2[3];
    double T1f = 2.0; double T2f = MT[5];
    double R1obs = 1/MT[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;
    double R2d = 1.0/T2d;    double dwd=wd*2.0*M_PI;
    double R2n = 1.0/T2n;    double dwn=wn*2.0*M_PI;
    double R2m = 1.0/T2m;    double dwm=wm*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);


            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));

            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // Amine

            double Ld = 2.0 * pow( ((kd+R2d)/kd) * pow(w1,2) + pow(kd+R2d,2) ,0.5);
            double AC = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwd,2) + (pow(dw,2)+pow(w1,2)) * R2d/kd + R2d * (kd+R2d));
            double Rexd = M0d * kd * AC / (pow(Ld/2,2)+pow(dw-dwd,2));

            // NOE -3.5ppm

            double Ln = 2.0 * pow( ((kn+R2n)/kn) * pow(w1,2) + pow(kn+R2n,2) ,0.5);
            double AD = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwn,2) + (pow(dw,2)+pow(w1,2)) * R2n/kn + R2n * (kn+R2n));
            double Rexn = M0n * kn * AD / (pow(Ln/2,2)+pow(dw-dwn,2));

            // NOE -1.7ppm

            double Lm = 2.0 * pow( ((km+R2m)/km) * pow(w1,2) + pow(km+R2m,2) ,0.5);
            double AE = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwm,2) + (pow(dw,2)+pow(w1,2)) * R2m/km + R2m * (km+R2m));
            double Rexm = M0m * km * AE / (pow(Lm/2,2)+pow(dw-dwm,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc + Rexd + Rexn + Rexm;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
      int wx; int wy;
      for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
            ffweighting[wx*63+wy] = ffweights[wy];
        }
      }


        int j;
        for(j = 0; j < 5*63; j++){
            fftemp[j] = Zarray[j]-specs[j];
            if(fftemp[j] < 0){
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }else{
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }
        }


    double ff = 0;
       int k;
       for (k=0; k<5*63; k++)
       {
             ff = ff + ffvec[k];
       }


    return ff;
}
double BMsim_b(double *x, int dim, void *params, double *specs, double *B1inh, double *MTi, double *APTi, double *NOEi, double *NOE2i, double *Cri) {

    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = x[0];     double kb = x[1];    double T2b = x[2];    double wb = x[3];
    double M0c = APTi[0];  double kc = APTi[1]; double T2c = APTi[2]; double wc = APTi[3];
    double M0d = Cri[0];   double kd = Cri[1];  double T2d = Cri[2];  double wd = Cri[3];
    double M0n = NOEi[0];  double kn = NOEi[1]; double T2n = NOEi[2]; double wn = NOEi[3];
    double M0m = NOE2i[0]; double km = NOE2i[1];double T2m = NOE2i[2];double wm = NOE2i[3];
    double T1f = 2.0; double T2f = MTi[5];
    double R1obs = 1/MTi[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;
    double R2d = 1.0/T2d;    double dwd=wd*2.0*M_PI;
    double R2n = 1.0/T2n;    double dwn=wn*2.0*M_PI;
    double R2m = 1.0/T2m;    double dwm=wm*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);


            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));

            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // Amine

            double Ld = 2.0 * pow( ((kd+R2d)/kd) * pow(w1,2) + pow(kd+R2d,2) ,0.5);
            double AC = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwd,2) + (pow(dw,2)+pow(w1,2)) * R2d/kd + R2d * (kd+R2d));
            double Rexd = M0d * kd * AC / (pow(Ld/2,2)+pow(dw-dwd,2));

            // NOE -3.5ppm

            double Ln = 2.0 * pow( ((kn+R2n)/kn) * pow(w1,2) + pow(kn+R2n,2) ,0.5);
            double AD = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwn,2) + (pow(dw,2)+pow(w1,2)) * R2n/kn + R2n * (kn+R2n));
            double Rexn = M0n * kn * AD / (pow(Ln/2,2)+pow(dw-dwn,2));

            // NOE -1.7ppm

            double Lm = 2.0 * pow( ((km+R2m)/km) * pow(w1,2) + pow(km+R2m,2) ,0.5);
            double AE = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwm,2) + (pow(dw,2)+pow(w1,2)) * R2m/km + R2m * (km+R2m));
            double Rexm = M0m * km * AE / (pow(Lm/2,2)+pow(dw-dwm,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc + Rexd + Rexn + Rexm;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
      int wx; int wy;
      for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
            ffweighting[wx*63+wy] = ffweights[wy];
        }
      }


        int j;
        for(j = 0; j < 5*63; j++){
            fftemp[j] = Zarray[j]-specs[j];
            if(fftemp[j] < 0){
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }else{
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }
        }


    double ff = 0;
    int k;
    for (k=0; k<5*63; k++)
    {
            ff = ff + ffvec[k];
    }


    return ff;
}
double BMsim_c(double *x, int dim, void *params, double *specs, double *B1inh, double *MT, double *APTi, double *NOEi, double *NOE2i, double *Cri) {

    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = MT[0];    double kb = MT[1];   double T2b = MT[2];   double wb = MT[3];
    double M0c = x[0];     double kc = x[1];    double T2c = x[2];    double wc = x[3];
    double M0d = Cri[0];   double kd = Cri[1];  double T2d = Cri[2];  double wd = Cri[3];
    double M0n = NOEi[0];  double kn = NOEi[1]; double T2n = NOEi[2]; double wn = NOEi[3];
    double M0m = NOE2i[0]; double km = NOE2i[1];double T2m = NOE2i[2];double wm = NOE2i[3];
    double T1f = 2.0; double T2f = MT[5];
    double R1obs = 1/MT[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;
    double R2d = 1.0/T2d;    double dwd=wd*2.0*M_PI;
    double R2n = 1.0/T2n;    double dwn=wn*2.0*M_PI;
    double R2m = 1.0/T2m;    double dwm=wm*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);


            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));

            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // Amine

            double Ld = 2.0 * pow( ((kd+R2d)/kd) * pow(w1,2) + pow(kd+R2d,2) ,0.5);
            double AC = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwd,2) + (pow(dw,2)+pow(w1,2)) * R2d/kd + R2d * (kd+R2d));
            double Rexd = M0d * kd * AC / (pow(Ld/2,2)+pow(dw-dwd,2));

            // NOE -3.5ppm

            double Ln = 2.0 * pow( ((kn+R2n)/kn) * pow(w1,2) + pow(kn+R2n,2) ,0.5);
            double AD = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwn,2) + (pow(dw,2)+pow(w1,2)) * R2n/kn + R2n * (kn+R2n));
            double Rexn = M0n * kn * AD / (pow(Ln/2,2)+pow(dw-dwn,2));

            // NOE -1.7ppm

            double Lm = 2.0 * pow( ((km+R2m)/km) * pow(w1,2) + pow(km+R2m,2) ,0.5);
            double AE = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwm,2) + (pow(dw,2)+pow(w1,2)) * R2m/km + R2m * (km+R2m));
            double Rexm = M0m * km * AE / (pow(Lm/2,2)+pow(dw-dwm,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc + Rexd + Rexn + Rexm;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
      int wx; int wy;
      for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
            ffweighting[wx*63+wy] = ffweights[wy];
        }
      }


        int j;
        for(j = 0; j < 5*63; j++){
            fftemp[j] = Zarray[j]-specs[j];
            if(fftemp[j] < 0){
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }else{
                ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
            }
        }


    double ff = 0;
       int k;
       for (k=0; k<5*63; k++)
       {
             ff = ff + ffvec[k];
       }


    return ff;
}
double BMsim_n(double *x, int dim, void *params, double *specs, double *B1inh, double *MT, double *APT, double *NOEi, double *NOE2i, double *Cri) {

    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = MT[0];    double kb = MT[1];   double T2b = MT[2];   double wb = MT[3];
    double M0c = APT[0];   double kc = APT[1];  double T2c = APT[2];  double wc = APT[3];
    double M0d = Cri[0];   double kd = Cri[1];  double T2d = Cri[2];  double wd = Cri[3];
    double M0n = x[0];     double kn = x[1];    double T2n = x[2];    double wn = x[3];
    double M0m = NOE2i[0]; double km = NOE2i[1];double T2m = NOE2i[2];double wm = NOE2i[3];
    double T1f = 2.0; double T2f = MT[5];
    double R1obs = 1/MT[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;
    double R2d = 1.0/T2d;    double dwd=wd*2.0*M_PI;
    double R2n = 1.0/T2n;    double dwn=wn*2.0*M_PI;
    double R2m = 1.0/T2m;    double dwm=wm*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);


            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));

            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // Amine

            double Ld = 2.0 * pow( ((kd+R2d)/kd) * pow(w1,2) + pow(kd+R2d,2) ,0.5);
            double AC = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwd,2) + (pow(dw,2)+pow(w1,2)) * R2d/kd + R2d * (kd+R2d));
            double Rexd = M0d * kd * AC / (pow(Ld/2,2)+pow(dw-dwd,2));

            // NOE -3.5ppm

            double Ln = 2.0 * pow( ((kn+R2n)/kn) * pow(w1,2) + pow(kn+R2n,2) ,0.5);
            double AD = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwn,2) + (pow(dw,2)+pow(w1,2)) * R2n/kn + R2n * (kn+R2n));
            double Rexn = M0n * kn * AD / (pow(Ln/2,2)+pow(dw-dwn,2));

            // NOE -1.7ppm

            double Lm = 2.0 * pow( ((km+R2m)/km) * pow(w1,2) + pow(km+R2m,2) ,0.5);
            double AE = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwm,2) + (pow(dw,2)+pow(w1,2)) * R2m/km + R2m * (km+R2m));
            double Rexm = M0m * km * AE / (pow(Lm/2,2)+pow(dw-dwm,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc + Rexd + Rexn + Rexm;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
      int wx; int wy;
      for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
            ffweighting[wx*63+wy] = ffweights[wy];
        }
      }


    int j;
    for(j = 0; j < 5*63; j++){
        fftemp[j] = Zarray[j]-specs[j];
        if(fftemp[j] < 0){
            ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
        }else{
            ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
        }
    }


    double ff = 0;
       int k;
       for (k=0; k<5*63; k++)
       {
             ff = ff + ffvec[k];
       }


    return ff;
}
double BMsim_m(double *x, int dim, void *params, double *specs, double *B1inh, double *MT, double *APT, double *NOE, double *NOE2i, double *Cri) {

    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = MT[0];    double kb = MT[1];   double T2b = MT[2];   double wb = MT[3];
    double M0c = APT[0];   double kc = APT[1];  double T2c = APT[2];  double wc = APT[3];
    double M0d = Cri[0];   double kd = Cri[1];  double T2d = Cri[2];  double wd = Cri[3];
    double M0n = NOE[0];   double kn = NOE[1];  double T2n = NOE[2];  double wn = NOE[3];
    double M0m = x[0];     double km = x[1];    double T2m = x[2];    double wm = x[3];
    double T1f = 2.0; double T2f = MT[5];
    double R1obs = 1/MT[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;
    double R2d = 1.0/T2d;    double dwd=wd*2.0*M_PI;
    double R2n = 1.0/T2n;    double dwn=wn*2.0*M_PI;
    double R2m = 1.0/T2m;    double dwm=wm*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);


            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));

            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // Amine

            double Ld = 2.0 * pow( ((kd+R2d)/kd) * pow(w1,2) + pow(kd+R2d,2) ,0.5);
            double AC = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwd,2) + (pow(dw,2)+pow(w1,2)) * R2d/kd + R2d * (kd+R2d));
            double Rexd = M0d * kd * AC / (pow(Ld/2,2)+pow(dw-dwd,2));

            // NOE -3.5ppm

            double Ln = 2.0 * pow( ((kn+R2n)/kn) * pow(w1,2) + pow(kn+R2n,2) ,0.5);
            double AD = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwn,2) + (pow(dw,2)+pow(w1,2)) * R2n/kn + R2n * (kn+R2n));
            double Rexn = M0n * kn * AD / (pow(Ln/2,2)+pow(dw-dwn,2));

            // NOE -1.7ppm

            double Lm = 2.0 * pow( ((km+R2m)/km) * pow(w1,2) + pow(km+R2m,2) ,0.5);
            double AE = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwm,2) + (pow(dw,2)+pow(w1,2)) * R2m/km + R2m * (km+R2m));
            double Rexm = M0m * km * AE / (pow(Lm/2,2)+pow(dw-dwm,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc + Rexd + Rexn + Rexm;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
      for(row=0;row<63;row++) {
        Zarray[ind] = Z[row][col];
        ind++;
      }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
    int wx; int wy;
    for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
            ffweighting[wx*63+wy] = ffweights[wy];
        }
    }


    int j;
    for(j = 0; j < 5*63; j++){
        fftemp[j] = Zarray[j]-specs[j];
        if(fftemp[j] < 0){
            ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
        }else{
            ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
        }
    }


    double ff = 0;
    int k;
    for (k=0; k<5*63; k++)
    {
        ff = ff + ffvec[k];
    }


    return ff;
}
double BMsim_d(double *x, int dim, void *params, double *specs, double *B1inh, double *MT, double *APT, double *NOE, double *NOE2, double *Cri) {

    // define vars
    double freq[63] = { // frequency array, in order
        -100000,   -20000,    -10000,    -9000,     -7000,     -5000,     -3000,     -2000,
        -1700,     -1500,     -1400,     -1300,     -1200,     -1150,     -1100,     -1050,
        -1000,     -950,      -900,      -800,      -700,      -650,      -600,      -550,
        -500,      -450,      -400,      -350,      -250,      -100,      -50,       0,
        50,        100,       250,       350,       400,       450,       500,       550,
        600,       650,       700,       800,       900,       950,       1000,      1050,
        1100,      1150,      1200,      1300,      1400,      1500,      1700,      2000,
        3000,      5000,      7000,      9000,      10000,     20000,     100000};

    double B1[5] = {B1inh[0]*0.33, B1inh[0]*0.67, B1inh[0]*1.00, B1inh[0]*1.33, B1inh[0]*1.67};

    double pd = 3.0;
    double y = 267.5153;

    double M0b = MT[0];    double kb = MT[1];   double T2b = MT[2];   double wb = MT[3];
    double M0c = APT[0];   double kc = APT[1];  double T2c = APT[2];  double wc = APT[3];
    double M0d = x[0];     double kd = x[1];    double T2d = x[2];    double wd = x[3];
    double M0n = NOE[0];   double kn = NOE[1];  double T2n = NOE[2];  double wn = NOE[3];
    double M0m = NOE2[0];  double km = NOE2[1]; double T2m = NOE2[2]; double wm = NOE2[3];
    double T1f = 2.0; double T2f = MT[5];
    double R1obs = 1/MT[4];

    double R2f = 1.0/T2f;    double R1f = 1/T1f;
    double R2b = 1.0/T2b;    double dwb=wb*2.0*M_PI;
    double R2c = 1.0/T2c;    double dwc=wc*2.0*M_PI;
    double R2d = 1.0/T2d;    double dwd=wd*2.0*M_PI;
    double R2n = 1.0/T2n;    double dwn=wn*2.0*M_PI;
    double R2m = 1.0/T2m;    double dwm=wm*2.0*M_PI;

    double wref=2.0*M_PI;
    double Z[63][5];

    int RF;
    for(RF=0; RF<5; RF++) { // no of powers
        double w1 = y*B1[RF];

        int f;
        for(f=0; f<63; f++) { // BM sim

            double dw = freq[f] * wref;
            double theta = atan(w1/dw);

            double Reff = R1f * pow(cos(theta),2) + R2f * pow(sin(theta),2);


            // MT

            double Lb = 2.0 * pow( ((kb+R2b)/kb) * pow(w1,2) + pow(kb+R2b,2) ,0.5);
            double AA = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwb,2) + (pow(dw,2) + pow(w1,2)) * R2b/kb + R2b * (kb+R2b));
            double Rexb = M0b * kb * AA / (pow(Lb/2,2) + pow(dw-dwb,2));

            // Amide

            double Lc = 2.0 * pow( ((kc+R2c)/kc) * pow(w1,2) + pow(kc+R2c,2) ,0.5);
            double AB = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwc,2) + (pow(dw,2)+pow(w1,2)) * R2c/kc + R2c * (kc+R2c));
            double Rexc = M0c * kc * AB / (pow(Lc/2,2)+pow(dw-dwc,2));

            // Amine

            double Ld = 2.0 * pow( ((kd+R2d)/kd) * pow(w1,2) + pow(kd+R2d,2) ,0.5);
            double AC = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwd,2) + (pow(dw,2)+pow(w1,2)) * R2d/kd + R2d * (kd+R2d));
            double Rexd = M0d * kd * AC / (pow(Ld/2,2)+pow(dw-dwd,2));

            // NOE -3.5ppm

            double Ln = 2.0 * pow( ((kn+R2n)/kn) * pow(w1,2) + pow(kn+R2n,2) ,0.5);
            double AD = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwn,2) + (pow(dw,2)+pow(w1,2)) * R2n/kn + R2n * (kn+R2n));
            double Rexn = M0n * kn * AD / (pow(Ln/2,2)+pow(dw-dwn,2));

            // NOE -1.7ppm

            double Lm = 2.0 * pow( ((km+R2m)/km) * pow(w1,2) + pow(km+R2m,2) ,0.5);
            double AE = pow(w1,2) / (pow(dw,2) + pow(w1,2)) * (pow(dwm,2) + (pow(dw,2)+pow(w1,2)) * R2m/km + R2m * (km+R2m));
            double Rexm = M0m * km * AE / (pow(Lm/2,2)+pow(dw-dwm,2));

            // Calculate spectrum

            double R1p = Reff + Rexb + Rexc + Rexd + Rexn + Rexm;

            double Zss = pow(cos(theta),2) * (R1obs / R1p);
            Z[f][RF] = (pow(cos(theta),2) - Zss) * exp(-R1p*pd) + Zss;
        }
    }



    // reshape Z
    double Zarray[5*63];
    int row; int col;int ind=0;
    for(col=0;col<5;col++) {
        for(row=0;row<63;row++) {
            Zarray[ind] = Z[row][col];
            ind++;
        }
    }


    // SSq difference
    double ffvec[5*63];
    double fftemp[5*63];
                        // -333, -66.7,-33.3, -30, -23.3,-16.7, -10, -6.7, -5.7,  -5,  -4.7, -4.3,  -4,  -3.8, -3.7, -3.5,
    double ffweights[63] ={ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        // -3.3, -3.2,  -3,  -2.7, -2.3, -2.2,  -2,  -1.8, -1.7, -1.5, -1.3, -1.2, -0.8, -0.3, -0.2,   0,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  0.2,  0.3,  0.8,  1.2,  1.3,  1.5,  1.7,  1.8,   2,   2.2,  2.3,  2.7,   3,   3.2,  3.3,  3.5,
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,
                        //  3.7,  3.8,   4,   4.3,  4.7,   5,   5.7,  6.7,  10,  16.7, 23.3,  30,  33.3, 66.7,  333
                            1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0};

    double ffweighting[5*63];
    int wx; int wy;
    for(wx=0;wx<5;wx++){
        for(wy=0;wy<63;wy++){
            ffweighting[wx*63+wy] = ffweights[wy];
        }
    }


    int j;
    for(j = 0; j < 5*63; j++){
        fftemp[j] = Zarray[j]-specs[j];
        if(fftemp[j] < 0){
            ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
        }else{
            ffvec[j] = ffweighting[j]*pow(fftemp[j],2);
        }
    }


    double ff = 0;
    int k;
    for (k=0; k<5*63; k++)
    {
        ff = ff + ffvec[k];
    }


    return ff;
}


int main(int argc, char *argv[]) {
  int runs;
  for(runs=1;runs<11;runs++){
  

      // read B1inh
      double B1inh[1];
      B1inh[0] = atof(argv[2]);



      // read specs
      double temp;
      FILE *specstxt;
      char* readprefix = "spec_";
      char* filetoread;
      filetoread = malloc(27 * sizeof(char));
      strcpy(filetoread, readprefix);
      strcat(filetoread, argv[1]);

      specstxt = fopen(filetoread, "r");
        double *specsclean;
        specsclean = malloc(5 * 63 * sizeof(double));
        int i;
        for (i=0;i<63*5;i++) {
                fscanf(specstxt,"%lf,",&temp);
                specsclean[i]= temp;
            }
      fclose(specstxt);
      free(filetoread);

      // Add noise to specs
      int nn;
      double noise[5*63];
      double specs[5*63];

      for (nn=0; nn<5*63; nn++){
          noise[nn] = 1.0;//rand() / 2147483647.0 * 0.005 + 0.9975;//1.0;// 
          specs[nn] = noise[nn] * specsclean[nn];
          // printf("%f\n", specs[nn]);
      }



      // define MT objective function
      pso_obj_fun_t_bi obj_fun_bi = BMsim_bi;
      // initialize MT pso settings
      pso_settings_t_bi settings_bi;
      // set the default settings
      pso_set_default_settings_bi(&settings_bi);
      // set the problem specific settings
      pso_set_BM_settings_bi(&settings_bi);
      // initialize MT GBEST solution
      pso_result_t_bi solution_bi;
      // allocate memory for the best position buffer
      solution_bi.gbest_bi = malloc(settings_bi.dim * sizeof(double));


      // run MT optimization algorithm
      pso_solve_bi(obj_fun_bi, NULL, &solution_bi, &settings_bi, specs, B1inh);

      double MTi[7]; int AA;
      for (AA=0; AA<7; AA++){
        MTi[AA] = solution_bi.gbest_bi[AA];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MTi[0]*100, MTi[1], MTi[2]*1000, MTi[3]/300, solution_bi.error_bi);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MTi[4], MTi[5]*1000);



      // define amide objective function
      pso_obj_fun_t_ci obj_fun_ci = BMsim_ci;
      // initialize amide pso settings
      pso_settings_t_ci settings_ci;
      // set the default settings
      pso_set_default_settings_ci(&settings_ci);
      // set the problem specific settings
      pso_set_BM_settings_ci(&settings_ci);
      // initialize amide GBEST solution
      pso_result_t_ci solution_ci;
      // allocate memory for the best position buffer
      solution_ci.gbest_ci = malloc(settings_ci.dim * sizeof(double));


      // run amide optimization algorithm
      pso_solve_ci(obj_fun_ci, NULL, &solution_ci, &settings_ci, specs, B1inh, MTi);

      double APTi[4]; int BB;
      for (BB=0; BB<4; BB++){
        APTi[BB] = solution_ci.gbest_ci[BB];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MTi[0]*100, MTi[1], MTi[2]*1000, MTi[3]/300, solution_bi.error_bi);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APTi[0]*100, APTi[1], APTi[2]*1000, APTi[3]/300, solution_ci.error_ci);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MTi[4], MTi[5]*1000);



      // define NOE objective function
      pso_obj_fun_t_ni obj_fun_ni = BMsim_ni;
      // initialize NOE pso settings
      pso_settings_t_ni settings_ni;
      // set the default settings
      pso_set_default_settings_ni(&settings_ni);
      // set the problem specific settings
      pso_set_BM_settings_ni(&settings_ni);
      // initialize NOE GBEST solution
      pso_result_t_ni solution_ni;
      // allocate memory for the best position buffer
      solution_ni.gbest_ni = malloc(settings_ni.dim * sizeof(double));


      // run NOE optimization algorithm
      pso_solve_ni(obj_fun_ni, NULL, &solution_ni, &settings_ni, specs, B1inh, MTi, APTi);

      double NOEi[4]; int CC;
      for (CC=0; CC<4; CC++){
        NOEi[CC] = solution_ni.gbest_ni[CC];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MTi[0]*100, MTi[1], MTi[2]*1000, MTi[3]/300, solution_bi.error_bi);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APTi[0]*100, APTi[1], APTi[2]*1000, APTi[3]/300, solution_ci.error_ci);
      printf("M0n = %.3f; kn = %.3fHz; T2n = %.3fms; wn = %.3fppm; NOE error = %.10f \n", NOEi[0]*100, NOEi[1], NOEi[2]*1000, NOEi[3]/300, solution_ni.error_ni);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MTi[4], MTi[5]*1000);


      // define NOE17 objective function
      pso_obj_fun_t_mi obj_fun_mi = BMsim_mi;
      // initialize NOE17 pso settings
      pso_settings_t_mi settings_mi;
      // set the default settings
      pso_set_default_settings_mi(&settings_mi);
      // set the problem specific settings
      pso_set_BM_settings_mi(&settings_mi);
      // initialize NOE17 GBEST solution
      pso_result_t_mi solution_mi;
      // allocate memory for the best position buffer
      solution_mi.gbest_mi = malloc(settings_mi.dim * sizeof(double));


      // run NOE 1.7 optimization algorithm
      pso_solve_mi(obj_fun_mi, NULL, &solution_mi, &settings_mi, specs, B1inh, MTi, APTi, NOEi);

      double NOE2i[4]; int DD;
      for (DD=0; DD<4; DD++){
        NOE2i[DD] = solution_mi.gbest_mi[DD];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MTi[0]*100, MTi[1], MTi[2]*1000, MTi[3]/300, solution_bi.error_bi);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APTi[0]*100, APTi[1], APTi[2]*1000, APTi[3]/300, solution_ci.error_ci);
      printf("M0n = %.3f; kn = %.3fHz; T2n = %.3fms; wn = %.3fppm; NOE error = %.10f \n", NOEi[0]*100, NOEi[1], NOEi[2]*1000, NOEi[3]/300, solution_ni.error_ni);
      printf("M0m = %.3f; km = %.3fHz; T2m = %.3fms; wm = %.3fppm; NOE -1.7 error = %.10f \n", NOE2i[0]*100, NOE2i[1], NOE2i[2]*1000, NOE2i[3]/300, solution_mi.error_mi);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MTi[4], MTi[5]*1000);




      // define Cr objective function
      pso_obj_fun_t_di obj_fun_di = BMsim_di;
      // initialize Cr pso settings
      pso_settings_t_di settings_di;
      // set the default settings
      pso_set_default_settings_di(&settings_di);
      // set the problem specific settings
      pso_set_BM_settings_di(&settings_di);
      // initialize Cr GBEST solution
      pso_result_t_di solution_di;
      // allocate memory for the best position buffer
      solution_di.gbest_di = malloc(settings_di.dim * sizeof(double));


      // run creatine optimization algorithm
      pso_solve_di(obj_fun_di, NULL, &solution_di, &settings_di, specs, B1inh, MTi, APTi, NOEi, NOE2i);


      double Cri[4]; int EE;
      for (EE=0; EE<4; EE++){
        Cri[EE] = solution_di.gbest_di[EE];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MTi[0]*100, MTi[1], MTi[2]*1000, MTi[3]/300, solution_bi.error_bi);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APTi[0]*100, APTi[1], APTi[2]*1000, APTi[3]/300, solution_ci.error_ci);
      printf("M0d = %.3f; kd = %.3fHz; T2d = %.3fms; wd = %.3fppm; Cr error = %.10f \n", Cri[0]*100, Cri[1], Cri[2]*1000, Cri[3]/300, solution_di.error_di);
      printf("M0n = %.3f; kn = %.3fHz; T2n = %.3fms; wn = %.3fppm; NOE error = %.10f \n", NOEi[0]*100, NOEi[1], NOEi[2]*1000, NOEi[3]/300, solution_ni.error_ni);
      printf("M0m = %.3f; km = %.3fHz; T2m = %.3fms; wm = %.3fppm; NOE -1.7 error = %.10f \n", NOE2i[0]*100, NOE2i[1], NOE2i[2]*1000, NOE2i[3]/300, solution_mi.error_mi);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MTi[4], MTi[5]*1000);




      // define MT objective function
      pso_obj_fun_t_b obj_fun_b = BMsim_b;
      // initialize MT pso settings
      pso_settings_t_b settings_b;
      // set the default settings
      pso_set_default_settings_b(&settings_b);
      // set the problem specific settings
      pso_set_BM_settings_b(&settings_b);
      // initialize MT GBEST solution
      pso_result_t_b solution_b;
      // allocate memory for the best position buffer
      solution_b.gbest_b = malloc(settings_b.dim * sizeof(double));


      // RErun MT optimization algorithm
      pso_solve_b(obj_fun_b, NULL, &solution_b, &settings_b, specs, B1inh, MTi, APTi, NOEi, NOE2i, Cri);


      double MT[5]; int FF;
      for (FF=0; FF<5; FF++){
        MT[FF] = solution_b.gbest_b[FF];
      }
      MT[4]=MTi[4];
      MT[5]=MTi[5];


      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MT[0]*100, MT[1], MT[2]*1000, MT[3]/300, solution_b.error_b);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APTi[0]*100, APTi[1], APTi[2]*1000, APTi[3]/300, solution_ci.error_ci);
      printf("M0d = %.3f; kd = %.3fHz; T2d = %.3fms; wd = %.3fppm; Cr error = %.10f \n", Cri[0]*100, Cri[1], Cri[2]*1000, Cri[3]/300, solution_di.error_di);
      printf("M0n = %.3f; kn = %.3fHz; T2n = %.3fms; wn = %.3fppm; NOE error = %.10f \n", NOEi[0]*100, NOEi[1], NOEi[2]*1000, NOEi[3]/300, solution_ni.error_ni);
      printf("M0m = %.3f; km = %.3fHz; T2m = %.3fms; wm = %.3fppm; NOE -1.7 error = %.10f \n", NOE2i[0]*100, NOE2i[1], NOE2i[2]*1000, NOE2i[3]/300, solution_mi.error_mi);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MT[4], MT[5]*1000);



      // define amide objective function
      pso_obj_fun_t_c obj_fun_c = BMsim_c;
      // initialize amide pso settings
      pso_settings_t_c settings_c;
      // set the default settings
      pso_set_default_settings_c(&settings_c);
      // set the problem specific settings
      pso_set_BM_settings_c(&settings_c);
      // initialize amide GBEST solution
      pso_result_t_c solution_c;
      // allocate memory for the best position buffer
      solution_c.gbest_c = malloc(settings_c.dim * sizeof(double));


      // run amide optimization algorithm
      pso_solve_c(obj_fun_c, NULL, &solution_c, &settings_c, specs, B1inh, MT, APTi, NOEi, NOE2i, Cri);

      double APT[4]; int GG;
      for (GG=0; GG<4; GG++){
        APT[GG] = solution_c.gbest_c[GG];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MT[0]*100, MT[1], MT[2]*1000, MT[3]/300, solution_b.error_b);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APT[0]*100, APT[1], APT[2]*1000, APT[3]/300, solution_c.error_c);
      printf("M0d = %.3f; kd = %.3fHz; T2d = %.3fms; wd = %.3fppm; Cr error = %.10f \n", Cri[0]*100, Cri[1], Cri[2]*1000, Cri[3]/300, solution_di.error_di);
      printf("M0n = %.3f; kn = %.3fHz; T2n = %.3fms; wn = %.3fppm; NOE error = %.10f \n", NOEi[0]*100, NOEi[1], NOEi[2]*1000, NOEi[3]/300, solution_ni.error_ni);
      printf("M0m = %.3f; km = %.3fHz; T2m = %.3fms; wm = %.3fppm; NOE -1.7 error = %.10f \n", NOE2i[0]*100, NOE2i[1], NOE2i[2]*1000, NOE2i[3]/300, solution_mi.error_mi);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MT[4], MT[5]*1000);

      // define NOE objective function
      pso_obj_fun_t_n obj_fun_n = BMsim_n;
      // initialize NOE pso settings
      pso_settings_t_n settings_n;
      // set the default settings
      pso_set_default_settings_n(&settings_n);
      // set the problem specific settings
      pso_set_BM_settings_n(&settings_n);
      // initialize NOE GBEST solution
      pso_result_t_n solution_n;
      // allocate memory for the best position buffer
      solution_n.gbest_n = malloc(settings_n.dim * sizeof(double));


      // run NOE optimization algorithm
      pso_solve_n(obj_fun_n, NULL, &solution_n, &settings_n, specs, B1inh, MT, APT, NOEi, NOE2i, Cri);

      double NOE[4]; int HH;
      for (HH=0; HH<4; HH++){
        NOE[HH] = solution_n.gbest_n[HH];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MT[0]*100, MT[1], MT[2]*1000, MT[3]/300, solution_b.error_b);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APT[0]*100, APT[1], APT[2]*1000, APT[3]/300, solution_c.error_c);
      printf("M0d = %.3f; kd = %.3fHz; T2d = %.3fms; wd = %.3fppm; Cr error = %.10f \n", Cri[0]*100, Cri[1], Cri[2]*1000, Cri[3]/300, solution_di.error_di);
      printf("M0n = %.3f; kn = %.3fHz; T2n = %.3fms; wn = %.3fppm; NOE error = %.10f \n", NOE[0]*100, NOE[1], NOE[2]*1000, NOE[3]/300, solution_n.error_n);
      printf("M0m = %.3f; km = %.3fHz; T2m = %.3fms; wm = %.3fppm; NOE -1.7 error = %.10f \n", NOE2i[0]*100, NOE2i[1], NOE2i[2]*1000, NOE2i[3]/300, solution_mi.error_mi);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MT[4], MT[5]*1000);

      // define NOE17 objective function
      pso_obj_fun_t_m obj_fun_m = BMsim_m;
      // initialize NOE17 pso settings
      pso_settings_t_m settings_m;
      // set the default settings
      pso_set_default_settings_m(&settings_m);
      // set the problem specific settings
      pso_set_BM_settings_m(&settings_m);
      // initialize NOE17 GBEST solution
      pso_result_t_m solution_m;
      // allocate memory for the best position buffer
      solution_m.gbest_m = malloc(settings_m.dim * sizeof(double));


      // run NOE 1.7 optimization algorithm
      pso_solve_m(obj_fun_m, NULL, &solution_m, &settings_m, specs, B1inh, MT, APT, NOE, NOE2i, Cri);

      double NOE2[4]; int JJ;
      for (JJ=0; JJ<4; JJ++){
        NOE2[JJ] = solution_m.gbest_m[JJ];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MT[0]*100, MT[1], MT[2]*1000, MT[3]/300, solution_b.error_b);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APT[0]*100, APT[1], APT[2]*1000, APT[3]/300, solution_c.error_c);
      printf("M0d = %.3f; kd = %.3fHz; T2d = %.3fms; wd = %.3fppm; Cr error = %.10f \n", Cri[0]*100, Cri[1], Cri[2]*1000, Cri[3]/300, solution_di.error_di);
      printf("M0n = %.3f; kn = %.3fHz; T2n = %.3fms; wn = %.3fppm; NOE error = %.10f \n", NOE[0]*100, NOE[1], NOE[2]*1000, NOE[3]/300, solution_n.error_n);
      printf("M0m = %.3f; km = %.3fHz; T2m = %.3fms; wm = %.3fppm; NOE -1.7 error = %.10f \n", NOE2[0]*100, NOE2[1], NOE2[2]*1000, NOE2[3]/300, solution_m.error_m);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MT[4], MT[5]*1000);

      // define Cr objective function
      pso_obj_fun_t_d obj_fun_d = BMsim_d;
      // initialize Cr pso settings
      pso_settings_t_d settings_d;
      // set the default settings
      pso_set_default_settings_d(&settings_d);
      // set the problem specific settings
      pso_set_BM_settings_d(&settings_d);
      // initialize Cr GBEST solution
      pso_result_t_d solution_d;
      // allocate memory for the best position buffer
      solution_d.gbest_d = malloc(settings_d.dim * sizeof(double));


      // run creatine optimization algorithm
      pso_solve_d(obj_fun_d, NULL, &solution_d, &settings_d, specs, B1inh, MT, APT, NOE, NOE2, Cri);


      double Cr[4]; int KK;
      for (KK=0; KK<4; KK++){
        Cr[KK] = solution_d.gbest_d[KK];
      }

      printf("M0b = %.3f; kb = %.3fHz; T2b = %.3fms; wb = %.3fppm; MT error = %.10f \n", MT[0]*100, MT[1], MT[2]*1000, MT[3]/300, solution_b.error_b);
      printf("M0c = %.3f; kc = %.3fHz; T2c = %.3fms; wc = %.3fppm; Amide error = %.10f \n", APT[0]*100, APT[1], APT[2]*1000, APT[3]/300, solution_c.error_c);
      printf("M0d = %.3f; kd = %.3fHz; T2d = %.3fms; wd = %.3fppm; Cr error = %.10f \n", Cr[0]*100, Cr[1], Cr[2]*1000, Cr[3]/300, solution_d.error_d);
      printf("M0n = %.3f; kn = %.3fHz; T2n = %.3fms; wn = %.3fppm; NOE error = %.10f \n", NOE[0]*100, NOE[1], NOE[2]*1000, NOE[3]/300, solution_n.error_n);
      printf("M0m = %.3f; km = %.3fHz; T2m = %.3fms; wm = %.3fppm; NOE -1.7 error = %.10f \n", NOE2[0]*100, NOE2[1], NOE2[2]*1000, NOE2[3]/300, solution_m.error_m);
      printf("T1f = %.3fs; T2f = %.3fms; \n", MT[4], MT[5]*1000);


      // save some results
      FILE *BMresults;
      char* writeprefix = "results//results1kP100R_";
      char* filetowrite;
      char* runno;
      char* run;
      runno = malloc(5 * sizeof(char));
      run = malloc(3 * sizeof(char));
      sprintf(runno,"%.0f",atof(argv[3]));
      sprintf(run,"%.d",runs);
      filetowrite = malloc(100 * sizeof(char));
      strcpy(filetowrite, writeprefix);
      strcat(filetowrite, runno);
      strcat(filetowrite, "_run");
      strcat(filetowrite, run);
      strcat(filetowrite, "_");
      strcat(filetowrite, argv[1]);
      strcat(filetowrite, ".csv");
      free(runno);

      BMresults = fopen(filetowrite,"w+");


      fprintf(BMresults, "%.5f, %.3f, %.6f, %.3f, ",MT[0],MT[1],MT[2],MT[3]);
      fprintf(BMresults, "%.5f, %.3f, %.6f, %.3f, ",APT[0],APT[1],APT[2],APT[3]);
      fprintf(BMresults, "%.5f, %.3f, %.6f,%.3f, ",Cr[0],Cr[1],Cr[2],Cr[3]);
      fprintf(BMresults, "%.5f, %.3f, %.6f, %.3f, ",NOE[0],NOE[1],NOE[2],NOE[3]);
      fprintf(BMresults, "%.5f, %.3f, %.6f, %.3f, ",NOE2[0],NOE2[1],NOE2[2],NOE2[3]);
      fprintf(BMresults, "%.3f, %.3f, ",MT[4],MT[5]);

      fprintf(BMresults, "%.10f\n",solution_d.error_d);


      fclose(BMresults);
      free(filetowrite);

  

      free(solution_bi.gbest_bi);
      free(solution_ci.gbest_ci);
      free(solution_ni.gbest_ni);
      free(solution_mi.gbest_mi);
      free(solution_di.gbest_di);
      free(solution_b.gbest_b);
      free(solution_c.gbest_c);
      free(solution_n.gbest_n);
      free(solution_m.gbest_m);
      free(solution_d.gbest_d);
      free(specsclean);
  }
  return 0;
}


