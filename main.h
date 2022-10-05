//**************************************************************************/
/*                                                                         */
/*  Program: main.h                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
//**************************************************************************/

// $Id: main.h,v 1.7 2005/07/07 22:19:51 heine Exp $
// Last updated on October 23, 2014 by M. Umar Qureshi
// coronary version
#ifndef _MAIN_H
#define _MAIN_H
#include <cmath>
#include "tools.h"
int nbrves, N_aorta;               // Number of vessels in the tree.
//**************************************************************************/
//    int    n = 15, tmstps = pow(2, n);              // The number of timesteps per period.
//    const char*  CO_filename    = "MynardSmolich.dat";   // Input flow file at the heart. //
//    const char*  small_filename    = "ExtP.dat";   // Input flow file at the heart. //
//**************************************************************************/
// data for smaller time steps
int    n = 16, tmstps = pow(2, n);              // The number of timesteps per period.

// CO_filename, small_filename;


// jay
const char*  CO_filename    = "./Data/JA.dat";   // Input flow file at the heart. //
const char*  small_filename = "./Data/JE.dat";   // External pressure profile

// scott
// const char*  CO_filename    = "./Data/SA.dat";   // Input flow file at the heart. //
// const char*  small_filename = "./Data/SE.dat";   // External pressure profile


// namshad
// const char*  CO_filename    = "./Data/NA.dat";   // Input flow file at the heart. //
// const char*  small_filename = "./Data/NE.dat";   // External pressure profile


        // const char*  CO_filename    = "./Data/corPress_scott.dat";   // Input flow file at the heart. //
        // const char*  small_filename    = "./Data/pres_sub3.dat";   // Input flow file at the heart. //
        // const char*  CO_filename    = "./Data/NamshadPAo.dat";   // Input flow file at the heart. //
        // const char*  small_filename    = "./Data/NamshadPExt.dat";   // Input flow file at the heart. //
//**************************************************************************/
// global parameters
double Tper = 0.70; // The period of one heart beat [s].
int plts = 64; // Number of plots per period.
double asym = 0.41, expo = 2.76, lrrA = 36, lrrV = 36; // structured tree params.
bool verbosity = 0;
double Lr  = 0.25;                  // characteristic vessel radius [cm].
//**************************************************************************/
//  NB: [g/cm/s]/s = 0.1 Pa = 7.5E-3 mmHg
double conv   = 1333.220,             // Conversion from mmHg to SI-units.
       rho    = 1.055,                // Density of blood [g/cm^3].
       mu     = 0.04053,              // Viscosity of blood [g/cm/s].
       mu_pl  = mu,                   // Viscosity of blood [g/cm/s].
       nu     = mu/rho,               // Kinematic viscosity
       Fcst   = 50,                   // Determines friction damping coeff.
       Lr2    = sq(Lr),               // The squared radius [cm2].
       Lr3    = cu(Lr),               // The radius to the third power [cm^3].
       g      = 981.0,                // The gravitational force [cm/s^2].
       q      = 10.0*Lr2,             // The characteristic flow [cm^3/s].
       Fr2    = sq(q)/g/pow(Lr,5),    // The squared Froudes number.
       Re     = q*rho/mu/Lr,          // Reynolds number.
       Period = Tper*q/Lr3,           // The dimension-less period.
       k      = Period/tmstps,        // Length of a timestep.
       Deltat = Period/plts,          // Interval between each point plottet.
       p0     = 0.0,
       bound_thick = sqrt((Period*nu)/(2.0*M_PI))/Lr;  // boundary layer thickness

double *fjac[18],
       *fmono[12],
       *ftrif[24],// Work space used by bound_bif.
       *fj[8];      // Work space used by bound_match.

//void looper (Tube *Arteries[], double threshold, int plts, double k);

#endif
