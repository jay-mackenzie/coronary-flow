//  PRESSURE DRIVEN
//*************************************************************//
//  Program: arteries.C
//  Version: 2.?
//  By: Mette Olufsen, Math-Tech, M. Umar Qureshi, UofG
//  Date: 14. Jan. 1997
//  This version: Aug 2021. By J. A. Mackenzie, UofG
//
//  This module can predict the flow and pressure in a tree of elastic vessels as described in IMFUFATEKST NR 297, and D2.1-4. The dependencies of the vessels in the tree must be specified in the main module according to the tree in question (for further details see documentation in problem.default.cxx). This module includes all the functions needed to solve the system of equations. That is the description of all functions in the class containing the vessel (for further details see arteries.h), and in particular the functions needed to solve the system of equations numerically.
//  The module is dependent on the utilities in tools.C, and their corresponding h-files, and also arteries.h that includes the declaration of the vessel-object.
//  This version uses the linear B tube law, takes spatially HOMOGENEOUS external pressures applied to large vessels and strucutred trees, and admits trifurcations and vessel stenoses.
//  JAM's final thesis version. Aug 2021.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include "tools.h"
#include "arteries.h"
#include "junction.h"
#include "nr3.h"
#include "ludcmp.h"
#include "qrdcmp.h"
#include "roots_multidim.h"
using namespace std;

extern int nbrves, tmstps, id;
extern char * CO_filename; // is the in data
extern char * large_filename; // is the in data // JAM. 05/2020
extern char * small_filename; // is the in data // JAM. 05/2020

//  A f90 subroutine that determines the impedance at the root of a structured tree. Used as BC for the tree consisting of the large vessels

extern "C"
void impedance_driver_(int * tmstps, double * Period, double * rho,
    double * mu,
    double * r_root, double * rmin,
    double * y11, double * y12,
    double * y21, double * y22,
    double * Lr, double * q, double * g,
    double * fa1, double * fa2, double * fa3,
    double * fv1, double * fv2, double * fv3,
    double * asym, double * expo,
    double * lrrA, double * lrrV);

// VecDoub vecfunc(VecDoub_I xb, VecDoub_I kB) {
//     //  from bound match:
//     //  kB[38] = IPP[qLnb];
//     //  kB[39] = RD->IPP[qLnb];
//     //  kB[40] = LVPA * IPP[qLnb-1];
//     //  kB[43] = RD->IPP[qLnb-1];

//     Doub PA = kB[17] * (sqrt(xb[0] / kB[18]) - 1.0) + kB[38]; // FlP
//     Doub PV = kB[19] * (sqrt(xb[4] / kB[20]) - 1.0) + kB[39]; // FlP
//     Doub PAh = kB[17] * (sqrt((0.5 * (kB[12] + xb[2])) / kB[18]) - 1.0) + (kB[38] + kB[40]) * 0.5; // FlP
//     Doub PVh = kB[19] * (sqrt((0.5 * (kB[13] + xb[6])) / kB[20]) - 1.0) + (kB[39] + kB[43]) * 0.5; // FlP

//     Doub BhNxb2 = kB[21] *
//         (cu(sqrt(xb[2])) - cu(sqrt(kB[22]))) / (3.0 * sqrt(kB[22]) * kB[23]);

//     Doub Fxb3xb2 = kB[24] * xb[3] / xb[2];

//     Doub dBdx1hNxb2 = (
//         (kB[25] / sqrt(kB[22]) - sqrt(kB[28]) * kB[21] / kB[22]) *
//         (cu(sqrt(xb[2])) - cu(sqrt(kB[22]))) / 3.0 +
//         kB[21] *
//         (kB[28] * kB[27] * (sqrt(cu(xb[2] / kB[22])) - 1.0) -
//             xb[2] * (sqrt(xb[2] / kB[22]) - 1) * kB[25] / kB[21])) * kB[26] / kB[23];

//     Doub RDBhMxb6 = kB[29] *
//         (cu(sqrt(xb[4])) - cu(sqrt(kB[30]))) / (3.0 * sqrt(kB[30]) * kB[23]);

//     Doub RDFxb7xb6 = kB[24] * xb[7] / xb[6];

//     Doub RDdBdx1hMxb6 = (
//         (kB[31] / sqrt(kB[30]) - sqrt(kB[28]) * kB[29] / kB[30]) *
//         (cu(sqrt(xb[4])) - cu(sqrt(kB[30]))) / 3.0 +
//         kB[29] * (kB[28] * kB[33] * (sqrt(cu(xb[4] / kB[30])) - 1.0) -
//             xb[4] * (sqrt(xb[4] / kB[30]) - 1.0) * kB[31] / kB[29])) * kB[32] / kB[23];

//     // residuals
//     VecDoub fvec(8);
//     fvec[0] = kB[1] - xb[0] - kB[0] * xb[3];
//     fvec[1] = kB[2] - xb[1] - kB[0] * (sq(xb[3]) / xb[2] + BhNxb2) + kB[15] * (Fxb3xb2 + dBdx1hNxb2);
//     fvec[2] = kB[3] - xb[1] + kB[4] * PA + kB[5] * PV - (kB[34] + kB[35]) * kB[41]; //  JAM
//     fvec[3] = kB[6] - xb[4] - kB[16] * xb[7];
//     fvec[4] = kB[7] - xb[5] - kB[16] * (sq(xb[7]) / xb[6] + RDBhMxb6) +
//         kB[15] * (RDFxb7xb6 + RDdBdx1hMxb6);
//     fvec[5] = kB[8] - xb[5] + kB[9] * PA + kB[10] * PV - (kB[36] + kB[37]) * kB[41]; //  JAM
//     fvec[6] = kB[11] - 0.5 * xb[3] + kB[4] * PAh + kB[5] * PVh - (kB[34] + kB[35]) * 0.5 * (kB[41] + kB[42]); //  JAM
//     fvec[7] = kB[14] - 0.5 * xb[7] + kB[9] * PAh + kB[10] * PVh - (kB[36] + kB[37]) * 0.5 * (kB[41] + kB[42]); //  JAM

//     return fvec;
// }

//  Methods of class Tube, see arteries.h for description of this.
//
//  The constructor. When an object is made this function will initialize all the attributes of the specific tube. The parameters for the length of the specific vessel, the top and bottom radii, and if applicable pointers to the daughter arteries will be initialized according to the actual parameters passed in the call. If the tube is terminal then the peripheral resistance must be set, and the daughter vessels should be NIL. Otherwise the pointers to the daughter vessels must be given. Further all the work arrays are declared and initialized, and the initial condition for the system equations is applied.

Tube::Tube(double Length,
        double topradius, double botradius,
        Tube * LeftDaughter, Tube * MidDaughter, Tube * RightDaughter,
        double rmin, double points, int init, //double K,
        double f1, double f2, double f3,
        double fa1, double fa2, double fa3,
        double fv1, double fv2, double fv3,
        double asym, double expo,
        double lrrA, double lrrV, double SVPA, double LVPA):
    L(Length), rtop(topradius), rbot(botradius),
    LD(LeftDaughter), MD(MidDaughter), RD(RightDaughter),
    rm(rmin), pts(points), init(init), // K_loss(K),
    ff1(f1), ff2(f2), ff3(f3), SVPA(SVPA), LVPA(LVPA) {

        // Initialization of the basic parameters
        N = int(pts * L);
        h = 1.0 / pts / Lr;

        // Declaration and Initialization of the needed intermediate arrays.
        Qnew = new double[N + 1];
        Anew = new double[N + 1];
        Qold = new double[N + 1];
        Aold = new double[N + 1];
        Qprv = new double[N + 1];
        Aprv = new double[N + 1];
        R1 = new double[N + 1];
        R2 = new double[N + 1];
        S1 = new double[N + 1];
        S2 = new double[N + 1];
        r0 = new double[N + 1];
        r0h = new double[N + 2];
        dr0dx = new double[N + 1];
        dr0dxh = new double[N + 2];
        wom = new double[N + 1];
        A0 = new double[N + 1];
        A0h = new double[N + 2];
        fr = new double[N + 1];
        frh = new double[N + 2];
        dfrdr0 = new double[N + 1];
        dfrdr0h = new double[N + 2];
        p1 = new double[N + 1];
        p1h = new double[N + 2];
        dp1dr0 = new double[N + 1];
        dp1dr0h = new double[N + 2];
        Ah = new double[N];
        Qh = new double[N];
        R1h = new double[N];
        R2h = new double[N];
        S1h = new double[N];
        S2h = new double[N];
        pL = new double[tmstps];
        y11 = new double[tmstps];
        y12 = new double[tmstps];
        y21 = new double[tmstps];
        y22 = new double[tmstps];
        double rgLr = 4.0 / 3.0 / rho / g / Lr;
        double rgLr2 = 4.0 / 3.0 / rho / g / Lr2;

        // Vessel geometry is tabulated and initial conditions are applied
        for (int i = 0; i <= N; i++) {
            r0[i] = rtop * exp(i * log(rbot / rtop) / N) / Lr;
            r0h[i] = rtop * exp((i - 0.5) * log(rbot / rtop) / N) / Lr;
            dr0dx[i] = (log(rbot / rtop) / h / N) * r0[i];
            dr0dxh[i] = (log(rbot / rtop) / h / N) * r0h[i];
            wom[i] = r0[i] * sqrt(2.0 * M_PI / Period / nu);
            A0[i] = M_PI * sq(r0[i]);
            A0h[i] = M_PI * sq(r0h[i]);
            fr[i] = (ff1 * exp(ff2 * r0[i]) + ff3) * rgLr;
            frh[i] = (ff1 * exp(ff2 * r0h[i]) + ff3) * rgLr;
            dfrdr0[i] = ff1 * ff2 * exp(ff2 * r0[i]) * rgLr2;
            dfrdr0h[i] = ff1 * ff2 * exp(ff2 * r0h[i]) * rgLr2;
            p1[i] = fr[i] / M_PI;
            p1h[i] = frh[i] / M_PI;
            dp1dr0[i] = dfrdr0[i] / M_PI;
            dp1dr0h[i] = dfrdr0h[i] / M_PI;
            Qnew[i] = 1.0;
            Anew[i] = A0[i];
        }
        r0h[N + 1] = rtop * exp((N + 0.5) * log(rbot / rtop) / N) / Lr;
        dr0dxh[N + 1] = log(rbot / rtop) / h / N * r0h[N + 1];
        A0h[N + 1] = M_PI * sq(r0h[N + 1]);
        frh[N + 1] = (ff1 * exp(ff2 * r0h[N + 1]) + ff3) * rgLr;
        dfrdr0h[N + 1] = ff1 * ff2 * exp(ff2 * r0h[N + 1]) * rgLr2;
        p1h[N + 1] = frh[N + 1] / M_PI;
        dp1dr0h[N + 1] = dfrdr0h[N + 1] / M_PI;
        Ah05 = A0[N];
        Qh05 = 1.0;

        //  =================  PRESSURE AS INLET CONDITION  =================
        if (init == 1) {
            Ps = new double[tmstps + 1];
            FILE * fi = fopen("corPress_scott.dat", "r"); // Reading in pressure
            for (int i = 0; i <= tmstps; i++) {
                fscanf(fi, "%lf", & Ps[i]);
                Ps[i] = (Ps[i]) * conv / rho / g / Lr; // non-dimensionalise if required
                //              Ps[i]  =80*conv/rho/g/Lr;
            }
            fclose(fi);
        }
        //  ==============  PRESSURE AT OUTLET OF LARGE VEINS  ==============
        if (init == 2) {
            Pout = new double[tmstps + 1];
            for (int i = 0; i <= tmstps; i++) {
                Pout[i] = 5.0 * conv / rho / g / Lr;
            }
        }
        //  ==================================================================

        //  sinusoidal pressure profile applied to large vessels
        IPP = new double[tmstps + 1];
        for (int i = 0; i <= tmstps; i++) {
            IPP[i] = LVPA * sin(i * 2 * M_PI / tmstps) * conv / rho / g / Lr; // sinusoidal
        }

        //  In case of an end-tube evaluate the impedances for the boundary condition. This is done by calling the f90 routine root_imp which calculates the impedance at the root of a structured tree. The underscores is sensitive to the compiler but can be seen at the bottom of the file root_imp.o.

        if (init == 3) {
            //  Pressure applied to SMALL vessels
            imp = new double[tmstps + 1];
            RD -> imp = new double[tmstps + 1];
            FILE * fii = fopen(small_filename, "r");
            for (int i = 0; i <= tmstps; i++) {
                fscanf(fii, "%lf", & imp[i]);
                imp[i] = SVPA * imp[i] * conv / rho / g / Lr; // from file
                //                imp[i] = SVPA*sin(i*2*M_PI/tmstps)* conv / rho / g / Lr; // sinusoidal
            }
            RD -> imp = imp;
            fclose(fii);

            //      fprintf(stdout,"\n Calling f90 subroutines\n"); // JAM suppressed
            impedance_driver_( & tmstps, & Period, & rho, & mu_pl, & rbot, & rmin, y11, y12, y21, y22, & Lr, & q, & g, & fa1, & fa2, & fa3, & fv1, & fv2, & fv3, & asym, & expo, & lrrA, & lrrV);

            //  Initialize the array pL used when determining the convolution in the right boundary condition (see the subroutine bound_right).

            for (int j = 0; j < tmstps; j++) {
                pL[j] = p1[N];
            }
        }
    }

//  The destructor. When the tube-objects terminates, all arrays are deleted, in order to free the memory occupied by the object.

Tube::~Tube() {
    delete[] Anew;
    delete[] Qnew;
    delete[] Aold;
    delete[] Qold;
    delete[] Aprv;
    delete[] Qprv;
    delete[] Ah;
    delete[] Qh;
    delete[] y11;
    delete[] y12;
    delete[] y21;
    delete[] y22;
    delete[] pL;
    delete[] R1h;
    delete[] R2h;
    delete[] S1h;
    delete[] S2h;
    delete[] R1;
    delete[] R2;
    delete[] S1;
    delete[] S2;
    delete[] r0;
    delete[] r0h;
    delete[] dr0dx;
    delete[] dr0dxh;
    delete[] A0;
    delete[] A0h;
    delete[] fr;
    delete[] frh;
    delete[] dfrdr0;
    delete[] dfrdr0h;
    delete[] p1;
    delete[] p1h;
    delete[] dp1dr0;
    delete[] dp1dr0h;
}

//  =================  PLOTTING ROUTINES WITH DIMENSIONS  =================
//  The following functions prints P, Q, A, and F as functions of (x, t). This is done in terms of the re-dimensionalized variables. In this case the functions is plotted for a fixed time, but for all x along the vessel in question. Since the doesn't have to be the first vessel in the tree, it would have some offset from the heart. Which determines the position for x. Therefore there are two arguments passed to this function the time and the offset. Umar modified the print routines to write important parameters in one file. Convenient files for WIA.

void Tube::printPxt(FILE * fd, double t, int offset, int sim_time, int vessel_id) {
    if (offset == 0) fprintf(fd, "\n");
    for (int i = 0; i <= N; i++) {
        /* fprintf (fd, "%13.10f %13.10f %15.10f %15.10f %15.10f %17.10f %15.10f\n",
                    t*Lr3/q, (i+offset)*h*Lr, (P(sim_time, i,Anew[i])+p0)*rho*g*Lr/conv,\
         Qnew[i]*q, Anew[i]*Lr2, \
                    c(i, Anew[i])*Fr2, (P(sim_time, i,Anew[i])+p0-IPP[i])*rho*g*Lr/conv); */

        // fprintf(fd, "%i, %9.8f,%9.8f,%9.8f,%9.8f,%9.8f,%9.8f,%9.8f,%9.8f\n",
        //     vessel_id, t * Lr3 / q, // vessel, time
        //     (i + offset) * h * Lr, // space
        //     (P(sim_time, i, Anew[i]) + p0) * rho * g * Lr / conv, //flp
        //     Qnew[i] * q, // flow
        //     Anew[i] * Lr2, //area
        //     (Qnew[i] * q) / (Anew[i] * Lr2), //fluid velocity
        //     (P(sim_time, i, Anew[i]) + p0 - IPP[sim_time]) * rho * g * Lr / conv, // tmp
        //     c(i, Anew[i]) * Fr2);
        fprintf(fd, "%i, %9.8f,%9.8f,%9.8f,%9.8f,%9.8f\n",
            vessel_id, // vessel 
            t * Lr3 / q, // time
            (i + offset) * h * Lr, // space
            (P(sim_time, i, Anew[i]) + p0) * rho * g * Lr / conv, //flp
            Qnew[i] * q, // flow
            (P(sim_time, i, Anew[i]) + p0 - IPP[sim_time]) * rho * g * Lr / conv // tmp
            );



        //  Time,
        //  Vessel Length,
        //  Fluid Pressure [mmHg],
        //  Flow rate [cm^3/sec], Cross-sectional Area [cm^2],
        //  fluid velocity
        //  transmural pressure [mmHg], Pulse wave velocity [m/sec]
    }
}

//  Further print functions can be added, and they would look similar to the two functions above!  The next function returns the pressure p as a function of a fixed x, and the corresponding cross-sectional area A. The pressure is defined according to the  mathematical model, described in IMFUFATEKST no 297, and D2.1-4.

double Tube::P(int qLnb, int i, double A) {
    double pold = fr[i] * (sqrt(A / A0[i]) - 1.0) + IPP[qLnb];
    // FLP iff A(TMP) Nov. 19.

    return pold;
}


double Tube::dPdA(int i, double A) {
    double pold = fr[i] / (sqrt(A0[i] * A) * 2.0); // Changed by JAM. FEB 2019.

    return pold;
}

double Tube::dPdx1(int i, double A)
// DERIVATIVE OF P WRT x IF WE USE A PROPERLY DEFINED AREA
{
    double pold = (dfrdr0[i] * (sqrt(A / A0[i]) - 1.0) - fr[i] * sqrt(M_PI * A) / A0[i]) *
        dr0dx[i];
    // Changed by JAM. FEB 2019.

    return pold;
}

double Tube::B(int i, double A) {
    //double pold = fr[i]*(sqrt(A0[i]*A)-A0[i])/Fr2;
    double pold = fr[i] * (sqrt(cu(A)) - sqrt(cu(A0[i]))) / (3.0 * sqrt(A0[i]) * Fr2);
    // correct assuming A(TMP)

    return pold;
}

double Tube::Bh(int i, double A) {
    int ip1 = i + 1;
    // double pold =  frh[ip1]*(sqrt(A0h[ip1]*A)-A0h[ip1])/Fr2;
    double pold = frh[ip1] * (sqrt(cu(A)) - sqrt(cu(A0h[ip1]))) /
        (3.0 * sqrt(A0h[ip1]) * Fr2);
    // correct assuming A(TMP)
    return pold;
}

double Tube::dBdx1(int i, double A) {
    double dfr = dfrdr0[i];
    double pold = ((dfr / sqrt(A0[i]) - sqrt(M_PI) * fr[i] / A0[i]) * (sqrt(cu(A)) -
            sqrt(cu(A0[i]))) / 3.0 +
        fr[i] *
        (M_PI * r0[i] * (sqrt(cu(A / A0[i])) - 1.0) - A *
            (sqrt(A / A0[i]) - 1.0) * dfr / fr[i])
    ) * dr0dx[i] / Fr2;

    // NEW TUBE LAW. JAM. JAN 2019. // LINE BREAKS. JAM. APR 2020.

    return pold;
}

double Tube::dBdx1h(int i, double A) {
    int ip1 = i + 1;

    double dfr = dfrdr0h[ip1];

    double pold = (
        (dfr / sqrt(A0h[ip1]) - sqrt(M_PI) * frh[ip1] / A0h[ip1]) *
        (sqrt(cu(A)) - sqrt(cu(A0h[ip1]))) / 3.0 +
        frh[ip1] *
        (M_PI * r0h[ip1] * (sqrt(cu(A / A0h[ip1])) - 1.0) -
            A * (sqrt(A / A0h[ip1]) - 1.0) * dfr / frh[ip1])
    ) * dr0dxh[ip1] / Fr2;

    // NEW TUBE LAW. JAM. JAN 2019.

    return pold;
}

double Tube::dBdAh(int i, double A) {
    int ip1 = i + 1;
    //  double pold = 0.5*frh[ip1]*sqrt(A0h[ip1]/A)/Fr2;
    double pold = frh[ip1] * sqrt(A / A0h[ip1]) / (2.0 * Fr2);
    // NEW TUBE LAW BY JAM. FEB 2019.

    return pold;
}

double Tube::d2BdAdxh(int i, double A) {
    int ip1 = i + 1;
    double dfr = dfrdr0h[ip1];
    //  double pold = (-dfr+1.0/sqrt(A)*(sqrt(M_PI)*frh[ip1]+
    //           sqrt(A0h[ip1])*dfr))*dr0dxh[ip1]/Fr2;
    double pold = (dfr +
        sqrt(A * M_PI) * frh[ip1] / A0h[ip1] -
        sqrt(A / A0h[ip1]) * dfr) * dr0dxh[ip1] / Fr2;

    // NEW TUBE LAW. JAM. FEB 2019. // LINE BREAKS. JAM. APR 2020.
    // CHECK AGAIN
    return pold;
}

//  When determining or checking the step-size (k) the CFL-condition is applied. This is determined according to the result reached from the analysis made using the method of characteristics (See IMFUFATEKST no 297). In this function the minimal step-size fulfilling this condition for this tube is returned.

double Tube::CFL() // The CFL-condition
{
    double minimum = 64000000.0;
    for (int i = 0; i <= N; i++) {
        double c_tmp = c(i, Anew[i]);
        double Vnew = Qnew[i] / Anew[i];
        double temp = min(h / fabs(Vnew - c_tmp),
            h / fabs(Vnew + c_tmp));
        if (temp < minimum) minimum = temp;
    }
    return (minimum);
}

//  When taking a Lax-Wendroff step, the flux of the system must be determined. This is evaluated at i + j/2, and the prediction is given as described in IMFUFATEKST no 297 and D2.1-4. The integer k determines whether we deal with the first or the second component of the vector.

double Tube::Rvec(int k, int i, int j, double Q, double A) {
    if (k == 1) return (Q);
    else
    if (k == 2) return (sq(Q) / A + ((j == 0) ? B(i, A) : Bh(i, A)));
    else error("arteries.cxx", "Call of non-existing vector-component of R");
    return (0);
}

//  Similarly the right hand side of the system of equations must be determined at i + j/2. Also in this case the function is given as stated in the mathematical model, and also in this case k states the needed component of the vector.

double Tube::Svec(int k, int i, int j, double Q, double A) {
    if (k == 1) return (0.0);
    else
    if (k == 2) return (F(Q, A) + ((j == 0) ? dBdx1(i, A) : dBdx1h(i, A)));
    else error("arteries.cxx", "Call of non-existing vector-component of S");
    return (0);
}

//  The solutions of Anew and Qnew are found for all interior points of the vessel at (t+k), where k is the length of the current time-step. This function saves the results in the arrays Anew and Qnew, and the function  is made according to Lax-Wendroff's method as described in IMFUFATEKST  no 297 and D2.1-4.

// void Tube::step(double k, int completedPeriods) {
//     double theta = k / h; // Theta is determined.
//     double gamma = 0.5 * k; // Gamma is determined.

//     for (int i = 0; i <= N; i++) // Remember the values at this time level.
//     {
//         Qold[i] = Qnew[i];
//         Aold[i] = Anew[i];
//     }

//     // Anew and Qnew are predicted at the new time level (t+k).
//     for (int i = 0; i <= N; i++) {
//         R1[i] = Rvec(1, i, 0, Qold[i], Aold[i]);
//         R2[i] = Rvec(2, i, 0, Qold[i], Aold[i]);
//         S1[i] = Svec(1, i, 0, Qold[i], Aold[i]);
//         S2[i] = Svec(2, i, 0, Qold[i], Aold[i]);
//     }

//     for (int i = 0; i < N; i++) {
//         Ah[i] = 0.5 * (Aold[i + 1] + Aold[i]) - 0.5 * theta * (R1[i + 1] - R1[i]) +
//             0.5 * gamma * (S1[i + 1] + S1[i]);
//         Qh[i] = 0.5 * (Qold[i + 1] + Qold[i]) - 0.5 * theta * (R2[i + 1] - R2[i]) +
//             0.5 * gamma * (S2[i + 1] + S2[i]);
//         R1h[i] = Rvec(1, i, 1, Qh[i], Ah[i]);
//         R2h[i] = Rvec(2, i, 1, Qh[i], Ah[i]);
//         S1h[i] = Svec(1, i, 1, Qh[i], Ah[i]);
//         S2h[i] = Svec(2, i, 1, Qh[i], Ah[i]);
//     }
//     for (int i = 1; i < N; i++) {
//         Anew[i] = Aold[i] - theta * (R1h[i] - R1h[i - 1]) + gamma * (S1h[i] + S1h[i - 1]);
//         Qnew[i] = Qold[i] - theta * (R2h[i] - R2h[i - 1]) + gamma * (S2h[i] + S2h[i - 1]);
//     }
// }

void Tube::step(double k, int ii, double t, double t0, char * file_name) {
    double theta = k / h; // Theta is determined.
    double gamma = 0.5 * k; // Gamma is determined.

    if (t0 != 0 && t == t0) {
        if (ii == 0) {printf("restarting\n");}
    } else { //if not restarting
        if (t == t0 + k && ii == 0) {printf("Second time step\n");}

        for (int i = 0; i <= N; i++) {
            // Remember the values at this time level.
            Qold[i] = Qnew[i];
            Aold[i] = Anew[i];

            //  Anew and Qnew are predicted at the new time level (t+k).
            R1[i] = Rvec(1, i, 0, Qold[i], Aold[i]);
            R2[i] = Rvec(2, i, 0, Qold[i], Aold[i]);
            S1[i] = Svec(1, i, 0, Qold[i], Aold[i]);
            S2[i] = Svec(2, i, 0, Qold[i], Aold[i]);
        }

        for (int i = 0; i < N; i++) {
            Ah[i] = 0.5 * (Aold[i + 1] + Aold[i]) - 0.5 * theta * (R1[i + 1] - R1[i]) + 0.5 * gamma * (S1[i + 1] + S1[i]);
            Qh[i] = 0.5 * (Qold[i + 1] + Qold[i]) - 0.5 * theta * (R2[i + 1] - R2[i]) + 0.5 * gamma * (S2[i + 1] + S2[i]);
            R1h[i] = Rvec(1, i, 1, Qh[i], Ah[i]);
            R2h[i] = Rvec(2, i, 1, Qh[i], Ah[i]);
            S1h[i] = Svec(1, i, 1, Qh[i], Ah[i]);
            S2h[i] = Svec(2, i, 1, Qh[i], Ah[i]);
        }

        for (int i = 1; i < N; i++) {
            Anew[i] = Aold[i] - theta * (R1h[i] - R1h[i - 1]) + gamma * (S1h[i] + S1h[i - 1]);
            Qnew[i] = Qold[i] - theta * (R2h[i] - R2h[i - 1]) + gamma * (S2h[i] + S2h[i - 1]);
        }
    }
}

//  ==========================  INPUT CONDITION  ==========================
//  The left boundary (x=0) uses this function to model an inflow into the system. The actual parameter given to the function is the model time. As stated in the mathematical model the constants of the function are chosen in order to ensure a certain CO (specified in main.h). Hence we have the specified value of b. Further the period (dimension-less) is assumed to be Period.
//  ---------------------------  FOR PRESSURE  ----------------------------
double Tube::Ps_init(double t, double k, double Period) // FLUID PRESSURE
{
    if (t < Period) return (Ps[int(t / k)]);
    else // FLUID PRESSURE
        if (t >= Period) return (Ps_init((t - Period), k, Period)); // FLUID PRESSURE
        //     if (t <= Period) return (Ps[int(t/k)]); else                // FLUID PRESSURE
        //     if (t >  Period) return (Ps_init((t-Period),k,Period));     // FLUID PRESSURE
        else return (0);
}
//  ==========================  LEFT BOUNDARY  ==========================
//  Uate of the left boundary at time t. This function uses Q0 to determine the flow rate at the next time-step. From this the value of A is predicted using Lax-Wendroff's numerical scheme. This function is only relevant when the tube is an inlet vessel.

void Tube::bound_left(int qLnb, double t, double k, double Period) {
    double Pnew = Ps_init(t, k, Period); // FLUID PRESSURE
    double P05 = Ps_init(t - k, k, Period); // FLUID PRESSURE

    if (int(t / k) < 0) {
        printf("t/k negative in bound_right\n");
    }

    Anew[0] = A0[0] * (sq(((Pnew - IPP[qLnb]) / fr[0]) + 1.0));

    double A05;
    if (qLnb == 0) {
        A05 = A0[0] * (sq(((P05 - IPP[qLnb]) / fr[0]) + 1.0));
    } else if (qLnb != 0) {
        A05 = A0[0] * (sq(((P05 - 0.5 * (IPP[qLnb] + IPP[qLnb - 1])) / fr[0]) + 1.0));
    } //  if statements ensure that qLnb in the range [0, tmstps]

    double Ahm05 = 2.0 * A05 - Ah[0];
    double Qhm05 = R1h[0] + h * (Anew[0] - Aold[0]) / k;
    double R2hm05 = Rvec(2, -1, 1, Qhm05, Ahm05);
    double S2hm05 = Svec(2, -1, 1, Qhm05, Ahm05);
    Qnew[0] = Qold[0] - k * (R2h[0] - R2hm05) / h + k * (S2h[0] + S2hm05) / 2.0;
}
//  -----------------------------------------------------------------------

//  The value at the right boundary at time t is predicted. NB: This should only be used with terminal vessels, i.e. for vessels that don't bifurcate into further branches. In that situation the bifurcation boundary function should be called instead. Again the procedure specified is given according to the mathematical theory presented in IMFUFATEKST no 297 and D2.1-4.

double Tube::c(int i, double A) // The wave speed through aorta.
{
    //double cnst =  0.5*fr[i]*sqrt(A0[i]/A)/Fr2;
    double cnst = 0.5 * (fr[i] / sqrt(A0[i] * A)) / Fr2;
    return sqrt(cnst);
}

// The matching boundary. // All P quantites should FLP
void Tube::bound_match(int qLnb, double t, double k, double theta, double gamma, int ID, double t0, char * file_name) {
        double PA, PV, PAh, PVh; // FLP. Defined ln. 77.

        // int j = 1, ok = false, ntrial = 50000, qLnb_1 = qLnb - 1;
        int j = 1, ok = false, ntrial = 5000, qLnb_1 = qLnb + 1;
 
        // Make sure that qLnb_1 runs in the interval [0:tmstps-1].
        if (qLnb_1 == (int) tmstps) {qLnb_1 = 0;}

        //  In order to make a relation between P(x_L, t+dt) and Q(x_L, t+dt), and P(x_L, t+dt/2) and Q(x_L,t+dt/2) we need to extract the term involving y[0] (see mathematical derivation). This term corresponds to the peripheral. The remaining terms in the convolution present at the boundary, see mathematical derivation.

        double paterms = 0.0, pvterms = 0.0, paterms_half = 0.0, pvterms_half = 0.0;
        int M = RD -> N;
        double RDtheta = k / RD -> h;

        if (t > Period) {
            //convolution integral
            for (int m = 1; m < tmstps; m++) {
                // Corrected time-stepping from NAH.
                int pindex = (qLnb_1 + tmstps - m) % tmstps;
                int pindex1 = pindex - 1;
                if (pindex == 0) {
                    pindex1 = tmstps - 1;
                }

                double small_vessel_curr = imp[pindex];
                double small_vessel_prev = imp[pindex1];

                paterms = paterms +
                    (pL[pindex] - small_vessel_curr) * y11[m] +
                    ((RD -> pL[pindex]) - small_vessel_curr) * y12[m];

                pvterms = pvterms +
                    (pL[pindex] - small_vessel_curr) * y21[m] +
                    ((RD -> pL[pindex]) - small_vessel_curr) * y22[m]; //YJ

                paterms_half = paterms_half +
                    0.5 * (pL[pindex] + pL[pindex1]) * y11[m] +
                    0.5 * ((RD -> pL[pindex]) + (RD -> pL[pindex1])) * y12[m] -
                    0.5 * (small_vessel_curr + small_vessel_prev) * y11[m] -
                    0.5 * (small_vessel_curr + small_vessel_prev) * y12[m];

                pvterms_half = pvterms_half +
                    0.5 * (pL[pindex] + pL[pindex1]) * y21[m] +
                    0.5 * ((RD -> pL[pindex]) + (RD -> pL[pindex1])) * y22[m] -
                    0.5 * (small_vessel_curr + small_vessel_prev) * y21[m]-
                    0.5 * (small_vessel_curr + small_vessel_prev) * y22[m];
            }
            paterms = k * paterms;
            pvterms = k * pvterms;
            paterms_half = k * paterms_half;
            pvterms_half = k * pvterms_half;
        }

        double k1 = Aold[N] + theta * R1h[N - 1];
        double k2 = Qold[N] + theta * R2h[N - 1] + gamma * S2h[N - 1];
        double k3 = paterms;
        double k4 = k * y11[0];
        double k5 = k * y12[0];
        double k6 = RD -> Aold[M] + RDtheta * (RD -> R1h[M - 1]);
        double k7 = RD -> Qold[M] + RDtheta * (RD -> R2h[M - 1]) + gamma * (RD -> S2h[M - 1]);
        double k8 = pvterms;
        double k9 = k * y21[0];
        double k10 = k * y22[0];
        double k11 = paterms_half - 0.5 * Qh[N - 1];
        double k12 = Ah[N - 1];
        double k13 = RD -> Ah[M - 1];
        double k14 = pvterms_half - 0.5 * (RD -> Qh[M - 1]);
        double k15 = k4; // JAM   //715;//
        double k16 = k5; // JAM   //716;//
        double k17 = k9; // JAM   //717;//
        double k18 = k10; // JAM   //718;//
        //  Unknowns declared, and initial guesses applied

        double xb[8];

        //  xb is a vector used to find the residuals for the boundary match
        if (t <= Period) {
            xb[0] = Ah[N - 1];
            xb[1] = Qh[N - 1];
            xb[2] = Aold[N];
            xb[3] = Qold[N];
            xb[4] = RD -> Ah[M - 1];
            xb[5] = RD -> Qh[M - 1];
            xb[6] = RD -> Aold[M];
            xb[7] = RD -> Qold[M];
        }

        if (t > Period) {
            xb[0] = Aold[N];
            xb[1] = Qold[N];
            xb[2] = Ah05; //Aold[N];
            xb[3] = Qh05; //Qold[N];
            xb[4] = RD -> Aold[M];
            xb[5] = RD -> Qold[M];
            xb[6] = RD -> Ah05; //RD->Aold[M];
            xb[7] = RD -> Qh05; //RD->Qold[M];
        }

        while (j <= ntrial && ok == false) {
            double fvec[8];
            PA = P(qLnb, N, xb[0]); // FLP iff P(i, A[i]) FLP.
            PV = RD -> P(qLnb, M, xb[4]); // FLP iff P(i, A[i]) FLP.
            PAh = P(qLnb, N, 0.5 * (k12 + xb[2])); // FLP iff P(i, A[i]) FLP.
            PVh = RD -> P(qLnb, M, 0.5 * (k13 + xb[6])); // FLP iff P(i, A[i]) FLP.
            double IMP = imp[qLnb]; // Ext P.

            //  Vector of residuals
            fvec[0] = k1 - xb[0] - theta * xb[3];

            fvec[1] = k2 - xb[1] - theta * (sq(xb[3]) / xb[2] + Bh(N, xb[2])) +
                gamma * (F(xb[3], xb[2]) + dBdx1h(N, xb[2]));

            fvec[2] = k3 - xb[1] + k4 * PA + k5 * PV - k15 * IMP - k16 * IMP;

            fvec[3] = k6 - xb[4] - RDtheta * xb[7];

            fvec[4] = k7 - xb[5] - RDtheta * (sq(xb[7]) / xb[6] + RD -> Bh(M, xb[6])) +
                gamma * (RD -> F(xb[7], xb[6]) + RD -> dBdx1h(M, xb[6]));

            fvec[5] = k8 - xb[5] + k9 * PA + k10 * PV - k17 * IMP - k18 * IMP;

            fvec[6] = k11 - xb[3] / 2.0 + k4 * PAh + k5 * PVh - k15 * IMP - k16 * IMP;

            fvec[7] = k14 - xb[7] / 2.0 + k9 * PAh + k10 * PVh - k17 * IMP - k18 * IMP;

            //  The Jacobian
            for (int row = 0; row < 8; row++) {
                for (int col = 0; col < 8; col++) {
                    fj[row][col] = 0.0;
                }
            }

            fj[0][0] = -1.0;
            fj[0][3] = -theta;

            fj[1][1] = -1.0;
            fj[1][2] = theta * (sq(xb[3] / xb[2]) - dBdAh(N, xb[2])) +
                gamma * (dFdA(xb[3], xb[2]) + d2BdAdxh(N, xb[2]));
            fj[1][3] = -2.0 * theta * xb[3] / xb[2] + gamma * dFdQ(xb[2]);

            fj[2][0] = k4 * dPdA(N, xb[0]);
            fj[2][1] = -1.0;
            fj[2][4] = k5 * RD -> dPdA(M, xb[4]);

            fj[3][4] = -1.0;
            fj[3][7] = -RDtheta;

            fj[4][5] = -1.0;
            fj[4][6] = RDtheta * (sq(xb[7] / xb[6]) - RD -> dBdAh(M, xb[6])) +
                gamma * (dFdA(xb[7], xb[6]) + RD -> d2BdAdxh(M, xb[6]));
            fj[4][7] = -2.0 * RDtheta * xb[7] / xb[6] + gamma * dFdQ(xb[6]);

            fj[5][0] = k9 * dPdA(N, xb[0]);
            fj[5][4] = k10 * RD -> dPdA(M, xb[4]);
            fj[5][5] = -1.0;

            fj[6][2] = k4 * dPdA(N, (k12 + xb[2]) / 2.0);
            fj[6][3] = -0.5;
            fj[6][6] = k5 * RD -> dPdA(M, (k13 + xb[6]) / 2.0);

            fj[7][2] = k9 * dPdA(N, (k12 + xb[2]) / 2.0);
            fj[7][6] = k10 * RD -> dPdA(M, (k13 + xb[6]) / 2.0);
            fj[7][7] = -0.5;

            // apply the zero finding function
            double tol = 1.0e-8;
            int ch = zero(xb, 8, tol, tol, fvec, fj);
            if (ch == 1) ok = true;
            j = j + 1;
        }

        Anew[N] = xb[0];
        Qnew[N] = xb[1];
        Ah05 = xb[2];
        Qh05 = xb[3];
        RD -> Anew[M] = xb[4];
        RD -> Qnew[M] = xb[5];
        RD -> Ah05 = xb[6];
        RD -> Qh05 = xb[7];
        pL[qLnb_1] = P(qLnb_1, N, Anew[N]); // FLP as P(N, Anew[N]) is FLP
        RD -> pL[qLnb_1] = RD -> P(qLnb, M, xb[4]); // FLP as P(N, Anew[N]) is FLP

        //      if (j%10 == 0)
        //      {printf("%d: AM=%10.15f, QM=%10.15f \n",j, xb[0],xb[1]);}

        if (j >= ntrial) {

            error("arteries.C", "Root not found in the matching");
            cout << "t = "<< t * Lr3 / q << "t/(k*tmstps) = " << t / (k * tmstps) << "\n";
            printf("AM=%10.15f, QM=%10.15f \n", xb[0], xb[1]);
            printf("AL=%10.15f, QL=%10.15f \n", xb[4], xb[5]);
            printf("rbot=%10.15f \n", rbot);
            printf("Arteries[%i] \n", ID);
            exit(1);
        }
}

//  The value at the junction point at time t is predicted. NB: This should only be done for tubes that do bifurcate into further branches. If this is not the case we have a terminal vessel and bound_right should be called instead. The procedure operates according to the specifications in the mathematical model as a link between this tube and its daughters. Therefore there will be three tubes involved in this function. One problem is however, that the rather complicated system of equations does not converge for all choices of parameters (the peripheral resistance, the top radius, and the bottom radius).

void Tube::call_junc(int qLnb, double theta, double gamma, Tube * Arteries[], int parent) {

    Tube * D1 = Arteries[parent] -> LD;
    Tube * D2 = Arteries[parent] -> MD;
    Tube * D3 = Arteries[parent] -> RD;

    int j = 1;
    int ok = false;
    const int ntrial = 4000;
    double tol = 1.0e-6;
    //  check how many daughter vessels are given for each junction to match, and call the appropriate function from junction.c

    if (D1 != 0 && D2 == 0 && D3 == 0) {
        monofurcation(qLnb, theta, gamma, Arteries, parent, ntrial, tol, j, ok);
    } else if (D1 != 0 && D2 == 0 && D3 != 0) {
        bifurcation(qLnb, theta, gamma, Arteries, parent, ntrial, tol, j, ok);
    } else if (D1 != 0 && D2 != 0 && D3 != 0) {
        trifurcation(qLnb, theta, gamma, Arteries, parent, ntrial, tol, j, ok);
    }
}

//  The right boundary (x=L) uses this function to model an inflow into the system. The actual parameter given to the function is the model time. As stated in the mathematical model the constants of the function are chosen in order to ensure a certain CO (specified in main.h). Hence we have the specified value of b. Further the period (dimension-less) is assumed to be Period.

double Tube::PL_init(double t, double k, double Period) {
    if (t <= Period) return (Pout[int(t / k)]);
    else
    if (t > Period) return (PL_init((t - Period), k, Period));
    else return (0);
}

//  Uate of the right boundary at time t. This function uses PL to determine A at the next time-step. From this the value of Q is predicted using Lax-Wendroff's numerical scheme. This function is only relevant when the tube is an outlet vessel.

void Tube::bound_right(int qLnb, double t, double k, double Period) {
    double Pnew = PL_init(t, k, Period);
    double P05 = PL_init(t - k, k, Period);

    if (int(t / k) < 0)
        printf("t/k negative in bound_right\n");

    Anew[0] = A0[0] * (sq(((Pnew - IPP[qLnb]) / fr[0]) + 1.0));

    //    double A05 = A0[0] * (sq(((P05 - IPP[qLnb]) / fr[0]) + 1.0));
    double A05;
    if (qLnb == 0) {
        A05 = A0[0] * (sq(((P05 - 0.5 * (IPP[qLnb] + IPP[qLnb])) / fr[0]) + 1.0));
    } else if (qLnb != 0) {
        A05 = A0[0] * (sq(((P05 - 0.5 * (IPP[qLnb] + IPP[qLnb - 1])) / fr[0]) + 1.0));
    }

    double Ahm05 = 2.0 * A05 - Ah[0];
    double Qhm05 = R1h[0] + h * (Anew[0] - Aold[0]) / k;
    double R2hm05 = Rvec(2, -1, 1, Qhm05, Ahm05);
    double S2hm05 = Svec(2, -1, 1, Qhm05, Ahm05);

    Qnew[0] = Qold[0] - k * (R2h[0] - R2hm05) / h + k * (S2h[0] + S2hm05) / 2.0;

}

//  calls solver and enforces convergence conditions

void looper(Tube * Arteries[], double threshold, int plts, double k, int max_its, int num_to_save, char * file_name, double t0) {

    int period_counter = 1, sol_ID = 0;
    double tstart = t0, tend = tstart + Deltat, pct[2], max_pct[2], matrix_1[nbrves][plts][2], matrix_2[nbrves][plts][2];
    for (int i = 0; i < 2; i++) {
        max_pct[i] = 1e+6;
    }
    // FILE * file1;char name1[100]; FILE * file2;char name2[100];

    // make file to which we print data
    const char * nameVessel; FILE * fVessel; char v[100]; 
    // sprintf(v, "./%s/%s.2d", file_name, file_name); 
    sprintf(v, "./Outputs/%s_%d.2d", file_name, int(Arteries[0]->SVPA)); 
    nameVessel = v; fVessel = fopen(nameVessel, "w");

    //====================================================================================================

    if (t0 == 0) {
        cout << "start first solve\n";


        while (tend <= period_counter * Period) {
            solver(Arteries, tstart, tend, k, t0, file_name);
            for (int ii = 0; ii < nbrves; ii++) {
                matrix_1[ii][sol_ID][0] = Arteries[ii] -> P((int) fmod(tstart / k, tmstps), 0, Arteries[ii] -> Anew[Arteries[ii] -> N]);
                matrix_1[ii][sol_ID][1] = Arteries[ii] -> Qnew[Arteries[ii] -> N];
            }
            tstart = tend;
            tend = tend + Deltat;
            sol_ID++; // advance one step
        }
        cout << "end first solve\n";



        //====================================================================================================

        if (max_its == 1) {
            max_pct[0] = 0.0;
            max_pct[1] = 0.0;
        } // if only running one cycle, skip to save loop

        while (abs(max_pct[0]) >= threshold || abs(max_pct[1]) >= threshold) {
            sol_ID = 0;
            max_pct[0] = 0.0;
            max_pct[1] = 0.0;
            period_counter++;

            while (tend <= period_counter * Period) {
                solver(Arteries, tstart, tend, k, t0, file_name);

                for (int ii = 0; ii < nbrves; ii++) {
                    matrix_2[ii][sol_ID][0] = Arteries[ii] -> P((int) fmod(tstart / k, tmstps), 0, Arteries[ii] -> Anew[Arteries[ii] -> N]);
                    matrix_2[ii][sol_ID][1] = Arteries[ii] -> Qnew[Arteries[ii] -> N];

                    // if (int(tstart / Period) == tstart / Period) {
                    //     //  save flow and area data from current and previous time steps in two files // higher res than other save
                    //     sprintf(name1, "./%s/sv_1_%i.2d", file_name, ii);
                    //     file1 = fopen(name1, "w");
                    //     Arteries[ii] -> save1(file1);
                    //     fclose(file1);
                    //     sprintf(name2, "./%s/sv_2_%i.2d", file_name, ii);
                    //     file2 = fopen(name2, "w");
                    //     Arteries[ii] -> save2(file1);
                    //     fclose(file2);

                    //     //  save pressure data for terminal art/vein pairs for the convolution integral
                    //     if (Arteries[ii] -> init == 3) {
                    //         char namepL[100];
                    //         sprintf(namepL, "./%s/pL_%d.2d", file_name, ii);
                    //         FILE * fpL = fopen(namepL, "w");
                    //         char nameRDpL[100];
                    //         sprintf(nameRDpL, "./%s/RDpL_%d.2d", file_name, ii);
                    //         FILE * fRDpL = fopen(nameRDpL, "w");
                    //         for (int m = 1; m < tmstps; m++) {
                    //             int qLnb = (int) fmod(tstart / k, tmstps);
                    //             int qLnb_1 = qLnb + 1;
                    //             if (qLnb_1 == (int) tmstps) {
                    //                 qLnb_1 = 0;
                    //             }
                    //             int pindex = (qLnb_1 + tmstps - m) % tmstps;
                    //             fprintf(fpL, "%20.15f\n", Arteries[ii] -> pL[pindex]);
                    //             fprintf(fRDpL, "%20.15f\n", (Arteries[ii] -> RD) -> pL[pindex]);
                    //         }
                    //         fclose(fpL);
                    //         fclose(fRDpL);
                    //     }
                    // }
                }
                tstart = tend;
                tend = tend + Deltat;
                sol_ID++; // advance one step
            }

            // find max pct diffs
            for (int ii = 0; ii < nbrves; ii++) {
                for (int jj = 0; jj < plts; jj++) {
                    for (int kk = 0; kk < 2; kk++) {
                        pct[kk] = 100 * (matrix_1[ii][jj][kk] - matrix_2[ii][jj][kk]) / matrix_2[ii][jj][kk];
                        if (abs(pct[kk]) > abs(max_pct[kk])) {
                            max_pct[kk] = pct[kk];
                        }
                    }
                }
            }

            memcpy(matrix_1, matrix_2, sizeof(matrix_2)); // replace previous period's matrix with the current period's

            cout << "completed periods\t\t: " << period_counter << "\n";
            cout << "Pressure:\tpct max diff\t: " << max_pct[0] << "\n";
            cout << "Flow:\t\tpct max diff\t: " << max_pct[1] << "\n\n\n";

            if (period_counter == max_its) {
                max_pct[0] = 0.0;
                max_pct[1] = 0.0;
                num_to_save = 1;
                cout << "Maximum number of permissable itterations exceeded\n";
            }
        }

        cout << "num cycles to convergence: " << period_counter << "\n";

        //====================================================================================================

        //  after convergence, run a given number of extra times, and save
        cout << "Start saving\nNow exiting\n";
        // exit(1);
        while (tend <= (period_counter + num_to_save) * Period) {
            solver(Arteries, tstart, tend, k, t0, file_name);
            for (int j = 0; j < nbrves; j++) {
                Arteries[j] -> printPxt(fVessel, tstart, 0, fmod(tstart / k, tmstps), j);
            }
            tstart = tend;
            tend = tend + Deltat; // move to the next time interval of length Deltat.
        }

    } 
    else if (t0 != 0) {

        // first restart period
        printf("restart %20.15f\n", 1.0 / 7.0);
        tstart = t0;
        tend = tstart + Deltat;

        while (tend <= t0 + Period) {
            solver(Arteries, tstart, tend, k, t0, file_name);
            for (int j = 0; j < nbrves; j++) {Arteries[j] -> printPxt(fVessel, tstart, 0, fmod(tstart / k, tmstps), j);}

            tstart = tend;
            tend = tend + Deltat;

        }
    }

}

//  Solves the non-linear E's (momentum and continuity eqn's.) from t = tstart to t= tend. This function checks the maximal possible size of the next time-step, reduces it to make sure that we don't walk to far, and takes the next step. This is done by executing the step routine, then uating the left boundary and finally uating bifurcation points and the right boundaries. This is carried out as long as the time hasn't passed the desired ending time (tend) which is passed to the function as a parameter.

void solver(Tube * Arteries[], double tstart, double tend, double k, double t0, char * file_name) {
    // The following definitions only used with variable time-stepping
    double t = tstart;
    int qLnb = (int) fmod(t / k, tmstps);

    // As long as we haven't passed the desired ending time do:
    while (t < tend) {
        // Check that the step we take is valid. If this error occurs when forcing a constant step-size the program must be terminated.
        if (t + k > tend) {
            double kold = k;
            k = tend - t;
            printf("ERROR (arteries.C): Step-size changed:\n");
            printf("t+k=%10.15f, tend=%10.15f k=%10.15f kold=%10.15f\n",
                t + kold, tend, k, kold);
        }

        // Check that the CFL-condition applies.
        for (int i = 0; i < nbrves; i++) {
            if (k > Arteries[i] -> CFL()) {
                cout << "Arteries[" << i << "]" << endl;
                cout << "time = " << (t * Lr3 / q) << "s" << endl;
                error("arteries.C", "Time-step size too large\n CFL-condition violated\n Now exiting.\n");
                //                exit(2); // suppress exit
            }
        }

        // solve for interior points, by calling step.
        // for (int i = 0; i < nbrves; i++) {Arteries[i] -> step(k, t0);}
        for (int i = 0; i < nbrves; i++) {
            Arteries[i] -> step(k, i, t, t0, file_name);
        }

        // Update left and right boundaries, and the bifurcation points.
        for (int i = 0; i < nbrves; i++) {
            if (Arteries[i] -> init == 1) {
                Arteries[i] -> bound_left(qLnb, t + k, k, Period);                
            };
            if (Arteries[i] -> init == 2) {
                Arteries[i] -> bound_right(qLnb, t + k, k, Period);
            };
            if (Arteries[i] -> init == 3) {
                double theta = k / Arteries[i] -> h, gamma = k / 2.0;
                Arteries[i] -> bound_match(qLnb, t, k, theta, gamma, i, t0, file_name);
            }
            
            if (Arteries[i] -> rm == 0) {
                double theta = k / Arteries[i] -> h, gamma = k / 2.0;
                Arteries[i] -> call_junc(qLnb, theta, gamma, Arteries, i);
            }
        }
        t = t + k; // Uate the time and position within one period.
        qLnb = (qLnb + 1) % tmstps;
    }
}