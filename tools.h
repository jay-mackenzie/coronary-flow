//***************************************************************/
/*                                                              */
/*  Program: tools.h                                            */
/*  Version: 2.0                                                */
/*  By: Mette Olufsen, Math-Tech                                */
/*  Date: 14. Jan. 1997                                         */
/*                                                              */
/*  A header file for the auxilliary functions in the module tools.C. Several of these are taken from Numerical Recipes in C, namely the routines used for finding zeros of a function. Either in several dimensions, namely ludcmp, lubksb, zeroludcmp, or in one dimension, zero_1d.                   */
/*                                                              */
//***************************************************************/

// $Id: tools.h,v 1.4 2005/06/27 19:48:11 heine Exp $

#ifndef _TOOLS_H
#define _TOOLS_H

void   error        (const char *s, const char *u);
void   ludcmp       (double *a[], int n, int indx[], double *d);
void   lubksb       (double *a[], int n, int indx[], double b[]);

void   fwdsub       (double *a[], int n, int indx[], double b[]);

int    zero         (double  *x, int n, double tolx, double tolf,
                     double fvec[], double *fjac[]);


// find zeros, with no matrix inversion
int    zeroNoDcmp   (double  *x, int n, double tolx, double tolf,
                     double fvec[], double *fjac[]);
// find det and exit if zero
void   fjacDet      (double *a[], int n, int indx[], double *d);


bool   zero_1d      (double *x, double f, double df, double tolx);

inline double sq(double x) { return x*x; }
inline double cu(double x) { return x*x*x; }

#endif
