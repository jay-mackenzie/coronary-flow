
#include "main.h"
#include "tools.h"
#include "arteries.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

using namespace std;
int main(int argc, char * argv[]) {
    double tstart, tend, finaltime;

    
    // Workspace used by bound_match
    for (int i = 0; i < 8; i++) fj[i] = new double[8];

    // Workspaces used in junction
    for (int i = 0; i < 12; i++) fmono[i] = new double[12];
    for (int i = 0; i < 18; i++) fjac[i] = new double[18];
    for (int i = 0; i < 24; i++) ftrif[i] = new double[24];

    clock_t c1 = clock(); // Only used when timing the program.

    char nameA0 [20]; sprintf(nameA0, "_A0_%d.2d", NHB); FILE *fA0 = fopen (nameA0, "w");
    char nameA1 [20]; sprintf(nameA1, "_A1_%d.2d", NHB); FILE *fA1 = fopen (nameA1, "w");
    char nameA2 [20]; sprintf(nameA2, "_A2_%d.2d", NHB); FILE *fA2 = fopen (nameA2, "w");
    char nameA3 [20]; sprintf(nameA3, "_A3_%d.2d", NHB); FILE *fA3 = fopen (nameA3, "w");
    char nameA4 [20]; sprintf(nameA4, "_A4_%d.2d", NHB); FILE *fA4 = fopen (nameA4, "w");
    char nameA5 [20]; sprintf(nameA5, "_A5_%d.2d", NHB); FILE *fA5 = fopen (nameA5, "w");
    char nameA6 [20]; sprintf(nameA6, "_A6_%d.2d", NHB); FILE *fA6 = fopen (nameA6, "w");
    char nameA7 [20]; sprintf(nameA7, "_A7_%d.2d", NHB); FILE *fA7 = fopen (nameA7, "w");
    char nameA8 [20]; sprintf(nameA8, "_A8_%d.2d", NHB); FILE *fA8 = fopen (nameA8, "w");
    char nameA9 [20]; sprintf(nameA9, "_A9_%d.2d", NHB); FILE *fA9 = fopen (nameA9, "w");
    char nameA10 [20]; sprintf(nameA10, "_A10_%d.2d", NHB); FILE *fA10 = fopen (nameA10, "w");
    char nameA11 [20]; sprintf(nameA11, "_A11_%d.2d", NHB); FILE *fA11 = fopen (nameA11, "w");
    char nameA12 [20]; sprintf(nameA12, "_A12_%d.2d", NHB); FILE *fA12 = fopen (nameA12, "w");
    char nameA13 [20]; sprintf(nameA13, "_A13_%d.2d", NHB); FILE *fA13 = fopen (nameA13, "w");
    char nameA14 [20]; sprintf(nameA14, "_A14_%d.2d", NHB); FILE *fA14 = fopen (nameA14, "w");
    char nameA15 [20]; sprintf(nameA15, "_A15_%d.2d", NHB); FILE *fA15 = fopen (nameA15, "w");

    nbrves    = 16;             // Total number of large arteries and veins in the network -- (temp)
    tstart    = 0.0;            // Starting time.
    finaltime = NHB*Period;     // Final end-time during a simulation.
    tend      = (NHB-1)*Period; // Timestep before the first plot-point
                                  // is reached.

    Tube   *Arteries[nbrves];                    // Array of blood vessels.
    
    // Parameters required to initiate class Tube (Length,topradius,botradius,LeftDaughter,RightDaughter,rmin, points,
                                               //init,K,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);
    
    Arteries[15] = new Tube(
                                            1.21, 0.08, 0.08,
               //                             1.21, 0.08, 0.08,
                            0, 0, 0, rm, 6, 2, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[14] = new Tube(
                                            1.65, 0.1, 0.1,
               //                             1.65, 0.1, 0.1,
                            0, 0, 0, rm, 6, 2, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[13] = new Tube(
                                            1.13, 0.07, 0.07,
               //                             1.13, 0.07, 0.07,
                            0, 0, 0, rm, 6, 2, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[12] = new Tube(
                                            1.25, 0.08, 0.08,
               //                             1.25, 0.08, 0.08,
                            0, 0, 0, rm, 6, 2, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[11] = new Tube(
                                            1.13, 0.07, 0.07,
               //                             1.13, 0.07, 0.07,
                            0, 0, 0, rm, 6, 2, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[10] = new Tube(
                                            0.87, 0.05, 0.05,
               //                             0.87, 0.05, 0.05,
                            0, 0, 0, rm, 7, 2, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[9] = new Tube(
                                            1.5, 0.09, 0.08,
               //                             1.5, 0.09, 0.08,
                            0, 0, Arteries[15], rm, 6, 3, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[8] = new Tube(
                                            0.79, 0.1, 0.1,
               //                             0.79, 0.1, 0.1,
                            0, 0, Arteries[14], rm, 8, 3, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[7] = new Tube(
                                            0.79, 0.15, 0.09,
               //                             0.79, 0.15, 0.09,
                            Arteries[8], 0, Arteries[9], 0, 8, 0, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[6] = new Tube(
                                            1.34, 0.15, 0.07,
               //                             1.34, 0.15, 0.07,
                            0, 0, Arteries[13], rm, 6, 3, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[5] = new Tube(
                                            3.17, 0.15, 0.15,
               //                             3.17, 0.15, 0.15,
                            Arteries[6], 0, Arteries[7], 0, 6, 0, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[4] = new Tube(
                                            3.57, 0.11, 0.08,
               //                             3.57, 0.11, 0.08,
                            0, 0, Arteries[12], rm, 6, 3, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[3] = new Tube(
                                            8.72, 0.14, 0.07,
               //                             8.72, 0.14, 0.07,
                            0, 0, Arteries[11], rm, 6, 3, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[2] = new Tube(
                                            5.48, 0.14, 0.05,
               //                             5.48, 0.14, 0.05,
                            0, 0, Arteries[10], rm, 6, 3, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[1] = new Tube(
                                            0.75, 0.14, 0.11,
               //                             0.75, 0.14, 0.11,
                            Arteries[4], 0, Arteries[5], 0, 8, 0, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    Arteries[0] = new Tube(
                                            0.83, 0.18, 0.14,
               //                             0.83, 0.18, 0.14,
                            Arteries[1], Arteries[2], Arteries[3], 0, 8, 1, 0,
                            f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);

    
    tstart = 0;
    tend = Deltat;

    if (verbosity){fprintf (stdout,"plots start\n");}

    while (tend <= finaltime)
    {
        for (int j=0; j<nbrves; j++)
        {
            int ArtjN = Arteries[j]->N;
            for (int i=0; i<ArtjN; i++)
            {
                Arteries[j]->Qprv[i+1] = Arteries[j]->Qnew[i+1];
                Arteries[j]->Aprv[i+1] = Arteries[j]->Anew[i+1];
            }
        }
        solver (Arteries, tstart, tend, k);
        if (verbosity){fprintf (stdout,".");}
        
        Arteries[0] -> printPxt (fA0, tend, 0,fmod(tstart/k, tmstps));
        Arteries[1] -> printPxt (fA1, tend, 0,fmod(tstart/k, tmstps));
        Arteries[2] -> printPxt (fA2, tend, 0,fmod(tstart/k, tmstps));
        Arteries[3] -> printPxt (fA3, tend, 0,fmod(tstart/k, tmstps));
        Arteries[4] -> printPxt (fA4, tend, 0,fmod(tstart/k, tmstps));
        Arteries[5] -> printPxt (fA5, tend, 0,fmod(tstart/k, tmstps));
        Arteries[6] -> printPxt (fA6, tend, 0,fmod(tstart/k, tmstps));
        Arteries[7] -> printPxt (fA7, tend, 0,fmod(tstart/k, tmstps));
        Arteries[8] -> printPxt (fA8, tend, 0,fmod(tstart/k, tmstps));
        Arteries[9] -> printPxt (fA9, tend, 0,fmod(tstart/k, tmstps));
        Arteries[10] -> printPxt (fA10, tend, 0,fmod(tstart/k, tmstps));
        Arteries[11] -> printPxt (fA11, tend, 0,fmod(tstart/k, tmstps));
        Arteries[12] -> printPxt (fA12, tend, 0,fmod(tstart/k, tmstps));
        Arteries[13] -> printPxt (fA13, tend, 0,fmod(tstart/k, tmstps));
        Arteries[14] -> printPxt (fA14, tend, 0,fmod(tstart/k, tmstps));
        Arteries[15] -> printPxt (fA15, tend, 0,fmod(tstart/k, tmstps));

        tstart = tend;
        tend   = tend + Deltat; // The current ending time is increased by Deltat.
    }

    if (verbosity){fprintf(stdout,"\n");}
    
    fprintf(stdout,"nbrves = %d, Lax, ", nbrves);
    clock_t c2 = clock(); // FIXME clock() may wrap after about 72 min.
    int tsec = (int) ((double) (c2-c1)/CLOCKS_PER_SEC + 0.5);
    fprintf(stdout,"cpu-time %d:%02d\n", tsec / 60, tsec % 60);
    
    if (verbosity) {fprintf(stdout,"\n");}
    
    for (int i=0; i<nbrves; i++) delete Arteries[i];
    for (int i=0; i<12; i++) delete[] fmono[i];
    for (int i=0; i<18; i++) delete[] fjac[i];
    for (int i=0; i<24; i++) delete[] ftrif[i];
    for (int i=0; i<8;  i++) delete[] fj[i];
    return 0;
}
