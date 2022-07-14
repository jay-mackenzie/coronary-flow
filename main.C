// makes own tubes from summary table

// best yet convergence test: compares all points in all tubes between current and prev. periods
// NOT RESTARTABLE
#include "coronary.h"
#include "pulmonary.h"
#include "tools.h"
#include "arteries.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;
int main(int argc, char * argv[]) {
    int gridPts, LVPA = 0, max_its, num_to_save;
    double SVPA;
    // int id;
    
    
    double rm, f3, fa3, threshold;
    
//    double rm = 0.002, f1 = 0, f2 = 1, f3 = 2.53*pow(10, 6), fa1 = 0, fa2 = 1, fa3 = pow(10, 6), 
//    asym = 0.41, expo = 2.76, lrrA = 36, lrrV = 36;
//    bool verbosity = 0; int gridPts = 8;
//    double SVPA=0; double LVPA=0;
        
    // gridPts = atoi(argv[1]), rm = atof(argv[2]), f3 = atof(argv[3])*pow(10, 6), fa3 = atof(argv[4])*pow(10, 6),
    // threshold =  atof(argv[5]),max_its = atoi(argv[6]),num_to_save = atoi(argv[7]);
    
    // double t0;
    // char * tree;

    // tree = argv[1]; // select the tree to make
    // t0 = atof(argv[2]); // the number of periods already run
    // bool inflow; // set inflow == 1 for pressure, else for flow
    
    double t0 =0.0;
    // int gridPts;
    char * tree;

    if (argc != 3) //argv[0] is the name of the program, here sor06
    {
        printf("Not enough input arguments, noargc %d and they are %s\n", argc, argv[0]);
        return 1;
    }
    // t0 = atoi(argv[1]);
    // inflow = atoi(argv[2]);
    tree = argv[1]; // select the tree to make
    SVPA = atof(argv[2]);
    // id = atoi(argv[2]);


    // if (t0 < 2.0) {t0 = 0.0;} // if fewer than two periods have been run, just start from 0
    // t0 = t0*Period;

    // can also be set in the header file or via command line args
    gridPts = 10, f3 = 2.53*pow(10, 6), fa3 = 1.03*pow(10, 6), threshold =  100.0 ,max_its = 100,num_to_save = 1; 
    rm = 0.002;
    
    // set the rest of the stiffnesses.
    double f1 = 0, f2 = 1, fa1 = 0, fa2 = 1, fv1 = fa1, fv2 = fa2, fv3 = fa3;
    
    // Workspace used by bound_match
    for (int i = 0; i < 8; i++) fj[i] = new double[8];

    // Workspaces used in junction
    for (int i = 0; i < 12; i++) fmono[i] = new double[12];
    for (int i = 0; i < 18; i++) fjac[i] = new double[18];
    for (int i = 0; i < 24; i++) ftrif[i] = new double[24];

    clock_t c1 = clock(); // Only used when timing the program.
    // mkdir(tree, S_IRWXU); // make a dir to keep the outputs for this tree
    
    // ==================================================================================================
    // tell the user smth useful.
    cout << "Loading the tree in " << tree <<".\n\n"
    << "gridPts = " << gridPts <<  ", f3 = " << f3 << ", fa3 = "<< fa3 << ", threshold = " << threshold << "\n"
    << "SVPA =  "<< SVPA << ".\n";
    // ==================================================================================================
    // Make the objects of class Tube from a file // JAM
    char str[20]; int char_count = 0; FILE * treefile;

    
    
    char fName [200]; 
    
    bool exclamationCheck = 0;
    if(strchr(tree, '/') != NULL)
    {
        sprintf(fName, "%s", tree);
    }
    else
    {
        sprintf(fName, "./NewSumTab/%s.tsv", tree);
    }
    
    
    treefile = fopen(fName, "r");
    
    // char fName [200];  treefile = fopen(fName, "r");
    

    // set nbrves and make enough tubes
    if (treefile) {while (fscanf(treefile, "%s", str) != EOF) { char_count++; }}; fclose(treefile);
    int line_len = 8; nbrves = int(char_count / line_len);
    Tube * Arteries[nbrves];

    double * array[nbrves];
    for (int i = 0; i < nbrves; i++) {array[i] = new double[line_len];}
    
    for (int row = 0; row < nbrves; row++) {
        for (int col = 0; col < line_len; col++) {
            array[row][col] = 0.0;
        }
    }

    treefile = fopen(fName, "r");
    for (int i = 0; i < char_count; i++) {
        int row = float(i) / float(line_len), col = i % line_len;
        fscanf(treefile, "%s", str); array[row][col] = atof(str);
    }
    fclose(treefile); printf("loaded tube data from tree file\n");

    for (int ii = 0; ii < nbrves; ii++) {


        for (int j = 0; j <= 6; j++)
        {
            if (int(array[ii][j]) == array[ii][j])
            {cout << int(array[ii][j]) << "\t";}
            else    {cout << array[ii][j] << "\t";}
        } // print vessels specs to terminal

        //  Parameters required to initiate class Tube (Length,topradius,botradius,LeftDaughter,MidDaughter,RightDaughter,rmin, grid points,
        //                                              init,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV,SVPA,LVPA);
        
        //  the summary table tsv file has columns
        //  index, length, toprad, botrad, LD, MD, RD, init
        //  only arterial inlet and venous outlet need to be given init values -- the rest are automated
            double rm_val;
            if (int(array[ii][4]) != 0){rm_val =0;}else {rm_val = rm;}
            double init_val;
            if (int(array[ii][4]) == 0 && int(array[ii][6]) != 0){init_val =3;} else{init_val =  array[ii][7];}
            
            // if (rm_val == 0){cout << "rm = " << Arteries[int(array[ii][0])]->rm << "\t\t";}
            // else{cout << "rm = " << Arteries[int(array[ii][0])]->rm << "\t";}
            
            // cout << "init = "<<Arteries[int(array[ii][0])]->init<<"\n";

            // cout << "rm = "<<rm_val<<"\n";
            // cout << "rm = "<<rm<<"\n";



            //  in order, handle trifs
            if (int(array[ii][4]) != 0 && int(array[ii][5]) != 0 && int(array[ii][6])!=0)
            {
                Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], Arteries[int(array[ii][4])],Arteries[int(array[ii][5])], Arteries[int(array[ii][6])],rm_val,  gridPts,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, SVPA, LVPA);
            }
            else if (int(array[ii][4]) != 0 && int(array[ii][5]) == 0 && int(array[ii][6])!=0)
            {
                Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], Arteries[int(array[ii][4])],0, Arteries[int(array[ii][6])],rm_val,  gridPts,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, SVPA, LVPA);
            }
            else if (int(array[ii][4]) != 0 && int(array[ii][5]) == 0 && int(array[ii][6])==0)
            {
                Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], Arteries[int(array[ii][4])],0, 0,rm_val,  gridPts,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, SVPA, LVPA);
            }
            else if (int(array[ii][4]) == 0 && int(array[ii][5]) == 0 && int(array[ii][6])!=0)
            {
                Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], 0,0, Arteries[int(array[ii][6])],rm_val,  gridPts,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, SVPA, LVPA);
            }
            else if (int(array[ii][4]) == 0 && int(array[ii][5]) == 0 && int(array[ii][6])==0)
            {
                Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], 0,0, 0,rm_val,  gridPts,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, SVPA, LVPA);
            }
            
            // Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], Arteries[int(array[ii][4])], Arteries[int(array[ii][5])], Arteries[int(array[ii][6])],rm_val,  gridPts,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, SVPA, LVPA);
            
            if (rm_val == 0){cout << "rm = " << Arteries[int(array[ii][0])]->rm << "\t\t";}
            else{cout << "rm = " << Arteries[int(array[ii][0])]->rm << "\t";}
            
            cout << "init = "<<Arteries[int(array[ii][0])]->init<<"\n";
    }
    printf("\nVessels are made\n\n");

    // ===================================================================================================================================================
    
    // printf("rm = %f\n",Arteries[nbrves-1]->rm);

    // exit(1);
    

    looper(Arteries, threshold, plts, k, max_its, num_to_save, tree, t0); // run the model and loop until convergence occurs
    

    fprintf(stdout, "nbrves = %d, Lax, ", nbrves);
    clock_t c2 = clock(); // FIXME clock() may wrap after about 72 min.
    int tsec = (int)((double)(c2 - c1) / CLOCKS_PER_SEC + 0.5);
    fprintf(stdout, "cpu-time %d:%02d\n", tsec / 60, tsec % 60);
    
    for (int i = 0; i < nbrves; i++) delete Arteries[i];
    for (int i = 0; i < 12; i++) delete[] fmono[i];
    for (int i = 0; i < 18; i++) delete[] fjac[i];
    for (int i = 0; i < 24; i++) delete[] ftrif[i];
    for (int i = 0; i < 8; i++) delete[] fj[i];
    
    return 0;
}
