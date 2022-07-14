# coronary-flow
J. A. Mackenzie.
July 2022.
jay.mackenzie@glasgow.ac.uk

This version of the coronary flow code is pressure driven, and includes external pressure imposed upon vascular beds.

Call like `$ ./main <tree> <SVPA>`

The trees are built from tab-seperated summary tables. These have columns
- vessel index
- length
- proximal radius
- distal radius
- left daughter
- middle daughter
- right daughter
- r_min
- init

Daughters are specified with a positive, integer vessel index
- For one-to-one branching, only the LD is given; RD and MD are 0
- for one-to-two, give LD and RD; MD is 0
- for one-to-three give LD, MD, RD
- for a terminal artery, give RD; LD and MD are 0
- otherwise, set all to 0

The tree can be given via it's path as a command line argument or, if it lives in `./NewSumTab`, just its name.


After the tree is built, a function called *looper* is called. It calls solver until the simulations have converged.
Convergence is assesed by checking the percentage difference between the flow and pressure at the ends of all arteries between sucessive cycles.
It is deemed to have converged if both fall below some given threshold or a maximum number of iterations (max_its) have been exceeded. These are set in main.C. The model is run once more and the output of this is saved.
Threshold and max_its can be set as command line args, as can any other parameter, but remember to update the function that checks their count.

Outputs are saves in the directory Outputs. The format is `[treeName]_[SVPA].2d`

There is one output per run
The matlab file plt.m loops

The columns of the output file are
- vessel id
- time (s)
- space (cm)
- fluid pressure (mmHg)
- flow (ml/s)
- transmural pressure (mmHg)
The matlab file basically works by looping through the files you load, and for each one seperates out the flp, flow, and tmp for each spatial point in each vessel

All trees in NewSumTab can be run using ./runAll.sh
If permission is denied try "chmod +x runAll.sh" and go again.
runAll.sh looks in NewSumTab and runs all the trees with SVPA = 100.
The outputs are saved as normal.