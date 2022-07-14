# coronary-flow
J. A. Mackenzie.
July 2022.
jay.mackenzie@glasgow.ac.uk

This version of the coronary flow code is pressure driven, and includes external pressure imposed upon vascular beds.

Call like
>> ./main <tree> <SVPA>

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

The tree can be given via it's path as a command line argument or, if it lives in ./NewSumTab, just its name.

runAll.sh looks in NewSumTab and runs all the trees with SVPA = 100.

Outputs are saves in the directory Outputs. The format is
> [treeName]_[SVPA].2d
