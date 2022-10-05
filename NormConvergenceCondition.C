/* ADDED BY Mitchel J. Colebank
   * Rather than specifying the number of cycles as an input to the function,
   * we want to test to see if the solution has converged. If so, we should exit.*/
  
  int period_counter = 1; // Count the number of periods you have solved for
  double norm_sol = 1e+6;
  double sol_tol  = 1e-5;
  printf("NORM_SOL: %f\n",norm_sol);
  double sol_p1[tmstps],sol_p2[tmstps];
  tend      = Deltat;
  
  
  // SOLVE THE MODEL ONCE
  // Note: Only want to test the pressure at the inlet
  int sol_ID = 0;
  while (tend<=period_counter*Period)
  {
  solver (Arteries, tstart, tend, k, WALL_MODEL);
  sol_p1[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0], WALL_MODEL); // for printing
  sol_p1[sol_ID] *= rho*g*Lr/conv;
  tstart = tend;
  tend   = tend + Deltat; // The current ending time is increased by Deltat.
  sol_ID++;
  }
  
  
  // LOOP FOR CONVERGENCE
  double sse;
  while (norm_sol>=sol_tol)
  {
      sol_ID = 0;
      sse    = 0;
      period_counter++;
      while (tend<=period_counter*Period)
      {
          solver (Arteries, tstart, tend, k, WALL_MODEL);
          sol_p2[sol_ID] = Arteries[0]->P(0,Arteries[0]->Anew[0], WALL_MODEL); // for printing
          sol_p2[sol_ID] *= rho*g*Lr/conv;
          sse = sse+ sq(sol_p1[sol_ID]-sol_p2[sol_ID]);
          tstart = tend;
          tend   = tend + Deltat; // The current ending time is increased by Deltat.
          sol_ID++;
      }
      norm_sol = sse;
      memcpy (sol_p1, sol_p2, sizeof(sol_p2));
      printf("NORM_SOL:%f\n",norm_sol);
  }
  printf("num_cylces:%d\n",period_counter);
