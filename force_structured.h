#ifndef RIGID_FORCE_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define RIGID_FORCE_H

# include "defs.h"

void forceUpdate(int step,int pairs[2][100][3] , int *pairs_now, int ptr_new, int ptr_old, int MaxPairs, vector<SubData>& particle,  double *p_energy);

#endif
