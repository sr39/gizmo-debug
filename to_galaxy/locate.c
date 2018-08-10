
/* 
 * Numerical Recipes 
 */

#include "math.h"
#include "locate.h"

/* 
 * Given an array xx[1..n] and given a value x returns 
 * a value j such that x is between xx[j] and xx[j+1]. 
 * xx must be monotonic. j=0 or j=n is returned to 
 * indicate that x is out of range. 
 * 
 * "A unit-offset array xx is assumed. To use locate with 
 *  a zero-offset array, remember to subtract 1 from the 
 *  address of xx, and also from the returned value j". 
 *
 */

void locate(float *xx, int n, float x, int *j){

  int ju, jm, jl;
  int ascnd;

  /* initialize lower */
  jl = 0; 
  /* and upper limits. */
  ju = n + 1;

  ascnd = (xx[n] >= xx[1]);

  while(ju-jl > 1){       /* if wer are not yet done, */ 
    jm = (ju + jl) >> 1;  /* compute a mid point      */
    if((x >= xx[jm]) == ascnd)
      jl = jm;            /* and replace either the lower limit */
    else 
      ju = jm;            /* or the upper limit, as appropriate */
  }

  if( x == xx[1] ) *j = 1; 
  else if(x == xx[n]) *j = n-1;
  else *j = jl;

  return;
}

void dlocate(double *xx, int n, double x, int *j){

  int ju, jm, jl;
  int ascnd;

  /* initialize lower */
  jl = 0; 
  /* and upper limits. */
  ju = n + 1;

  ascnd = (xx[n] >= xx[1]);

  while(ju-jl > 1){       /* if wer are not yet done, */ 
    jm = (ju + jl) >> 1;  /* compute a mid point      */
    if((x >= xx[jm]) == ascnd)
      jl = jm;            /* and replace either the lower limit */
    else 
      ju = jm;            /* or the upper limit, as appropriate */
  }

  if( x == xx[1] ) *j = 1; 
  else if(x == xx[n]) *j = n-1;
  else *j = jl;

  return;
}
