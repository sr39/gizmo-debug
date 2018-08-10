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

void locate(float *xx, int n, float x, int *j);
void dlocate(double *xx, int n, double x, int *j);
