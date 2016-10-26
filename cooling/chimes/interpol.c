#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "proto.h"

/*
 * ----------------------------------------------------------------------
 * This routine returns the position i of a value x in a 1D table and the
 * displacement dx needed for the interpolation.  The table is assumed to
 * be evenly spaced.
 * ----------------------------------------------------------------------
 */

void get_index_1d_mydbl(double *table, int ntable, double x, int *i, double *dx)
{
  double dxm1;

  dxm1 = (double) (ntable - 1) / (table[ntable - 1] - table[0]);

  if((double) x <= table[0])
    {
      *i = 0;
      *dx = 0;
    }
  else if((double) x >= table[ntable - 1])
    {
      *i = ntable - 2;
      *dx = 1;
    }
  else
    {
      *i = (int) floor(((double) x - table[0]) * dxm1);
      *dx = ((double) x - table[*i]) * dxm1;
    }
}

/*
 * ----------------------------------------------------------------------
 * This routine returns the position i of a value x in a 1D table and the
 * displacement dx needed for the interpolation, but without assuming 
 * that the table is evenly spaced
 * ----------------------------------------------------------------------
 */
 
void get_index_1d_irregular(double *table, int ntable, double x, int *i, double *dx)
{
	if ((double) x <= table[0])
	{
		*i = 0;
		*dx = 0;
	}
	else if ((double) x >= table[ntable - 1])
	{
		*i = ntable - 2;
		*dx = 1;
	}
	else
	{
		*i = 0;
		while (table[*i] < x)
			*i += 1;
		*i -= 1;
		*dx = ((double) x - table[*i]) / (table[(*i) + 1] - table[*i]);
	}
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_1d_mydbl(double *table, int i, double dx)
{
  double result;

  if (dx <= 0.0)
    result = table[i];
  else if (dx >= 1.0)
    result = table[i + 1];
  else
    result = (1 - dx) * table[i] + dx * table[i + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a linear interpolation. It is 
 * to be used when the table is in single precision, 
 * but you want to cast the result to double precision.
 * ----------------------------------------------------------------------
 */

double interpol_1d_fltdbl(float *table, int i, double dx)
{
  double result;

  if (dx <= 0.0)
    result = (double) table[i];
  else if (dx >= 1.0)
    result = (double) table[i + 1];
  else
    result = (1 - dx) * ((double) table[i]) + dx * ((double) table[i + 1]);

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_2d_mydbl(double **table, int i, int j, double dx, double dy)
{
  double result;

  result =
    (1 - dx) * (1 - dy) * table[i][j] +
    (1 - dx) * dy * table[i][j + 1] +
    dx * (1 - dy) * table[i + 1][j] +
    dx * dy * table[i + 1][j + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a bi-linear interpolation. It is 
 * to be used when the table is in single precision, 
 * but you want to cast the result to double precision.
 * ----------------------------------------------------------------------
 */

double interpol_2d_fltdbl(float **table, int i, int j, double dx, double dy)
{
  double result;

  result =
    (1 - dx) * (1 - dy) * ((float) table[i][j]) +
    (1 - dx) * dy * ((float) table[i][j + 1]) +
    dx * (1 - dy) * ((float) table[i + 1][j]) +
    dx * dy * ((float) table[i + 1][j + 1]);

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a tri-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_3d_mydbl(double ***table, int i, int j, int k, double dx, double dy, double dz)
{
  double result;

  result =
    (1 - dx) * (1 - dy) * (1 - dz) * table[i][j][k] +
    (1 - dx) * (1 - dy) * dz * table[i][j][k + 1] +
    (1 - dx) * dy * (1 - dz) * table[i][j + 1][k] +
    (1 - dx) * dy * dz * table[i][j + 1][k + 1] +
    dx * (1 - dy) * (1 - dz) * table[i + 1][j][k] +
    dx * (1 - dy) * dz * table[i + 1][j][k + 1] +
    dx * dy * (1 - dz) * table[i + 1][j + 1][k] + dx * dy * dz * table[i + 1][j + 1][k + 1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a tri-linear interpolation, but for a 4D table
 * (Note that the last index is kept fixed). This is used for the
 * equilibrium abundance tables.
 * ----------------------------------------------------------------------
 */

double interpol_3d_special(double ****table, int i, int j, int k, int l, double dx, double dy, double dz)
{
  double result;

  result =
    (1 - dx) * (1 - dy) * (1 - dz) * table[i][j][k][l] +
    (1 - dx) * (1 - dy) * dz * table[i][j][k+1][l] +
    (1 - dx) * dy * (1 - dz) * table[i][j+1][k][l] +
    (1 - dx) * dy * dz * table[i][j+1][k+1][l] +
    dx * (1 - dy) * (1 - dz) * table[i+1][j][k][l] +
    dx * (1 - dy) * dz * table[i+1][j][k+1][l] +
    dx * dy * (1 - dz) * table[i+1][j+1][k][l] +
    dx * dy * dz * table[i+1][j+1][k+1][l];

  return result;
}
  

/*
 * ----------------------------------------------------------------------
 * This routine performs a quadri-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_4d_mydbl(double ****table, int i, int j, int k,
		  int l, double dx, double dy, double dz, double dw)
{
  double result;
  
  result =
    (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * table[i][j][k][l] +
    (1 - dx) * (1 - dy) * (1 - dz) * dw * table[i][j][k][l+1] +
    (1 - dx) * (1 - dy) * dz * (1 - dw) * table[i][j][k+1][l] +
    (1 - dx) * (1 - dy) * dz * dw * table[i][j][k+1][l+1] +
    (1 - dx) * dy * (1 - dz) * (1 - dw) * table[i][j+1][k][l] +
    (1 - dx) * dy * (1 - dz) * dw * table[i][j+1][k][l+1] +
    (1 - dx) * dy * dz * (1 - dw) * table[i][j+1][k+1][l] +
    (1 - dx) * dy * dz * dw * table[i][j+1][k+1][l+1] +
    dx * (1 - dy) * (1 - dz) * (1 - dw) * table[i+1][j][k][l] +
    dx * (1 - dy) * (1 - dz) * dw * table[i+1][j][k][l+1] +
    dx * (1 - dy) * dz * (1 - dw) * table[i+1][j][k+1][l] +
    dx * (1 - dy) * dz * dw * table[i+1][j][k+1][l+1] +
    dx * dy * (1 - dz) * (1 - dw) * table[i+1][j+1][k][l] +
    dx * dy * (1 - dz) * dw * table[i+1][j+1][k][l+1] +
    dx * dy * dz * (1 - dw) * table[i+1][j+1][k+1][l] + 
    dx * dy * dz * dw * table[i+1][j+1][k+1][l+1];

  return result;
}

/*
 * ----------------------------------------------------------------------
 * This routine performs a quinti-linear interpolation
 * ----------------------------------------------------------------------
 */

double interpol_5d_mydbl(double *****table, int i, int j, int k,
		  int l, int m, double dx, double dy, double dz, double dw, double dv)
{
  double result;
  
  result =
    (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * (1 - dv) * table[i][j][k][l][m] +
    (1 - dx) * (1 - dy) * (1 - dz) * dw * (1 - dv) * table[i][j][k][l+1][m] +
    (1 - dx) * (1 - dy) * dz * (1 - dw) * (1 - dv) * table[i][j][k+1][l][m] +
    (1 - dx) * (1 - dy) * dz * dw * (1 - dv) * table[i][j][k+1][l+1][m] +
    (1 - dx) * dy * (1 - dz) * (1 - dw) * (1 - dv) * table[i][j+1][k][l][m] +
    (1 - dx) * dy * (1 - dz) * dw * (1 - dv) * table[i][j+1][k][l+1][m] +
    (1 - dx) * dy * dz * (1 - dw) * (1 - dv) * table[i][j+1][k+1][l][m] +
    (1 - dx) * dy * dz * dw * (1 - dv) * table[i][j+1][k+1][l+1][m] +
    dx * (1 - dy) * (1 - dz) * (1 - dw) * (1 - dv) * table[i+1][j][k][l][m] +
    dx * (1 - dy) * (1 - dz) * dw * (1 - dv) * table[i+1][j][k][l+1][m] +
    dx * (1 - dy) * dz * (1 - dw) * (1 - dv) * table[i+1][j][k+1][l][m] +
    dx * (1 - dy) * dz * dw * (1 - dv) * table[i+1][j][k+1][l+1][m] +
    dx * dy * (1 - dz) * (1 - dw) * (1 - dv) * table[i+1][j+1][k][l][m] +
    dx * dy * (1 - dz) * dw * (1 - dv) * table[i+1][j+1][k][l+1][m] +
    dx * dy * dz * (1 - dw) * (1 - dv) * table[i+1][j+1][k+1][l][m] + 
    dx * dy * dz * dw * (1 - dv) * table[i+1][j+1][k+1][l+1][m] + 
    (1 - dx) * (1 - dy) * (1 - dz) * (1 - dw) * dv * table[i][j][k][l][m+1] +
    (1 - dx) * (1 - dy) * (1 - dz) * dw * dv * table[i][j][k][l+1][m+1] +
    (1 - dx) * (1 - dy) * dz * (1 - dw) * dv * table[i][j][k+1][l][m+1] +
    (1 - dx) * (1 - dy) * dz * dw * dv * table[i][j][k+1][l+1][m+1] +
    (1 - dx) * dy * (1 - dz) * (1 - dw) * dv * table[i][j+1][k][l][m+1] +
    (1 - dx) * dy * (1 - dz) * dw * dv * table[i][j+1][k][l+1][m+1] +
    (1 - dx) * dy * dz * (1 - dw) * dv * table[i][j+1][k+1][l][m+1] +
    (1 - dx) * dy * dz * dw * dv * table[i][j+1][k+1][l+1][m+1] +
    dx * (1 - dy) * (1 - dz) * (1 - dw) * dv * table[i+1][j][k][l][m+1] +
    dx * (1 - dy) * (1 - dz) * dw * dv * table[i+1][j][k][l+1][m+1] +
    dx * (1 - dy) * dz * (1 - dw) * dv * table[i+1][j][k+1][l][m+1] +
    dx * (1 - dy) * dz * dw * dv * table[i+1][j][k+1][l+1][m+1] +
    dx * dy * (1 - dz) * (1 - dw) * dv * table[i+1][j+1][k][l][m+1] +
    dx * dy * (1 - dz) * dw * dv * table[i+1][j+1][k][l+1][m+1] +
    dx * dy * dz * (1 - dw) * dv * table[i+1][j+1][k+1][l][m+1] + 
    dx * dy * dz * dw * dv * table[i+1][j+1][k+1][l+1][m+1];

  return result;
}

