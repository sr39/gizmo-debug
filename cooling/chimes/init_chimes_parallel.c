#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <hdf5.h>
#include "../../allvars.h" 

#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif

void ReadPhotoIonTables_parallel(struct globalVariables *myGlobalVars, char *photoIonTablePath, struct NonEq_Ionization *myNonEqIon, int N_Elements_in_Bens_tables, double *dustG_arr, double *H2_dissocJ_arr, int current_spectrum); 
void GetPhotoIonTables_parallel(struct globalVariables *myGlobalVars, int N_Elements_in_Bens_tables, struct All_rate_variables_structure **this_all_rates, double *dustG_arr, double *H2_dissocJ_arr);
void initialise_bens_tables_parallel(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates, double *dustG_arr, double *H2_dissocJ_arr); 

/* The following routine reads a photoionisation table into the given NonEqIon structure. */
void ReadPhotoIonTables_parallel(struct globalVariables *myGlobalVars, char *photoIonTablePath, struct NonEq_Ionization *myNonEqIon, int N_Elements_in_Bens_tables, double *dustG_arr, double *H2_dissocJ_arr, int current_spectrum)
{
  hid_t file_id, dataset, dataspace_id, memspace_id;
  herr_t status;
  hsize_t dims[2];
  hsize_t dims2D[2], count2D[2], offset2D[2];
  hsize_t dims3D[3], count3D[3], offset3D[3];
  hsize_t dims4D[4], count4D[4], offset4D[4];
  int rank, i, j, k, l, m, ns, r, included_index, N_arrayCells;
  char set_name[500];
  double *sigmaphot, *E_thresh, *epsilon, *cosmicrays;
  double *H2CO_N, *COself_N, *CO_S, *x_ion, *n_ion_HI, *n_ion_HeI, *shieldColumn;
  double H2_dissocJ, dust_G_parameter;
  float *shieldFac_1D, *shieldFac_2D; 
  dims[1] = 1;
	
  // Each element only needs to be read in by 1 task. 
  if (ThisTask < chimesRateTables.NonEqIon->N_Elements) 
    file_id = H5Fopen(photoIonTablePath, H5F_ACC_RDONLY, H5P_DEFAULT); 

  if (current_spectrum == 0)
    {
      /* First time this routine is called, read in tables and 
       * parameters that are only needed once, not for each 
       * UV spectrum. */ 
      int N_Auger[N_Elements_in_Bens_tables];
      myNonEqIon->N_Auger = (int *) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(int));
  
      if (ThisTask == 0)
	{
	  sprintf(set_name, "/N_Auger");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_Auger);
	  status = H5Dclose(dataset);
  
	  sprintf(set_name, "/Column_density_dimension");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, myNonEqIon->shieldingColumnDimensions);
	  status = H5Dclose(dataset);

	  /* Tables for the shielding of CO */
	  sprintf(set_name, "/shielding_dimensions");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, myNonEqIon->shielding_dimensions);
	  status = H5Dclose(dataset);

	  /* Tables for secondary CR ionisation of HI and HeI */
	  sprintf(set_name, "/secondary_ionisation/bin_dimensions");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, myNonEqIon->secondary_ionisation_dims);
	  status = H5Dclose(dataset);
	}

      MPI_Bcast(N_Auger, N_Elements_in_Bens_tables, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(myNonEqIon->shieldingColumnDimensions, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(myNonEqIon->shielding_dimensions, 2, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(myNonEqIon->secondary_ionisation_dims, 1, MPI_INT, 0, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);

      included_index = 0;
      for (i = 0; i < N_Elements_in_Bens_tables; i++)
	{
	  if (i < 2 || myGlobalVars->element_included[max(i - 2, 0)] == 1)
	    {
	      myNonEqIon->N_Auger[included_index] = N_Auger[i];
	      included_index += 1;
	    }
	}

      myNonEqIon->shieldingColumnDensities = (double *) malloc(myNonEqIon->shieldingColumnDimensions[0] * sizeof(double));
      shieldColumn = (double *) malloc(myNonEqIon->shieldingColumnDimensions[0] * sizeof(double));
      
      myNonEqIon->COself_shielding_N = (double *) malloc(myNonEqIon->shielding_dimensions[0] * sizeof(double));
      COself_N = (double *) malloc(myNonEqIon->shielding_dimensions[0] * sizeof(double));

      myNonEqIon->H2CO_shielding_N = (double *) malloc(myNonEqIon->shielding_dimensions[1] * sizeof(double));
      H2CO_N = (double *) malloc(myNonEqIon->shielding_dimensions[1] * sizeof(double));

      myNonEqIon->CO_shielding_S = (double **) malloc(myNonEqIon->shielding_dimensions[0] * sizeof(double *));
      for (i = 0; i < myNonEqIon->shielding_dimensions[0]; i++)
	myNonEqIon->CO_shielding_S[i] = (double *) malloc(myNonEqIon->shielding_dimensions[1] * sizeof(double));
      
      myNonEqIon->x_ion_fraction = (double *) malloc(myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      x_ion = (double *) malloc(myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      
      myNonEqIon->n_ion_HI = (double *) malloc(myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      n_ion_HI = (double *) malloc(myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      
      myNonEqIon->n_ion_HeI = (double *) malloc(myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));
      n_ion_HeI = (double *) malloc(myNonEqIon->secondary_ionisation_dims[0] * sizeof(double));

      if (ThisTask == 0) 
	{
	  sprintf(set_name, "/Column_density_bins");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, shieldColumn);
	  status = H5Dclose(dataset);
      
	  sprintf(set_name, "/COself_shielding_N");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, COself_N);
	  status = H5Dclose(dataset);
    
	  sprintf(set_name, "/H2CO_shielding_N");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, H2CO_N);
	  status = H5Dclose(dataset);

	  CO_S = (double *) malloc(myNonEqIon->shielding_dimensions[0] * sizeof(double));

	  sprintf(set_name, "/CO_shielding_S");
	  dataset = H5Dopen(file_id, set_name);
	  dataspace_id = H5Dget_space(dataset);
	  status = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
	  dims[0] = myNonEqIon->shielding_dimensions[0];
	  rank = 1;
	  memspace_id = H5Screate_simple(rank, dims, NULL);
      
	  for (j = 0; j < myNonEqIon->shielding_dimensions[1]; j++)
	    {
	      offset2D[0] = 0;
	      offset2D[1] = j;
	      count2D[0] = myNonEqIon->shielding_dimensions[0];
	      count2D[1] = 1;
	      status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);
	      
	      status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, CO_S);
	      
	      for (i = 0; i < myNonEqIon->shielding_dimensions[0]; i++)
		myNonEqIon->CO_shielding_S[i][j] = CO_S[i];
	    }
	  H5Sclose(memspace_id);
	  H5Sclose(dataspace_id);
	  status = H5Dclose(dataset);
	  free(CO_S);
      
	  sprintf(set_name, "/secondary_ionisation/x_ion");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x_ion);
	  status = H5Dclose(dataset);
      
	  sprintf(set_name, "/secondary_ionisation/n_ion_HI");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_ion_HI);
	  status = H5Dclose(dataset);

	  sprintf(set_name, "/secondary_ionisation/n_ion_HeI");
	  dataset = H5Dopen(file_id, set_name);
	  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_ion_HeI);
	  status = H5Dclose(dataset);
	}

      MPI_Bcast(shieldColumn, myNonEqIon->shieldingColumnDimensions[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(COself_N, myNonEqIon->shielding_dimensions[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(H2CO_N, myNonEqIon->shielding_dimensions[1], MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(x_ion, myNonEqIon->secondary_ionisation_dims[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(n_ion_HI, myNonEqIon->secondary_ionisation_dims[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(n_ion_HeI, myNonEqIon->secondary_ionisation_dims[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);

      /* 2D array - broadcast into 1D buffer first. */
      CO_S = (double *) malloc(myNonEqIon->shielding_dimensions[0] * myNonEqIon->shielding_dimensions[1] * sizeof(double));
      if (ThisTask == 0)
        {
          for (i = 0; i < myNonEqIon->shielding_dimensions[0]; i++)
	    {
	      for (j = 0; j < myNonEqIon->shielding_dimensions[1]; j++)
		CO_S[i * myNonEqIon->shielding_dimensions[1] + j] = myNonEqIon->CO_shielding_S[i][j];
	    }
        }
      MPI_Bcast(CO_S, myNonEqIon->shielding_dimensions[0] * myNonEqIon->shielding_dimensions[1], MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
	  
      if (ThisTask != 0)
        {
          for (i = 0; i < myNonEqIon->shielding_dimensions[0]; i++)
	    {
	      for (j = 0; j < myNonEqIon->shielding_dimensions[1]; j++)
		myNonEqIon->CO_shielding_S[i][j] = CO_S[i * myNonEqIon->shielding_dimensions[1] + j];
	    }
        }

      for (i = 0; i < myNonEqIon->shieldingColumnDimensions[0]; i++)
	myNonEqIon->shieldingColumnDensities[i] = log10(shieldColumn[i]);
  
      for (i = 0; i < myNonEqIon->shielding_dimensions[0]; i++)
	myNonEqIon->COself_shielding_N[i] = COself_N[i];
  
      for (i = 0; i < myNonEqIon->shielding_dimensions[1]; i++)
	myNonEqIon->H2CO_shielding_N[i] = H2CO_N[i];
      
      for (i = 0; i < myNonEqIon->secondary_ionisation_dims[0]; i++)
	myNonEqIon->x_ion_fraction[i] = log10(x_ion[i]);
      
      for (i = 0; i < myNonEqIon->secondary_ionisation_dims[0]; i++)
	myNonEqIon->n_ion_HI[i] = log10(n_ion_HI[i]);
      
      for (i = 0; i < myNonEqIon->secondary_ionisation_dims[0]; i++)
	myNonEqIon->n_ion_HeI[i] = log10(n_ion_HeI[i]);

      free(shieldColumn); 
      free(H2CO_N);
      free(COself_N);
      free(CO_S);
      free(x_ion);
      free(n_ion_HI);
      free(n_ion_HeI);
      
      /* Also allocate memory for tables first time this routine is called. */ 
      ns = 0;
      for (l = 0; l < N_Elements_in_Bens_tables; l++)
	{  
	  if (l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
	    {
	      if (chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
		N_arrayCells = 3;   /* Hydrogen photoionisation also includes H- & H2*/
	      else
		N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

	      myNonEqIon->NonEqRates[ns].sigmaphot = (double ***) malloc(myGlobalVars->N_spectra * sizeof(double **));
	      for (j = 0; j < myGlobalVars->N_spectra; j++)
		{
		  myNonEqIon->NonEqRates[ns].sigmaphot[j] = (double **) malloc(N_arrayCells * sizeof(double *));  
		  for (k = 0; k < N_arrayCells; k++)
		    myNonEqIon->NonEqRates[ns].sigmaphot[j][k] = (double *) malloc(myNonEqIon->N_Auger[ns] * sizeof(double)); 
		}
	      
	      myNonEqIon->NonEqRates[ns].shieldFactor1D = (float ****) malloc(myGlobalVars->N_spectra * sizeof(float ***));
	      for (j = 0; j < myGlobalVars->N_spectra; j++)
		{
		  myNonEqIon->NonEqRates[ns].shieldFactor1D[j] = (float ***) malloc(3 * sizeof(float **));
		  for (k = 0; k < 3; k++)
		    {
		      myNonEqIon->NonEqRates[ns].shieldFactor1D[j][k] = (float **) malloc(N_arrayCells * sizeof(float *));
		      for (m = 0; m < N_arrayCells; m++)
			myNonEqIon->NonEqRates[ns].shieldFactor1D[j][k][m] = (float *) malloc(myNonEqIon->shieldingColumnDimensions[0] * sizeof(float)); 
		    }
		}

	      myNonEqIon->NonEqRates[ns].shieldFactor2D = (float *****) malloc(myGlobalVars->N_spectra * sizeof(float ****));
	      for (m = 0; m < myGlobalVars->N_spectra; m++)
		{
		  myNonEqIon->NonEqRates[ns].shieldFactor2D[m] = (float ****) malloc(6 * sizeof(float ***));
		  for (j = 0; j < 6; j++)
		    {
		      myNonEqIon->NonEqRates[ns].shieldFactor2D[m][j] = (float ***) malloc(N_arrayCells * sizeof(float **));
		      for (k = 0; k < N_arrayCells; k++)
			{
			  myNonEqIon->NonEqRates[ns].shieldFactor2D[m][j][k] = (float **) malloc(myNonEqIon->shieldingColumnDimensions[0] * sizeof(float *));
			  for (r = 0; r < myNonEqIon->shieldingColumnDimensions[0]; r++) 
			    myNonEqIon->NonEqRates[ns].shieldFactor2D[m][j][k][r] = (float *) malloc(myNonEqIon->shieldingColumnDimensions[0] * sizeof(float)); 
			}
		    }
		}  

	      myNonEqIon->NonEqRates[ns].epsilon = (double **) malloc(myGlobalVars->N_spectra * sizeof(double *));
	      for (j = 0; j < myGlobalVars->N_spectra; j++)
		myNonEqIon->NonEqRates[ns].epsilon[j] = (double *) malloc(N_arrayCells * sizeof(double)); 
	      
	      /* Read in E_thresh and cosmicrays from table. */
	      E_thresh = (double *) malloc(N_arrayCells * sizeof(double));
	      cosmicrays = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[ns] - 1) * sizeof(double));  /* For hydrogen, H- & H2 cr rates are stored elsewhere */

	      if (ThisTask == ns)
		{ 
		  sprintf(set_name, "/%s/E_thresh", chimesRateTables.NonEqIon->ElementName[ns]);
		  dataset = H5Dopen(file_id, set_name);
		  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, E_thresh);
		  status = H5Dclose(dataset);
	      
		  sprintf(set_name, "/%s/cosmicRays", chimesRateTables.NonEqIon->ElementName[ns]);
		  dataset = H5Dopen(file_id, set_name);
		  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, cosmicrays);
		  status = H5Dclose(dataset);
		}

	      MPI_Bcast(E_thresh, N_arrayCells, MPI_DOUBLE, ns, MPI_COMM_WORLD);
	      MPI_Bcast(cosmicrays, N_arrayCells, MPI_DOUBLE, ns, MPI_COMM_WORLD);
	      MPI_Barrier(MPI_COMM_WORLD);

	      for (i = 0; i < N_arrayCells; i++)
		{
		  myNonEqIon->NonEqRates[ns].E_thresh[i] = E_thresh[i];
		  myNonEqIon->NonEqRates[ns].cosmicRays[i] = cosmicrays[i];
		}
	      free(E_thresh);
	      free(cosmicrays);
	      
	      ns += 1;
	    }
	}
    }
  
  // The following need to be read in for each UV spectrum. 
  ns = 0;
  for (l = 0; l < N_Elements_in_Bens_tables; l++)
    {  
      if (l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
	{
	  if (ThisTask == ns) 
	    {
	      if (chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
		N_arrayCells = 3;   /* Hydrogen photoionisation also includes H- & H2*/
	      else
		N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

	      sigmaphot = (double *) malloc(N_arrayCells * sizeof(double));
	      epsilon = (double *) malloc(N_arrayCells * sizeof(double));
	      shieldFac_1D = (float *) malloc(3 * sizeof(float));
	      shieldFac_2D = (float *) malloc(6 * sizeof(float));
	  
	      sprintf(set_name, "/%s/sigmaPhot", chimesRateTables.NonEqIon->ElementName[ns]);
	      dataset = H5Dopen(file_id, set_name);  
	  
	      if (ns == 0 || ns == 1)
		{
		  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sigmaphot);
		  for (i = 0; i < N_arrayCells; i++)
		    myNonEqIon->NonEqRates[ns].sigmaphot[current_spectrum][i][0] = sigmaphot[i];
		  status = H5Dclose(dataset);
		}
	      else
		{		
		  dataspace_id = H5Dget_space(dataset);
		  status = H5Sget_simple_extent_dims(dataspace_id, dims2D, NULL);
		  dims[0] = N_arrayCells;
		  rank = 1;
		  memspace_id = H5Screate_simple(rank, dims, NULL);
		  
		  for (j = 0; j < myNonEqIon->N_Auger[ns]; j++)
		    {
		      offset2D[0] = 0;
		      offset2D[1] = j;
		      count2D[0] = N_arrayCells;
		      count2D[1] = 1;
		      status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset2D, NULL, count2D, NULL);
		      
		      status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, sigmaphot);
		      
		      for (i = 0; i < N_arrayCells; i++)
			myNonEqIon->NonEqRates[ns].sigmaphot[current_spectrum][i][j] = sigmaphot[i];
		    }
		  H5Sclose(memspace_id);
		  H5Sclose(dataspace_id);
		  status = H5Dclose(dataset);
		}
	      
	      sprintf(set_name, "/%s/epsilonPhot", chimesRateTables.NonEqIon->ElementName[ns]);
	      dataset = H5Dopen(file_id, set_name);
	      status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, epsilon);
	      status = H5Dclose(dataset);    
	      
	      sprintf(set_name, "/%s/shieldFactor_1D", chimesRateTables.NonEqIon->ElementName[ns]);
	      dataset = H5Dopen(file_id, set_name);
	      dataspace_id = H5Dget_space(dataset);
	      status = H5Sget_simple_extent_dims(dataspace_id, dims3D, NULL);
	      dims[0] = 3;
	      rank = 1;
	      memspace_id = H5Screate_simple(rank, dims, NULL);
	      
	      for (i = 0; i < N_arrayCells; i++)
		{
		  for (j = 0; j < myNonEqIon->shieldingColumnDimensions[0]; j++)
		    {
		      offset3D[0] = 0;
		      offset3D[1] = i;
		      offset3D[2] = j;
		      count3D[0] = 3;
		      count3D[1] = 1;
		      count3D[2] = 1;
		      status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset3D, NULL, count3D, NULL);
		      status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, shieldFac_1D);
		      
		      for (k = 0; k < 3; k++)
			myNonEqIon->NonEqRates[ns].shieldFactor1D[current_spectrum][k][i][j] = (float) max(shieldFac_1D[k], -300.0);
		    }
		}
	      H5Sclose(memspace_id);
	      H5Sclose(dataspace_id);
	      H5Dclose(dataset);
	      
	      sprintf(set_name, "/%s/shieldFactor_2D", chimesRateTables.NonEqIon->ElementName[ns]);
	      dataset = H5Dopen(file_id, set_name);
	      dataspace_id = H5Dget_space(dataset);
	      status = H5Sget_simple_extent_dims(dataspace_id, dims4D, NULL);
	      dims[0] = 6;
	      rank = 1;
	      memspace_id = H5Screate_simple(rank, dims, NULL);
	      
	      for (m = 0; m < N_arrayCells; m++)
		{
		  for (i = 0; i < myNonEqIon->shieldingColumnDimensions[0]; i++)
		    {
		      for (j = 0; j < myNonEqIon->shieldingColumnDimensions[0]; j++)
			{
			  offset4D[0] = 0;
			  offset4D[1] = m;
			  offset4D[2] = i;
			  offset4D[3] = j;
			  count4D[0] = 6;
			  count4D[1] = 1;
			  count4D[2] = 1;
			  count4D[3] = 1;
			  status = H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset4D, NULL, count4D, NULL);
			  status = H5Dread(dataset, H5T_NATIVE_FLOAT, memspace_id, dataspace_id, H5P_DEFAULT, shieldFac_2D);
			  
			  for (k = 0; k < 6; k++)
			    myNonEqIon->NonEqRates[ns].shieldFactor2D[current_spectrum][k][m][i][j] = (float) max(shieldFac_2D[k], -300.0);
			}
		    }
		}
	      H5Sclose(memspace_id);
	      H5Sclose(dataspace_id);
	      H5Dclose(dataset);
	      
	      for (i = 0; i < N_arrayCells; i++)
		{
		  myNonEqIon->NonEqRates[ns].epsilon[current_spectrum][i] = epsilon[i];
		}
	      
	      free(sigmaphot);
	      free(epsilon);
	      free(shieldFac_1D);
	      free(shieldFac_2D);

	      ns += 1;
	    }
	  else 
	    ns += 1; 
	}
    }

  /* Now broadcast each element to all tasks. */
  ns = 0;
  for (l = 0; l < N_Elements_in_Bens_tables; l++)
    {  
      if (l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
	{
	  if (chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
	    N_arrayCells = 3;   /* Hydrogen photoionisation also includes H- & H2*/
	  else
	    N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;
	  
	  /* For multi-D arrays, use a 1D buffer. */
	  epsilon = (double *) malloc(N_arrayCells * sizeof(double)); 
	  if (ThisTask == ns) 
	    {
	      for (i = 0; i < N_arrayCells; i++)
		epsilon[i] = myNonEqIon->NonEqRates[ns].epsilon[current_spectrum][i]; 
	    }
	  
	  MPI_Bcast(epsilon, N_arrayCells, MPI_DOUBLE, ns, MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD); 

	  if (ThisTask != ns) 
	    {
	      for (i = 0; i < N_arrayCells; i++)
		myNonEqIon->NonEqRates[ns].epsilon[current_spectrum][i] = epsilon[i]; 
	    }
	      
	  if (ns == 0 || ns == 1)
	    sigmaphot = (double *) malloc(N_arrayCells * sizeof(double));
	  else
	    sigmaphot = (double *) malloc(N_arrayCells * myNonEqIon->N_Auger[ns] * sizeof(double));
	  
	  shieldFac_1D = (float *) malloc(3 * N_arrayCells * myNonEqIon->shieldingColumnDimensions[0] * sizeof(float));
	  shieldFac_2D = (float *) malloc(6 * N_arrayCells * myNonEqIon->shieldingColumnDimensions[0] * myNonEqIon->shieldingColumnDimensions[0] * sizeof(float));
		  
	  if (ThisTask == ns)
	    {
	      if (ns == 0 || ns == 1)
		{
		  for (i = 0; i < N_arrayCells; i++)
		    sigmaphot[i] = myNonEqIon->NonEqRates[ns].sigmaphot[current_spectrum][i][0];
		}
	      else
		{
		  for (i = 0; i < N_arrayCells; i++)
		    {
		      for (j = 0; j < myNonEqIon->N_Auger[ns]; j++)
			sigmaphot[(i * myNonEqIon->N_Auger[ns]) + j] = myNonEqIon->NonEqRates[ns].sigmaphot[current_spectrum][i][j];
		    }
		}
	      
	      for (i = 0; i < 3; i++)
		{
		  for(j = 0; j < N_arrayCells; j++)
		    {
		      for (k = 0; k < myNonEqIon->shieldingColumnDimensions[0]; k++)
			shieldFac_1D[(i * myNonEqIon->shieldingColumnDimensions[0] * N_arrayCells) + (j * myNonEqIon->shieldingColumnDimensions[0]) + k] = myNonEqIon->NonEqRates[ns].shieldFactor1D[current_spectrum][i][j][k];
		    }
		}
	      
	      for (i = 0; i < 6; i++)
		{
		  for(j = 0; j < N_arrayCells; j++)
		    {
		      for (k = 0; k < myNonEqIon->shieldingColumnDimensions[0]; k++)
			{
			  for (m = 0; m < myNonEqIon->shieldingColumnDimensions[0]; m++)
			    shieldFac_2D[(i * myNonEqIon->shieldingColumnDimensions[0] * myNonEqIon->shieldingColumnDimensions[0] * N_arrayCells) + (j * myNonEqIon->shieldingColumnDimensions[0] * myNonEqIon->shieldingColumnDimensions[0]) + (k * myNonEqIon->shieldingColumnDimensions[0]) + m] = myNonEqIon->NonEqRates[ns].shieldFactor2D[current_spectrum][i][j][k][m];
			}
		    }
		}
	    }
	  
	  if (ns == 0 || ns == 1)
	    MPI_Bcast(sigmaphot, N_arrayCells, MPI_DOUBLE, ns, MPI_COMM_WORLD);
	  else
	    MPI_Bcast(sigmaphot, N_arrayCells * myNonEqIon->N_Auger[ns], MPI_DOUBLE, ns, MPI_COMM_WORLD);
	  
	  MPI_Bcast(shieldFac_1D, 3 * N_arrayCells * myNonEqIon->shieldingColumnDimensions[0], MPI_FLOAT, ns, MPI_COMM_WORLD);
	  MPI_Bcast(shieldFac_2D, 6 * N_arrayCells * myNonEqIon->shieldingColumnDimensions[0] * myNonEqIon->shieldingColumnDimensions[0], MPI_FLOAT, ns, MPI_COMM_WORLD);
	  
	  MPI_Barrier(MPI_COMM_WORLD);
	  
	  if (ThisTask != ns)
	    {
	      if (ns == 0 || ns == 1)
		{
		  for (i = 0; i < N_arrayCells; i++)
		    myNonEqIon->NonEqRates[ns].sigmaphot[current_spectrum][i][0] = sigmaphot[i];
		}
	      else
		{
		  for (i = 0; i < N_arrayCells; i++)
		    {
		      for (j = 0; j < myNonEqIon->N_Auger[ns]; j++)
			myNonEqIon->NonEqRates[ns].sigmaphot[current_spectrum][i][j] = sigmaphot[(i * myNonEqIon->N_Auger[ns]) + j];
		    }
		}
	      
	      for (i = 0; i < 3; i++)
		{
		  for(j = 0; j < N_arrayCells; j++)
		    {
		      for (k = 0; k < myNonEqIon->shieldingColumnDimensions[0]; k++)
			myNonEqIon->NonEqRates[ns].shieldFactor1D[current_spectrum][i][j][k] = shieldFac_1D[(i * myNonEqIon->shieldingColumnDimensions[0] * N_arrayCells) + (j * myNonEqIon->shieldingColumnDimensions[0]) + k];
		    }
		}
	      
	      for (i = 0; i < 6; i++)
		{
		  for(j = 0; j < N_arrayCells; j++)
		    {
		      for (k = 0; k < myNonEqIon->shieldingColumnDimensions[0]; k++)
			{
			  for (m = 0; m < myNonEqIon->shieldingColumnDimensions[0]; m++)
			    myNonEqIon->NonEqRates[ns].shieldFactor2D[current_spectrum][i][j][k][m] = shieldFac_2D[(i * myNonEqIon->shieldingColumnDimensions[0] * myNonEqIon->shieldingColumnDimensions[0] * N_arrayCells) + (j * myNonEqIon->shieldingColumnDimensions[0] * myNonEqIon->shieldingColumnDimensions[0]) + (k * myNonEqIon->shieldingColumnDimensions[0]) + m];
			}
		    }
		}
	    } 
	  
	  MPI_Barrier(MPI_COMM_WORLD);

	  free(epsilon); 
	  free(sigmaphot);
	  free(shieldFac_1D);
	  free(shieldFac_2D);
	  
	  ns += 1;
	}
    }
  
  if (ThisTask == 0) 
    {
      /* Read in H2_dissocJ and dust_G_parameter for each UV spectrum. */ 
      sprintf(set_name, "/H2_dissocJ");
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &H2_dissocJ);
      status = H5Dclose(dataset);    
      
      sprintf(set_name, "/dust_G_parameter");
      dataset = H5Dopen(file_id, set_name);
      status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dust_G_parameter);
      status = H5Dclose(dataset);    
    }
  
  MPI_Bcast(&H2_dissocJ, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&dust_G_parameter, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
      
  H2_dissocJ_arr[current_spectrum] = H2_dissocJ;
  dustG_arr[current_spectrum] = dust_G_parameter; 
  
  if (ThisTask < chimesRateTables.NonEqIon->N_Elements) 
    status = H5Fclose(file_id);
}

void GetPhotoIonTables_parallel(struct globalVariables *myGlobalVars, int N_Elements_in_Bens_tables, struct All_rate_variables_structure **this_all_rates, double *dustG_arr, double *H2_dissocJ_arr)
{
  int j, l, ns, N_arrayCells;

  /* Read in list of paths to the PhotoIon tables, 
   * and read in each table in turn. */
  FILE *fdin; 
  char buffer[500]; 
  char table_path[500]; 
  
  if ((fdin = fopen(myGlobalVars->PhotoIonTablePath, "r"))) 
    {
      int current_spectrum = 0; 
      while (fgets(buffer, 500, fdin)) 
	{
	  if (current_spectrum > myGlobalVars->N_spectra - 1) 
	    {
	      printf("ERROR: too many UV spectra specified in PhotoIon table list. Aborting. \n"); 
	      exit(-1); 
	    }
	  sscanf(buffer, "%s", table_path);
	  printf("Reading PhotoIon table: %s \n", table_path); 
	  ReadPhotoIonTables_parallel(myGlobalVars, table_path, chimesRateTables.NonEqIon, N_Elements_in_Bens_tables, dustG_arr, H2_dissocJ_arr, current_spectrum);
	  current_spectrum += 1;
	}
      if (current_spectrum < myGlobalVars->N_spectra)
	{
	  printf("ERROR: too few UV spectra specified in PhotoIon table list. Aborting. \n"); 
	  exit(-1); 
	}
    }
  else
    {
      printf("ERROR: PhotoIon table list %s not found. \n", myGlobalVars->PhotoIonTablePath);
      exit(-1);
    }
  

  /* Allocate memory for the photoionisation rates
   * in the AllRates structure. */
  ns = 0;
  for (l = 0; l < N_Elements_in_Bens_tables; l++)
    {  
      if (l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
	{
	  if (chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
	    N_arrayCells = 3;   /* Hydrogen photoionisation also includes H- & H2*/
	  else
	    N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

	  (*this_all_rates)->BensRates[ns].PhotoIon = (double **) malloc(N_arrayCells * sizeof(double *));
	  for (j = 0; j < N_arrayCells; j++)
	    (*this_all_rates)->BensRates[ns].PhotoIon[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Auger[ns] * sizeof(double));
	  
	  ns += 1;
	}
    }
}

void initialise_bens_tables_parallel(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates, double *dustG_arr, double *H2_dissocJ_arr)
{
  
  int i,j,k, l, m; 
  hid_t file_id, dataset, datatype;
  hsize_t el_name_length = EL_NAME_LENGTH;
  hsize_t file_name_length = NONEQION_NAME_LENGTH; 
  int N_Elements_in_Bens_tables; 
  int included_index;
  char fname[256];

  chimesRateTables.NonEqIon = (struct NonEq_Ionization *) malloc(sizeof(struct NonEq_Ionization));

  if (ThisTask == 0) 
    {
      /* Task0 reads in header info. */
      sprintf(fname, "%sIonRates_Head.HM01.hdf5", myGlobalVars->BenTablesPath);		/* Note that we ignore everything that depends on the UV spectrum, so we can pick any of them */

      file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

      if(file_id < 0)
	printf("[read_noneq_tables()]: unable to open file %s\n", fname);

      dataset = H5Dopen(file_id, "N_Elements");
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_Elements_in_Bens_tables);
      H5Dclose(dataset);
      
      dataset = H5Dopen(file_id, "N_Temperatures");
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.NonEqIon->N_Temperatures);
      H5Dclose(dataset);

      dataset = H5Dopen(file_id, "/nei_headers/BinSizes");
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.nei_cooling_table_dimensions[0]);
      H5Dclose(dataset);

      dataset = H5Dopen(file_id, "/chianti_headers/BinSizes");
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.chianti_cooling_table_dimensions[0]);
      H5Dclose(dataset);
    }

  /* Broadcast header info to all tasks. */
  MPI_Bcast(&N_Elements_in_Bens_tables, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&chimesRateTables.NonEqIon->N_Temperatures, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&chimesRateTables.nei_cooling_table_dimensions[0], 4, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&chimesRateTables.chianti_cooling_table_dimensions[0], 2, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  char element_names[N_Elements_in_Bens_tables][EL_NAME_LENGTH];
  char file_names[N_Elements_in_Bens_tables][NONEQION_NAME_LENGTH];

  if (ThisTask == 0) 
    {
      /* Read in the element names and filenames to temporary arrays. */
      datatype = H5Tcopy(H5T_C_S1);
      H5Tset_size(datatype, el_name_length);
      dataset = H5Dopen(file_id, "Element_Names");
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names);
      H5Dclose(dataset);
      H5Tclose(datatype);
      
      datatype = H5Tcopy(H5T_C_S1);
      H5Tset_size(datatype, file_name_length);
      dataset = H5Dopen(file_id, "File_Names");
      H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, file_names);
      H5Dclose(dataset);
      H5Tclose(datatype);
    }

  /* Broadcast element names and file names to all tasks. */
  MPI_Bcast(element_names, N_Elements_in_Bens_tables * el_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(file_names, N_Elements_in_Bens_tables * file_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
      
  /* Now we need to determine how many elements are actually
   * included in our network. */
  chimesRateTables.NonEqIon->N_Elements = 2;  /* H & He */
  for (i = 0; i < 9; i++)
    if (myGlobalVars->element_included[i] == 1)
      chimesRateTables.NonEqIon->N_Elements += 1;
  
  /* Finally, set the element names and file names
   * in the arrays for bens tables. */
  chimesRateTables.NonEqIon->ElementName = (char **) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(char *));
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    chimesRateTables.NonEqIon->ElementName[i] = (char *) malloc(EL_NAME_LENGTH * sizeof(char));

  included_index = 0;
  for(i = 0; i < N_Elements_in_Bens_tables; i++)
    {
      if (i < 2 || myGlobalVars->element_included[max(i - 2, 0)] == 1)
	{
	  strcpy(chimesRateTables.NonEqIon->ElementName[included_index], element_names[i]);
	  included_index += 1;
	}
    }

  chimesRateTables.NonEqIon->FileName = (char **) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(char *));
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    chimesRateTables.NonEqIon->FileName[i] = (char *) malloc(NONEQION_NAME_LENGTH * sizeof(char));

  included_index = 0;
  for(i = 0; i < N_Elements_in_Bens_tables; i++)
    {
      if (i < 2 || myGlobalVars->element_included[max(i - 2, 0)] == 1)
	{
	  strcpy(chimesRateTables.NonEqIon->FileName[included_index], file_names[i]);
	  included_index += 1;
	}
    }

  /* Allocate memory to remaining header arrays. */
  int atomic_number[N_Elements_in_Bens_tables];
  int N_Ions[N_Elements_in_Bens_tables];
  double atomic_weights[N_Elements_in_Bens_tables];

  chimesRateTables.NonEqIon->AtomicNumber = (int *) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(int));
  chimesRateTables.NonEqIon->N_Ions = (int *) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(int));
  chimesRateTables.NonEqIon->AtomicWeights = (double *) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(double));
  chimesRateTables.NonEqIon->Temperatures = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
  chimesRateTables.nei_cooling_temperature = (double *) malloc(chimesRateTables.nei_cooling_table_dimensions[0] * sizeof(double));
  chimesRateTables.nei_cooling_HIAbundance = (double *) malloc(chimesRateTables.nei_cooling_table_dimensions[1] * sizeof(double));
  chimesRateTables.nei_cooling_ElectronAbundance = (double *) malloc(chimesRateTables.nei_cooling_table_dimensions[2] * sizeof(double));
  chimesRateTables.nei_cooling_HIIAbundance = (double *) malloc(chimesRateTables.nei_cooling_table_dimensions[3] * sizeof(double));
  chimesRateTables.chianti_cooling_temperature = (double *) malloc(chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double));
  chimesRateTables.chianti_cooling_ElectronDensity = (double *) malloc(chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));
  chimesRateTables.NonEqIon->NonEqRates = (struct NonEq_Rates *) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(struct NonEq_Rates));

  if (ThisTask == 0)
    {
      /* Task0 reads in remaining header tables. */

      dataset = H5Dopen(file_id, "Temperatures");
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &chimesRateTables.NonEqIon->Temperatures[0]);
      H5Dclose(dataset);

      for(i = 0; i < chimesRateTables.NonEqIon->N_Temperatures; i++)
	chimesRateTables.NonEqIon->Temperatures[i] = log10(chimesRateTables.NonEqIon->Temperatures[i]);

      dataset = H5Dopen(file_id, "AtomicNumbers");
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomic_number);
      H5Dclose(dataset);

      dataset = H5Dopen(file_id, "N_Ions");
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, N_Ions);
      H5Dclose(dataset);

      dataset = H5Dopen(file_id, "AtomicWeights");
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, atomic_weights);
      H5Dclose(dataset);
  
      dataset = H5Dopen(file_id, "/nei_headers/Temperatures");
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.nei_cooling_temperature);
      H5Dclose(dataset);
      
      dataset = H5Dopen(file_id, "/nei_headers/HIDensities");
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.nei_cooling_HIAbundance);
      H5Dclose(dataset);
      
      dataset = H5Dopen(file_id, "/nei_headers/ElectronDensities");
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.nei_cooling_ElectronAbundance);
      H5Dclose(dataset);
  
      dataset = H5Dopen(file_id, "/nei_headers/HIIDensities");
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.nei_cooling_HIIAbundance);
      H5Dclose(dataset);
  
      dataset = H5Dopen(file_id, "/chianti_headers/Temperatures");
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.chianti_cooling_temperature);
      H5Dclose(dataset);
  
      dataset = H5Dopen(file_id, "/chianti_headers/ElectronDensity");
      H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, chimesRateTables.chianti_cooling_ElectronDensity);
      H5Dclose(dataset);

      H5Fclose(file_id);
    }

  /* Broadcast these arrays to all tasks. */
  MPI_Bcast(&chimesRateTables.NonEqIon->Temperatures[0], chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(atomic_number, N_Elements_in_Bens_tables, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(N_Ions, N_Elements_in_Bens_tables, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(atomic_weights, N_Elements_in_Bens_tables, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(chimesRateTables.nei_cooling_temperature, chimesRateTables.nei_cooling_table_dimensions[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(chimesRateTables.nei_cooling_HIAbundance, chimesRateTables.nei_cooling_table_dimensions[1], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(chimesRateTables.nei_cooling_ElectronAbundance, chimesRateTables.nei_cooling_table_dimensions[2], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(chimesRateTables.nei_cooling_HIIAbundance, chimesRateTables.nei_cooling_table_dimensions[3], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(chimesRateTables.chianti_cooling_temperature, chimesRateTables.chianti_cooling_table_dimensions[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(chimesRateTables.chianti_cooling_ElectronDensity, chimesRateTables.chianti_cooling_table_dimensions[1], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  
  included_index = 0;
  for(i = 0; i < N_Elements_in_Bens_tables; i++)
    {
      if (i < 2 || myGlobalVars->element_included[max(i - 2, 0)] == 1)
	{
	  chimesRateTables.NonEqIon->AtomicNumber[included_index] = atomic_number[i];
	  chimesRateTables.NonEqIon->N_Ions[included_index] = N_Ions[i];
	  chimesRateTables.NonEqIon->AtomicWeights[included_index] = atomic_weights[i];
	  included_index += 1;
	}
    }
  
  chimesRateTables.nei_cooling_CI = (double ****) malloc(chimesRateTables.nei_cooling_table_dimensions[0] * sizeof(double ***));
  for (i = 0; i < chimesRateTables.nei_cooling_table_dimensions[0]; i++)
    {
      chimesRateTables.nei_cooling_CI[i] = (double ***) malloc(chimesRateTables.nei_cooling_table_dimensions[1] * sizeof(double **));
      for (j = 0; j < chimesRateTables.nei_cooling_table_dimensions[1]; j++)
	{
	  chimesRateTables.nei_cooling_CI[i][j] = (double **) malloc(chimesRateTables.nei_cooling_table_dimensions[2] * sizeof(double *));
	  for (k = 0; k < chimesRateTables.nei_cooling_table_dimensions[2]; k++)
	    {	
	      chimesRateTables.nei_cooling_CI[i][j][k] = (double *) malloc(chimesRateTables.nei_cooling_table_dimensions[3] * sizeof(double));
	    }
	}
    }

  chimesRateTables.nei_cooling_OI = (double ****) malloc(chimesRateTables.nei_cooling_table_dimensions[0] * sizeof(double ***));
  for (i = 0; i < chimesRateTables.nei_cooling_table_dimensions[0]; i++)
    {
      chimesRateTables.nei_cooling_OI[i] = (double ***) malloc(chimesRateTables.nei_cooling_table_dimensions[1] * sizeof(double **));
      for (j = 0; j < chimesRateTables.nei_cooling_table_dimensions[1]; j++)
	{
	  chimesRateTables.nei_cooling_OI[i][j] = (double **) malloc(chimesRateTables.nei_cooling_table_dimensions[2] * sizeof(double *));
	  for (k = 0; k < chimesRateTables.nei_cooling_table_dimensions[2]; k++)
	    {
	      chimesRateTables.nei_cooling_OI[i][j][k] = (double *) malloc(chimesRateTables.nei_cooling_table_dimensions[3] * sizeof(double));
	    }
	}
    }
    
  chimesRateTables.chianti_cooling_CII = (double **) malloc(chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double *));
  for (i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[0]; i++)
    chimesRateTables.chianti_cooling_CII[i] = (double *) malloc(chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));
    
  chimesRateTables.chianti_cooling_NII = (double **) malloc(chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double *));
  for (i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[0]; i++)
    chimesRateTables.chianti_cooling_NII[i] = (double *) malloc(chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));  
    
  chimesRateTables.chianti_cooling_SiII = (double **) malloc(chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double *));
  for (i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[0]; i++)
    chimesRateTables.chianti_cooling_SiII[i] = (double *) malloc(chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double)); 
    
  chimesRateTables.chianti_cooling_FeII = (double **) malloc(chimesRateTables.chianti_cooling_table_dimensions[0] * sizeof(double *));
  for (i = 0; i < chimesRateTables.chianti_cooling_table_dimensions[0]; i++)
    chimesRateTables.chianti_cooling_FeII[i] = (double *) malloc(chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));
 
  *this_all_rates = (struct All_rate_variables_structure *) malloc(sizeof(struct All_rate_variables_structure)); 

  (*this_all_rates)->BensRates = (struct Bens_rate_structure *) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(struct Bens_rate_structure));
  
  chimesRateTables.NonEqIon->IonIndexBegin = (int *) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(int));  
  
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {
      if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Hydrogen") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = HI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Helium") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = HeI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Carbon") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = CI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Nitrogen") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = NI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Oxygen") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = OI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Neon") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = NeI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Magnesium") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = MgI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Silicon") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = SiI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Sulphur") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = SI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Calcium") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = CaI;
      else if (strcmp(chimesRateTables.NonEqIon->ElementName[i], "Iron") == 0)
	chimesRateTables.NonEqIon->IonIndexBegin[i] = FeI;
    
      chimesRateTables.NonEqIon->NonEqRates[i].N_Ions = chimesRateTables.NonEqIon->N_Ions[i];
      chimesRateTables.NonEqIon->NonEqRates[i].alpharad = (double **) malloc(chimesRateTables.NonEqIon->N_Ions[i]* sizeof(double *));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	chimesRateTables.NonEqIon->NonEqRates[i].alpharad[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      chimesRateTables.NonEqIon->NonEqRates[i].alphadi = (double **) malloc(chimesRateTables.NonEqIon->N_Ions[i]* sizeof(double *));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	chimesRateTables.NonEqIon->NonEqRates[i].alphadi[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      chimesRateTables.NonEqIon->NonEqRates[i].betacoll = (double **) malloc(chimesRateTables.NonEqIon->N_Ions[i]* sizeof(double *));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	chimesRateTables.NonEqIon->NonEqRates[i].betacoll[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));

      chimesRateTables.NonEqIon->NonEqRates[i].cool = (double **) malloc(chimesRateTables.NonEqIon->N_Ions[i]* sizeof(double *));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	chimesRateTables.NonEqIon->NonEqRates[i].cool[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      
      chimesRateTables.NonEqIon->NonEqRates[i].E_thresh = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));

      chimesRateTables.NonEqIon->NonEqRates[i].cosmicRays = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));

      chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof = (double **) malloc(chimesRateTables.NonEqIon->N_Ions[i]* sizeof(double *));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHionof = (double **) malloc(chimesRateTables.NonEqIon->N_Ions[i]* sizeof(double *));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	chimesRateTables.NonEqIon->NonEqRates[i].CTHionof[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof = (double **) malloc(chimesRateTables.NonEqIon->N_Ions[i]* sizeof(double *));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof = (double **) malloc(chimesRateTables.NonEqIon->N_Ions[i]* sizeof(double *));
      for(j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      
      chimesRateTables.NonEqIon->NonEqRates[i].CTHion_mask = (int *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * sizeof(int));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHrec_mask = (int *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * sizeof(int));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHeion_mask = (int *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * sizeof(int));
      chimesRateTables.NonEqIon->NonEqRates[i].CTHerec_mask = (int *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * sizeof(int));
    
      /* Each element is read in by only 1 task. These                                                                                                                                    
       * will be broadcast to the other tasks later. */
      if (ThisTask == i)
	GetNonEqTables(i, myGlobalVars);
    
      (*this_all_rates)->BensRates[i].CollisIon = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].Recomb = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHrec = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHion = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHerec = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHeion = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].cosmicRays = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
    } 

  MPI_Barrier(MPI_COMM_WORLD);

  /* Array buffers. */
  double *alpharad, *alphadi, *betacoll, *cool;
  double *CTHrecof, *CTHerecof, *CTHionof, *CTHeionof;
  double *twoDcool, *fourDcool;

  /* Broadcast NonEq tables for each element. */
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {
      MPI_Bcast(chimesRateTables.NonEqIon->NonEqRates[i].CTHion_mask, chimesRateTables.NonEqIon->N_Ions[i], MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(chimesRateTables.NonEqIon->NonEqRates[i].CTHrec_mask, chimesRateTables.NonEqIon->N_Ions[i], MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(chimesRateTables.NonEqIon->NonEqRates[i].CTHeion_mask, chimesRateTables.NonEqIon->N_Ions[i], MPI_INT, i, MPI_COMM_WORLD);
      MPI_Bcast(chimesRateTables.NonEqIon->NonEqRates[i].CTHerec_mask, chimesRateTables.NonEqIon->N_Ions[i], MPI_INT, i, MPI_COMM_WORLD);

      /* Read multidimensional arrays into buffers first. */
      alpharad = (double *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      alphadi = (double *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      betacoll = (double *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      cool = (double *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      CTHrecof = (double *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      CTHerecof = (double *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      CTHionof = (double *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));
      CTHeionof = (double *) malloc(chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures * sizeof(double));

      if (chimesRateTables.NonEqIon->IonIndexBegin[i] == CI || chimesRateTables.NonEqIon->IonIndexBegin[i] == NI || chimesRateTables.NonEqIon->IonIndexBegin[i] == SiI || chimesRateTables.NonEqIon->IonIndexBegin[i] == FeI) 
	twoDcool = (double *) malloc(chimesRateTables.chianti_cooling_table_dimensions[0] * chimesRateTables.chianti_cooling_table_dimensions[1] * sizeof(double));

      if (chimesRateTables.NonEqIon->IonIndexBegin[i] == CI || chimesRateTables.NonEqIon->IonIndexBegin[i] == OI)
	fourDcool = (double *) malloc(chimesRateTables.nei_cooling_table_dimensions[0] * chimesRateTables.nei_cooling_table_dimensions[1] * chimesRateTables.nei_cooling_table_dimensions[2] * chimesRateTables.nei_cooling_table_dimensions[3] * sizeof(double));

      if (i == ThisTask)
	{
	  for (j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	    {
	      for (k = 0; k < chimesRateTables.NonEqIon->N_Temperatures; k++)
		{
		  alpharad[j * chimesRateTables.NonEqIon->N_Temperatures + k] = chimesRateTables.NonEqIon->NonEqRates[i].alpharad[j][k];
		  alphadi[j * chimesRateTables.NonEqIon->N_Temperatures + k] = chimesRateTables.NonEqIon->NonEqRates[i].alphadi[j][k];
		  betacoll[j * chimesRateTables.NonEqIon->N_Temperatures + k] = chimesRateTables.NonEqIon->NonEqRates[i].betacoll[j][k];
		  cool[j * chimesRateTables.NonEqIon->N_Temperatures + k] = chimesRateTables.NonEqIon->NonEqRates[i].cool[j][k];
		  CTHrecof[j * chimesRateTables.NonEqIon->N_Temperatures + k] = chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof[j][k];
		  CTHerecof[j * chimesRateTables.NonEqIon->N_Temperatures + k] = chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof[j][k];
		  CTHionof[j * chimesRateTables.NonEqIon->N_Temperatures + k] = chimesRateTables.NonEqIon->NonEqRates[i].CTHionof[j][k];
		  CTHeionof[j * chimesRateTables.NonEqIon->N_Temperatures + k] = chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof[j][k];
		}
	    }
	  if (chimesRateTables.NonEqIon->IonIndexBegin[i] == CI)
	    {
	      for (j = 0; j < chimesRateTables.chianti_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.chianti_cooling_table_dimensions[1]; k++)
		    twoDcool[j * chimesRateTables.chianti_cooling_table_dimensions[1] + k] = chimesRateTables.chianti_cooling_CII[j][k];
		}
	      
	      for (j = 0; j < chimesRateTables.nei_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.nei_cooling_table_dimensions[1]; k++)
		    {
		      for (l = 0; l < chimesRateTables.nei_cooling_table_dimensions[2]; l++)
			{
			  for (m = 0; m < chimesRateTables.nei_cooling_table_dimensions[3]; m++)
			    fourDcool[(j * chimesRateTables.nei_cooling_table_dimensions[3] * chimesRateTables.nei_cooling_table_dimensions[2] * chimesRateTables.nei_cooling_table_dimensions[1]) + (k * chimesRateTables.nei_cooling_table_dimensions[3] * chimesRateTables.nei_cooling_table_dimensions[2]) + (l * chimesRateTables.nei_cooling_table_dimensions[3]) + m] = chimesRateTables.nei_cooling_CI[j][k][l][m];
			}
		    }
		}
	    }
	  else if (chimesRateTables.NonEqIon->IonIndexBegin[i] == OI)
	    {
	      for (j = 0; j < chimesRateTables.nei_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.nei_cooling_table_dimensions[1]; k++)
		    {
		      for (l = 0; l < chimesRateTables.nei_cooling_table_dimensions[2]; l++)
			{
			  for (m = 0; m < chimesRateTables.nei_cooling_table_dimensions[3]; m++)
			    fourDcool[(j * chimesRateTables.nei_cooling_table_dimensions[3] * chimesRateTables.nei_cooling_table_dimensions[2] * chimesRateTables.nei_cooling_table_dimensions[1]) + (k * chimesRateTables.nei_cooling_table_dimensions[3] * chimesRateTables.nei_cooling_table_dimensions[2]) + (l * chimesRateTables.nei_cooling_table_dimensions[3]) + m] = chimesRateTables.nei_cooling_OI[j][k][l][m];
			}
		    }
		}
	    }
	  else if (chimesRateTables.NonEqIon->IonIndexBegin[i] == NI)
	    {
	      for (j = 0; j < chimesRateTables.chianti_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.chianti_cooling_table_dimensions[1]; k++)
		    twoDcool[j * chimesRateTables.chianti_cooling_table_dimensions[1] + k] = chimesRateTables.chianti_cooling_NII[j][k];
		}
	    }
	  else if (chimesRateTables.NonEqIon->IonIndexBegin[i] == SiI)
	    {
	      for (j = 0; j < chimesRateTables.chianti_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.chianti_cooling_table_dimensions[1]; k++)
		    twoDcool[j * chimesRateTables.chianti_cooling_table_dimensions[1] + k] = chimesRateTables.chianti_cooling_SiII[j][k];
		}
	    }
	  else if (chimesRateTables.NonEqIon->IonIndexBegin[i] == FeI)
	    {
	      for (j = 0; j < chimesRateTables.chianti_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.chianti_cooling_table_dimensions[1]; k++)
		    twoDcool[j * chimesRateTables.chianti_cooling_table_dimensions[1] + k] = chimesRateTables.chianti_cooling_FeII[j][k];
		}
	    }
	}
      
      MPI_Bcast(alpharad, chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(alphadi, chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(betacoll, chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(cool, chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(CTHrecof, chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(CTHerecof, chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(CTHionof, chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, i, MPI_COMM_WORLD);
      MPI_Bcast(CTHeionof, chimesRateTables.NonEqIon->N_Ions[i] * chimesRateTables.NonEqIon->N_Temperatures, MPI_DOUBLE, i, MPI_COMM_WORLD);
      
      if (chimesRateTables.NonEqIon->IonIndexBegin[i] == CI || chimesRateTables.NonEqIon->IonIndexBegin[i] == NI || chimesRateTables.NonEqIon->IonIndexBegin[i] == SiI || chimesRateTables.NonEqIon->IonIndexBegin[i] == FeI) 
	MPI_Bcast(twoDcool, chimesRateTables.chianti_cooling_table_dimensions[0] * chimesRateTables.chianti_cooling_table_dimensions[1], MPI_DOUBLE, i, MPI_COMM_WORLD);
      
      if (chimesRateTables.NonEqIon->IonIndexBegin[i] == CI || chimesRateTables.NonEqIon->IonIndexBegin[i] == OI)
	MPI_Bcast(fourDcool, chimesRateTables.nei_cooling_table_dimensions[0] * chimesRateTables.nei_cooling_table_dimensions[1] * chimesRateTables.nei_cooling_table_dimensions[2] * chimesRateTables.nei_cooling_table_dimensions[3], MPI_DOUBLE, i, MPI_COMM_WORLD);
      
      if (i != ThisTask)
	{
	  for (j = 0; j < chimesRateTables.NonEqIon->N_Ions[i]; j++)
	    {
	      for (k = 0; k < chimesRateTables.NonEqIon->N_Temperatures; k++)
		{
		  chimesRateTables.NonEqIon->NonEqRates[i].alpharad[j][k] = alpharad[j * chimesRateTables.NonEqIon->N_Temperatures + k];
		  chimesRateTables.NonEqIon->NonEqRates[i].alphadi[j][k] = alphadi[j * chimesRateTables.NonEqIon->N_Temperatures + k];
		  chimesRateTables.NonEqIon->NonEqRates[i].betacoll[j][k] = betacoll[j * chimesRateTables.NonEqIon->N_Temperatures + k];
		  chimesRateTables.NonEqIon->NonEqRates[i].cool[j][k] = cool[j * chimesRateTables.NonEqIon->N_Temperatures + k];
		  chimesRateTables.NonEqIon->NonEqRates[i].CTHrecof[j][k] = CTHrecof[j * chimesRateTables.NonEqIon->N_Temperatures + k];
		  chimesRateTables.NonEqIon->NonEqRates[i].CTHerecof[j][k] = CTHerecof[j * chimesRateTables.NonEqIon->N_Temperatures + k];
		  chimesRateTables.NonEqIon->NonEqRates[i].CTHionof[j][k] = CTHionof[j * chimesRateTables.NonEqIon->N_Temperatures + k];
		  chimesRateTables.NonEqIon->NonEqRates[i].CTHeionof[j][k] = CTHeionof[j * chimesRateTables.NonEqIon->N_Temperatures + k];
		}
	    }
	  if (chimesRateTables.NonEqIon->IonIndexBegin[i] == CI)
	    {
	      for (j = 0; j < chimesRateTables.chianti_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.chianti_cooling_table_dimensions[1]; k++)
		    chimesRateTables.chianti_cooling_CII[j][k] = twoDcool[j * chimesRateTables.chianti_cooling_table_dimensions[1] + k];
		}
	      
	      for (j = 0; j < chimesRateTables.nei_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.nei_cooling_table_dimensions[1]; k++)
		    {
		      for (l = 0; l < chimesRateTables.nei_cooling_table_dimensions[2]; l++)
			{
			  for (m = 0; m < chimesRateTables.nei_cooling_table_dimensions[3]; m++)
			    chimesRateTables.nei_cooling_CI[j][k][l][m] = fourDcool[(j * chimesRateTables.nei_cooling_table_dimensions[3] * chimesRateTables.nei_cooling_table_dimensions[2] * chimesRateTables.nei_cooling_table_dimensions[1]) + (k * chimesRateTables.nei_cooling_table_dimensions[3] * chimesRateTables.nei_cooling_table_dimensions[2]) + (l * chimesRateTables.nei_cooling_table_dimensions[3]) + m];
			}
		    }
		}
	    }
	  else if (chimesRateTables.NonEqIon->IonIndexBegin[i] == OI)
	    {
	      for (j = 0; j < chimesRateTables.nei_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.nei_cooling_table_dimensions[1]; k++)
		    {
		      for (l = 0; l < chimesRateTables.nei_cooling_table_dimensions[2]; l++)
			{
			  for (m = 0; m < chimesRateTables.nei_cooling_table_dimensions[3]; m++)
			    chimesRateTables.nei_cooling_OI[j][k][l][m] = fourDcool[(j * chimesRateTables.nei_cooling_table_dimensions[3] * chimesRateTables.nei_cooling_table_dimensions[2] * chimesRateTables.nei_cooling_table_dimensions[1]) + (k * chimesRateTables.nei_cooling_table_dimensions[3] * chimesRateTables.nei_cooling_table_dimensions[2]) + (l * chimesRateTables.nei_cooling_table_dimensions[3]) + m];
			}
		    }
		}
	    }
	  else if (chimesRateTables.NonEqIon->IonIndexBegin[i] == NI)
	    {
	      for (j = 0; j < chimesRateTables.chianti_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.chianti_cooling_table_dimensions[1]; k++)
		    chimesRateTables.chianti_cooling_NII[j][k] = twoDcool[j * chimesRateTables.chianti_cooling_table_dimensions[1] + k];
		}
	    }
	  else if (chimesRateTables.NonEqIon->IonIndexBegin[i] == SiI)
	    {
	      for (j = 0; j < chimesRateTables.chianti_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.chianti_cooling_table_dimensions[1]; k++)
		    chimesRateTables.chianti_cooling_SiII[j][k] = twoDcool[j * chimesRateTables.chianti_cooling_table_dimensions[1] + k];
		}
	    }
	  else if (chimesRateTables.NonEqIon->IonIndexBegin[i] == FeI)
	    {
	      for (j = 0; j < chimesRateTables.chianti_cooling_table_dimensions[0]; j++)
		{
		  for (k = 0; k < chimesRateTables.chianti_cooling_table_dimensions[1]; k++)
		    chimesRateTables.chianti_cooling_FeII[j][k] = twoDcool[j * chimesRateTables.chianti_cooling_table_dimensions[1] + k];
		}
	    }
	}
      
      free(alpharad);
      free(alphadi);
      free(betacoll);
      free(cool);
      free(CTHrecof);
      free(CTHerecof);
      free(CTHionof);
      free(CTHeionof);
      if (chimesRateTables.NonEqIon->IonIndexBegin[i] == CI || chimesRateTables.NonEqIon->IonIndexBegin[i] == NI || chimesRateTables.NonEqIon->IonIndexBegin[i] == SiI || chimesRateTables.NonEqIon->IonIndexBegin[i] == FeI) 
	free(twoDcool);
      if (chimesRateTables.NonEqIon->IonIndexBegin[i] == CI || chimesRateTables.NonEqIon->IonIndexBegin[i] == OI)
	free(fourDcool);
      MPI_Barrier(MPI_COMM_WORLD);
    }

  GetPhotoIonTables_parallel(myGlobalVars, N_Elements_in_Bens_tables, this_all_rates, dustG_arr, H2_dissocJ_arr);
}

void init_chimes_parallel(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates, struct Reactions_Structure **this_all_reactions_root, struct Reactions_Structure **this_nonmolecular_reactions_root, double *dustG_arr, double *H2_dissocJ_arr)
{
  int i;

  myGlobalVars->totalNumberOfSpecies = set_species_index_array(myGlobalVars);

  /* Read in Bens tables and the photoionisation
   * tables. */
  initialise_bens_tables_parallel(myGlobalVars, this_all_rates, dustG_arr, H2_dissocJ_arr);

  /* Read in the tables of the additional rates. */
  initialise_additional_rates_tables(myGlobalVars);

  /* Read in tables of equilibrium abundances */
  GetEqAbundancesTables(myGlobalVars);
  
  /* Read in tables for the molecular, CI and OI
   * cooling. */
  initialise_cooling(myGlobalVars);
  
  /* Allocate memory to the reaction lists
   * and initialise them.*/
  *this_all_reactions_root = (struct Reactions_Structure *) malloc(sizeof(struct Reactions_Structure));
  initialise_reactions(*this_all_reactions_root, 1, myGlobalVars, *this_all_rates);

  *this_nonmolecular_reactions_root = (struct Reactions_Structure *) malloc(sizeof(struct Reactions_Structure));
  initialise_reactions(*this_nonmolecular_reactions_root, 0, myGlobalVars, *this_all_rates);
}

void free_all_rates_structure(struct All_rate_variables_structure *this_all_rates, struct globalVariables *myGlobalVars)
{
  int ns, i, j, l, N_arrayCells; 

  /* Allocate memory for the photoionisation rates
   * in the AllRates structure. */
  ns = 0;
  for (l = 0; l < chimesRateTables.NonEqIon->N_Elements; l++)
    {  
      if (l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
	{
	  if (chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
	    N_arrayCells = 3;   /* Hydrogen photoionisation also includes H- & H2*/
	  else
	    N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

	  for (j = 0; j < N_arrayCells; j++)
	    free(this_all_rates->BensRates[ns].PhotoIon[j]); 
	  free(this_all_rates->BensRates[ns].PhotoIon); 
	  
	  ns += 1;
	}
    }

  for (i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {
      free(this_all_rates->BensRates[i].CollisIon); 
      free(this_all_rates->BensRates[i].Recomb);
      free(this_all_rates->BensRates[i].CTHrec);
      free(this_all_rates->BensRates[i].CTHion);
      free(this_all_rates->BensRates[i].CTHerec);
      free(this_all_rates->BensRates[i].CTHeion);
      free(this_all_rates->BensRates[i].cosmicRays);
    }

  free(this_all_rates->BensRates); 

  free(this_all_rates); 

  return; 
} 

void free_reactions_list(struct Reactions_Structure *root_node) 
{
  struct Reactions_Structure *current_reaction, *prev_reaction; 
  
  current_reaction = root_node; 
  while (current_reaction->next_reaction != NULL) 
    {
      prev_reaction = current_reaction; 
      current_reaction = current_reaction->next_reaction; 
      free(prev_reaction); 
    }
  free(current_reaction); 
  root_node = NULL; 

  return; 
} 

void init_chimes_omp(struct globalVariables *myGlobalVars, struct All_rate_variables_structure **this_all_rates, struct Reactions_Structure **this_all_reactions_root, struct Reactions_Structure **this_nonmolecular_reactions_root)
{

  int i,j, l; 
  hid_t file_id, dataset, datatype;
  int N_Elements_in_Bens_tables; 
  char fname[256];

  if (ThisTask == 0) 
    {
      /* Task0 reads in header info. */
      sprintf(fname, "%sIonRates_Head.HM01.hdf5", myGlobalVars->BenTablesPath);		/* Note that we ignore everything that depends on the UV spectrum, so we can pick any of them */

      file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

      if(file_id < 0)
	printf("[read_noneq_tables()]: unable to open file %s\n", fname);

      dataset = H5Dopen(file_id, "N_Elements");
      H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &N_Elements_in_Bens_tables);
      H5Dclose(dataset);
    }

  /* Broadcast header info to all tasks. */
  MPI_Bcast(&N_Elements_in_Bens_tables, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (ThisTask == 0)
    {
      /* Task0 closes header table. */
      H5Fclose(file_id);
    }
 
  *this_all_rates = (struct All_rate_variables_structure *) malloc(sizeof(struct All_rate_variables_structure)); 

  (*this_all_rates)->BensRates = (struct Bens_rate_structure *) malloc(chimesRateTables.NonEqIon->N_Elements * sizeof(struct Bens_rate_structure));
  
  for(i = 0; i < chimesRateTables.NonEqIon->N_Elements; i++)
    {    
      (*this_all_rates)->BensRates[i].CollisIon = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].Recomb = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHrec = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHion = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHerec = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].CTHeion = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
      (*this_all_rates)->BensRates[i].cosmicRays = (double *) malloc((chimesRateTables.NonEqIon->N_Ions[i] - 1) * sizeof(double));
    } 

  /* Allocate memory for the photoionisation rates
   * in the AllRates structure. */
  int ns, N_arrayCells; 
  ns = 0;
  for (l = 0; l < N_Elements_in_Bens_tables; l++)
    {  
      if (l < 2 || myGlobalVars->element_included[max(l - 2, 0)] == 1)
	{
	  if (chimesRateTables.NonEqIon->IonIndexBegin[ns] == HI)
	    N_arrayCells = 3;   /* Hydrogen photoionisation also includes H- & H2*/
	  else
	    N_arrayCells = chimesRateTables.NonEqIon->N_Ions[ns] - 1;

	  (*this_all_rates)->BensRates[ns].PhotoIon = (double **) malloc(N_arrayCells * sizeof(double *));
	  for (j = 0; j < N_arrayCells; j++)
	    (*this_all_rates)->BensRates[ns].PhotoIon[j] = (double *) malloc(chimesRateTables.NonEqIon->N_Auger[ns] * sizeof(double));
	  
	  ns += 1;
	}
    }

  *this_all_reactions_root = (struct Reactions_Structure *) malloc(sizeof(struct Reactions_Structure));
  initialise_reactions(*this_all_reactions_root, 1, myGlobalVars, *this_all_rates);

  *this_nonmolecular_reactions_root = (struct Reactions_Structure *) malloc(sizeof(struct Reactions_Structure));
  initialise_reactions(*this_nonmolecular_reactions_root, 0, myGlobalVars, *this_all_rates);
}





