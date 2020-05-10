/*
    This program is used to generate morphological evolution pattern
    based on the Allen-Cahn equation (which is used for a non-conserved
    order parameter).

    Copyright (C) 2020  Abhinav Roy, M.P. Gururajan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/
//Including all the requisite headre files.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include "fftw3.h"
#include<gsl/gsl_rng.h>
#include "../header_file/headers.h"

//The main program begins.
int main(void)
{
  FILE *fr, *fw;
  int Nx, Ny;
  double dx,dy;
  double A,kappa;
  int time_step;
  double dt;
  double phi_zero, phi_noise;
  double average_phi;
  fftw_complex *phi;
  gsl_rng * ran_num;
  const gsl_rng_type * Taus;
  int i1,i2;

  //Clearing the output directory in order to print fresh result from every run of the code.

  (void) system("rm -rf output/*");

  //Opening the file simulation_data where all the physical constants and system information related to the simulation is stored.

  if ((fw = fopen("output/simulation_data","w")) == NULL)
  {
    printf("Unable to open output/simulation_data. Exiting.");
    exit(0);
  }
  else
  {
    fw = fopen("output/simulation_data","w");
  }

  {
    if ((fr = fopen("input/constants","r")) == NULL)
    {
      printf("Unable to open input/constants. Exiting.");
      exit(0);
    }
    else
    {
      fr = fopen("input/constants","r");
    }
    (void)fscanf(fr,"%le%le",&kappa,&A);
    (void)fclose(fr);
    fprintf(fw,"kappa = %le\nA = %le\n",kappa,A);
  }

  //Taking the input for time information.

  if ((fr = fopen("input/time_info","r")) == NULL)
  {
    printf("Unable to open input/time_info. Exiting.");
    exit(0);
  }
  else
  {
    fr = fopen("input/time_info","r");
  }
  (void)fscanf(fr,"%d%le",&time_step,&dt);
  (void)fclose(fr);
  fprintf(fw,"time_step = %d\ndt = %le\n",time_step,dt);

  //Taking the input for system information, i.e.,grid points and grid spacing.

  if ((fr = fopen("input/system_info","r")) == NULL)
  {
    printf("Unable to open input/system_info. Exiting.");
    exit(0);
  }
  else
  {
    fr = fopen("input/system_info","r");
  }
  (void)fscanf(fr,"%d%d%le%le",&Nx,&Ny,&dx,&dy);
  fclose(fr);
  (void)fprintf(fw,"Nx = %d\nNy = %d\n", Nx, Ny);
  (void)fprintf(fw,"dx = %le\ndy = %le\n", dx, dy);

  //Taking the input for order parameter profile, i.e., nominal order parameter and order parameter noise strength.

  if ((fr = fopen("input/order_parameter_profile","r")) == NULL)
  {
    printf("Unable to open input/order_parameter_profile. Exiting.");
    exit(0);
  }
  else
  {
    fr = fopen("input/order_parameter_profile","r");
  }
  (void)fscanf(fr,"%le%le", &phi_zero, &phi_noise);
  (void)fclose(fr);
  fprintf(fw,"phi_zero = %le\nphi_noise = %le\n", phi_zero, phi_noise);

  //Closing the file pointer for writing the simulation data.

  (void)fclose(fw);
  fflush(fw);
  //Performing the input data test to make sure all the simulation data provided are within acceptable range.

  input_data_test(kappa, A, time_step, dt, Nx, Ny, dx, dy, phi_zero, phi_noise);

  /*Defining the complex composition variable of fftw_complex type.
  For more information, kindly refer to FFTW manual (version 3.3.4).*/

  phi = fftw_malloc(Nx*Ny* sizeof(fftw_complex));

  /*Defining the Tausworthe uniform random number generator.
  For more information, kindly refer to the GNU Scientific Library (GSL) documentation.*/

  (void) gsl_rng_env_setup();
  Taus = gsl_rng_taus;
  ran_num = gsl_rng_alloc (Taus);

  //Setting the initial order parameter profile.

  average_phi = 0.0;
  for(i1=0; i1 < Nx; i1++)
  {
    for(i2=0; i2 < Ny; i2++)
    {
      __real__(phi[i2+Ny*i1]) = phi_zero + phi_noise*(0.5 - gsl_rng_uniform_pos(ran_num));
      __imag__(phi[i2+Ny*i1]) = 0.0;
      average_phi = average_phi + __real__(phi[i2+Ny*i1]);
    }
  }
  average_phi = average_phi/(Nx*Ny);

  //Writing the average order parameter as part of simulation data.

  if ((fw = fopen("output/simulation_data","a")) == NULL)
  {
    printf("Unable to open output/simulation_data. Exiting.");
    exit(0);
  }
  else
  {
    fw = fopen("output/simulation_data","a");
  }
  fprintf(fw,"average_phi = %le\n", average_phi);
  (void)fclose(fw);
  //Performing the evolution using the Allen-Cahn equation.
  allen_cahn_evolution(kappa, A, time_step, dt, Nx, Ny, dx, dy, phi);
  //free the composition and the random number variable.
  fftw_free(phi);
  gsl_rng_free(ran_num);

  return 0;
}
/*------------------------------------------------End of CODE---------------------------------------------------------*/
