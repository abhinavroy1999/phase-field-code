/*  
    Code to perform test on the input data to determing whether the given
    simulation data is in the acceptable range.

    This program is used to generate morphological evolution pattern
    based on the Cahn-Hilliard equation (which is used for a conserved
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

#include<stdio.h>
#include<stdlib.h>

void input_data_test(double kappa, double A, int time_step, double dt, int Nx, int Ny, double dx, double dy, double c_zero, double c_noise)
{
  FILE *fw;
  int index = 0;
  //Check the gradient free energy term.
  if (kappa<=0)
  {
    printf("The value of gradient free energy coefficient term(kappa) cannot be 0 or negative. Kindly change the value accordingly.\n");
    index=1;
  }
  //Check the value of A (constant determining the height of the free energy curve).
  if (A<=0)
  {
    printf("The value of A cannot be 0 or negative. Kindly change the value accordingly.\n");
    index=1;
  }
  //Check the value of total time of evolution.
  if (time_step<=0)
  {
    printf("The value of total time of evolution cannot be 0 or negative. Kindly change the value accordingly.\n");
    index=1;
  }
  //Check the value of time-step.
  if (dt<=0)
  {
    printf("The value of time steps cannot be 0 or negative. Kindly change the value accordingly.\n");
    index=1;
  }
  //Check the value of grid points.
  if (Nx<=0||Ny<=0)
  {
    printf("The value of grid points cannot be 0 or negative. Kindly change the value accordingly.\n");
    index=1;
  }
  //Check the value of grid spacings.
  if (dx<=0||dy<=0)
  {
    printf("The value of grid spacing cannot be 0 or negative. Kindly change the value accordingly.\n");
    index+=1;
  }
  //Check the value of nominal alloy composition.
  if (c_zero<=0||c_zero>=1)
  {
    printf("The average alloy composition should lie between 0 and 1 (excluding the extreme points). Kindly change the value accordingly.\n");
    index=1;
  }
  //Check the value of composition noise strength.
  if (c_noise>=1)
  {
    printf("The strength of noise in the composition should not be greater than 1. Kindly change the value accordingly.\n");
    index=1;
  }
  if (index == 1)
  {
    printf("change the respective value accordingly. Exiting.\n");
    exit(0);
  }
  else //The case where all the tests have passed.
  {
    if( (fw = fopen("output/input_test_code_result","w")) == NULL)
    {
		    printf("Unable to open output/input_test_code_result. Exiting\n");
		    exit(0);
	  }
	else
    {
	     fw = fopen("output/input_test_code_result","w");
    }
    (void) fprintf(fw, "The input data has passed the test.\n");
    (void) fprintf(fw, "The input data are within correct range to perform simulation.\n");
    fclose(fw);
  }

}
/*---------------------------------------------End of CODE--------------------------------------------------------------------*/
