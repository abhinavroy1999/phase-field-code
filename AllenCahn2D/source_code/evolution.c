/*
    This code performs evolution of the order parameter guided by the Allen-Cahn equaiton.

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
#include<math.h>
#include<complex.h>
#include "fftw3.h"
#include<gsl/gsl_math.h>
#include "../header_file/headers.h"

void allen_cahn_evolution(double kappa, double A, int time_step, double dt, int Nx, int Ny, double dx, double dy, fftw_complex *phi)
{
  FILE *fw;
  char file_name[50];
  char ps_file_name[50];
  double *c;
  int halfNx, halfNy;
  float kx, ky, delkx, delky;
  double k2;
  double denominator;
  int i1, i2;
  int temp = 0;
  fftw_complex *g;
  //double pointer variable for writing the order parameter profile to a file.
 c = (double *)malloc((size_t) Nx*Ny*sizeof(double));
 g = fftw_malloc(Nx*Ny*sizeof(fftw_complex));
  for (i1 = 0; i1<Nx; ++i1)
  {
    for (i2 = 0; i2<Ny; ++i2)
    {
      c[i2 + i1*Ny] = __real__(phi[i2 + i1*Ny]);
    }
  }
  sprintf(file_name, "output/time%d.dat", temp);
  fw = fopen(file_name, "wb");
  fwrite(&c[0], sizeof(double),(size_t) Nx*Ny,fw);
  (void) fclose (fw);
  fflush(fw);
  sprintf(ps_file_name,"output/time%d.ps",temp);
  ps_file(file_name, ps_file_name, Nx, Ny, Nx);


  //Defining the plans for fourier transforms.
  fftw_plan plan1, plan2, plan3;
  plan1 = fftw_plan_dft_2d(Nx, Ny, phi, phi, FFTW_FORWARD, FFTW_ESTIMATE);
  plan2 = fftw_plan_dft_2d(Nx, Ny, g, g, FFTW_FORWARD, FFTW_ESTIMATE);
  plan3 = fftw_plan_dft_2d(Nx, Ny, phi, phi, FFTW_BACKWARD, FFTW_ESTIMATE);



  //half of the grid points required for Periodic Boundary Condition (PBC) to get rid of the surface effects during simulation.
  halfNx = (int) Nx/2;
  halfNy = (int) Ny/2;

  //delta Kx and delta Ky - Fourier space vectors.
  delkx = (2*M_PI)/(Nx*dx);
  delky = (2*M_PI)/(Ny*dy);

  //Loop for temporal evolution.
  for (temp = 1; temp <time_step+1; ++temp)
  {
    //Defining the derivative of local free energy densty with respect to the composition.
    for(i1=0; i1<Nx; ++i1)
    {
      for(i2=0; i2<Ny; ++i2)
      {
        g[i2+Ny*i1] = 2.0*A*((phi[i2 + i1*Ny]))*(1.0-(phi[i2 + i1*Ny]))*(1.0-2.0*(phi[i2 + i1*Ny]));
      }
    }

    fftw_execute(plan1);
    fftw_execute(plan2);
    //Loop for implementing the PBC and subsequent evolution of the composition profile.
    for(i1=0; i1 < Nx; ++i1)
    {
	     if(i1 < halfNx)
       { kx = i1*delkx; }
	      else
       { kx = (i1-Nx)*delkx; }
       for(i2=0; i2 < Ny; ++i2)
       {
	        if(i2 < halfNy)
          { ky = i2*delky; }
	        else
          { ky = (i2-Ny)*delky; }
	         k2 = kx*kx + ky*ky;
           denominator = (1.0 + 2.0*kappa*k2*dt);
          phi[i2 + Ny*i1] = (phi[i2+Ny*i1] - dt*g[i2 + Ny*i1])/denominator;
        }
    }

    fftw_execute(plan3);
    //Determining the normalized order parameter (required for the fft algorithm to work).
    for(i1=0; i1<Nx; ++i1)
    {
          for(i2=0; i2<Ny; ++i2)
          {
           phi[i2+Ny*i1] = phi[i2+Ny*i1]/(Nx*Ny);
          }
    }

    //Taking the output every 1000 units interval of time.
      if (temp%1000 == 0)
      {
        for (i1 = 0; i1<Nx; ++i1)
        {
          for (i2 = 0; i2<Ny; ++i2)
          {
            c[i2 + i1*Ny] = __real__(phi[i2 + i1*Ny]);
          }
        }

        sprintf(file_name,"output/time%d.dat",temp);
        fw = fopen(file_name,"wb");
        fwrite(&c[0], sizeof(double),(size_t)Nx*Ny, fw);
	      sprintf(ps_file_name,"output/time%d.ps",temp);
        ps_file(file_name, ps_file_name, Nx, Ny, Nx);
        (void) fclose (fw);
	      fflush(fw);
      }
  }
  fftw_free(g);
  free(c);
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  fftw_destroy_plan(plan3);

}
/*-------------------------------------------End of CODE------------------------------------------------------------------*/
