/*Header file for the source code*/

extern void input_data_test(double kappa, double A, int time_step, double dt, int Nx, int Ny, double dx, double dy, double c_zero, double c_noise);

extern void cahn_hilliard_evolution(double kappa, double A, int time_step, double dt, int Nx, int Ny, double dx, double dy, fftw_complex *comp);

extern void ps_file(char* file_name, char* ps_file_name, int Nx, int Ny, int N);