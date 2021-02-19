#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<unistd.h>
#include<sys/wait.h>
#include<sys/types.h>
#include "mt19937ar.h"

#define Kb 1.3806503E-23 //[J/K]
#define h 6.626068E-34 //[J*s]
#define c 2.9979E+10 //[cm/s]
#define Rgas 8.314472 //[J/mol/K]
#define Na 6.0221E+23 //Avogadro number [at/mol]

unsigned long csSeed;
int E;
int J;
double T,s;
char quantum_method[20];
char basis_set[20];
char restart[20];
char start_log[20];
int ATOMS, step;
double Internal_Energy;
double Kinetic_Energy;
double EBarr_kcalmol;
double Delta_Energy_rea;
double Delta_Energy_pro;
double Tinitial;
double Tsteps;
double Tdelta;

int MAXSTEP;
double deltas;//=0.0104;
int E_dim=150000;
int NpointsInt;
int rearrange;
int intcoord;
int dsmethod;
int saddlepbkf;
int saddlep;
int MAXSTEP2TS;
double Maxtdev;
int stepmin;
int isctvtst;
int izct;
int brows,bcolumns;
double redmu;
double iminfreq;

int nrot;
double angular_velocity_x , angular_velocity_y, angular_velocity_z;  
int onlyrotors, numrotors;
int  *pivotA, *pivotB, *numatomsintopA;
int **atomsintopA;
int **atomsintopB;
int **igroupA;
int **igroupB;


double dt;
double redmass_TS;
double **incoox,**incooy,**incooz;
double **distance;
double **force_constants;
double **Bmat;
double **Amat;
double **Gmat;
double **ABLmat;
double ***Cmat;
double **Lcartint;
double *Freqint;
double *lambdaint;

int bonds_number;
int bonds_angle_number;
int dihedral_angle_number;
double *TR1, *TR2, *TR3, *TR4, *TR5, *TR6;
double **D;
double *Rx_coord;
//int dim_truong_mueff;
double *mueff_read, *x_mueff_read;
int dim_newpes;
double *pes_read, *x_pes_read;
double *mass;
double proj_rea_coo;
int an_freq;

double     pi;
double     c_light;
double     c_light_cm_s;
// Magnetic constant or permeability of vacuum [N/A**2]
double     mu_perm ; 

// Electric constant or permittivity of vacuum [F/m]
double     permittivity; 
//_VERSION == 2006
// Recommended fundamental constants of physics
// and chemistry based on the 2006 adjustment
 
// Planck constant [J*s]
double     h_planck;
double     h_bar; 
//Elementary charge [C]
double     e_charge; 

//     ! Electron mass [kg]
double     e_mass; 
 
//     ! Proton mass [kg]
double     p_mass;    
 
//     ! Electron g factor [ ]
double     e_gfactor; 
 
//     ! Fine-structure constant
// !MK a_fine = 0.5 *mu_perm*c_light*e_charge**2/h_planck
double     a_fine ;

//     ! Rydberg constant [1/m]
// !MK rydberg = 0.5 *e_mass*c_light*a_fine**2/h_planck
double     rydberg ;
 
//     ! Avogradro constant [1/mol]
double     n_avogadro; 
 
//     ! Boltzmann constant [J/K]
double     boltzmann  ;

//     ! Atomic mass unit [kg]; conversion factor [u] -> [kg]
double     a_mass  ;
 
//    ! Bohr radius [m]
// !MK a_bohr = a_fine/(4.0 *pi*rydberg)
double     a_bohr ; 

 
//     ! Conversion factors
 
//     ! [u] -> [a.u.]
double     massunit;
//     ! [Bohr] -> [Angstrom]
double angstrom;
//     ! [Angstrom] -> [Bohr]
double bohr ;
//     ! [a.u.] -> [s]
double seconds;
//     ! [a.u.] -> [fs]
double femtoseconds;
//     ! [a.u.] -> [ps]
double picoseconds;
//     ! [a.u.] -> [J]
double joule;
//     ! [a.u.] -> [K]
double kelvin;
//     ! [a.u.] -> [kJ/mol]
double kjmol ;
//     ! [a.u.] -> [kcal/mol]
double kcalmol ;
//     ! [a.u.] -> [Pa]
double pascal ;
//     ! [a.u.] -> [bar]
double bar  ;
//     ! [a.u.] -> [atm]
double atm; 
//     ! [a.u.] -> [eV]
double evolt;
//     ! [a.u.] -> [Hz]
double hertz;
//     ! [a.u./Bohr**2] -> [1/cm] (wave numbers)
double vibfac;
//     ! [a.u.] -> [1/cm] (wave numbers)
double wavenumbers;
 
double ***Hessian;
double **positions, **position_derivative;
double **gradient, **gradient_derivative;

double *CM_x,*CM_y, *CM_z;
double ***I;
double ***v, **eigen;
double *lambda;
double *Energy;
double *mueff_mu;

//double **L;
//double **D1;
typedef struct {
  char *atom_name;
  double atomic_mass;
  int at_numb;
  double *ext_coord_x;
  double *ext_coord_y;
  double *ext_coord_z;
  double ext_coord_x_bohr;
  double ext_coord_y_bohr;
  double ext_coord_z_bohr;
  double *grad_x;
  double *grad_y;
  double *grad_z;
  double force_x_old;
  double force_y_old;
  double force_z_old;
  double velocity_x;
  double velocity_y;
  double velocity_z;
  double vib_velocity_x;
  double vib_velocity_y;
  double vib_velocity_z;
  double LRot1x;
  double LRot1y;
  double LRot1z;
  double LRot2x;
  double LRot2y;
  double LRot2z;
  double LRot3x;
  double LRot3y;
  double LRot3z;
  double LTr1x;
  double LTr1y;
  double LTr1z;
  double LTr2x;
  double LTr2y;
  double LTr2z;
  double LTr3x;
  double LTr3y;
  double LTr3z;


}atom_struct;

atom_struct* atoms_data;

typedef struct {
  double distance;
  int atom1;
  int atom2;
}bond_struct;

bond_struct* bonds_data;

typedef struct {
  double angle_value;
  int atom1;
  int atom2;
  int atom3;

}bond_angle_struct;

bond_angle_struct* bond_angle;

typedef struct {
  double angle_value;
  int atom1;
  int atom2;
  int atom3;
  int atom4;

}dihedral_angle_struct;

dihedral_angle_struct* dihedral_angle;

typedef struct {
  double positions;
  double velocities;
  char name;

}internal_coordinates_struct;

internal_coordinates_struct* internal_coordinates;

void read_inputfile ();
void read_Bmat_Cmat ();
void calc_Amat ();
void intcoord_hessian(int step, double **FC, double *grad);
void determine_top_atoms ();
void read_VaG_mueff ();
void allocation( int dim );
void convert_coordinates_to_bohr(int step);
void calculate_center_of_mass(int step);
void set_origin_to_cofm(int step );
void calculate_inertia(int step);
void coordinates_in_principal_axes(int step);
void jacobi(double **a, int n, double *d, double **v, int nrot);
void eigsrt(double d[], double **v, int n);
void cross_prod(double ax, double ay, double az, double bx, double by, double bz, double *px, double *py, double *pz);
double ** force_constants_mass_weight(double **FC);
double **diagonalize_hessian(int dim, double**FC);
double *calc_freq(double *lam, double *freq);
double calc_var_corr(int dim, double T, double **frequencies_save);
double **invert_inertia_matrix( double **matrix, int dim, double **inverse_matrix);
double **projector_matrix(double **P, int step);
void projector_matrix_Rot(int step);
void mass_weight_coordinates(int step);
void save_displ(double *eval,double **evect,int ist,int dim,int step,FILE *fname);
double **prod_mat(double **A, double **B, int dim);
double **prod_mat2(double **A, int dim1, int dim2, double **B, int dim3, int dim4);
double **invert_matrix(double **A, int dim);
double **transpose_double(double**A, int l_side, int h_side,double **B);
double dot_prod(double *v1,double *v2, int dim);
void calc_reaction_coordinate_direction_vector(int step);
double **calc_deriv_matrix(int step,double **Matrix_back,double **Matrix_for , double **Deriv_Matrix ,int dim);
double *calc_BkF(int dim,int step, double **L_int, double *BkF);
int compare(const void *v1 ,const void *v2);
void eigen_QL(double **a, int n, double d[], double e[], double **z);
void tred2(double **a, int n, double d[], double e[]);
void tqli(double d[], double e[], int n, double **z);
double **inverse_matrix(double **a, int n, double **inv);
double ** force_constants_mass_unweight(double **FC);
void read_pesrxfile ();
void read_pesfile ();
void read_rxfile ();
//void reading_newpes_file ();
void write_traj();
void smooth_hess();
void grad_hess();
double ***Hessian_deriv;
int calc_VaG(double **frequencies_save);
double *calc_transmission_coeff(double *mueff,int Emax, double *Transmission_coeff );
double find_sl_match(double E, double dstep);
double find_sr_match(double E, double dstep);
double integrate_function(int Efirst, int Elast, double *function);
//void read_mueff_file (FILE *fp);
void deriv_position(void);
double trapzd(double *x,double *function,double *function2, int a,int b,int n, double olds, double s);
int max(int a, int b);
double splint(double *xa, double *ya, double *y2a, int n, double x, double y);
double *spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
double interpolate_grad(int index, int step);
void deriv_gradient(void);
double ***check_L_sign(int step, double ***L_int_save);
int sign(double q);
void check_L_order(int step, double ***L_int_save);
void interpolate_freqs(double **freqs);
void interpolate_hess(double **freqs);
