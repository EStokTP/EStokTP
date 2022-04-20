/******************************************************************************
 * January 2021
 *
 * Daniela Polino
 * Carlo Cavallotti
 * projection of hindered rotors from Hessian
 * tunneling and coupling constants through RPH - SCT theories
 * frequencies in internal coordinates 
 *
 *****************************************************************************/
#include "RPHt.h"
#include "nrutil.h"

int main(void) {
  
// Constants

  pi=3.141592653589793;          // greek pi
  c_light = 299792458.0;         //Speed of light in vacuum [m/s]
  c_light_cm_s = 29979245800.0;  //Speed of light in vacuum [cm/s]
  mu_perm = 4.0 *pi*1.0e-7;      // Magnetic constant or permeability of vacuum [N/A**2]
  permittivity = 1.0 /(mu_perm*c_light*c_light);  // Electric constant or permittivity of vacuum [F/m]

// CODATA_VERSION == 2006
// Recommended fundamental constants of physics and chemistry based on the 2006 adjustment
 
  h_planck = 6.62606896e-34;     // Planck constant [J*s]
  h_bar = h_planck/(2.0*pi);     // Planck constant/2/pi [J*s]
  e_charge = 1.602176487e-19;    // Elementary charge [C]
  e_mass = 9.10938215e-31;       // Electron mass [kg]
  p_mass = 1.672621637e-27;      // Proton mass [kg]
  e_gfactor = -2.0023193043622;  // Electron g factor [ ]
  a_fine = 7.2973525376e-3;      // Fine-structure constant: MK a_fine = 0.5 *mu_perm*c_light*e_charge**2/h_planck
  rydberg = 10973731.568527 ;    // Rydberg constant [1/m] MK rydberg = 0.5 *e_mass*c_light*a_fine**2/h_planck
  n_avogadro = 6.02214179e+23;   // Avogadro constant [1/mol]
  boltzmann = 1.3806504e-23;     // Boltzmann constant [J/K]
  a_mass = 1.660538782e-27;      // Atomic mass unit [kg]; conversion factor [u] -> [kg]
  a_bohr = 0.52917720859e-10;    // Bohr radius [m] MK a_bohr = a_fine/(4.0 *pi*rydberg)

 //     ! Conversion factors
 
  massunit = a_mass/e_mass;                       //  [u] -> [a.u.]
  angstrom = 1.0e+10 *a_bohr;                     //  [Bohr] -> [Angstrom]
  bohr = 1.0 /angstrom;                           //  [Angstrom] -> [Bohr]
  seconds = 1.0 /(4.0 *pi*rydberg*c_light);       //  [a.u.] -> [s] time in a.u. = 2.418884e-17 s
  femtoseconds = 1.0e+15 *seconds;                //  [a.u.] -> [fs]
  picoseconds = 1.0e+12 *seconds;                 //  [a.u.] -> [ps]
  joule = 2.0 *rydberg*h_planck*c_light;          //  [a.u.] -> [J]
  kelvin = joule/boltzmann;                       //  [a.u.] -> [K]
  kjmol = 0.001 *joule*n_avogadro;                //  [a.u.] -> [kJ/mol]
  kcalmol = kjmol/4.184;                          //  [a.u.] -> [kcal/mol]
  pascal = joule/(a_bohr*a_bohr*a_bohr);          //  [a.u.] -> [Pa]
  bar = pascal/1.0e+5;                            //  [a.u.] -> [bar]
  atm = pascal/1.013250e+5;                       //  [a.u.] -> [atm]
  evolt = joule/e_charge;                         //  [a.u.] -> [eV]
  hertz = joule/h_planck;                         //  [a.u.] -> [Hz]
  vibfac = 5.0 *sqrt(kjmol)/(pi*a_bohr*c_light);  //  [a.u./Bohr**2] -> [1/cm] (wave numbers)
  wavenumbers = 0.02 *rydberg;                    //     ! [a.u.] -> [1/cm] (wave numbers)

 /* Code Structure
  *
  * (1) Reading input file
  *
  */

  //  FILE *freq_results,*freq_results1, *freq_orig_results, *Km_res, *P_E;
  //  FILE *mueff_results,*L_save,*L_orig_save,*Cart_displ,*anim_freq; 
  //  FILE *Cart_displ; 
  int step;
  int i,j;
  double **L_int;
  double ***L_int_save;
  double *BkF;
  double *k_ome;
  double **Km,*K,*t, *a, *t_dev;
  FILE *Cart_displ; 
  //  double **Km,*K,*t, *a, *t_dev,*mueff_mu;

  //Reading input files: RPHt_input.dat, RPHt_PES.dat, Rx_coord.txt
  //the hessian and gradient inputs must be in Hartree/Bohr^2 and Hartree/Bohr (not mass weighted)

  read_inputfile();
  // choice for computation of average centrifigual component
  // 0 = SCT
  // 1 = RPH or similar (to be checked)
  dsmethod=0;

  if(onlyrotors==1){
    read_pesrxfile();    
    write_traj();
  }

  if(intcoord==1){
    read_Bmat_Cmat();
    calc_Amat();
  }
  
  //  read_newpes_file();

  // save trajectories in xyz format redable by VMD and Molden

  // allocate vectors

  int dim=3*ATOMS;
  double *frequencies;
 
  L_int_save=(double ***) malloc( MAXSTEP*sizeof(double **) );
  for(step=0;step<MAXSTEP;step++){ L_int_save[step]=(double **) malloc( (int)(dim)*sizeof(double *) );}
  for(step=0;step<MAXSTEP;step++){ for(j=0; j <(int)(dim); j++) { L_int_save[step][j]= (double*) malloc((int)(dim)*sizeof(double) );} }

  Hessian_deriv=(double ***) malloc( MAXSTEP*sizeof(double **) );for(step=0;step<MAXSTEP;step++){ Hessian_deriv[step]=(double **) malloc( (int)(dim)*sizeof(double *) );}

  for(step=0;step<MAXSTEP;step++){ for(j=0; j <(int)(dim); j++) { Hessian_deriv[step][j]= (double*) malloc((int)(dim)*sizeof(double) );} }

  double **frequencies_save;
  frequencies_save =  (double **) malloc((int)(MAXSTEP)* sizeof(double*));
  for(step=0; step<MAXSTEP; step++) { frequencies_save[step]=  (double *) malloc((int)(3*ATOMS)* sizeof(double));}
 
  for(step=0; step<MAXSTEP; step++) {
    for(i=0; i<dim; i++) {frequencies_save[step][i]=0. ;}
  }
  
  Km=(double **) malloc( MAXSTEP*sizeof(double *) );
  for(j=0; j < MAXSTEP; j++) {    Km[j]= (double*) malloc( (int)(dim-7)*sizeof(double) ); }
  K= (double*) malloc( MAXSTEP*sizeof(double) );
  t= (double*) malloc( MAXSTEP*sizeof(double) );
  t_dev= (double*) malloc( MAXSTEP*sizeof(double) );
  a= (double*) malloc( MAXSTEP*sizeof(double) );
  k_ome= (double*) malloc( MAXSTEP*sizeof(double) );
  mueff_mu= (double*) malloc( MAXSTEP*sizeof(double) );

  for(step=0; step<MAXSTEP; step++) {
    k_ome[step]=0 ;
    K[step]=0;
  }




  //**** beginning of cycle over different input structures (steps)****

  // open output files for the step cycle

  //  if( onlyrotors==1 || MAXSTEP==1) {

  //  if((freq_results1=fopen("freqs1.txt","w"))==NULL) {
  //    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","freqs1.txt");
  //    exit(1);
  //  }
  
  //  if((L_orig_save=fopen("L_orig_save.txt","w"))==NULL) {
  //    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","L_orig_save.txt");
  //    exit(1);
  //  }

  /*
  FILE *Cart_displ; 
  if((Cart_displ=fopen("Cart_displ.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","Cart_displ.txt");
    exit(1);
  }
  */
  //  if((freq_orig_results=fopen("freqs_orig.txt","w"))==NULL) {
  //    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","freqs_orig.txt");
  //    exit(1);
  //  }

  //  if((freq_results=fopen("freqs.txt","w"))==NULL) {
  //    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","freqs.txt");
  //    exit(1);
  //  }
  //  if((Km_res=fopen("Km_res.txt","w"))==NULL) {
  //    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","Km_res.txt");
  //    exit(1);
  //  }
  //  if((anim_freq=fopen("anim_freq.xyz","w"))==NULL) {
  //    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","anim_freq.txt");
  //    exit(1);
  //  }

    //  }

  //the purpose of the cycle is the diagonalization of the Hessian and its projection

  if( (onlyrotors==1) || (onlyrotors==0)){
    for(step=0;step<MAXSTEP;step++){

      printf("Step Number %d\n",step);

    //Converto le coordinate da Angstrom a Bohr
      convert_coordinates_to_bohr(step);

    //calcolo il centro di massa
      calculate_center_of_mass(step);
  
    //pongo l'origine degli assi pari al centro di massa
      set_origin_to_cofm(step);

    //calculate inertia moments in amu x bohr^2 and eigenvectors 
      calculate_inertia(step); 

    // don't use, this function must be checked; it has been check it works now
    // be careful, the result is sensitive to the initial orientation.
    //coordinates_in_principal_axes(step);


    //mass weight the Hessian and initialize mass vector
      Hessian[step]=force_constants_mass_weight(Hessian[step]);

      //this is the call to internal coords used for MDtunnel
      //for the only rotor call it is redundant
      if(intcoord==1 && onlyrotors!=0){
	intcoord_hessian(step,Hessian[step],gradient[step]);
      } 

    //diagonalize mass weighted hessian

      if(onlyrotors!=0){
	double **L_mwc;
	L_mwc=diagonalize_hessian((int)(3*ATOMS),Hessian[step]);
	/*    
	fprintf(L_orig_save,"step=%d\n\n",step);
	for(j=1;j<(int)(dim+1);j++){
	  for(i=1;i<(int)(dim+1);i++){	
	    fprintf(L_orig_save,"%lf\t",L_mwc[j][i]);fflush(0);
	  }  
	  fprintf(L_orig_save,"\n");
	}
	fprintf(L_orig_save,"\n");
	*/
	for(i=0;i<(dim+1);i++){ free (L_mwc[i]);}free(L_mwc);
    
	frequencies=  (double *) malloc((int)(3*ATOMS)* sizeof(double));

    //lambda, the mass weighted Hessian eigenvalues, are converted in frequencies
    //the 6 lowest frequencies correspond to translations and internal rotations

	frequencies=calc_freq(lambda,frequencies);
  
	free(lambda);

	//	fprintf(freq_orig_results, "%d\t",step);
	//	for(i=0;i<dim;i++){
	//	  fprintf(freq_orig_results, "%lf\t",frequencies[i]);
	//	}
	//	fprintf(freq_orig_results, "\n");
      }

    //    here we have two alternatives RPH of MHA or Rotor Projection as implemented by Green et al. in Cantherm
    // the code terminates if we use the Rotor projection option

      if(onlyrotors==0){
	determine_top_atoms(step);
	if(intcoord==1){
	  intcoord_hessian(step,Hessian[step],gradient[step]);
	  for(i=0;i<(dim);i++){ free (Lcartint[i]);}free(Lcartint);
	  free(Freqint);
	  free(lambdaint);
	  projector_matrix_Rot(step);
	} else {
	  projector_matrix_Rot(step);
	} 
	//	fclose(anim_freq);
	//	fclose(L_orig_save);
	//	fclose(freq_results);
	//	fclose(freq_orig_results);

	//	fclose(Cart_displ);
	//	fclose(Km_res);
	exit(0);
      } 

    //generation of projector matrix P to project out 7 dof: 3 transl + 3 rotation + reaction coordinate 
    //P used as Project=I-P as defined by MHA in equation 1.5b
    //Miller,Handy, and Adams J. Chem. Phys. vol.72 p90 (1980) 
    

      double **Project;
      Project=(double **) malloc( dim*sizeof(double *) );
      for(j=0; j < dim; j++) {
	Project[j]= (double *) malloc( dim*sizeof(double) );
      } 

      //      printf("Bef proj mat %d\n",step);

      Project=projector_matrix(Project,step);
    
      //      printf("Aft proj mat %d\n",step);

    // NB from now on coordinates are mass weighted

      double **force_constants_int, **temp1;
  
    //calculation of projected Hessian (mass weighted force constant matrix)

      temp1=prod_mat(Hessian[step],Project,dim);
      force_constants_int=prod_mat(Project,temp1,dim);

    // force_constants_int is the projected force constant matrix KP
    // diagonalization of projected Hessian

      L_int=diagonalize_hessian(dim,force_constants_int);
  
    // save Hessian

      // Cartesian Coordinates
      for(j=0;j<(int)(dim);j++){
	for(i=0;i<(int)(dim);i++){
	  L_int_save[step][i][j]=L_int[i+1][j+1];      
	}
      } 
      //if internal coordinates then a part of the Eigenvector is replaced for 3N-6 coordinates
      //with the internal Eigenvectors
      //NB there is no check that the order is maintained

      if(intcoord==1){
	for(j=0;j<(int)(dim-6);j++){
	  for(i=0;i<(int)(dim);i++){
	    L_int_save[step][i][j]=Lcartint[i][j];      
	  }
	} 
      for(i=0;i<(dim);i++){ free (Lcartint[i]);}free(Lcartint);
      } 

      if((Cart_displ=fopen("Cart_displ.txt","w"))==NULL) {
	printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","Cart_displ.txt");
	exit(1);
      }

      save_displ(lambda,L_int_save[step],0,3*ATOMS-6,step,Cart_displ);
      if(step==saddlep){
	if(redmu<0.){
	  redmu=redmass_TS;
	  printf("setting as reduced mass of TS %6.2f \n",redmass_TS);
	}
      }

      fclose(Cart_displ);
    //compute and print cartesian displacements

      double **cart_displ;
      double **mass_matr;
      double **temp2;

      mass_matr=(double **) malloc( dim*sizeof(double *) );
      for(j=0; j < dim; j++) {
	mass_matr[j]= (double *) malloc( dim*sizeof(double) );
      } 

      for(j=0;j<(int)(dim);j++){
	for(i=0;i<(int)(dim);i++){
	  mass_matr[i][j]=0.;
	  if(i==j){
	    mass_matr[i][j]=1/sqrt(mass[i]);
	  }
	}
      } 

      temp2=prod_mat(Project,L_int_save[step],dim);
    //temp2=prod_mat(L_int_save[step],Project,dim);
    //    cart_displ=prod_mat(mass_matr,temp2,dim);
      cart_displ=prod_mat(temp2,mass_matr,dim);
  
      //      if(onlyrotors!=0){
      /*
      fprintf(Cart_displ,"step=%d\n\n",step);
      for(j=0;j<(int)(dim);j++){
	for(i=0;i<(int)(dim);i++){	
	    //	fprintf(Cart_displ,"%lf\t",cart_displ[j][i]);fflush(0);
	    //	    	    fprintf(Cart_displ,"%lf\t",L_int_save[step][j][i]);
	    //	    fprintf(Cart_displ,"%8.3g ",temp2[j][i]);fflush(0);
	  fprintf(Cart_displ,"%8.3g ",cart_displ[j][i]);fflush(0);
	}  
	fprintf(Cart_displ,"\n");
      }
      fprintf(Cart_displ,"\n");
      */
	//      }


    //evaluation of vibr freqs. 7 should be zero or near 0 
    //( 3 rotational + 3 translation + projected reaction path)
    //they are excluded from the calculation of VaG

      frequencies=calc_freq(lambda,frequencies);

    //    fprintf(freq_results,"step\tfreqs\n");
      
      //      fprintf(freq_results, "%d\t",step);
      for(i=0;i<dim;i++){     
	frequencies_save[step][i]=frequencies[i];
	//	fprintf(freq_results, "%lf\t",frequencies_save[step][i]);  
      }
      //      fprintf(freq_results, "\n");

      if(intcoord==1){
	for(i=0;i<dim-6;i++){     
	  frequencies_save[step][i]=Freqint[i];
	  //	  if(step==MAXSTEP2TS/2){
	  //  printf("freq_check %d %lf\n",i,frequencies_save[step][i]);   
	  // }
	}
	/*
	for(i=0;i<dim;i++){     
	  fprintf(freq_results, "%lf\t",Freqint[i]);  
	}
	*/
	//	if(step==MAXSTEP2TS/2){
	//	  exit(0);
	//	}
	free(Freqint);
	free(lambdaint);
      }


      free(lambda);

      //      fprintf(freq_results, "\n");


      if(MAXSTEP==1){


	double *start_coord,*arr_pcoord,*arr_ncoord;
	start_coord=  (double *) malloc((int)(3*ATOMS)* sizeof(double));
	arr_pcoord=  (double *) malloc((int)(3*ATOMS)* sizeof(double));
	arr_ncoord=  (double *) malloc((int)(3*ATOMS)* sizeof(double));
     
	int ij=0;
	for (i=0;i<ATOMS;i++){
	  start_coord[ij]=atoms_data[i].ext_coord_x[step]/bohr/sqrt(mass[ij]);
	  start_coord[ij+1]=atoms_data[i].ext_coord_y[step]/bohr/sqrt(mass[ij+1]);
	  start_coord[ij+2]=atoms_data[i].ext_coord_z[step]/bohr/sqrt(mass[ij+2]);
	  
	  arr_pcoord[ij]=atoms_data[i].ext_coord_x[step]/bohr/sqrt(mass[ij])+cart_displ[ij][an_freq]/bohr;
	  arr_pcoord[ij+1]=atoms_data[i].ext_coord_y[step]/bohr/sqrt(mass[ij+1])+cart_displ[ij+1][an_freq]/bohr;
	  arr_pcoord[ij+2]=atoms_data[i].ext_coord_z[step]/bohr/sqrt(mass[ij+2])+cart_displ[ij+2][an_freq]/bohr;
	  
	  arr_ncoord[ij]=atoms_data[i].ext_coord_x[step]/bohr/sqrt(mass[ij])-cart_displ[ij][an_freq]/bohr;
	  arr_ncoord[ij+1]=atoms_data[i].ext_coord_y[step]/bohr/sqrt(mass[ij+1])-cart_displ[ij+1][an_freq]/bohr;
	  arr_ncoord[ij+2]=atoms_data[i].ext_coord_z[step]/bohr/sqrt(mass[ij+2])-cart_displ[ij+2][an_freq]/bohr;

	  ij=ij+3;
	}
	/*
	for(i=0;i<100;i++) {
	  fprintf(anim_freq,"    %u  \n",ATOMS);
	  fprintf(anim_freq,"step_number%d\n",i);
	  ij=0;
	  for(j=0;j<ATOMS;j++){
	    fprintf(anim_freq,"%s\t%lf\t%lf\t%lf\n",atoms_data[j].atom_name,start_coord[ij],start_coord[ij+1],start_coord[ij+2]);
	    ij=ij+3;
	  }    
	  fprintf(anim_freq,"    %u  \n",ATOMS);
	  fprintf(anim_freq,"step_pos_displ_number%d\n",i);
	  ij=0;
	  for(j=0;j<ATOMS;j++){
	    fprintf(anim_freq,"%s\t%lf\t%lf\t%lf\n",atoms_data[j].atom_name,arr_pcoord[ij],arr_pcoord[ij+1],arr_pcoord[ij+2]);
	    ij=ij+3;
	  }    	
	  fprintf(anim_freq,"    %u  \n",ATOMS);
	  fprintf(anim_freq,"step_number%d\n",i);
	  ij=0;
	  for(j=0;j<ATOMS;j++){
	    fprintf(anim_freq,"%s\t%lf\t%lf\t%lf\n",atoms_data[j].atom_name,start_coord[ij],start_coord[ij+1],start_coord[ij+2]);
	    ij=ij+3;
	  }    
	  fprintf(anim_freq,"    %u  \n",ATOMS);
	  fprintf(anim_freq,"step_neg_displ_number%d\n",i);
	  ij=0;
	  for(j=0;j<ATOMS;j++){
	    fprintf(anim_freq,"%s\t%lf\t%lf\t%lf\n",atoms_data[j].atom_name,arr_ncoord[ij],arr_ncoord[ij+1],arr_ncoord[ij+2]);
	    ij=ij+3;
	  }    	
	}    
	*/
	free(start_coord);
	free(arr_pcoord);
	free(arr_ncoord);
	//	fclose(anim_freq);
      }

    //free vectors

      for(i=0;i<(dim+1);i++){ free (L_int[i]); } free(L_int);
      for(i=0;i<(dim);i++){ free (Project[i]);}free(Project);
      for(i=0;i<(dim);i++){ free (force_constants_int[i]);}free(force_constants_int);
      for(i=0;i<(dim);i++){ free (temp1[i]);}free(temp1);
      for(i=0;i<(dim);i++){ free (temp2[i]);}free(temp2);
      for(i=0;i<(dim);i++){ free (cart_displ[i]);}free(cart_displ);
      for(i=0;i<(dim);i++){ free (mass_matr[i]);}free(mass_matr);
      free(frequencies);
   
    }




    //    if(onlyrotors!=0){
      //      fclose(L_orig_save);
      //      fclose(freq_results);
      //      fclose(freq_orig_results);
      //      fclose(Cart_displ);
    //    }

  //end of cycle over different structures
  //now all projected frequencies are available


  // print L matrix for n steps in file L_save.txt
    /*
    if((L_save=fopen("L_save.txt","w"))==NULL) {
      printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","L_save.txt");
      exit(1);
    }

    for(step=0;step<MAXSTEP;step++){
      fprintf(L_save,"step=%d\n\n",step);
      for(j=0;j<(int)(dim);j++){
	for(i=0;i<(int)(dim);i++){
	  fprintf(L_save,"%lf\t",L_int_save[step][j][i]);
	}  
	fprintf(L_save,"\n");
      }
      fprintf(L_save,"\n");
    }
    fclose(L_save);
    */

    if(MAXSTEP == 1){
      printf(" single point analysis\n");
      printf(" correct end of calculations \n");
      exit(1);
    }

  //  interpolate_freqs(frequencies_save);


  //Beginning of calculations for RPH/SCT theories


 //********************************
 //Calculation of reduced masses
 //********************************

  //Evaluation of Hessian derivative

  // Hessian interpolation in the saddle point region of (+/- 9 points)

  //  smooth_hess();

  // Hessian gradient (finite differences)

    grad_hess();
  
  //internal checks - not particularly influent on final results  

  /*
  for(step=0;step<MAXSTEP;step++){
    L_int_save=check_L_sign(step, L_int_save);
  }

  for(step=1;step<MAXSTEP;step++){
    check_L_order(step, L_int_save);
  }
 
  */
  //calcolo BKf l'ultimo vettore della nonadiabatic coupling matrix Bkl per ogni punto della coordianta di reazione per cui vale  Km[s]=-BkF[s]


    //    fprintf(Km_res,"BkF:\n");
    
    for(step=0;step<MAXSTEP;step++){
        
      BkF=(double *) malloc( dim*sizeof(double ) );
      for(i=0;i<dim;i++) BkF[i]=0.;

      BkF=calc_BkF(dim,step,L_int_save[step],BkF);
    //    BkF=calc_BkF(dim-7,step,L_int_save[step],BkF);
    
      //      fprintf(Km_res,"%d\t",step);

      for(i=0;i<dim-7;i++){
	Km[step][i]= -BkF[i];  
      //cc test
	if(iminfreq!=0){
	  if(frequencies_save[step][i]<iminfreq) Km[step][i]=0.;
	}    
	//	fprintf(Km_res,"%le\t",Km[step][i]);      
      }    
      //      fprintf(Km_res,"\n");
      free(BkF);
    }

  // if not enough points before of after the saddle point the programs computes the VaG and stops

    if(((MAXSTEP-saddlep)<3)||(saddlep<3) ){
      calc_VaG(frequencies_save);
      printf("not enough points to coompute the tunneling contribution\n"); 
      printf("the program stops now\n"); 
      exit(0);
    }    
    
 //evaluation of K = magnitude of the curvature, it is expressed in 1/Bohr/uma^1/2

    for(step=0;step<MAXSTEP;step++){
      for(i=0;i<dim-7;i++){      
	K[step]=K[step]+ Km[step][i]*Km[step][i];
      }  
      K[step]=sqrt(K[step]);
      //      printf("step %d curvature %le\n",step,K[step]);        
    }
    //    exit(0);
 // interpolation in the saddle point region of Km (+/- N points)

 // interpolate_freqs(Km);

 //Calculate t_average = the maximum concave-side vibrational displacement along the curvature direction and k_average
 // k_ome is the sum of k^2i x omega^2i over i=1,N-7
 // frequencies in cm-1  are converted  in s-1 (freq x c_light) and then to omega (x 2pi), the eigenvalues of the hessian

    for(step=0;step<MAXSTEP;step++){
      k_ome[step]=0.;
      for(i=0;i<dim-7;i++){
	if(dsmethod!=1) k_ome[step]=k_ome[step]+(Km[step][i]*Km[step][i]*(frequencies_save[step][i]*c_light_cm_s*2.*pi)*(frequencies_save[step][i]*c_light_cm_s*2.*pi));
    //    if(dsmethod==1) k_ome[step]=k_ome[step]+(Km[step][i]*Km[step][i]*(frequencies_save[step][i]*c_light_cm_s*2.*pi)*(frequencies_save[step][i]*c_light_cm_s*2.*pi));
    // if activated must include square option 
    //        if(dsmethod==1) k_ome[step]=k_ome[step]+fabs(Km[step][i])*fabs(Km[step][i])*fabs((frequencies_save[step][i]*c_light_cm_s*2.*pi));
	if(dsmethod==1) k_ome[step]=k_ome[step]+fabs(Km[step][i])*sqrt(fabs((frequencies_save[step][i]*c_light_cm_s*2.*pi)));
    //    if(dsmethod==1) k_ome[step]=k_ome[step]+fabs(Km[step][i])*fabs(Km[step][i]);
    //    if(dsmethod==1)   k_ome[step]=k_ome[step]+((frequencies_save[step][i]*c_light_cm_s*2.*pi)*(frequencies_save[step][i]*c_light_cm_s*2.*pi));
    //    if(dsmethod==1)   k_ome[step]=k_ome[step]+(fabs(frequencies_save[step][i]*c_light_cm_s*2.*pi));
   }
 // it goes with the second option
      if(dsmethod==1) k_ome[step]=k_ome[step]*k_ome[step]*k_ome[step]*k_ome[step];
      printf("step %d t_average %le\n",step,k_ome[step]);        
    }

   
 // linear procedure to avoid oscillations in K. Active if Npoint different from 0
 // replace npoints of k_ome before and after tstate using linear average between maxstep/2-npoint and tstate

    int N_point=NpointsInt;

    for(i=0;i<N_point;i++){

   // left side
      k_ome[MAXSTEP2TS/2-N_point+i]= k_ome[MAXSTEP2TS/2-N_point]-(k_ome[MAXSTEP2TS/2-N_point]-k_ome[MAXSTEP2TS/2])/(N_point)*(i);
      K[MAXSTEP2TS/2-N_point+i]= K[MAXSTEP2TS/2-N_point]-(K[MAXSTEP2TS/2-N_point]-K[MAXSTEP2TS/2])/(N_point)*(i);

   // right side
      k_ome[MAXSTEP2TS/2+i]= k_ome[MAXSTEP2TS/2+N_point]-(k_ome[MAXSTEP2TS/2+N_point]-k_ome[MAXSTEP2TS/2])/(N_point)*(N_point-i);
      K[MAXSTEP2TS/2+i]= K[MAXSTEP2TS/2+N_point]-(K[MAXSTEP2TS/2+N_point]-K[MAXSTEP2TS/2])/(N_point)*(N_point-i);

    }

 // evaluate t_average and a_average
 // t has the dimensions of a length, it is calculated in meters and converted to Bohr

    for(step=0;step<MAXSTEP;step++){
      t[step]=sqrt(K[step]*h_bar/(a_mass*redmu))/sqrt(sqrt(k_ome[step]))/(a_bohr);
   // using this approximation we do not "normalize" the curvature by K, as it is not requested in RPH
      if(dsmethod==1)   t[step]=sqrt(h_bar/(a_mass*redmu))/sqrt(sqrt(k_ome[step]))/(a_bohr)/K[step];
   // if(dsmethod==1)   t[step]=sqrt(h_bar/(a_mass*redmu))/sqrt(k_ome[step])/(a_bohr);
      a[step]=fabs(K[step]*t[step]);
   //     if(dsmethod==1)    a[step]=fabs(K[step]);
   //fabs necessary since the product must be positive (concave side) 
   //   a[step]=(K[step]*t[step]);
      printf("step %d newcurvature %le\n",step,K[step]);        
      printf("step %d tparam %le\n",step,t[step]);        
      printf("step %d aparam %le\n",step,a[step]);        
    }

 //calculate  d(t_average)/ds

 
    double *t_spline;
    double *x_spline;
    double temp0=0.;
    double temp3=0.;
    double temp1=0, temp2=0;

    t_spline=(double *) malloc( MAXSTEP*sizeof(double ) );
    x_spline=(double *) malloc( MAXSTEP*sizeof(double ) );
    for(i=0;i<MAXSTEP;i++)x_spline[i]=i;
    t_spline=spline(x_spline,t ,MAXSTEP-1, 1e30, 1e30,t_spline );
 

    deltas=Rx_coord[1]-Rx_coord[0];
    t_dev[0]=(t[1]-t[0])/(deltas);
    deltas=Rx_coord[2]-Rx_coord[0];
    t_dev[1]=(t[2]-t[0])/(deltas);
    deltas=Rx_coord[MAXSTEP-1]-Rx_coord[MAXSTEP-3];
    t_dev[MAXSTEP-2]=(t[MAXSTEP-1]-t[MAXSTEP-3])/(deltas);
    deltas=Rx_coord[MAXSTEP-1]-Rx_coord[MAXSTEP-2];
    t_dev[MAXSTEP-1]=(t[MAXSTEP-1]-t[MAXSTEP-2])/(deltas);

 // finite diff from  4 points estimation

    for(step=2;step<MAXSTEP-2;step++){ 
      temp0=splint(x_spline,t ,t_spline , MAXSTEP-1, step-2, temp0);
      temp1=splint(x_spline,t ,t_spline , MAXSTEP-1, step-1, temp1);
      temp2=splint(x_spline,t ,t_spline , MAXSTEP-1, step+1, temp2);
      temp3=splint(x_spline,t ,t_spline , MAXSTEP-1, step+2, temp3);
      deltas=Rx_coord[step-2]-Rx_coord[step+2];
      t_dev[step]=(-temp3+8*temp2-8*temp1+temp0)/(3*deltas);
   //   t_dev[step]=(-temp3+8*temp2-8*temp1+temp0)/(3*deltas);
   //   deltas=Rx_coord[step-1]-Rx_coord[step+1];
   //   t_dev[step]=(temp2-temp1)/(deltas);
   //   t_dev[step]=(-t[step+2]+8*t[step+1]-8*t[step-1]+t[step-2])/(3*deltas);
   //   t_dev[step]=(t[step]-t[step-1])/(Rx_coord[step]-Rx_coord[step-1]);
      printf("step is %d tdev %le\n",step,t_dev[step]);        
      if(t_dev[step]>Maxtdev)t_dev[step]=Maxtdev;
      if(t_dev[step]<-Maxtdev)t_dev[step]=-Maxtdev;
    }

    free(t_spline);
    free(x_spline);
  
    double *mueff_mu1;
    mueff_mu1=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );


    //    if((mueff_results=fopen("mueff.txt","w"))==NULL) {
    //      printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","mueff.txt");
    //      exit(1);
    //    }

    //    fprintf(mueff_results," step pow(k_ome[step],-0.25) t[step]      t_dev[step]      a[step]         K[step]        mueff_mu1[step]  mueff_mu[step]\n");  

    for(step=0;step<MAXSTEP;step++){


   // SCT  approximation, mu_eff is expressed in atomic units
      mueff_mu1[step]=-2.*a[step]-(a[step]*a[step])+(t_dev[step]*t_dev[step]);
      mueff_mu[step]=exp(mueff_mu1[step])*massunit;

   // cc alternative to SCT  be used in alternative to exponential: shows spurious behavior if a>1
   // with these modifications the theory is equal to the simplest level of the RPH of Miller
   // check also above for the different weighting of t
      if(dsmethod==1)   t_dev[step]=0.;
      if(dsmethod==1) mueff_mu[step]=(pow(1-a[step],2)+pow(t_dev[step],2))*massunit;

      if(mueff_mu[step]>massunit){ mueff_mu[step]=massunit;}
      printf("step %d mueff %le\n",step,mueff_mu[step]);        

    }

 // write output 
 
    for(step=0;step<MAXSTEP;step++){
      if (N_point != 0){
	if(step == (MAXSTEP2TS/2-N_point+1)){mueff_mu[step]=(mueff_mu[step-2]+mueff_mu[step+2])/2.;}
	if(step == (MAXSTEP2TS/2-N_point)){mueff_mu[step]=(mueff_mu[step-2]+mueff_mu[step+2])/2.;}
	if(step == (MAXSTEP2TS/2+N_point)){mueff_mu[step]=(mueff_mu[step-2]+mueff_mu[step+2])/2.;}
	if(step == (MAXSTEP2TS/2+N_point-1)){mueff_mu[step]=(mueff_mu[step-2]+mueff_mu[step+2])/2.;}
      }
      //      fprintf(mueff_results,"%d\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",step,pow(k_ome[step],-0.25),t[step],t_dev[step],a[step],K[step],mueff_mu1[step],mueff_mu[step]);  
    }
 

//    fclose(mueff_results);

    free(K);free(a);free(t);free(t_dev);free(k_ome);

    for(i=0;i<MAXSTEP;i++){free(Km[i]);} free(Km);
  }


 // end of evaluation of reduced mass. 

 //********************************
 //Calculation of tunneling factors
 //********************************


 // evaluation of vibrationally adiabatic potential energies
  int Emax=0.;

  if(onlyrotors==0||onlyrotors==1){
     Emax=calc_VaG(frequencies_save);
  }
  if(onlyrotors==2){
     read_VaG_mueff();
  }

 int Emaxvar=Energy[MAXSTEP2TS/2];
 if(isctvtst==1){
   for(i=0;i<MAXSTEP;i++){
     if(Energy[i]>Emaxvar) Emaxvar=Energy[i];
   } 
   Emax=Emaxvar;
 } 
   printf(" Emax is : %d \n",Emax);  

 // evaluation of transmission coefficient using  imaginary action integral
 
 double dE=1,den; 
 double inv;
 double En;
 double *Transmission_coeff, *Transmission_coeff_spline ;

 Transmission_coeff=(double *) malloc( E_dim*sizeof(double ) );
 
 Transmission_coeff=calc_transmission_coeff(mueff_mu,Emax,Transmission_coeff);

 double *x,y=0;
 x=(double *) malloc( E_dim*sizeof(double ) );

  
 Transmission_coeff_spline=(double *) malloc( E_dim*sizeof(double ) ); 
 
 for(j=0;j<E_dim;j++){   
   x[j]=j;
 }

 Transmission_coeff_spline=spline(x, Transmission_coeff ,E_dim-1, 1e30, 1e30,Transmission_coeff_spline );
 
 // output file for E resolved  transmission coefficients
 /*
 if((P_E=fopen("Transmission_coefficient.txt","w"))==NULL) {
   printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","Transmission_coefficient.txt");
   exit(1);
 }
 */
 // fprintf(P_E,"E[cm-1]\t P(E)\n");
 for(En=0;En<E_dim;En=En+dE){
   y=splint(x,Transmission_coeff ,Transmission_coeff_spline , E_dim-1, En, y);
   //   fprintf(P_E, "%lf\t %le\n",En, y);
   //   fprintf(P_E, "%lf\t %le\n",En, Transmission_coeff[E]);
 }
 

 /*
 if((P_E=fopen("Transmission_coefficient.txt","w"))==NULL) {
   printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","Transmission_coefficient.txt");
   exit(1);
 }

 for(i=0;i<E_dim;i=i+1){
   fprintf(P_E, "%d %le\n",i, Transmission_coeff[i]);
 }
 */
 //  fclose(P_E);

 // exit(0);

 //Evaluation of tunneling correction factor kSCT

  dE=0.1;
  // dE=10.0;
 double Tfinal;
 double var_corr;
 // double den2;

 Tfinal=Tinitial+(Tsteps+1.)*Tdelta;

 for(T=Tinitial;T<Tfinal;T=T+Tdelta){   
   var_corr=calc_var_corr(dim,T,frequencies_save);

   //   if((isctvtst==1) && (Energy[stepmin]>Energy[MAXSTEP/2])) Emaxvar=Energy[stepmin];

   double kT,num=0.;
   den=0.;
   inv= h_planck*c_light_cm_s/boltzmann/T; 
   for(En=0;En<E_dim;En=En+dE){
     y=splint(x,Transmission_coeff ,Transmission_coeff_spline , E_dim-1, En, y);
     num=num+y*exp(-En*inv);
   } 
   for(En=Emaxvar;En<E_dim;En=En+dE){
     den=den+ exp(-En*inv);
   }
   //   den2=exp(-Emaxvar*h_planck*c_light_cm_s/boltzmann/T)*boltzmann*T/h_planck/c_light_cm_s/dE;
   kT=num/den;
   printf("T=%.1f  kSCT=%le/%le=%lf var corr %lf stepmin %d Emax %d \n",T,num,den,kT,var_corr,stepmin,Emaxvar);  
 }

 if(onlyrotors==2){
   return 0;
 }

 //libero la memoria allocata

 free(x);
 free(Transmission_coeff_spline);
 free(Transmission_coeff);
 free(mueff_mu);

 for(i=0;i<MAXSTEP;i++){for(j=0;j<3*ATOMS;j++){free(L_int_save[i][j]);} free(L_int_save[i]);} free(L_int_save);
 for(i=0;i<MAXSTEP;i++){free(frequencies_save[i]);}free(frequencies_save);
 for(i=0;i<MAXSTEP;i++){free(gradient[i]);}free(gradient);
 // for(i=0;i<3*ATOMS;i++){free(gradient_derivative[i]);}free(gradient_derivative);
 for(i=0;i<MAXSTEP;i++){for(j=0;j<3*ATOMS;j++){free(Hessian[i][j]);} free(Hessian[i]);} free(Hessian);
 free(Energy);
 free(Rx_coord);
 for(i=0;i<MAXSTEP;i++){for(j=0;j<4;j++){free(I[i][j]);} free(I[i]);} free(I);
 for(i=0;i<MAXSTEP;i++){for(j=0;j<4;j++){free(v[i][j]);} free(v[i]);} free(v);
 for(i=0;i<MAXSTEP;i++){free(eigen[i]);}free(eigen);
 free(CM_x);
 free(CM_y);
 free(CM_z);
 free(incoox);
 free(incooy);
 free(incooz);
 for(i=0;i<ATOMS;i++){
   free(atoms_data[i].ext_coord_x);
   free(atoms_data[i].ext_coord_y);
   free(atoms_data[i].ext_coord_z);

   free(atoms_data[i].grad_x);
   free(atoms_data[i].grad_y);
   free(atoms_data[i].grad_z);
 }
 free(atoms_data);
 

 // fclose(freq_results1);
 // fclose(Km_res);
 // fclose(P_E);


 //***************************************************************************************************************/ 

  return 0;
}


void allocation( int dim ) {

  

  int i,j,k;

  for(i=0;i<ATOMS;i++) {
    atoms_data[i].ext_coord_x =(double*) malloc(dim * sizeof(double));
    atoms_data[i].ext_coord_y =(double*) malloc(dim * sizeof(double));
    atoms_data[i].ext_coord_z =(double*) malloc(dim * sizeof(double));

    atoms_data[i].grad_x =(double*) malloc(dim * sizeof(double));
    atoms_data[i].grad_y =(double*) malloc(dim * sizeof(double));
    atoms_data[i].grad_z =(double*) malloc(dim * sizeof(double));

  }

  gradient=(double**)malloc((dim) * sizeof(double*));
  for(i=0;i<dim;i++) {
    gradient[i] = (double*)malloc((int)(3*ATOMS)*sizeof(double));
  }

  
  Hessian = (double***)malloc((dim) * sizeof(double**));
  for(i=0;i<dim;i++) {
    Hessian[i] = (double**)malloc((int)(3*ATOMS) * sizeof(double*));
  }

  for(i=0;i<dim;i++) {
    for(j=0;j<((int)(3*ATOMS));j++) {
      Hessian[i][j] = (double*)malloc((int)(3*ATOMS)*sizeof(double));
    }
  }

  for(k=0;k<(int)(dim);k++) {
    for(i=0;i<((int)(3*ATOMS));i++) {
      for(j=0;j<((int)(3*ATOMS));j++) {
	
	Hessian[k][i][j]=0.0;
	//	printf("%lf\t",Hessian[i][j][k]);
      }
    }
  }

  CM_x= (double*) malloc( dim*sizeof(double) );
  CM_y= (double*) malloc( dim*sizeof(double) );
  CM_z= (double*) malloc( dim*sizeof(double) );

   I = (double***)malloc((int)(dim) * sizeof(double**));
  
  for(i=0;i<dim;i++) {
    I[i] = (double**)malloc((int)(4) * sizeof(double*));
  }

  for(i=0;i<dim;i++) {
    for(j=0;j<4;j++) {
      I[i][j] = (double*)malloc((int)(4)*sizeof(double));
    }
  }

   v = (double***)malloc((int)(dim) * sizeof(double**));
  
  for(i=0;i<dim;i++) {
    v[i] = (double**)malloc((int)(4) * sizeof(double*));
  }

  for(i=0;i<dim;i++) {
    for(j=0;j<4;j++) {
      v[i][j] = (double*)malloc((int)(4)*sizeof(double));
    }
  }

  eigen=(double**)malloc((dim) * sizeof(double*));
  for(i=0;i<dim;i++) {
    eigen[i] = (double*)malloc((int)(4)*sizeof(double));
  }
 
  Energy= (double*) malloc( MAXSTEP*sizeof(double) );

 
  Rx_coord= (double*) malloc( MAXSTEP*sizeof(double) );

}



void read_inputfile () {

  char buffer[200];
  char filename[20];
  int j,i,ij;
  int natread;
  int s=0;
  FILE *fp;


  //open file RPHt_input_data.dat

  sprintf(filename,"RPHt_input_data.dat");

  if((fp=fopen(filename,"r"))==NULL){
    printf("main: IMPOSSIBILE APRIRE IL FILE: %s\nProgram now exits\n",filename);
    exit(1);
  }



  // Leggo il numero di atomi
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&ATOMS);
    //    printf("ATOMS NUMBER:      %d\n", ATOMS);
    if(ATOMS<0) {
      printf("the number of atoms cannot be negative\n"
	     "Change the value and restart the program\n"
	     "the program will be stopped now\n");
      exit(0);
    }
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }
  
 // Leggo l'energia SCF+ZPE dei reagenti
 //Updated, it is the energy barrier in kcal/mol wrt reacts without ZPE
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&EBarr_kcalmol);
    //    printf("Energy of the reactants:      %lf\n", Ereactants);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }


  // Read Initial Temperature for T scan
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&Tinitial);
    //    printf("INITIAL TEMPERATURE:      %lf\n", Tinitial);
    if(Tinitial<0) {
      printf("the temperature cannot be negative\n"
	     "Change the value and restart the program\n"
	     "the program will be stopped now\n");
      exit(0);
    }
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read number of Temperature steps
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&Tsteps);
    //    printf("TEMPERATURE step:      %lf\n", Tsteps);
    if(Tsteps<0) {
      printf("the temperature steps must be positive\n"
	     "Change the value and restart the program\n"
	     "the program will be stopped now\n");
      exit(0);
    }
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read Temperature increment
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&Tdelta);
    //    printf("TEMPERATURE increment:      %lf\n", Tdelta);
    if(Tdelta<0) {
      printf("the temperature increment must be positive\n"
	     "Change the value and restart the program\n"
	     "the program will be stopped now\n");
      exit(0);
    }
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }


  // Read Delta Energy for reactant side: energy to shift the maximum and rescale linearly
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&Delta_Energy_rea);
    //    printf("Delta Energy:      %lf\n", Delta_Energy);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read Delta Energy for product side: energy to shift the maximum and rescale linearly
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&Delta_Energy_pro);
    //    printf("Delta Energy:      %lf\n", Delta_Energy);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }


  // Read Maxstep: total number of points of the PES
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&MAXSTEP);
    //    printf("Delta Energy:      %d\n", MAXSTEP);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read NpointsInt: number of points around the saddle point to fit
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&NpointsInt);
    //    printf("Npoints interpolated:      %d\n", NpointsInt);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read Maxtdev: parameter to limit the maximum value of change of tdev
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&Maxtdev);
    //    printf("Maxtedev:      %lf\n", Maxtdev);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read rearrange: parameter to decide if rearrange with absolute value of eigenvalues
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&rearrange);
    //    printf("Rearrange:      %d\n", rearrange);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read saddlep: point at which we have the TS
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&saddlep);
    //    printf("Saddle Point :      %d\n", saddlep);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }
  //initialize saddle point parameter (saddlepoint*2)
  MAXSTEP2TS=saddlep*2+1;
  //intialize

  // Read coordinates type
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&intcoord);
    //    printf("Cartesian0/Iternal1 coord:      %d\n", intcoord);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read index for variation/sct or sct only theory
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&isctvtst);
    //    printf("using variational tst(1):      %d\n", isctvtst);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read index for ZCT or sct  theory
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&izct);
    //    printf("using zero curvature(1):      %d\n", izct);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read reduced mass (multiple of atomic mass)
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&redmu);
    printf("using reduced mass:      %lf\n", redmu);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }


  // Read iminfreq (minimum frequency accepted in calculations)
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&iminfreq);
    //    printf("setting minimum frequency:      %lf\n", iminfreq);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }


  // Read frequency to animate
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&an_freq);
     printf("frequency to animate:      %d\n", an_freq);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

    // Read flag for rotor projection code
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&onlyrotors);
    //    printf("rotor flag:      %d\n", onlyrotors);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  // Read if project along reaction coordinate
  if(!feof(fp)) {
    fscanf(fp,"%s %lf",buffer,&proj_rea_coo);
    //    printf("projection along react coord:      %lf\n", proj_rea_coo);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

    // Read number of rotors
  if(!feof(fp)) {
    fscanf(fp,"%s %d",buffer,&numrotors);
    printf("numrotors:      %d\n", numrotors);
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  pivotA= (int *) malloc(numrotors * sizeof(int) );
  pivotB= (int *) malloc(numrotors * sizeof(int) );
  numatomsintopA= (int *) malloc(numrotors * sizeof(int) );

  atomsintopA =  (int **) malloc((int)(numrotors)* sizeof(int*));
  for(j=0; j < numrotors; j++) {
    atomsintopA[j]= (int *) malloc(ATOMS * sizeof(int) );
  }
  atomsintopB =  (int **) malloc((int)(numrotors)* sizeof(int*));
  for(j=0; j < numrotors; j++) {
    atomsintopB[j]= (int *) malloc(ATOMS * sizeof(int) );
  }

  for (i=0;i<numrotors;i++) {
    for (j=0;j<ATOMS;j++) {
      atomsintopA[i][j]=0;
      atomsintopB[i][j]=0;
    }
  }


  for (i=0;i<numrotors;i++) {

    // Read index of pivotA atom
    if(!feof(fp)) {
      fscanf(fp,"%s %d",buffer,&pivotA[i]);
      //          printf("pivotA:      %d\n", pivotA[i]);
    } else {
      printf("Data file error\nthe program will be stopped now\n");
      exit(0);
    }
 
    // Read index of pivotB atom
    if(!feof(fp)) {
      fscanf(fp,"%s %d",buffer,&pivotB[i]);
      //    printf("pivotB:      %d\n", pivotB);
    } else {
      printf("Data file error\nthe program will be stopped now\n");
      exit(0);
    }

    // Read number of atoms in group A
    if(!feof(fp)) {
      fscanf(fp,"%s %d",buffer,&numatomsintopA[i]);
      //      printf("numatomsintopA:      %d\n", numatomsintopA[i]);
    } else {
      printf("Data file error\nthe program will be stopped now\n");
      exit(0);
    }


    //  }
    
    
    // Read atoms in group A
    fscanf(fp,"%s ",buffer);  
    //    printf("buffer is:      %s\n", buffer);
    for(j=0;j<numatomsintopA[i];j++){
      if(!feof(fp)) {
	//	fscanf(fp,"%d",&atomsintopA[i][j]);
	fscanf(fp,"%d",&natread);
	//	printf("atomsintopA:      %d\n", natread);
	atomsintopA[i][natread-1]=1;
      } else {
	printf("Data file error\nthe program will be stopped now\n");
	exit(0);
      }
    }
    for(j=0;j<ATOMS;j++){
      if((atomsintopA[i][j]!=1) && (j!=pivotA[i]-1) && (j!=pivotB[i]-1)){
	atomsintopB[i][j]=1;
      }
    }
  }

  /*
  for (i=0;i<numrotors;i++) {
    for (j=0;j<ATOMS;j++) {
      printf("at A i % d j %d is:  %d\n",i,j,atomsintopA[i][j]);
    }
  }
  for (i=0;i<numrotors;i++) {
    for (j=0;j<ATOMS;j++) {
      printf("at B i % d j %d is:  %d\n",i,j,atomsintopB[i][j]);
    }
  }

  exit(0);
  */

  /*  
  if(!feof(fp)) {
    fscanf(fp,"%s",buffer);    
  } else {
    printf("Data file error: a line in the input  is missing\nthe program will be stopped now\n");
    exit(0);
  }
  */
  atoms_data = malloc(sizeof(atom_struct) * ATOMS); 
  allocation(MAXSTEP);

  //  printf("MAXSTEP is %d \n",MAXSTEP);        

  //Leggo le masse atomiche 
  // not anymore necessary
  /*
  int i1;
  for(j=0;j<ATOMS;j++){
    if(!feof(fp)) {
      fscanf(fp,"%s",&atoms_data[j].atom_name);  
      //      printf("Atom_name: %s\t",&atoms_data[j].atom_name);

      fscanf(fp,"%lf",&atoms_data[j].atomic_mass);
      //      printf("Atom_mass: %lf\n",atoms_data[j].atomic_mass);

      if(atoms_data[j].atomic_mass<0) {
	printf("atomic mass cannot be neagtive\n"
	       "Change the value and restart the program\n"
	       "the program will be stopped now\n");
	exit(0);
      }
    } else {
      printf("Data file error\nthe program will be stopped now\n");
      exit(0);
    }
  }
  */ 

  incoox =  (double **) malloc((int)(MAXSTEP)* sizeof(double*));
  for(j=0; j < MAXSTEP; j++) {
    incoox[j]= (double *) malloc(ATOMS * sizeof(double) );
  }
  incooy =  (double **) malloc((int)(MAXSTEP)* sizeof(double*));
  for(j=0; j < MAXSTEP; j++) {
    incooy[j]= (double *) malloc(ATOMS * sizeof(double) );
  }
  incooz =  (double **) malloc((int)(MAXSTEP)* sizeof(double*));
  for(j=0; j < MAXSTEP; j++) {
    incooz[j]= (double *) malloc(ATOMS * sizeof(double) );
  }

      
  if(onlyrotors!=2){
    while(s<MAXSTEP){

    // Leggo lo step
      if(!feof(fp)) {
      //      fscanf(fp,"%s %d %s",buffer,&step, buffer);
	fscanf(fp,"%s %d ",buffer,&step);
      //      printf("bufferread1:      %s\n", buffer);
      //      step=step-1;
	printf("STEP NUMBERread:      %d\n", step);
	step=s;
      //inverted numbering to be consistent with estoktp format
	//	step=MAXSTEP-1-s;

	printf("STEP NUMBER:      %d\n", step);
	if(step<0) {
	  printf("the number of step cannot be negative\n"
		 "Change the value and restart the program\n"
		 "the program will be stopped now\n");
	  exit(0);
	}
      } else {
	printf("Data file error\nthe program will be stopped now\n");
	exit(0);
      }

    //    printf("\nCoordinates: \n");

      fscanf(fp,"%s",buffer);  
    //    printf("bufferread2:      %s\n", buffer);
    

    // read coordinates
      int atomic_number;
      i=0; 
      for(j=0;j<ATOMS;j++){
	if(!feof(fp)) {
	
	  fscanf(fp,"%s",buffer);  	
	  fscanf(fp,"%d",&atomic_number);
	  fscanf(fp,"%s",buffer);  	
	  fscanf(fp,"%lf",&atoms_data[j].ext_coord_x[step]);
	  fscanf(fp,"%lf",&atoms_data[j].ext_coord_y[step]);
	  fscanf(fp,"%lf",&atoms_data[j].ext_coord_z[step]);
	//	printf("z coo:      %lf\n",atoms_data[j].ext_coord_z[step] );
	  incoox[step][j]=atoms_data[j].ext_coord_x[step];
	  incooy[step][j]=atoms_data[j].ext_coord_y[step];
	  incooz[step][j]=atoms_data[j].ext_coord_z[step];

	  if(s==0){
	
	    if (atomic_number==1){
	      atoms_data[j].atom_name="H";	  
	      atoms_data[j].at_numb=1;	  
	    // gaussian and IUPAC coincide
	      atoms_data[j].atomic_mass=1.00783;
	    }
	    else if(atomic_number==6){
	      atoms_data[j].atom_name="C";	  
	      atoms_data[j].at_numb=6;	  
	    // gaussian mass
	      atoms_data[j].atomic_mass=12.0;
	      	    // mass from IUPAC atomic weight of the elements 1979 Pure & Appl. Chem.., Vol.52, pp.2349—2384
	    //atoms_data[j].atomic_mass=12.0107;
	    }
	    else if (atomic_number==7){
	      atoms_data[j].atom_name="N";	  
	      atoms_data[j].at_numb=7;	  
	    //gaussian mass
	      atoms_data[j].atomic_mass=14.003074;
	    // mass from IUPAC atomic weight of the elements 1979 Pure & Appl. Chem.., Vol.52, pp.2349—2384
	    //atoms_data[j].atomic_mass=14.0067;
	    }
	    else if (atomic_number==8){
	      atoms_data[j].atom_name="O";	  
	      atoms_data[j].at_numb=8;	  
	    // gaussian mass
	      atoms_data[j].atomic_mass=15.99491;
	      // mass from IUPAC atomic weight of the elements 1979 Pure & Appl. Chem.., Vol.52, pp.2349—2384
	    //atoms_data[j].atomic_mass=15.9994;
	    }
	    else if (atomic_number==9){
	      atoms_data[j].atom_name="F";	  
	      atoms_data[j].at_numb=9;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=18.9984;
	    }
	    else if (atomic_number==10){
	      atoms_data[j].atom_name="Ne";	  
	      atoms_data[j].at_numb=10;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=19.9924391;
	    }
	    else if (atomic_number==13){
	      atoms_data[j].atom_name="Al";	  
	      atoms_data[j].at_numb=13;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=26.981538;
	    }
	    else if (atomic_number==14){
	      atoms_data[j].atom_name="Si";	  
	      atoms_data[j].at_numb=14;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=27.976928;
	    }
	    else if (atomic_number==15){
	      atoms_data[j].atom_name="P";	  
	      atoms_data[j].at_numb=15;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=30.9737634;
	    }
	    else if (atomic_number==16){
	      atoms_data[j].atom_name="S";	  
	      atoms_data[j].at_numb=16;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=31.9720718;
	    }
	    else if (atomic_number==17){
	      atoms_data[j].atom_name="Cl";	  
	      atoms_data[j].at_numb=17;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=34.968852729;
	    }
	    else if (atomic_number==18){
	      atoms_data[j].atom_name="Ar";	  
	      atoms_data[j].at_numb=18;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=39.9623831;
	    }
	    else if (atomic_number==31){
	      atoms_data[j].atom_name="Ga";	  
	      atoms_data[j].at_numb=31;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=68.92558;
	    }
	    else if (atomic_number==32){
	      atoms_data[j].atom_name="Ge";	  
	      atoms_data[j].at_numb=32;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=73.92117788;
	    }
	    else if (atomic_number==33){
	      atoms_data[j].atom_name="As";	  
	      atoms_data[j].at_numb=33;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=74.9215955;
	    }
	    else if (atomic_number==34){
	      atoms_data[j].atom_name="Se";	  
	      atoms_data[j].at_numb=34;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=79.91652;
	    }
	    else if (atomic_number==35){
	      atoms_data[j].atom_name="Br";	  
	      atoms_data[j].at_numb=35;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=78.9183361;
	    }
	    else if (atomic_number==36){
	      atoms_data[j].atom_name="Kr";	  
	      atoms_data[j].at_numb=36;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=83.911506;
	    }
	    else if (atomic_number==49){
	      atoms_data[j].atom_name="In";	  
	      atoms_data[j].at_numb=49;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=114.9041;
	    }
	    else if (atomic_number==50){
	      atoms_data[j].atom_name="Sn";	  
	      atoms_data[j].at_numb=50;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=117.9018;
	    }
	    else if (atomic_number==51){
	      atoms_data[j].atom_name="Sb";	  
	      atoms_data[j].at_numb=51;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=120.9038;
	    }
	    else if (atomic_number==52){
	      atoms_data[j].atom_name="Te";	  
	      atoms_data[j].at_numb=52;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=129.9067;
	    }
	    else if (atomic_number==53){
	      atoms_data[j].atom_name="I";	  
	      atoms_data[j].at_numb=53;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=126.9004;
	    }
	    else if (atomic_number==54){
	      atoms_data[j].atom_name="Xe";	  
	      atoms_data[j].at_numb=54;	  
	      // gaussian mass
	      atoms_data[j].atomic_mass=131.905429;
	    }
	    else {
	      printf("The mass of this element is not in the code\n");
	      printf("Add it and recompile\n");
	      printf("The atomic number is %d \n",atomic_number);
	      exit(0);
	    }
	  }
      

	//	printf("%lf\t",atoms_data[j].ext_coord_x[step]);
	//	printf("%lf\t",atoms_data[j].ext_coord_y[step]);
	//	printf("%lf\n",atoms_data[j].ext_coord_z[step]);
	  
	  i=i+3;
		
	} else {
	  printf("Data file error\nthe program will be stopped now\n");
	  exit(0);
	}
      }
    
      if(!feof(fp)) {
	fscanf(fp,"%s",buffer);
      
      } else {
	printf("Data file error\nthe program will be stopped now\n");
	exit(0);
      }
    
      //    printf("\nGradient: \n");

      // Read the gradient
      int i1;
      i=0;
      for(j=0;j<ATOMS;j++){
	if(!feof(fp)) {
	
	  fscanf(fp,"%d",&i1);  
	  fscanf(fp,"%d",&i1);
	  fscanf(fp,"%lf",&atoms_data[j].grad_x[step]);
	  fscanf(fp,"%lf",&atoms_data[j].grad_y[step]);
	  fscanf(fp,"%lf",&atoms_data[j].grad_z[step]);
	  
	// mass weight the gradient
	// the gradient is in Hartree/Bohr/uma
	
	  gradient[step][i]=atoms_data[j].grad_x[step]/sqrt(atoms_data[j].atomic_mass);
	  gradient[step][i+1]=atoms_data[j].grad_y[step]/sqrt(atoms_data[j].atomic_mass);
	  gradient[step][i+2]=atoms_data[j].grad_z[step]/sqrt(atoms_data[j].atomic_mass);
	  
	//	printf("%lf\t",gradient[step][i]);
	//	printf("%lf\t",gradient[step][i+1]);
	//	printf("%lf\n",gradient[step][i+2]);
	  i=i+3;
	
	}
      
      }
    
      if(!feof(fp)) {
	fscanf(fp,"%s",buffer);
      } else {
	printf("Data file error\nthe program will be stopped now\n");
	exit(0);
      }
        
    //Read Hessian in cartesian coordinates
    
    //    printf("\nHessian: \n");
    
      int dim=3*ATOMS;
    
      if(!feof(fp)) {
	int start=0;      
	for(i=0;i<(int)(3*ATOMS/5);i++){	
	  for(ij=0;ij<5;ij++){
	    fscanf(fp,"%s",buffer);	  
	  }	
	  for(ij=start;ij<3*ATOMS;ij++){
	    fscanf(fp,"%s",buffer);
	    for(j=start;(j<5+start && j<=ij);j++ ){
	      fscanf(fp,"%le  ",&Hessian[step][ij][j]);
	    }	
	  }
	  start=start +5;
	}
	int rest=0;
	rest=dim-(int)(dim/5)*5;
	if(rest!=0){
	  for(ij=0;ij<rest;ij++){
	    fscanf(fp,"%s",buffer);
	  }           
	  for(ij=start;ij<(int)(dim);ij++){	
	    fscanf(fp,"%s  ",buffer);
	    for(j=start;(j<rest+start && j<=ij);j++ ){
	      fscanf(fp,"%lf  ",&Hessian[step][ij][j]);
	    }
	  }
	}        
      } else {
	printf("Data file error\nthe program will be stopped now\n");
	exit(0);
      }    
      for(j=0;j<(int)(3*ATOMS);j++){            
	for(i=0;i<=j;i++){
	  //	printf("%lf\t",Hessian[step][j][i]);	
	  //	  Hessian[step][j][i]=Hessian[step][j][i]*1.1;	
	}
	//      printf("\n");
      }
      s=s+1;
    }  
  }  

  // end of cycle on steps

  fclose(fp);  
}

void read_VaG_mueff () {

  int i;
  char buffer[100];
  char filename[25];
  FILE *fp;

  sprintf(filename,"RPHt_step_en_mueff.dat");
  if((fp=fopen(filename,"r"))==NULL){
    printf("main: IMPOSSIBILE APRIRE IL FILE: %s\nProgram now exits\n",filename);
    exit(1);
  }

  if(!feof(fp)) {
    fscanf(fp,"%s %s %s ",buffer,buffer,buffer);
    for(i=0;i<MAXSTEP;i++){
      //inverted numbering to be consistent with estoktp format
      //      step=MAXSTEP-1-i;
      step=i;
      fscanf(fp,"%lf %lf %lf ",&Rx_coord[step],&Energy[step],&mueff_mu[step]);      
      printf("%lf %lf %lf \n",Rx_coord[step],Energy[step],mueff_mu[step]);      
    }
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }

  for(i=0;i<MAXSTEP;i++){
    Energy[i]=Energy[i]*349.75;
  }

  fclose(fp);



}


void read_pesrxfile () {

  int i,step;
  char buffer[100];
  char filename[20];
  FILE *fp;

  sprintf(filename,"RPHt_coord_en.dat");
  if((fp=fopen(filename,"r"))==NULL){
    printf("main: IMPOSSIBILE APRIRE IL FILE: %s\nProgram now exits\n",filename);
    exit(1);
  }

  if(!feof(fp)) {
    fscanf(fp,"%s %s %s %s %s",buffer,buffer,buffer,buffer,buffer);
    for(i=0;i<MAXSTEP;i++){
      //inverted numbering to be consistent with estoktp format
      //      step=MAXSTEP-1-i;
      step=i;
      fscanf(fp,"%s %lf %lf %s  %s ",buffer,&Rx_coord[step],&Energy[step],buffer,buffer);      
    }
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }
  fclose(fp);

  //reaction energy rescaling

      
  double deltar=Delta_Energy_rea; 
  double deltap=Delta_Energy_pro; 
  for(i=0;i<MAXSTEP;i++){
    if(i <= (MAXSTEP2TS-1)/2) Energy[i]=Energy[i]+deltar/627.5*((i-((MAXSTEP2TS-1)/2.))/((MAXSTEP2TS-1)/2.));
    //    if(i > (MAXSTEP-1)/2) Energy[i]=Energy[i]+deltap/627.5*((MAXSTEP-1)/2.-i)/((MAXSTEP-1)/2.);
    if(i > (MAXSTEP2TS-1)/2) Energy[i]=Energy[i]+deltap/627.5*((MAXSTEP2TS-1)/2.-i+1)/((MAXSTEP2TS-1)/2.-1);

    //    if(i <= MAXSTEP/2) Energy[i]=Energy[i]+deltar/627.5*(i/(MAXSTEP/2.));
    //    if(i > MAXSTEP/2) Energy[i]=Energy[i]+deltap/627.5*(MAXSTEP-i)/(MAXSTEP/2.);
  }
  
  
}



void read_pesfile () {

  int i;
  char buffer[100];
  char filename[20];
  FILE *fp;

  sprintf(filename,"RPHt_PES.dat");
  if((fp=fopen(filename,"r"))==NULL){
    printf("main: IMPOSSIBILE APRIRE IL FILE: %s\nProgram now exits\n",filename);
    exit(1);
  }

  if(!feof(fp)) {
    for(i=0;i<MAXSTEP;i++){
      fscanf(fp,"%s %s %s %s",buffer,buffer,buffer,buffer);
      fscanf(fp,"%lf %s %s %s %s ",&Energy[i],buffer,buffer,buffer,buffer);      
    }
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }
  fclose(fp);

  //reaction energy rescaling

      
  double delta=Delta_Energy_rea; 
  for(i=0;i<MAXSTEP;i++){
    if(i <= MAXSTEP2TS/2) Energy[i]=Energy[i]+delta/627.5*((i-(MAXSTEP2TS/2.))/(MAXSTEP2TS/2.));

    if(i > MAXSTEP2TS/2) Energy[i]=Energy[i]+delta/627.5*(MAXSTEP2TS/2.-i)/(MAXSTEP2TS/2.);
  }
  
  
}


void read_Bmat_Cmat () {

  int i,j,k;
  //  char buffer[100];
  char filename[20];
  FILE *fp;
  int browsc;
  int bcolumnsc;

  sprintf(filename,"bmat.dat");
  if((fp=fopen(filename,"r"))==NULL){
    printf("main: IMPOSSIBILE APRIRE IL FILE: %s\nProgram now exits\n",filename);
    exit(1);
  }
  fscanf(fp,"%d %d ",&brows,&bcolumns);

  // allocate bmat
  printf("B rows and columns %d %d\n",brows,bcolumns);

  Bmat =  (double **) malloc((int)(brows)* sizeof(double*));
  for(step=0; step<brows; step++) { Bmat[step]=  (double *) malloc((int)(bcolumns)* sizeof(double));}

  Cmat =  (double ***) malloc((int)(brows)* sizeof(double**));
  for(step=0; step<brows; step++) { Cmat[step]=  (double **) malloc((int)(bcolumns)* sizeof(double*));}
  for(step=0; step<brows; step++){ for(j=0; j <(int)(bcolumns); j++) { Cmat[step][j]=  (double *) malloc((int)(bcolumns)* sizeof(double));}}


  if(!feof(fp)) {
    for(i=0;i<brows;i++){
      for(j=0;j<bcolumns;j++){
	fscanf(fp,"%lf",&Bmat[i][j]);
      }
    }
  } else {
    printf("Data file error - Bmat\nthe program will be stopped now\n");
    exit(0);
  }

fclose(fp);
  //  printf("last value %lf\n",Bmat[brows-1][bcolumns-1]);
  //  exit(0);

  /*
  L_int_save=(double ***) malloc( MAXSTEP*sizeof(double **) );
  for(step=0;step<MAXSTEP;step++){ L_int_save[step]=(double **) malloc( (int)(dim)*sizeof(double *) );}
  for(step=0;step<MAXSTEP;step++){ for(j=0; j <(int)(dim); j++) { L_int_save[step][j]= (double*) malloc((int)(dim)*sizeof(double) );} }

  Hessian_deriv=(double ***) malloc( MAXSTEP*sizeof(double **) );for(step=0;step<MAXSTEP;step++){ Hessian_deriv[step]=(double **) malloc( (int)(dim)*sizeof(double *) );}

  for(step=0;step<MAXSTEP;step++){ for(j=0; j <(int)(dim); j++) { Hessian_deriv[step][j]= (double*) malloc((int)(dim)*sizeof(double) );} }

  */

  sprintf(filename,"cmat.dat");
  if((fp=fopen(filename,"r"))==NULL){
    printf("main: IMPOSSIBILE APRIRE IL FILE: %s\nProgram now exits\n",filename);
    exit(1);
  }

 fscanf(fp,"%d %d ",&browsc,&bcolumnsc);
 if(browsc!=brows || bcolumnsc != bcolumns){
   printf("main: Bmat and Cmat not compatible: rows %d %d\n",brows,browsc);
   printf("main: Bmat and Cmat not compatible: columns %d %d\n",bcolumns,bcolumnsc);
   printf("Program now exits\n");
   exit(1);
 }

 int indexi;

  if(!feof(fp)) {
    for(i=0;i<browsc;i++){
      fscanf(fp,"%d",&indexi);
      for(j=0;j<bcolumnsc;j++){
	for(k=0;k<bcolumnsc;k++){
	fscanf(fp,"%lf",&Cmat[i][j][k]);
	}
      }
      //      printf("reading Cmat component %d \n",indexi);
    }
  } else {
    printf("Data file error - Bmat\nthe program will be stopped now\n");
    exit(0);
  }


  fclose(fp);
  //  exit(0);

  /*
  if(!feof(fp)) {
    for(i=0;i<MAXSTEP;i++){
      fscanf(fp,"%s %s %s %s",buffer,buffer,buffer,buffer);
      fscanf(fp,"%lf %s %s %s %s ",&Energy[i],buffer,buffer,buffer,buffer);      
    }
  } else {
    printf("Data file error\nthe program will be stopped now\n");
    exit(0);
  }
  */

        
  
}

void calc_Amat () {

  int i,j;
  int kind,aind;
  //  char buffer[100];
  double ** Umat;
  double ** BmatT;
  double ** Gtemp;
  //  double ** Gmat;
  double ** GmatM1;
  double ** Atemp;

  // allocate bmat
  //  printf("B rows and columns %d %d\n",brows,bcolumns);

  Umat =  (double **) malloc((int)(bcolumns)* sizeof(double*));
  for(step=0; step<bcolumns; step++) { Umat[step]=  (double *) malloc((int)(bcolumns)* sizeof(double));}
  BmatT =  (double **) malloc((int)(bcolumns)* sizeof(double*));
  for(step=0; step<bcolumns; step++) { BmatT[step]=  (double *) malloc((int)(brows)* sizeof(double));}
  Gmat =  (double **) malloc((int)(brows)* sizeof(double*));
  for(step=0; step<brows; step++) { Gmat[step]=  (double *) malloc((int)(brows)* sizeof(double));}
  //  GmatM1 =  (double **) malloc((int)(brows)* sizeof(double*));
  //  for(step=0; step<brows; step++) { GmatM1[step]=  (double *) malloc((int)(brows)* sizeof(double));}
  Amat =  (double **) malloc((int)(bcolumns)* sizeof(double*));
  for(step=0; step<bcolumns; step++) { Amat[step]=  (double *) malloc((int)(brows)* sizeof(double));}

  kind=0;
  aind=0;
  for(i=0;i<bcolumns;i++){
    for(j=0;j<bcolumns;j++){
      if(i==j){
		Umat[i][j]=1./atoms_data[aind].atomic_mass;
		//	Umat[i][j]=1.;
	//	printf("Umat is %lf \n",Umat[i][j]);
	kind=kind+1;
	if(kind==3){
	  aind=aind+1;
	  kind=0;
	}
      } else
	Umat[i][j]=0.;
    }
  }

  /*  
  for(i=0;i<bcolumns;i++){
    for(j=0;j<bcolumns;j++){
	printf(" %lf ",Umat[i][j]);
    }
	printf("\n");
  }
  exit(0);
  
  */
  //  printf("Aind is %d \n",aind);
  //  exit(0);


  BmatT=transpose_double(Bmat,bcolumns,brows,BmatT);

  //  printf("last value transp B %lf\n",BmatT[bcolumns-1][brows-1]);

  Gtemp=prod_mat2(Umat,bcolumns,bcolumns,BmatT,bcolumns,brows);
  //  printf("last value Gtemp %lf\n",Gtemp[bcolumns-1][brows-1]);

  Gmat=prod_mat2(Bmat,brows,bcolumns,Gtemp,bcolumns,brows);
  //  printf("last value Gmat %lf\n",Gmat[brows-1][brows-1]);

  //  exit(0);
  GmatM1=invert_matrix(Gmat,brows);

  //  printf("last value GmatM1 %lf\n",GmatM1[brows-1][brows-1]);

  Atemp=prod_mat2(BmatT,bcolumns,brows,GmatM1,brows,brows);
  //  printf("last value Atemp %lf\n",Atemp[bcolumns-1][brows-1]);

  Amat=prod_mat2(Umat,bcolumns,bcolumns,Atemp,bcolumns,brows);
  //  printf("last value Amatrix %lf\n",Amat[bcolumns-1][brows-1]);

  // check consistency of B and A matrixes

  /*
  double **testmat;
  testmat=prod_mat2(Bmat,brows,bcolumns,Amat,bcolumns,brows);
  for (i=0; i<brows; i++) {
    for (j=0; j<brows; j++) {
    printf(" %lf",testmat[i][j]);
    }
    printf(" \n");
  }
  exit(0);
  */

  // now diagonalize and get eigenvalues and eigen freqs for Gmatrix
  // here we use the approach of Miyazawa (JCP 29 (1958) 246).

  int dim;

  dim=brows+1;

  double *Glambda;
  double **GL_temp;
  double **GL;
  double **BL;
  double **G_temp;
  
  Glambda = (double *) malloc(sizeof(double) * dim);
  GL_temp = (double **) malloc((int)dim*sizeof(double *) );
  G_temp = (double **) malloc((int)dim*sizeof(double *) );
  GL = (double **) malloc((int)(dim-1)*sizeof(double *) );
  BL = (double **) malloc((int)(dim-1)*sizeof(double *) );
 
  for(j=0; j<(int)(dim); j++) {
    G_temp[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
    GL_temp[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
  }

  for(j=0; j<(int)(dim-1); j++) {
    GL[j]=(double *) malloc((int)(dim-1)*sizeof(double)) ;
    BL[j]=(double *) malloc((int)(dim-1)*sizeof(double)) ;
  }


  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      if(j==i) GL_temp[i][j]=1;
      else GL_temp[i][j]=0.;
      G_temp[i][j]=0.;
    }
  }
  
  for (i=1; i<dim; i++) {
    for (j=1; j<dim; j++) {
      G_temp[i][j]= Gmat[i-1][j-1];
    }
  }
  
  jacobi(G_temp, dim-1, Glambda, GL_temp, nrot);
  eigsrt(Glambda, GL_temp, dim-1);

  for (i=1; i<dim; i++) {
    for (j=1; j<dim; j++) {
      GL[i-1][j-1]= GL_temp[i][j];
    }
  }

  //  for(j=1;j<(int)(dim);j++){
  //    printf("Glambda %6.2f  \n",Glambda[j]);
  //  }

  //  for(j=0;j<(int)(dim-1);j++){
  //    printf("GL  %8.4f  \n",GL[0][j]);
  //  }

  // now fill in the diagonal matrix BL with sqrt of eigenvalues 


  for (i=0; i<dim-1; i++) {
    for (j=0; j<dim-1; j++) {
      if(j==i) BL[i][j]=sqrt(Glambda[i+1]);
       else BL[i][j]=0.;
    }
  }
  
  // now compute the AB matrix of Miyazawa


  //  double **ABL;

  ABLmat=prod_mat2(GL,brows,brows,BL,brows,brows);

  // free vectors

  for(i=0;i<bcolumns;i++){ free (Umat[i]); } free(Umat);
  for(i=0;i<bcolumns;i++){ free (BmatT[i]); } free(BmatT);
  for(i=0;i<bcolumns;i++){ free (Atemp[i]); } free(Atemp);
  for(i=0;i<brows;i++){ free (GmatM1[i]); } free(GmatM1);
  //  for(i=0;i<brows;i++){ free (Gmat[i]); } free(Gmat);
  for(i=0;i<bcolumns;i++){ free (Gtemp[i]); } free(Gtemp);

  for(i=0;i<brows;i++){ free (BL[i]); } free(BL);
  for(i=0;i<brows;i++){ free (GL[i]); } free(GL);
  for(i=0;i<brows+1;i++){ free (GL_temp[i]); } free(GL_temp);

  free(Glambda);

  //  free(e_temp);


  //  exit(0);

}


void intcoord_hessian (int step, double **FC, double *grad) {

  int i,j,k,istep;
  //  int kind,aind;
  //  char buffer[100];
  double ** FCtemp;
  double ** Ftemp;
  double ** Ftemp2;
  double ** Ftemp3;
  double ** Ftemp4;
  double ** AmatT;
  double ** ABLmatT;
  double ** Fmat_int;
  double ** Fmat_int2;
  double ** Gradmat;
  double ** Pmat1;
  double ** Pmat2;
  double * Grad_int;
  double * Grad_corr;
  //  double ** ACmat_temp1;
  double *** ACmat_temp2;
  double ** ACmat_temp3;
  //  double ** ACmat_temp4;
  double ** ACmat_tempa;
  double ** ACmat_tempb;

  //  double * frequencies_int;

  // convert Hessian to Angstron

  double conv1;

  //  conv1=bohr;
  FCtemp =  (double **) malloc((int)(bcolumns)* sizeof(double*));
  for(istep=0; istep<bcolumns;istep++) { FCtemp[istep]= (double *) malloc((int)(bcolumns)* sizeof(double));}
  AmatT =  (double **) malloc((int)(brows)* sizeof(double*));
  for(istep=0; istep<brows; istep++) { AmatT[istep]=  (double *) malloc((int)(bcolumns)* sizeof(double));}
  ABLmatT =  (double **) malloc((int)(brows)* sizeof(double*));
  for(istep=0; istep<brows; istep++) { ABLmatT[istep]=  (double *) malloc((int)(brows)* sizeof(double));}
  Grad_int=(double *) malloc((int)(brows)* sizeof(double));
  Grad_corr=(double *) malloc((int)(bcolumns)* sizeof(double));
  Gradmat =  (double **) malloc((int)(brows)* sizeof(double*));
  for(istep=0; istep<brows; istep++) { Gradmat[istep]=  (double *) malloc((int)(brows)* sizeof(double));}

  //  ACmat_temp1 =  (double ***) malloc((int)(brows)* sizeof(double**));
  //  for(step=0; step<brows; step++) { ACmat_temp1[step]=  (double **) malloc((int)(bcolumns)* sizeof(double*));}
  //  for(step=0; step<brows; step++){ for(j=0; j <(int)(bcolumns); j++) { ACmat_temp1[step][j]=  (double *) malloc((int)(brows)* sizeof(double));}}

  ACmat_temp2 =  (double ***) malloc((int)(brows)* sizeof(double**));
  for(istep=0; istep<brows; istep++) { ACmat_temp2[istep]=  (double **) malloc((int)(brows)* sizeof(double*));}
  for(istep=0; istep<brows; istep++){ for(j=0; j <(int)(brows); j++) { ACmat_temp2[istep][j]=  (double *) malloc((int)(brows)* sizeof(double));}}

  ACmat_temp3 =  (double **) malloc((int)(brows)* sizeof(double*));
  for(istep=0; istep<brows; istep++) { ACmat_temp3[istep]=  (double *) malloc((int)(brows)* sizeof(double));}

  //  ACmat_temp1 =  (double **) malloc((int)(bcolumns)* sizeof(double*));
  //  for(step=0; step<bcolumns; step++) { ACmat_temp1[step]=  (double *) malloc((int)(bcolumns)* sizeof(double));}

  //  ACmat_temp4 =  (double **) malloc((int)(brows)* sizeof(double*));
  //  for(step=0; step<brows; step++) { ACmat_temp4[step]=  (double *) malloc((int)(brows)* sizeof(double));}

    
  //mass unweight the Hessian and go to Angstrom

  conv1=bohr*bohr;

  for (i=0;i<bcolumns;i++){
    for (j=0;j<bcolumns;j++){
      FCtemp[i][j]=FC[i][j]*conv1; 
    }
  }


  // set gradient to 0 for saddle point
  for (i=0;i<bcolumns;i++){
    Grad_corr[i]=1.;
    if(step==saddlep){
      Grad_corr[i]=0.;
    }
  }

  FCtemp=force_constants_mass_unweight(FCtemp);

  /*  
  printf(" unproj FC matrix \n"); 
  for (i=0;i<bcolumns;i++){
    for (j=0;j<bcolumns;j++){
      printf(" %lf ",FCtemp[i][j]); 
    }
      printf("\n"); 
  }
  printf(" A matrix \n"); 
  for (i=0;i<bcolumns;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",Amat[i][j]); 
    }
      printf("\n"); 
  }
  printf(" G matrix \n"); 
  for (i=0;i<brows;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",Gmat[i][j]); 
    }
      printf("\n"); 
  }
  */

  AmatT=transpose_double(Amat,brows,bcolumns,AmatT);

  Ftemp=prod_mat2(FCtemp,3*ATOMS,3*ATOMS,Amat,bcolumns,brows);

  Ftemp2=prod_mat2(AmatT,brows,bcolumns,Ftemp,3*ATOMS,brows);

  /*  
  printf(" A matrix \n"); 
  for (i=0;i<brows;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",Amat[i][j]); 
    }
      printf("\n"); 
  }
  */
  // Ftemp2 is the Hessian in internal coordinates. 

  // now we compute the corretion for the gradient

  //first we compute the gradient in Hartree/ang
  // we also put back weights 

  int iatom=0;
  int iind=0;
  double gcheck=0.;
  for (i=0;i<brows;i++){
    Grad_int[i]=0.;
    //    printf(" %d\n ",i); 
    //    printf(" grad %d %lf \n ",i,grad[i]); 
    iind=0;
    iatom=0;
    for (j=0;j<bcolumns;j++){
      if(iind==3){
	iatom=iatom+1;
	iind=0;
      }
      iind=iind+1;
      //      printf(" %d\n ",iatom); 
      //      Grad_int[i]+=grad[j]*sqrt(atoms_data[iatom].atomic_mass)*AmatT[i][j];
      Grad_int[i]+=Grad_corr[j]*grad[j]*sqrt(atoms_data[iatom].atomic_mass)*AmatT[i][j]*bohr;
      gcheck+=grad[j]*Grad_corr[j];
      //      Grad_int[i]+=grad[j]*sqrt(atoms_data[iatom].atomic_mass)*AmatT[i][j];
    }
  }

  for (j=0;i<brows;j++){
    for (i=0;i<brows;i++){
      for (k=0;k<brows;k++){
	ACmat_temp2[j][i][k]=0.;
      }
    }
  }

  /*  
    for (i=0;i<bcolumns;i++){
      for (k=0;k<bcolumns;k++){
	ACmat_temp1[i][k]=Cmat[0][i][k];
      }
    }

    for (i=0;i<brows;i++){
      for (k=0;k<brows;k++){
	ACmat_temp4[i][k]=0.;
      }
    }
  
  */

    for (j=0;j<brows;j++){
      ACmat_tempa=prod_mat2(Cmat[j],bcolumns,bcolumns,Amat,bcolumns,brows);
    //  ACmat_tempa=prod_mat2(ACmat_temp1,bcolumns,bcolumns,Amat,bcolumns,brows);
      ACmat_tempb=prod_mat2(AmatT,brows,bcolumns,ACmat_tempa,bcolumns,brows);
      for (i=0;i<brows;i++){
	for (k=0;k<brows;k++){
	  //	  ACmat_temp2[j][i][k]=ACmat_tempb[i][k]*Grad_int[j];
	  ACmat_temp2[j][i][k]=ACmat_tempb[i][k]*Grad_int[j];
	  //	  ACmat_temp4[i][k]=ACmat_tempb[i][k];
	}
      }
      for(i=0;i<bcolumns;i++){free(ACmat_tempa[i]);}free(ACmat_tempa);
      for(i=0;i<brows;i++){free(ACmat_tempb[i]);}free(ACmat_tempb);
      //      printf("Gr %d %lf\n ",j,Grad_int[j]); 
    }


  /*
  for (j=0;j<brows;j++){
    ACmat_temp1[j]=prod_mat2(Cmat[j],bcolumns,bcolumns,Amat,bcolumns,brows);
    ACmat_temp2[j]=prod_mat2(AmatT,brows,bcolumns,ACmat_temp1[j],bcolumns,brows);
    for (i=0;i<brows;i++){
      for (k=0;k<brows;k++){
	ACmat_temp2[j][i][k]=ACmat_temp2[j][i][k]*Grad_int[j];
      }
    }
   }

  */


  for (i=0;i<brows;i++){
    for (k=0;k<brows;k++){
      ACmat_temp3[i][k]=0.;
    }
  }


  for (j=0;j<brows;j++){
    for (i=0;i<brows;i++){
      for (k=0;k<brows;k++){
	ACmat_temp3[i][k]=ACmat_temp3[i][k]+ACmat_temp2[j][i][k];
      }
    }
  }


  /*         
  printf("ACM_3\n"); 
  for (i=0;i<brows;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",ACmat_temp4[i][j]); 
    }
    printf("\n"); 
  }
  printf("\n"); 

  printf("ACtemp_2\n"); 
  for (i=0;i<brows;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",ACmat_temp2[0][i][j]); 
    }
    printf("\n"); 
  }
  printf("\n"); 


      printf("Cmat\n"); 
  for (i=0;i<bcolumns;i++){
    for (j=0;j<bcolumns;j++){
      printf(" %lf ",Cmat[0][i][j]); 
    }
      printf("\n"); 
  }
      printf("\n"); 
  */
  


    //        double inc_fact=1.;

  for (j=0;j<brows;j++){
    for (k=0;k<brows;k++){
      //with C matrix contribution
      // negative (should be the correct one!)
      //            Ftemp2[j][k]=Ftemp2[j][k]-ACmat_temp3[j][k]; 
      // positive
            Ftemp2[j][k]=Ftemp2[j][k]+ACmat_temp3[j][k]; 
      //without C matrix contribution
      //            Ftemp2[j][k]=Ftemp2[j][k]+ACmat_temp3[j][k]*10.; 
      //	  Ftemp2[j][k]=Ftemp2[j][k];
    }
  }


  double grnorm=0.;

  for (i=0;i<brows;i++){
    for (k=0;k<brows;k++){
      //grnorm+=Grad_int[k]*Grad_int[k]*Gmat[i][k];
      grnorm+=Grad_int[i]*Grad_int[k]*Gmat[i][k];
    }
  }

  // check modification
    grnorm=grnorm;

  if(grnorm != 0.){
    for (i=0;i<brows;i++){
      for (k=0;k<brows;k++){
	Gradmat[i][k]=Grad_int[i]*Grad_int[k]/grnorm;
      }
    }
  }

  Pmat1=prod_mat2(Gmat,brows,brows,Gradmat,brows,brows);
  Pmat2=prod_mat2(Gradmat,brows,brows,Gmat,brows,brows);

  for (i=0;i<brows;i++){
    for (k=0;k<brows;k++){
      if(i==k){
	Pmat1[i][k]=1.-Pmat1[i][k];
	Pmat2[i][k]=1.-Pmat2[i][k];
      } else{
	Pmat1[i][k]=-Pmat1[i][k];
	Pmat2[i][k]=-Pmat2[i][k];
      }
    }
  }


  /*       
      printf("Grad Mat\n"); 
  for (i=0;i<brows;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",Gradmat[i][j]); 
    }
      printf("\n"); 
  }
      printf("\n"); 
  */
      //      exit(0);

  Ftemp4=prod_mat2(Ftemp2,brows,brows,Pmat1,brows,brows);
  Fmat_int2=prod_mat2(Pmat2,brows,brows,Ftemp4,brows,brows);

  //Now we can determine the eigenvalues we use the AB mat of Miyazawa

  // with reaction path projection (but not at the TS)

  if(gcheck!=0){
    Ftemp3=prod_mat2(Fmat_int2,brows,brows,ABLmat,brows,brows);
  } else {
        Ftemp3=prod_mat2(Ftemp2,brows,brows,ABLmat,brows,brows);
  }


  // without reaction path projection
  //        Ftemp3=prod_mat2(Ftemp2,brows,brows,ABLmat,brows,brows);

  // we transpose ABL

  ABLmatT=transpose_double(ABLmat,brows,brows,ABLmatT);

  //we compute the symmetrized F internal matrix

  Fmat_int=prod_mat2(ABLmatT,brows,brows,Ftemp3,brows,brows);

  // now convert the Hessian to J/m^2*kg

  conv1=joule/a_mass*1.0e10*1.0e10;

  for (j=0;j<brows;j++){
    for (k=0;k<brows;k++){
      Fmat_int[j][k]=Fmat_int[j][k]*conv1; 
    }
  }


  // now calculate eigenvalues

  int dim;

  dim=brows+1;

  double *lambda;
  double **L;
  double **FC_temp;
  
  lambda = (double *) malloc(sizeof(double) * dim);
  lambdaint = (double *) malloc(sizeof(double) * dim);
  Freqint = (double *) malloc(sizeof(double) * brows);
  L = (double **) malloc((int)dim*sizeof(double *) );
  FC_temp = (double **) malloc((int)dim*sizeof(double *) );
 
  for(j=0; j<(int)(dim); j++) {
    FC_temp[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
    L[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
  }

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      if(j==i) L[i][j]=1;
      else L[i][j]=0.;
      FC_temp[i][j]=0.;
    }
  }
  
  for (i=1; i<dim; i++) {
    for (j=1; j<dim; j++) {
      FC_temp[i][j]= Fmat_int[i-1][j-1];
    }
  }
  
  jacobi(FC_temp, dim-1, lambda, L, nrot);
  eigsrt(lambda, L, dim-1);

  for(j=1;j<(int)(dim);j++){
    printf(" %6.2f %6.2f \n",sqrt(fabs(lambda[j]))/(2 * pi * c_light_cm_s ),lambda[j]/1.0E20);
  }
  lambdaint[0]=0.;
  for(j=1;j<(int)(dim);j++){
    Freqint[j-1]=sqrt(fabs(lambda[j]))/(2 * pi * c_light_cm_s )*lambda[j]/fabs(lambda[j]);
    lambdaint[j]=lambda[j];
  }
  
  FILE *hrproj_freq;
  if((hrproj_freq=fopen("hrprojint_freq.dat","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","hrproj_freq.dat");
    exit(1);
  }

  for(j=1;j<brows+1;j++){
    if(lambda[j]>0.) {
      fprintf(hrproj_freq," %6.2f %6.2f \n",sqrt(fabs(lambda[j]))/(2 * pi * c_light_cm_s ),lambda[j]/1.0E20);
    }
    else {
      fprintf(hrproj_freq," %6.2f %6.2f \n",-sqrt(fabs(lambda[j]))/(2 * pi * c_light_cm_s ),lambda[j]/1.0E20);
    }
  }

  fclose(hrproj_freq);
  
  // now renumber eigenvectors

  double **L_temp; 

  L_temp = (double **) malloc((int)brows*sizeof(double *) );
  for(j=0; j<(int)(brows); j++) {
    L_temp[j]=(double *) malloc((int)(brows)*sizeof(double)) ;
  }

  for(j=0; j<(int)(brows); j++) {
    for(i=0; i<(int)(brows); i++) {
      //      L_temp[i][j]=L[i+1][j+1]/L[1][j+1];
      L_temp[i][j]=L[i+1][j+1];
    }
  }

  double **L_int; 
  double **L_intM1; 
  double **L_intM1T; 
  double **Ctemp; 
  double **Ctemp1; 
  double **Cdiag; 
  double **Lcart_temp; 
  //  double **Lcart; 
  double *normfact;

  // the internal coordinates Eigenvectors are recovered by multiplication with the AB mat

  L_int=prod_mat2(ABLmat,brows,brows,L_temp,brows,brows);

  // the inverse of L_int is determined

  L_intM1=invert_matrix(L_int,brows);

  L_intM1T =  (double **) malloc((int)(brows)* sizeof(double*));
  for(istep=0; istep<brows; istep++) { L_intM1T[istep]=  (double *) malloc((int)(brows)* sizeof(double));}

  Cdiag =  (double **) malloc((int)(brows)* sizeof(double*));
  for(istep=0; istep<brows; istep++) { Cdiag[istep]=  (double *) malloc((int)(brows)* sizeof(double));}

  normfact =  (double *) malloc((int)(brows)* sizeof(double));

  L_intM1T=transpose_double(L_intM1,brows,brows,L_intM1T);

  //  Ctemp1=prod_mat2(L_int,brows,brows,L_intT,brows,brows);

  Ctemp1=prod_mat2(Gmat,brows,brows,L_intM1T,brows,brows);
  Ctemp=prod_mat2(L_intM1,brows,brows,Ctemp1,brows,brows);
  

  for (i=0; i<brows; i++) {
    for (j=0; j<brows; j++) {
      if(j==i) Cdiag[i][j]=sqrt(Ctemp[i][j]);
      else Cdiag[i][j]=0.;
    }
  }

  // now we can compute the Cartesian Eigenvectors
  // please note, Lcart is globally defined, but the memery is allocted here and freed in the main

  Lcart_temp=prod_mat2(L_int,brows,brows,Cdiag,brows,brows);
  Lcartint=prod_mat2(Amat,bcolumns,brows,Lcart_temp,brows,brows);

  //  Lcart=prod_mat2(Amat,bcolumns,brows,L_int,brows,brows);

  // check Eigenvectors
  /*
  printf("L_intM\n"); 
  for (i=0;i<brows;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",Ctemp1[i][j]); 
    }
      printf("\n"); 
  }
      printf("\n"); 

  printf("G Mat\n"); 
  //  for (i=0;i<bcolumns;i++){
  for (i=0;i<brows;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",Gmat[i][j]); 
    }
      printf("\n"); 
  }
      printf("\n"); 
  */
         
      // now normalize over masses


  iatom=0;
  iind=0;
  for (j=0;j<brows;j++){
    normfact[j]=0.;
    iind=0;
    iatom=0;
    for (i=0;i<bcolumns;i++){
      if(iind==3){
	iatom=iatom+1;
	iind=0;
      }
      iind=iind+1;
      Lcartint[i][j]=Lcartint[i][j]*sqrt(atoms_data[iatom].atomic_mass);
      normfact[j]+=Lcartint[i][j]*Lcartint[i][j]*(atoms_data[iatom].atomic_mass);
    }
  }

  // here normalization condition suppressed to match Cartesian Eigen vector definition 
  // in SCT Cartesian implementation.

  for (j=0;j<brows;j++){
    for (i=0;i<bcolumns;i++){
      //      Lcart[i][j]=Lcart[i][j]/sqrt(normfact[j]);
      Lcartint[i][j]=Lcartint[i][j];
    }
  }

  /*
  printf("Cart Eigen Mat\n"); 
  for (i=0;i<bcolumns;i++){
    for (j=0;j<brows;j++){
      printf(" %lf ",Lcartint[i][j]); 
    }
      printf("\n"); 
  }
      printf("\n"); 
  */

  //  if(step==saddlep){
  //    exit(0);
  //  }

  for(i=0;i<dim;i++){free(L[i]);}free(L);
  for(i=0;i<dim;i++){free(FC_temp[i]);}free(FC_temp);
  for(i=0;i<bcolumns;i++){free(FCtemp[i]);}free(FCtemp);
  for(i=0;i<brows;i++){free(Ftemp2[i]);}free(Ftemp2);
  for(i=0;i<brows;i++){free(Ftemp3[i]);}free(Ftemp3);
  for(i=0;i<brows;i++){free(Ftemp4[i]);}free(Ftemp4);

  for(i=0;i<brows;i++){free(AmatT[i]);}free(AmatT);
  for(i=0;i<brows;i++){free(ABLmatT[i]);}free(ABLmatT);
  for(i=0;i<brows;i++){free(Gradmat[i]);}free(Gradmat);
  for(i=0;i<brows;i++){for(j=0;j<brows;j++){free(ACmat_temp2[i][j]);}free(ACmat_temp2[i]);}free(ACmat_temp2);
  for(i=0;i<brows;i++){free(ACmat_temp3[i]);}free(ACmat_temp3);
  for(i=0;i<brows;i++){free(Pmat1[i]);}free(Pmat1);
  for(i=0;i<brows;i++){free(Pmat2[i]);}free(Pmat2);

  for(i=0;i<brows;i++){free(Fmat_int2[i]);}free(Fmat_int2);
  for(i=0;i<brows;i++){free(Fmat_int[i]);}free(Fmat_int);
  for(i=0;i<brows;i++){free(L_temp[i]);}free(L_temp);
  for(i=0;i<brows;i++){free(Lcart_temp[i]);}free(Lcart_temp);
  for(i=0;i<brows;i++){free(L_int[i]);}free(L_int);
  for(i=0;i<brows;i++){free(L_intM1[i]);}free(L_intM1);
  for(i=0;i<brows;i++){free(L_intM1T[i]);}free(L_intM1T);

  for(i=0;i<brows;i++){free(Cdiag[i]);}free(Cdiag);
  for(i=0;i<brows;i++){free(Ctemp[i]);}free(Ctemp);
  for(i=0;i<brows;i++){free(Ctemp1[i]);}free(Ctemp1);

  free(lambda);
  free(normfact);
  free(Grad_int);
  free(Grad_corr);

  //  exit(0);
  //Hessian[step]

}


void read_rxfile () {
  int i,i1;
  double i2;
  char buffer[100];
  FILE *fp;
  char filename[20];

  sprintf(filename,"Rx_coord.txt");
  if((fp=fopen(filename,"r"))==NULL){
    printf("main: IMPOSSIBILE APRIRE IL FILE: %s\nProgram now exits\n",filename);
    exit(1);
  }
  
  if(!feof(fp)) {
    fscanf(fp,"%s %s %s",buffer,buffer,buffer);
    for(i=0;i<MAXSTEP;i++){  
      fscanf(fp,"%d %lf %lf ",&i1,&i2,&Rx_coord[i]);
      //      Rx_coord[i]=Rx_coord[i]/1.07;
      printf("%d\t%lf\n",i1,Rx_coord[i]);
    }
  } else {
      printf("Data file error\nthe program will be stopped now\n");
      exit(0);
  }
  fclose(fp);
}

/*
void read_mueff_file (FILE *fp) {
  int i;
  
  char buffer[100];
 
  if(!feof(fp)) {
    
    fscanf(fp,"%d",&dim_truong_mueff);
    mueff_read=(double *) malloc( dim_truong_mueff*sizeof(double ) );
    x_mueff_read=(double *) malloc( dim_truong_mueff*sizeof(double ) );
  
    fscanf(fp,"%s %s",buffer,buffer);
    for(i=0;i<dim_truong_mueff;i++){  
      fscanf(fp,"%lf %lf ",&x_mueff_read[i],&mueff_read[i]);
      printf("%lf\t%lf\n",x_mueff_read[i],mueff_read[i]);
    }
  } else {
      printf("Data file error\nthe program will be stopped now\n");
      exit(0);
  }

}
*/

void read_newpes_file () {

  int i;
  char buffer[100];
  FILE *fp;
  char filename[20];

  sprintf(filename,"read_pes.txt");
  if((fp=fopen(filename,"r"))==NULL){
    printf("main: IMPOSSIBILE APRIRE IL FILE: %s\nProgram now exits\n",filename);
    exit(1);
  }

  if(!feof(fp)) {
    fscanf(fp,"%d",&dim_newpes);

    pes_read=(double *) malloc( dim_newpes*sizeof(double ) );
    x_pes_read=(double *) malloc( dim_newpes*sizeof(double ) );

    fscanf(fp,"%s %s",buffer,buffer);
    for(i=0;i<dim_newpes;i++){  
      fscanf(fp,"%lf %lf ",&x_pes_read[i],&pes_read[i]);
      printf("%lf\t%lf\n",x_pes_read[i],pes_read[i]);
    }
  } else {
      printf("Data file error\nthe program will be stopped now\n");
      exit(0);
  }

  fclose(fp);

}

void write_traj () {

  FILE *trajec_file;
  //  char filename[20];
  int j;

  if((trajec_file=fopen("trajec.xyz","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","trajec.txt");
    exit(1);
  }

  for(step=0;step<MAXSTEP;step++){
    fprintf(trajec_file,"    %u  \n",ATOMS);
    fprintf(trajec_file,"step number %d\n",step);
    for(j=0;j<ATOMS;j++){
      fprintf(trajec_file,"%s\t%lf\t%lf\t%lf\n",atoms_data[j].atom_name,atoms_data[j].ext_coord_x[step],atoms_data[j].ext_coord_y[step],atoms_data[j].ext_coord_z[step]);
    }    
  }
  fclose(trajec_file);

}

void convert_coordinates_to_bohr(int step){

  int i;

  //  printf("Coordinates in atomic units\n");
  //  printf("Atom\t     x\t     y\t     z\n");

  for(i=0;i<ATOMS;i++){
    
    atoms_data[i].ext_coord_x[step]= atoms_data[i].ext_coord_x[step]*bohr;
    atoms_data[i].ext_coord_y[step]= atoms_data[i].ext_coord_y[step]*bohr;   
    atoms_data[i].ext_coord_z[step]= atoms_data[i].ext_coord_z[step]*bohr;
   
    //    printf("%d\t%lf\t%lf\t%lf\n",i,atoms_data[i].ext_coord_x[step],atoms_data[i].ext_coord_y[step],atoms_data[i].ext_coord_z[step]);
  }
}

void calculate_center_of_mass(int step){

  int i;
  double Total_mass=0;

  CM_x[step] =0.0;
  CM_y[step] =0.0;
  CM_z[step] =0.0;

  for(i=0;i<ATOMS;i++){

    CM_x[step] =   CM_x[step] + atoms_data[i].atomic_mass *  atoms_data[i].ext_coord_x[step];
    CM_y[step] =   CM_y[step] + atoms_data[i].atomic_mass *  atoms_data[i].ext_coord_y[step];
    CM_z[step] =   CM_z[step] + atoms_data[i].atomic_mass *  atoms_data[i].ext_coord_z[step];

    Total_mass =  Total_mass + atoms_data[i].atomic_mass;
  }

  CM_x[step]=CM_x[step]/Total_mass;
  CM_y[step]=CM_y[step]/Total_mass;
  CM_z[step]=CM_z[step]/Total_mass;

}


void set_origin_to_cofm(int step ){
  
  int i;
  //  printf("Coordinates with origin set to Center of Mass\n");
  //  printf("Atom\t     x\t     y\t     z\n");

  for(i=0;i<ATOMS;i++){
    
    atoms_data[i].ext_coord_x[step]= atoms_data[i].ext_coord_x[step]-CM_x[step];    
    atoms_data[i].ext_coord_y[step]= atoms_data[i].ext_coord_y[step]-CM_y[step];   
    atoms_data[i].ext_coord_z[step]= atoms_data[i].ext_coord_z[step]-CM_z[step];
   
    //    printf("%d\t%lf\t%lf\t%lf\n",i, atoms_data[i].ext_coord_x[step], atoms_data[i].ext_coord_y[step], atoms_data[i].ext_coord_z[step]);
  }
}

void calculate_inertia(int step){

  int i,j;
  int dim=4;
  double temp=0,temp1=0,temp2=0,temp3=0,temp4=0,temp5=0;
  
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      I[step][i][j]=0.;
    }
  }

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      v[step][i][j]=0.;
    }
  }

 
  for(i=0;i<ATOMS;i++){

    //I(1,1)=Ixx=sum(mi x yi^2 + mi x zi^2)
    temp =  temp  + atoms_data[i].atomic_mass * ((atoms_data[i].ext_coord_y[step]*atoms_data[i].ext_coord_y[step])+(atoms_data[i].ext_coord_z[step]*atoms_data[i].ext_coord_z[step]));      

    //I(2,2)=Iyy=sum(mi x xi^2 + mi x zi^2)
    temp1 = temp1 + atoms_data[i].atomic_mass * ((atoms_data[i].ext_coord_x[step]*atoms_data[i].ext_coord_x[step])+(atoms_data[i].ext_coord_z[step]*atoms_data[i].ext_coord_z[step]));
   
    //I(3,3)=Izz=sum(mi x yi^2 + mi x zi^2)
    temp2 = temp2 + atoms_data[i].atomic_mass * ((atoms_data[i].ext_coord_y[step]*atoms_data[i].ext_coord_y[step])+(atoms_data[i].ext_coord_x[step]*atoms_data[i].ext_coord_x[step]));

    //I(1,2)=Ixy=-sum(mi x yi x xi)
    temp3 = temp3 - atoms_data[i].atomic_mass * (atoms_data[i].ext_coord_y[step]*atoms_data[i].ext_coord_x[step]);

    //I(1,3)=Ixz=-sum(mi x xi x zi)
    temp4 = temp4 - atoms_data[i].atomic_mass * (atoms_data[i].ext_coord_x[step]*atoms_data[i].ext_coord_z[step]);

    //I(2,3)=Iyz=-sum(mi x yi x zi)
    temp5 = temp5 - atoms_data[i].atomic_mass * (atoms_data[i].ext_coord_y[step]*atoms_data[i].ext_coord_z[step]);	
    //    printf("%le\n", atoms_data[i].atomic_mass);

  }

  I[step][1][1]=temp; 
  I[step][2][2]=temp1;
  I[step][3][3]=temp2;
  I[step][1][2]=temp3;
  I[step][1][3]=temp4;
  I[step][2][3]=temp5;
  I[step][2][1]=I[step][1][2];
  I[step][3][2]=I[step][2][3];
  I[step][3][1]=I[step][1][3];

  //    printf("Moment of Inertia tensor expressed in uma*bohr^2\n");
  //    printf("%le\t%le\t%le\n", I[step][1][1],I[step][1][2],I[step][1][3]);
  //  printf("%le\t%le\t%le\n", I[step][2][1],I[step][2][2],I[step][2][3]);
  //  printf("%le\t%le\t%le\n", I[step][3][1],I[step][3][2],I[step][3][3]);
    //    exit(0);

  //calculate eigenvectors and eigenvalues
  jacobi(I[step], 3, eigen[step], v[step], nrot);   

  //order eigenvectors and eigenvalues
  eigsrt(eigen[step], v[step], 3);

  //    printf("Principal moments of inertia\n");
  //    printf("%lf\t%lf\t%lf\n", eigen[step][1],eigen[step][2],eigen[step][3]);  

  //  printf("Eigenvectors\n");
  //  printf("%lf\t%lf\t%lf\n", v[step][1][1],v[step][1][2],v[step][1][3]); 
  //  printf("%lf\t%lf\t%lf\n", v[step][2][1],v[step][2][2],v[step][2][3]); 
  //  printf("%lf\t%lf\t%lf\n", v[step][3][1],v[step][3][2],v[step][3][3]); 

  // refill inertia tensor to correct for jacobi modification


  I[step][1][2]=I[step][2][1];
  I[step][2][3]=I[step][3][2];
  I[step][1][3]=I[step][3][1];


  //    printf("Moment of Inertia tensor expressed in uma*bohr^2\n");
  //    printf("%le\t%le\t%le\n", I[step][1][1],I[step][1][2],I[step][1][3]);
  //  printf("%le\t%le\t%le\n", I[step][2][1],I[step][2][2],I[step][2][3]);
  //  printf("%le\t%le\t%le\n", I[step][3][1],I[step][3][2],I[step][3][3]);
  //  exit(0);

}


void coordinates_in_principal_axes(int step){
  int i;
  double temp_x,temp_y,temp_z;
  
  for(i=0;i<ATOMS;i++){
    
     temp_x= atoms_data[i].ext_coord_x[step] * v[step][1][1] + atoms_data[i].ext_coord_y[step]* v[step][2][1] + atoms_data[i].ext_coord_z[step] * v[step][3][1];

     temp_y= atoms_data[i].ext_coord_x[step] * v[step][1][2] + atoms_data[i].ext_coord_y[step] * v[step][2][2] + atoms_data[i].ext_coord_z[step] * v[step][3][2];

     temp_z= atoms_data[i].ext_coord_x[step] * v[step][1][3] + atoms_data[i].ext_coord_y[step] * v[step][2][3] + atoms_data[i].ext_coord_z[step] * v[step][3][3];

     atoms_data[i].ext_coord_x[step] =  temp_x;
     atoms_data[i].ext_coord_y[step] =  temp_y;
     atoms_data[i].ext_coord_z[step] =  temp_z;

  }

  /*
 printf("\nCoordinates in the principal axes system\n");
 for(i=0;i<ATOMS;i++){
  
   printf("%d\t%lf\t%lf\t%lf\n",i, atoms_data[i].ext_coord_x[step], atoms_data[i].ext_coord_y[step], atoms_data[i].ext_coord_z[step]);
 }
  */
}

double ** force_constants_mass_weight(double **FC){

  
  int i, j,k;
  //  double *mass;
  
  mass=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );
    
  k=0;  
  for(i=0;i<3*ATOMS;i=i+3){
    for(j=0;j<3;j++){
      mass[i+j]=atoms_data[k].atomic_mass;
    }
    k=k+1;
  }
  
  //  printf("\n");
        
  for(i=0;i<3*ATOMS;i++){
    for(j=0;j<3*ATOMS;j++){
      FC[i][j]= FC[i][j]/sqrt(mass[i]*mass[j]);
    }
  }

  //  free(mass);
 
  //the hessian missing elements are filled in, since the input is triangular

  for(i=0;i<(int)(3*ATOMS);i++){
    //    for(j=0;j<(int)(3*ATOMS);j++){
    for(j=i;j<(int)(3*ATOMS);j++){
      FC[i][j]=FC[j][i];
    }  
  }

  return FC;
}

double ** force_constants_mass_unweight(double **FC){

  
  int i, j,k;
  double *mass;
  
  mass=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );
 
  k=0;
  
    for(i=0;i<3*ATOMS;i=i+3){
      for(j=0;j<3;j++){
	mass[i+j]=atoms_data[k].atomic_mass;
	  }
      k=k+1;
  }
  
    //   printf("mass vector\n");
    //    for(i=0;i<3*ATOMS;i++){
      //  printf("%lf   ",mass[i]);
      //    }
    
    //    printf("\n");
    
    
    for(i=0;i<3*ATOMS;i++){
      for(j=0;j<3*ATOMS;j++){
	FC[i][j]= FC[i][j]*sqrt(mass[i]*mass[j]);
      }
    }
    free(mass);
 
    //  printf("\nHessian in mass weighted cartesian coordinates\n");
  for(i=0;i<(int)(3*ATOMS);i++){
    for(j=0;j<(int)(3*ATOMS);j++){
      FC[i][j]=FC[j][i];
      //printf("%lf\t",FC[i][j]);
    }  
    //  printf("\n");
  }
  //printf("\n");  
    return FC;

}

double **diagonalize_hessian(int dim, double**FC){

  //int dim=(3*ATOMS-6)+1;
 
  dim=dim+1;
  int j,i;
  double **FC_temp,**L, *e_temp;

  lambda = (double *) malloc(sizeof(double) * dim);
  e_temp = (double *) malloc(sizeof(double) * dim);
  L = (double **) malloc((int)dim*sizeof(double *) );
  FC_temp = (double **) malloc((int)dim*sizeof(double *) );
 
  for(j=0; j<(int)(dim); j++) {
    FC_temp[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
    L[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
  }

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      if(j==i) L[i][j]=1;
      else L[i][j]=0.;
      FC_temp[i][j]=0.;
    }
  }
  
  for (i=1; i<dim; i++) {
    for (j=1; j<dim; j++) {
      FC_temp[i][j]= FC[i-1][j-1];
    }
  }
  
  jacobi(FC_temp, dim-1, lambda, L, nrot);
  eigsrt(lambda, L, dim-1);

  printf("Hessian diagonalization\n Eigenvalues --- ");
  for(j=1; j<(int)(dim); j++) {    
    printf("%lf\t",lambda[j]);
  }

  printf("\n");	

  if(rearrange == 1) {

    double *lambda_ord;
    double **Lint_temp;
    int k;

    lambda_ord = (double *) malloc(sizeof(double) * dim);
    Lint_temp = (double **) malloc((int)dim*sizeof(double *) );
    for(j=0; j<(int)(dim); j++) {
      Lint_temp[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
    }

    for (i=0; i<dim; i++) {lambda_ord[i]=0.;}
    //    lambda_ord[1]=0.;

    // find the largest eigenvalue and eigenvector and put it in position 1
    for (j=1; j<dim; j++) {
      if(fabs(lambda[j]) > fabs(lambda_ord[1])) {
	lambda_ord[1]=lambda[j];
	for (i=1; i<dim; i++) {
	  //	  Lint_temp[1][i]=L[j][i];
	  Lint_temp[i][1]=L[i][j];
	}
      }
    }

    //    one by one, determines the order of the other eigenvalues and eigenstates
    for (i=2; i<dim; i++) {
      //      lambda_ord[i]=0.;
      for (j=1; j<dim; j++) {
	if((fabs(lambda[j]) > fabs(lambda_ord[i])) && (fabs(lambda[j])<fabs(lambda_ord[i-1])))  {
	  lambda_ord[i]=lambda[j];
	  for (k=1; k<dim; k++) {
	    //	  Lint_temp[i][k]=L[j][k];
	    Lint_temp[k][i]=L[k][j];
	  }
	}
      }
    } 

    // reassignes new order of eigenvalues and eigenstates

    for(i=1; i<(int)(dim); i++) {    
      lambda[i]=lambda_ord[i];
       //    printf("%lf\t",lambda_ord[i]);
      for(j=1; j<(int)(dim); j++) {    
	L[i][j]=Lint_temp[i][j];
      }
    }

    for(i=0;i<dim;i++){ free (Lint_temp[i]); } free(Lint_temp);
    free(lambda_ord);

    printf("rearranging order of eigenvalues and eigenstates\n");	
    printf("\n");	

  } 
  else if (rearrange == 0) {
    printf("order of eigenvalues and eigenstates not rearranged\n");	
  } 
  else {
    printf("rearrange parameter must be 0 or 1\n");	
    printf("change and restart\n");	
    exit(-1);
  } 

  //  exit(-1);

  for(i=0;i<dim;i++){ free (FC_temp[i]); } free(FC_temp);
  free(e_temp);
	       			
  return L;
}

double *calc_freq(double *lam, double *freq){

  int i; 
  double constant=2*pi*c_light_cm_s*seconds;

  //  printf("\n");

  //printf("\n u.t.a --> seconds= %le \n",seconds);
  // printf("\n u.m.a./me = %le \n",massunit);

  for(i=0;i<3*ATOMS;i++){

    if(lam[i+1]>0)

      freq[i]=sqrt(fabs(lam[i+1])/massunit)/constant;
    
    else freq[i]=-sqrt(fabs(lam[i+1])/massunit)/constant;
  }
  /*
  printf("frequencies in cm-1\n");
  for(i=0;i<3*ATOMS;i++){
    printf("%lf   ",freq[i]);
  }
    
  printf("\n");
  */    
  return freq;
}



double **invert_inertia_matrix( double **matrix, int dim, double **inv_matrix){

  int i, j;
  double det, **b;
  double **cc;

  //calcolo il determinante della matrice
  det=matrix[1][1]*(matrix[2][2]*matrix[3][3]-matrix[2][3]*matrix[3][2]);
  det=det - matrix[1][2]*(matrix[2][1]*matrix[3][3]-matrix[2][3]*matrix[3][1]);
  det=det + matrix[1][3]*(matrix[2][1]*matrix[3][2]-matrix[2][2]*matrix[3][1]);


  printf("det is  %lf  \n ",det);

  
   b=(double **) malloc( dim*sizeof(double *) );
  for(j=0; j < dim; j++) {
    b[j]= (double*) malloc( dim*sizeof(double) );
  }
   cc=(double **) malloc( dim*sizeof(double *) );
  for(j=0; j < dim; j++) {
    cc[j]= (double*) malloc( dim*sizeof(double) );
  }

  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      cc[i][j]=matrix[j][i];
    }
  }


  b[1][1]=(cc[2][2]*cc[3][3]-cc[2][3]*cc[3][2]);

  b[1][2]=-(cc[2][1]*cc[3][3]-cc[2][3]*cc[3][1]);

  b[1][3]=(cc[2][1]*cc[3][2]-cc[2][2]*cc[3][1]);
  

  b[2][1]=-(cc[1][2]*cc[3][3]-cc[1][3]*cc[3][2]);

  b[2][2]=(cc[1][1]*cc[3][3]-cc[1][3]*cc[3][1]);

  b[2][3]=-(cc[1][1]*cc[3][2]-cc[1][2]*cc[3][1]);

  
  b[3][1]=(cc[1][2]*cc[2][3]-cc[1][3]*cc[2][2]);

  b[3][2]=-(cc[1][1]*cc[2][3]-cc[1][3]*cc[2][1]);

  b[3][3]=(cc[1][1]*cc[2][2]-cc[1][2]*cc[2][1]);


  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){

      inv_matrix[i][j]=((b[i][j])/det);

    }
  }
  for(i=0;i<dim;i++){ free (b[i]); } free(b);
  for(i=0;i<dim;i++){ free (cc[i]); } free(cc);

  //  for(i=0;i<dim;i++)free(b[i]);free(b);
  return inv_matrix;
}


void determine_top_atoms(int step){

  int i,j;
  int i1,i2;

  //  int **igroupB;
  double **dist_matrix;
  //  int ncountA,ncountB;
  int nvar;
  double *ItopA,*ItopB;
  int *natomsA,*natomsB;
  double *massA,*massB;
  double masstot;
  double dist_CH,dist_HC,dist_CO,dist_CN;
  double dist_OH,dist_HO;
  double dist_OO;
  //  double dist_NO;
  double dist_CC;
  double dist_SH,dist_NH;
  double dist_HH;

  dist_CH=1.5;
  //  dist_CH=1.5;
  //  dist_HC=1.2;
  //  dist_CC=2.4;
  dist_CC=2.4;
  dist_CN=1.8;

  dist_CO=1.6;
  dist_OH=1.6;
  //dist_HO=1.2;
  dist_HO=1.5;
  dist_OO=1.6;
  dist_SH=1.5;
  dist_NH=1.3;
  //  dist_NO=1.6;
  //cc value
  dist_HC=1.25;
  //automech value
  // dist_HC=1.5;
  dist_HH=1.2;

  // check if distance file is present and open it

  FILE *dist_rotpr;
  FILE *rotmass;
  char buffer[200];
  char aname1[2];
  char aname2[2];
  int iread=1;
  int numdata=1000;
  double atdist=0.;
  //  char filename[20];

  if((dist_rotpr=fopen("dist_rotpr.dat","r"))==NULL) {
    printf("** using default distances to determine rotating moieties\n");
    iread=0;
  } else{
    printf("** updating distances of rotating moieties from file dist_rotpr.dat\n");
  }

  if((rotmass=fopen("rotmass.dat","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","trajec.txt");
    exit(1);
  }


  if(iread==1){

    fgets(buffer,sizeof(buffer),dist_rotpr);
    fgets(buffer,sizeof(buffer),dist_rotpr);
    //    fscanf(dist_rotpr,"%s %d",buffer,&numdata);
    for(j=0; j < numdata; j++) {
      fscanf(dist_rotpr,"%s %s %lf",aname1,aname2,&atdist);
      if(ferror(dist_rotpr)  ||  feof(dist_rotpr)) break;
      if((strcmp(aname1,"C")==0)&&(strcmp(aname2,"H")==0)){dist_CH=atdist;}
      if((strcmp(aname1,"H")==0)&&(strcmp(aname2,"C")==0)){dist_HC=atdist;}
      if((strcmp(aname1,"C")==0)&&(strcmp(aname2,"C")==0)){dist_CC=atdist;}
      if((strcmp(aname1,"C")==0)&&(strcmp(aname2,"N")==0)){dist_CN=atdist;}
      if((strcmp(aname1,"C")==0)&&(strcmp(aname2,"O")==0)){dist_CO=atdist;}
      if((strcmp(aname1,"O")==0)&&(strcmp(aname2,"H")==0)){dist_OH=atdist;}
      if((strcmp(aname1,"H")==0)&&(strcmp(aname2,"O")==0)){dist_HO=atdist;}
      if((strcmp(aname1,"O")==0)&&(strcmp(aname2,"O")==0)){dist_OO=atdist;}
      if((strcmp(aname1,"S")==0)&&(strcmp(aname2,"H")==0)){dist_SH=atdist;}
      if((strcmp(aname1,"N")==0)&&(strcmp(aname2,"H")==0)){dist_NH=atdist;}
      if((strcmp(aname1,"H")==0)&&(strcmp(aname2,"H")==0)){dist_HH=atdist;}
    }

    //    printf("%d numdata\n",numdata);
    //    printf("%s aname2\n",aname2);
    //    printf("%lf dist_NH\n",dist_NH);
    //    exit(0);
  }


  ItopA = (double *) malloc(numrotors*sizeof(double));
  ItopB = (double *) malloc(numrotors*sizeof(double));

  natomsA = (int *) malloc(numrotors*sizeof(int));
  natomsB = (int *) malloc(numrotors*sizeof(int));

  massA = (double *) malloc(numrotors*sizeof(double));
  massB = (double *) malloc(numrotors*sizeof(double));

  igroupB = (int **) malloc(numrotors*sizeof(int *));
  for(j=0; j < numrotors; j++) {
    igroupB[j]= (int*) malloc(ATOMS*sizeof(int));
  }

  //  int **igroupA;
  igroupA = (int **) malloc(numrotors*sizeof(int *));
  for(j=0; j < numrotors; j++) {
    igroupA[j]= (int*) malloc(ATOMS*sizeof(int));
  }


  dist_matrix = (double **) malloc(ATOMS*sizeof(double *));
  for(j=0; j < ATOMS; j++) {
    dist_matrix[j]= (double*) malloc(ATOMS*sizeof(double));
  }


  for (i=0;i<ATOMS;i++) {
    for (j=0;j<ATOMS;j++) {
      dist_matrix[i][j]=sqrt(pow((atoms_data[i].ext_coord_x[step]-atoms_data[j].ext_coord_x[step]),2)+
			     pow((atoms_data[i].ext_coord_y[step]-atoms_data[j].ext_coord_y[step]),2)+
			     pow((atoms_data[i].ext_coord_z[step]-atoms_data[j].ext_coord_z[step]),2));
    }
  }

  // convert to Angstrom


  for (i=0;i<ATOMS;i++){
    for (j=0;j<ATOMS;j++){
      dist_matrix[i][j]=dist_matrix[i][j]/bohr;
	//      printf("distance between at %d and atom %d is %lf  \n",i,j,dist_matrix[i][j]); 
    }
  }

  ///bohr


  for (i=0;i<numrotors;i++) {
    for (j=0;j<ATOMS;j++) {
      igroupB[i][j]=0;
      if (j==pivotB[i]-1){
	igroupB[i][j]=1;
      }
    }
  }
  
  for (i=0;i<numrotors;i++) {

    //    ncountB=1;
    nvar=1;
    int ncycle=0;
    while (nvar!=0 && ncycle < 1000){
      nvar=0;
      for (i1=0;i1<ATOMS;i1++) {
	if(igroupB[i][i1]==1){
	  for (i2=0;i2<ATOMS;i2++) {
	    if(i2!=pivotA[i]-1){
	      //check if C-H bond
	      // printf("at i1 %d is %s and atom i2 %d is  %s and dist is  %lf  \n",i1,atoms_data[i1].atom_name,i2,atoms_data[i2].atom_name,dist_matrix[i1][i2]); 
	      
	      if((strcmp(atoms_data[i1].atom_name,"C")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_CH)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		  igroupB[i][i2]=1;
		  nvar=nvar+1;
		}
	      //check if H-C bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"C")==0)&&(dist_matrix[i1][i2]<dist_HC)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		  igroupB[i][i2]=1;
		  nvar=nvar+1;
		}
	      //check if C-C bond
	      if((strcmp(atoms_data[i1].atom_name,"C")==0)&&(strcmp(atoms_data[i2].atom_name,"C")==0)&&(dist_matrix[i1][i2]<dist_CC)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		  igroupB[i][i2]=1;
		  nvar=nvar+1;
		}
	      //check if C-O bond
	      if((strcmp(atoms_data[i1].atom_name,"C")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_CO)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		  igroupB[i][i2]=1;
		  nvar=nvar+1;
		}
	      //check if O-C bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"C")==0)&&(dist_matrix[i1][i2]<dist_CO)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if C-N bond
	      if((strcmp(atoms_data[i1].atom_name,"C")==0)&&(strcmp(atoms_data[i2].atom_name,"N")==0)&&(dist_matrix[i1][i2]<dist_CN)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		  igroupB[i][i2]=1;
		  nvar=nvar+1;
		}
	      //check if N-C bond
	      if((strcmp(atoms_data[i1].atom_name,"N")==0)&&(strcmp(atoms_data[i2].atom_name,"C")==0)&&(dist_matrix[i1][i2]<dist_CN)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if O-O bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_OO)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if O-O bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_OO)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if H-O bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_HO)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if O-H bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_OH)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if H-S bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"S")==0)&&(dist_matrix[i1][i2]<dist_SH)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if S-H bond
	      if((strcmp(atoms_data[i1].atom_name,"S")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_SH)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if H-N bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"N")==0)&&(dist_matrix[i1][i2]<dist_NH)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if N-H bond
	      if((strcmp(atoms_data[i1].atom_name,"N")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_NH)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if N-O bond
	      if((strcmp(atoms_data[i1].atom_name,"N")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_CO)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		  igroupB[i][i2]=1;
		  nvar=nvar+1;
		}
	      //check if O-N bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"N")==0)&&(dist_matrix[i1][i2]<dist_CO)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if H-H bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_HH)&&(igroupB[i][i2]!=1)&&(atomsintopB[i][i2]==1)){
		igroupB[i][i2]=1;
		nvar=nvar+1;
	      }
	    }
	  }
	}
      }
      ncycle+=1;
    }
    //    exit(0);
  }  


  for (i=0;i<numrotors;i++) {
    for (j=0;j<ATOMS;j++) {
      igroupA[i][j]=0;
      if (j==pivotA[i]-1){
	igroupA[i][j]=1;
      }
    }
  }
  
  for (i=0;i<numrotors;i++) {

    //    ncountA=1;
    nvar=1;
    int ncycle=0;
    while (nvar!=0 && ncycle < 1000){
      nvar=0;
      for (i1=0;i1<ATOMS;i1++) {
	if(igroupA[i][i1]==1){
	  for (i2=0;i2<ATOMS;i2++) {
	    if(i2!=pivotB[i]-1){
	      //check if C-H bond
	      // printf("at i1 %d is %s and atom i2 %d is  %s and dist is  %lf  \n",i1,atoms_data[i1].atom_name,i2,atoms_data[i2].atom_name,dist_matrix[i1][i2]); 
	      
	      if((strcmp(atoms_data[i1].atom_name,"C")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_CH)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if H-C bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"C")==0)&&(dist_matrix[i1][i2]<dist_HC)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if C-C bond
	      if((strcmp(atoms_data[i1].atom_name,"C")==0)&&(strcmp(atoms_data[i2].atom_name,"C")==0)&&(dist_matrix[i1][i2]<dist_CC)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		  igroupA[i][i2]=1;
		  nvar=nvar+1;
		}
	      //check if C-O bond
	      if((strcmp(atoms_data[i1].atom_name,"C")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_CO)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if O-C bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"C")==0)&&(dist_matrix[i1][i2]<dist_CO)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if C-N bond
	      if((strcmp(atoms_data[i1].atom_name,"C")==0)&&(strcmp(atoms_data[i2].atom_name,"N")==0)&&(dist_matrix[i1][i2]<dist_CN)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if N-C bond
	      if((strcmp(atoms_data[i1].atom_name,"N")==0)&&(strcmp(atoms_data[i2].atom_name,"C")==0)&&(dist_matrix[i1][i2]<dist_CN)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if O-O bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_OO)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if O-O bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_OO)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if H-O bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_HO)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if O-H bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_OH)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if H-S bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"S")==0)&&(dist_matrix[i1][i2]<dist_SH)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if S-H bond
	      if((strcmp(atoms_data[i1].atom_name,"S")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_SH)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if H-N bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"N")==0)&&(dist_matrix[i1][i2]<dist_NH)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if N-H bond
	      if((strcmp(atoms_data[i1].atom_name,"N")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_NH)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if O-N bond
	      if((strcmp(atoms_data[i1].atom_name,"N")==0)&&(strcmp(atoms_data[i2].atom_name,"O")==0)&&(dist_matrix[i1][i2]<dist_CO)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
		}
	      //check if N-O bond
	      if((strcmp(atoms_data[i1].atom_name,"O")==0)&&(strcmp(atoms_data[i2].atom_name,"N")==0)&&(dist_matrix[i1][i2]<dist_CO)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
	      }
	      //check if H-H bond
	      if((strcmp(atoms_data[i1].atom_name,"H")==0)&&(strcmp(atoms_data[i2].atom_name,"H")==0)&&(dist_matrix[i1][i2]<dist_HH)&&(igroupA[i][i2]!=1)&&(atomsintopA[i][i2]==1)){
		igroupA[i][i2]=1;
		nvar=nvar+1;
	      }
	    }
	  }
	}
      }
      ncycle+=1;
    }
    //  for (i=0;i<numrotors;i++){
    //    for (j=0;j<ATOMS;j++){
    //        printf("igroup A valueb for rot %d natom %d is %d  \n",i,j,igroupA[i][j]); 
	    //	    printf("igroup B value for rot %d natom %d is %d  \n",i,j,igroupB[i][j]); 
    //}
    //  }

    //    exit(0);
  }  

  /*
  for (i=0;i<ATOMS;i++){
    for (j=0;j<ATOMS;j++){
      printf("distance between at %d and atom %d is %lf  \n",i,j,dist_matrix[i][j]); 
    }
  }
  */

  // determine inertia momentts for topA and topB

  double IxxT,IyyT,IzzT;
  double rotaxx,rotaxy,rotaxz;
  double rotatx,rotaty,rotatz;
  double distx,disty,distz;
  double norm;
  //  double distest;


  for (i=0; i<numrotors; i++) {
    rotaxx=(atoms_data[pivotA[i]-1].ext_coord_x[step]-atoms_data[pivotB[i]-1].ext_coord_x[step])/bohr;  
    rotaxy=(atoms_data[pivotA[i]-1].ext_coord_y[step]-atoms_data[pivotB[i]-1].ext_coord_y[step])/bohr;  
    rotaxz=(atoms_data[pivotA[i]-1].ext_coord_z[step]-atoms_data[pivotB[i]-1].ext_coord_z[step])/bohr;  
    norm=dist_matrix[pivotA[i]-1][pivotB[i]-1];
    IxxT=0;
    IyyT=0;
    IzzT=0;
    for (j=0; j<ATOMS; j++) {
      if(igroupA[i][j]==1){
	rotatx=(atoms_data[j].ext_coord_x[step]-atoms_data[pivotA[i]-1].ext_coord_x[step])/bohr;  
	rotaty=(atoms_data[j].ext_coord_y[step]-atoms_data[pivotA[i]-1].ext_coord_y[step])/bohr;  
	rotatz=(atoms_data[j].ext_coord_z[step]-atoms_data[pivotA[i]-1].ext_coord_z[step])/bohr;  

	cross_prod(rotatx,rotaty,rotatz,rotaxx,rotaxy,rotaxz,&distx,&disty,&distz);
	IxxT=IxxT+distx*distx*(atoms_data[j].atomic_mass)/norm/norm;
	IyyT=IyyT+disty*disty*(atoms_data[j].atomic_mass)/norm/norm;
	IzzT=IzzT+distz*distz*(atoms_data[j].atomic_mass)/norm/norm;
	//	distest=sqrt((distx*distx+disty*disty+distz*distz))/norm;
	//	printf("j is %d norm is %lf and dist is %lf \n",j,norm,distest);
      }
    }
    ItopA[i]=IxxT+IyyT+IzzT;
    //    printf("ItopA of rotor %d is  %lf \n",i,ItopA[i]);
    //    exit(0);
  }

  // now compute the same for the B top


  for (i=0; i<numrotors; i++) {
    rotaxx=(atoms_data[pivotB[i]-1].ext_coord_x[step]-atoms_data[pivotA[i]-1].ext_coord_x[step])/bohr;  
    rotaxy=(atoms_data[pivotB[i]-1].ext_coord_y[step]-atoms_data[pivotA[i]-1].ext_coord_y[step])/bohr;  
    rotaxz=(atoms_data[pivotB[i]-1].ext_coord_z[step]-atoms_data[pivotA[i]-1].ext_coord_z[step])/bohr;  
    norm=dist_matrix[pivotB[i]-1][pivotA[i]-1];
    IxxT=0;
    IyyT=0;
    IzzT=0;
    for (j=0; j<ATOMS; j++) {
      if(igroupB[i][j]==1){
	rotatx=(atoms_data[j].ext_coord_x[step]-atoms_data[pivotB[i]-1].ext_coord_x[step])/bohr;  
	rotaty=(atoms_data[j].ext_coord_y[step]-atoms_data[pivotB[i]-1].ext_coord_y[step])/bohr;  
	rotatz=(atoms_data[j].ext_coord_z[step]-atoms_data[pivotB[i]-1].ext_coord_z[step])/bohr;  

	cross_prod(rotatx,rotaty,rotatz,rotaxx,rotaxy,rotaxz,&distx,&disty,&distz);
	IxxT=IxxT+distx*distx*(atoms_data[j].atomic_mass)/norm/norm;
	IyyT=IyyT+disty*disty*(atoms_data[j].atomic_mass)/norm/norm;
	IzzT=IzzT+distz*distz*(atoms_data[j].atomic_mass)/norm/norm;
	//	distest=sqrt((distx*distx+disty*disty+distz*distz))/norm;
	//	printf("j is %d norm is %lf and dist is %lf \n",j,norm,distest);
      }
    }
    ItopB[i]=IxxT+IyyT+IzzT;
    //    printf("ItopB of rotor %d is  %lf \n",i,ItopB[i]);
    //    exit(0);
  }

  /*  
  for (i=0;i<numrotors;i++){
    for (j=0;j<ATOMS;j++){
            printf("igroup A value for rot %d natom %d is %d  \n",i,j,igroupA[i][j]); 
	    //	    printf("igroup B value for rot %d natom %d is %d  \n",i,j,igroupB[i][j]); 
    }
  }
  for (i=0;i<numrotors;i++){
    for (j=0;j<ATOMS;j++){
	    printf("igroup B value for rot %d natom %d is %d  \n",i,j,igroupB[i][j]); 
    }
  }
  */
  //  exit(0);
  

  for (i=0;i<numrotors;i++){
    natomsA[i]=0;
    natomsB[i]=0;
    massA[i]=0.;
    massB[i]=0.;
    masstot=0.;
    for (j=0;j<ATOMS;j++){
      if(igroupA[i][j] == 1 ){
	natomsA[i]+=1;
	massA[i]+=atoms_data[j].atomic_mass;

      } 
      if(igroupB[i][j] == 1 ){
	natomsB[i]+=1;
	massB[i]+=atoms_data[j].atomic_mass;
      } 
      masstot+=atoms_data[j].atomic_mass;
      // fprintf(rotmass," %d %lf  %lf  %lf  \n",i,masstot,massA[i],massB[i]);
    }
    fprintf(rotmass," %d %lf  %lf  %lf  \n",i,masstot,massA[i],massB[i]);
  }

  fclose(rotmass);

  // of the two tops identified I take the one with the smallest inertia moment    
  int indred;

  indred=0;
  
  for (i=0;i<numrotors;i++){
    if(natomsA[i] < numatomsintopA[i]){indred=1;}
    if(natomsB[i] < ATOMS-1-numatomsintopA[i]){indred=2;}
    indred=0;
    for (j=0;j<ATOMS;j++){
      //      printf("igroup B value for rot %d natom %d is %d  \n",i,j,igroupB[i][j]); 
      if(ItopB[i]==0 ){
	igroupB[i][j]=0;
      } 
      else if(ItopA[i] < ItopB[i] && indred == 0){
	igroupB[i][j]=0;
      } 
      else if(ItopA[i] < ItopB[i] && indred == 1){
	igroupB[i][j]=0;
      } 
      else if(ItopA[i] < ItopB[i] && indred == 2){
	igroupA[i][j]=0;
      } 
      else {
	igroupA[i][j]=0;
      }
    }
  }

  /*          
  for (i=0;i<numrotors;i++) {
    for (j=0;j<ATOMS;j++) {
      igroupA[i][j]=0;
      if (j==pivotA[i]-1){
	igroupA[i][j]=1;
      }
      if (atomsintopA[i][j]==1){
	igroupA[i][j]=1;
      }
    }
  }

  for (i=0;i<numrotors;i++) {
    for (j=0;j<ATOMS;j++) {
      igroupB[i][j]=0;
      if (j==pivotB[i]-1){
	igroupB[i][j]=1;
      }
      if (atomsintopB[i][j]==1){
	igroupB[i][j]=1;
      }
    }
  }

  */
 
        
  for (i=0;i<numrotors;i++){
    for (j=0;j<ATOMS;j++){
            printf("igroup A value for rot %d natom %d is %d  \n",i,j,igroupA[i][j]); 
	    //	    printf("igroup B value for rot %d natom %d is %d  \n",i,j,igroupB[i][j]); 
    }
  }
  for (i=0;i<numrotors;i++){
    for (j=0;j<ATOMS;j++){
	    printf("igroup B value for rot %d natom %d is %d  \n",i,j,igroupB[i][j]); 
    }
  } 
  
  free(natomsA);
  free(natomsB);
  free(massA);
  free(massB);
  //    exit(0);
  

}


void projector_matrix_Rot(int step){

  //Projector matrix implemented as described by Green in Cantherm to project rotors

  int i,j,k,ik;
  double **I_temp;

  // set coordinates to principal axes NB working in Bohr

  //reverse v (inertia eigenvector) in order to be consistent with cantherm implementation (from smallest to highest eigen)

  double tempx,tempy,tempz;
  
  tempx=v[step][1][1];
  tempy=v[step][2][1];
  tempz=v[step][3][1];
  v[step][1][1]=v[step][1][3];
  v[step][2][1]=v[step][2][3];
  v[step][3][1]=v[step][3][3];
  v[step][1][3]=tempx;
  v[step][2][3]=tempy;
  v[step][3][3]=tempz;

  // save coordinates and convert to angs

  double **coords;
    
  coords=(double **) malloc( (ATOMS)*sizeof(double *) );
  for(j=0; j < ATOMS; j++) {
    coords[j]= (double*) malloc( (3)*sizeof(double) );
  }
  for(i=0; i < ATOMS; i++) {
    coords[i][0]=atoms_data[i].ext_coord_x[step]/bohr;
    coords[i][1]=atoms_data[i].ext_coord_y[step]/bohr;
    coords[i][2]=atoms_data[i].ext_coord_z[step]/bohr;
  }

  
    coordinates_in_principal_axes(step);


  // convert coord to ang

  for (j=0;j<ATOMS;j++){
    atoms_data[j].ext_coord_x[step]=atoms_data[j].ext_coord_x[step]/bohr;  
    atoms_data[j].ext_coord_y[step]=atoms_data[j].ext_coord_y[step]/bohr;  
    atoms_data[j].ext_coord_z[step]=atoms_data[j].ext_coord_z[step]/bohr;  
  }

  /*
  for (j=0;j<ATOMS;j++){
    printf("atomic coo are: %lf %lf %lf \n",atoms_data[j].ext_coord_x[step],atoms_data[j].ext_coord_y[step],atoms_data[j].ext_coord_z[step]); 
  }
  */

  //construct projection matrix
  int external=6;
  // cc tests

  if(proj_rea_coo == 1){ external=6;}
  if(proj_rea_coo == 0){ external=7;}

  //  int external=6;
  //  int external=7;

  double **Dmat;
    
  Dmat=(double **) malloc( (3*ATOMS)*sizeof(double *) );
  for(j=0; j < 3*ATOMS; j++) {
    Dmat[j]= (double*) malloc( (external)*sizeof(double) );
  }
  
  for(i=0; i < 3*ATOMS; i++) {
    for(j=0; j < external; j++) {
      Dmat[i][j]=0.;
    }
  }


  //    Dmat=(double **) malloc( (external-1)*sizeof(double *) );
  //    for(j=0; j < external-1; j++) {
  //      Dmat[j]= (double*) malloc( (3*ATOMS-1)*sizeof(double) );
  //    }

  for (i=0;i<ATOMS;i++){
    Dmat[3*i+0][0]=sqrt(atoms_data[i].atomic_mass);
    Dmat[3*i+1][1]=sqrt(atoms_data[i].atomic_mass);
    Dmat[3*i+2][2]=sqrt(atoms_data[i].atomic_mass);
  }
  
  for (i=0;i<ATOMS;i++){

    //Dmat[3*i][3]=(atoms_data[i].ext_coord_z[step]*v[step][1][1]-atoms_data[i].ext_coord_x[step]*v[step][1][3])*sqrt(atoms_data[i].atomic_mass);
    //Dmat[3*i+1][3]=(atoms_data[i].ext_coord_z[step]*v[step][2][1]-atoms_data[i].ext_coord_x[step]*v[step][2][3])*sqrt(atoms_data[i].atomic_mass);
    //Dmat[3*i+2][3]=(atoms_data[i].ext_coord_z[step]*v[step][3][1]-atoms_data[i].ext_coord_x[step]*v[step][3][3])*sqrt(atoms_data[i].atomic_mass);

    //Dmat[3*i][4]=(atoms_data[i].ext_coord_x[step]*v[step][1][2]-atoms_data[i].ext_coord_y[step]*v[step][1][1])*sqrt(atoms_data[i].atomic_mass);
    //Dmat[3*i+1][4]=(atoms_data[i].ext_coord_x[step]*v[step][2][2]-atoms_data[i].ext_coord_y[step]*v[step][2][1])*sqrt(atoms_data[i].atomic_mass);
    //Dmat[3*i+2][4]=(atoms_data[i].ext_coord_x[step]*v[step][3][2]-atoms_data[i].ext_coord_y[step]*v[step][3][1])*sqrt(atoms_data[i].atomic_mass);


    Dmat[3*i][3]=(atoms_data[i].ext_coord_y[step]*v[step][1][3]-atoms_data[i].ext_coord_z[step]*v[step][1][2])*sqrt(atoms_data[i].atomic_mass);
    Dmat[3*i+1][3]=(atoms_data[i].ext_coord_y[step]*v[step][2][3]-atoms_data[i].ext_coord_z[step]*v[step][2][2])*sqrt(atoms_data[i].atomic_mass);
    Dmat[3*i+2][3]=(atoms_data[i].ext_coord_y[step]*v[step][3][3]-atoms_data[i].ext_coord_z[step]*v[step][3][2])*sqrt(atoms_data[i].atomic_mass);

    Dmat[3*i][4]=(atoms_data[i].ext_coord_z[step]*v[step][1][1]-atoms_data[i].ext_coord_x[step]*v[step][1][3])*sqrt(atoms_data[i].atomic_mass);
    Dmat[3*i+1][4]=(atoms_data[i].ext_coord_z[step]*v[step][2][1]-atoms_data[i].ext_coord_x[step]*v[step][2][3])*sqrt(atoms_data[i].atomic_mass);
    Dmat[3*i+2][4]=(atoms_data[i].ext_coord_z[step]*v[step][3][1]-atoms_data[i].ext_coord_x[step]*v[step][3][3])*sqrt(atoms_data[i].atomic_mass);

    Dmat[3*i][5]=(atoms_data[i].ext_coord_x[step]*v[step][1][2]-atoms_data[i].ext_coord_y[step]*v[step][1][1])*sqrt(atoms_data[i].atomic_mass);
    Dmat[3*i+1][5]=(atoms_data[i].ext_coord_x[step]*v[step][2][2]-atoms_data[i].ext_coord_y[step]*v[step][2][1])*sqrt(atoms_data[i].atomic_mass);
    Dmat[3*i+2][5]=(atoms_data[i].ext_coord_x[step]*v[step][3][2]-atoms_data[i].ext_coord_y[step]*v[step][3][1])*sqrt(atoms_data[i].atomic_mass);

    //cc test projection along reaction coordinate
    if(external ==7){
      Dmat[3*i][6]=atoms_data[i].grad_x[step]/sqrt(atoms_data[i].atomic_mass);
      Dmat[3*i+1][6]=atoms_data[i].grad_y[step]/sqrt(atoms_data[i].atomic_mass);
      Dmat[3*i+2][6]=atoms_data[i].grad_z[step]/sqrt(atoms_data[i].atomic_mass);
    }
    //

  }

  // check rotational projection vectors module for degeneracy
  // in case of degeneracy remove projection and renumber external

  double check3=0.;
  double check4=0.;
  double check5=0.;
  int skiprot=0;

  for (i=0;i<ATOMS;i++){
    check3=check3+Dmat[3*i][3]*Dmat[3*i][3]+Dmat[3*i+1][3]*Dmat[3*i+1][3]+Dmat[3*i+2][3]*Dmat[3*i+2][3];
    check4=check4+Dmat[3*i][4]*Dmat[3*i][4]+Dmat[3*i+1][4]*Dmat[3*i+1][4]+Dmat[3*i+2][4]*Dmat[3*i+2][4];
    check5=check5+Dmat[3*i][5]*Dmat[3*i][5]+Dmat[3*i+1][5]*Dmat[3*i+1][5]+Dmat[3*i+2][5]*Dmat[3*i+2][5];
  }

  //  printf("check 3 is %lf\n",check3);
  //  printf("check 4 is %lf\n",check4);
  //  printf("check 5 is %lf\n",check5);

  if (check3<1.0E-10){
    for (i=0;i<ATOMS;i++){
      Dmat[3*i][3]= Dmat[3*i][4];
      Dmat[3*i+1][3]= Dmat[3*i+1][4];
      Dmat[3*i+2][3]= Dmat[3*i+2][4];
      Dmat[3*i][4]= Dmat[3*i][5];
      Dmat[3*i+1][4]= Dmat[3*i+1][5];
      Dmat[3*i+2][4]= Dmat[3*i+2][5];
      skiprot=1;
    }
    external=5;
  }

  if (check4<1.0E-10){
    for (i=0;i<ATOMS;i++){
      Dmat[3*i][4]= Dmat[3*i][5];
      Dmat[3*i+1][4]= Dmat[3*i+1][5];
      Dmat[3*i+2][4]= Dmat[3*i+2][5];
    }
    external=5;
    skiprot=1;
  }

  if (check5<1.0E-10){
    external=5;
    skiprot=1;
  }

  // now construct the projection matrix

  double **Pmat;
    
  Pmat=(double **) malloc( (3*ATOMS)*sizeof(double *) );
  for(j=0; j < 3*ATOMS; j++) {
    Pmat[j]= (double*) malloc( (external+3*ATOMS)*sizeof(double) );
  }
  
  for(i=0; i < 3*ATOMS; i++) {
    for(j=0; j < external+3*ATOMS; j++) {
      Pmat[i][j]=0.;
    }
  }


  double **Identity;
 
  Identity=(double **) malloc( (int)(3*ATOMS)*sizeof(double *) );
  for(j=0; j < 3*ATOMS; j++) {
    Identity[j]= (double*) malloc( (int)(3*ATOMS)*sizeof(double) );
  }

  for(i=0;i<3*ATOMS;i++){
    for(j=0;j<3*ATOMS;j++){
      if(i==j) Identity[i][j]=1.0;
      else Identity[i][j]=0.0;
    }
  }
  
  for(i=0; i < 3*ATOMS; i++) {
    for(j=0; j < external+3*ATOMS; j++) {
      if(j<external){
	Pmat[i][j]=Dmat[i][j];
      }
      if(j>=external){
	Pmat[i][j]=Identity[i][j-external];
      }
    }
  }

  //now orthonormalize matrix


  
  double norm;
  double proj;
  
  for (i=0;i<3*ATOMS+external; i++) {
    norm=0.;
    for (j=0;j<3*ATOMS; j++) {
      norm+= Pmat[j][i]*Pmat[j][i];
    }
    for (j=0;j< 3*ATOMS;j++) {
      if (norm>1.0E-15){
	Pmat[j][i]=Pmat[j][i]/sqrt(norm);
      }
      else{
	Pmat[j][i]=0.;
      }
    }

    for (j=i+1;j<3*ATOMS+external;j++){
      proj=0.0;
      for (k=0;k<3*ATOMS;k++){
	proj+=Pmat[k][i]*Pmat[k][j];
      }
      for (k=0;k<3*ATOMS;k++){
	Pmat[k][j]-=proj*Pmat[k][i]; 
      }
    }
  }


  //now order vectors

  i=0;
  while(i<3*ATOMS){
    norm=0.;
    for (j=0;j< 3*ATOMS;j++) {
      norm+= Pmat[j][i]*Pmat[j][i];
    }
    if (norm<0.5){
      for(ik=i;ik<3*ATOMS+external;ik++){
	for(k=0;k<3*ATOMS;k++){
	  Pmat[k][ik]=Pmat[k][ik+1];
	}
      }
    }
    else{
      i=i+1;
    }
  }
  
  // now extract the projection matrix

  double **Tmat;
    
  Tmat=(double **) malloc( (3*ATOMS)*sizeof(double *) );
  for(j=0; j < 3*ATOMS; j++) {
    Tmat[j]= (double*) malloc( (3*ATOMS-external)*sizeof(double) );
  }
  
  for(i=0; i < 3*ATOMS; i++) {
    for(j=0; j < 3*ATOMS-external; j++) {
      Tmat[i][j]=Pmat[i][j+external];
    }
  }


  // now convert the Hessian to J/m^2*kg

  double conv;

  conv=joule*bohr*bohr/a_mass*1.0e10*1.0e10;
  
  for (j=0;j<3*ATOMS;j++){
    for (k=0;k<3*ATOMS;k++){
      Hessian[step][j][k]=Hessian[step][j][k]*conv; 
    }
  }


  // now project the Hessian and construct the internal force constant matrix

  //double **I_temp;
  
  //  I_temp=(double **) malloc( (3*ATOMS)*sizeof(double *) );
  //  for(j=0; j < 3*ATOMS; j++) {
  //    I_temp[j]= (double*) malloc( (3*ATOMS-external)*sizeof(double) );
  //  }
  
  I_temp=prod_mat2(Hessian[step],3*ATOMS,3*ATOMS,Tmat,3*ATOMS,3*ATOMS-external);

  double **TmatT;
    
  TmatT=(double **) malloc( (3*ATOMS-external)*sizeof(double *) );
  for(j=0; j < 3*ATOMS-external; j++) {
    TmatT[j]= (double*) malloc( (3*ATOMS)*sizeof(double) );
  }

  TmatT=transpose_double(Tmat,3*ATOMS-external,3*ATOMS,TmatT);
  
  double **Fc_int;
  
  Fc_int=prod_mat2(TmatT,3*ATOMS-external,3*ATOMS,I_temp,3*ATOMS,3*ATOMS-external);

  
  // now calculate eigenvalues

  int dim;

  dim=3*ATOMS-external+1;

  double *lambda;
  double **L;
  double **Lt;
  double **FC_temp;
  
  lambda = (double *) malloc(sizeof(double) * dim);
  L = (double **) malloc((int)dim*sizeof(double *) );
  Lt = (double **) malloc((int)dim*sizeof(double *) );
  FC_temp = (double **) malloc((int)dim*sizeof(double *) );
 
  for(j=0; j<(int)(dim); j++) {
    FC_temp[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
    L[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
    Lt[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
  }

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      if(j==i) L[i][j]=1;
      else L[i][j]=0.;
      FC_temp[i][j]=0.;
    }
  }
  
  for (i=1; i<dim; i++) {
    for (j=1; j<dim; j++) {
      FC_temp[i][j]= Fc_int[i-1][j-1];
    }
  }

  //  exit(0);
  
  //  printf("Hessian ibefore diag 1,1 %lf\n",FC_temp[1][1]);

  jacobi(FC_temp, dim-1, lambda, L, nrot);
  eigsrt(lambda, L, dim-1);

  for (i=0; i<dim-1; i++) {
    for (j=0; j<dim-1; j++) {
      Lt[i][j]=L[i+1][j+1];
    }
  }


  double **TLmat;
  TLmat=prod_mat2(Tmat,3*ATOMS,3*ATOMS-external,Lt,3*ATOMS-external,3*ATOMS-external);


  //  printf("6D Projected FC diagonalization\n Eigenvalues --- \n");
  //  for(j=1; j<(int)(dim); j++) {    
  //    printf("%6.2f\n",sqrt(fabs(lambda[j]))/(2 * pi * c_light_cm_s ));
  //  }

  FILE *RTproj_freq;
  FILE *RTproj_cd;
  //  char filename[20];

  if((RTproj_freq=fopen("RTproj_freq.dat","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","RTproj_freq.dat");
    exit(1);
  }
  if((RTproj_cd=fopen("RTproj_cd.dat","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","RTproj_freq.dat");
    exit(1);
  }

  for(j=1;j<(int)(dim);j++){
    fprintf(RTproj_freq," %6.2f  \n",sqrt(fabs(lambda[j]))/(2 * pi * c_light_cm_s ));
  }

  fclose(RTproj_freq);

  save_displ(lambda,TLmat,0,dim-1,step,RTproj_cd);

  fclose(RTproj_cd);

  //  exit(0);

  //now project rotors
  //here we work with cartesian coordinates

  double **Dmatrot;
    
  Dmatrot=(double **) malloc( (3*ATOMS)*sizeof(double *) );
  for(j=0; j < 3*ATOMS; j++) {
    Dmatrot[j]= (double*) malloc( (numrotors)*sizeof(double) );
  }
  
  for(i=0; i < 3*ATOMS; i++) {
    for(j=0; j < numrotors; j++) {
      Dmatrot[i][j]=0.;
    }
  }

  /*
  int **igroupA;
  igroupA = (int **) malloc(numrotors*sizeof(int *));
  for(j=0; j < numrotors; j++) {
    igroupA[j]= (int*) malloc(ATOMS*sizeof(int));
  }

  for(i=0; i < numrotors; i++) {
    for(j=0; j < ATOMS; j++) {
      igroupA[i][j]=0.;
    }
  }
  
  for (i=0; i<numrotors; i++) {
    for(j=0; j < ATOMS; j++) {
      igroupA[i][j]=0;
      for(k=0; k< numatomsintopA[i]; k++){
	if(j==(atomsintopA[i][k]-1)||j==(pivotA[i]-1)){
	  igroupA[i][j]=1;
	}
      }
    }
  }
  */
  /*
  double *cmAx,*cmAy,*cmAz;
  double *cmBx,*cmBy,*cmBz;
  double masstopA;
  double masstopB;
  cmAx = (double *) malloc(sizeof(double) * numrotors);
  cmAy = (double *) malloc(sizeof(double) * numrotors);
  cmAz = (double *) malloc(sizeof(double) * numrotors);
  cmBx = (double *) malloc(sizeof(double) * numrotors);
  cmBy = (double *) malloc(sizeof(double) * numrotors);
  cmBz = (double *) malloc(sizeof(double) * numrotors);

  
  for (i=0; i<numrotors; i++) {
    cmAx[i]=0;
    cmAy[i]=0;
    cmAz[i]=0;
    cmBx[i]=0;
    cmBy[i]=0;
    cmBz[i]=0;
    masstopA=0;
    masstopB=0;

    for(j=0; j < ATOMS; j++) {
      if(igroupA[i][j]==1){
	cmAx[i]+=coords[j][0]*atoms_data[j].atomic_mass;
	cmAy[i]+=coords[j][1]*atoms_data[j].atomic_mass;
	cmAz[i]+=coords[j][2]*atoms_data[j].atomic_mass;
	masstopA+=atoms_data[j].atomic_mass;
	printf("i top A  is %d  \n",j);
      }
      else{
	cmBx[i]+=coords[j][0]*atoms_data[j].atomic_mass;
	cmBy[i]+=coords[j][1]*atoms_data[j].atomic_mass;
	cmBz[i]+=coords[j][2]*atoms_data[j].atomic_mass;
	masstopB+=atoms_data[j].atomic_mass;
	printf("i top B  is %d  \n",j);
      }
    }

    cmAx[i]=cmAx[i]/masstopA;
    cmAy[i]=cmAy[i]/masstopA;
    cmAz[i]=cmAz[i]/masstopA;
    cmBx[i]=cmBx[i]/masstopB;
    cmBy[i]=cmBy[i]/masstopB;
    cmBz[i]=cmBz[i]/masstopB;

    printf("CmA x is %lf pivot x is %lf \n",cmAx[i],coords[pivotA[i]-1][0]);
    printf("CmA y is %lf pivot y is %lf \n",cmAy[i],coords[pivotA[i]-1][1]);
    printf("CmA z is %lf pivot z is %lf \n",cmAz[i],coords[pivotA[i]-1][2]);
    printf("CmB x is %lf pivot x is %lf \n",cmBx[i],coords[pivotB[i]-1][0]);
    printf("CmB y is %lf pivot y is %lf \n",cmBy[i],coords[pivotB[i]-1][1]);
    printf("CmB z is %lf pivot z is %lf \n",cmBz[i],coords[pivotB[i]-1][2]);
    
  }
  */

  //  exit(0);
  
  // now change Pivot A and B to centers of mass of A and B
  
  double e12x,e12y,e12z;
  double e31x,e31y,e31z;
  double e32x,e32y,e32z;
  double p31x,p31y,p31z;
  double p32x,p32y,p32z;

  //  double *redA, *redB;
  //  redA = (double *) malloc(sizeof(double) * numrotors);
  //  redB = (double *) malloc(sizeof(double) * numrotors);

  
  for (i=0; i<numrotors; i++) {
    // e12x=cmAx[i]-cmBx[i];
    //e12y=cmAy[i]-cmBy[i];
    //e12z=cmAz[i]-cmBz[i];
    e12x=coords[pivotA[i]-1][0]-coords[pivotB[i]-1][0];
    e12y=coords[pivotA[i]-1][1]-coords[pivotB[i]-1][1];
    e12z=coords[pivotA[i]-1][2]-coords[pivotB[i]-1][2];
    //    redA[i]=0;
    //    redB[i]=0;
    
    for (j=0; j<ATOMS; j++) {
      if(igroupA[i][j]==1){
       	e31x=coords[j][0]-coords[pivotA[i]-1][0];
       	e31y=coords[j][1]-coords[pivotA[i]-1][1];
       	e31z=coords[j][2]-coords[pivotA[i]-1][2];

	//e31x=coords[j][0]-cmAx[i];
	//e31y=coords[j][1]-cmAy[i];
	//e31z=coords[j][2]-cmAz[i];
	cross_prod(e31x,e31y,e31z,e12x,e12y,e12z,&p31x,&p31y,&p31z);
	Dmatrot[3*j][i]=p31x*sqrt(atoms_data[j].atomic_mass);
	Dmatrot[3*j+1][i]=p31y*sqrt(atoms_data[j].atomic_mass);
	Dmatrot[3*j+2][i]=p31z*sqrt(atoms_data[j].atomic_mass);
	//	redA[i]+=Dmatrot[3*j][i]+Dmatrot[3*j+1][i]+Dmatrot[3*j+2][i];
      }
      else if(igroupB[i][j]==1) {
	e32x=coords[j][0]-coords[pivotB[i]-1][0];
       	e32y=coords[j][1]-coords[pivotB[i]-1][1];
	e32z=coords[j][2]-coords[pivotB[i]-1][2];
	//e32x=coords[j][0]-cmBx[i];
	//e32y=coords[j][1]-cmBy[i];
	//e32z=coords[j][2]-cmBz[i];
	cross_prod(e32x,e32y,e32z,e12x,e12y,e12z,&p32x,&p32y,&p32z);
	Dmatrot[3*j][i]=-p32x*sqrt(atoms_data[j].atomic_mass);
	Dmatrot[3*j+1][i]=-p32y*sqrt(atoms_data[j].atomic_mass);
	Dmatrot[3*j+2][i]=-p32z*sqrt(atoms_data[j].atomic_mass);
	//	redB[i]+=Dmatrot[3*j][i]+Dmatrot[3*j+1][i]+Dmatrot[3*j+2][i];
      }
    }
  }

  /*
  for (i=0; i<numrotors; i++) {
    for (j=0; j<ATOMS; j++) {
      if(igroupA[i][j]==1){
	Dmatrot[3*j][i]=Dmatrot[3*j][i]*redB[i]/(redA[i]+redB[i]);
	Dmatrot[3*j+1][i]=Dmatrot[3*j+1][i]*redB[i]/(redA[i]+redB[i]);
	Dmatrot[3*j+2][i]=Dmatrot[3*j+2][i]*redB[i]/(redA[i]+redB[i]);
      }
      else{
	Dmatrot[3*j][i]=Dmatrot[3*j][i]*redA[i]/(redA[i]+redB[i]);
	Dmatrot[3*j+1][i]=Dmatrot[3*j+1][i]*redA[i]/(redA[i]+redB[i]);
	Dmatrot[3*j+2][i]=Dmatrot[3*j+2][i]*redA[i]/(redA[i]+redB[i]);
      }
    }
  }
  */
  /*
  for(i=0; i<(int)(3*ATOMS); i++) {    
    for(j=0; j<(int)(numrotors); j++) {    
      printf("%le ",Dmatrot[i][j]);
    }
      printf("\n");
  }
  exit(0);
  */
  // now convert the normal modes back to cartesian coord. Use the rototransl projection matrix


  double **Vmwmat;
  double **VmwmatT;

  // first reverse the L and lambda order
  double *lambda_rev;
  double **L_rev;
  int irev;
  
  dim=3*ATOMS-external;

  lambda_rev = (double *) malloc(sizeof(double) * dim);
  L_rev = (double **) malloc((int)dim*sizeof(double *) );
  for(j=0; j<(int)(dim); j++) {
    L_rev[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
  }

  irev=dim;
  for(j=0; j<(int)(dim); j++) {
    lambda_rev[j]=lambda[irev];
    for(i=0; i<(int)(dim); i++) {
      L_rev[i][j]=L[i+1][irev];
    }
    irev=irev-1;
  }
  
  Vmwmat=prod_mat2(Tmat,3*ATOMS,3*ATOMS-external,L_rev,3*ATOMS-external,3*ATOMS-external);

  VmwmatT=(double **) calloc( (3*ATOMS-external) , sizeof(double *) );
  for(j=0; j < 3*ATOMS-external; j++) {
    VmwmatT[j]= (double*) calloc( 3*ATOMS, sizeof(double) );
  }

  VmwmatT=transpose_double(Vmwmat,3*ATOMS-external,3*ATOMS,VmwmatT);

  double **eigM;
  eigM=(double **) calloc( (3*ATOMS-external) , sizeof(double *) );
  for(j=0; j < 3*ATOMS-external; j++) {
    eigM[j]= (double*) calloc( 3*ATOMS-external, sizeof(double) );
  }

  for (i=0;i<3*ATOMS-external;i++){
    eigM[i][i]=lambda_rev[i];
  }

  double**Int1;

  Int1=prod_mat2(eigM,3*ATOMS-external,3*ATOMS-external,VmwmatT,3*ATOMS-external,3*ATOMS);
    
  double **FMI;

  FMI=prod_mat2(Vmwmat,3*ATOMS,3*ATOMS-external,Int1,3*ATOMS-external,3*ATOMS);


    //  projection of internal rotors on internal cartesion FC matrix

  double **Dprojmat;

  Dprojmat=prod_mat2(VmwmatT,3*ATOMS-external,3*ATOMS,Dmatrot,3*ATOMS,numrotors);
  
  for(i=0; i < numrotors; i++) {
    for(j=0; j < 3*ATOMS; j++) {
      Dmatrot[j][i]=0.;
      for(k=0; k < 3*ATOMS-external; k++) {
	Dmatrot[j][i]+=Dprojmat[k][i]*Vmwmat[j][k];
      }
    }
  }
    
  //now orthonormalize matrix

  for(i=0; i < numrotors; i++) {
    norm=0.;
    for(j=0; j < 3*ATOMS; j++) {
      norm+=Dmatrot[j][i]*Dmatrot[j][i];
    }
    for(j=0; j < 3*ATOMS; j++) {
      Dmatrot[j][i]=Dmatrot[j][i]/sqrt(norm);
    }
    for(j=i+1; j < numrotors; j++) {
      proj=0.;
      for(k=0; k < 3*ATOMS; k++) {
	proj+=Dmatrot[k][i]*Dmatrot[k][j];
      }
      for(k=0; k < 3*ATOMS; k++) {
	Dmatrot[k][j]-=proj*Dmatrot[k][i];
      }
    }
  }

  // now project FC matrix and calculate eigenvalues


  double **DmatrotT;


  DmatrotT=(double **) malloc( (numrotors)*sizeof(double *) );
  for(j=0; j < numrotors; j++) {
    DmatrotT[j]= (double*) malloc( (3*ATOMS)*sizeof(double) );
  }
  
  for(i=0; i < numrotors; i++) {
    for(j=0; j < 3*ATOMS; j++) {
      DmatrotT[i][j]=0.;
    }
  }


  /*
  DmatrotT=(double **) calloc( numrotors , sizeof(double *) );
  for(j=0; j < numrotors; j++) {
    DmatrotT[j]= (double*) calloc( 3*ATOMS, sizeof(double) );
  }
  */

  DmatrotT=transpose_double(Dmatrot,numrotors,3*ATOMS,DmatrotT);

  double **Pmatrot;

  Pmatrot=prod_mat2(Dmatrot,3*ATOMS,numrotors,DmatrotT,numrotors,3*ATOMS);
  
  for(i=0; i < 3*ATOMS; i++) {
    for(j=0; j < 3*ATOMS; j++) {
      Pmatrot[i][j]=Identity[i][j]-Pmatrot[i][j];
    }
  }

  
  double**Int2;

  Int2=prod_mat2(FMI,3*ATOMS,3*ATOMS,Pmatrot,3*ATOMS,3*ATOMS);

  double**FMIpro;

  FMIpro=prod_mat2(Pmatrot,3*ATOMS,3*ATOMS,Int2,3*ATOMS,3*ATOMS);
  
  dim=3*ATOMS+1;

  double *lambdarot;
  double **Lrot;
  double **FC_temprot;
  
  lambdarot = (double *) malloc(sizeof(double) * dim);
  Lrot = (double **) malloc((int)dim*sizeof(double *) );
  FC_temprot = (double **) malloc((int)dim*sizeof(double *) );
 
  for(j=0; j<(int)(dim); j++) {
    FC_temprot[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
    Lrot[j]=(double *) malloc((int)(dim)*sizeof(double)) ;
  }

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      if(j==i) Lrot[i][j]=1;
      else Lrot[i][j]=0.;
      FC_temprot[i][j]=0.;
    }
  }
  
  for (i=1; i<dim; i++) {
    for (j=1; j<dim; j++) {
      FC_temprot[i][j]= FMIpro[i-1][j-1];
    }
  }
  /*
  for(i=1; i<(int)(dim); i++) {    
    for(j=1; j<(int)(dim); j++) {    
      printf("%lf ",FC_temprot[i][j]);
    }
      printf("\n");
  }
  */

  if(skiprot==0){  
    jacobi(FC_temprot, dim-1, lambdarot, Lrot, nrot);
    eigsrt(lambdarot, Lrot, dim-1);
  //  printf("6D + Nrot  Projected FC diagonalization\n Eigenvalues --- \n");
  //  for(j=1; j<(int)(dim); j++) {    
  //    printf("%6.2f \n",sqrt(fabs(lambdarot[j]))/(2 * pi * c_light_cm_s ));
  //  }

    FILE *hrproj_freq,*hrproj_cd;
    //FILE *hrproj_freq;
    //  char filename[20];

    if((hrproj_freq=fopen("hrproj_freq.dat","w"))==NULL) {
      printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","hrproj_freq.dat");
      exit(1);
    }
    if((hrproj_cd=fopen("hrproj_cd.dat","w"))==NULL) {
      printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","hrproj_cd.dat");
      exit(1);
    }

    for(j=1;j<3*ATOMS+1;j++){
      if(lambdarot[j]>0.){
	fprintf(hrproj_freq," %6.2f  \n",sqrt(fabs(lambdarot[j]))/(2 * pi * c_light_cm_s ));
      } 
      else {	
	fprintf(hrproj_freq," %6.2f  \n",-sqrt(fabs(lambdarot[j]))/(2 * pi * c_light_cm_s ));
	//fprintf(hrproj_freq," %6.2f  \n",-1.0);
      } 
    }
 
    //save displacements and evaluate reduced mass
    save_displ(lambdarot,Lrot,1,3*ATOMS,step,hrproj_cd);

    fclose(hrproj_freq);
    fclose(hrproj_cd);
  }


  if(intcoord==1){

    
    FILE *hrprojint_cd;

    if((hrprojint_cd=fopen("hrprojint_cd.dat","w"))==NULL) {
      printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","hrprojint_cd.dat");
      exit(1);
    }
    
    // convert back to atomic units

    conv=joule*bohr*bohr/a_mass*1.0e10*1.0e10;
  
    //    printf("Hessian in proj funct 1,1 %lf\n",Hessian[step][0][0]);

    for (j=0;j<3*ATOMS;j++){
      for (k=0;k<3*ATOMS;k++){
	FMIpro[j][k]=FMIpro[j][k]/conv; 
      }
    }

    intcoord_hessian(step,FMIpro,gradient[step]);
    /*
    fprintf(cartdispl_rot_int," %6.2f  \n",Freqint[35]);

    for(j=0;j<3*ATOMS;j++){
      fprintf(cartdispl_rot_int," %6.2f  \n",Lcartint[j][35]);
    } 
    */
    //    fclose(cartdispl_rot_int);

    //save displacements and evaluate reduced mass
    save_displ(lambdaint,Lcartint,0,3*ATOMS-6,step,hrprojint_cd);

    for(i=0;i<(dim-6);i++){ free (Lcartint[i]);}free(Lcartint);
    free(Freqint);
    free(lambdaint);
    //  exit(0);
  } 


  //  todo free all vectors and matrixes

  for(i=0;i<dim;i++){free(Lrot[i]);}free(Lrot);
  for(i=0;i<dim;i++){free(FC_temprot[i]);}free(FC_temprot);
  free(lambdarot);
  for(i=0;i<3*ATOMS;i++){free(Int2[i]);}free(Int2);
  for(i=0;i<3*ATOMS;i++){free(FMIpro[i]);}free(FMIpro);
  for(i=0;i<3*ATOMS;i++){free(Pmatrot[i]);}free(Pmatrot);
  for(i=0;i<numrotors;i++){free(DmatrotT[i]);}free(DmatrotT);
  for(i=0;i<3*ATOMS-external;i++){free(Dprojmat[i]);}free(Dprojmat);
  for(i=0;i<3*ATOMS;i++){free(FMI[i]);}free(FMI);
  for(i=0;i<3*ATOMS-external;i++){free(Int1[i]);}free(Int1);
  for(i=0;i<3*ATOMS-external;i++){free(eigM[i]);}free(eigM);

  free(lambda_rev);
  for(i=0;i<3*ATOMS;i++){free(Vmwmat[i]);}free(Vmwmat);
  for(i=0;i<numrotors;i++){free(igroupA[i]);}free(igroupA);
  for(i=0;i<numrotors;i++){free(igroupB[i]);}free(igroupB);
  for(i=0;i<3*ATOMS;i++){free(Dmatrot[i]);}free(Dmatrot);
  free(lambda);
  for(i=0;i<3*ATOMS-external+1;i++){free(L[i]);}free(L);
  for(i=0;i<3*ATOMS-external+1;i++){free(Lt[i]);}free(Lt);
  for(i=0;i<3*ATOMS-external+1;i++){free(FC_temp[i]);}free(FC_temp);

  for(i=0;i<3*ATOMS-external;i++){free(Fc_int[i]);}free(Fc_int);
  for(i=0;i<3*ATOMS-external;i++){free(TmatT[i]);}free(TmatT);
  for(i=0;i<3*ATOMS;i++){free(Tmat[i]);}free(Tmat);
  for(i=0;i<3*ATOMS;i++){free(TLmat[i]);}free(TLmat);
  for(i=0;i<3*ATOMS;i++){free(Identity[i]);}free(Identity);
  for(i=0;i<3*ATOMS;i++){free(Pmat[i]);}free(Pmat);
  for(i=0;i<3*ATOMS;i++){free(Dmat[i]);}free(Dmat);
  for(i=0;i<ATOMS;i++){free(coords[i]);}free(coords);
  for(i=0;i<3*ATOMS;i++){free(I_temp[i]);}free(I_temp);

  
  return;
}

void save_displ(double *eval,double **evect,int ist,int dim,int step,FILE *fname){

  int i,j,is,iat;
  //  int nlines,ncols;
  int nlines;
  int jw;
  double la,lb,lc;
  double ren;
  double *redmass;

  i=1;
  nlines=dim/3;
  //  ncols=dim%3;

  redmass = (double*)malloc((int)(dim+1)*sizeof(double));

  //mass unweight and renormalize. Compute reduced masses
  ren=0.;

  for(j=ist;j<dim+ist;j++){
    ren=0.;
    for(i=ist;i<3*ATOMS+ist;i++){
      iat=(i-ist)/3;
      evect[i][j]=evect[i][j]/sqrt(atoms_data[iat].atomic_mass);
      redmass[j]+=evect[i][j]*evect[i][j];
      ren+=evect[i][j]*evect[i][j];
    } 
    ren=sqrt(ren);
    redmass[j]=1./redmass[j];
    for(i=ist;i<3*ATOMS+ist;i++){
      evect[i][j]=evect[i][j]/ren;
    } 
  } 
  // assigned reduced mass of last eigenvector to TS

  redmass_TS=redmass[dim+ist-1];

  /*
  for(j=ist;j<dim+ist;j++){
    for(i=ist;i<3*ATOMS+ist;i++){
      redmass[j]+=evect[i][j]*evect[i][j];
    } 
    redmass[j]=1./redmass[j];
  } 
  */
  //  fprintf(fname," %6.2f  \n",evect[1][1]);

  //  fprintf(fname," %d %d  \n",nlines,ncols);
  //  fprintf(fname," %d %d  \n",nlines,ncols);
  //  fprintf(fname," %d  \n",dim);
  fprintf(fname," Entering Gaussian System, \n");
  fprintf(fname,"                         Standard orientation:\n");
  fprintf(fname," ---------------------------------------------------------------------\n");
  fprintf(fname," Center     Atomic      Atomic             Coordinates (Angstroms)\n");
  fprintf(fname," Number     Number       Type             X           Y           Z\n");
  fprintf(fname," ---------------------------------------------------------------------\n");
  for(i=1;i<ATOMS+1;i++){
    fprintf(fname,"     %2d         %2d           0      %10.6f  %10.6f  %10.6f \n",i,atoms_data[i-1].at_numb,incoox[step][i-1],incooy[step][i-1],incooz[step][i-1]);
  } 
  fprintf(fname," ---------------------------------------------------------------------\n");

  fprintf(fname,"\n");

  fprintf(fname,"    99 basis functions,   99 primitive gaussians,   9 cartesian basis functions\n");
  fprintf(fname,"    99 alpha electrons       99 beta electrons\n");
  fprintf(fname,"\n");
  fprintf(fname," **********************************************************************\n");
  fprintf(fname,"\n");
  fprintf(fname," Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering\n");
  fprintf(fname," activities (A**4/AMU), depolarization ratios for plane and unpolarized\n");
  fprintf(fname," incident light, reduced masses (AMU), force constants (mDyne/A),\n");
  fprintf(fname," and normal coordinates: \n");
  //  fprintf(fname," \n");

  la=0.;
  lb=0.;
  lc=0.;
  int jcorr;
  jcorr=0;
  if(ist==0){jcorr=1;};
  
  for(j=1;j<nlines+1;j++){
    is=ist;
    jw=dim+ist-1-(j-1)*3;
    la=sqrt(fabs(eval[jw+jcorr]))/(2 * pi * c_light_cm_s );
    lb=sqrt(fabs(eval[jw-1+jcorr]))/(2 * pi * c_light_cm_s );
    lc=sqrt(fabs(eval[jw-2+jcorr]))/(2 * pi * c_light_cm_s );
   
    if(eval[jw+jcorr]<0.){
      la=-la;
    } 
    if(eval[jw-1+jcorr]<0.){
      lb=-lb;
    } 
    if(eval[jw-2+jcorr]<0.){
      lc=-lc;
    } 
    fprintf(fname,"                    %3d                    %3d                    %3d \n",1+3*(j-1),2+3*(j-1),3+3*(j-1));
    fprintf(fname,"                      A                      A                      A \n");
    fprintf(fname," Frequencies --   %9.4f              %9.4f              %9.4f \n",la,lb,lc);
    fprintf(fname," Red. masses --   %9.4f              %9.4f              %9.4f \n",redmass[jw],redmass[jw-1],redmass[jw-2]);
    fprintf(fname," Frc consts  --      0.1000                 0.1000                 0.1000\n");
    fprintf(fname," IR Inten    --      1.0000                 1.0000                 1.0000\n");
    fprintf(fname,"  Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z \n");
    //    fprintf(fname,"dim  %d  \n",dim);
    //    for(i=1;i<ATOMS+1;i++){      
    for(i=1;i<ATOMS+1;i++){
    //    for(i=1;i<ATOMS-6;i++){
      fprintf(fname,"    %2d  %2d   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f   %6.2f %6.2f %6.2f \n",i,atoms_data[i-1].at_numb,evect[is][jw],evect[is+1][jw],evect[is+2][jw],evect[is][jw-1],evect[is+1][jw-1],evect[is+2][jw-1],evect[is][jw-2],evect[is+1][jw-2],evect[is+2][jw-2]);
      is=is+3;
    } 
  } 

  return;
}


double **projector_matrix(double **P, int step){

  //Projector matrix implemented as described by Miller,Handy, and Adams J. Chem. Phys. vol.72 p90 (1980) 

  int i,j;

  /*
  int i,j;
  double **I_1,**I_temp;
 
  I_1=(double **) malloc( 4*sizeof(double *) );
  I_temp=(double **) malloc( 4*sizeof(double *) );
  for(j=0; j < 4; j++) {
    I_1[j]= (double*) malloc( 4*sizeof(double) );
    I_temp[j]= (double*) malloc( 4*sizeof(double) );
  }

  for(j=0; j < 4; j++) {
    for(i=0; i < 4; i++) {
      I_temp[i][j]=I[step][i][j];
      I_1[i][j]=0.;
    }
  }

  for(i=0;i<4;i++){
      printf("i is : %d  \n",i);
    for(j=0;j<4;j++){
      //      printf("i %d %d %lf   ",i,j,I_1[i][j]);
      printf("i %d %d %lf   ",i,j,I_temp[i][j]);
    }
  }
  */


  double **Identity;
 
  Identity=(double **) malloc( (int)(3*ATOMS)*sizeof(double *) );
  for(j=0; j < 3*ATOMS; j++) {
    Identity[j]= (double*) malloc( (int)(3*ATOMS)*sizeof(double) );
  }
  
  //Calcolo la matrice Io^-1/2

  // eliminated invert_inertia routine - not necessary anymore
  // check check check !!!


  double **I_1,**I_temp;

  I_1=(double **) malloc( 4*sizeof(double *) );
  I_temp=(double **) malloc( 4*sizeof(double *) );
  for(j=0; j < 4; j++) {
    I_1[j]= (double*) malloc( 4*sizeof(double) );
    I_temp[j]= (double*) malloc( 4*sizeof(double) );
  }

  for(j=0; j < 4; j++) {
    for(i=0; i < 4; i++) {
      I_temp[i][j]=I[step][i][j];
      I_1[i][j]=0.;
    }
  }

  //      I_1 = invert_inertia_matrix( I[step],4,I_1);
  //I_1 = invert_inertia_matrix( I_temp,4,I_1);
       	 I_1 = inverse_matrix(I_temp,3,I_1);

  for(i=0;i<4;i++){free(I_temp[i]);}free(I_temp);

  // compute square root of I_1 (diagonal matrix)
  
  double ** v_temp;
  double ** d_temp;
  double ** v_tempt;
  double * eigen_temp;

  v_temp = (double**)malloc((int)(4) * sizeof(double*));
  for(j=0;j<4;j++) {
    v_temp[j] = (double*)malloc((int)(4)*sizeof(double));
  }
  d_temp = (double**)malloc((int)(4) * sizeof(double*));
  for(j=0;j<4;j++) {
    d_temp[j] = (double*)malloc((int)(4)*sizeof(double));
  }
  v_tempt = (double**)malloc((int)(4) * sizeof(double*));
  for(j=0;j<4;j++) {
    v_tempt[j] = (double*)malloc((int)(4)*sizeof(double));
  }
  eigen_temp = (double*)malloc((int)(4)*sizeof(double));
 
  //calculate eigenvectors and eigenvalues of I_1
  jacobi(I_1, 3, eigen_temp, v_temp, nrot);   

  //calculate transpose of eigenvalue matrix

  for(i=1;i<4;i++){
    for(j=1;j<4;j++){
      v_tempt[i][j]=v_temp[j][i];
    }
  }

  //initialize diagonal matrix with square root eigenvalues 

  for(i=1;i<4;i++){
    for(j=1;j<4;j++){
      d_temp[i][j]=0;
      if(i == j)  d_temp[i][j]=sqrt(eigen_temp[i]);
    }
  }


  //calculate sqrt of I_1 matrix

  I_1=prod_mat(d_temp,v_tempt,4);
  I_1=prod_mat(v_temp,I_1,4);

  //free temporary vectors

  for(i=0;i<4;i++){ free (v_temp[i]); } free(v_temp);
  for(i=0;i<4;i++){ free (d_temp[i]); } free(d_temp);
  for(i=0;i<4;i++){ free (v_tempt[i]); } free(v_tempt);
  free(eigen_temp);

  //Calcolo gli Liy legati alle traslazioni

  double Total_mass=0.0;
  for(i=0;i<ATOMS;i++){
    Total_mass =  Total_mass + atoms_data[i].atomic_mass;
  }
  for(i=0;i<ATOMS;i++){
    // lungo x:
    atoms_data[i].LTr1x=sqrt( atoms_data[i].atomic_mass/Total_mass);
    atoms_data[i].LTr1y=0.0;
    atoms_data[i].LTr1z=0.0;
    // lungo y:
    atoms_data[i].LTr2x=0.0;
    atoms_data[i].LTr2y=sqrt( atoms_data[i].atomic_mass/Total_mass);
    atoms_data[i].LTr2z=0.0;
    // lungo z:
    atoms_data[i].LTr3x=0.0;
    atoms_data[i].LTr3y=0.0;
    atoms_data[i].LTr3z=sqrt( atoms_data[i].atomic_mass/Total_mass);
  }

  // here coordinates are re-initialized as mass weighted
 
  mass_weight_coordinates(step); 
  
  //calcolo gli Liy legati alle rotazioni
  
  for(i=0;i<ATOMS;i++){
    //cc double checked, it is ok
    
    //intorno a x
    atoms_data[i].LRot1x = (atoms_data[i].ext_coord_z[step]*I_1[1][2]- atoms_data[i].ext_coord_y[step]*I_1[1][3]);
    atoms_data[i].LRot1y = (atoms_data[i].ext_coord_x[step]*I_1[1][3]- atoms_data[i].ext_coord_z[step]*I_1[1][1]);
    atoms_data[i].LRot1z = (atoms_data[i].ext_coord_y[step]*I_1[1][1]- atoms_data[i].ext_coord_x[step]*I_1[1][2]);
    //intorno a y
    atoms_data[i].LRot2x = (atoms_data[i].ext_coord_z[step]*I_1[2][2]- atoms_data[i].ext_coord_y[step]*I_1[2][3]);
    atoms_data[i].LRot2y = (atoms_data[i].ext_coord_x[step]*I_1[2][3]- atoms_data[i].ext_coord_z[step]*I_1[2][1]);
    atoms_data[i].LRot2z = (atoms_data[i].ext_coord_y[step]*I_1[2][1]- atoms_data[i].ext_coord_x[step]*I_1[2][2]);
    //intorno a z
    atoms_data[i].LRot3x = (atoms_data[i].ext_coord_z[step]*I_1[3][2]- atoms_data[i].ext_coord_y[step]*I_1[3][3]);
    atoms_data[i].LRot3y = (atoms_data[i].ext_coord_x[step]*I_1[3][3]- atoms_data[i].ext_coord_z[step]*I_1[3][1]);
    atoms_data[i].LRot3z = (atoms_data[i].ext_coord_y[step]*I_1[3][1]- atoms_data[i].ext_coord_x[step]*I_1[3][2]);

  }
  
  //genero i 6 vettori Liy lunghi 3*N 

  double *Liy_rot1,*Liy_rot2,*Liy_rot3,*Liy_tr1,*Liy_tr2,*Liy_tr3;
 
  Liy_rot1=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );
  Liy_rot2=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );
  Liy_rot3=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );
  Liy_tr1=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );
  Liy_tr2=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );
  Liy_tr3=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );


  j=0;
  for(i=0;i<3*ATOMS;i=i+3){
    Liy_rot1[i]=atoms_data[j].LRot1x;
    Liy_rot1[i+1]=atoms_data[j].LRot1y;
    Liy_rot1[i+2]=atoms_data[j].LRot1z;
    
    Liy_rot2[i]=atoms_data[j].LRot2x;
    Liy_rot2[i+1]=atoms_data[j].LRot2y;
    Liy_rot2[i+2]=atoms_data[j].LRot2z;

    Liy_rot3[i]=atoms_data[j].LRot3x;
    Liy_rot3[i+1]=atoms_data[j].LRot3y;
    Liy_rot3[i+2]=atoms_data[j].LRot3z;

    Liy_tr1[i]=atoms_data[j].LTr1x;
    Liy_tr1[i+1]=atoms_data[j].LTr1y;
    Liy_tr1[i+2]=atoms_data[j].LTr1z;

    Liy_tr2[i]=atoms_data[j].LTr2x;
    Liy_tr2[i+1]=atoms_data[j].LTr2y;
    Liy_tr2[i+2]=atoms_data[j].LTr2z;
     
    Liy_tr3[i]=atoms_data[j].LTr3x;
    Liy_tr3[i+1]=atoms_data[j].LTr3y;
    Liy_tr3[i+2]=atoms_data[j].LTr3z;
     
    j=j+1;

  }
  
  // printf("\n\n6 Liy vectors that corresponds to infinitesimal translations and rotations of the N-Atoms system\n");
  // printf ("\nLiy_tr1\t Liy_tr2\t Liy_tr3\t Liy_rot1\t Liy_rot2\t Liy_rot3\n");
 
  // genero i vettori Liy  per  dof=reaction coordinate

  double *Liy_rc, tot=0;
  Liy_rc=(double *) malloc( (int)(3*ATOMS)*sizeof(double ) );

  j=0;

  for(i=0;i<3*ATOMS;i=i+3){
    //gradient in hartree/bohr is converted to mass weighted coordinates
    Liy_rc[i]=atoms_data[j].grad_x[step]/sqrt(atoms_data[j].atomic_mass);
    Liy_rc[i+1]=atoms_data[j].grad_y[step]/sqrt(atoms_data[j].atomic_mass);
    Liy_rc[i+2]=atoms_data[j].grad_z[step]/sqrt(atoms_data[j].atomic_mass);

    j=j+1;    
  }


  for(i=0;i<3*ATOMS;i++){
    tot = tot +  Liy_rc[i]* Liy_rc[i]; 
  }

  tot=sqrt(tot);
 
  if(tot!=0){
    for(i=0;i<3*ATOMS;i++){
      Liy_rc[i]=-Liy_rc[i]/tot; 
      //    Liy_rc[i]=0.;
      printf("%lf\t %lf\n", Liy_rc[i],gradient[step][i]);    
    }
  }

  // genero la matrice P e la matrice I
  
  for(i=0;i<3*ATOMS;i++){
    for(j=0;j<3*ATOMS;j++){

      P[i][j]=Liy_tr1[i]*Liy_tr1[j]+Liy_tr2[i]*Liy_tr2[j]+Liy_tr3[i]*Liy_tr3[j]+ Liy_rot1[i]*Liy_rot1[j]+ Liy_rot2[i]*Liy_rot2[j]+ Liy_rot3[i]*Liy_rot3[j]+ Liy_rc[i]* Liy_rc[j];
      //cc test: do not project reaction coordinate
      if(proj_rea_coo==1){
	P[i][j]=Liy_tr1[i]*Liy_tr1[j]+Liy_tr2[i]*Liy_tr2[j]+Liy_tr3[i]*Liy_tr3[j]+ Liy_rot1[i]*Liy_rot1[j]+ Liy_rot2[i]*Liy_rot2[j]+ Liy_rot3[i]*Liy_rot3[j];
      }
      if(i==j) Identity[i][j]=1.0;
      else Identity[i][j]=0.0;
 
    }
  }

  free(Liy_rot1); free(Liy_rot2); free(Liy_rot3); free(Liy_tr1);free(Liy_tr2);free(Liy_tr3);free(Liy_rc);
  for(i=0;i<4;i++){free(I_1[i]);}free(I_1);

  
  //  printf("ATOM is i %d \n",ATOMS);
  /*
  for(i=0;i<3*ATOMS;i++){
    for(j=0;j<3*ATOMS;j++){
      printf("proj value at i %d j %d is %lf ",i,j,P[i][j]);
    }
      printf("\n");
  }
    exit(-1);
  */
  
  for(i=0;i<3*ATOMS;i++){
    for(j=0;j<3*ATOMS;j++){
      P[i][j]=Identity[i][j]-P[i][j];
    }
  }

  for(i=0;i<3*ATOMS;i++){free(Identity[i]);}free(Identity);

  return P;
}


void mass_weight_coordinates(int step){
  int i;

  
  printf("Coordinates mass weighted\n");
  printf("Atom\t     x\t     y\t     z\n");

  for(i=0;i<ATOMS;i++){
    
    atoms_data[i].ext_coord_x[step]= atoms_data[i].ext_coord_x[step]*sqrt(atoms_data[i].atomic_mass);    
    atoms_data[i].ext_coord_y[step]= atoms_data[i].ext_coord_y[step]*sqrt(atoms_data[i].atomic_mass);   
    atoms_data[i].ext_coord_z[step]= atoms_data[i].ext_coord_z[step]*sqrt(atoms_data[i].atomic_mass);
   
    printf("%d\t%lf\t%lf\t%lf\n",i, atoms_data[i].ext_coord_x[step], atoms_data[i].ext_coord_y[step], atoms_data[i].ext_coord_z[step]);
  }
}


double **prod_mat(double **A, double **B, int dim){
  

  int i, j,k;
  double **Prod;

  Prod=(double **) malloc( dim*sizeof(double *) );
  for(j=0; j < dim; j++) {
    Prod[j]= (double*) malloc( dim*sizeof(double) );
  }

  for(i=0;i<dim;i++){

    for(j=0;j<dim;j++){

      Prod[i][j]=0.;

      for(k=0;k<dim;k++){

	Prod[i][j]= Prod[i][j]+  A[i][k]*B[k][j] ;
	//printf("Prod[%d][%d]=%lf\n ",i,j,Prod[i][j]);
	
      }
    }
  }

  return Prod;

}

double **prod_mat2(double **A, int dim1, int dim2, double **B, int dim3, int dim4){
  

  int i,j,k;
  double **Prod;

  Prod=(double **) malloc( dim1*sizeof(double *) );
  for(j=0; j < dim1; j++) {
    Prod[j]= (double*) malloc( dim4*sizeof(double) );
  }

  if(dim2!=dim3){
    printf("Error columns and rows of matixes A and B must be equal\n");
    exit(-1);
  }
    
  for(i=0;i<dim1;i++){

    for(j=0;j<dim4;j++){

      Prod[i][j]=0.;

      for(k=0;k<dim2;k++){

	Prod[i][j]= Prod[i][j]+  A[i][k]*B[k][j] ;
	//printf("Prod[%d][%d]=%lf\n ",i,j,Prod[i][j]);
	
      }
    }
  }

  return Prod;

}


double **invert_matrix(double **A, int dim){
  
  int i,j;
  double **AM1;
  double **Atemp;
  double **AM1temp;


  AM1=(double **) malloc( dim*sizeof(double *) );
  for(j=0; j < dim; j++) {
    AM1[j]= (double*) malloc( dim*sizeof(double) );
  }

  Atemp=(double **) malloc( (dim+1)*sizeof(double *) );
  for(j=0; j < dim+1; j++) {
    Atemp[j]= (double*) malloc( (dim+1)*sizeof(double) );
  }
  AM1temp=(double **) malloc( (dim+1)*sizeof(double *) );
  for(j=0; j < dim+1; j++) {
    AM1temp[j]= (double*) malloc( (dim+1)*sizeof(double) );
  }

  for(i=0;i<dim+1;i++){
    for(j=0;j<dim+1;j++){
      Atemp[i][j]=0.;
      AM1temp[i][j]=0.;
    }
  }
    
  for(i=0;i<dim;i++){
    for(j=0;j<dim;j++){
      Atemp[i+1][j+1]=A[i][j];
    }
  }

  AM1temp=inverse_matrix(Atemp,dim,AM1temp);

  for(i=1;i<dim+1;i++){
    for(j=1;j<dim+1;j++){
      AM1[i-1][j-1]=AM1temp[i][j];
    }
  }

  for(i=0;i<(dim);i++){ free (Atemp[i]); } free(Atemp);
  for(i=0;i<(dim);i++){ free (AM1temp[i]); } free(AM1temp);
  
  return AM1;

}


void cross_prod(double ax, double ay, double az, double bx, double by, double bz, double *px, double *py, double *pz)
{

  *px=ay*bz-az*by;
  *py=az*bx-ax*bz;
  *pz=ax*by-ay*bx;

  return;

}
 

double **transpose_double(double **A, int l_side, int h_side, double **B)
{
  int i,j;

  for(j=0; j < l_side; j++) {
    for (i=0; i<h_side; i++) {
      B[j][i]=A[i][j];      
    }
  }
 return B;
}


double dot_prod(double *v1,double *v2, int dim){

  int i;
  double prod=0;

  for(i=0;i<dim;i++){
    prod = prod + v1[i]*v2[i];
  }
  return prod;
}


double **calc_deriv_matrix(int step,double **Matrix_back,double **Matrix_for, double **Deriv_Matrix ,int dim){

  int i,j;
  dim=dim+1;

  Deriv_Matrix=(double **) malloc( dim*sizeof(double *) );
  for(j=0; j < dim; j++) {
    Deriv_Matrix[j]= (double*) malloc( dim*sizeof(double) );
  }

  for(j=1;j<(dim);j++){
    for(i=1;i<dim;i++){
      printf("%lf\n",Matrix_for[i][j]);
    }
  }
  
  for(j=0;j<(dim);j++){
    for(i=0;i<dim;i++){
      Deriv_Matrix[i][j]=(Matrix_for[i][j]-Matrix_back[i][j])/(2*deltas);	   
    }
  }
  
  return Deriv_Matrix;

}



double *calc_BkF(int dim,int step, double **L_int, double *BkF){

  //calcolo BkF secondo Page, McIver J. Chem. Phys. vol 88, pp 922-935 (1988) 
  // Bkf=Li_t*H*v/c

  int i,k,j;
    
  for(k=0;k<dim;k++){
    BkF[k]=0;    
  }

  double **temp, cost=0;
  double *bkfcorr;
  //Hessian[step]=force_constants_mass_unweight(Hessian[step]);

  temp=(double **) malloc( (int)(dim)*sizeof(double *) );
  bkfcorr=(double *) malloc( (int)(dim)*sizeof(double) );
  for(i=0; i < (int)(dim); i++) {
    temp[i]= (double*) malloc( (int)(dim)*sizeof(double) );
  }
        
  for(j=0;j<dim;j++){
    bkfcorr[j]=0.;
    for(k=0;k<dim;k++){
      temp[j][k]=0.;
    }
  }

  // product Li_t*H  N.B. Li transposed

  for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){	
      for(i=0;i<dim;i++){
	//cc	temp[j][k]=temp[j][k] + L_int[i][k]*Hessian[step][i][j];
	temp[k][j]=temp[k][j] + L_int[i][k]*Hessian[step][i][j];
      }
    }
  }

  // Normalization constant C

  for(j=0;j<dim;j++){    
    cost=cost+gradient[step][j]*gradient[step][j];
  }
  cost=sqrt(cost);     
  if((step == saddlep) && (cost==0.)){
    cost=1.;        
  }  

  // product gradt*H*grad. its a number

  double ffc=0;
  double *ffv_tempc;
  ffv_tempc=(double *) malloc( (int)(dim)*sizeof(double) );

  for(j=0;j<dim;j++){
    ffv_tempc[j]=0;
  }

  for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      ffv_tempc[j]=ffv_tempc[j]+Hessian[step][j][i]*(gradient[step][i]/cost);
    }
  }

  ffc=0.;
  for(j=0;j<dim;j++){
    ffc=ffc+ffv_tempc[j]*gradient[step][j]/cost;
  }

  free(ffv_tempc);

  // product Li_t*grad*ffv  N.B. Li transposed


  for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
	bkfcorr[k]=bkfcorr[k]+ffc*L_int[j][k]*(gradient[step][j])/cost;
    }
  }



  for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
      //      if(step<=MAXSTEP/2)
	//cc	BkF_temp[k]=BkF_temp[k]+temp[j][k]*(gradient[step][j])/cost;
      BkF[k]=BkF[k]+temp[k][j]*(gradient[step][j])/cost;
	//else
	//cc	BkF_temp[k]=BkF_temp[k]+temp[j][k]*(gradient[step][j])/cost;
	//cc NB changed signed as steepest descent, but no effect 
	//	BkF_temp[k]=BkF_temp[k]-temp[k][j]*(gradient[step][j])/cost;
    }
    //    BkF[k]=BkF[k]-bkfcorr[k];
  }

  //cc check sign + or - ??? !!!!!

  for(k=0;k<dim;k++){
    BkF[k]=-BkF[k]/cost;
    //  BkF[k]=BkF[k]/cost;
    //    if(step == 100 ){
    //  printf("bkf %le\n",BkF[k]);        
    //}  
  }
  //  if(step == 100 ){
  //    exit(0);
  //  }

  for(i=0;i<(dim);i++){ free (temp[i]); } free(temp);
  free(bkfcorr);


  // calculation of BkF in the TS region with eq 49 of M. Page and J.W. McIver J.Chem.Phys. vol 88,p922 (1988)
  // npoints on the left and right of the saddle point will be evaluated
  // however at present only the saddle point is used in the calculations
  
  //  int npoints_rem=2;
  int npoints_rem;
  saddlepbkf=1;
  // option disabled as long as it has been debugged
  saddlepbkf=0;
  if(saddlepbkf==0) npoints_rem=0;
  if(saddlepbkf!=0) npoints_rem=saddlepbkf;
  
  if(step>((int)(MAXSTEP2TS/2-npoints_rem))&& step<((int)(MAXSTEP2TS/2+npoints_rem))){

    //calcolo fv=vt*F*v è uno scalare!

    double ffv=0;
    double *ffv_tempv;
    ffv_tempv=(double *) malloc( (int)(dim)*sizeof(double) );

    for(j=0;j<dim;j++){
      ffv_tempv[j]=0;
    }

    for(j=0;j<dim;j++){
      for(i=0;i<dim;i++){
	ffv_tempv[j]=ffv_tempv[j]+Hessian[step][j][i]*(gradient[step][i]/cost);
      }
    }

    ffv=0.;
    for(j=0;j<dim;j++){
      ffv=ffv+ffv_tempv[j]*gradient[step][j]/cost;
    }

    free(ffv_tempv);

    //    printf("step %d ffv=%lf\n",step,ffv);
    

    //calcolo Delta=(2*vt*F*v*I-F)

    double **Delta;
    Delta=(double **) malloc( (int)(dim+1)*sizeof(double *) );
    for(i=0; i < (int)(dim+1); i++) {
      Delta[i]= (double*) malloc( (int)(dim+1)*sizeof(double) );
    }

    for(j=0;j<dim+1;j++){
      for(i=0;i<dim+1;i++){
	Delta[i][j]=0.;
      }
    }


    for(i=1;i<=dim;i++){
      for(j=1;j<=dim;j++){
     	if(i==j) Delta[i][j]=2*ffv-Hessian[step][i-1][j-1];
	//	else Delta[i][j]=Hessian[step][i-1][j-1];		   
	else Delta[i][j]=-Hessian[step][i-1][j-1];		   
      }
    }

    //calcolo la matrice inversa di Delta --> Delta_inv

    double **Delta_inv;
    Delta_inv=(double **) malloc( (int)(dim+1)*sizeof(double *) );
    for(i=0; i < (int)(dim+1); i++) {
      Delta_inv[i]= (double*) malloc( (int)(dim+1)*sizeof(double) );
    }

    for(j=0;j<dim;j++){
      for(i=0;i<dim;i++){
	Delta_inv[i][j]=0.;
      }
    }

    Delta_inv=inverse_matrix(Delta,dim,Delta_inv);
    
    for(i=0;i<(dim);i++){ free (Delta[i]); } free(Delta);

    //Calcolo il prodotto Lt*Delta_inv

    double **LDelta_inv;
    LDelta_inv=(double **) malloc( (int)(dim)*sizeof(double*) );
    for(i=0; i < (int)(dim); i++) {
      LDelta_inv[i]= (double*) malloc( (int)(dim)*sizeof(double) );
    }

    for(k=0;k<dim;k++){
      for(j=0;j<dim;j++){
	LDelta_inv[k][j]=0.;
      }
    }

    //NB L_int transposed
    for(k=0;k<dim;k++){
      for(j=0;j<dim;j++){
	for(i=0;i<dim;i++){
	  //cc_err	  LDelta_inv[j][k]= LDelta_inv[j][k]+ L_int[i][k]*Delta_inv[i+1][j+1];
	  LDelta_inv[k][j]= LDelta_inv[k][j]+ L_int[i][k]*Delta_inv[i+1][j+1];
	}
      }
    }

    for(i=0;i<(dim+1);i++){ free (Delta_inv[i]); } free(Delta_inv);
      
    //calcolo dffv=vt*dF/ds*v è uno scalare!

    double dffv;
    double *dffv_tempv;
    dffv_tempv=(double *) malloc( (int)(dim)*sizeof(double) );

    for(j=0;j<dim;j++){
      dffv_tempv[j]=0;
    }

    dffv=0;

    for(j=0;j<dim;j++){
      for(i=0;i<dim;i++){
	dffv_tempv[j]=dffv_tempv[j]+Hessian_deriv[step][j][i]*(gradient[step][i]/cost);	
      }
    }

    for(i=0;i<dim;i++){	
      dffv=dffv+gradient[step][i]/cost*dffv_tempv[i];
    }

    free(dffv_tempv);

    printf("step %d dffv_n= %lf\n",step,dffv);


   //calcolo delta_piccolo=(dF/ds-vt*dF/ds*v*I)

    double **delta_piccolo;
    delta_piccolo=(double **) malloc( (int)(dim)*sizeof(double *) );
    for(i=0; i < (int)(dim); i++) {
      delta_piccolo[i]= (double*) malloc( (int)(dim)*sizeof(double) );
    }

    for(j=0;j<dim;j++){
      for(i=0;i<dim;i++){
	delta_piccolo[i][j]=0.;
      }
    }


    for(i=0;i<dim;i++){
      for(j=0;j<dim;j++){
     	if(i==j) delta_piccolo[i][j]=Hessian_deriv[step][i][j]-dffv;
	else delta_piccolo[i][j]=Hessian_deriv[step][i][j];
      }
    }

    //calcolo BkF nella regione del TS

    double *ld_temp;
    ld_temp=(double *) malloc( (int)(dim)*sizeof(double) );

    for(j=0;j<dim;j++){
      ld_temp[j]=0;
      BkF[j]=0;
    }

    for(k=0;k<dim;k++){
      for(j=0;j<dim;j++){
	ld_temp[k]= ld_temp[k]+delta_piccolo[k][j]*(gradient[step][j]/cost);
      }
    }


    for(k=0;k<dim;k++){
      for(i=0;i<dim;i++){
	BkF[k]=BkF[k]+LDelta_inv[k][i]*ld_temp[i];
	//if(step<=MAXSTEP/2)BkF_temp[k]=BkF_temp[k]+LDelta_inv[i][k]*ld_temp[i];
	//	else  BkF_temp[k]=BkF_temp[k]-LDelta_inv[i][k]*ld_temp[i];

      }
    }

    for(i=0;i<(dim);i++){ free (delta_piccolo[i]); } free(delta_piccolo);
    free(ld_temp);

    //    if(step<MAXSTEP/2){
    for(k=0;k<dim;k++){
      BkF[k]=-BkF[k];
    }
      //    }
  }

  //    k=0;
  //    printf("step %d bkf= %le \n",step, BkF[k]);

  return BkF;
}

double ***check_L_sign(int step, double ***L_int_save){

  double *temp,*temp2;
  int k,i;

  temp=(double *) malloc((int)(3*ATOMS-7)*sizeof(double ) );
  temp2=(double *) malloc((int)(3*ATOMS-7)*sizeof(double ) );

  for(k=0;k<(3*ATOMS-7);k++){
    
    temp[k]=0;
    temp2[k]=0;
    
  }


  for(k=0;k<(3*ATOMS-7);k++){
    for(i=0;i<3*ATOMS;i++){
      if(step>0){
      temp[k]= temp[k]+((L_int_save[step][i][k]));
      temp2[k]= temp2[k]+((L_int_save[step-1][i][k]));
    }
    }
    if((fabs(fabs(temp[k])-fabs(temp2[k]))<fabs(0.3*temp2[k]))&& (sign(temp[k])!=sign(temp2[k]))){
      printf("warning! changing L sing at step=%d for the degree of freedom n° %d\n",step,k);
      
      for(i=0;i<3*ATOMS;i++){

	L_int_save[step][i][k]=-L_int_save[step][i][k];
      }

    }

  } 

  return L_int_save;
  
}

void check_L_order(int step, double ***L_int_save){

  int k;

  for(k=0;k<(3*ATOMS-7);k++){
    if(fabs(fabs(L_int_save[step][0][k])-fabs(L_int_save[step-1][0][k]))>fabs(0.3*L_int_save[step-1][0][k])){

      printf("warning changing L order at step %d for the degree of freedom n°%d\n",step,k);     

    }
  }
}



int sign(double q){

  if(q<0)return -1;
  else if(q>0)return +1;
  else return 0;

}


int compare(const void *v1, const void *v2) {

  if(*(double *)v1>*(double *)v2) { return 1;}
  else if(*(double *)v1<*(double *)v2) {return -1;}
  else {return 0;}

}


int calc_VaG(double **frequencies_save){

  double *E_temp;
  int step;
  int i,j;

  E_temp=(double *) malloc( MAXSTEP*sizeof(double ) );
  
  double *xE, *Energy_temp, *Energy_spline, *E_to_sort, *E_int_temp, *E_int_spline, *Ezpe;
  int N_spline=0;
  xE=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
  Energy_temp=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
  Energy_spline=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
  E_to_sort=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );
  Ezpe=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );

  E_int_temp=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
  E_int_spline=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );

  FILE *VaG;
  if((VaG=fopen("VaG.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","kuni.txt");
    exit(1);
  }

  for (step=0;step<MAXSTEP;step++){    
    E_temp[step]=0;
    Ezpe[step]=0;
  }

  // add zero point energy to MEP, using projected frequencies

  for(step=0;step<MAXSTEP;step++){    
    for(i=0;i<(int)(3*ATOMS-7);i++){        
      //      E_temp[step]=E_temp[step]+0.5*(frequencies_save[step][i]); //Energy in cm-1      
      // watch out! this condition means that negative frequencies are not used to calculat the ZPE correction
      // this is consistent since it is expected that neg. freq. are fictitious.
      if(frequencies_save[step][i]>0) E_temp[step]=E_temp[step]+0.5*fabs(frequencies_save[step][i]); //Energy in cm-1      
      //cctest
      if(iminfreq!=0){
	if(frequencies_save[step][i]<iminfreq) E_temp[step]=E_temp[step]+0.5*iminfreq;
      }
      Ezpe[step]+=0.5*fabs(frequencies_save[step][i]);
      if(step==MAXSTEP2TS/2){
	//	printf("freq_check %d %lf\n",i,frequencies_save[step][i]);   
      }
    }
    E_temp[step]=E_temp[step];//*0.5*n_avogadro*h_planck*c_light_cm_s*4187*349.75;    
  }

  //Interpolo E_int//

  E_int_temp=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
  E_int_spline=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );

  for(j=0;j<(int)(MAXSTEP2TS/2-N_spline/2);j++){  
    E_int_temp[j]= E_temp[j];
  }
   
   
  for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<MAXSTEP-N_spline;j++){  
    E_int_temp[j]= E_temp[j+N_spline];
  }
      
  for(j=0;j<MAXSTEP-N_spline;j++){ 
    xE[j]=j;
    if(j>=(int)(MAXSTEP2TS/2-N_spline/2)) xE[j]=j+N_spline;
  }
     
  for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<(int)(MAXSTEP2TS/2+N_spline/2);j++){       

    E_int_spline= spline(xE,E_int_temp, MAXSTEP-N_spline-1, 1e30, 1e30,E_int_spline);
    E_temp[j]=splint(xE, E_int_temp,E_int_spline, MAXSTEP-N_spline-1,j,E_temp[j]);    
  }
     
  //  double E_n=fabs(Energy[0]);
  
  for(step=0;step<MAXSTEP;step++){    
    //updated to new format of the input
    //    Energy[step]=((Energy[step]-Ereactants)*627.5)*349.75; //Energy in cm-1
    //    printf("Energy_check[%d]=%lf\n",step,Energy[step]);   
    Energy[step]=((Energy[step])*627.5+EBarr_kcalmol)*349.75; //Energy in cm-1
    Energy[step]=(Energy[step]+E_temp[step]); //Energy in cm-1         
  }


  ///////////////INTERPOLO VaG///////////////////////
  
  for(j=0;j<(int)(MAXSTEP2TS/2-N_spline/2);j++){
    Energy_temp[j]= Energy[j];     
  }
    
  for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<MAXSTEP-N_spline;j++){ 
    Energy_temp[j]= Energy[j+N_spline];
  }
        
  for(j=0;j<MAXSTEP-N_spline;j++){   
    xE[j]=j;
    if(j>=(int)(MAXSTEP2TS/2-N_spline/2)) xE[j]=j+N_spline;
  }
   
  for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<(int)(MAXSTEP2TS/2+N_spline/2);j++){   
    Energy_spline= spline(xE,Energy_temp, MAXSTEP-N_spline-1, 1e30, 1e30,Energy_spline);
    Energy[j]=splint(xE, Energy_temp,Energy_spline, MAXSTEP-N_spline-1,j,Energy[j]);
  }
   
 ////////////////////////////////////////

  //commento se uso Energia calcolata trmaite IRC, decommento se uso Vag letta da file (es Eckart) 

  /*

  double *pes_spline;
  pes_spline=(double *) malloc( dim_newpes*sizeof(double ) ); 
  pes_spline= spline(x_pes_read, pes_read, dim_newpes-1, 1e30, 1e30,pes_spline);

  for(i=0;i<MAXSTEP;i++){
    Energy[i]= (splint(x_pes_read, pes_read, pes_spline, dim_newpes-1, i, Energy[i]))*349.75;
  }
  */

  //////////////////////////////////

  for(i=0;i<MAXSTEP;i++){
    //    Energy[i]= (Energy[i]);
    E_to_sort[i]=Energy[i];     
    printf("Energy[%d]=%lf\n",i,Energy[i]);   
  }

  for(i=0;i<MAXSTEP;i++){
    fprintf(VaG,"%d\t%lf\t%lf\t%lf\n",i,(Energy[i]-E_temp[i])/349.75,Ezpe[i]/349.75,((Energy[i]-E_temp[i]+Ezpe[i])/349.75));
  }

  fclose(VaG); 
 
  qsort(E_to_sort,MAXSTEP, sizeof (double), compare);
  qsort(E_to_sort,MAXSTEP, sizeof (double), compare);

  // int  Emax=E_to_sort[MAXSTEP-1];
  int Emax=Energy[(int)(MAXSTEP2TS/2)];

  //  int Emax=E_to_sort[MAXSTEP-1]+0.5;

  printf("Emax=%d\n",Emax); 

  free(Ezpe);
  free(E_temp);
  
  return Emax;
}


double *calc_transmission_coeff(double *mueff,int Emax, double *Transmission_coeff ){
 
  int E, E0;
  int step;
  int i,j;
  double sl,sr;
  double *theta;
  double **dtheta;
  double *xE, *Energy_spline;
  FILE *dthetaE, *imactint ;

 
  //open output file

  if((dthetaE=fopen("dthetaE.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","dthetaE.txt");
    exit(1);
  }

  if((imactint=fopen("imactint.dat","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","imactint.dat");
    exit(1);
  }

  // initialize vectors

  xE=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );
  for(j=0;j<MAXSTEP;j++){   
    xE[j]=j;
  }
  
  Energy_spline=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );
  Energy_spline= spline(xE,Energy, MAXSTEP-1, 1e30, 1e30,Energy_spline);

  for (E=0;E<E_dim;E++){    
    Transmission_coeff[E]=0;
  }
   
  theta=(double *) malloc( (int)(Emax)*sizeof(double ) );
  dtheta=(double **) malloc( (int)(Emax)*sizeof(double *) );

  for(i=0; i < (int)(Emax); i++) {
    dtheta[i]= (double*) malloc( (int)(MAXSTEP)*sizeof(double) );
  }
  
  for(i=0; i < (int)(Emax); i++) { 
    theta[i]=0.0;
    for(j=0;j<(int)(MAXSTEP);j++){
      dtheta[i][j]=0.0;
    }    
  }

  // define E0 (it is an integer): the minimum value for which the imaginary action integral is defined

  if( Energy[0]==Energy[MAXSTEP-1] || Energy[0]>Energy[MAXSTEP-1] ) E0=(Energy[0]+0.5);
  else if(Energy[MAXSTEP-1]>Energy[0]) E0=(Energy[MAXSTEP-1]+0.5);
  else {printf("Error in defining E0"); exit(-1);}

  printf("Energy[0]=%lf\tEnergy[MAXSTEP-1]=%lf\tE0=%d\n",Energy[0],Energy[MAXSTEP-1],E0);

  //calculate dtheta, the function to be integrated

  /*
  double *x, *mueff_spline,*dtheta_spline,y=0;
  FILE *mueff_spline_result, *thetaE;

  if((mueff_spline_result=fopen("mueff_spline_result.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","kuni.txt");
    exit(1);
  }
  
  if((thetaE=fopen("theta.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","kuni.txt");
    exit(1);
  }
  
  //  mueff_spline=(double *) malloc( MAXSTEP*sizeof(double ) ); 
  
  */

  
  /*
  //mueff_spline= spline(x, mueff, MAXSTEP-1, 1e30, 1e30,mueff_spline);


  mueff_spline=(double *) malloc( dim_truong_mueff*sizeof(double ) ); 
  dtheta_spline=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) ); 

  mueff_spline= spline(x_mueff_read, mueff_read, dim_truong_mueff-1, 1e30, 1e30,mueff_spline);
  
  double stepn;

  for(stepn=0;stepn<MAXSTEP;stepn=stepn+0.01){    
    y=splint(x_mueff_read, mueff_read, mueff_spline, dim_truong_mueff-1, stepn, y);
    fprintf(mueff_spline_result,"%lf\t%le\n",stepn,y);
  }
  
  y=0;
  */

  double Energy_test=0;
     
  for(E=0;E<(int)((Emax-E0));E++){
    for(step=0;step<(int)(MAXSTEP);step++){

      //for SCT calculation NB: d_theta[E][step] è adimensionale
            
      //      Energy_test=splint(xE,Energy, Energy_spline, MAXSTEP-1, (step), Energy_test);
      Energy_test=Energy[step];
      dtheta[E][step]=sqrt(2*mueff[step]*e_mass*redmu*(fabs(Energy_test-((E+E0)))*c_light_cm_s*h_planck));//*a_bohr;

      //for ZCT calculation
      if(izct==1){
	dtheta[E][step]=sqrt(2*redmu*(fabs(Energy_test-((E+E0)))*c_light_cm_s*h_planck));//*a_bohr;          
      }
    }
  }
  
  

  /*
  double *dtheta_dfirst, *dtheta_dlast;
  dtheta_dfirst=(double *) malloc( (int)(Emax)*sizeof(double ) );
  dtheta_dlast=(double *) malloc( (int)(Emax)*sizeof(double ) );


  for(E=0;E<(int)((Emax-E0));E++){

    dtheta_dfirst[E]= (dtheta[E][1]-dtheta[E][0]);
    dtheta_dlast[E]= (dtheta[E][MAXSTEP-1]-dtheta[E][MAXSTEP-2]);
  }

  for(step=0;step<MAXSTEP;step++){    
    fprintf(dthetaE,"%d\t",step);
  }

  */

  


  /*
  for(step=0;step<MAXSTEP;step++){
    
    fprintf(dthetaE,"%d\t",step);
    
    for(E=0;E<(int)((Emax-E0));E++){
      
      dtheta_spline= spline(x,dtheta[E], (int)(MAXSTEP-1), dtheta_dfirst[E] , dtheta_dlast[E] ,dtheta_spline);
      dtheta_temp2=splint(x,dtheta[E] , dtheta_spline, (int)(MAXSTEP-1), step, dtheta_temp2);
      fprintf(dthetaE,"%le\t",dtheta[E][step]);

    }
    fprintf(dthetaE,"\n");
  }
  
  */

  // calculate imaginary action integral
  
  double *x, *dtheta_spline;

  x=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );
  dtheta_spline=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) ); 

  for(j=0;j<(int)(MAXSTEP);j++){    
    x[j]=j;
  }

  for(step=0;step<MAXSTEP;step++){
    printf("at step = %lf Rx_coord = %lf Energy[%d]=%lf\n",x[step],Rx_coord[step],step,Energy[step]);   
  }


  double dstep=1.0e-5;   
  double dtheta_temp1=0.0, dtheta_temp2=0.0;
  double *Rx_coord_spline;
  Rx_coord_spline=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) ); 
  double Rx1=0, Rx2=0;
  Rx_coord_spline= spline(x, Rx_coord, MAXSTEP-1, 1e30, 1e30,Rx_coord_spline);
 
  // for(E=(int)(E0/ratio);E<(int)(Emax/ratio);E++){
  
  for(E=0;E<(int)((Emax-E0));E++){
        
    sl=find_sl_match(((E+E0)),dstep);
    sr=find_sr_match(((E+E0)),dstep);
        
    double s=0;
    double Energy_sl=0.0;
    double Energy_sr=0.0;

    //theta[E]=integrate_function(sl,sr,dtheta[E]);

    dstep=(sr-sl)*1.0e-3;
    
    for(s=sl;s<sr;s=s+(dstep)){

      dtheta_spline= spline(x,dtheta[E], (int)(MAXSTEP-1), 1e30, 1e30,dtheta_spline);      
      dtheta_temp2=splint(x,dtheta[E] , dtheta_spline, (int)(MAXSTEP-1), s+dstep, dtheta_temp2);
      dtheta_temp1=splint(x,dtheta[E] , dtheta_spline, (int)(MAXSTEP-1), s, dtheta_temp1);
      Rx2=splint(x,Rx_coord,Rx_coord_spline,(int) (MAXSTEP-1), s+dstep,Rx2);
      Rx1=splint(x,Rx_coord,Rx_coord_spline,(int) (MAXSTEP-1), s,Rx1);
      //      deltas=(Rx2-Rx1);
      deltas=fabs(Rx2-Rx1);
      //deltas=Rx_coord[(int)(s+dstep)]-Rx_coord[s];
      //      theta[E]+=(2*pi/h_planck)*dtheta_temp1*deltas*a_bohr; //integrazione con metodo dei rettangoli
      theta[E]+=(2*pi/h_planck)*(dtheta_temp2+dtheta_temp1)*deltas/2*a_bohr;//1.07;  //integrazione con metodo dei trapezi
      //theta[E]+=(2*pi/h_planck)*dtheta[E][s]*deltas*a_bohr;

    }

    fprintf(imactint,"%d %lf\n",E-Emax+E0,2*theta[E]);
    //    fprintf(imactint,"%d %lf\n",E-Emax+E0,theta[E]);

    
    Energy_sl=splint(xE,Energy, Energy_spline, MAXSTEP-1, (sl), Energy_sl);
    Energy_sr=splint(xE,Energy, Energy_spline, MAXSTEP-1, (sr), Energy_sr);
    
    printf("sl=%.7f sr=%.7f \tE_sl=%.3f\tE_sr=%.3f\t%d\t%le\t%le\n",sl,sr, Energy_sl, Energy_sr,(E+E0),(theta[E]/1.0e25), 1/(1+exp(2*theta[E])));
  }


  for (E=0;E<E0;E++){    
    Transmission_coeff[E]=0.;
    //    printf("E%d trans= %le \n",E,Transmission_coeff[E]);
   }

  for (E=E0;E<=Emax;E++){    
    Transmission_coeff[E]=1.0/(1+exp(2*theta[(int)(E-E0)]));
    //    printf("E%d trans= %le \n",E,Transmission_coeff[E]);
   }
  
  for (E=Emax;E<(2*Emax-E0);E++){
    Transmission_coeff[E]=1.0 - Transmission_coeff[(int)(2*Emax-E)];
    if( Transmission_coeff[E] > 1.0)  Transmission_coeff[E]=1;
    //     printf("E%d trans= %le \n",E,Transmission_coeff[E]);
   }
    
  for (E=(2*Emax-E0);E<E_dim;E++){
    Transmission_coeff[E]=1;
    //        printf("E%d trans= %le \n",E,Transmission_coeff[E]);
  }

  /*
  FILE *P_Eout;

  if((P_Eout=fopen("Transmission_coefficient.dat","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","Transmission_coefficient.dat");
    exit(1);
  }

  for(E=0;E<E_dim;E++){
    fprintf(P_Eout, "%d %le\n",E, Transmission_coeff[E]);
  }
  //  fclose(P_Eout);
  */
  fclose(dthetaE);
  fclose(imactint);

  //  exit(0);

  //  fclose(mueff_spline_result);
 
  //  for(i=0;i<Emax;i++){free(dtheta[i]);}free(dtheta);

  //  for(i=0;i<MAXSTEP;i++){free(dtheta[i]);}free(dtheta);
  //  free(theta);
  //  free(x);

  //  free(mueff_spline);

  
  return Transmission_coeff;
    
}



double find_sl_match(double E, double dstep){

  double step=0;int j;
  double sl=0;
  double *xE, Energy_temp=0.0, *Energy_spline;
 
  xE=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );
  
  Energy_spline=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );

  for(j=0;j<MAXSTEP;j++){
    xE[j]=j;
  }

  //  dstep=1.0e-5;

  Energy_spline= spline(xE,Energy, MAXSTEP-1, 1e30, 1e30,Energy_spline);
  
  Energy_temp=splint(xE, Energy,Energy_spline, MAXSTEP-1,step,Energy_temp);

  while((Energy_temp-E)<1.0e-6){

    step=step+dstep;
    sl=0;
    Energy_temp=splint(xE, Energy,Energy_spline, MAXSTEP-1,step,Energy_temp);
  }

  sl=step;

  return sl;
}


double find_sr_match(double E,double dstep){

  double step=MAXSTEP-1;int j;
  double sr=0;
  double *xE, Energy_temp=0.0, *Energy_spline;
 
  xE=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );

  Energy_spline=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );

  for(j=0;j<MAXSTEP;j++){
   
    xE[j]=j;
   
  }
 
  Energy_spline= spline(xE,Energy, MAXSTEP-1, 1e30, 1e30,Energy_spline);
  
  Energy_temp=splint(xE, Energy,Energy_spline, MAXSTEP-1,step,Energy_temp);

  //  dstep=1.0e-5;

  while((Energy_temp-E)<1.0e-6 ){

    step=step-dstep;      
    sr=0;
    Energy_temp=splint(xE, Energy,Energy_spline, MAXSTEP-1,step,Energy_temp);      
  }

  sr=step; 
  return sr;
}


double integrate_function(int Efirst, int Elast, double *function){
  void nrerror(char error_text[]);

  double EPS=1.0e-7;
  int JMAX=20;

  int j;
  double s=0,olds=0;

  double *x, *function2;
  x=(double *) malloc( MAXSTEP*sizeof(double ) );
  function2=(double *) malloc( MAXSTEP*sizeof(double ) ); 

  for(j=0;j<MAXSTEP;j++){

    x[j]=(Rx_coord[j]-Rx_coord[0]);
  }

  function2= spline(x, function, MAXSTEP-1, 1e30, 1e30,function2);

  /*  for(j=1;j<=JMAX;j++){
    st=trapzd(function,Efirst,Elast,j);
    // printf("s=%lf\n",st);
    s=(4*st-ost)/3.0;
    if(j>5)
      if(fabs(s-os)<EPS*fabs(os)||(s==0.0 && os==0.0)) return s;
    os=s;
    ost=st;
  }
  nrerror("Too many steps in routine qsimp");
  return 0.0;
  */

  Efirst=x[Efirst];
  Elast=x[Elast];

  for(j=1;j<=JMAX;j++){
    s=trapzd(x,function,function2,Efirst,Elast,j,olds,s);
    if(j>5)
      if(fabs(s-olds)<EPS*fabs(olds) ||(s==0.0 && olds==0.0))
	{ free(x);
	free(function2);
	return s;}
    olds=s;
  }

  nrerror("Too many steps in routine qsimp");



  return 0.0;

}

double trapzd(double *x,double *function,double *function2, int a,int b,int n, double olds, double s){

  double sum,z,del,tnm;
 
  int j,it;

  double ya=0,yb=0,y=0;

  

    if(n==1){
      ya=splint(x, function, function2, MAXSTEP-1, a, ya);
      yb=splint(x, function, function2, MAXSTEP-1, b, yb);
      return (s=0.5*(b-a)*(ya+yb));
    }else{
      it=pow(2,(n-2));
      tnm=it;
      del=(b-a)/tnm;
      z=a+0.5*del;
      sum=0.0;
      for(j=1;j<=it;j++){
	
	y=splint(x, function, function2, MAXSTEP-1, z, y);
	//	printf("spline[%lf]=%lf\n",z,y);
	sum=sum+y;
	z=z+del;
      }
      s=0.5*(olds+(b-a)*sum/tnm);
      return s;
    }
    
}

int max(int a, int b){

  if(a>=b)return a;
  else return b;

}


void deriv_gradient(void){

  FILE *grad_d,*grad;

  if((grad_d=fopen("gradient_derivative.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","kuni.txt");
    exit(1);
  }
 
  int N_spline=18;

  double **gradient_transposed, *tot1;
  double *x,*gradient_spline,grad1=0,grad2=0,grad0=0,grad3=0,*x1,*gradient_derivative_spline ;
  int j,i;
  double **gradient_derivative_temp,**gradient_transposed_temp ;
 
  x=(double *) malloc( MAXSTEP*sizeof(double ) );
  gradient_spline=(double *) malloc( MAXSTEP*sizeof(double ) );
 
  x1=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
 
  gradient_derivative_spline=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
 
  tot1=(double *) malloc(MAXSTEP*sizeof(double ) );
 
 
  gradient_transposed=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    gradient_transposed[j] = (double*)malloc((int)(MAXSTEP)*sizeof(double));
  }
 
  gradient_derivative=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    gradient_derivative[j] = (double*)malloc((int)(MAXSTEP)*sizeof(double));
  }
 
  gradient_derivative_temp=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    gradient_derivative_temp[j] = (double*)malloc((int)(MAXSTEP-N_spline)*sizeof(double));
  }

  gradient_transposed_temp=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    gradient_transposed_temp[j] = (double*)malloc((int)(MAXSTEP-N_spline)*sizeof(double));
  }
 
  gradient_transposed=transpose_double(gradient,3*ATOMS,MAXSTEP, gradient_transposed);
 

  if((grad=fopen("gradient.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","gradient.txt");
    exit(1);
  }

  for(i=0;i<MAXSTEP;i++){
    fprintf(grad,"%d\t",i);
    for(j=0;j<3*ATOMS;j++){
      fprintf(grad,"%lf\t",gradient[i][j]);
    }
    fprintf(grad,"\n");
  }

  fclose(grad);



  for(j=0;j<MAXSTEP;j++){
    x[j]=j;
  }
 
 
  for(j=0;j<MAXSTEP;j++){
    tot1[j]=0;
  }
 
  for(i=0;i<MAXSTEP;i++){
    for(j=0;j<3*ATOMS;j++){
      tot1[i]=tot1[i]+(gradient_transposed[j][i]*gradient_transposed[j][i]);
    }   
  }
 
  for(j=0;j<3*ATOMS;j++){
    for(i=0;i<MAXSTEP;i++){
      gradient_transposed[j][i]= gradient_transposed[j][i]/sqrt(tot1[i]);
    }    
  }
 
  printf("gradient transposed\n");
 
  for(i=0;i<MAXSTEP;i++){
    printf("%d\t",i);
    for(j=0;j<3*ATOMS;j++){
      printf("%lf\t",gradient_transposed[j][i]);
    }
    printf("\n");
  }
 
 
 
  for(i=0;i<3*ATOMS;i++){

    gradient_spline= spline(x,gradient_transposed[i], MAXSTEP-1, 1e30, 1e30,gradient_spline);
     
    grad1=-splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,0, grad1);
    grad2=-splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,1, grad2);
    deltas=Rx_coord[1]-Rx_coord[0];
    gradient_derivative[i][0]= (grad2-grad1)/deltas;
    
    grad2=-splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,2, grad2);
    deltas=Rx_coord[2]-Rx_coord[0];
    gradient_derivative[i][1]= (grad2-grad1)/deltas;

    grad1=splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,MAXSTEP-2, grad1);
    grad2=splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,MAXSTEP-1, grad2);
    deltas=Rx_coord[MAXSTEP-1]-Rx_coord[MAXSTEP-2];
    gradient_derivative[i][MAXSTEP-1]= (grad2-grad1)/deltas;
  
    grad1=splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,MAXSTEP-3, grad1);
    deltas=Rx_coord[MAXSTEP-1]-Rx_coord[MAXSTEP-3];
    gradient_derivative[i][MAXSTEP-2]= (grad2-grad1)/deltas;

    for(j=2;j<MAXSTEP2TS/2;j++){

      
      grad1=-splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,j-1, grad1);
      grad2=-splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,j+1, grad2);
   
      //prova a usare derivate a 5 punti
      
      grad0=-splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,j-2, grad0);
      grad3=-splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,j+2, grad3);
      deltas=Rx_coord[j-2]-Rx_coord[j+2];
      gradient_derivative[i][j]= (-grad3+8*grad2-8*grad1+grad0)/3/deltas;
     
      //gradient_derivative[i][j]= (grad2-grad1)/2/deltas;

    }

    for(j=MAXSTEP2TS/2;j<MAXSTEP-2;j++){

      grad1=splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,j-1, grad1);
      grad2=splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,j+1, grad2);      
      grad0=splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,j-2, grad0);
      grad3=splint(x, gradient_transposed[i],gradient_spline, MAXSTEP-1,j+2, grad3);
      deltas=Rx_coord[j-2]-Rx_coord[j+2];
      gradient_derivative[i][j]= (-grad3+8*grad2-8*grad1+grad0)/3/deltas;
     
      //gradient_derivative[i][j]= (grad2-grad1)/2/deltas;

    }
  }

  for(j=0;j<3*ATOMS;j++){
    for(i=0;i<MAXSTEP;i++){
      gradient_transposed[j][i]= gradient_transposed[j][i]*sqrt(tot1[i]);
    }   
  }
  // eliminating N_spline points from gr_der_temp vectors

  for(i=0;i<3*ATOMS;i++){
    for(j=0;j<(int)(MAXSTEP2TS/2-N_spline/2);j++){   
      gradient_derivative_temp[i][j]= gradient_derivative[i][j];
      gradient_transposed_temp[i][j]= gradient_transposed[i][j];
    }   
    for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<MAXSTEP-N_spline;j++){
      gradient_derivative_temp[i][j]= gradient_derivative[i][j+N_spline];
      gradient_transposed_temp[i][j]= gradient_transposed[i][j+N_spline];
    }   
  }

  for(j=0;j<MAXSTEP-N_spline;j++){
   
    x1[j]=j;
    if(j>=(int)(MAXSTEP2TS/2-N_spline/2)) x1[j]=j+N_spline;
  }
 
  //substituting central N_spline points with spline interpolations
  
 for(i=0;i<3*ATOMS;i++){  
   for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<(int)(MAXSTEP2TS/2+N_spline/2);j++){

     gradient_derivative_spline= spline(x1,gradient_derivative_temp[i], MAXSTEP-N_spline-1, 1e30, 1e30,gradient_derivative_spline);
     gradient_derivative[i][j]=splint(x1,gradient_derivative_temp[i],gradient_derivative_spline, MAXSTEP-N_spline-1,j,gradient_derivative[i][j]);
    
     gradient_spline= spline(x1,gradient_transposed_temp[i], MAXSTEP-N_spline-1, 1e30, 1e30,gradient_spline);
     gradient_transposed[i][j]=splint(x1, gradient_transposed_temp[i],gradient_spline, MAXSTEP-N_spline-1,j,gradient_transposed[i][j]);
   }      
 }
  
 
  for(i=0;i<MAXSTEP;i++){
    fprintf(grad_d,"%d\t",i);
    for(j=0;j<3*ATOMS;j++){
      fprintf(grad_d,"%lf\t",gradient_derivative[j][i]);
    }
    fprintf(grad_d,"\n");
  }

  fclose(grad_d);

  gradient=transpose_double(gradient_transposed,MAXSTEP,3*ATOMS ,gradient);


  free(x);free(x1);
  free(gradient_spline);free(gradient_derivative_spline);
  free(tot1);
  for(i=0;i<3*ATOMS;i++){free(gradient_transposed[i]);}free(gradient_transposed);
  for(i=0;i<3*ATOMS;i++){free(gradient_derivative_temp[i]);}free(gradient_derivative_temp);
  for(i=0;i<3*ATOMS;i++){free(gradient_transposed_temp[i]);}free(gradient_transposed_temp);
}


void deriv_position(void){

  FILE *pos_d;

  if((pos_d=fopen("positions_derivative.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","kuni.txt");
    exit(1);
  }
 
  int N_spline=18;

  double **position_transposed;
  double *x,*position_spline,pos0=0,pos1=0,pos2=0,pos3=0,*x1,*position_derivative_spline ;
  int j,i;
  double **position_derivative_temp,**position_transposed_temp ;
 
  x=(double *) malloc( MAXSTEP*sizeof(double ) );
  position_spline=(double *) malloc( MAXSTEP*sizeof(double ) );
 
  x1=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
 
  position_derivative_spline=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
 
 // tot1=(double *) malloc(MAXSTEP*sizeof(double ) );
 
 
  position_transposed=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    position_transposed[j] = (double*)malloc((int)(MAXSTEP)*sizeof(double));
  }
 
  position_derivative=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    position_derivative[j] = (double*)malloc((int)(MAXSTEP)*sizeof(double));
  }
 
  position_derivative_temp=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    position_derivative_temp[j] = (double*)malloc((int)(MAXSTEP-N_spline)*sizeof(double));
  }

  position_transposed_temp=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    position_transposed_temp[j] = (double*)malloc((int)(MAXSTEP-N_spline)*sizeof(double));
  }
 
  position_transposed=transpose_double(positions,3*ATOMS,MAXSTEP, position_transposed);
 
  for(j=0;j<MAXSTEP;j++){   
    x[j]=j;
  }
 

  for(i=0;i<3*ATOMS;i++){

    position_spline= spline(x,position_transposed[i], MAXSTEP-1, 1e30, 1e30,position_spline);
   
   
    pos1=splint(x, position_transposed[i],position_spline, MAXSTEP-1,0, pos1);
    pos2=splint(x, position_transposed[i],position_spline, MAXSTEP-1,1, pos2);
    
    position_derivative[i][0]= (pos2-pos1)/deltas;

    
    pos2=splint(x, position_transposed[i],position_spline, MAXSTEP-1,2, pos2);
    
    position_derivative[i][1]= (pos2-pos1)/2/deltas;

    pos1=splint(x, position_transposed[i],position_spline, MAXSTEP-1,MAXSTEP-2, pos1);
    pos2=splint(x, position_transposed[i],position_spline, MAXSTEP-1,MAXSTEP-1, pos2);
    
    position_derivative[i][MAXSTEP-1]= (pos2-pos1)/deltas;

    pos1=splint(x, position_transposed[i],position_spline, MAXSTEP-1,MAXSTEP-3, pos1);
    position_derivative[i][MAXSTEP-2]= (pos2-pos1)/2/deltas;

    for(j=2;j<MAXSTEP-2;j++){

      
      pos1=splint(x, position_transposed[i],position_spline, MAXSTEP-1,j-1, pos1);
      pos2=splint(x, position_transposed[i],position_spline, MAXSTEP-1,j+1, pos2);      
      pos0=splint(x, position_transposed[i],position_spline, MAXSTEP-1,j-2, pos0);
      pos3=splint(x, position_transposed[i],position_spline, MAXSTEP-1,j+2, pos3);
      
      position_derivative[i][j]= (-pos3+8*pos2-8*pos1+pos0)/12/deltas;
     
      //position_derivative[i][j]= (pos2-pos1)/2/deltas;

    }

 for(i=0;i<3*ATOMS;i++){

   for(j=0;j<(int)(MAXSTEP2TS/2-N_spline/2);j++){
   
     position_derivative_temp[i][j]= position_derivative[i][j];
     position_transposed_temp[i][j]= position_transposed[i][j];
   }
   

   for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<MAXSTEP-N_spline;j++){
   
     position_derivative_temp[i][j]= position_derivative[i][j+N_spline];
     position_transposed_temp[i][j]= position_transposed[i][j+N_spline];
   }
   
   
 }

 for(j=0;j<MAXSTEP-N_spline;j++){
   
   x1[j]=j;
   if(j>=(int)(MAXSTEP2TS/2-N_spline/2)) x1[j]=j+N_spline;
 }
 
 
 for(i=0;i<3*ATOMS;i++){
  
  
   for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<(int)(MAXSTEP2TS/2+N_spline/2);j++){
     position_derivative_spline= spline(x1,position_derivative_temp[i], MAXSTEP-N_spline-1, 1e30, 1e30,position_derivative_spline);
     position_derivative[i][j]=splint(x1, position_derivative_temp[i],position_derivative_spline, MAXSTEP-N_spline-1,j,position_derivative[i][j]);
    
     position_spline= spline(x1,position_transposed_temp[i], MAXSTEP-N_spline-1, 1e30, 1e30,position_spline);
     position_transposed[i][j]=splint(x1, position_transposed_temp[i],position_spline, MAXSTEP-N_spline-1,j,position_transposed[i][j]);
   }

      
 }

  }

  for(i=0;i<MAXSTEP;i++){
    fprintf(pos_d,"%d\t",i);
    for(j=0;j<3*ATOMS;j++){
      fprintf(pos_d,"%lf\t",position_derivative[j][i]);
    }
    fprintf(pos_d,"\n");
  }

  fclose(pos_d);

}

void interpolate_freqs(double **freqs){

double **freqs_transposed, *freq_spline, *x, **freqs_temp;
int N_spline=0;
int j,i;
x=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
freq_spline=(double *) malloc((int)(MAXSTEP-N_spline)*sizeof(double ) );

freqs_transposed=(double**)malloc((3*ATOMS) * sizeof(double*));
for(j=0;j<3*ATOMS;j++) {
freqs_transposed[j] = (double*)malloc((int)(MAXSTEP)*sizeof(double));
}

 freqs_temp=(double**)malloc((3*ATOMS) * sizeof(double*));
 for(j=0;j<3*ATOMS;j++) {
 freqs_temp[j] = (double*)malloc((int)(MAXSTEP-N_spline)*sizeof(double));
 }
 
 freqs_transposed=transpose_double(freqs,3*ATOMS,MAXSTEP, freqs_transposed);
 

for(i=0;i<3*ATOMS;i++){

   for(j=0;j<(int)(MAXSTEP2TS/2-N_spline/2);j++){
   
     freqs_temp[i][j]= freqs_transposed[i][j];
     
   }
   
   for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<MAXSTEP-N_spline;j++){
      freqs_temp[i][j]= freqs_transposed[i][j+N_spline];    
   }   
 }

 for(j=0;j<MAXSTEP-N_spline;j++){
   
   x[j]=j;
   if(j>=(int)(MAXSTEP2TS/2-N_spline/2)) x[j]=j+N_spline;
 }
 
 
 for(i=0;i<3*ATOMS;i++){
  
  
   for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<(int)(MAXSTEP2TS/2+N_spline/2);j++){
     freq_spline = spline(x,freqs_temp[i], MAXSTEP-N_spline-1, 1e30, 1e30,freq_spline);
     freqs_transposed[i][j]=splint(x, freqs_temp[i],freq_spline, MAXSTEP-N_spline-1,j,freqs_transposed[i][j]);
    
   
   }

      
 }
   
 freqs=transpose_double(freqs_transposed,MAXSTEP,3*ATOMS ,freqs);


}


void interpolate_hess(double **freqs){

  double **freqs_transposed, *freq_spline, *x, **freqs_temp;
  int N_spline=18;
  int j,i;
  x=(double *) malloc( (int)(MAXSTEP-N_spline)*sizeof(double ) );
  freq_spline=(double *) malloc((int)(MAXSTEP-N_spline)*sizeof(double ) );

  freqs_transposed=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    freqs_transposed[j] = (double*)malloc((int)(MAXSTEP)*sizeof(double));
  }

  freqs_temp=(double**)malloc((3*ATOMS) * sizeof(double*));
  for(j=0;j<3*ATOMS;j++) {
    freqs_temp[j] = (double*)malloc((int)(MAXSTEP-N_spline)*sizeof(double));
  }
 
  freqs_transposed=transpose_double(freqs,3*ATOMS,MAXSTEP, freqs_transposed);
 

  for(i=0;i<3*ATOMS;i++){
    for(j=0;j<(int)(MAXSTEP2TS/2-N_spline/2);j++){   
      freqs_temp[i][j]= freqs_transposed[i][j];
    }   
    for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<MAXSTEP-N_spline;j++){
      freqs_temp[i][j]= freqs_transposed[i][j+N_spline];
    }   
  }

  for(j=0;j<MAXSTEP-N_spline;j++){
    x[j]=j;
    if(j>=(int)(MAXSTEP2TS/2-N_spline/2)) x[j]=j+N_spline;
  }
  
  for(i=0;i<3*ATOMS;i++){  
    for(j=(int)(MAXSTEP2TS/2-N_spline/2);j<(int)(MAXSTEP2TS/2+N_spline/2);j++){
      freq_spline = spline(x,freqs_temp[i], MAXSTEP-N_spline-1, 1e30, 1e30,freq_spline);
      freqs_transposed[i][j]=splint(x, freqs_temp[i],freq_spline, MAXSTEP-N_spline-1,j,freqs_transposed[i][j]);
    }      
  }   
  freqs=transpose_double(freqs_transposed,MAXSTEP,3*ATOMS ,freqs);
}


void smooth_hess(void){
  
  int i,j,k;
  int dim;
  double ***hess_step_col;

  dim=3*ATOMS;
  
  hess_step_col=(double ***) malloc( dim*sizeof(double **) );for(step=0;step<dim;step++){  hess_step_col[step]=(double **) malloc( (int)(MAXSTEP)*sizeof(double *) );}
  for(step=0;step<dim;step++){ for(j=0; j <(int)(MAXSTEP); j++) {  hess_step_col[step][j]= (double*) malloc((int)(dim)*sizeof(double) );} }

  for(k=0;k<dim;k++){
    for(i=0;i<MAXSTEP;i++){
      for(j=0;j<dim;j++){       
	hess_step_col[k][i][j]=Hessian[i][k][j];       
      }
    }
    interpolate_hess(hess_step_col[k]);
  }
  
  for(k=0;k<dim;k++){
    for(i=0;i<MAXSTEP;i++){
      for(j=0;j<dim;j++){       
	Hessian[i][k][j]=hess_step_col[k][i][j];
      }
    }   
  }
 for(i=0;i<dim;i++){for(j=0;j<dim;j++){free(hess_step_col[i][j]);} free(hess_step_col[i]);} free(hess_step_col);

}
  

void grad_hess(void){

 // Hessian gradient (finite differences)

  int i,j;
  int dim;
  int step;

  dim=3*ATOMS;

  for(step=0;step<MAXSTEP;step++){     
    if(step==0){
      deltas=Rx_coord[step+1]-Rx_coord[step];
      for(j=0;j<(dim);j++){
	for(i=0;i<dim;i++){
	  Hessian_deriv[step][i][j]=-(Hessian[step+1][i][j]-Hessian[step][i][j])/(deltas);
	}
      }
    } else if(step==MAXSTEP-1){
      deltas=Rx_coord[step]-Rx_coord[step-1];
      for(j=0;j<(dim);j++){
	for(i=0;i<dim;i++){	 
	  Hessian_deriv[step][i][j]=(Hessian[step][i][j]-Hessian[step-1][i][j])/(deltas);
	}
      }
    } else if(step<MAXSTEP2TS/2){
     deltas=Rx_coord[step+1]-Rx_coord[step-1];
     for(j=0;j<(dim);j++){
       for(i=0;i<dim;i++){	 
	 Hessian_deriv[step][i][j]=(Hessian[step+1][i][j]-Hessian[step-1][i][j])/(2*deltas);
       }
     }
    } else if(step>MAXSTEP2TS/2){
      deltas=Rx_coord[step+1]-Rx_coord[step-1];
      for(j=0;j<(dim);j++){
	for(i=0;i<dim;i++){
	  Hessian_deriv[step][i][j]=(Hessian[step+1][i][j]-Hessian[step-1][i][j])/(2*deltas);	 
	}
      }
    } else if(step==MAXSTEP2TS/2){
      deltas=(Rx_coord[step+1]-Rx_coord[step]);
     for(j=0;j<(dim);j++){
       for(i=0;i<dim;i++){	 
	 // five points 
	 //	 Hessian_deriv[step][i][j]=(-Hessian[step+2][i][j]+8*Hessian[step+1][i][j]-8*Hessian[step-1][i][j]+Hessian[step-2][i][j])/(12*deltas);
	 // four points forward
	 //	 Hessian_deriv[step][i][j]=(-Hessian[step+2][i][j]+6*Hessian[step+1][i][j]-3*Hessian[step][i][j]-2*Hessian[step-1][i][j])/(6*deltas);
	 // centered
	 Hessian_deriv[step][i][j]=(Hessian[step+1][i][j]-Hessian[step-1][i][j])/(4*deltas);
       }
     }
    }
  }


  /*
  FILE *Lderiv; 
  if((Lderiv=fopen("Lderiv.txt","w"))==NULL) {
    printf("********IMPOSSIBILE APRIRE IL FILE %s*************\n","Lderiv.txt");
    exit(1);
  }
  for(step=0;step<MAXSTEP;step++){
    fprintf(Lderiv,"step=%d\n\n",step);
    for(j=0;j<(int)(dim);j++){
      for(i=0;i<(int)(dim);i++){
	fprintf(Lderiv,"%lf\t",Hessian_deriv[step][j][i]);
      }  
      fprintf(Lderiv,"\n");
    }
    fprintf(Lderiv,"\n");
  }
  fclose(Lderiv);
  */
}


double calc_var_corr(int dim, double Trif, double **frequencies_save){

  // evaluation of variational minimum of reactant flux
  // at assigned T

  double *qvib,*qrot;
  double varc=1.;
  double qrfact;
  double minflux;
  double stepflux;
  double tsflux;
  int i;
  //  int stepmin;

  qvib=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );
  qrot=(double *) malloc( (int)(MAXSTEP)*sizeof(double ) );

  for(step=0;step<MAXSTEP;step++){
    qvib[step]=1.;
    qrot[step]=1.;
  }

  //  qrfact=sqrt(pow((8*pi*boltzmann*Trif/h_planck/h_planck),3));
  qrfact=8*pow(pi,2)*sqrt(pow((2*pi*boltzmann*Trif/h_planck/h_planck),3));

  for(step=0;step<MAXSTEP;step++){    
    for(i=0;i<dim-7;i++){
      if(frequencies_save[step][i]>0.){
	qvib[step]=qvib[step]/(1.-exp(-h_planck*frequencies_save[step][i]*c_light_cm_s/boltzmann/Trif));
      }
      qrot[step]=qrfact*sqrt(eigen[step][1]*eigen[step][2]*eigen[step][3]*pow((a_mass*a_bohr*a_bohr),3));
    }
    //    if(step==MAXSTEP/2) printf("at %lf K qvib of tstate is %le\n",Trif,qvib[step]);
    //    if(step==MAXSTEP/2) printf("at %lf K qrot of tstate is %le\n",Trif,qrot[step]);
  }

  tsflux=qvib[MAXSTEP2TS/2]*qrot[MAXSTEP2TS/2]*exp(-Energy[MAXSTEP2TS/2]*1000./349.75/1.987/Trif);
  minflux=tsflux;
  stepmin=MAXSTEP2TS/2;

  for(step=0;step<MAXSTEP;step++){    
    stepflux=qvib[step]*qrot[step]*exp(-Energy[step]*1000./349.75/1.987/Trif);
    if (stepflux < minflux) {
      minflux=stepflux;
      stepmin=step;
    }
  }

  varc=tsflux/minflux;

  //   printf("at %lf K vcorr is %lf at step %d exp is: %le energy is: %lf  \n",Trif,varc,stepmin,exp(-Energy[stepmin]*1000./349.75/1.987/Trif),Energy[stepmin]*1000./349.75);

  free(qvib);
  free(qrot);

  return varc;
}
