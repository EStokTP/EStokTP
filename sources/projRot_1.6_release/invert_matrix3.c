#include <math.h>
#include "nrutil.h"
#include <stdio.h>
#include <stdlib.h>
#define TINY 1.0e-20                       // A small number.

void ludcmp(double **a, int n, int *indx, double d);
void lubksb(double **a, int n, int *indx, double *b);

double **inverse_matrix(double **a, int N, double **y){
  
  double d=1,*col;
  int i,j,*indx;

  col=(double *) malloc( (int)(N+1)*sizeof(double ) );
  indx=(int *) malloc( (int)(N+1)*sizeof(int ) );
 
  ludcmp(a,N,indx,d);                  //Decompose the matrix just once.


  for(j=1;j<=N;j++) {                   //Find inverse by columns.
    for(i=1;i<=N;i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(a,N,indx,col);
    for(i=1;i<=N;i++) y[i][j]=col[i];
  }


  return y;
}

void ludcmp(double **a, int n, int *indx, double d)
//Given a matrix a[1..n][1..n] , this routine replaces it by the LU decomposition of a rowwise
//permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
//indx[1..n] is an output vector that records the row permutation effected by the partial
//pivoting; d is output as ±1 depending on whether the number of row interchanges was even
//or odd, respectively. This routine is used in combination with lubksb to solve linear equations
//or invert a matrix.
{
     int i,imax=0,j,k;
     double big,dum,sum,temp;

     double *vv;   

     vv=(double *) malloc( (int)(n+1)*sizeof(double ) );          //vv stores the implicit scaling of each row.
  //cc-replaced     vv=(double *) malloc( (int)(n)*sizeof(double ) );          //vv stores the implicit scaling of each row.
     // vv=vector(0,n);

     d=1.0;                                //No row interchanges yet.
     for (i=1;i<=n;i++) {                   //Loop over rows to get the implicit scaling informa-
       big=0.0;                               // tion.
         for (j=1;j<=n;j++)
              if ((temp=fabs(a[i][j])) > big) big=temp;
         if (big == 0.0) nrerror("Singular matrix in routine ludcmp");//No nonzero largest element.
         vv[i]=1.0/big;                     //Save the scaling.
     }
     for (j=1;j<=n;j++) {                   //This is the loop over columns of Crout’s method.
       for (i=1;i<j;i++) {                //This is equation (2.3.12) except for i = j .
	 sum=a[i][j];
	 for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
	 a[i][j]=sum;
       }
       big=0.0;                           //Initialize for the search for largest pivot element.
       for (i=j;i<=n;i++) {               //This is i = j of equation (2.3.12) and i = j + 1 . . . N
	 sum=a[i][j];                     //     of equation (2.3.13).
	 for (k=1;k<j;k++) sum -= a[i][k]*a[k][j];
	 a[i][j]=sum;
	 if ( (dum=vv[i]*fabs(sum)) >= big) {  //Is the figure of merit for the pivot better than the best so far?
	   big=dum;
	   imax=i;
	 }
       }
       if (j != imax) {                  //Do we need to interchange rows?
	 for (k=1;k<=n;k++) {         //Yes, do so...
	   dum=a[imax][k];
	   a[imax][k]=a[j][k];
	   a[j][k]=dum;
	 }
	 d = -(d);                  //...and change the parity of d.
	 vv[imax]=vv[j];              //Also interchange the scale factor.
									  }
      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=TINY;
      // If the pivot element is zero the matrix is singular (at least to the precision of the
      //algorithm). For some applications on singular matrices, it is desirable to substitute
      //TINY for zero.
      if (j != n) {                     //Now, finally, divide by the pivot element.
           dum=1.0/(a[j][j]);
           for (i=j+1;i<=n;i++) a[i][j] *= dum;
      }
     }                                     //Go back for the next column in the reduction.


  free_vector(vv,1,n);
}


void lubksb(double **a, int n, int *indx, double *b)
//Solves the set of n linear equations A · X = B . Here a[1..n][1..n] is input, not as the matrix
//A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
//as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
//B , and returns with the solution vector X . a, n, and indx are not modiﬁed by this routine
//and can be left in place for successive calls with diﬀerent right-hand sides b. This routine takes
//into account the possibility that b will begin with many zero elements, so it is eﬃcient for use
//in matrix inversion.
{
     int i,ii=0,ip,j;
     double sum;
     for (i=1;i<=n;i++) {                //When ii is set to a positive value, it will become the
       ip=indx[i];                       //index of the ﬁrst nonvanishing element of b. We now
       sum=b[ip];                        //do the forward substitution, equation (2.3.6). The
       b[ip]=b[i];                       //only new wrinkle is to unscramble the permutation
       if (ii)                           //as we go.
	 for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
       else if (sum) ii=i;               //A nonzero element was encountered, so from now on we
       b[i]=sum;                         // will have to do the sums in the loop above.
     }
     for (i=n;i>=1;i--) {                //Now we do the backsubstitution, equation (2.3.7).
         sum=b[i];
         for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
         b[i]=sum/a[i][i];               //Store a component of the solution vector X .
     }                                   //All done!
}


