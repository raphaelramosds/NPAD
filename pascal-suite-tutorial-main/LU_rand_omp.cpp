// LU decomposition with Doolittle's algorithm from 
// Philip Wallstedt http://www.sci.utah.edu/~wallstedt
// with random initialization, partial pivoting, and OpenMP parallelization by Samuel Xavier-de-Souza

#include <time.h>

#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
using namespace std;

// These functions print matrices and vectors in a nice format
void coutMatrix(int d,double*m){
   cout<<'\n';
   for(int i=0;i<d;++i){
      for(int j=0;j<d;++j)cout<<setw(14)<<m[i*d+j];
      cout<<'\n';
   }
}

void coutVectorUnsignInt(int d,unsigned int*v){
   cout<<'\n';
   for(int j=0;j<d;++j)cout<<setw(14)<<v[j];
   cout<<'\n';
}

void coutVector(int d,double*v){
   cout<<'\n';
   for(int j=0;j<d;++j)cout<<setw(14)<<v[j];
   cout<<'\n';
}

// The following compact LU factorization schemes are described
// in Dahlquist, Bjorck, Anderson 1974 "Numerical Methods".
//
// S and D are d by d matrices.  However, they are stored in
// memory as 1D arrays of length d*d.  Array indices are in
// the C style such that the first element of array A is A[0]
// and the last element is A[d*d-1].
//
// These routines are written with separate source S and
// destination D matrices so the source matrix can be retained
// if desired.  However, the compact schemes were designed to
// perform in-place computations to save memory.  In
// other words, S and D can be the SAME matrix.  

// Doolittle's algorithm with partial pivoting - pivot is the permutation vector | by SXdS 
void Doolittle_pivot(int d,double*S,double*D,unsigned int*pivot){
   int k,i,j,p,max;
   double sum;
#pragma omp parallel private(sum,k,i,j,p) default(none) shared(d,S,D,pivot,max,cout)
{
   #pragma omp master
   cout<<"Number of active threads is "<<omp_get_num_threads()<<"\n";
   #pragma omp for
   for(k=0;k<d;++k)
      pivot[k]=k;
   for(k=0;k<d;++k){
      #pragma omp single
      {            
      for(max=k,p=k;p<d;++p)
         if(abs(D[pivot[p]*d+k])>abs(D[pivot[max]*d+k]))
            max=p;
      unsigned int tmp = pivot[k];
      pivot[k]=pivot[max];
      pivot[max]=tmp;
      }
      #pragma omp for schedule(static,10)
      for(j=k;j<d;++j){ // Compute and store U's row k in place
         for(sum=0., p=0;p<k;++p)
            sum+=D[pivot[k]*d+p]*D[pivot[p]*d+j];
         D[pivot[k]*d+j]=(S[pivot[k]*d+j]-sum); // not dividing by diagonals
      }
      #pragma omp for schedule(static,10)
      for(i=k+1;i<d;++i){ // Compute and store L's column k in place
         for(sum=0, p=0;p<k;++p)
            sum+=D[pivot[i]*d+p]*D[pivot[p]*d+k];
         D[pivot[i]*d+k]=(S[pivot[i]*d+k]-sum)/D[pivot[k]*d+k];
      }
   }
} // end parallel region
} // enf function

void solveDoolittle_pivot(int d,double*LU,double*b,double*x,unsigned int*pivot){
   double y[d];
   double sum;
   int k,i;
   for(i=0;i<d;++i){
      for(sum=0., k=0;k<i;++k)
         sum+=LU[pivot[i]*d+k]*y[k];
      y[i]=(b[pivot[i]]-sum); // not dividing by diagonals
   }
   for(i=d-1;i>=0;--i){
      for(sum=0., k=i+1;k<d;++k)
         sum+=LU[pivot[i]*d+k]*x[k];
      x[i]=(y[i]-sum)/LU[pivot[i]*d+i];
   }
}

int main(int argc, char **argv){

   double start,finish,enlapsed;
   unsigned int n = strtol(argv[1],NULL,10);
   double *A = (double *) malloc(n*n*sizeof(double));
   double *LU = (double *) malloc(n*n*sizeof(double)); 
   double *b = (double *) malloc(n*sizeof(double));
   double *x = (double *) malloc(n*sizeof(double));
   unsigned int *pivot = (unsigned int *) malloc(n*sizeof(unsigned int));
   
   // Arbitrary size and random inicialization
   if(argc==3)
      srand48(strtol(argv[2],NULL,10)); // same seed for the pseudorandom number generation
   else
      srand48(time(NULL)); // different seed for the pseudorandom number generation
   for(int i=0;i<n;i++){
      for(int j=0;j<n;j++)
         A[i*n+j] = drand48()*drand48(); // symetric for use with Cholesky
      b[i] = drand48();
   }

   // cout<<"A";
   // coutMatrix(n,A);
   start=omp_get_wtime();
   Doolittle_pivot(n,A,LU,pivot);
   finish=omp_get_wtime();
   cout<<"LU decomposition took "<<finish-start<<" seconds\n";
   //cout<<"Doolittle_pivot";
   //coutVectorUnsignInt(n,pivot);
   //coutMatrix(n,LU);
   solveDoolittle_pivot(n,LU,b,x,pivot);
   //cout<<"solveDoolittle_pivot";
   // coutVector(n,x);

   // Computes the mean square roots of Ax - b for validation
   double sum = 0.;
   for(int i=0;i<n;i++){
      for(int j=0;j<n;j++)
         b[i] -= A[i*n+j]*x[j];
      sum+=b[i]*b[i];
   }
   cout<<"MSR of (Ax - b) is "<<setw(14)<<sum/n<<'\n';

}



