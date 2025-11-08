#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image.h"
#include "../c_libs/stb_image_write.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

//finding norm of *x
double norm(int n,double *x){
  double s=0;
  for(int i=0;i<n;i++)
  s = s + x[i]*x[i];
  return sqrt(s);
  }
  
  //finding a new matrix by multiplying two matrices A and v
  void mad(int m,int n,double *A,double *v,double *y){
  for(int i=0;i<m;i++){
  double s=0;
  double *r=A+i*n;
  for(int j=0;j<n;j++)
  s = s + r[j]*v[j];
  y[i]=s;
  }}
  
  //finding a new matrix by multiplying two matrices A^t and u
  void man(int m,int n,double *A,double *u,double *z){
  for(int j=0;j<n;j++){
  double s=0;
  for(int i=0;i<m;i++)
  s = s + A[i*n+j]*u[i];
  z[j]=s;
  }
  }
  
  //taking a random vector v and applying power iteration to get eigen vector corresponding to highest singular value and finding singular value and singular vectors from that.
  double SVD(int m,int n,double *A,double *u,double *v,int p,double tol){
   double *y=calloc(m,sizeof(double));
   double *z=calloc(n,sizeof(double));
    for(int j=0;j<n;j++)
    v[j]=(rand()%1000)/1000.0;
    double ps=0;
    double s=0;
    for(int k=0;k<p;k++){
        mad(m,n,A,v,y);
        s=norm(m,y);
        if(s==0)
        break;
        man(m,n,A,y,z);
        double nz=norm(n,z);
        if(nz==0)
        break;
        for(int j=0;j<n;j++)
        v[j]=z[j]/nz;
        if(fabs(s-ps)<tol*(1+fabs(ps)))
        break;
        ps=s;
    }
    mad(m,n,A,v,y);
    for(int i=0;i<m;i++)
    u[i]=y[i]/s;  //computing u from y and s
    free(y);
    free(z); 
    return s;
}

//subtracting matrix corresponding to k singular values from complete matrix to get a new matrix.This process is called as deflation.
void submatrix(int m,int n,double *A,double *u,double *v,double s){
    for(int i=0;i<m;i++){
        double t=u[i]*s;
        double *r=A+i*n;
        for(int j=0;j<n;j++)
        r[j] = r[j] - t*v[j];
    }
}


void addmatrix(int m,int n,double *C,double *u,double *v,double s){
    for(int i=0;i<m;i++){
        double t=u[i]*s;
        double *r=C+i*n;
        for(int j=0;j<n;j++)
        r[j] = r[j] + t*v[j];
    }
}

//Finding Frobenius absolute error
double FrobeniusError(double *A,double *C,int m,int n){
    double sum = 0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            double d = A[i*n+j] - C[i*n+j];
            sum = sum + d*d;
        }
    }
    return 255.0 * sqrt(sum);
}

// storing eigen vectors through **U and **Vv 
void saveCompressed(double *A,int m,int n,int K,char *fname){

clock_t start = clock();

double **U=malloc(K*sizeof(double*));
double **Vv=malloc(K*sizeof(double*));
double *S=malloc(K*sizeof(double)); //storing singular values through *S
    double *Anew=malloc(m*n*sizeof(double)); 
    memcpy(Anew,A,m*n*sizeof(double));


    for(int k=0;k<K;k++){
    U[k]=calloc(m,sizeof(double)); //storing array of size m*1 in U
    Vv[k]=calloc(n,sizeof(double)); //storing array of size n*1 in Vv
    }

//if S[k] is very small i.e approaching to 0 then break
    for(int k=0;k<K;k++){
        S[k]=SVD(m,n,Anew,U[k],Vv[k],400,1e-6);
        if(S[k]<1e-8){   
        K=k;
        break;
        }
        submatrix(m,n,Anew,U[k],Vv[k],S[k]);
    }

//C is the matrix for reconstructed image
    double *C=calloc(m*n,sizeof(double));
    for(int k=0;k<K;k++)
    addmatrix(m,n,C,U[k],Vv[k],S[k]);

    double err = FrobeniusError(A,C,m,n);
    printf("Frobenius Error for K=%d : %f\n",K,err);

//forming *out to store m*n elements in it
    unsigned char *out=malloc(m*n);
    for(int i=0;i<m;i++)
    for(int j=0;j<n;j++){
        double x=C[i*n+j];
        if(x<0)
        x=0;
        if(x>1)
        x=1; 
        out[i*n+j]=(unsigned char)(x*255); 
    }

//input jpg given and we will get output jpg 
    stbi_write_jpg(fname,n,m,1,out,90);
    
    clock_t end = clock();
    double time_taken = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time taken for SVD + reconstruction for K=%d : %f seconds\n",K,time_taken);

    free(out);
    free(C);
    free(Anew);
    for(int k=0;k<K;k++){
    free(U[k]);
    free(Vv[k]);
    }
    free(U);
    free(Vv);
    free(S);
}

//using argc and **argv to get output jpg
int main(int argc,char **argv){
    if(argc<3){
        printf("Usage: %s inputImage outputPrefix\n",argv[0]);
        return 0;
    }

    int w,h,ch;
    unsigned char *img=stbi_load(argv[1],&w,&h,&ch,0);
    if(!img){
    printf("Load failed\n");
    return 0;
    }

    int m=h,n=w;
    double *A=malloc(m*n*sizeof(double)); //A is a m*n matrix

    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++){
            int idx=(i*n+j)*ch;
            double r=img[idx]/255.0;
            double g=(ch>1?img[idx+1]/255.0:r);
            double b=(ch>2?img[idx+2]/255.0:r);
            A[i*n+j]=0.299*r+0.587*g+0.114*b;
        }

    stbi_image_free(img);

//getting 4 different output jpg for 4 different values of k.
    char f1[100],f2[100],f3[100],f4[100];
    sprintf(f1,"%s_k5.jpg",argv[2]);
    sprintf(f2,"%s_k20.jpg",argv[2]);
    sprintf(f3,"%s_k50.jpg",argv[2]);
    sprintf(f4,"%s_k100.jpg",argv[2]);

    saveCompressed(A,m,n,5,f1);
    saveCompressed(A,m,n,20,f2);
    saveCompressed(A,m,n,50,f3);
    saveCompressed(A,m,n,100,f4);

    free(A);
    return 0;
}


 
  
