#include <Accelerate/Accelerate.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main()
{

    double Mr,Mi,alpha,k,kend,kabs; 
    int i,j,n1,n2,n3,ii,jj;
    FILE *fp;
    char jobz,uplo;
    int lda,lwork,info;
    double *M,*rwork,*work,*w,kx,ky,kz;
    double Rx,Ry,Rz,Mrr,Mii,theta;
    double *MM0,*Tx,*Txy1,*Txy2,*Tz;
	double *xx,*yy,*zz,dx,dy,dz,dx2,dy2,dz2,dd,dd2;
	int iu,ju,ib,jb; 


	const double aa = 4.1355/sqrt(3.0);
	const double cc = 28.6152/3.0;
	const double AA1 = 2.26/cc;
	const double AA2 = 3.33*(2.0/3.0)/aa;
	const double BB1 = 6.86/(cc*cc);
	const double BB2 = 44.5/(aa*aa);
	const double CC1 = 5.74/(cc*cc);
	const double CC2 = 30.4/(aa*aa);
	const double M0 = -0.28 -4*BB2 -2*BB1;

	const int Nx = 2;
	const int Ny = 10;
	const int Nz = 10;
	const int NSyz = Ny*Nz;
	const int Nl = NSyz*4;
	const int ND = Nx*Ny*Nz*4;

	const int Nm = 4;

        const double pi=3.1415926536;
//printf("no error!"); 	
    

	MM0=(double *)malloc(Nm*Nm*2*sizeof(double));
	Tx=(double *)malloc(Nm*Nm*2*sizeof(double));
	Tz=(double *)malloc(Nm*Nm*2*sizeof(double));
	Txy1=(double *)malloc(Nm*Nm*2*sizeof(double));
	Txy2=(double *)malloc(Nm*Nm*2*sizeof(double));
	


	#include "create_TT_matrices.c"



	xx=(double *)malloc(ND*sizeof(double));
	yy=(double *)malloc(ND*sizeof(double));
	zz=(double *)malloc(ND*sizeof(double));

	for(i=0;i<ND;i++)
	{
		iu = i/4;
		xx[i] = (iu/(Ny*Nz))*aa+((iu%Ny+1)%2)*0.5*aa;
		yy[i] = (iu%Ny)*sqrt(3.0)*0.5*aa;
		zz[i] = (iu%(Ny*Nz))/Ny*cc;
	//if(zz[i]<1.01*2*cc && zz[i]>0.99*2*cc)   printf("%.3f\t%.3f\n",xx[i],yy[i]);
	}

  
    lwork=2*Nl-1;
   
    M=(double*)malloc(Nl*Nl*2*sizeof(double));
    
    for(kx=0.0;kx<2*pi/aa;kx+=pi/(aa*100))
    {         

		#include "MM.c"

               w=(double*)malloc(Nl*sizeof(double));
               work=(double*)malloc(2*lwork*sizeof(double));
               rwork=(double*)malloc((3*Nl-2)*sizeof(double));
 
               jobz='N';   
               uplo='U'; 
               lda=Nl;

//               ZHEEV(&jobz,&uplo,&Nl,M,&lda,w,work,&lwork,rwork,&info);
               zheev_(&jobz,&uplo,&Nl,M,&lda,w,work,&lwork,rwork,&info);
             
               free(work);
               free(rwork);
              
               printf("%5f\t",kx);
               for(i=0;i<Nl;i++)
               {
                     printf("%5f\t", w[i]);
               }
               printf("\n"); 
             
               free(w); 

     }

   
     free(M); 


	free(xx);
	free(yy);
	free(zz);

	free(MM0);
	free(Tx);
	free(Tz);
	free(Txy1);
	free(Txy2);



     return 0;
   
       
}

