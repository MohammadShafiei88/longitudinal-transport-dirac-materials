/* #include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
main()
{
	int i,j,i1,j1,k,Nl,ND,Ly,l,v;
	int flag,info1,lda1,M1,N1,*ipiv1,lwork1;             j:for site l:number of sites
	double   t,Lso,e,eta,Lv;
	double   *d,*Dd,*A,*B,*green,*green_up_up,*green_down_up,*green_up_down,*green_down_down,*work1;
	double   *sigma_leftlead,*sigma_rightlead,*T1_leftlead,*T1_rightlead,*T2_leftlead,*T2_rightlead;
	v=ND-Ly;
	Nl=8;
	t=1.0;
	Lso=((0.3)*t)/(3*(3^0.5));
	Lv=(0.2)*t;
	eta=0.00001;
	e=-4.0;
	l=Nl;
	green_up_up=(double*)malloc(ND*ND*2*sizeof(double));
	green_up_down=(double*)malloc(ND*ND*2*sizeof(double));
	green_down_up=(double*)malloc(ND*ND*2*sizeof(double));
	green_down_down=(double*)malloc(ND*ND*2*sizeof(double));
	Dd=(double *)malloc(N1*N1*2*sizeof(double));
  initialize_green function  */

for(i=0;i<ND;i++)
{
	for(j=0;j<ND;j++)
	{
		green[2*(j*ND+i)]=0.0;
		green[2*(j*ND+i)+1]=0.0;
	}
}


//d =(double *)malloc(Nl*Nl*2*sizeof(double));
for(i=0;i<ND;i++)
{
	for(j=0;j<ND;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j];
	dy = yy[i] - yy[j];
	dz = zz[i] - zz[j];

	dx2 = dx*dx;
	dy2 = dy*dy;
	dz2 = dz*dz;

	dd2 = dx2+dy2+dz2;
	dd = sqrt(dd2);
	

	if(i==j)
		{  green[2*(j*(ND)+i)] += e;    green[2*(j*(ND)+i)+1] += eta;  }

	if(dd<0.01*aa)
	{
		{      green[2*(j*(ND)+i)] += -MM0[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -MM0[2*(jb*(4)+ib)+1];     }
	}
    if(j>i)
    {
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  green[2*(j*(ND)+i)] += -Tx[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Tx[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  green[2*(j*(ND)+i)] += -Txy1[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Txy1[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  green[2*(j*(ND)+i)] += -Txy2[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Txy2[2*(jb*(4)+ib)+1];  }
		 }


	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      green[2*(j*(ND)+i)] += -Tz[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Tz[2*(jb*(4)+ib)+1];     }
	}

}	
}

    if(j<i)
    {
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  green[2*(j*(ND)+i)] += -Txd[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Txd[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  green[2*(j*(ND)+i)] += -Txy1d[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Txy1d[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  green[2*(j*(ND)+i)] += -Txy2d[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Txy2d[2*(jb*(4)+ib)+1];  }
		 }


	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      green[2*(j*(ND)+i)] += -Tzd[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Tzd[2*(jb*(4)+ib)+1];     }
	}

}	
}

}
}

/* end of matrix green  */

/********************disordrer?*******************/

/* ALREADY COMMENTED OUT
	  for(i=0;i<ND;i++)
		 { green[2*((3*i)*ND+(3*i))]+=Er[i];
		  green[2*((3*i+1)*ND+(3*i+1))]+=Er[i];
		  green[2*((3*i+2)*ND+(3*i+2))]+=Er[i]; }
*/
/******************************** end **********************/

#include"add_sigma_to_green.c"

/* ALREADY COMMENTED OUT 
  printf("HDG matrix:\n");
	for(i=0;i<ND_T;i++)
	  {for(k=0;k<ND_T;k++)
		 { printf("%.4f  %u %u   %.5f\t%.5f\t   %.5f\t%.5f\n",e,i,k,green_up_up[2*(k*ND_T+i)],green_up_up[2*(k*ND_T+i)+1],green_down_down[2*(k*ND_T+i)],green_down_down[2*(k*ND_T+i)+1]);
		  }
	   }

printf("\n\n");
*/

lda1=ND;
N1=ND;
M1=ND;
lwork1=33*ND;
work1=(double *)malloc(lwork1*2*sizeof(double));
ipiv1=(int *)malloc(ND*sizeof(int));
//ZGETRF(&M1,&N1,green,&lda1,ipiv1,&info1);
//ZGETRI(&N1,green,&lda1,ipiv1,work1,&lwork1,&info1);
zgetrf_(&M1,&N1,green,&lda1,ipiv1,&info1);
zgetri_(&N1,green,&lda1,ipiv1,work1,&lwork1,&info1);
free(work1);
free(ipiv1);

/* ALREADY COMMENTED OUT
	lda1=ND_T;
	N1=ND_T;
	M1=ND_T;
	lwork1=33*ND_T;
	work1=(double *)malloc(lwork1*2*sizeof(double));
	ipiv1=(int *)malloc(ND_T*sizeof(int));
	ZGETRF(&M1,&N1,green_down_down,&lda1,ipiv1,&info1);
	ZGETRI(&N1,green_down_down,&lda1,ipiv1,work1,&lwork1,&info1);
	free(work1);
	free(ipiv1);

*/

/* ALREADY COMMENTED OUT
 printf("HD matrix:\n");
	for(i=0;i<ND_T;i++)
		{for(k=0;k<ND_T;k++)
		   { printf("%.4f  %u %u   %.5f\t%.5f\t   %.5f\t%.5f\n",e,i,k,green_up_up[2*(k*ND_T+i)],green_up_up[2*(k*ND_T+i)+1],green_down_down[2*(k*ND_T+i)],green_down_down[2*(k*ND_T+i)+1]);
			}
		}

   printf("\n\n");

*/

/* ALREADY COMMENTED OUT
   return 0;
  }
*/
