/* #include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
main()
{
	int i,j,i1,j1,k,(Nl+2),ND,Ly,l,v;
	int flag,info1,lda1,M1,N1,*ipiv1,lwork1;             j:for site l:number of sites
	double   t,Lso,e,eta,Lv;
	double   *d,*Dd,*A,*B,*green,*green,*green_down_up,*green_up_down,*green_down_down,*work1;
	double   *sigma_leftlead,*sigma_rightlead,*T1_leftlead,*T1_rightlead,*T2_leftlead,*T2_rightlead;
	v=ND-Ly;
	(Nl+2)=8;
	t=1.0;
	Lso=((0.3)*t)/(3*(3^0.5));
	Lv=(0.2)*t;
	eta=0.00001;
	e=-4.0;
	l=(Nl+2);
	green=(double*)malloc(ND*ND*2*sizeof(double));
	green_up_down=(double*)malloc(ND*ND*2*sizeof(double));
	green_down_up=(double*)malloc(ND*ND*2*sizeof(double));
	green_down_down=(double*)malloc(ND*ND*2*sizeof(double));
	Dd=(double *)malloc(N1*N1*2*sizeof(double));
  initialize_green function

	printf("HD matrix:\n");
	for(i=0;i<ND;i++)
	  {for(k=0;k<ND;k++)
		 { printf("%u %u   %.4f\t%.4f\t   %.4f\t%.4f\n",i,k,green[2*(k*ND+i)],green[2*(k*ND+i)+1],green_down_down[2*(k*ND+i)],green_down_down[2*(k*ND+i)+1]);
		  }
	   }

printf("\n\n");*/

#include "create_TT_matrices.c"
for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_rightlead[2*(j*(Nl)+i)]=0.0;
		T1_rightlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T2_rightlead[2*(j*(Nl)+i)]=0.0;
		T2_rightlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		sigma_rightlead[2*(j*(Nl)+i)]=0.0;
		sigma_rightlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


/************ matrix A of left lead using for T2_rightlead **************/
for(i=0;i<Nl;i++)
	for(j=0;j<Nl;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j] + aa;	 /*  for left lead  */
	dy = yy[i] - yy[j];
	dz = zz[i] - zz[j];

	dx2 = dx*dx;
	dy2 = dy*dy;
	dz2 = dz*dz;

	dd2 = dx2+dy2+dz2;
	dd = sqrt(dd2);
	/*
			  if(dd<0.01*aa)
		  {
				 {      A[2*(j*(Nl)+i)] += MM0[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=MM0[2*(jb*(4)+ib)+1] + eta;     }
		  }
	*/
if(j>i)
{
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  T2_rightlead[2*(j*(Nl)+i)] += Tx[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1]+=Tx[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  T2_rightlead[2*(j*(Nl)+i)] += Txy1[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1]+=Txy1[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  T2_rightlead[2*(j*(Nl)+i)] += Txy2[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1]+=Txy2[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      T2_rightlead[2*(j*(Nl)+i)] += Tz[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1]+=Tz[2*(jb*(4)+ib)+1];     }
	}

}

if(j<i)
{
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  T2_rightlead[2*(j*(Nl)+i)] += Txd[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1]+=Txd[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  T2_rightlead[2*(j*(Nl)+i)] += Txy1d[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1]+=Txy1d[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  T2_rightlead[2*(j*(Nl)+i)] += Txy2d[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1]+=Txy2d[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      T2_rightlead[2*(j*(Nl)+i)] += Tzd[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1]+=Tzd[2*(jb*(4)+ib)+1];     }
	}

}
}							 /* end of matrix A using for T2_rightlead  */


/************************/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_rightlead[2*(j*(Nl)+i)]=T2_rightlead[2*(i*(Nl)+j)];
		T1_rightlead[2*(j*(Nl)+i)+1]=-T2_rightlead[2*(i*(Nl)+j)+1];
	}
}


/*
	printf(" T1_rightlead:\n");
   for(i=0;i<(Nl+2);i++)
	 {for(j=0;j<(Nl+2);j++)
		{    printf("%.4f   %u %u  %.5f\t%.5f\n",e,i,j,T1_rightlead[2*(j*(Nl+2)+i)],T1_rightlead[2*(j*(Nl+2)+i)+1]);
		}
	 }
*//*
   printf("T2_rightlead:\n");
   for(i=0;i<(Nl+2);i++)
	 {for(j=0;j<ND;j++)
		{    printf("%.4f\t%u\t%u\t%.5f\t%.5f\n",e,i,j,T2_rightlead[2*(j*(Nl+2)+i)],T2_rightlead[2*(j*(Nl+2)+i)+1]);
		}
	 }
*/

T1_dr=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
T1_dl=(double*)malloc((Nl)*(Nl)*2*sizeof(double));

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_dr[2*(j*(Nl)+i)]=0.0;
		T1_dr[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		for(k=0;k<(Nl);k++)
		{
			T1_dr[2*(j*(Nl)+i)]+=T1_rightlead[2*(k*(Nl)+i)]*dr[2*(j*(Nl)+k)]-T1_rightlead[2*(k*(Nl)+i)+1]*dr[2*(j*(Nl)+k)+1];
			T1_dr[2*(j*(Nl)+i)+1]+=T1_rightlead[2*(k*(Nl)+i)]*dr[2*(j*(Nl)+k)+1]+T1_rightlead[2*(k*(Nl)+i)+1]*dr[2*(j*(Nl)+k)];
		}
	}
}


/*   printf(" T1_dr:\n");
  for(i=0;i<(Nl+2);i++)
	{for(j=0;j<(Nl+2);j++)
	   {    printf("%.4f   %u %u  %.5f\t%.5f\n",e,i,j,T1_rightlead[2*(j*(Nl+2)+i)],T1_rightlead[2*(j*(Nl+2)+i)+1]);
	   }
	}

  printf("dr:\n");
  for(i=0;i<(Nl+2);i++)
	{for(j=0;j<(Nl);j++)
	   {    printf("%.4f   %u %u  %.5f\t%.5f\n",e,i,j,T2_rightlead[2*(j*(Nl+2)+i)],T2_rightlead[2*(j*(Nl+2)+i)+1]);
	   }
	}

*/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		for(k=0;k<(Nl);k++)
		{
			sigma_rightlead[2*(j*(Nl)+i)]+=T1_dr[2*(k*(Nl)+i)]*T2_rightlead[2*(j*(Nl)+k)]-T1_dr[2*(k*(Nl)+i)+1]*T2_rightlead[2*(j*(Nl)+k)+1];
			sigma_rightlead[2*(j*(Nl)+i)+1]+=T1_dr[2*(k*(Nl)+i)]*T2_rightlead[2*(j*(Nl)+k)+1]+T1_dr[2*(k*(Nl)+i)+1]*T2_rightlead[2*(j*(Nl)+k)];
		}
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_leftlead[2*(j*(Nl)+i)]=0.0;
		T1_leftlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T2_leftlead[2*(j*(Nl)+i)]=0.0;
		T2_leftlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		sigma_leftlead[2*(j*(Nl)+i)]=0.0;
		sigma_leftlead[2*(j*(Nl)+i)+1]=0.0;
	}

}


/************ matrix A right using as matrix T2_leftlead **************/
for(i=0;i<Nl;i++)
	for(j=0;j<Nl;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j] - aa;	 /*  for right lead  */
	dy = yy[i] - yy[j];
	dz = zz[i] - zz[j];

	dx2 = dx*dx;
	dy2 = dy*dy;
	dz2 = dz*dz;

	dd2 = dx2+dy2+dz2;
	dd = sqrt(dd2);
	/*
			  if(dd<0.01*aa)
		  {
				 {      A[2*(j*(Nl)+i)] += MM0[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=MM0[2*(jb*(4)+ib)+1] + eta;     }
		  }
	*/
	
	if(j>i)
	{
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  T2_leftlead[2*(j*(Nl)+i)] += Tx[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1]+=Tx[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  T2_leftlead[2*(j*(Nl)+i)] += Txy1[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1]+=Txy1[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  T2_leftlead[2*(j*(Nl)+i)] += Txy2[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1]+=Txy2[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      T2_leftlead[2*(j*(Nl)+i)] += Tz[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1]+=Tz[2*(jb*(4)+ib)+1];     }
	}

}

if(j<i)
	{
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  T2_leftlead[2*(j*(Nl)+i)] += Txd[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1]+=Txd[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  T2_leftlead[2*(j*(Nl)+i)] += Txy1d[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1]+=Txy1d[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  T2_leftlead[2*(j*(Nl)+i)] += Txy2d[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1]+=Txy2d[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      T2_leftlead[2*(j*(Nl)+i)] += Tzd[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1]+=Tzd[2*(jb*(4)+ib)+1];     }
	}

}
}							 /* end of matrix A right using as T2_leftlead */


/************************/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_leftlead[2*(j*(Nl)+i)]=T2_leftlead[2*(i*(Nl)+j)];
		T1_leftlead[2*(j*(Nl)+i)+1]=-T2_leftlead[2*(i*(Nl)+j)+1];
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_dl[2*(j*(Nl)+i)]=0.0;
		T1_dl[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		for(k=0;k<(Nl);k++)
		{
			T1_dl[2*(j*(Nl)+i)]+=T1_leftlead[2*(k*(Nl)+i)]*dl[2*(j*(Nl)+k)]-T1_leftlead[2*(k*(Nl)+i)+1]*dl[2*(j*(Nl)+k)+1];
			T1_dl[2*(j*(Nl)+i)+1]+=T1_leftlead[2*(k*(Nl)+i)]*dl[2*(j*(Nl)+k)+1]+T1_leftlead[2*(k*(Nl)+i)+1]*dl[2*(j*(Nl)+k)];
		}

	}
}


/*   printf(" T1_leftlead:\n");
  for(i=v;i<ND;i++)
	{for(j=0;j<(Nl+2);j++)
	   {    printf("%.4f   %u %u  %.5f\t%.5f\n",e,i,j,T1_leftlead[2*(j*ND+i)],T1_leftlead[2*(j*ND+i)+1]);
	   }
	}
*//*
  printf("T2_leftlead:\n");
  for(i=0;i<(Nl+2);i++)
	{for(j=v;j<ND;j++)
	   {    printf("%.4f   %u %u  %.5f\t%.5f\n",e,i,j,T2_leftlead[2*(j*(Nl+2)+i)],T2_leftlead[2*(j*Nl+i)+1]);
	   }
	}
*/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		for(k=0;k<(Nl);k++)
		{
			sigma_leftlead[2*(j*(Nl)+i)]+=T1_dl[2*(k*(Nl)+i)]*T2_leftlead[2*(j*(Nl)+k)]-T1_dl[2*(k*(Nl)+i)+1]*T2_leftlead[2*(j*(Nl)+k)+1];
			sigma_leftlead[2*(j*(Nl)+i)+1]+=T1_dl[2*(k*(Nl)+i)]*T2_leftlead[2*(j*(Nl)+k)+1]+T1_dl[2*(k*(Nl)+i)+1]*T2_leftlead[2*(j*(Nl)+k)];
		}
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		green[2*(j*ND+i)]+=(-sigma_rightlead[2*(j*(Nl)+i)]);
		green[2*(j*ND+i)+1]+=(-sigma_rightlead[2*(j*(Nl)+i)+1]);
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		green[2*((j+ND-Nl)*ND+(i+ND-Nl))]+=(-sigma_leftlead[2*(j*(Nl)+i)]);
		green[2*((j+ND-Nl)*ND+(i+ND-Nl))+1]+=(-sigma_leftlead[2*(j*(Nl)+i)+1]);
	}
}


free(T1_dr);
free(T1_dl);

/*   printf(" sigma_right&leftlead:\n");
  for(i=0;i<ND;i++)
	{for(j=0;j<ND;j++)
	   {    printf("%.4f   %u %u  %.5f\t%.5f\t  %.5f\t%.5f\n",e,i,j,sigma_rightlead[2*(j*ND+i)],sigma_rightlead[2*(j*ND+i)+1],sigma_leftlead[2*(j*ND+i)],sigma_leftlead[2*(j*ND+i)+1]);
	   }
	}
*/

/*
   return 0;
  }
*/
