/*#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
main()
{
	int i,j,N1,l;               j:for site l:number of sites
	double   t,e,eps;
	double   *d,*Dd,*A,*B;
	N1=7;
	t=1.0;
	eps=0.01;
	e=-4.0;
	l=N1;
	d=(double *)malloc(N1*N1*2*sizeof(double));
	Dd=(double *)malloc(N1*N1*2*sizeof(double));
	A=(double *)malloc(N1*N1*2*sizeof(double));
	B=(double *)malloc(N1*N1*2*sizeof(double));
  initialize_green function*/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		d[2*(j*(Nl)+i)]=0.0;
		d[2*(j*(Nl)+i)+1]=0.0;
		Dd[2*(j*(Nl)+i)]=0.0;
		Dd[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<Nl;i++)
	for(j=0;j<Nl;j++)
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
		{  d[2*(j*(Nl)+i)] += e;    d[2*(j*(Nl)+i)+1] += eta;  }
		
		

		if(dd<0.01*aa)
	{
		{      d[2*(j*(Nl)+i)] += -MM0[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -MM0[2*(jb*(4)+ib)+1];     }
	}
	
	if(j>i)
		{

	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  d[2*(j*(Nl)+i)] += -Tx[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Tx[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  d[2*(j*(Nl)+i)] += -Txy1[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Txy1[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  d[2*(j*(Nl)+i)] += -Txy2[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Txy2[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      d[2*(j*(Nl)+i)] += -Tz[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Tz[2*(jb*(4)+ib)+1];     }
	}
}

	if(j<i)
		{

	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  d[2*(j*(Nl)+i)] += -Txd[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Txd[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  d[2*(j*(Nl)+i)] += -Txy1d[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Txy1d[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  d[2*(j*(Nl)+i)] += -Txy2d[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Txy2d[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      d[2*(j*(Nl)+i)] += -Tzd[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Tzd[2*(jb*(4)+ib)+1];     }
	}
}
}							 /* end of matrix d  */


for(i=0;i<Nl;i++)
{
	for(j=0;j<Nl;j++)
	{
		Dd[2*(j*Nl+i)]=d[2*(j*Nl+i)];
		Dd[2*(j*Nl+i)+1]=d[2*(j*Nl+i)+1];
	}
}


for(i=0;i<Nl;i++)
{
	for(j=0;j<Nl;j++)
	{
		A[2*(j*Nl+i)]=0.0;
		A[2*(j*Nl+i)+1]=0.0;
		B[2*(j*Nl+i)]=0.0;
		B[2*(j*Nl+i)+1]=0.0;
		C[2*(j*Nl+i)]=0.0;
		C[2*(j*Nl+i)+1]=0.0;
	}
}


/************ matrix A **************/
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
				{  A[2*(j*(Nl)+i)] += Tx[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=Tx[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  A[2*(j*(Nl)+i)] += Txy1[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=Txy1[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  A[2*(j*(Nl)+i)] += Txy2[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=Txy2[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      A[2*(j*(Nl)+i)] += Tz[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=Tz[2*(jb*(4)+ib)+1];     }
	}

}

if(j<i)
	{
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  A[2*(j*(Nl)+i)] += Txd[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=Txd[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  A[2*(j*(Nl)+i)] += Txy1d[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=Txy1d[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  A[2*(j*(Nl)+i)] += Txy2d[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=Txy2d[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      A[2*(j*(Nl)+i)] += Tzd[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1]+=Tzd[2*(jb*(4)+ib)+1];     }
	}

}
}							 /* end of matrix A  */


/************************/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)

	{
		B[2*(j*(Nl)+i)]=A[2*(i*(Nl)+j)];
		B[2*(j*(Nl)+i)+1]=-A[2*(i*(Nl)+j)+1];
	}
}


/*

   for(i=0;i<(Nl+2);i++)
	{for(j=0;j<(Nl+2);j++)
	  printf("i=%d  , j=%d , dr[2*(j*(Nl+2)+i)]=%.6f ,   dr[2*(j*(Nl+2)+i)+1]=%.6f ,  A[2*(j*(Nl+2)+i)]=%.6f  ,  B[2*(j*(Nl+2)+i)]=%.6f\n",i,j,d[2*(j*(Nl+2)+i)],d[2*(j*(Nl+2)+i)+1],A[2*(j*(Nl+2)+i)],B[2*(j*(Nl+2)+i)]);
	 }

  */
