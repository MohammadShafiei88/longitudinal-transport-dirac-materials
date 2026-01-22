for(i=0;i<Nl;i++)
for(j=0;j<Nl;j++)
	    {	
		M[2*(Nl*j+i)]=0.0;
		M[2*(Nl*j+i)+1]=0.0;
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
		{  M[2*(j*(Nl)+i)] += 0.0;    M[2*(j*(Nl)+i)+1] += 0.0;  }

		if(dd<0.01*aa)
	{
		{      M[2*(j*(Nl)+i)] += -MM0[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += -MM0[2*(jb*(4)+ib)+1];     }
	}

	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  M[2*(j*(Nl)+i)] += -Tx[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += -Tx[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  M[2*(j*(Nl)+i)] += -Txy1[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += -Txy1[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  M[2*(j*(Nl)+i)] += -Txy2[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += -Txy2[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      M[2*(j*(Nl)+i)] += -Tz[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += -Tz[2*(jb*(4)+ib)+1];     }
	}

}							 /* end of matrix d  */


	theta = kx*aa;

/************ matrix with negative phase **************/
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

	Mr = 0.0;
	Mi = 0.0;
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  Mr += Tx[2*(jb*(4)+ib)];    Mi+=Tx[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  Mr += Txy1[2*(jb*(4)+ib)];    Mi+=Txy1[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  Mr += Txy2[2*(jb*(4)+ib)];    Mi+=Txy2[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      Mr += Tz[2*(jb*(4)+ib)];    Mi+=Tz[2*(jb*(4)+ib)+1];     }
	}
	

	Mrr = Mr*cos(theta) - Mi*sin(theta);	
	Mii = Mr*sin(theta) + Mi*cos(theta);
	
	M[2*(Nl*j+i)]+=Mrr;
	M[2*(Nl*j+i)+1]+=Mii;
}							 /* end of negative phase  */


/************************/



	theta = -kx*aa;

/************ matrix with positive phase **************/
for(i=0;i<Nl;i++)
	for(j=0;j<Nl;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j] - aa;	 /*  for left lead  */
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

	Mr = 0.0;
	Mi = 0.0;
	if((dd<1.01*aa)&&(0.99*aa<dd))
	{
		if(dz2<0.01*cc*cc)		 /*  x-bond or xy1-bond or xy2-bond  */
		{
			if(dy2<0.01*aa*aa)	 /*  x-bond  */
				{  Mr += Tx[2*(jb*(4)+ib)];    Mi+=Tx[2*(jb*(4)+ib)+1];  }
								 /*  xy1-bond  */
					else if(dy<-0.8*aa)
					{  Mr += Txy1[2*(jb*(4)+ib)];    Mi+=Txy1[2*(jb*(4)+ib)+1];  }
								 /*  xy2-bond  */
						else if(dy>0.8*aa)
						{  Mr += Txy2[2*(jb*(4)+ib)];    Mi+=Txy2[2*(jb*(4)+ib)+1];  }
		}
	}

	if((dd<1.01*cc)&&(0.99*cc<dd))
	{
								 /*  z-bond  */
		if((dx2<0.01*aa*aa) && (dy2<0.01*aa*aa))
			{      Mr += Tz[2*(jb*(4)+ib)];    Mi+=Tz[2*(jb*(4)+ib)+1];     }
	}
	

	Mrr = Mr*cos(theta) - Mi*sin(theta);	
	Mii = Mr*sin(theta) + Mi*cos(theta);
	
	M[2*(Nl*j+i)]+=Mrr;
	M[2*(Nl*j+i)+1]+=Mii;
}							 /* end of positive phase  */


/************************/

