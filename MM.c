for(i=0;i<Nl;i++)
for(j=0;j<Nl;j++)
	    {	
		M[2*(Nl*j+i)]=0.0;
		M[2*(Nl*j+i)+1]=0.0;
            }


for(i=0;i<Nl;i++)
{
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
	{
		M[2*(j*(Nl)+i)] += MM0[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += MM0[2*(jb*(4)+ib)+1]; } 
		
		if(zz[i]<0.01*aa && zz[i]>-0.01*aa)
			{ M[2*(j*(Nl)+i)] += zeeman_bottom[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += zeeman_bottom[2*(jb*(4)+ib)+1];}
		
		if(zz[i]<1.01*(Nz-1)*aa && zz[i]>0.99*(Nz-1)*aa)
		 	{M[2*(j*(Nl)+i)] += zeeman_top[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += zeeman_top[2*(jb*(4)+ib)+1];}
	
		

    if((dd<1.01*aa)&&(dd>0.99*aa))
	{
	if(j<i)
	{
	   
	    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  M[2*(j*(Nl)+i)] += Txd[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += Txd[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  M[2*(j*(Nl)+i)] += Tyd[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += Tyd[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  M[2*(j*(Nl)+i)] += Tzd[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += Tzd[2*(jb*(4)+ib)+1];  }
		   }
		   if(j>i)
		   {
		    
		if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  M[2*(j*(Nl)+i)] += Tx[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += Tx[2*(jb*(4)+ib)+1];  }
		
	    	if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  M[2*(j*(Nl)+i)] += Ty[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += Ty[2*(jb*(4)+ib)+1];  }
		
	    	if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  M[2*(j*(Nl)+i)] += Tz[2*(jb*(4)+ib)];    M[2*(j*(Nl)+i)+1] += Tz[2*(jb*(4)+ib)+1];  }
		
		}
	
	}

}							 /* end of matrix d  */
}

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
	
	if((dd<1.01*aa)&&(dd>0.99*aa))
	{
	
   	    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  Mr += Tx[2*(jb*(4)+ib)];    Mi += Tx[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  Mr += Ty[2*(jb*(4)+ib)];    Mi += Ty[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  Mr += Tz[2*(jb*(4)+ib)];    Mi += Tz[2*(jb*(4)+ib)+1];  }
	
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
	
	
	if((dd<1.01*aa)&&(dd>0.99*aa))
	{
	
	
	    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  Mr += Txd[2*(jb*(4)+ib)];    Mi += Txd[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  Mr += Tyd[2*(jb*(4)+ib)];   Mi += Tyd[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  Mr += Tzd[2*(jb*(4)+ib)];    Mi += Tzd[2*(jb*(4)+ib)+1];  }
		}
	

	Mrr = Mr*cos(theta) - Mi*sin(theta);	
	Mii = Mr*sin(theta) + Mi*cos(theta);
	
	M[2*(Nl*j+i)]+=Mrr;
	M[2*(Nl*j+i)+1]+=Mii;
}							 /* end of positive phase  */


/************************/

	theta = ky*aa;

/************ matrix with negative phase **************/
for(i=0;i<Nl;i++)
	for(j=0;j<Nl;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j];	 /*  for left lead  */
	dy = yy[i] - yy[j] + aa;
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
	
	if((dd<1.01*aa)&&(dd>0.99*aa))
	{
	
   	    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  Mr += Tx[2*(jb*(4)+ib)];    Mi += Tx[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  Mr += Ty[2*(jb*(4)+ib)];    Mi += Ty[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  Mr += Tz[2*(jb*(4)+ib)];    Mi += Tz[2*(jb*(4)+ib)+1];  }
	
}
	

	Mrr = Mr*cos(theta) - Mi*sin(theta);	
	Mii = Mr*sin(theta) + Mi*cos(theta);
	
	M[2*(Nl*j+i)]+=Mrr;
	M[2*(Nl*j+i)+1]+=Mii;
}							 /* end of negative phase  */


/************************/



	theta = -ky*aa;

/************ matrix with positive phase **************/
for(i=0;i<Nl;i++)
	for(j=0;j<Nl;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j];	 
	dy = yy[i] - yy[j] - aa;
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
	
	
	if((dd<1.01*aa)&&(dd>0.99*aa))
	{
	
	
	    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  Mr += Txd[2*(jb*(4)+ib)];    Mi += Txd[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  Mr += Tyd[2*(jb*(4)+ib)];   Mi += Tyd[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  Mr += Tzd[2*(jb*(4)+ib)];    Mi += Tzd[2*(jb*(4)+ib)+1];  }
		}
	

	Mrr = Mr*cos(theta) - Mi*sin(theta);	
	Mii = Mr*sin(theta) + Mi*cos(theta);
	
	M[2*(Nl*j+i)]+=Mrr;
	M[2*(Nl*j+i)+1]+=Mii;
}							 /* end of positive phase  */

