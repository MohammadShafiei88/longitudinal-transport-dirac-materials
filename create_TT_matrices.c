for(i=0;i<Nm;i++)
{
	for(j=0;j<Nm;j++)
	{
		Tx[2*(Nm*j+i)]    =0.0;
		Tx[2*(Nm*j+i)+1]  =0.0;
		Txy1[2*(Nm*j+i)]    =0.0;
		Txy1[2*(Nm*j+i)+1]  =0.0;
		Txy2[2*(Nm*j+i)]    =0.0;
		Txy2[2*(Nm*j+i)+1]  =0.0;
		Tz[2*(Nm*j+i)]    =0.0;
		Tz[2*(Nm*j+i)+1]  =0.0;
		Txd[2*(Nm*j+i)]    =0.0;
		Txd[2*(Nm*j+i)+1]  =0.0;
		Txy1d[2*(Nm*j+i)]    =0.0;
		Txy1d[2*(Nm*j+i)+1]  =0.0;
		Txy2d[2*(Nm*j+i)]    =0.0;
		Txy2d[2*(Nm*j+i)+1]  =0.0;
		Tzd[2*(Nm*j+i)]    =0.0;
		Tzd[2*(Nm*j+i)+1]  =0.0;
		MM0[2*(Nm*j+i)]    =0.0;
		MM0[2*(Nm*j+i)+1]  =0.0;
	}
}


Tx[0]  = CC2+BB2;
Tx[7]  = -0.5*AA2;
Tx[10] = CC2-BB2;
Tx[13] = -0.5*AA2;
Tx[19] = -0.5*AA2;
Tx[20] = CC2+BB2;
Tx[25] = -0.5*AA2;
Tx[30] = CC2-BB2;

Txy1[0]  = CC2+BB2;
Txy1[6]  = 0.25*sqrt(3.0)*AA2;
Txy1[7]  = -0.25*AA2;
Txy1[10] = CC2-BB2;
Txy1[12] = 0.25*sqrt(3.0)*AA2;
Txy1[13] = -0.25*AA2;
Txy1[18] = -0.25*sqrt(3.0)*AA2;
Txy1[19] = -0.25*AA2;
Txy1[20] = CC2+BB2;
Txy1[24] = -0.25*sqrt(3.0)*AA2;
Txy1[25] = -0.25*AA2;
Txy1[30] = CC2-BB2;

Txy2[0]  = CC2+BB2;
Txy2[6]  = -0.25*sqrt(3.0)*AA2;
Txy2[7]  =-0.25*AA2;
Txy2[10] = CC2-BB2;
Txy2[12] = -0.25*sqrt(3.0)*AA2;
Txy2[13] =-0.25*AA2;
Txy2[18] = 0.25*sqrt(3.0)*AA2;
Txy2[19] = -0.25*AA2;
Txy2[20] = CC2+BB2;
Txy2[24] = 0.25*sqrt(3.0)*AA2;
Txy2[25] = -0.25*AA2;
Txy2[30] = CC2-BB2;

Tz[0]  = CC1+BB1;
Tz[3]  =-0.5*AA1;
Tz[9]  =-0.5*AA1;
Tz[10] = CC1-BB1;
Tz[20] = CC1+BB1;
Tz[23] = 0.5*AA1;
Tz[29] = 0.5*AA1;
Tz[30] = CC1-BB1;

Txd[0]  = CC2+BB2;
Txd[7]  = 0.5*AA2;
Txd[10] = CC2-BB2;
Txd[13] = 0.5*AA2;
Txd[19] = 0.5*AA2;
Txd[20] = CC2+BB2;
Txd[25] = 0.5*AA2;
Txd[30] = CC2-BB2;

Txy1d[0]  = CC2+BB2;
Txy1d[6]  = -0.25*sqrt(3.0)*AA2;
Txy1d[7]  = 0.25*AA2;
Txy1d[10] = CC2-BB2;
Txy1d[12] = -0.25*sqrt(3.0)*AA2;
Txy1d[13] = 0.25*AA2;
Txy1d[18] = 0.25*sqrt(3.0)*AA2;
Txy1d[19] = 0.25*AA2;
Txy1d[20] = CC2+BB2;
Txy1d[24] = 0.25*sqrt(3.0)*AA2;
Txy1d[25] = 0.25*AA2;
Txy1d[30] = CC2-BB2;

Txy2d[0]  = CC2+BB2;
Txy2d[6]  = 0.25*sqrt(3.0)*AA2;
Txy2d[7]  =-0.25*AA2;
Txy2d[10] = CC2-BB2;
Txy2d[12] = 0.25*sqrt(3.0)*AA2;
Txy2d[13] = 0.25*AA2;
Txy2d[18] = -0.25*sqrt(3.0)*AA2;
Txy2d[19] = 0.25*AA2;
Txy2d[20] = CC2+BB2;
Txy2d[24] = -0.25*sqrt(3.0)*AA2;
Txy2d[25] = 0.25*AA2;
Txy2d[30] = CC2-BB2;

Tzd[0]  = CC1+BB1;
Tzd[3]  = 0.5*AA1;
Tzd[9]  = 0.5*AA1;
Tzd[10] = CC1-BB1;
Tzd[20] = CC1+BB1;
Tzd[23] =-0.5*AA1;
Tzd[29] =-0.5*AA1;
Tzd[30] = CC1-BB1;

MM0[0]  = M0;
MM0[10] =-M0;
MM0[20] = M0;
MM0[30] =-M0;
