
for(i=0;i<Nm;i++)
{
	for(j=0;j<Nm;j++)
	{
		Tx[2*(Nm*j+i)]    =0.0;
		Tx[2*(Nm*j+i)+1]  =0.0;
		Ty[2*(Nm*j+i)]    =0.0;
		Ty[2*(Nm*j+i)+1]  =0.0;
		Tz[2*(Nm*j+i)]    =0.0;
		Tz[2*(Nm*j+i)+1]  =0.0;	

		zeeman_top[2*(Nm*j+i)]    =0.0;
		zeeman_top[2*(Nm*j+i)+1]  =0.0;
	        zeeman_bottom[2*(Nm*j+i)]    =0.0;
		zeeman_bottom[2*(Nm*j+i)+1]  =0.0;

	        Txd[2*(Nm*j+i)]    =0.0;
		Txd[2*(Nm*j+i)+1]  =0.0;
		Tyd[2*(Nm*j+i)]    =0.0;
		Tyd[2*(Nm*j+i)+1]  =0.0;
		Tzd[2*(Nm*j+i)]    =0.0;
		Tzd[2*(Nm*j+i)+1]  =0.0;
	
		
		MM0[2*(Nm*j+i)]    =0.0;
		MM0[2*(Nm*j+i)+1]  =0.0;
	}
}



Tx[0]  = BB;
Tx[25]  = -0.5*AA;
Tx[10] = -BB;
Tx[19] = -0.5*AA;
Tx[13] = -0.5*AA;
Tx[20] = BB;
Tx[7] = -0.5*AA;
Tx[30] = -BB;

Ty[0]  = BB;
Ty[24]  = -0.5*AA;
Ty[10] = -BB;
Ty[18] = -0.5*AA;
Ty[12] = 0.5*AA;
Ty[20] = BB;
Ty[6] = 0.5*AA;
Ty[30] = -BB;


Tz[0]  = BB;
Tz[9]  =-0.5*AA;
Tz[3]  = -0.5*AA;
Tz[10] = -BB;
Tz[20] = BB;
Tz[29] = 0.5*AA;
Tz[23] = 0.5*AA;
Tz[30] = -BB;



Txd[0]  = BB;
Txd[25]  = 0.5*AA;
Txd[10] = -BB;
Txd[19] = 0.5*AA;
Txd[13] = 0.5*AA;
Txd[20] = BB;
Txd[7] = 0.5*AA;
Txd[30] = -BB;

Tyd[0]  = BB;
Tyd[24]  = 0.5*AA;
Tyd[10] = -BB;
Tyd[18] = 0.5*AA;
Tyd[12] = -0.5*AA;
Tyd[20] = BB;
Tyd[6] = -0.5*AA;
Tyd[30] = -BB;


Tzd[0]  = BB;
Tzd[9]  = 0.5*AA;
Tzd[3]  = 0.5*AA;
Tzd[10] = -BB;
Tzd[20] = BB;
Tzd[29] = -0.5*AA;
Tzd[23] = -0.5*AA;
Tzd[30] = -BB;




MM0[0]  = MM-6*BB;
MM0[10] =-(MM-6*BB);
MM0[20] = MM-6*BB;
MM0[30] =-(MM-6*BB);

zeeman_top[0]  = delta_top ;
//zeeman_top[2]  = Mx_top;
//zeeman_top[3]  = My_top;
//zeeman_top[8]  = Mx_top;
//zeeman_top[9]  =-My_top ;
zeeman_top[10] = delta_top;
zeeman_top[20] = -delta_top;
//zeeman_top[22] = Mx_top;
//zeeman_top[23] = My_top ;
//zeeman_top[28] = Mx_top;
//zeeman_top[29] =-My_top;
zeeman_top[30] =-delta_top;
 
zeeman_bottom[0]  = delta_bottom ;
//zeeman_bottom[2]  = Mx_bottom;
//zeeman_bottom[3]  = My_bottom;
//zeeman_bottom[8]  = Mx_bottom;
//zeeman_bottom[9]  =-My_bottom ;
zeeman_bottom[10] = delta_bottom;
zeeman_bottom[20] = -delta_bottom;
//zeeman_bottom[22] = Mx_bottom;
//zeeman_bottom[23] = My_bottom ;
//zeeman_bottom[28] = Mx_bottom;
//zeeman_bottom[29] =-My_bottom;
zeeman_bottom[30] =-delta_bottom;

