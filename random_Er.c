//#include <string.h>
//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
//main()
//{
//  float    *Er;
//  int       i,j,N;
//   long       random_seed2;

//   random_seed2=0;
//   N=20;
//   Er=(float *)malloc(N*sizeof(float));

for(i=0;i<ND_1;i++)
	Er[i]=0.0;

srand48(random_seed);
for(i=0;i<ND_1;i++)
	Er[i]=drand48();

for(i=0;i<ND_1;i++)
	Er[i]=((Er[i]-0.5)*We);

//n_site=0;
//Er[248]=10000;
//Er[280]=10000;
//Er[249]=10000;
//Er[279]=10000;
//Er[247]=10000;
//Er[264]=10000;
//Er[255]=10000;

//for(i=0;i<N;i++)
//printf("Er[%d]=%f\n",i,Er[i]);

//free(Er);
//return 0;
//}
