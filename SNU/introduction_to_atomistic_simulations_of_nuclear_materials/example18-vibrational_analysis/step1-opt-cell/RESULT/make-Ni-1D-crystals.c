#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define PROJ "Ni_1D"
#define MS1 58.7	// mass, Fe

#define NA 1		// number of atoms in unit cell, Fe
#define LC 2.78   // lattice constant (A)
#define SA 1		// supercell size for a-axis
#define SB 1		//   , b-axis
#define SC 1		//   , c-axis


int main(void)
{
	int i,j,k,a,b,c,n,m;
	int s;
	char fname1[300];
        int p,q;
	double lat0[3][3]={
		{1.0, 0.0, 0.0},
		{0.0, 10.0, 0.0},
		{0.0, 0.0, 10.0} 
	};
	double lat1[3][3];
	int sc[3]={SA,SB,SC};
	// position of iron atoms
        double pos1[NA][3]={
        	{0.0, 0.0, 0.0},
        };

	double tpos1[3],tpos2[3];
	FILE *fw[2],*fr[2];

	
	///////	
	// obtain the cell vectors of the supercell
	//////
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			lat1[i][j]=lat0[i][j]*LC*sc[i];


	/////
	// lammps output
	/////
	for(s=0;s<1;s++)
	{
		sprintf(fname1,"data.%s",PROJ);
		if( (fw[s]=fopen(fname1,"w"))==NULL )
		{
			printf("error: cannot open the file at %d\n",s);
			exit(1);	
		}
       	 	fprintf(fw[s],"#  #\n");  
		fprintf(fw[s],"\n");
		fprintf(fw[s],"%d    atoms\n",SA*SB*SC*NA);
		fprintf(fw[s],"%d  atom types\n",1);
		fprintf(fw[s],"\n");

		fprintf(fw[s],"%.10lf  %.10lf   xlo xhi\n",0.0,lat1[0][0]);
		fprintf(fw[s],"%.10lf  %.10lf   ylo yhi\n",0.0,lat1[1][1]);
		fprintf(fw[s],"%.10lf  %.10lf   zlo zhi\n",0.0,lat1[2][2]);
		fprintf(fw[s],"\n");

		fprintf(fw[s],"Masses\n");
		fprintf(fw[s],"\n");
		fprintf(fw[s],"1  %.6lf\n",MS1);	
		fprintf(fw[s],"\n");

		fprintf(fw[s],"Atoms\n");
		fprintf(fw[s],"\n");


	        n=0;
       	 	for(a=0;a<SA;a++)
                	for(b=0;b<SB;b++)
                        	for(c=0;c<SC;c++)
					for(m=0;m<NA;m++)
					{
						tpos1[0]=(pos1[m][0]+a)/SA;
						tpos1[1]=(pos1[m][1]+b)/SB;
						tpos1[2]=(pos1[m][2]+c)/SC;
						for(i=0;i<3;i++)
						{
							tpos2[i]=0.0;
							for(j=0;j<3;j++)	tpos2[i]+=(tpos1[j]*lat1[j][i]);
						}
						n++;
                                       		fprintf(fw[s],"%-7d 1  %13.10lf  %13.10lf  %13.10lf\n",
                                        		n,tpos2[0],tpos2[1],tpos2[2]);
					}
		fclose(fw[s]);
	}


	return 0;

}
