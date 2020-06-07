#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#define N0 8   		// the number of atoms in the supercell for force evaluation 
#define N1 1 		// the number of atoms that will be independently treated in vibration analysis (condition: N1<=N0, N0%N1==0)
#define LCON 2.752	// lattice constant (A)
#define NK 11 		// the number of sampling points in Brillouin zone
#define TOL1 1.0E-13	// if fabs(d)<TOL1 is achieved, d is set to be 0.0.
#define CUTD 10.0	// cutoff distance in interaction

//#define _DEBUG


int main (void)
{
	if(N0%N1)	{  printf("error: N0 should be dividable with N1\n");  exit(1);  }
	if(N0<N1)	{  printf("error: N1 should be smaller than (or equal to) N0\n");  exit(1);  }
	int i,j,k,t;
	const double bvec=2.0*M_PI/(LCON*N1);		// reciprocal lattice vector
	double db;
	if(NK!=1)	db=1.0*bvec/((NK-1));	// interval between sampling points (over Brillouin zone)
	else	db=0.0;
	double kvec; 		// to store k-vector
	double dif_phase;	// to store phase difference
        double dd[N0]={		// data for dynamical matrix
		0.063811641,
		-0.032496518,
		0.000578577,
		1.21202E-05,
		0.0,
		1.21202E-05,
		0.000578577,
		-0.032496518
	};
	double data0[N0][N0];	// data for dynamical matrix
	double pos0[N0];	// equilibrium position (force=0) of atoms in the unit cell
	for(i=0;i<N0;i++)	pos0[i]=i*LCON;	
	double pos0sc[3][N0];	// position of atoms in the supercell ([0][]: left cell, [1][]: original cell, [2][0]: right cell
	for(i=0;i<3;i++)        for(j=0;j<N0;j++) 	pos0sc[i][j]=(i-1)*LCON*N0+pos0[j];

  	double dynm1r[N1][N1],dynm1i[N1][N1];	// dynamical matrix components real part (dynm1r) imaginary part (dynm1i)
	double dynm1_lin[N1*N1*2];	// to convert dynm1r and dynm1i to gsl-complex: [i*2] is real par and [i*2+1] is imaginary part of i-th complex
	const double conv=(1.6E-19/(1.0E-20))/(0.001/6.02E23);	// unti conversion factor from [eV, A, mass-number] to [J, m, kg]
        static double w2_res[NK][N1];  	// [i][j]: to store the w^2 value (eigenvalue) for j-th mode of i-th k-point
        static double evec_res_r[NK][N1][N1],evec_res_i[NK][N1][N1]; 	// [i][j][k]: to store eigenvector of k-th atom for j-th mode at i-th k-point


	// to avoid artificial force due to the periodic boundary, this condition is needed
	if(CUTD>=(N0*LCON/2.0))	{  printf("error: CUTD should be smaller than a half of cell size\n");  exit(1);  }

	// constructing data for dynamical matrix
	printf("prepartaion-for-raw-dasta-for-dynamical-matrix\n");
	for(i=0;i<N0;i++)
       	{
         	for(j=0;j<N0;j++)
               	{
                   	data0[i][j]=dd[(j-i+N0)%N0];
			printf("%lf  ",data0[i][j]);
		}
		printf("\n");
	}
	for(i=0;i<3;i++)	for(j=0;j<N0;j++)	{  pos0sc[i][j]=(i-1)*LCON*N0+pos0[j];  printf("%d\t%lf\n",i*N0+j,pos0sc[i][j]);  }

	// sampling over NK k-points
	for(t=0;t<NK;t++)
	{
		if(NK!=1)	kvec=-0.5*bvec+t*db;	// k-vector for Brullian zone sampling
		else	kvec=0.0;

		// zeroing the dynamical matrix elements
		for(i=0;i<N1;i++)	
			for(j=0;j<N1;j++)	
			{
				dynm1r[i][j]=0.0;
				dynm1i[i][j]=0.0;
			}

		// constructing dynamical matrix with given k-vector
		for(i=0;i<N1;i++)
		{
			for(j=0;j<N0;j++)
			{
				for(k=0;k<3;k++)
				{
					if(sqrt(pow(pos0sc[k][j]-pos0[i],2.0))<CUTD )
					{
						dif_phase=kvec*(pos0sc[k][N1*(j/N1)]-pos0[0]); 
						//dif_phase=kvec*(pos0sc[k][j]-pos0[i]);
						dynm1r[i][j%N1]+=(data0[i][j]*cos(-dif_phase));
						dynm1i[i][j%N1]+=(data0[i][j]*sin(-dif_phase));
					}
				}
			}
		}
               
#ifdef _DEBUG
             	printf("***real-part****\n");
              	for(i=0;i<N1;i++)
               	{
                   	for(j=0;j<N1;j++)	printf("%.3e ",dynm1r[i][j]);
                       	printf("\n");
              	}
              	printf("***imaginary-part****\n");
              	for(i=0;i<N1;i++)
              	{
                 	for(j=0;j<N1;j++)	printf("%.3e ",dynm1i[i][j]);
                     	printf("\n");
                }
#endif

		// checking whether the matrix is symmetric
		for(i=0;i<N1;i++)
		{
			for(j=i+1;j<N1;j++)
			{
				if(dynm1r[i][j]!=dynm1r[j][i])
				{
					// if matrix is not symmetric but almost symmetric 
					if( fabs(dynm1r[i][j])<TOL1 && fabs(dynm1r[j][i])<TOL1)
					{
                                                printf("warning: symmetrized (zeroing) at t=%d i=%d j=%d : %e %e\n",
                                                        t,i,j,dynm1r[i][j],dynm1r[j][i]);
						sleep(0.5);
						dynm1r[i][j]=0.0;
						dynm1r[j][i]=0.0;
					}
					else if( dynm1r[i][j]!=0.0 && dynm1r[j][i]!=0.0 && fabs(fabs(dynm1r[j][i]/dynm1r[i][j])-1.0)<0.01)
					{
						printf("warning: symmetrized at t=%d i=%d j=%d : %e %e\n",
							t,i,j,dynm1r[i][j],dynm1r[j][i]);
						sleep(0.5);
						dynm1r[i][j]=0.5*(dynm1r[i][j]+dynm1r[j][i]);
						dynm1r[j][i]=dynm1r[i][j];
					}	
					else	// if matrix is not symmetric, exit with an error message
					{
						printf("error: not symmetric-r, %e\t%e at t=%d i=%d j=%d\n",dynm1r[i][j],dynm1r[j][i],t,i,j);
						exit(1);
					}
				}
                                if(dynm1i[i][j]!=(-dynm1i[j][i]))
                                {

					// if matrix is not symmetric but almost symmetric 
                                        if( fabs(dynm1i[i][j])<TOL1 && fabs(dynm1i[j][i])<TOL1)
                                        {
						printf("warning: symmetrized at t=%d i=%d j=%d : %e %e\n",
                                                        t,i,j,dynm1i[i][j],dynm1i[j][i]);
						sleep(0.5);
                                                dynm1i[i][j]=0.0;
                                                dynm1i[j][i]=0.0;
                                        }
                                        else if( dynm1i[i][j]!=0.0 && dynm1i[j][i]!=0.0 && fabs(fabs(dynm1i[j][i]/dynm1i[i][j])-1.0)<0.01)
                                        {
                                                printf("warning: symmetrized at t=%d i=%d j=%d : %e %e\n",
                                                        t,i,j,dynm1i[i][j],dynm1i[j][i]);
						sleep(0.5);
                                                dynm1i[i][j]=0.5*(dynm1i[i][j]-dynm1i[j][i]);
                                                dynm1i[j][i]=-dynm1i[i][j];
                                        }
                                        else	// if matrix is not symmetric, exit with an error message
                                        {
                                                printf("error: not symmetric-i, %e\t%e at t=%d i=%d j=%d\n",dynm1i[i][j],dynm1i[j][i],t,i,j);
                                                exit(1);
                                        }
                                }
			}
		}				
	
		// to copy the data in dynm1r and dynm1i into dynm1_lin
		for(i=0;i<N1;i++)	
			for(j=0;j<N1;j++)	
			{
				dynm1_lin[i*N1*2+j*2]=dynm1r[i][j];
				dynm1_lin[i*N1*2+j*2+1]=dynm1i[i][j];
			}	

#ifdef _DEBUG
		printf("***real-part****\n");
		for(i=0;i<N1;i++)
		{
			for(j=0;j<N1;j++)	printf("%.3e ",dynm1r[i][j]);
			printf("\n");
		}
             	printf("***imaginary-part****\n");
              	for(i=0;i<N1;i++)
             	{
                   	for(j=0;j<N1;j++)	printf("%.3e ",dynm1i[i][j]);
                   	printf("\n");
		}
#endif

		// solving eigen value problem of Hermitian matrix using gsl
  		gsl_matrix_complex_view m = gsl_matrix_complex_view_array (dynm1_lin, N1, N1);
  		gsl_vector *eval = gsl_vector_alloc (N1);
  		gsl_matrix_complex *evec = gsl_matrix_complex_alloc (N1, N1);

		gsl_eigen_hermv_workspace * w = gsl_eigen_hermv_alloc (N1);
		gsl_eigen_hermv (&m.matrix, eval, evec, w);
		gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
 		gsl_eigen_hermv_free (w);
		gsl_complex zz;

		// store eigenvalues and eigenvectors in w2_res and "evec_res_r and evec_res_i", respectively.
 		for (i=0;i<N1;i++)
		{
			gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec, i);
			w2_res[t][i]=conv*(gsl_vector_get(eval, i));
			for(j=0;j<N1;j++)
			{
				zz=gsl_vector_complex_get(&evec_i.vector,j);
				evec_res_r[t][i][j]=GSL_REAL(zz);
				evec_res_i[t][i][j]=GSL_IMAG(zz);
			}
		}
 		gsl_vector_free (eval);
  		gsl_matrix_complex_free (evec);

/*
		// This is to test cosine-only version.
		// Cosine only version is only correct at some k-points, such as gamma-point.
		// So, in general, cosine-only version is wrong. We need to use complex number version
                gsl_matrix_view m_cos
                = gsl_matrix_view_array (dynm1_lin, N1, N1);
                gsl_vector *eval_cos = gsl_vector_alloc (N1);
                gsl_matrix *evec_cos = gsl_matrix_alloc (N1, N1);

                gsl_eigen_symmv_workspace * w_cos =  gsl_eigen_symmv_alloc (N1);
                gsl_eigen_symmv (&m_cos.matrix, eval_cos, evec_cos, w_cos);
                gsl_eigen_symmv_sort (eval_cos, evec_cos, GSL_EIGEN_SORT_ABS_ASC);
		gsl_eigen_symmv_free (w_cos);


                for (i=0;i<N1;i++)
                {
                        gsl_vector_view evec_cos_i
                                = gsl_matrix_column (evec_cos, i);
                        w2_res[t][i]=conv*(gsl_vector_get(eval_cos, i));
                        for(j=0;j<N1;j++)       evec_res_r[t][i][j]=gsl_vector_get(&evec_cos_i.vector,j);
                }
                gsl_vector_free (eval_cos);
                gsl_matrix_free (evec_cos);
*/

	}

	// print out normal mode information
	printf("k-vector(1/A);  frequency(cm-1), eigenvector\n");
	for(i=0;i<NK;i++)
	{
		printf("mode-%d\n",i);
		kvec=-0.5*bvec+i*db;
		for(j=0;j<N1;j++)		
		{
			// if imaginary frequency (w2<0), minus sign is put on frequency (w)
			// the unit of frequency is converted to cm-1
			if(w2_res[i][j]<0.0)	printf("%10.5lf\t%10.2lf\n",kvec,-33.35641*sqrt(fabs(w2_res[i][j]))/1.0E12);
			else printf("%10.5lf\t%10.2lf\n",kvec,33.35641*sqrt(w2_res[i][j])/1.0E12);
			// print out eigenvectors 
			// since we have an conjugate part, the imaginary part will always disapear.
			for(k=0;k<N1;k++)	printf("    %7.3lf + %.3lf i \n",evec_res_r[i][j][k],evec_res_i[i][j][k]);
			printf("\n");
		}
		printf("\n");
	}

	// print out  only frequency information
	printf("frequency-summary\n");
        printf("k-vector(1/A)\tfrequency(cm-1)....\n");
	
        for(i=0;i<NK;i++)
        {
		kvec=-0.5*bvec+i*db;
		printf("%lf\t",kvec);
		for(j=0;j<N1;j++)
		{
			// if imaginary frequency (w2<0), minus sign is put on frequency (w)
			// the unit of frequency is converted to cm-1
                        if(w2_res[i][j]<0.0) 	printf("%8.2lf\t",-33.35641*sqrt(fabs(w2_res[i][j]))/1.0E12);
                        else 	printf("%8.2lf\t",33.35641*sqrt(w2_res[i][j])/1.0E12);
                }
		printf("\n");

	}

  	return 0;
}
