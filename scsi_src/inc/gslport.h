//Converts gsl matrix to double precision nr matrix and creates new pointer to pointer
#define GSL2NRDM1(d,m,p,s) 	double* d[m->size1+1];\
	{\
	int port_i;\
	for(port_i=s;port_i<=m->size1;port_i++) {\
		d[port_i+1] = m->data + port_i*m->size2 - 1;\
		}\
	}\
	double** p = d;
	
//Converts gsl matrix to double precision nr matrix and assigns to already existing pointer to pointer
#define GSL2NRDM2(d,m,p,s) 	double* d[m->size1+1];\
	{\
	int port_i;\
	for(port_i=s;port_i<=m->size1;port_i++) {\
		d[port_i+1] = m->data + port_i*m->size2 - 1;\
		}\
	}\
	p = d;

//Sets values of gsl matrix from an array of data 
#define GSLSET(data,mat,n,m,nn,mm) {\
int port_i, port_j;\
for (port_i=n; port_i<nn; port_i++)\
	for (port_j=m; port_j <mm; port_j++)\
		gsl_matrix_set (mat, port_i, port_j, data[port_i][port_j]);\
}	

//Sets values of gsl matrix from nr matrix
#define GSLSET2(nr,mat,n,m,nn,mm) {\
int port_i, port_j;\
for (port_i=n; port_i<nn; port_i++)\
	for (port_j=m; port_j <mm; port_j++)\
		gsl_matrix_set (mat, port_i, port_j, nr[port_i+1][port_j+1]);\
}

//Prints gsl matrix to standard output
#define GSLPRINT(x,n,m,nn,mm) {\
	int port_i, port_j;\
	for ( port_i=n; port_i<nn; port_i++ ) {\
      for ( port_j=m; port_j<mm; port_j++ ) {\
         printf("%g ", gsl_matrix_get(x, port_i, port_j));\
      }\
	  printf("\n");\
	}\
}
	
//Prints nr matrix to standard output	
#define NRPRINT(x,n,m,nn,mm) {\
int port_i, port_j;\
	for ( port_i=n; port_i<=nn; port_i++ ) {\
      for ( port_j=m; port_j<=mm; port_j++ ) {\
         printf("%g ", x[port_i][port_j]);\
      }\
	  printf("\n");\
	  }\
	 }
	 
//Converts gsl vector to double nr vector and creates pointer to it
#define GSL2NRDV1(v,p) double *p = v->data-1;
	
//Converts gsl vector to double nr vector 
#define GSL2NRDV2(v,p) p = v->data-1;


double ***gslport_tensor(int t1, int t2, int r1, int r2, int s1, int s2);

long unsigned int *gslport_lvector(int s, int e);

void gslport_matrixdump( FILE *outf, char *text, gsl_matrix *m, char *format);
