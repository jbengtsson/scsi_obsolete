double ***gslport_tensor(int t1, int t2, int r1, int r2, int s1, int s2) {

	double ***array = (double***)malloc(t2*sizeof(double**));
    int i, j;

	for (i = t1; i <=t2; i++) {
		array[i] = (double**)malloc(r2*sizeof(double*));
		
		for (j = r1; j < r2; j++) {
			array[i][j] = (double*)malloc(s2*sizeof(double));
		}
	}
	
	return array;

}


long unsigned int *gslport_lvector(int s, int e) {
	
	long unsigned int *data =
	  (long unsigned int*)malloc((e+1)*sizeof(long unsigned int));
	
	return data;
}


void gslport_matrixdump( FILE *outf, char *text, gsl_matrix *m, char *format)
{

   //fprintf( outf, "%s", text);
   printf("%s\n", text);
   //gsl_matrix_fprintf(outf, m, format);
   GSLPRINT(m,0,0,m->size1,m->size2);
   printf("\n", text);
   //fprintf( outf, "\n");
}
