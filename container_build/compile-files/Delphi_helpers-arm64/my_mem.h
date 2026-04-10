
char  *alloc_1D_char     (int dim1) ;
int   *alloc_1D_int      (int dim1) ; 
float *alloc_1D_float    (int dim1) ; 

char  **alloc_2D_char    (int dim1, int dim2) ;
int   **alloc_2D_int     (int dim1, int dim2) ;
float **alloc_2D_float   (int dim1, int dim2) ;

char  ***alloc_3D_char   (int dim1, int dim2, int dim3) ;
int   ***alloc_3D_int    (int dim1, int dim2, int dim3) ;
int   ***alloc_3D_uint   (int dim1, unsigned short int *dim2, int dim3) ;
float ***alloc_3D_float  (int dim1, int dim2, int dim3) ;
float ***alloc_3D_ufloat (int dim1, unsigned short int *dim2, int dim3) ;

char  ****alloc_4D_char  (int dim1, int dim2, int dim3, int dim4) ; 
char  ****alloc_4D_uchar (int dim1, unsigned short int *dim2,
			  int dim3, int dim4) ;
int   ****alloc_4D_int   (int dim1, int dim2, int dim3, int dim4) ;
float ****alloc_4D_float (int dim1, int dim2, int dim3, int dim4) ;


float **mat_alloc_2       (int row, unsigned short int *col) ;


void clean_1D_char  (char *dim) ;
void clean_1D_int   (int *dim) ;
void clean_1D_float (float *dim) ;

void clean_2D_char  (char **mat, int dim2) ;
void clean_2D_int   (int **mat, int dim2) ;
void clean_2D_float (float **mat, int dim2) ;

void clean_3D_char  (char ***mat, int dim3, int dim2) ;
void clean_3D_int   (int ***mat, int dim3, int dim2) ;
void clean_3D_uint  (int ***mat, int dim3, unsigned short int *dim2) ;
void clean_3D_float (float ***mat, int dim3, int dim2) ;
void clean_3D_ufloat(float ***mat, int dim3, unsigned short int *dim2) ;

void clean_4D_char  (char ****mat, int dim4, int dim3, int dim2) ;
void clean_4D_uchar (char ****mat, int dim4, int dim3, 
		     unsigned short int *dim2) ;
void clean_4D_int   (int ****mat, int dim4, int dim3, int dim2) ;
void clean_4D_float (float ****mat, int dim4, int dim3, int dim2) ;



