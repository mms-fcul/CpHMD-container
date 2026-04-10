#include <stdlib.h>

/* Memory allocation and clean */

char *alloc_1D_char (int dim1) {
  char *vec;
  vec = calloc(dim1, sizeof(char)) ;
  return vec;
}

int *alloc_1D_int (int dim1) {
  int *vec;
  vec = calloc(dim1, sizeof(int)) ;
  return vec;
}

float *alloc_1D_float (int dim1) {
  float *vec;
  vec = calloc(dim1, sizeof(float)) ;
  return vec;
}


char **alloc_2D_char (int dim1, int dim2) {
  int i;
  char **mat;
  mat = calloc(dim1, sizeof(char *)) ;
  for (i = 0 ; i < dim1 ; i++)
    mat[i] = calloc(dim2, sizeof(char)) ;
  return mat;
}

int **alloc_2D_int (int dim1, int dim2) {
  int i;
  int **mat;
  mat = calloc(dim1, sizeof(int *)) ;
  for (i = 0 ; i < dim1 ; i++)
    mat[i] = calloc(dim2, sizeof(int)) ;
  return mat;
}

float **alloc_2D_float (int dim1, int dim2) {
  int i;
  float **mat;
  mat = calloc(dim1, sizeof(float *)) ;
  for (i = 0 ; i < dim1 ; i++)
    mat[i] = calloc(dim2, sizeof(float)) ;
  return mat;
}


char ***alloc_3D_char (int dim1, int dim2, int dim3) {
  int i, j ;
  char ***mat ;
  mat = calloc(dim1, sizeof(char **)) ;
  for (i = 0; i < dim1; i++) {
    mat[i] = calloc(dim2, sizeof(char *)) ;
    for (j = 0; j < dim2; j++) 
      mat[i][j] = calloc(dim3, sizeof(char )) ;
  }
  return mat ;
}

int ***alloc_3D_int (int dim1, int dim2, int dim3) {
  int i, j;
  int ***mat;
  mat = calloc(dim1, sizeof(int **)) ;
  for (i = 0 ; i < dim1 ; i++) {
    mat[i] = calloc(dim2, sizeof(int *)) ;
    for (j = 0; j < dim2; j++) 
      mat[i][j] = calloc(dim3, sizeof(int )) ;
  }
  return mat;
}

int ***alloc_3D_uint (int dim1, unsigned short int *dim2, int dim3) {
  int i, j;
  int ***mat;
  mat = calloc(dim1, sizeof(int **)) ;
  for (i = 0 ; i < dim1 ; i++) {
    mat[i] = calloc(*dim2, sizeof(int *)) ;
    for (j = 0; j < *dim2; j++) 
      mat[i][j] = calloc(dim3, sizeof(int )) ;
  }
  return mat;
}

float ***alloc_3D_float (int dim1, int dim2, int dim3) {
  int i, j;
  float ***mat;
  mat = calloc(dim1, sizeof(float **)) ;
  for (i = 0 ; i < dim1 ; i++) {
    mat[i] = calloc(dim2, sizeof(float *)) ;
    for (j = 0; j < dim2; j++) 
      mat[i][j] = calloc(dim3, sizeof(float )) ;
  }
  return mat;
}

float ***alloc_3D_ufloat (int dim1, unsigned short int *dim2, int dim3) {
  int i, j;
  float ***mat;
  mat = calloc(dim1, sizeof(float **)) ;
  for (i = 0 ; i < dim1 ; i++) {
    mat[i] = calloc(*dim2, sizeof(float *)) ;
    for (j = 0; j < *dim2; j++) 
      mat[i][j] = calloc(dim3, sizeof(float )) ;
  }
  return mat;
}


char ****alloc_4D_char (int dim1, int dim2, int dim3, int dim4) {
  int i, j, k ;
  char ****mat ;
  mat = calloc(dim1, sizeof(char ***)) ;
  for (i = 0; i < dim1; i++) {
    mat[i] = calloc(dim2, sizeof(char **)) ;
    for (j = 0; j < dim2; j++) {
      mat[i][j] = calloc(dim3, sizeof(char *)) ;
      for (k = 0; k < dim3; k++) 
	mat[i][j][k] = calloc(dim4, sizeof(char )) ;
    }
  }
  return mat ;
}

char ****alloc_4D_uchar (int dim1, unsigned short int *dim2,
			 int dim3, int dim4) {
  int i, j, k ;
  char ****mat ;
  mat = calloc(dim1, sizeof(char ***)) ;
  for (i = 0; i < dim1; i++) {
    mat[i] = calloc(*dim2, sizeof(char **)) ;
    for (j = 0; j < *dim2; j++) {
      mat[i][j] = calloc(dim3, sizeof(char *)) ;
      for (k = 0; k < dim3; k++) 
	mat[i][j][k] = calloc(dim4, sizeof(char )) ;
    }
  }
  return mat ;
}

int ****alloc_4D_int (int dim1, int dim2, int dim3, int dim4) {
  int i, j, k ;
  int ****mat ;
  mat = calloc(dim1, sizeof(int ***)) ;
  for (i = 0; i < dim1; i++) {
    mat[i] = calloc(dim2, sizeof(int **)) ;
    for (j = 0; j < dim2; j++) {
      mat[i][j] = calloc(dim3, sizeof(int *)) ;
      for (k = 0; k < dim3; k++) 
	mat[i][j][k] = calloc(dim4, sizeof(int )) ;
    }
  }
  return mat ;
}

float ****alloc_4D_float (int dim1, int dim2, int dim3, int dim4) {
  int i, j, k ;
  float ****mat ;
  mat = calloc(dim1, sizeof(float ***)) ;
  for (i = 0; i < dim1; i++) {
    mat[i] = calloc(dim2, sizeof(float **)) ;
    for (j = 0; j < dim2; j++) {
      mat[i][j] = calloc(dim3, sizeof(float *)) ;
      for (k = 0; k < dim3; k++) 
	mat[i][j][k] = calloc(dim4, sizeof(float )) ;
    }
  }
  return mat ;
}


/*------------------------------------------*/


float **mat_alloc_2 (int row, unsigned short int *col) {
  int i ;
  float **mat ;
  mat = calloc(*col, sizeof(float *)) ;
  for (i = 0; i < *col; i++) 
    mat[i] = calloc(row, sizeof(float)) ;
  return mat ;
}


/*------------------------------------------*/


void clean_1D_char (char *vec) {
  free (vec);
}

void clean_1D_int (int *vec) {
  free (vec);
}

void clean_1D_float (float *vec) {
  free (vec);
}


void clean_2D_char (char **mat, int dim2) {
  int i;
  for (i = 0 ; i < dim2 ; i++) free(mat[i]) ;
  free(mat) ;
}

void clean_2D_int (int **mat, int dim2) {
  int i;
  for (i = 0 ; i < dim2 ; i++) free(mat[i]) ;
  free(mat) ;
}

void clean_2D_float (float **mat, int dim2) {
  int i;
  for (i = 0 ; i < dim2 ; i++) free(mat[i]) ;
  free(mat) ;
}


void clean_3D_char (char ***mat, int dim3, int dim2) {
  int i, j ;
  for (i = 0; i < dim3; i++) {
    free(mat[i]);
      for (j = 0; j < dim2; j++) 
	free(mat[i][j]) ;
  }
  free(mat) ;
}

void clean_3D_int (int ***mat, int dim3, int dim2) {
  int i, j;
  for (i = 0 ; i < dim3; i++) {
    free(mat[i]) ;
    for (j = 0; j < dim2; j++) 
      free(mat[i][j]) ;
  }
  free(mat) ;
}

void clean_3D_uint (int ***mat, int dim3, unsigned short int *dim2) {
  int i, j;
  for (i = 0 ; i < dim3; i++) {
    free(mat[i]) ;
    for (j = 0; j < *dim2; j++) 
      free(mat[i][j]) ;
  }
  free(mat) ;
}

void clean_3D_float (float ***mat, int dim3, int dim2) {
  int i, j;
  for (i = 0 ; i < dim3; i++) {
    free(mat[i]) ;
    for (j = 0; j < dim2; j++) 
      free(mat[i][j]) ;
  }
  free(mat) ;
}

void clean_3D_ufloat (float ***mat, int dim3, unsigned short int *dim2) {
  int i, j;
  for (i = 0 ; i < dim3; i++) {
    free(mat[i]) ;
    for (j = 0; j < *dim2; j++) 
      free(mat[i][j]) ;
  }
  free(mat) ;
}

void clean_4D_char (char ****mat, int dim4, int dim3, int dim2) {
  int i, j, k ;
  for (i = 0; i < dim4; i++) {
    free(mat[i]);
    for (j = 0; j < dim3; j++) {
      free(mat[i][j]) ;
      for (k = 0; k < dim2; k++) 
	free(mat[i][j][k]) ;
    }
  }
  free(mat) ;
}

void clean_4D_uchar (char ****mat, int dim4, int dim3, 
		     unsigned short int *dim2) {
  int i, j, k ;
  for (i = 0; i < dim4; i++) {
    free(mat[i]);
    for (j = 0; j < dim3; j++) {
      free(mat[i][j]) ;
      for (k = 0; k < *dim2; k++) 
	free(mat[i][j][k]) ;
    }
  }
  free(mat) ;
}

void clean_4D_int (int ****mat, int dim4, int dim3, int dim2) {
  int i, j, k ;
  for (i = 0; i < dim4; i++) {
    free(mat[i]);
    for (j = 0; j < dim3; j++) {
      free(mat[i][j]) ;
      for (k = 0; k < dim2; k++) 
	free(mat[i][j][k]) ;
    }
  }
  free(mat) ;
}

void clean_4D_float (float ****mat, int dim4, int dim3, int dim2) {
  int i, j, k ;
  for (i = 0; i < dim4; i++) {
    free(mat[i]);
    for (j = 0; j < dim3; j++) {
      free(mat[i][j]) ;
      for (k = 0; k < dim2; k++) 
	free(mat[i][j][k]) ;
    }
  }
  free(mat) ;
}










