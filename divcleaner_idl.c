/* idl_bicgstab.c
 *
 * IDL wrapper over bicgstab.c
 *
 */
 
#include "divcleaner_idl.h"

void poisson(int argc, void* argv[]) {
	/* Input arguments */
	int bX, bY, bZ;
	double thr;
	int maxiter;

	/* Keyword arguments */
	double hx=1., hy=1., hz=1.;

	/* Local variables */
	int a1, a2, a3, N, nghost=1;
	void* stream = NULL;

	double *b  = NULL, \
	       *x0 = NULL, \
	       *S  = NULL, \
	       *x  = NULL, \
	       *g  = NULL, \
        *b_arr = NULL, \
        *d_arr = NULL, \
	       *x_out = NULL, \
	       *x0_arr = NULL;

 int *dims  = NULL;
	
	Lambda L = (Lambda)&stencil3DL;

	/* Parse input arguments */
 if (argc != 9) {
  fprintf(stderr, \
    "divcleaner_poisson: wrong number of arguments.\n");
  return;
 }
 x_out   =   (double*) argv[0];
 x0_arr  =   (double*) argv[1];
 dims    =   (int*)    argv[2];
 bX      = ( (int*)    argv[3])[0];
 bY      = ( (int*)    argv[3])[1];
 bZ      = ( (int*)    argv[3])[2];
 thr     =  *(double*) argv[4];
 maxiter =  *(int*)    argv[5];
 /* Optional arguments */
 b_arr   =   (double*) argv[6];
 d_arr   =   (double*) argv[7];
 if (argv[8] != NULL) {
  hx     = ( (double*) argv[8])[0];
  hy     = ( (double*) argv[8])[1];
  hz     = ( (double*) argv[8])[2];
 }

 if (d_arr == NULL && \
			(bX == BC_DIRICHLET || bY == BC_DIRICHLET \
			 || bZ == BC_DIRICHLET)) {
  fprintf(stderr, \
    "divcleaner_poisson: With Dirichlet boundary \
    conditions the keyword argument `boundary` \
    must be provided.\n");
  return;
 }

	/* Get dimensions of the computational box */
	a1 = (int)dims[0]+2*nghost;
	a2 = (int)dims[1]+2*nghost;
	a3 = (int)dims[2]+2*nghost;
	N = a1*a2*a3;

	/* Create the stencil */
	S = stencil2OLaplacian3D(hx, hy, hz);

	/* Create the stream for the laplacian lambda operator */
	if (d_arr != NULL)
		g = saveGhostCells3D(d_arr, a1, a2, a3, nghost);
	stream = streamForLaplacian3DL(a1, a2, a3, nghost, \
			bX, bY, bZ, S, g);

	/* Allocate C-arrays */
	b = arrmalloc(N);
	x = arrmalloc(N);

	/* Set b */
	if (argv[6] != NULL) {
 	nullifyGhostFrame3D(b, a1, a2, a3, nghost);
		addGhostFrame3D(b, b_arr, a1, a2, a3, nghost);
	}
	else nullify(b, N);

	/* Update b taking into account constant boundary conditions */
	nonLinearPart(x, L, a1, a2, a3, nghost, stream);
	add(b, b, -1., x, N);

	/* Update stream */
	if (bX == BC_DIRICHLET) bX = BC_ZERO;
	if (bY == BC_DIRICHLET) bY = BC_ZERO;
	if (bZ == BC_DIRICHLET) bZ = BC_ZERO;
	free(g);
	free(stream);
	stream = streamForLaplacian3DL(a1, a2, a3, nghost, \
			bX, bY, bZ, S, NULL);

	/* Set x0 */
 x0 = arrmalloc(N);
	nullifyGhostFrame3D(x0, a1, a2, a3, nghost);
	addGhostFrame3D(x0, x0_arr, a1, a2, a3, nghost);

	/* Solve the Poisson equation */
	solveL(x, L, stream, b, x0, thr, maxiter, N);

	/* Remove ghost cells */
	removeGhostFrame3D(x_out, x, a1, a2, a3, nghost);

	/* Free memory */	
	free(b);
 free(x0);
	free(x);
	free(stream);
	free(S);
 
 return;
}

void divergence(int argc, void* argv[]) {
	/* Keyword arguments */
	double hx=1., hy=1., hz=1.;

	/* Local variables */
	int a1, a2, a3, i, j, k;
	double *Bb1, *Bb2, *Bb3, *div_out;
	int *dims = NULL;

	/* Parse input arguments */
 if (argc != 6) {
  fprintf(stderr, \
    "divcleaner_div: wrong number of arguments.\n");
  return;
 }
 div_out =   (double*) argv[0];
 Bb1     =   (double*) argv[1];
 Bb2     =   (double*) argv[2];
 Bb3     =   (double*) argv[3];
 dims    =   (int*)    argv[4];
 if (argv[5] != NULL) {
  hx     = ( (double*) argv[5])[0];
  hy     = ( (double*) argv[5])[1];
  hz     = ( (double*) argv[5])[2];
 }

	a1 = dims[0];
	a2 = dims[1];
	a3 = dims[2];
 
	/* Compute divergence */
	for (k=a3-1; k>=0; k--) {
	  for (j=a2-1; j>=0; j--) {
	    for (i=a1-1; i>=0; i--) {
	      div_out[(k*a2+j)*a1+i] = - Bb1[(k*a2+j)*(a1+1)+i]/hx \
	                               + Bb1[(k*a2+j)*(a1+1)+i+1]/hx \
	                               - Bb2[(k*(a2+1)+j)*a1+i]/hy \
	                               + Bb2[(k*(a2+1)+j+1)*a1+i]/hy \
	                               - Bb3[(k*a2+j)*a1+i]/hz \
	                               + Bb3[((k+1)*a2+j)*a1+i]/hz;
	    }
	  }
	}

	return;
}

void fixTopBottomBoundary(int argc, void* argv[]) {
	/* Keyword arguments */
	double hx=1., hy=1., hz=1.;

	/* Local variables */
	int a1, a2, a3, i, j;
	double *Bb1, *Bb2, *Bb3, *Bb3_out;
	int *dims = NULL;

	/* Parse input arguments */
 if (argc != 6) {
  fprintf(stderr, \
    "divcleaner_fixTopBottomBoundary: wrong number of arguments.\n");
  return;
 }
 Bb3_out =   (double*) argv[0];
 Bb1     =   (double*) argv[1];
 Bb2     =   (double*) argv[2];
 Bb3     =   (double*) argv[3];
 dims    =   (int*)    argv[4];
 if (argv[5] != NULL) {
  hx     = ( (double*) argv[5])[0];
  hy     = ( (double*) argv[5])[1];
  hz     = ( (double*) argv[5])[2];
 }

	a1 = dims[0];
	a2 = dims[1];
	a3 = dims[2];

	/* Fix top and bottom boundaries */
	arrcpy(Bb3_out, Bb3, a1*a2*(a3+1));
	for (j=a2-1; j>=0; j--) {
	  for (i=a1-1; i>=0; i--) {
	    Bb3_out[j*a1+i] =         - hy*hz*Bb1[j*(a1+1)+i] \
	                              + hy*hz*Bb1[j*(a1+1)+i+1] \
	                              - hz*hx*Bb2[j*a1+i] \
	                              + hz*hx*Bb2[(j+1)*a1+i] \
	                              + hx*hy*Bb3[(a2+j)*a1+i];
	    Bb3_out[j*a1+i] /= hx*hy;
	    Bb3_out[(a3*a2+j)*a1+i] = + hy*hz*Bb1[((a3-1)*a2+j)*(a1+1)+i] \
	                              - hy*hz*Bb1[((a3-1)*a2+j)*(a1+1)+i+1] \
	                              + hz*hx*Bb2[((a3-1)*(a2+1)+j)*a1+i] \
	                              - hz*hx*Bb2[((a3-1)*(a2+1)+j+1)*a1+i] \
	                              + hx*hy*Bb3[((a3-1)*a2+j)*a1+i];
	    Bb3_out[(a3*a2+j)*a1+i] /= hx*hy;
	  }
	}

	return;
}

void phiToBbc(int argc, void* argv[]) {
	/* Keyword arguments */
	double hx=1., hy=1., hz=1.;

	/* Local variables */
	int a1, a2, a3, i, j, k;
	double *phi, *Bb1, *Bb2, *Bb3;
	int *dims = NULL;

	/* Parse input arguments */
 if (argc != 6) {
  fprintf(stderr, \
    "divcleaner_phiToBbc: wrong number of arguments.\n");
  return;
 }
 Bb1     =   (double*) argv[0];
 Bb2     =   (double*) argv[1];
 Bb3     =   (double*) argv[2];
 phi     =   (double*) argv[3];
 dims    =   (int*)    argv[4];
 if (argv[5] != NULL) {
  hx     = ( (double*) argv[5])[0];
  hy     = ( (double*) argv[5])[1];
  hz     = ( (double*) argv[5])[2];
 }

	a1 = (int)dims[0];
	a2 = (int)dims[1];
	a3 = (int)dims[2];

	/* Compute Bb */
	for (k=a3-1; k>=0; k--) {
	  for (j=a2-1; j>=0; j--) {
	    for (i=a1-1; i>=1; i--) {
	      Bb1[(k*a2+j)*(a1+1)+i] = \
	           (phi[(k*a2+j)*a1+i] - phi[(k*a2+j)*a1+i-1])/hx;
	    }
	    Bb1[(k*a2+j)*(a1+1)] = \
	         (phi[(k*a2+j)*a1] - phi[(k*a2+j)*a1+a1-1])/hx;
	    Bb1[(k*a2+j)*(a1+1)+a1] = Bb1[(k*a2+j)*(a1+1)];
	  }
	}
	for (k=a3-1; k>=0; k--) {
	  for (i=a1-1; i>=0; i--) {
	    for (j=a2-1; j>=1; j--) {
	      Bb2[(k*(a2+1)+j)*a1+i] = \
	           (phi[(k*a2+j)*a1+i] - phi[(k*a2+j-1)*a1+i])/hy;
	    }
	    Bb2[(k*(a2+1))*a1+i] = \
	         (phi[(k*a2)*a1+i] - phi[(k*a2+a2-1)*a1+i])/hy;
	    Bb2[(k*(a2+1)+a2)*a1+i] = Bb2[(k*(a2+1))*a1+i];
	  }
	}
	for (j=a2-1; j>=0; j--) {
	  for (i=a1-1; i>=0; i--) {
	    for (k=a3-1; k>=1; k--) {
	      Bb3[(k*a2+j)*a1+i] = \
	           (phi[(k*a2+j)*a1+i] - phi[((k-1)*a2+j)*a1+i])/hz;
	    }
	    Bb3[j*a1+i] = \
	         (phi[j*a1+i] - phi[((a3-1)*a2+j)*a1+i])/hz;
	    Bb3[(a3*a2+j)*a1+i] = Bb3[(k*a2+j)*a1+i];
	  }
	}

	return;
}

void killdiv(int argc, void* argv[]) {
	/* Input arguments */
	double *Bb1, *Bb2, *Bb3;

	/* Keyword arguments */
	int bX=BC_PERIODIC, bY=BC_PERIODIC, bZ=BC_ZERO;
 int *dims;
	double *x0;
 double thr=1.e-5;
	int maxiter=20;

	/* Output variables */
	double *Bb1_out, *Bb2_out, *Bb3_out;

	/* Local variables */
 void *argv_in[9];
 double *div_out, *phi, *Bb1_irrot, *Bb2_irrot, *Bb3_irrot, *Bb3_tmp;
	int x0alloc=0;
 int i, bc[3];

	/* Parse input arguments */
 if (argc != 13) {
  fprintf(stderr, \
    "divcleaner_killdiv: wrong number of arguments.\n");
  return;
 }
 Bb1_out =   (double*) argv[0];
 Bb2_out =   (double*) argv[1];
 Bb3_out =   (double*) argv[2];
 Bb1     =   (double*) argv[3];
 Bb2     =   (double*) argv[4];
 Bb3     =   (double*) argv[5];
 dims    =   (int*)    argv[6];
 /* Optional arguments */
 if (argv[8] != NULL) {
  bX     = ( (int*)    argv[8])[0];
  bY     = ( (int*)    argv[8])[1];
  bZ     = ( (int*)    argv[8])[2];
 }
 x0      =   (double*) argv[10];
 if (argv[11] != NULL)
  thr    =  *(double*) argv[11];
 if (argv[12] != NULL)
  maxiter=  *(int*)    argv[12];

	/* Compute div(B) */
 div_out = (double*) malloc(dims[0]*dims[1]*dims[2]*sizeof(double));
 argv_in[0] = (void*)div_out;
 argv_in[1] = (void*)Bb1;
 argv_in[2] = (void*)Bb2;
 argv_in[3] = (void*)Bb3;
 argv_in[4] = (void*)dims;
 argv_in[5] = argv[7];
 divergence(6, argv_in);

	/* Build x0 if it was not provided */
 if (x0 == NULL) {
  x0 = (double*) malloc(dims[0]*dims[1]*dims[2]*sizeof(double));
  for (i=0;i<dims[0]*dims[1]*dims[2];i++) x0[i] = 1.0;
  x0alloc = 1;
 }

	/* Solve the poisson equation */
 phi = (double*) malloc(dims[0]*dims[1]*dims[2]*sizeof(double));
 bc[0] = bX;
 bc[1] = bY;
 bc[2] = bZ;
 argv_in[0] = (void*)phi;
 argv_in[1] = (void*)x0;
 argv_in[2] = (void*)dims;
 argv_in[3] = (void*)bc;
 argv_in[4] = (void*)&thr;
 argv_in[5] = (void*)&maxiter;
 argv_in[6] = (void*)div_out;
 argv_in[7] = argv[9];
 argv_in[8] = argv[7];
 poisson(9, argv_in);
 free(div_out);

	/* Compute irrotational field component */
 Bb1_irrot = (double*) malloc((dims[0]+1)*dims[1]*dims[2]*sizeof(double));
 Bb2_irrot = (double*) malloc(dims[0]*(dims[1]+1)*dims[2]*sizeof(double));
 Bb3_irrot = (double*) malloc(dims[0]*dims[1]*(dims[2]+1)*sizeof(double));
 argv_in[0] = (void*)Bb1_irrot;
 argv_in[1] = (void*)Bb2_irrot;
 argv_in[2] = (void*)Bb3_irrot;
 argv_in[3] = (void*)phi;
 argv_in[4] = argv[6];
 argv_in[5] = argv[7];
 phiToBbc(6, argv_in);
 free(phi);

	/* Fix top and bottom boundaries */
 Bb3_tmp = (double*) malloc(dims[0]*dims[1]*(dims[2]+1)*sizeof(double));
 for (i=0;i<(dims[0]+1)*dims[1]*dims[2];i++) Bb1_out[i] = Bb1[i]-Bb1_irrot[i];
 for (i=0;i<dims[0]*(dims[1]+1)*dims[2];i++) Bb2_out[i] = Bb2[i]-Bb2_irrot[i];
 for (i=0;i<dims[0]*dims[1]*(dims[2]+1);i++) Bb3_tmp[i] = Bb3[i]-Bb3_irrot[i];
 argv_in[0] = (void*)Bb3_out;
 argv_in[1] = (void*)Bb1_out;
 argv_in[2] = (void*)Bb2_out;
 argv_in[3] = (void*)Bb3_tmp;
 argv_in[4] = argv[6];
 argv_in[5] = argv[7];
 fixTopBottomBoundary(6, argv_in);
 free(Bb3_tmp);
 
 if(x0alloc == 1) free(x0);

	return;
}
