/* bicgstab.c
 *
 * Implementation of the Bi-CGSTAB (biconjugate gradient stabilized) method
 * for the solution of the Ax = b linear system.
 *
 * Conventions for all functions:
 *  - Output arrays are passed as (already allocated) first argument pointers
 *  - Output scalars are returned by the functions
 *  - Input arrays follow in the usual order, as pointers
 *  - In some cases memory for output arrays is allocated inside functions;
 *    then a pointer to the array is returned
 *  - Sizes of the arrays are at the end and are not pointers
 *  - Ordering of arrays is column-major in the linear algebra routines
 *  - Ordering of arrays is row-major in stencils and arrays with
 *    spatial coordinates
 *  - Arrays are (in general) uppercase and scalars are lowercase; vectors
 *    depend on context
 *
 */

#include "bicgstab.h"

/* solve
 *
 * Solve the system Ax = b, with intial guess x0, iterating till
 * threshold; shape(A) = (n,n)
 *
 * A has to be known explicitly
 *
 */
void solve(double* x, const double* A, const double* b, const double* x0, \
		double thr, int maxiter, int n) {
	void* stream;
	Lambda L = (Lambda)&explicitL;

	stream = malloc(sizeof(int) + n*sizeof(double));
	((int*)stream)[0] = n;
	*(const double**)&((int*)stream)[1] = A;
	solveL(x, L, stream, b, x0, thr, maxiter, n);
}


/* solveL
 *
 * Solve the system Ax = b, with intial guess x0, iterating till
 * threshold; shape(A) = (n,n)
 *
 * Ax is computed through the Lambda-operator L
 *
 */
void solveL(double* x, Lambda L, const void* stream, const double* b, \
		const double* x0, double thr, int maxiter, int n) {
	double rho0, rho1, alpha, beta, omega, err=0.;
	double *r, *r0hat, *v, *p, *s, *t;
	int i;

	arrcpy(x, x0, n);
	r = arrmalloc(n);
	arrcpy(r, b, n);
	v = arrmalloc(n);
	arrcpy(v, x0, n);
	(*L)(v, x0, n, stream);
	add(r, r, -1., v, n);
	r0hat = arrmalloc(n);
	arrcpy(r0hat, r, n);
	rho1 = 1.;
	alpha = 1.;
	omega = 1.;
	nullify(v, n);
	p = arrmalloc(n);
	nullify(p, n);
	s = arrmalloc(n);
	t = arrmalloc(n);
	arrcpy(t, x0, n);

	i = 0;
	do {
		rho0 = rho1;
		dot(&rho1, r0hat, r, 1, n, 1);
		beta = rho1/rho0 * alpha/omega;
		add(p, r, beta, p, n);
		add(p, p, -beta*omega, v, n);
		(*L)(v, p, n, stream);
		dot(&alpha, r0hat, v, 1, n, 1);
		alpha = rho1/alpha;
		add(s, r, -alpha, v, n);
		err = normmax(s, n);
		if (err <= thr) {
			add(x, x, alpha, p, n);
			fprintf(stdout, "Finished at iteration %d with err = %e.\n", i, err);
			break;
		}
		(*L)(t, s, n, stream);
		dot(&beta, t, t, 1, n, 1);
		dot(&omega, t, s, 1, n, 1);
		omega /= beta;
		add(x, x, alpha, p, n);
		add(x, x, omega, s, n);
		add(r, s, -omega, t, n);
		i += 1;
		fprintf(stdout, "Iteration %d (err = %e)...\n", i, err);
	} while (err > thr && (i < maxiter || maxiter < 0));
	if (maxiter >= 0 && i >= maxiter) {
		fprintf(stdout, "Reached maximum nomber of iterations with err = %e.\n", err);
	}

	free(r);
	free(v);
	free(r0hat);
	free(p);
	free(s);
	free(t);

	return;
}

/* Building blocks of the Bi-CGSTAB method */

/* arrmalloc
 *
 * Allocate memory for an array of size n
 *
 */
double* arrmalloc(int n) {
	double *p = (double*)malloc(n*sizeof(double));

	if (p == NULL) {
		printf("arrmalloc: allocation failed.\n");
	}

	return p;
}

/* nullify
 *
 * Set to zero all values of array
 *
 */
void nullify(double* A, int n) {
	int i;

	for (i=n-1; i>=0; i--) {
		A[i] = 0.;
	}

	return;
}

/* oneify
 *
 * Set to one all values of array
 *
 */
void oneify(double* A, int n) {
	int i;

	for (i=n-1; i>=0; i--) {
		A[i] = 1.;
	}

	return;
}

/* arrcpy
 *
 * Copy array
 *
 */
void arrcpy(double* A, const double* B, int n) {
	memcpy((void*)A, (void*)B, (size_t)(n*sizeof(double)));

	return;
}

/* dot
 *
 * Dot product of two 2D-arrays; shape(A) = (a,b) and shape(B) = (b,c)
 *
 */
void dot(double* C, const double* A, const double* B, int a, int b, int c) {
	int i, j, k;

	for (k=c-1 ; k>=0 ; k--) {
		for (i=a-1 ; i>=0 ; i--) {
			C[i*c+k] = 0.;
			for (j=b-1 ; j>=0 ; j--) {
				C[i*c+k] += A[i*b+j]*B[j*c+k];
			}
		}
	}

	return;
}

/* add
 *
 * Add two arrays:
 * 
 *   C = A + r B
 *
 */
void add(double* C, const double* A, double r, const double* B, int n) {
	int i;

	for (i=n-1 ; i>=0 ; i--) {
		C[i] = A[i] + r*B[i];
	}

	return;
}

/* norm2
 *
 * Returns the Euclidean norm of a vector
 *
 */
double norm2(const double* v, int n) {
	int i;
	double norm = 0.;

	for (i=n-1; i>=0; i--) {
		norm += pow(v[i], 2.0);
	}

	return sqrt(norm);
}

/* normmax
 *
 * Returns the maximum norm of a vector
 *
 */
double normmax(const double* v, int n) {
	int i;
	double norm = -INFINITY;

	for (i=n-1; i>=0; i--) {
		norm = max(norm, fabs(v[i]));
	}
	return norm;
}

/* explicitL
 *
 * Lambda operator with evaluations from explicit A matrix
 *
 * The matrix A is encoded in the stream variable, which consists of:
 *
 *   stream = (int n, double** p),
 *
 * with shape(A) = (n, n) and p is the address of A
 *
 */
void explicitL(double* b, const double* x, int n, const void* stream) {
	int m = SVAR(stream, int, 0);
	double *A = SVAR(stream, double*, sizeof(int));

	if (n != m) {
		fprintf(stderr, "explicitL: invalid vector size: sz=%d, shape(Lambda operator)=(%d, %d).\n", n, m, m);
		exit(-1);
	}

	dot(b, A, x, n, n, 1);

	return;
}

/* stencil3DL
 *
 * Lambda operator for a general stencil
 *
 * The (2*n+1)^3 stencil is encoded in stream. First integer is boudary
 * conditions. Next three integers in stream are the shape of the spatial
 * array, next integer is n (hence the size of the stencil is (2*n+1)^3),
 * then follows a pointer to the stencil and finally a pointer to the ghost
 * cells providing optional constant boundary conditions (can be NULL if
 * boundary conditions are not constant along any axis).
 *
 */
void stencil3DL(double* b, double* x, int n, const void* stream) {
	int bc = SVAR(stream,int,0);
	int bX = (bc >> BC_X) & BC;
	int bY = (bc >> BC_Y) & BC;
	int bZ = (bc >> BC_Z) & BC;
	int a1 = SVAR(stream, int, sizeof(int));
	int a2 = SVAR(stream, int, 2*sizeof(int));
	int a3 = SVAR(stream, int, 3*sizeof(int));
	int nghost = SVAR(stream, int, 4*sizeof(int));
	double *S = SVAR(stream, double*, 5*sizeof(int));
	double *g = SVAR(stream, double*, 5*sizeof(int)+sizeof(double*));
	double *g_save = NULL;

	if (n != a1*a2*a3) {
		fprintf(stderr, "stencilRegularBCL: invalid vector size: sz=%d, shape(Lambda operator)=(%d, %d).\n", n, a1*a2*a3, a1*a2*a3);
		exit(-1);
	}

	g_save = saveGhostCells3D(x, a1, a2, a3, nghost);
	restoreGhostCells3D(x, g, a1, a2, a3, nghost);
	nullifyGhostEdges3D(x, a1, a2, a3, nghost);
	setBoundaryConditions3D(x, a1, a2, a3, nghost, bX, bY, bZ);
	applyStencil3D(b, S, x, a1, a2, a3, nghost);
	popGhostCells3D(x, g_save, a1, a2, a3, nghost);
	nullifyGhostFrame3D(b, a1, a2, a3, nghost);

	return;
}

/* nonLinearPart
 *
 * Extract the non linear part from the given L operator
 *
 * Some Lambda operators (e.g. stencil operators) with specific boundary
 * conditions might contain non-linear contributions. This means that it
 * might happen that L(x) = Ax + b instead of just L(x) = Ax.
 *
 * This routine extracts the non-linear part b.
 *
 */
void nonLinearPart(double* b, Lambda L, \
		int a1, int a2, int a3, \
		int n, const void* stream) {
	int N = a1*a2*a3;
	double *x0 = arrmalloc(N);

	nullify(x0, N);
	oneifyGhostFrame3D(x0, a1, a2, a3, n);
	nullifyGhostEdges3D(x0, a1, a2, a3, n);
	(*L)(b, x0, N, stream);
	free(x0);
	nullifyGhostFrame3D(b, a1, a2, a3, n);

	return;
}

/* 
 * Other functions unused in the present implementation of the
 * Bi-CGSTAB method
 */

/* dotadd
 *
 * Adds the dot product of two arrays to a third array:
 *
 *   C <-- C + r < A ; B >
 *
 */
void dotadd(double* C, double r, const double* A, const double* B, \
		int a, int b, int c) {
	int i, j, k;

	for (k=c-1 ; k>=0 ; k--) {
		for (i=a-1 ; i>=0 ; i--) {
			for (j=b-1 ; j>=0 ; j--) {
				C[i*c+k] += r*A[i*b+j]*B[j*c+k];
			}
		}
	}

	return;
}

/* mult
 *
 * Multiply an array by a scalar
 *
 */
void mult(double* C, double s, const double* A, int n) {
	int i;

	for (i=n-1 ; i>=0 ; i--) {
		C[i] = s*A[i];
	}

	return;
}

/*
 * Useful functions in the construction of Lambda operators:
 * Stencil management
 *
 * /!\ Row-major ordering of arrays and stencils /!\
 *
 */

/* streammalloc
 *
 * Allocate memory for a stream of length n (bytes)
 *
 */
double* streammalloc(int n) {
	return (void*)malloc(n);
}

/* addGhostFrame3D
 *
 * Add ghost cells to an array
 *
 */
void addGhostFrame3D(double* C, const double* A, \
		int a1, int a2, int a3, int n) {
	int i, j, k;

	for (k=a3-n-1 ; k>=n ; k--) {
	  for (j=a2-n-1 ; j>=n ; j--) {
	    for (i=a1-n-1 ; i>=n ; i--) {
	      C[(k*a2+j)*a1+i] = A[((k-n)*(a2-2*n)+j-n)*(a1-2*n)+i-n];
	    }
	  }
	}
	return;
}

/* removeGhostFrame3D
 *
 * Remove ghost cells from an array
 *
 */
void removeGhostFrame3D(double* C, const double* A, \
		int a1, int a2, int a3, int n) {
	int i, j, k;

	for (k=a3-n-1 ; k>=n ; k--) {
	  for (j=a2-n-1 ; j>=n ; j--) {
	    for (i=a1-n-1 ; i>=n ; i--) {
	      C[((k-n)*(a2-2*n)+j-n)*(a1-2*n)+i-n] = A[(k*a2+j)*a1+i];
	    }
	  }
	}
	return;
}

/* nullifyGhostFrame3D
 *
 * Set to zero ghost cells of 3D array
 *
 */
void nullifyGhostFrame3D(double* C, int a1, int a2, int a3, int n) {
	int i, j, k;

	if (n == 0) return;
	for (k=a3-1 ; k>=0 ; k--) {
	  for (j=a2-1 ; j>=0 ; j--) {
	    for (i=a1-1 ; i>=0 ; i--) {
	      C[(k*a2+j)*a1+i] = 0.;
	      if (i == a1-n && (j>=n && j<a2-n) && (k>=n && k<a3-n)) i = n;
	    }
	  }
	}

	return;
}

/* oneifyGhostFrame3D
 *
 * Set to one ghost cells of 3D array
 *
 */
void oneifyGhostFrame3D(double* C, int a1, int a2, int a3, int n) {
	int i, j, k;

	if (n == 0) return;
	for (k=a3-1 ; k>=0 ; k--) {
	  for (j=a2-1 ; j>=0 ; j--) {
	    for (i=a1-1 ; i>=0 ; i--) {
	      C[(k*a2+j)*a1+i] = 1.;
	      if (i == a1-n && (j>=n && j<a2-n) && (k>=n && k<a3-n)) i = n;
	    }
	  }
	}

	return;
}

/* nullifyGhostEdges3D
 *
 * Set to zero ghost edges of 3D array
 *
 */
void nullifyGhostEdges3D(double* C, int a1, int a2, int a3, int n) {
	int i, j, k;

	if (n == 0) return;
	for (k=a3-1 ; k>=0 ; k--) {
	  for (j=a2-1 ; j>=0 ; j--) {
	    for (i=a1-1 ; i>=0 ; i--) {
	      C[(k*a2+j)*a1+i] = 0.;
	      if (i == a1-n && ((j>=n && j<a2-n) || (k>=n && k<a3-n))) i = n;
	    }
	    if (j == a2-n && (k>=n && k<a3-n)) j = n;
	  }
	}

	return;
}

/* nullifyGhostFace3D
 *
 * Set to zero specified ghost face of 3D array
 *
 */
void nullifyGhostFace3D(double* C, int a1, int a2, int a3, int n, int f) {
	int i, j, k, i0=n, j0=n, k0=n, i1=a1-n-1, j1=a2-n-1, k1=a3-n-1;

	if (n == 0) return;
	switch(f) {
	  case 0 :
	    i0 = 0;
	    i1 = n-1;
	    break;
	  case 1 :
	    j0 = 0;
	    j1 = n-1;
	    break;
	  case 2 :
	    k0 = 0;
	    k1 = n-1;
	    break;
	  default :
	    return;
	}
	for (k=k1 ; k>=k0 ; k--) {
	  for (j=j1 ; j>=j0 ; j--) {
	    for (i=i1 ; i>=i0 ; i--) {
	      C[(k*a2+j)*a1+i] = 0.;
	    }
	  }
	}
	switch(f) {
	  case 0 :
	    i0 = a1-n;
	    i1 = a1-1;
	    break;
	  case 1 :
	    j0 = a2-n;
	    j1 = a2-1;
	    break;
	  case 2 :
	    k0 = a3-n;
	    k1 = a3-1;
	    break;
	  default :
	    return;
	}
	for (k=k1 ; k>=k0 ; k--) {
	  for (j=j1 ; j>=j0 ; j--) {
	    for (i=i1 ; i>=i0 ; i--) {
	      C[(k*a2+j)*a1+i] = 0.;
	    }
	  }
	}

	return;
}

/* saveGhostCells3D
 *
 * Returns an array with ghost cells
 *
 */
double* saveGhostCells3D(const double* A, int a1, int a2, int a3, int n) {
	int i, j, k, l=0;
	int N = a1*a2*a3-(a1-2*n)*(a2-2*n)*(a3-2*n);
	double *G = NULL;

	if (n == 0) return G;
	G = arrmalloc(N);
	for (k=a3-1 ; k>=0 ; k--) {
	  for (j=a2-1 ; j>=0 ; j--) {
	    for (i=a1-1 ; i>=0 ; i--) {
	      G[l] = A[(k*a2+j)*a1+i];
	      l++;
	      if (i == a1-n && (j>=n && j<a2-n) && (k>=n && k<a3-n)) i = n;
	    }
	  }
	}

	if (l != N) {
		fprintf(stderr, "saveGhostCells3D: assertion error.");
		exit(-1);
	}

	return G;
}

/* popGhostCells3D
 *
 * Restore ghost cells to array A and free memory, or just free memory
 * if C = NULL
 *
 */
void popGhostCells3D(double* C, double* G, int a1, int a2, int a3, int n) {
	int i, j, k, l=0;
	int N = a1*a2*a3-(a1-2*n)*(a2-2*n)*(a3-2*n);

	if (G == NULL) return;
	if (C != NULL) {
	  for (k=a3-1 ; k>=0 ; k--) {
	    for (j=a2-1 ; j>=0 ; j--) {
	      for (i=a1-1 ; i>=0 ; i--) {
	        C[(k*a2+j)*a1+i] = G[l];
	        l++;
	        if (i == a1-n && (j>=n && j<a2-n) && (k>=n && k<a3-n)) i = n;
	      }
	    }
	  }
	  if (l != N) {
	    fprintf(stderr, "popGhostCells3D: assertion error.");
	    exit(-1);
	  }
	}
	free(G);
	G = NULL;

	return;
}

/* restoreGhostCells3D
 *
 * Restore ghost cells to array A and keep G (don't free memory)
 *
 */
void restoreGhostCells3D(double* C, const double* G, \
		int a1, int a2, int a3, int n) {
	int i, j, k, l=0;
	int N = a1*a2*a3-(a1-2*n)*(a2-2*n)*(a3-2*n);

	if (G == NULL) {
		return;
	}
	for (k=a3-1 ; k>=0 ; k--) {
	  for (j=a2-1 ; j>=0 ; j--) {
	    for (i=a1-1 ; i>=0 ; i--) {
	      C[(k*a2+j)*a1+i] = G[l];
	      l++;
	      if (i == a1-n && (j>=n && j<a2-n) && (k>=n && k<a3-n)) i = n;
	    }
	  }
	}
	if (l != N) {
	  fprintf(stderr, "restoreGhostCells3D: assertion error.");
	  exit(-1);
	}

	return;
}

/* setBoundaryConditions3D()
 *
 * Set the ghost cells of array to satisfy given boundary conditions
 *
 */
void setBoundaryConditions3D(double* C, int a1, int a2, int a3, int n, \
		int bX, int bY, int bZ) {
	int i, j, k;

	nullifyGhostEdges3D(C, a1, a2, a3, n);
	switch(bX) {
	  case BC_ZERO :
	    for (k=a3-n-1 ; k>=n ; k--) {
	      for (j=a2-n-1 ; j>=n ; j--) {
	        for (i=n-1 ; i>=0 ; i--) {
	          C[(k*a2+j)*a1+i] = 0.;
	          C[(k*a2+j+1)*a1-i-1] = 0.;
	        }
	      }
	    }
	    break;
	  case BC_DIRICHLET :
	    break;
	  case BC_PERIODIC :
	    for (k=a3-n-1 ; k>=n ; k--) {
	      for (j=a2-n-1 ; j>=n ; j--) {
	        for (i=n-1 ; i>=0 ; i--) {
	          C[(k*a2+j)*a1+i] = C[(k*a2+j+1)*a1-2*n+i];
	          C[(k*a2+j+1)*a1-i-1] = C[(k*a2+j)*a1+(2*n-i-1)];
	        }
	      }
	    }
	    break;
	  case BC_SYMMETRIC :
	    for (k=a3-n-1 ; k>=n ; k--) {
	      for (j=a2-n-1 ; j>=n ; j--) {
	        for (i=n-1 ; i>=0 ; i--) {
	          C[(k*a2+j)*a1+i] = C[(k*a2+j)*a1+2*n-1-i];
	          C[(k*a2+j+1)*a1-i-1] = C[(k*a2+j+1)*a1-2*n+i];
	        }
	      }
	    }
	    break;
	  case BC_ANTISYMMETRIC :
	    for (k=a3-n-1 ; k>=n ; k--) {
	      for (j=a2-n-1 ; j>=n ; j--) {
	        for (i=n-1 ; i>=0 ; i--) {
	          C[(k*a2+j)*a1+i] = -C[(k*a2+j)*a1+2*n-1-i];
	          C[(k*a2+j+1)*a1-i-1] = -C[(k*a2+j+1)*a1-2*n+i];
	        }
	      }
	    }
	    break;
	  default :
	    fprintf(stderr, "setBoundaryConditions3D: unknown boundary condition along X axis.");
	    exit(-1);
	}

	switch(bY) {
	  case BC_ZERO :
	    for (k=a3-n-1 ; k>=n ; k--) {
	      for (j=n-1 ; j>=0 ; j--) {
	        for (i=a1-n-1 ; i>=n ; i--) {
	          C[(k*a2+j)*a1+i] = 0.;
	          C[((k+1)*a2-j-1)*a1+i] = 0.;
	        }
	      }
	    }
	    break;
	  case BC_DIRICHLET :
	    break;
	  case BC_PERIODIC :
	    for (k=a3-n-1 ; k>=n ; k--) {
	      for (j=n-1 ; j>=0 ; j--) {
	        for (i=a1-n ; i>=n ; i--) {
	          C[(k*a2+j)*a1+i] = C[((k+1)*a2-2*n+j)*a1+i];
	          C[((k+1)*a2-j-1)*a1+i] = C[(k*a2+2*n-j-1)*a1+i];
	        }
	      }
	    }
	    break;
	  case BC_SYMMETRIC :
	    for (k=a3-n-1 ; k>=n ; k--) {
	      for (j=n-1 ; j>=0 ; j--) {
	        for (i=a1-n ; i>=n ; i--) {
	          C[(k*a2+j)*a1+i] = C[(k*a2+2*n-1-j)*a1+i];
	          C[((k+1)*a2-j-1)*a1+i] = C[((k+1)*a2-2*n+j)*a1+i] ;
	        }
	      }
	    }
	    break;
	  case BC_ANTISYMMETRIC :
	    for (k=a3-n-1 ; k>=n ; k--) {
	      for (j=n-1 ; j>=0 ; j--) {
	        for (i=a1-n ; i>=n ; i--) {
	          C[(k*a2+j)*a1+i] = -C[(k*a2+2*n-1-j)*a1+i];
	          C[((k+1)*a2-j-1)*a1+i] = -C[((k+1)*a2-2*n+j)*a1+i] ;
	        }
	      }
	    }
	    break;
	  default :
	    fprintf(stderr, "setBoundaryConditions3D: unknown boundary condition along Y axis.");
	    exit(-1);
	}

	switch(bZ) {
	  case BC_ZERO :
	    for (k=n-1 ; k>=0 ; k--) {
	      for (j=a2-n-1 ; j>=n ; j--) {
	        for (i=a1-n-1 ; i>=n ; i--) {
	          C[(k*a2+j)*a1+i] = 0.;
	          C[((a3-1-k)*a2+j)*a1+i] = 0.;
	        }
	      }
	    }
	    break;
	  case BC_DIRICHLET :
	    break;
	  case BC_PERIODIC :
	    for (k=n-1 ; k>=0 ; k--) {
	      for (j=a2-n-1 ; j>=n ; j--) {
	        for (i=a1-n-1 ; i>=n ; i--) {
	          C[(k*a2+j)*a1+i] = C[((k+a3-2*n)*a2+j)*a1+i];
	          C[((a3-1-k)*a2+j)*a1+i] = C[((2*n-1-k)*a2+j)*a1+i];
	        }
	      }
	    }
	    break;
	  case BC_SYMMETRIC :
	    for (k=n-1 ; k>=0 ; k--) {
	      for (j=a2-n-1 ; j>=n ; j--) {
	        for (i=a1-n-1 ; i>=n ; i--) {
	          C[(k*a2+j)*a1+i] = C[((2*n-1-k)*a2+j)*a1+i];
	          C[((a3-1-k)*a2+j)*a1+i] = C[((a3-2*n+k)*a2+j)*a1+i];
	        }
	      }
	    }
	    break;
	  case BC_ANTISYMMETRIC :
	    for (k=n-1 ; k>=0 ; k--) {
	      for (j=a2-n-1 ; j>=n ; j--) {
	        for (i=a1-n-1 ; i>=n ; i--) {
	          C[(k*a2+j)*a1+i] = -C[((2*n-1-k)*a2+j)*a1+i];
	          C[((a3-1-k)*a2+j)*a1+i] = -C[((a3-2*n+k)*a2+j)*a1+i];
	        }
	      }
	    }
	    break;
	  default :
	    fprintf(stderr, "setBoundaryConditions3D: unknown boundary condition along Z axis.");
	    exit(-1);
	}

	return;
}

/* applyStencil3D
 *
 * Apply a 3D stencil S of size (2n+1)^3 to a 3D array A
 *
 */
void applyStencil3D(double* C, const double* S, const double* A, \
		int a1, int a2, int a3, int n) {
	int i, j, k;
	int cross = 1;

	if (a1 < 2*n+1 || a2 < 2*n+1 || a3 < 2*n+1) {
		fprintf(stderr, "applyStencil3D: too small array (%d, %d, %d) or missing ghost cells.\n", a1, a2, a3);
		exit(-1);
	}
	for (k=2*n ; k>=0 ; k--) {
	  for (j=2*n ; j>=0 ; j--) {
	    for (i=2*n ; i>=0 ; i--) {
	      if (i!=j && j!=k && k!=i && S[(k*(2*n+1)+j)*(2*n+1)+i] != 0.) {
	        cross = 0;
		printf("applyStencil3D: non-zero value (%e) at (%d, %d, %d): applying general stencil\n", S[(k*(2*n+1)+j)*(2*n+1)+i], i, j, k);
		break;
	      }
	    }
	    if (!cross) break;
	  }
	  if (!cross) break;
	}
	
	if (cross) applyCrossStencil3D(C, S, A, a1, a2, a3, n);
	else applyCubeStencil3D(C, S, A, a1, a2, a3, n);

	return;
}

/* applyCrossStencil3D
 *
 * Apply a 3D cross stencil S of size (2n+1)^3 to a 3D array A
 *
 */
void applyCrossStencil3D(double* C, const double* S, const double* A, \
		int a1, int a2, int a3, int n) {
	int i, j, k, s;

	for (k=a3-n-1 ; k>=n ; k--) {
	  for (j=a2-n-1 ; j>=n ; j--) {
	    for (i=a1-n-1 ; i>=n ; i--) {
	      C[(k*a2+j)*a1+i] = 0.;
	      for (s=2*n ; s>=0 ; s--) {
	        C[(k*a2+j)*a1+i] += S[(n*(2*n+1)+n)*(2*n+1)+s] * \
				    A[(k*a2+j)*a1+i+s-n];
	        C[(k*a2+j)*a1+i] += S[(n*(2*n+1)+s)*(2*n+1)+n] * \
				    A[(k*a2+j+s-n)*a1+i];
	        C[(k*a2+j)*a1+i] += S[(s*(2*n+1)+n)*(2*n+1)+n] * \
				    A[((k+s-n)*a2+j)*a1+i];
	        continue;
	      }
	      C[(k*a2+j)*a1+i] -= 2.*S[(n*(2*n+1)+n)*(2*n+1)+n] * \
				  A[(k*a2+j)*a1+i];
	    }
	  }
	}
	return;
}

/* applyCubeStencil3D
 *
 * Apply a general 3D stencil S of size (2n+1)^3 to a 3D array A
 *
 */
void applyCubeStencil3D(double* C, const double* S, const double* A, \
		int a1, int a2, int a3, int n) {
	int i, j, k, s1, s2, s3;

	for (k=a3-n-1 ; k>=n ; k--) {
	  for (j=a2-n-1 ; j>=n ; j--) {
	    for (i=a1-n-1 ; i>=n ; i--) {
	      C[(k*a2+j)*a1+i] = 0.;
	      for (s3=2*n ; s3>=0 ; s3--) {
	        for (s2=2*n ; s2>=0 ; s2--) {
	          for (s1=2*n ; s1>=0 ; s1--) {
	            C[(k*a2+j)*a1+i] += S[(s3*(2*n+1)+s2)*(2*n+1)+s1] * \
					A[((k+s3-n)*a2+j+s2-n)*a1+i+s1-n];
		  }
		}
	      }
	    }
	  }
	}
	return;
}

/* stencil2OLaplacian3D
 *
 * Returns a pointer to the stencil of a 2th order 3D Laplacian
 * This stencil has n=1, hence it has 3x3x3 size
 *
 */
double* stencil2OLaplacian3D(double hx, double hy, double hz) {
	double* S = arrmalloc(27);

	nullify(S, 27);
	S[Card3dO(1)] = -2. * (1./(hx*hx) + 1./(hy*hy) + 1./(hz*hz));
	S[Card3dN(1,1)] = 1./(hz*hz);
	S[Card3dS(1,1)] = 1./(hz*hz);
	S[Card3dE(1,1)] = 1./(hx*hx);
	S[Card3dW(1,1)] = 1./(hx*hx);
	S[Card3dR(1,1)] = 1./(hy*hy);
	S[Card3dF(1,1)] = 1./(hy*hy);

	return S;
}

/* stencil4OLaplacian3D
 *
 * Returns a pointer to the stencil of a 4th order 3D Laplacian
 * This stencil has n=2, hence it has 5x5x5 size
 *
 */
double* stencil4OLaplacian3D(double hx, double hy, double hz) {
	double* S = arrmalloc(125);

	nullify(S, 125);
	S[Card3dO(2)] = -5./2. * (1./(hx*hx) + 1./(hy*hy) + 1./(hz*hz));
	S[Card3dN(2,1)] = 4./3. * 1./(hz*hz);
	S[Card3dN(2,2)] = -1./12. * 1./(hz*hz);
	S[Card3dS(2,1)] = 4./3. * 1./(hz*hz);
	S[Card3dS(2,2)] = -1./12. * 1./(hz*hz);
	S[Card3dE(2,1)] = 4./3. * 1./(hx*hx);
	S[Card3dE(2,2)] = -1./12. * 1./(hx*hx);
	S[Card3dW(2,1)] = 4./3. * 1./(hx*hx);
	S[Card3dW(2,2)] = -1./12. * 1./(hx*hx);
	S[Card3dR(2,1)] = 4./3. * 1./(hy*hy);
	S[Card3dR(2,2)] = -1./12. * 1./(hy*hy);
	S[Card3dF(2,1)] = 4./3. * 1./(hy*hy);
	S[Card3dF(2,2)] = -1./12. * 1./(hy*hy);

	return S;
}

/* streamForExplicitL()
 *
 * Returns a stream to feed explicitL
 *
 */
void* streamForExplicitL(double* A, int n) {
	void *stream = streammalloc(sizeof(int)+sizeof(double*));

	SVAR(stream,int,0) = n;
	SVAR(stream,double*,sizeof(int)) = A;

	return stream;
}

/* streamForRBC3DLaplacian
 *
 * Returns a stream suited for stencil3DL and implementing the 3D Laplacian
 * with regular boundary conditions (i.e. zero, constant, periodic, symmetric
 * or antisymmetric)
 *
 */
void* streamForLaplacian3DL(int a1, int a2, int a3, int nghost, \
		int bX, int bY, int bZ, double* S, double* g) {
	void *stream = streammalloc(5*sizeof(int)+2*sizeof(double*));
	int bc = (bX << BC_X) | (bY << BC_Y) | (bZ << BC_Z);
	int N = a1*a2*a3, n=(a1-2*nghost)*(a2-2*nghost)*(a3-2*nghost);
	double *x0, *g0;

	SVAR(stream, int, 0) = bc;
	SVAR(stream, int, sizeof(int)) = a1;
	SVAR(stream, int, 2*sizeof(int)) = a2;
	SVAR(stream, int, 3*sizeof(int)) = a3;
	SVAR(stream, int, 4*sizeof(int)) = nghost;
	SVAR(stream, double*, 5*sizeof(int)) = S;
	if (g != NULL) {
		x0 = arrmalloc(N);
		restoreGhostCells3D(x0, g, a1, a2, a3, nghost);
		nullifyGhostEdges3D(x0, a1, a2, a3, nghost);
		if (bX != BC_DIRICHLET)
			nullifyGhostFace3D(x0, a1, a2, a3, nghost, FACE_X);
		if (bY != BC_DIRICHLET)
			nullifyGhostFace3D(x0, a1, a2, a3, nghost, FACE_Y);
		if (bZ != BC_DIRICHLET)
			nullifyGhostFace3D(x0, a1, a2, a3, nghost, FACE_Z);
		g0 = saveGhostCells3D(x0, a1, a2, a3, nghost);
		arrcpy(g, g0, N-n);
		free(g0);
		free(x0);
	}
	SVAR(stream, double*, 5*sizeof(int)+sizeof(double*)) = g;

	return stream;
}
