/* bicgstab.c
 *
 * Implementation of the Bi-CGSTAB (biconjugate gradient stabilized) method
 * for the solution of the Ax = b linear system.
 *
 */

#ifndef _BICGSTAB_H
#define _BICGSTAB_H 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define max(x,y)       (((x) > (y)) ? (x) : (y))
#define min(x,y)       (((x) < (y)) ? (x) : (y))
#define SVAR(s,t,n)    (*((t*)(s+n)))

/*
 * Card3d* definitions provide the index of the central cell of a
 * (2*n+1)^3 array and of its (non-diagonal) neighbours located at
 * relative distance k
 *
 */
#define Card3dO(n)    n*((2*n+1)*(2*n+2)+1)               /* (x0, y0, z0) */
#define Card3dN(n,k)  (n-k)*(2*n+1)*(2*n+1)+n*(2*n+1)+n   /* (x0, y0, z0-1) */
#define Card3dS(n,k)  (n+k)*(2*n+1)*(2*n+1)+n*(2*n+1)+n   /* (x0, y0, z0+1) */
#define Card3dE(n,k)  n*(2*n+1)*(2*n+1)+n*(2*n+1)+n+k     /* (x0+1, y0, z0) */
#define Card3dW(n,k)  n*(2*n+1)*(2*n+1)+n*(2*n+1)+n-k     /* (x0-1, y0, z0) */
#define Card3dF(n,k)  n*(2*n+1)*(2*n+1)+(n-k)*(2*n+1)+n   /* (x0, y0-1, z0) */
#define Card3dR(n,k)  n*(2*n+1)*(2*n+1)+(n+k)*(2*n+1)+n   /* (x0, y0+1, z0) */

/*
 * Boundary conditions
 *
 * They are stored in one single integer. The first four LSB contain the
 * boundary condition in the X axis, the following four bits contain the
 * boundary condition in the Y axis, and finally follow the Z axis.
 *
 * Extra bits are unused.
 *
 */
#define BC_ZERO            0x00
#define BC_DIRICHLET       0x01
#define BC_PERIODIC        0x02
#define BC_SYMMETRIC       0x03
#define BC_ANTISYMMETRIC   0x04
#define BC                 0x0f
#define BC_X               0x00
#define BC_Y               0x04
#define BC_Z               0x08

/* Faces */
#define FACE_X             0x00
#define FACE_Y             0x01
#define FACE_Z             0x02

/* Contents
 *
 * Typedefs
 *    |
 *    \- Lambda(b, x, n, stream)
 *
 * Functions
 *    |
 *    |- solve(x, A, b, x0, thr, n)
 *    |- solveL(x, L, stream, b, x0, thr, n)
 *    |
 *    |- A = arrmalloc(n)
 *    |- nullify(A, n)
 *    |- oneify(A, n)
 *    |- arrcpy(A, B, n)
 *    |- dot(C, A, B, a, b, c)
 *    |- add(C, A, r, B, n)
 *    |- N = norm2(v, n)
 *    |
 *    |- explicitL(b, x, n, stream)
 *    |- stencil3DL(b, x, n, stream)
 *    |
 *    |- nonLinearPart(b, L, a1, a2, a3, n, stream)
 *    |
 *    |- dotadd(C, r, A, B, a, b, c)
 *    |- mult(C, s, A, n)
 *    |
 *    |- stream = streammalloc(n)
 *    |- addGhostFrame3D(C, A, a1, a2, a3, n)
 *    |- removeGhostFrame3D(C, A, a1, a2, a3, n)
 *    |- nullifyGhostFrame3D(C, a1, a2, a3, n)
 *    |- oneifyGhostFrame3D(C, a1, a2, a3, n)
 *    |- nullifyGhostEdges3D(C, a1, a2, a3, n)
 *    |- nullifyGhostFace3D(C, a1, a2, a3, n, f)
 *    |- G = saveGhostCells3D(A, a1, a2, a3, n)
 *    |- popGhostCells3D(A, G, a1, a2, a3, n)
 *    |- restoreGhostCells3D(A, G, a1, a2, a3, n)
 *    |- setBoundaryConditions3D(C, a1, a2, a3, n, bX, bY, bZ)
 *    |
 *    |- applyStencil3D(C, S, A, a1, a2, a3, n)
 *    |- applyCrossStencil3D(C, S, A, a1, a2, a3, n)
 *    |- applyCubeStencil3D(C, S, A, a1, a2, a3, n)
 *    |
 *    |- S = stencil2OLaplacian3D(hx, hy, hz)
 *    |- S = stencil4OLaplacian3D(hx, hy, hz)
 *    |- stream = streamForExplicitL(A, n)
 *    \- stream = streamForLaplacian3DL(a1, a2, a3, nghost, bX, bY, bZ, S, g)
 */

/* Lambda
 *
 * Lambda is the type for a function that evaluates the lambda operator
 * on a vector x, resulting on a vector b (vector sizes: n).
 *
 * Any additional parameters needed to construct the lambda operator can
 * be passed through the stream pointer.
 *
 * See explicitL for an example.
 *
 */
typedef void (*Lambda)(double* /* b */, const double* /* x */, int /* n */, \
		const void* /* stream */);

/* solve
 *
 * Solve the system Ax = b, with intial guess x0, iterating till
 * threshold; shape(A) = (n,n)
 *
 * A has to be know explicitly
 *
 */
void solve(double* /* x */, const double* /* A */, const double* /* b */, \
	       	const double* /* x0 */, double /* thr */, int /* maxiter */, \
		int /* n */);

/* solveL
 *
 * Solve the system Ax = b, with intial guess x0, iterating till
 * threshold; shape(A) = (n,n)
 *
 * Ax is computed through the Lambda-operator L
 *
 */
void solveL(double* /* x */, Lambda /* L */, const void* /* stream */, \
	       	const double* /* b */, const double* /* x0 */, \
		double /* thr */, int /* maxiter */, int /* n */);

/* arrmalloc
 *
 * Allocate memory for an array of size n
 *
 */
double* arrmalloc(int /* n */);

/* nullify
 *
 * Set to zero all values of array
 *
 */
void nullify(double* /* A */, int /* n */);

/* oneify
 *
 * Set to one all values of array
 *
 */
void oneify(double* /* A */, int /* n */);

/* arrcpy
 *
 * Copy array
 *
 */
void arrcpy(double* /* A */, const double* /* B */, int /* n */);

/* dot
 *
 * Dot product of two 2D-arrays; shape(A) = (a,b) and shape(B) = (b,c)
 *
 */
void dot(double* /* C */, const double* /* A */, const double* /* B */, \
		int a /* a */, int /* b */, int /* c */);

/* add
 *
 * Add two arrays:
 * 
 *   C = A + r B
 *
 */
void add(double* /* C */, const double* /* A */, double /* r */, \
		const double* /* B */, \
	       	int /* n */);

/* norm2
 *
 * Returns the Euclidean norm of a vector
 *
 */
double norm2(const double* /* v */, int /* n */);

/* normmax
 *
 * Returns the maximum norm of a vector
 *
 */
double normmax(const double* /* v */, int /* n */);

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
void explicitL(double* /* b */, const double* /* x */, int /* n */, \
		const void* /* stream */);

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
void stencil3DL(double* /* b */, double* /* x */, int /* n */, \
		const void* /* stream */);

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
void nonLinearPart(double* /* b*/, Lambda /* L */, \
		int /* a1 */, int /* a2 */, int /* a3 */, \
		int /* n */, const void* /* stream */);

/* dotadd
 *
 * Adds the dot product of two arrays to a third array:
 *
 *   C <-- C + r < A ; B >
 *
 */
void dotadd(double* /* C */, double /* r */, const double* /* A */, \
		const double* /* B */, int a /* a */, int /* b */, \
		int /* c */);

/* mult
 *
 * Multiply an array by a scalar
 *
 */
void mult(double* /* C */, double /* s */, const double* /* A */, \
		int /* n */);

/* streammalloc
 *
 * Allocate memory for a stream of length n (bytes)
 *
 */
double* streammalloc(int /* n */);

/* addGhostFrame3D
 *
 * Add ghost cells to an array
 *
 */
void addGhostFrame3D(double* /* C */, const double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* removeGhostFrame3D
 *
 * Remove ghost cells from an array
 *
 */
void removeGhostFrame3D(double* /* C */, const double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* nullifyGhostFrame3D
 *
 * Set to zero ghost cells of 3D array
 *
 */
void nullifyGhostFrame3D(double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* oneifyGhostFrame3D
 *
 * Set to one ghost cells of 3D array
 *
 */
void oneifyGhostFrame3D(double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* nullifyGhostEdges3D
 *
 * Set to zero ghost corners of 3D array
 *
 */
void nullifyGhostEdges3D(double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* nullifyGhostFace3D
 *
 * Set to zero specified ghost face of 3D array
 *
 */
void nullifyGhostFace3D(double* /* C */, \
		int /* a1 */, int /* a2 */, int /* a3 */, \
		int /* n */, int /* f */);


/* saveGhostCells3D
 *
 * Returns an array with ghost cells
 *
 */
double* saveGhostCells3D(const double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* popGhostCells3D
 *
 * Restore ghost cells to array A and free memory, or just free memory
 * if C = NULL
 *
 */
void popGhostCells3D(double* /* A */, double* /* G */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* restoreGhostCells3D
 *
 * Restore ghost cells to array A and keep G (don't free memory)
 *
 */
void restoreGhostCells3D(double* /* A */, const double* /* G */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* setBoundaryConditions3D()
 *
 * Set the ghost cells of array to satisfy given boundary conditions
 *
 */
void setBoundaryConditions3D(double* /* C */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */, \
		int /* bX */, int /* bY */, int /* bZ */);

/* applyStencil3D
 *
 * Apply a 3D stencil S of size (2n+1)^3 to a 3D array A
 *
 */
void applyStencil3D(double* /* C */, const double* /* S */,
		const double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* applyCrossStencil3D
 *
 * Apply a 3D cross stencil S of size (2n+1)^3 to a 3D array A
 *
 */
void applyCrossStencil3D(double* /* C */, const double* /* S */,
		const double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* applyCubeStencil3D
 *
 * Apply a general 3D stencil S of size (2n+1)^3 to a 3D array A
 *
 */
void applyCubeStencil3D(double* /* C */, const double* /* S */,
		const double* /* A */, \
		int /* a1 */, int /* a2 */, int /* a3 */, int /* n */);

/* stencil2OLaplacian3D
 *
 * Returns a pointer to the stencil of a 2th order 3D Laplacian
 * This stencil has n=2, hence it has 5x5x5 size
 *
 */
double* stencil2OLaplacian3D(double /* hx */, double /* hy */, \
		double /* hz */);

/* stencil4OLaplacian3D
 *
 * Returns a pointer to the stencil of a 4th order 3D Laplacian
 * This stencil has n=2, hence it has 5x5x5 size
 *
 */
double* stencil4OLaplacian3D(double /* hx */, double /* hy */, \
		double /* hz */);

/* streamForExplicitL()
 *
 * Returns a stream to feed explicitL
 *
 */
void* streamForExplicitL(double* /* A */, int /* n */);

/* streamForRBC3DLaplacian
 *
 * Returns a stream suited for stencil3DL and implementing the 3D Laplacian
 * with regular boundary conditions (i.e. zero, constant, periodic, symmetric
 * or antisymmetric)
 *
 */
void* streamForLaplacian3DL(int /* a1 */, int /* a2 */, int /* a3 */, \
		int /* nghost */, int /* bX */, int /* bY */, int /* bZ */, \
		double* /* S */, double* /* g */);

#endif /* bicgstab.h */
