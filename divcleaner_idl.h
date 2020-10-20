/* idl_bicgstab.c
 *
 * IDL wrapper over bicgstab.c
 *
 */
 
#include "bicgstab.h"

/* poisson
 * 
 * Arguments in argv (argc = 9)
 * 
 * (double*) x_out       < nx * ny * nz >   | Solution of the Poisson equation
 * 
 * (double*) x0_arr      < nx * ny * nz >   | Initial guess
 * (int*)    nx, ny, nz                     | Size of the box
 * (int*)    bx, by, bz                     | Boundary conditions (0: zero,
 *                                          | 1: Dirichlet, 2: periodic,
 *                                          | 3: symmetric, 4: antisymmetric)
 *                                          | (default: <2, 2, 0>)
 * (double*) thr                            | Iteration threshold
 * (int*)    maxiter                        | Max. number of iterations
 *
 Optional args. (NULL for default values)
 * 
 * (double*) b_arr       < nx * ny * nz >   | Right-hand-side of the equation
 *                                          | (default: zero)
 * (double*) d_arr <(nx+1)*(ny+1)*(nz+1)>   | Boundary (for Dirichlet conditions)
 * (int*)    hx, hy, hz                     | Cell size (default: <1, 1, 1>)
 */
void              poisson(int /* argc */, void** /* argv */);

/* divergence
 * 
 * Arguments in argv (argc = 6)
 * 
 * (double*) div_out     < nx * ny * nz >   | Divergence
 * 
 * (double*) Bb1     < (nx+1) * ny * nz >   | Vector field (X-component)
 * (double*) Bb2     < nx * (ny+1) * nz >   | Vector field (Y-component)
 * (double*) Bb3     < nx * ny * (nz+1) >   | Vector field (Z-component)
 * (int*)    nx, ny, nz                     | Size of the box
 *
 Optional args. (NULL for default values)
 * 
 * (int*)    hx, hy, hz                     | Cell size (default: <1, 1, 1>)
 */
void           divergence(int /* argc */, void** /* argv */);

/* fixTopBottomBoundary
 * 
 * Arguments in argv (argc = 6)
 * 
 * (double*) Bb3_out < nx * ny * (nz+1) >   | Corrected v.field (Z-component)
 * 
 * (double*) Bb1     < (nx+1) * ny * nz >   | Vector field (X-component)
 * (double*) Bb2     < nx * (ny+1) * nz >   | Vector field (Y-component)
 * (double*) Bb3     < nx * ny * (nz+1) >   | Vector field (Z-component)
 * (int*)    nx, ny, nz                     | Size of the box
 *
 Optional args. (NULL for default values)
 * 
 * (int*)    hx, hy, hz                     | Cell size (default: <1, 1, 1>)
 */
void fixTopBottomBoundary(int /* argc */, void** /* argv */);

/* phiToBbc
 * 
 * Arguments in argv (argc = 6)
 * 
 * (double*) Bb1_out < (nx+1) * ny * nz >   | Vector field (X-component)
 * (double*) Bb2_out < nx * (ny+1) * nz >   | Vector field (Y-component)
 * (double*) Bb3_out < nx * ny * (nz+1) >   | Vector field (Z-component)
 * 
 * (double*) phi     < nx * ny * nz >       | Phi (B = gradient phi)
 * (int*)    nx, ny, nz                     | Size of the box
 *
 Optional args. (NULL for default values)
 * 
 * (int*)    hx, hy, hz                     | Cell size (default: <1, 1, 1>)
 */
void             phiToBbc(int /* argc */, void** /* argv */);

/* killdiv
 * 
 * Arguments in argv (argc = 13)
 * 
 * (double*) Bb1_out < (nx+1) * ny * nz >   | Vector field (X-component)
 * (double*) Bb2_out < nx * (ny+1) * nz >   | Vector field (Y-component)
 * (double*) Bb3_out < nx * ny * (nz+1) >   | Vector field (Z-component)

 * (double*) Bb1     < (nx+1) * ny * nz >   | Vector field (X-component)
 * (double*) Bb2     < nx * (ny+1) * nz >   | Vector field (Y-component)
 * (double*) Bb3     < nx * ny * (nz+1) >   | Vector field (Z-component)
 * (int*)    nx, ny, nz                     | Size of the box
 *
 Optional args. (NULL for default values)
 * (int*)    hx, hy, hz                     | Cell size (default: <1, 1, 1>)
 * (int*)    bx, by, bz                     | Boundary conditions (0: zero,
 *                                          | 1: Dirichlet, 2: periodic,
 *                                          | 3: symmetric, 4: antisymmetric)
 *                                          | (default: <2, 2, 0>)
 * (double*) d_arr <(nx+1)*(ny+1)*(nz+1)>   | Boundary (for Dirichlet conditions)
 * (double*) x0_arr      < nx * ny * nz >   | Initial guess for phi
 * (double*) thr                            | Iteration threshold
 * (int*)    maxiter                        | Max. number of iterations
 */
void              killdiv(int /* argc */, void** /* argv */);
