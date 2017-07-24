/****************************************************************************
 *
 * Local Approximate Gaussian Process Regression
 * Copyright (C) 2013, The University of Chicago
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
 *
 ****************************************************************************/


#include "matrix.h"
#include "util.h"
#include "linalg.h"
#include "rhelp.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <Rmath.h>
#include "gp.h"


#define SDEPS sqrt(DOUBLE_EPS)


/*
 * distance:
 * 
 * C-side version of distance_R
 */

void distance(double **X1, const unsigned int n1, double **X2,
              const unsigned int n2, const unsigned int m,
              double **D)
{
  unsigned int i,j,k;

  /* for each row of X1 and X2 */
  for(i=0; i<n1; i++) {
    for(j=0; j<n2; j++) {

      /* sum the squared entries */
      D[i][j] = 0.0;
      for(k=0; k<m; k++) {
              D[i][j] += sq(X1[i][k] - X2[j][k]);
      }

    }
  }
}



/*
 * distance_R:
 *
 * function for calculating the distance matrix between
 * the rows of X1 and X2, with output in D_out -- using
 * a built-in R interface
 */

void distance_R(double *X1_in, int *n1_in, double *X2_in,
                int *n2_in, int *m_in, double *D_out)
{
  double **X1, **X2, **D;

  /* make matrix bones */
  X1 = new_matrix_bones(X1_in, *n1_in, *m_in);
  X2 = new_matrix_bones(X2_in, *n2_in, *m_in);
  D = new_matrix_bones(D_out, *n1_in, *n2_in);

  distance(X1, *n1_in, X2, *n2_in, *m_in, D);

  /* clean up */
  free(X1);
  free(X2);
  free(D);
}


/*
 * distance_symm_R:
 *
 * function for calculating the distance matrix between
 * the rows of X1 and itself, with output in the symmetric
 * D_out matrix -- using a built-in R interface
 */

void distance_symm_R(double *X_in, int *n_in, int *m_in, double *D_out)
{
  int n, m, i, j, k;
  double **X, **D;

  /* copy integers */
  n = *n_in;
  m = *m_in;

  /* make matrix bones */
  X = new_matrix_bones(X_in, n, m);
  D = new_matrix_bones(D_out, n, n);

  /* for each row of X and itself */
  for(i=0; i<n; i++) {
    D[i][i] = 0.0;
    for(j=i+1; j<n; j++) {
      D[i][j] = 0.0;
      for(k=0; k<m; k++)
        D[i][j] += sq(X[i][k] - X[j][k]);
      D[j][i] = D[i][j];
    }
  }

  /* clean up */
  free(X);
  free(D);
}




/*
 * covar_symm:
 *
 * calculate the correlation (K) between X1 and X2 with 
 * an sepseparable.g power exponential correlation function 
 * with range d and nugget g
 */

void covar_symm(const int col, double **X, const int n, 
		    double *d, double *g, int *gi, double *mult, 
        double **K)
{
  int i, j, k;
  double multi = 1.0;

  /* calculate the covariance */
  for(i=0; i<n; i++) {
    if(mult) multi = mult[i];
    K[i][i] = 1.0 + g[gi[i]]/multi;
    for(j=i+1; j<n; j++) {
      K[i][j] = 0.0;
      for(k=0; k<col; k++) K[i][j] += sq(X[i][k] - X[j][k])/d[k];
      K[i][j] = exp(0.0 - K[i][j]);
      K[j][i] = K[i][j];
    }
  }
}



/*
 * covar:
 *
 * calculate the correlation (K) between X1 and X2 with 
 * an isotropic power exponential correlation function 
 * with range d and nugget g
 */

void covar(const int col, double **X1, const int n1, double **X2,
     const int n2, double *d, double **K)
{
  int i, j, k;

  /* calculate the covariance */
  for(i=0; i<n1; i++)
    for(j=0; j<n2; j++) {
      K[i][j] = 0.0;
      for(k=0; k<col; k++) K[i][j] += sq(X1[i][k] - X2[j][k])/d[k];
      K[i][j] = exp(0.0 - K[i][j]);
    }
}


/*
 * diff_covar:
 *
 * calculate the first and 2nd derivative (wrt d) of the correlation (K)
 * between X1 and X2 with an separable.g power exponential 
 * correlation function with range d and nugget g (though g not
 * needed)
 */

void diff_covar(const int col, double **X1, const int n1, 
        double **X2, const int n2, double *d, double **K, double ***dK)
{
  int i, j, k;
  double d2k;

  /* first copy K into each dK */
  for(k=0; k<col; k++) {
    d2k = sq(d[k]);
    for(i=0; i<n1; i++) {
      for(j=0; j<n2; j++) {
        dK[k][i][j] = K[i][j]*sq(X1[i][k] - X2[j][k])/d2k;
      }
    }
  }
}


/*
 * diff_covar_symm:
 *
 * calculate the first and 2nd derivative (wrt d) of the correlation (K)
 * between X1 and X2 with an separable.g power exponential 
 * correlation function with range d and nugget g (though g not
 * needed) -- assumes symmetric matrix
 */

void diff_covar_symm(const int col, double **X, const int n, 
        double *d, double **K, double ***dK)
{
  int i, j, k;
  double d2k;

  /* first copy K into each dK */
  for(k=0; k<col; k++) {
    d2k = sq(d[k]);
    for(i=0; i<n; i++) {
      for(j=i+1; j<n; j++) {
        dK[k][j][i] = dK[k][i][j] = K[i][j]*sq(X[i][k] - X[j][k])/d2k;
      }
      dK[k][i][i] = 0.0;
    }
  }
}


/*
 * Global variables used to accumulate data on the C-side
 * as passed from initialization and updtating R calls 
 */

unsigned int NGP = 0;
GP **gps = NULL;


/*
 * get_gp:
 *
 * returns an integer reference to a free separable.g gp 
 */

unsigned int get_gp(void)
{
  unsigned int i;
  if(NGP == 0) {
    assert(gps == NULL);
    gps = (GP**) malloc(sizeof(GP*));
    gps[0] = NULL;
    NGP = 1;
    return 0;
  } else {
    for(i=0; i<NGP; i++) {
      if(gps[i] == NULL) return i;
    }
    gps = (GP**) realloc(gps, sizeof(GP*) * (2*NGP));
    for(i=NGP; i<2*NGP; i++) gps[i] = NULL;
    NGP *= 2;
    return NGP/2;
  }
}


/*
 * deletedKGP:
 *
 * delete the dK components of the gp 
 */

void deletedKGP(GP *gp)
{
  unsigned int k;
  if(gp->dK) {
    for(k=0; k<gp->m; k++) {
      assert(gp->dK[k]);
      delete_matrix(gp->dK[k]);
    }
    free(gp->dK);  
  }
}


/*
 * deletedKGP_R:
 *
 * R-interface to code destroying dK information 
 * so that they are not updated in future calculations
 */

void deletedKGP_R(/* inputs */
      int *gpi_in)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* check if needed */
  if(! gp->dK) error("derivative info not in gp");
  
  /* call real C routine */
  deletedKGP(gp);
}


/* 
 * deleteGP:
 *
 * free the memory allocated to a separable.g gp structure 
 */


void deleteGP(GP* gp)
{
  unsigned int k;

  assert(gp);
  assert(gp->X); delete_matrix(gp->X);
  assert(gp->Z); free(gp->Z);
  assert(gp->K); delete_matrix(gp->K);
  assert(gp->Ki); delete_matrix(gp->Ki);  
  assert(gp->KiZ); free(gp->KiZ);
  deletedKGP(gp);
  assert(gp->d); free(gp->d);
  assert(gp->g); free(gp->g);
  assert(gp->gi); free(gp->gi);
  assert(gp->mult); free(gp->mult);
  if(gp->onenug != NULL) {
    for(k=0; k<gp->glen; k++) {
      assert(gp->onenug[k]);
      deleteGP(gp->onenug[k]);
    }
  }
  free(gp);
}


/* 
 * deleteGP_index:
 *
 * delete the i-th gp
 */

void deleteGP_index(unsigned int i)
{
  if(!(gps == NULL || i >= NGP || gps[i] == NULL)) { 
    deleteGP(gps[i]);
    gps[i] = NULL;
  } else error("gp %d is not allocated\n", i);
}


/*
 * deleteGP_R:
 *
 * R-interface to deleteGP
 */

void deleteGP_R(int *gp)
{
  deleteGP_index(*gp);
}


/*
 * deleteGPs:
 *
 * delete all of the gps stored in
 * the array, and destroy the array
 */

void deleteGPs(void)
{
  int i;
  for(i=0; i<NGP; i++) {
    if(gps[i]) {
      MYprintf(MYstdout, "removing gp %d\n", i);
      deleteGP(gps[i]);
    }
  }
  if(gps) free(gps);
  gps = NULL;
  NGP = 0;
}


/*
 * deleteGPs_R:
 *
 * R interface to deleteGPs
 */

void deleteGPs_R(void)
{ 
  if(gps) deleteGPs();
}


/*
 * calc_ZtKiZ:
 *
 * re-calculates phi = ZtKiZ from Ki and Z stored in
 * the GP object; also update KiZ on which it depends
 */

void calc_ZtKiZ(GP *gp) 
{
  assert(gp);
  /* phi <- t(Z) %*% Ki %*% Z */
  if(gp->KiZ == NULL) gp->KiZ = new_vector(gp->n);
  linalg_dsymv(gp->n,1.0,gp->Ki,gp->n,gp->Z,1,0.0,gp->KiZ,1);
  gp->phi = linalg_ddot(gp->n, gp->Z, 1, gp->KiZ, 1);
}


/*
 * newdKGP:
 *
 * allocate new space for dK and d2K calculations, and 
 * cancluate derivatives 
 */

void newdKGP(GP *gp)
{
  unsigned int j;
  assert(gp->dK == NULL);
  gp->dK = (double ***) malloc(sizeof(double **) * gp->m);
  for(j=0; j<gp->m; j++) gp->dK[j] = new_matrix(gp->n, gp->n);
  diff_covar_symm(gp->m, gp->X, gp->n, gp->d, gp->K, gp->dK);
}


/*
 * buildKGP_R:
 *
 * R-interface to code allocating dK information 
 * for future calculations
 */

void buildKGP_R(/* inputs */
    int *gpi_in)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* check if needed */
  if(gp->dK) error("derivative info already in gp");
  
  /* call real C routine */
  newdKGP(gp);
}


/*
 * buildGP:
 *
 * intended for newly created separable GPs, e.g., via newGP
 * does all of the correlation calculations, etc., after data and
 * parameters are defined; automaticall builds derivative info too
 */

GP* buildGP(GP *gp)
{ 
  double **Kchol, **X;
  unsigned int n, m;
  int info;

  assert(gp && gp->K == NULL);
  n = gp->n;
  m = gp->m;
  X = gp->X;

  /* build covari ance matrix */
  gp->K = new_matrix(n, n);
  covar_symm(m, X, n, gp->d, gp->g, gp->gi, gp->mult, gp->K);
  
  /* invert covariance matrix */
  gp->Ki = new_id_matrix(n);
  Kchol = new_dup_matrix(gp->K, n, n);
  info = linalg_dposv(n, Kchol, gp->Ki);
  if(info) {
#ifdef UNDEBUG
    printMatrix(gp->K, n, n, stdout);
#endif
    MYprintf(MYstdout, "d = ");
    printVector(gp->d, m, MYstdout, HUMAN);
    MYprintf(MYstdout, "g = ");
    printVector(gp->g, gp->glen, stdout, HUMAN);
    MYprintf(MYstdout, "gi = ");
    printIVector(gp->gi, n, MYstdout);
    error("bad Cholesky decomp (info=%d)", info);
  }
  gp->ldetK = log_determinant_chol(Kchol, n);
  delete_matrix(Kchol);

  /* phi <- t(Z) %*% Ki %*% Z */
  gp->KiZ = NULL;
  calc_ZtKiZ(gp);

  /* calculate derivatives ? */
  gp->dK = NULL;
  newdKGP(gp);

  /* return new structure */
  return(gp);
}


/*
 * newGP:
 *
 * allocate a new separable GP structure using the data and parameters
 * provided
 */ 

GP* newGP(const unsigned int m, const unsigned int n, double **X,
	  double *Z, double *d, double *g, const int glen, int *gi, 
    double *mult)
{
  GP* gp;

  /* new gp structure */
  gp = (GP*) malloc(sizeof(GP));
  gp->m = m;
  gp->n = n;
  gp->X = new_dup_matrix(X, n, m);
  gp->Z = new_dup_vector(Z, n);
  gp->d = new_dup_vector(d, m);
  gp->g = new_dup_vector(g, glen);
  gp->glen = glen;
  gp->gi = new_dup_ivector(gi, n);
  gp->mult = new_dup_vector(mult, n);
  gp->dn = sumv(mult, n);
  gp->K = NULL;
  gp->dK = NULL;
  gp->onenug = NULL;

  return buildGP(gp);
}


/*
 * newGP_sub:
 *
 * allocate a new GP structure using the parameters
 * provided, and the subset (rows) of the data specified by p
 *
 * restricted to one nugget for use in mleGP_nug; maybe it would
 * be simpler if this just took a gp_old argument
 */

GP* newGP_sub(const unsigned int m, const unsigned int n, int *p,
              double **X, double *Z, double *d, const double g, 
              double *mult)
{
  unsigned int i;
  GP* gp;

  /* new gp structure */
  gp = (GP*) malloc(sizeof(GP));
  gp->m = m;
  gp->n = n;
  gp-> X = new_p_submatrix_rows(p, X, n, gp->m, 0);
  gp->Z = new_vector(n);
  gp->mult = new_vector(n);
  for(i=0; i<n; i++) {
    gp->Z[i] = Z[p[i]];
    gp->mult[i] = mult[p[i]];
  }
  gp->dn = sumv(gp->mult, n);
  gp->d = new_dup_vector(d, m);
  gp->g = new_vector(1);
  gp->g[0] = g;
  gp->glen = 1;
  gp->gi = new_zero_ivector(n);
  gp->K = NULL;
  gp->dK = NULL;
  gp->onenug = NULL;

  return buildGP(gp);
}



/*
 * newGP_R:
 *
 * R-interface initializing a new separable GP, allocating and
 * assigning values to the global variables, which are
 * written over if already in use
 */

void newGP_R(/* inputs */
       int *m_in,
       int *n_in,
       double *X_in,
       double *Z_in,
       double *d_in,
       double *g_in,
       int *glen_in,
       int *gi_in,
       double *mult_in,
       
       /* outputs */
       int *gp_index)
{
  double **X;

  /* assign a new gp index */
  *gp_index = get_gp();

  /* create a new GP; */
  X = new_matrix_bones(X_in, *n_in, *m_in);
  gps[*gp_index] = newGP(*m_in, *n_in, X, Z_in, d_in, g_in, *glen_in, gi_in, mult_in);
  free(X);
}


/*
 * llikGP:
 *
 * calculate and return the log marginal likelihood of a separably
 * parameterized GP (which only really matters for the prior)
 */

double llikGP(GP *gp, double *dab, double *gab)
{
  unsigned int k;
  double llik;

  /* proportional to likelihood calculation */
  llik = 0.0 - 0.5*(((double) gp->n) * log(0.5 * gp->phi) + gp->ldetK);
  // MYprintf(MYstdout, "d=%g, g=%g, phi=%g, llik=%g\n", gp->d, gp->g, gp->phi, llik); 
  /* llik += lgamma(0.5*((double) gp->n)) - ((double) gp->n)*M_LN_SQRT_2PI; */

  /* if priors are being used; for lengthscale */
  if(dab && dab[0] > 0 && dab[1] > 0) {
    for(k=0; k<gp->m; k++) {
      if(gp->d[k] > 0) llik += dgamma(gp->d[k], dab[0], 1.0/dab[1], 1);
    }
  }

  /* if priors are being used; for nugget */
  if(gab && gab[0] > 0 && gab[1] > 0) {
    for(k=0; k<gp->glen; k++) {
      if(gp->g[k] > 0) llik += dgamma(gp->g[k], gab[0], 1.0/gab[1], 1);
    }
  }

  return(llik);
}


/*
 * llikGP_R:
 *
 * R-interface to calculate the marginal likelihood of a 
 * separable GP
 */

void llikGP_R(/* inputs */
        int *gpi_in,
        double *dab_in,
        double *gab_in,

        /* outputs */
        double *llik_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* calculate log likelihood */
  *llik_out = llikGP(gp, dab_in, gab_in);
}


/*
 * dllikGP:
 *
 * batch calculation of the gradient of the log likelihood 
 * of a separable gp parameterization, with respect to the 
 * vector lengthscale parameter, d; requires that derivatives
 * for K and Ki be pre-calculated
 */

void dllikGP(GP *gp, double *ab, double *dllik)
{
  double *KiZtwo;
  unsigned int i, j, n, k;
  double phirat ;

  /* sanity check */
  assert(gp->dK);
  assert(dllik);

  /* copy dims for fast access */
  n = gp->n;

  KiZtwo = new_vector(n);
  for(k=0; k<gp->m; k++) {

    /* deal with possible prior */
    if(ab && ab[0] > 0 && ab[1] > 0) {
      dllik[k] = (ab[0] - 1.0)/gp->d[k] - ab[1];
    } else dllik[k] = 0.0;
  
    /* dllik = - 0.5 * tr(Ki %*% dK) */
    for(i=0; i<n; i++) {
      for(j=0; j<i; j++) /* off diagonal */
        dllik[k] -= gp->Ki[i][j] * gp->dK[k][i][j];

      /* on-diagonal */
      dllik[k] -= 0.5 * gp->Ki[i][i] * gp->dK[k][i][i];
    }

    /* now third part of the expression, re-using KiZtwo */
    /* KiZtwo = dK %*% KiZ */
    linalg_dsymv(n,1.0,gp->dK[k],n,gp->KiZ,1,0.0,KiZtwo,1);
    /* now t(KiZ) %*% dK %*% KiZ */
    phirat = linalg_ddot(n, gp->KiZ, 1, KiZtwo, 1) / gp->phi;
    dllik[k] += 0.5*(gp->dn)*phirat;
  }
  
  /* clean up */
  free(KiZtwo);
}


/*
 * grad_llikGP_nug:
 *
 * batch calculation of the first derivative
 * of the log likelihood of a gp, with respect to the 
 * NUGGET parameter, g
 */

void grad_llikGP_nug(GP *gp, double *ab, double *dllik)
{
  unsigned int i, n;

  /* sanity check */
  assert(dllik); 

  /* copy dims for fast access */
  n = gp->n;
  
  /* dllik = - 0.5 * tr(Ki dK) +  (n/2phi) * Yt Ki dK Ki Y*/
  zerov(dllik, gp->glen);
  for(i=0; i<n; i++) {
    dllik[gp->gi[i]] -= 0.5 * gp->Ki[i][i] / gp->mult[i];
    dllik[gp->gi[i]] += sq(gp->KiZ[i]) * 0.5 * (gp->dn) / (gp->mult[i] * gp->phi);
  }

  if(ab && ab[0] > 0 && ab[1] > 0) {
    for(i=0; i<gp->glen; i++) {
      dllik[i] += (ab[0] - 1.0)/gp->g[i] - ab[1];
    }
  }
}


/*
 * dllikGP_nug:
 *
 * batch calculation of the first derivative
 * of the log likelihood of a gp, with respect to the 
 * NUGGET parameter, g
 *
 * MULT ADJUSTMENT NOT MADE IN THIS FUNCTION -- MUST USE 
 * GRAD LLIKGP NUG ABOVE
 */

void dllikGP_nug(GP *gp, double *ab, double *dllik, double *d2llik)
{
  unsigned int i, j, n;
  double *KiZtwo;
  double **two, **dKKidK;
  double phirat, dlp, d2lp;

  /* sanity check */
  assert(dllik); 

  /* only made for single nugget; must first make sub-GP 
     with one nugget */
  assert(gp->glen == 1);

  /* deal with possible prior */
  if(ab && ab[0] > 0 && ab[1] > 0) {
    dlp = (ab[0] - 1.0)/gp->g[0] - ab[1];
    d2lp = 0.0 - (ab[0] - 1.0)/sq(gp->g[0]);
  } else dlp = d2lp = 0.0;

  /* copy dims for fast access */
  n = gp->n;

  if(d2llik) {
    two = new_matrix(n, n);
    dKKidK = gp->Ki; 
  } else two = dKKidK = NULL;
  
  /* d2llik = - 0.5 * tr(Ki %*% [0.0 - Ki]); the first expression */
  /* dllik = - 0.5 * tr(Ki) */
  if(d2llik) *d2llik = d2lp;
  *dllik = dlp; 
  for(i=0; i<n; i++) {
    if(d2llik) {
      for(j=0; j<i; j++) { /* off diagonal */
        *d2llik += gp->Ki[i][j] * dKKidK[i][j];
        two[i][j] = two[j][i] = 2.0*dKKidK[i][j];
      }
    }
    /* on-diagonal */
    *dllik -= 0.5 * gp->Ki[i][i];
    if(d2llik) {
      *d2llik += 0.5 * gp->Ki[i][i] * dKKidK[i][i];
      two[i][i] = 2.0*dKKidK[i][i];
    }
  }

  /* now the second part of the expression: */
  /* d2llik -= 0.5 * KiZ %*% two %*% KiZ */
  if(d2llik) {
    KiZtwo = new_vector(n);
    linalg_dsymv(n,1.0,two,n,gp->KiZ,1,0.0,KiZtwo,1);
    *d2llik -= 0.5*(gp->dn)*linalg_ddot(n, gp->KiZ, 1, KiZtwo, 1) / gp->phi;
    free(KiZtwo);
  }

  /* now third part of the expression, re-using KiZtwo */
  /* now t(KiZ) %*% dK %*% KiZ */
  phirat = linalg_ddot(n, gp->KiZ, 1, gp->KiZ, 1) / gp->phi;
  if(d2llik) *d2llik += 0.5*(gp->dn)*sq(phirat);
  *dllik += 0.5*(gp->dn)*phirat;

  /* clean up */
  if(two) delete_matrix(two);
}


/*
 * dllikGP_R:
 *
 * R-interface to calculate the derivatives of the
 * likelihood of a GP - wrt lengthscale
 */

void dllikGP_R(/* inputs */
               int *gpi_in,
               double *ab_in,

               /* outputs */
               double *dllik_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* double check that derivatives have been calculated */
  if(! gp->dK) 
    error("derivative info not in gp; use newGP with dK=TRUE");

  /* calculate derivative(s) log likelihood wrt theta*/
  dllikGP(gp, ab_in, dllik_out);
}


/*
 * grad_llikGP_nug_R:
 *
 * R-interface to calculate the derivatives of the
 * likelihood of a GP - wrt the NUGGET
 *
 * deals explicitly with gradient in multiple nugget setup
 */

void grad_llikGP_nug_R(/* inputs */
               int *gpi_in,
               double *ab_in,

               /* outputs */
               double *dllik_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* calculate gradient of the log likelihood wrt multiple gs */
  grad_llikGP_nug(gp, ab_in, dllik_out);
}


/*
 * grad_llikGP_R:
 *
 * R-interface to calculate the derivatives of the
 * likelihood of a GP
 *
 * joint gradient of lengthscale and multiple nugget
 */

void grad_llikGP_R(/* inputs */
               int *gpi_in,
               double *ab_in,

               /* outputs */
               double *dllik_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* calculate derivative(s) log likelihood wrt theta*/
  dllikGP(gp, ab_in, dllik_out);
  /* calculate gradient of the log likelihood wrt multiple gs */
  grad_llikGP_nug(gp, ab_in + 2, dllik_out + gp->m);
}


/*
 * dllikGP_nug_R:
 *
 * R-interface to calculate the derivatives of the
 * likelihood of a GP - wrt the NUGGET
 */

void dllikGP_nug_R(/* inputs */
               int *gpi_in,
               double *ab_in,

               /* outputs */
               double *dllik_out,
               double *d2llik_out)
{
  GP *gp;
  double *d2llik;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];


  /* check to see if we want 2nd derivative or not */
  if(d2llik_out[0] == 1) d2llik = d2llik_out;
  else d2llik = NULL;

  /* calculate derivative of the log likelihood wrt g alone*/
  dllikGP_nug(gp, ab_in, dllik_out, d2llik);
}


/*
 * getmGP_R:
 *
 * R-interface accessing the input dimension m
 * of a GP
 */

void getmGP_R(/* inputs */
              int *gpi_in,
               /* outputs */
              int *m_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  *m_out = gp->m;
}

/*
 * getphiGP_R:
 *
 * R-interface accessing the input dimension m
 * of a GP
 */

void getphiGP_R(/* inputs */
              int *gpi_in,
               /* outputs */
              double *phi_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  *phi_out = gp->phi;
}



/*
 * getgGP_R:
 *
 * R-interface accessing the nugget of a separable
 * of a GP
 */


void getgGP_R(/* inputs */
              int *gpi_in,
              int *glen_in,
               /* outputs */
              double *g_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* check glen */
  if(*glen_in != gp->glen) error("glen (%d) does not match gp->glen (%d)", 
                                  *glen_in, gp->glen);

  dupv(g_out, gp->g, gp->glen);
}


/*
 * getglenGP_R:
 *
 * R-interface accessing the number of nuggets of a separable
 * of a GP
 */


void getglenGP_R(/* inputs */
              int *gpi_in,
               /* outputs */
              int *glen_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  *glen_out = gp->glen;
}


/*
 * getdGP_R:
 *
 * R-interface accessing the separable.g lengthscale parameter
 * of a GP
 */

void getdGP_R(/* inputs */
              int *gpi_in,
               /* outputs */
              double *d_out)
{
  GP *gp;
  unsigned int gpi;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL)
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* double check that derivatives have been calculated */
  dupv(d_out, gp->d, gp->m);
}


/*
 * newparamsGP:
 *
 * change the lengthscale and nugget parameters to the gp
 *
 * SIMIAR to newparamsGP except vectorized d and always does
 * gradient
 */ 

void newparamsGP(GP* gp, double *d, double* g)
{
  int info, m, n;
  double **Kchol;

  /* sanity check */
  assert(g >= 0);

  /* build covariance matrix */
  m = gp->m; n = gp->n;
  dupv(gp->d, d, m);
  dupv(gp->g, g, gp->glen);
  covar_symm(m, gp->X, n, gp->d, gp->g, gp->gi, gp->mult, gp->K);
  
  /* invert covariance matrix */
  id(gp->Ki, n);
  Kchol = new_dup_matrix(gp->K, n, n);
  info = linalg_dposv(n, Kchol, gp->Ki);
  if(info) {
#ifdef UNDEBUG
    printMatrix(gp->K, n, n, stdout);
#endif
    MYprintf(MYstdout, "d =");
    printVector(gp->d, m, MYstdout, HUMAN);
    MYprintf(MYstdout, "g = ");
    printVector(gp->g, gp->glen, stdout, HUMAN);
    error("bad Cholesky decomp (info=%d)", info);
  }
  gp->ldetK = log_determinant_chol(Kchol, n);
  delete_matrix(Kchol);

  /* phi <- t(Z) %*% Ki %*% Z */
  calc_ZtKiZ(gp);

  /* calculate derivatives ? */
  if(gp->dK) 
    diff_covar_symm(gp->m, gp->X, gp->n, gp->d, gp->K, gp->dK);    
}


/*
 * newparamsGP_R:
 *
 * R-interface allowing the internal/global separable.g GP representation
 * to change its parameterization without destroying the
 * memory and then re-allocating it
 */

void newparamsGP_R(/* inputs */
    int *gpi_in,
    double *d_in,
    double *g_in,
    int *glen_in)
{
  GP *gp;
  unsigned int gpi, k;
  int dsame, gsame;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* check glen */
  if(*glen_in != gp->glen) 
    error("glen (%d) does not match gp->glen (%d)", 
                                  *glen_in, gp->glen);

  /* check if any are old */
  dsame = 1;
  for(k=0; k<gp->m; k++) {
    if(d_in[k] <= 0) d_in[k] = gp->d[k];
    else if(d_in[k] != gp->d[k]) dsame = 0;
  }
  gsame = 1;
  for(k=0; k<gp->glen; k++) {
    if(g_in[k] <= 0) g_in[k] = gp->g[k];
    else if(g_in[k] != gp->g[k]) gsame = 0;
  }

  /* check if there is nothing to do bc the params are the same */
  if(dsame && gsame) return;

  /* call real C routine */
  newparamsGP(gp, d_in, g_in);
}


/*
 * utility structure for fcnnllik and fcnndllik defined below
 * for use with lbfgsb (R's optim with "L-BFGS-B" method)
 * for optimization over the lengthscale parameter only
 */

struct callinfo {
  GP *gp;
  double *dab;
  double *gab;
  int its;  /* updated but not used since lbfgsb counts fmin and gr evals */
  int verb;
};


/*
 * fcnnllik:
 * 
 * a utility function for lbfgsb (R's optim with "L-BFGS-B" method) to 
 * evaluating the separable.g GP log likelihood after changes to the 
 * lengthscale parameter 
 */

static double fcnnllik(int n, double *p, struct callinfo *info)
{
  double llik;
  int psame, k, m;

  /* sanity check */
  m = info->gp->m;
  assert(n == m || n == m + info->gp->glen);

  /* check if parameters in p are new */
  psame = 1;
  for(k=0; k<n; k++) {
    if(k < m && p[k] != info->gp->d[k]) { psame = 0; break; }
    else if(k >= m && p[k] != info->gp->g[k-m]) { psame = 0; break; }
  }

  /* update GP with new parameters */
  if(!psame) {
    (info->its)++;
    if(n == m) newparamsGP(info->gp, p, info->gp->g);
    else newparamsGP(info->gp, p, p+m);
  }

  /* evaluate likelihood with potentially new paramterization */
  llik = llikGP(info->gp, info->dab, info->gab);

  /* progress meter */
  if(info->verb > 0) {
    if(n == m) MYprintf(MYstdout, "fmin it=%d, d=(%g", info->its, p[0]);
    else MYprintf(MYstdout, "fmin it=%d, theta=(%g", info->its, p[0]);
    for(k=1; k<n; k++) MYprintf(MYstdout, " %g", p[k]);
    MYprintf(MYstdout, "), llik=%g\n", llik);
  }

  /* done */
  return 0.0-llik;
}
 

/*
 * fcnngradllik:
 * 
 * a utility function for lbfgsb (R's optim with "L-BFGS-B" method) 
 * evaluating the derivative of separable GP log likelihood after 
 * changes to the lengthscale  and nugget parameter 
 */

static void fcnngradllik(int n, double *p, double *df, struct callinfo *info)
{
  int psame, k;
  unsigned int m;

  /* sanity check */
  m = info->gp->m;
  assert(n == m + info->gp->glen);

  /* check if parameters in p are new */
  psame = 1;
  for(k=0; k<m; k++) if(p[k] != info->gp->d[k]) { psame = 0; break; }
  if(psame) for(k=m; k<n; k++) if(p[k] != info->gp->g[k-m]) { psame = 0; break; }

  /* update GP with new parameters */
  if(!psame) {
    (info->its)++;
    newparamsGP(info->gp, p, p+m);
  }

  /* evaluate likelihood with potentially new paramterization */
  dllikGP(info->gp, info->dab, df);
  grad_llikGP_nug(info->gp, info->gab, df+m);

  /* negate values */
  for(k=0; k<n; k++) df[k] = 0.0-df[k];

  /* progress meter */
  if(info->verb > 1) {
    MYprintf(MYstdout, "grad it=%d, theta=(%g", info->its, p[0]);
    for(k=1; k<n; k++) MYprintf(MYstdout, " %g", p[k]);
    MYprintf(MYstdout, ")\n");
  }
}


/*
 * fcnndllik:
 * 
 * a utility function for lbfgsb (R's optim with "L-BFGS-B" method) 
 * evaluating the derivative of separable GP log likelihood after 
 * changes to the lengthscale parameter 
 */

static void fcnndllik(int n, double *p, double *df, struct callinfo *info)
{
  int dsame, k;

  /* sanity check */
  assert(n == info->gp->m);

  /* check if parameters in p are new */
  dsame = 1;
  for(k=0; k<n; k++) if(p[k] != info->gp->d[k]) { dsame = 0; break; }

  /* update GP with new parameters */
  if(!dsame) {
    (info->its)++;
    newparamsGP(info->gp, p, info->gp->g);
  }

  /* evaluate likelihood with potentially new paramterization */
  dllikGP(info->gp, info->dab, df);

  /* negate values */
  for(k=0; k<n; k++) df[k] = 0.0-df[k];

  /* progress meter */
  if(info->verb > 1) {
    MYprintf(MYstdout, "grad it=%d, d=(%g", info->its, info->gp->d[0]);
    for(k=1; k<n; k++) MYprintf(MYstdout, " %g", info->gp->d[k]);
    MYprintf(MYstdout, "), dd=(%g", df[0]);
    for(k=1; k<n; k++) MYprintf(MYstdout, " %g", df[k]);
    MYprintf(MYstdout, ")\n");
  }
} 


/*
 * mleGP_both:
 *
 * update the separable GP to use its MLE separable
 * lengthscale and multiple nugget parameterization using the current data,
 * via the lbfgsb function 
 */

void mleGP_both(GP* gp, double* tmin, double *tmax, double *ab, 
  const unsigned int maxit, int verb, double *p, int *its, char *msg, 
  int *conv, int fromR)
{
  double rmse;
  int k, lbfgs_verb;
  double *told;

  /* create structure for lbfgsb */
  struct callinfo info;
  info.gp = gp;
  info.dab = ab;
  info.gab = ab+2;
  info.its = 0;
  info.verb = verb-6;

  /* copy the starting value */
  dupv(p, gp->d, gp->m);
  dupv(p+(gp->m), gp->g, gp->glen);
  told = new_dup_vector(p, gp->m + gp->glen);

  if(verb > 0) {
    MYprintf(MYstdout, "(theta=[%g", p[0]);
    for(k=1; k<gp->m+gp->glen; k++) MYprintf(MYstdout, ",%g", p[k]);
    MYprintf(MYstdout, "], llik=%g) ", llikGP(gp, ab, ab+2));
  }

  /* set ifail argument and verb/trace arguments */
  *conv = 0;
  if(verb <= 1) lbfgs_verb = 0;
  else lbfgs_verb = verb - 1;

  /* call the C-routine behind R's optim function with method = "L-BFGS-B" */
  MYlbfgsb(gp->m+gp->glen, p, tmin, tmax, 
         (double (*)(int, double*, void*)) fcnnllik, 
         (void (*)(int, double *, double *, void *)) fcnngradllik,
         conv, &info, its, maxit, msg, lbfgs_verb, fromR);

  /* check if parameters in p are new */
  rmse = 0.0;
  for(k=0; k<gp->m; k++) rmse += sq(p[k] - gp->d[k]);
  if(sqrt(rmse/((double) gp->m)) > SDEPS) warning("stored d not same as d-hat");
  rmse = 0.0;
  for(k=0; k<gp->glen; k++) rmse += sq(p[k+gp->m] - gp->g[k]);
  if(sqrt(rmse/((double) gp->glen)) > SDEPS) warning("stored g not same as g-hat");
  rmse = 0.0;
  for(k=0; k<gp->m+gp->glen; k++) rmse += sq(p[k] - told[k]);
  if(sqrt(rmse/((double) (gp->m + gp->glen))) < SDEPS) {
    sprintf(msg, "lbfgs initialized at minima");
    *conv = 0;
    its[0] = its[1] = 0;
  }

  /* print progress */
  if(verb > 0) {
    MYprintf(MYstdout, "-> %d lbfgsb its -> (theta=[%g", its[1], p[0]);
    for(k=1; k<gp->m+gp->glen; k++) MYprintf(MYstdout, ",%g", p[k]);
    MYprintf(MYstdout, "], llik=%g)\n", llikGP(gp, ab, ab+2));
  }

  /* clean up */
  free(told);
}


/*
 * mleGP_both_R:
 *
 * R-interface to update the separable GP to use its MLE 
 * separable lengthscale and multiple nugget 
 * parameterization using the current data
 */

void mleGP_both_R(/* inputs */
       int *gpi_in,
       int *maxit_in,
       int *verb_in,
       double *tmin_in,
       double *tmax_in,
       double *ab_in,

       /* outputs */
       double *mle_out,
       int *its_out,
       char **msg_out,
       int *conv_out)
{
  GP *gp;
  unsigned int gpi, j;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* check d against tmax and tmin */
  for(j=0; j<gp->m; j++) { 
    if(tmin_in[j] <= 0) tmin_in[j] = SDEPS;
    if(tmax_in[j] <= 0) tmax_in[j] = sq((double) gp->m);
    if(gp->d[j] > tmax_in[j]) 
      error("d[%d]=%g > tmax[%d]=%g\n", j, gp->d[j], j, tmax_in[j]);
    else if(gp->d[j] < tmin_in[j]) 
      error("d[%d]=%g < tmin[%d]=%g\n", j, gp->d[j], j, tmin_in[j]);
  }

  /* check g and tmax */
  for(j=0; j<gp->glen; j++) { 
    if(tmin_in[j+gp->m] <= 0) tmin_in[j+gp->m] = SDEPS;
    if(gp->g[j] >= tmax_in[j+gp->m]) error("g=%g >= tmax=%g\n", gp->g[j], tmax_in[j+gp->m]);
    else if(gp->g[j] <= tmin_in[j+gp->m]) error("g=%g <= tmin=%g\n", gp->g[j], tmin_in[j+gp->m]);
  }

  /* check a & b */
  if(ab_in[0] < 0 || ab_in[1] < 0 || ab_in[2] < 0 || ab_in[3] < 0) 
    error("ab_in must be a positive 4-vector");

  /* double check that derivatives have been calculated */
  if(!gp->dK) 
    error("derivative info not in gp; use newGP with dK=TRUE");  

  /* call C-side MLE */
  /* dupv(mle_out, gp->d, gp->m);
  dupv(mle_out+gp->m, gp->g, gp->glen); */ /* already done inside mleGP_both */
  mleGP_both(gp, tmin_in, tmax_in, ab_in, *maxit_in, *verb_in, mle_out,
           its_out, *msg_out, conv_out, 1);
}


/*
 * mleGP:
 *
 * update the separable GP to use its MLE arable
 * lengthscale parameterization using the current data,
 * via the lbfgsb function 
 *
 */

void mleGP(GP* gp, double* dmin, double *dmax, double *ab, 
  const unsigned int maxit, int verb, double *p, int *its, char *msg, 
  int *conv, int fromR)
{
  double rmse;
  int k, lbfgs_verb;
  double *dold;

  /* create structure for Brent_fmin */
  struct callinfo info;
  info.gp = gp;
  info.dab = ab;
  info.gab = NULL;
  info.its = 0;
  info.verb = verb-6;

  /* copy the starting value */
  dupv(p, gp->d, gp->m);
  dold = new_dup_vector(gp->d, gp->m);

  if(verb > 0) {
    MYprintf(MYstdout, "(d=[%g", gp->d[0]);
    for(k=1; k<gp->m; k++) MYprintf(MYstdout, ",%g", gp->d[k]);
    MYprintf(MYstdout, "], llik=%g) ", llikGP(gp, ab, NULL));
  }

  /* set ifail argument and verb/trace arguments */
  *conv = 0;
  if(verb <= 1) lbfgs_verb = 0;
  else lbfgs_verb = verb - 1;

  /* call the C-routine behind R's optim function with method = "L-BFGS-B" */
  MYlbfgsb(gp->m, p, dmin, dmax, 
         (double (*)(int, double*, void*)) fcnnllik, 
         (void (*)(int, double *, double *, void *)) fcnndllik,
         conv, &info, its, maxit, msg, lbfgs_verb, fromR);

  /* check if parameters in p are new */
  rmse = 0.0;
  for(k=0; k<gp->m; k++) rmse += sq(p[k] - gp->d[k]);
  if(sqrt(rmse/((double) gp->m)) > SDEPS) warning("stored d not same as d-hat");
  rmse = 0.0;
  for(k=0; k<gp->m; k++) rmse += sq(p[k] - dold[k]);
  if(sqrt(rmse/((double) gp->m)) < SDEPS) {
    sprintf(msg, "lbfgs initialized at minima");
    *conv = 0;
    its[0] = its[1] = 0;
  }

  /* print progress */
  if(verb > 0) {
    MYprintf(MYstdout, "-> %d lbfgsb its -> (d=[%g", its[1], gp->d[0]);
    for(k=1; k<gp->m; k++) MYprintf(MYstdout, ",%g", gp->d[k]);
    MYprintf(MYstdout, "], llik=%g)\n", llikGP(gp, ab, NULL));
  }

  /* clean up */
  free(dold);
}


/*
 * mleGP_R:
 *
 * R-interface to update the separable GP to use its MLE 
 * separable lengthscale parameterization using the current data
 *
 * SIMPLIFIED compared to mleGP_R since only the (separable)
 * lengthscale is supported; for the nugget see mleGP_nug_R
 */

void mleGP_R(/* inputs */
       int *gpi_in,
       int *maxit_in,
       int *verb_in,
       double *dmin_in,
       double *dmax_in,
       double *ab_in,

       /* outputs */
       double *mle_out,
       int *its_out,
       char **msg_out,
       int *conv_out)
{
  GP *gp;
  unsigned int gpi, j;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* check d against dmax and dmin */
  for(j=0; j<gp->m; j++) { 
    if(dmin_in[j] <= 0) dmin_in[j] = SDEPS;
    if(dmax_in[j] <= 0) dmax_in[j] = sq((double) gp->m);
    if(gp->d[j] > dmax_in[j]) 
      error("d[%d]=%g > dmax[%d]=%g\n", j, gp->d[j], j, dmax_in[j]);
    else if(gp->d[j] < dmin_in[j]) 
      error("d[%d]=%g < dmin[%d]=%g\n", j, gp->d[j], j, dmin_in[j]);
  }

  /* check a & b */
  if(ab_in[0] < 0 || ab_in[1] < 0) error("ab_in must be a positive 2-vector");

  /* double check that derivatives have been calculated */
  if(!gp->dK) 
    error("derivative info not in gp; use newGP with dK=TRUE");  

  /* call C-side MLE */
  /* dupv(mle_out, gp->d, gp->m); */ /* already done inside mleGP */
  mleGP(gp, dmin_in, dmax_in, ab_in, *maxit_in, *verb_in, mle_out,
           its_out, *msg_out, conv_out, 1);
}


/*
 * utility structure for fcnnllik_nug defined below
 * for use with Brent_fmin (R's optimize) or uniroot
 *
 * SIMPLIFIED compared to callinfo in gp.c because it only does the nugget
 */

struct callinfo_nug {
  GP *gp;
  double *gab;
  int its;
  int verb;
};


/*
 * fcnnllik_nug:
 * 
 * a utility function for Brent_fmin (R's optimize) to apply to the separable.g 
 * GP log likelihood after changes to the nugget parameter 
 *
 * SIMPLIFIED compared to fcnnllik in gp.c since it only does the nugget
 */

static double fcnnllik_nug(double x, struct callinfo_nug *info)
{
  double llik;
  assert(infp->gp->glen == 0);
  (info->its)++;
  newparamsGP(info->gp, info->gp->d, &x);
  llik = llikGP(info->gp, NULL, info->gab);
  if(info->verb > 1)
    MYprintf(MYstdout, "fmin it=%d, g=%g, llik=%g\n", info->its, info->gp->g[0], llik);
  return 0.0-llik;
} 


/*
 * Ropt_nug:
 *
 * use R's Brent Fmin routine (from optimize) to optimize
 *
 * SIMPLIFIED compared to Ropt in GP because it only does the nugget
 */

double Ropt_nug(GP* gp, double tmin, double tmax, 
                   double *ab, char *msg, int *its, int verb)
{
  double tnew, th;
  double Tol = SDEPS;

  /* sanity check */
  assert(tmin < tmax);

  /* get parameter */
  assert(gp->glen == 0);
  th = gp->g[0];

  /* create structure for Brent_fmin */
  struct callinfo_nug info;
  info.gp = gp;
  info.gab = ab;
  info.its = 0;
  info.verb = verb;

  /* call the C-routine behind R's optimize function */
  while(1) { /* check to make sure solution is not on boundary */
   tnew = Brent_fmin(tmin, tmax, (double (*)(double, void*)) fcnnllik_nug, &info, Tol);  
   if(tnew > tmin && tnew < tmax) break;
   if(tnew == tmin) { /* left boundary found */
    tmin *= 2;
    if(verb > 0) MYprintf(MYstdout, "Ropt: tnew=tmin, increasing tmin=%g\n", tmin);
   } else { /* right boundary found */
    tmax /= 2.0;
    if(verb > 0) MYprintf(MYstdout, "Ropt: tnew=tmax, decreasing tmax=%g\n", tmax);
  }
  /* check that boundaries still valid */
  if(tmin >= tmax) error("unable to opimize in fmin()");
  } 

  /* check that last value agrees with GP parameterization */
  if(gp->g[0] != tnew) newparamsGP(gp, gp->d, &tnew);

  /* possible print message and return */
  if(verb > 0) MYprintf(MYstdout, "Ropt %s: told=%g -[%d]-> tnew=%g\n",
      msg, th, info.its, tnew);

  *its += info.its;
  return(tnew);
}


/*
 * mleGP_onenug:
 *
 * calculate the MLE with respect to the nugget parameter;
 * derivatives for the Newton method are calculated on the fly;
 *
 * assumes there is only one nugget parameter (in this multiple
 * nugget setup)
 */

double mleGP_onenug(GP* gp, double tmin, double tmax, double *ab, 
             int verb, int *its)
{
  double tnew, dllik, d2llik, llik_init, llik_new, adj, rat;
  double th; //, th_init;
  double *gab, *dab;
  int restoredKGP;

  /* set priors based on Theta */
  dab = NULL;
  gab = ab;

  /* initialization */
  *its = 0;
  restoredKGP = 0;
  /* th_init = */ th = gp->g[0];

  /* check how close we are to tmin */
  if(fabs(th - tmin) < SDEPS) {
    if(verb > 0) MYprintf(MYstdout, "(g=%g) -- starting too close to min (%g)\n", th, tmin);
    goto alldone;
  }

  /* initial likelihood calculation */
  llik_init = llikGP(gp, dab, gab);

  /* initial printing */
  if(verb > 0) 
      MYprintf(MYstdout, "(g=%g, llik=%g) ", gp->g[0], llik_init);
  if(verb > 1) MYprintf(MYstdout, "\n");

  while(1) { /* checking for improved llik */
    while(1) {  /* Newton step(s) */
      llik_new = 0.0-1e300*1e300;
      while(1) {  /* Newton proposal */

        /* calculate first and second derivatives */
        dllikGP_nug(gp, gab, &dllik, &d2llik);

        /* check for convergence by root */
        if(fabs(dllik) < SDEPS) {
          if(*its == 0) {
            if(verb > 0) MYprintf(MYstdout, "-- Newton not needed\n");
            goto alldone;
          } else goto newtondone;
        }

        /* Newton update */
        rat = dllik/d2llik; adj = 1.0; (*its)++;

        /* check if we're going the right way */
        if((dllik < 0 && rat < 0) || (dllik > 0 && rat > 0)) {
          if(!gp->dK && restoredKGP == 0) { 
            deletedKGP(gp); restoredKGP = 1; 
          }
          th = Ropt_nug(gp, tmin, tmax, ab, "[slip]", its, verb); goto mledone; 
        } else tnew = th - adj*rat;  /* right way: Newton: */

        /* check that we haven't proposed a tnew out of range */
        while((tnew <= tmin || tnew >= tmax) && adj > SDEPS) {
          adj /= 2.0; tnew = th - adj*rat;
        }

        /* if still out of range, restart? */
        if(tnew <= tmin || tnew >= tmax) { 
          if(!gp->dK && restoredKGP == 0) { 
            deletedKGP(gp); restoredKGP = 1; 
          }
          th = Ropt_nug(gp, tmin, tmax, ab, "[range]", its, verb);
          goto mledone;
        } else break;
      } /* end inner while -- Newton proposal */

      /* else, resets gp->g = tnew */
      if(!gp->dK && restoredKGP == 0) { 
        deletedKGP(gp); restoredKGP = 1; 
      }
      newparamsGP(gp, gp->d, &tnew);

      /* print progress */
      if(verb > 1) MYprintf(MYstdout, "\ti=%d g=%g, c(a,b)=(%g,%g)\n", 
                            *its, tnew, ab[0], ab[1]);      

      /* check for convergence, and break or update */
      if(fabs(tnew - th) < SDEPS) break;
      else th = tnew;

      /* check for max its */
      if(*its >= 100) {
        if(verb > 0) warning("Newton 100/max iterations");
        /* could also call Ropt here as last resort */
       goto alldone;
      }
    } /* end middle while -- Newton step */

    /* sanity check check that we went in the right direction */
newtondone:
    llik_new = llikGP(gp, dab, gab);
    if(llik_new < llik_init-SDEPS) { 
      if(verb > 0) MYprintf(MYstdout, "llik_new = %g\n", llik_new);
      llik_new = 0.0-1e300*1e300;
      if(!gp->dK && restoredKGP == 0) { 
        deletedKGP(gp); restoredKGP = 1; 
      }
      th = Ropt_nug(gp, tmin, tmax, ab, "[dir]", its, verb); 
      goto mledone;
    } else break;
  } /* outer improved llik check while(1) loop */

  /* capstone progress indicator */
mledone:
  if(!R_FINITE(llik_new)) llik_new = llikGP(gp, dab, gab);
  if(verb > 0) {
    MYprintf(MYstdout, "-> %d Newtons -> (g=%g, llik=%g)\n", //, tdiff*ldiff=%g)\n", 
            *its, gp->g[0], llik_new); //, fabs(llik_init - llik_new)*fabs(th_init - th));
  }
  /* to determine convergence we report a two its instead of the truth */
  if(fabs(llik_new - llik_init) < SDEPS) *its = 2;

  /* return theta-value found */
alldone:
  if(restoredKGP) newdKGP(gp);
  return th;
}


/*
 * mleGP_nug:
 *
 * calculate the MLE with respect to the multiple nugget parameter;
 * derivatives for the Newton method are calculated on the fly;
 *
 * if there is only one nugget, then mleGP_onenug is used instead
 */

void mleGP_nug(GP *gp, double tmin, double tmax, double *ab, 
            int verb, int *its, double *g)
{
  unsigned int k, nk, i, isum, init;
  int its_temp;
  int *p;

  if(gp->glen == 1) *g = mleGP_onenug(gp, tmin, tmax, ab, verb, its);
  else {

    /* possibly allocate new one-nugget GPs */
    if(gp->onenug == NULL) {
      p = new_ivector(gp->n);
      gp->onenug = malloc(sizeof(GP*) * gp->glen);
      for(k=0; k<gp->glen; k++) {
        nk = 0;
        for(i=0; i<gp->n; i++) {
          if(gp->gi[i] == k) p[nk++] = i;
        }
        gp->onenug[k] = newGP_sub(gp->m, nk, p, gp->X, gp->Z, gp->d, gp->g[k], gp->mult);
      }
      free(p);
      init = 1;
    } else init = 0;

    isum = 0;
    /* for each onenugget GP */
    #pragma omp parallel for private(its_temp) reduction(+:isum)
    for(k=0; k<gp->glen; k++) {

      /* only initialize if not just allocated above */
      if(!init) newparamsGP(gp->onenug[k], gp->d, gp->onenug[k]->g);
      
      /* now performed onenug nugget optimization */
      g[k] = mleGP_onenug(gp->onenug[k], tmin, tmax, ab, verb, &its_temp);
      g[k] *= gp->onenug[k]->phi / ((double) gp->onenug[k]->n); 
      g[k] /= gp->phi / ((double) gp->n);

      // #pragma omp atomic
      isum += its_temp;
    }
    *its = isum;

    /* don't forget to make assignment in big gp */
    newparamsGP(gp, gp->d, g);
  }
}


/*
 * mleGP_nug_R:
 *
 * R-interface to update the separable GP to use its MLE 
 * nugget parameterization using the current data
 *
 */

void mleGP_nug_R(/* inputs */
       int *gpi_in,
       int *verb_in,
       double *tmin_in,
       double *tmax_in,
       double *ab_in,

       /* outputs */
       double *mle_out,
       int *its_out)
{
  GP *gp;
  unsigned int gpi, k;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* check theta and tmax */
  if(*tmin_in <= 0) *tmin_in = SDEPS;
  for(k=0; k<gp->glen; k++) { 
    if(gp->g[k] >= *tmax_in) error("g=%g >= tmax=%g\n", gp->g[k], *tmax_in);
    else if(gp->g[k] <= *tmin_in) error("g=%g <= tmin=%g\n", gp->g[k], *tmin_in);
  }

  /* check a & b */
  if(ab_in[0] < 0 || ab_in[1] < 0) error("ab_in must be a positive 2-vector");

  /* call C-side MLE */
  mleGP_nug(gp, *tmin_in, *tmax_in, ab_in, *verb_in, its_out, mle_out);
}


/*
 * jmleGP:
 *
 * calculate joint mle for separable.g lengthscale (d) and nugget (g) 
 * by a coordinite-wise search, iterating over d and g searches via mleGP
 * and mleGP_nug
 */

void jmleGP(GP *gp, int maxit, double *dmin, double *dmax, 
               double *grange, double *dab, double *gab, int verb, 
               int *dits, int *gits, int *dconv, int fromR) 
  {
    unsigned int i;
    int dit[2], git;
    char msg[60];
    double *d, *g;

    /* sanity checks */
    assert(gab && dab);
    assert(dmin && dmax && grange);

    /* auxillary space for d-parameter values(s) */
    d = new_vector(gp->m);
    /* and g parameters */
    g = new_vector(gp->glen);

    /* loop over coordinate-wise iterations */
    *dits = *gits = 0;
    for(i=0; i<100; i++) {
      mleGP(gp, dmin, dmax, dab, maxit, verb, d, dit, msg, dconv, fromR);
      if(dit[1] > dit[0]) dit[0] = dit[1];
      *dits += dit[0];
      mleGP_nug(gp, grange[0], grange[1], gab, verb, &git, g);
      *gits += git;
      if((git <= 2*gp->glen && (dit[0] <= gp->m+1 && *dconv == 0)) || *dconv > 1) break;
    }
    if(i == 100 && verb > 0) warning("max outer its (N=100) reached");

    /* clean up */
    free(d);
    free(g);
  }


/*
 * jmleGP_R:
 *
 * R-interface to update the GP to use its joint MLE (lengthscale
 * and nugget) parameterization using the current data
 */

void jmleGP_R(/* inputs */
       int *gpi_in,
       int *maxit_in,
       int *verb_in,
       double *dmin_in,
       double *dmax_in,
       double *grange_in,
       double *dab_in,
       double *gab_in, 

       /* outputs */
       double *d_out,
       double *g_out,
       int *dits_out,
       int *gits_out,
       int *dconv_out)
{
  GP *gp;
  unsigned int gpi, k;

  /* get the cloud */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];

  /* check theta and tmax */
  assert(grange_in[0] >= 0 && grange_in[0] < grange_in[1]);
  for(k=0; k<gp->m; k++) {
    assert(dmin_in[k] >= 0 && dmin_in[k] < dmax_in[k]);
    if(gp->d[k] < dmin_in[k] || gp->d[k] > dmax_in[k])
      error("gp->d[%d]=%g outside drange[%d]=[%g,%g]", 
        k, gp->d[k], k, dmin_in[k], dmax_in[k]);
  }
  for(k=0; k<gp->glen; k++) {
    if(gp->g[k] < grange_in[0] || gp->g[k] > grange_in[1])
      error("gp->g[k]=%g outside grange=[%g,%g]", 
        k, gp->g[k], grange_in[0], grange_in[1]);
  }

  /* double check that derivatives have been calculated */
  if(! gp->dK) 
    error("derivative info not in gp; use newGP with dK=TRUE");

  /* call C-side MLE */
  jmleGP(gp, *maxit_in, dmin_in, dmax_in, grange_in, dab_in, gab_in, *verb_in, 
    dits_out, gits_out, dconv_out, 1);

  /* write back d and g */
  dupv(d_out, gp->d, gp->m);
  dupv(g_out, gp->g, gp->glen);
}


/* 
 * pred_generic:
 *
 * a function that captures the predictive equations without reference
 * struct gp objects.  Created so that code can be shared between GP
 * and GPsep objects, and beyond
 */

void pred_generic(const unsigned int n, const double phidf, double *Z, 
  double **Ki, const unsigned int nn, double **k, double *mean, double **Sigma)
{
  int i, j;
  double **ktKi, **ktKik;

  /* ktKi <- t(k) %*% util$Ki */
  ktKi = new_matrix(n, nn);
  linalg_dsymm(CblasRight,nn,n,1.0,Ki,n,k,nn,0.0,ktKi,nn);
  /* ktKik <- ktKi %*% k */
  ktKik = new_matrix(nn, nn);
  linalg_dgemm(CblasNoTrans,CblasTrans,nn,nn,n,
               1.0,k,nn,ktKi,nn,0.0,ktKik,nn);

  /* mean <- ktKi %*% Z */
  linalg_dgemv(CblasNoTrans,nn,n,1.0,ktKi,nn,Z,1,0.0,mean,1);

  /* Sigma <- phi*(Sigma - ktKik)/df */
  for(i=0; i<nn; i++) {
    Sigma[i][i] = phidf * (Sigma[i][i] - ktKik[i][i]);
    for(j=0; j<i; j++)
      Sigma[j][i] = Sigma[i][j] = phidf * (Sigma[i][j] - ktKik[i][j]);
  }

  /* clean up */
  delete_matrix(ktKi);
  delete_matrix(ktKik);
}


/*
 * predGP:
 *
 * return the student-t predictive equations,
 * i.e., parameters to a multivatiate t-distribution
 * for XX predictive locations of dimension (n*m)
 */

void predGP(GP* gp, unsigned int nn, double **XX, int *ggi, double *mean, 
      double **Sigma, double *df, double *llik)
{
  unsigned int m, n;
  double **k;
  double phidf;

  /* easier referencing for dims */
  n = gp->n;  m = gp->m;

  /* variance (s2) components */
  *df = gp-> dn; 
  phidf = gp->phi/(*df);

  /* calculate marginal likelihood (since we have the bits) */
  *llik = 0.0 - 0.5*((*df) * log(0.5 * gp->phi) + gp->ldetK);
  /* continuing: - (gp->dn)*M_LN_SQRT_2PI;*/

  /* k <- covar(X1=X, X2=XX, d=Zt$d, g=0) */
  k = new_matrix(n, nn);
  covar(m, gp->X, n, XX, nn, gp->d, k);
  /* Sigma <- covar(X1=XX, d=Zt$d, g=Zt$g) */
  covar_symm(m, XX, nn, gp->d, gp->g, ggi, NULL, Sigma);
  
  /* call generic function that would work for all GP covariance specs */
  pred_generic(n, phidf, gp->Z, gp->Ki, nn, k, mean, Sigma);

  /* clean up */
  delete_matrix(k);
}



/* 
 * new_predutil_generic_lite:
 *
 * a function allocates space and calculate portions of the GP predictive 
 * equations without reference struct gp objects.  Created so that code can 
 * be shared between GP and GPsep objects, and beyond
 */

void new_predutil_generic_lite(const unsigned int n, double ** Ki, 
  const unsigned int nn, double **k, double ***ktKi, double **ktKik)
{
  unsigned int i, j;

  /* ktKi <- t(k) %*% util$Ki */
  *ktKi = new_matrix(n, nn);
  linalg_dsymm(CblasRight,nn,n,1.0,Ki,n,k,nn,0.0,*ktKi,nn);
  /* ktKik <- diag(ktKi %*% k) */
  *ktKik = new_zero_vector(nn); 
  for(i=0; i<nn; i++) for(j=0; j<n; j++) (*ktKik)[i] += (*ktKi)[j][i]*k[j][i];
}


/*
 * new_predutilGP_lite:
 *
 * utility function that allocates and calculate useful vectors 
 * and matrices for prediction; used by predGP_lite and dmus2GP
 */

void new_predutilGP_lite(GP *gp, unsigned int nn, double **XX, 
  double ***k, double ***ktKi, double **ktKik)
{
  /* k <- covar(X1=X, X2=XX, d=Zt$d, g=0) */
  *k = new_matrix(gp->n, nn);
  covar(gp->m, gp->X, gp->n, XX, nn, gp->d, *k);
  
  /* call generic function that would work for all GP covariance specs */
  new_predutil_generic_lite(gp->n, gp->Ki, nn, *k, ktKi, ktKik);
}


/*
 * predGP_lite:
 *
 * return the student-t predictive equations,
 * i.e., parameters to a multivatiate t-distribution
 * for XX predictive locations of dimension (n*m);
 * lite because sigma2 not Sigma is calculated
 */

void predGP_lite(GP* gp, unsigned int nn, double **XX, int *ggi, 
  double *mean, double *sigma2, double *df, double *llik)
{
  unsigned int i;
  double **k, **ktKi;
  double *ktKik;
  double phidf;
  
  /* sanity checks */
  assert(df);
  *df = gp->dn; 

  /* utility calculations */
  new_predutilGP_lite(gp, nn, XX, &k, &ktKi, &ktKik);

  /* mean <- ktKi %*% Z */
  if(mean) linalg_dgemv(CblasNoTrans,nn,gp->n,1.0,ktKi,nn,gp->Z,
                        1,0.0,mean,1);

  /* Sigma <- phi*(Sigma - ktKik)/df */
  /* *df = n - m - 1.0; */  /* only if estimating beta */
  if(sigma2) {
    phidf = gp->phi/(*df);
    for(i=0; i<nn; i++) sigma2[i] = phidf * (1.0 + gp->g[ggi[i]] - ktKik[i]);
  }

  /* calculate marginal likelihood (since we have the bits) */
  /* might move to updateGP if we decide to move phi to updateGP */
  if(llik) *llik = 0.0 - 0.5*(((double) gp->n) * log(0.5* gp->phi) + 
    gp->ldetK);
  /* continuing: - (gp->dn)*M_LN_SQRT_2PI;*/

  /* clean up */
  delete_matrix(k);
  delete_matrix(ktKi);
  free(ktKik);
}


/*
 * predGP_R:
 *
 * R-interface to C-side function that 
 * returns the student-t predictive equations,
 * i.e., parameters to a multivatiate t-distribution
 * for XX predictive locations of dimension (n*m)
 * using the stored GP parameterization
 */

void predGP_R(/* inputs */
        int *gpi_in,
        int *m_in,
        int *nn_in,
        double *XX_in,
        int *ggi_in,
        int *lite_in,
        
        /* outputs */
        double *mean_out,
        double *Sigma_out,
        double *df_out,
        double *llik_out)
{
  GP* gp;
  unsigned int gpi;
  double **Sigma, **XX;

  /* get the gp */
  gpi = *gpi_in;
  if(gps == NULL || gpi >= NGP || gps[gpi] == NULL) 
    error("gp %d is not allocated\n", gpi);
  gp = gps[gpi];
  if((unsigned) *m_in != gp->m) 
    error("ncol(X)=%d does not match GP/C-side (%d)", *m_in, gp->m);

  /* sanity check and XX representation */
  XX = new_matrix_bones(XX_in, *nn_in, *m_in);
  if(! *lite_in) Sigma = new_matrix_bones(Sigma_out, *nn_in, *nn_in);
  else Sigma = NULL;

  /* call the C-only Predict function */
  if(*lite_in) predGP_lite(gp, *nn_in, XX, ggi_in, mean_out, Sigma_out, df_out,
                              llik_out);
  else predGP(gp, *nn_in, XX, ggi_in, mean_out, Sigma, df_out, llik_out);
  
  /* clean up */
  free(XX);
  if(Sigma) free(Sigma);
}

