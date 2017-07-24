
#ifndef __GP_H__
#define __GP_H__

typedef struct gp {
  double **X;       /* design matrix */
  double **K;       /* covariance between design points */
  double **Ki;      /* inverse of K */
  double ***dK;     /* gradient of K */
  double ldetK;     /* log of the determinant of K */
  double *Z;        /* response vector */
  double *KiZ;      /* Ki %*% Z */
  unsigned int m;   /* number of cols in X */
  unsigned int n;   /* number of rows in X; length of Z */
  double *d;        /* arable lengthscale parameter to correlation */
  double *g;        /* nugget parameter to correlation */
  int glen;         /* number of nuggets in g */
  int *gi;          /* index into g of length n */
  double *mult;     /* multiplicity parameter for each obs */
  double dn;        /* effective number of obs: dn = sum(mult) */
  double phi;       /* t(Z) %*% Ki %*% Z = t(Z) %*% KiZ, used for s2 */
  struct gp **onenug;      /* separate GPs for each nugget (subset) */
} GP;


GP* newGP(const unsigned int m, const unsigned int n, double **X,
    double *Z, double *d, double *g, const int glen, int *gi, double *mult);
void newGP_R(int *m_in, int *n_in, double *X_in, double *Z_in,
    double *d_in, double *g_in, int *glen_in, int *gi_in, double *mult_in, 
    int *gp_index);
GP* newGP_sub(const unsigned int m, const unsigned int n, int *p,
    double **X, double *Z, double *d, const double g, double *mult);
unsigned int get_gp(void);
void deletedKGP(GP *gp);
void deleteGP(GP* gp);
void deleteGP_index(unsigned int i);
void deleteGP_R(int *gp);
void deleteGPs(void);
void deleteGPs_R(void);
void calc_ZtKiZ_(GP *gp);
void newdKGP(GP *gp);
GP* buildGP(GP *gp);
double llikGP(GP *gp, double *dab, double *gab);
void llikGP_R(int *gpi_in, double *dab_in, double *gab_in,
        double *llik_out);
void dllikGP(GP *gp, double *ab, double *dllik);
void dllikGP_nug(GP *gp, double *ab, double *dllik, double *d2llik);
void dllikGP_R(int *gpi_in, double *ab_in, double *dllik_out);
void dllikGP_nug_R(int *gpi_in, double *ab_in, double *dllik_out,
     double *d2llik_out);
void getmGP_R(int *gpi_in, int *m_out);
void getgGP_R(int *gpi_in, int *gplen_in, double *g_out);
void getdGP_R(int *gpi_in, double *d_out);
void newparamsGP(GP* gp, double *d, double *g);
void newparamsGP_R(int *gpi_in, double *d_in, double *g_in, 
     int *gplen_in);
void jmleGP(GP *gp, int maxit, double *dmin, double *dmax, 
     double *grange, double *dab, double *gab, int verb, 
     int *dits, int *gits, int *dconv, int fromR);
void mleGP(GP* gp, double* dmin, double *dmax, double *ab, 
     const unsigned int maxit, int verb, double *p, int *its, 
     char *msg, int *conv, int fromR);
void mleGP_nug(GP* gp, double tmin, double tmax, double *ab, 
     int verb, int *its, double *g);
void mleGP_nug_R(int *gpi_in, int *verb_in, double *tmin_in,
     double *tmax_in, double *ab_in, double *mle_out, int *its_out);
void predGP(GP* gp, unsigned int nn, double **XX, int *ggi, 
     double *mean, double **Sigma, double *df, double *llik);
void new_predutilGP_lite(GP *gp, unsigned int nn, double **XX, 
     double ***k, double ***ktKi, double **ktKik);
void predGP_lite(GP* gp, unsigned int nn, double **XX, int *ggi,
     double *mean, double *sigma2, double *df, double *llik);
void predGP_R(int *gpi_in, int *m_in, int *nn_in, double *XX_in, 
     int *ggi_in, int *lite_in, double *mean_out, double *Sigma_out, 
     double *df_out, double *llik_out);


void covar_symm(const int col, double **X, const int n, 
     double *d, double *g, int *gi, double *mult, double **K);
void covar(const int col, double **X1, const int n1, double **X2,
     const int n2, double *d, double **K);
void diff_covar(const int col, double **X1, const int n1, 
     double **X2, const int n2, double *d, double **K, double ***dK);
void diff_covar_symm(const int col, double **X, const int n, 
     double *d, double **K, double ***dK);

#endif
