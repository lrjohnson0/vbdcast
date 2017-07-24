This sub-directory to the supplementary material for Johnson, et al (2016)
"Phenomenological forecasting of disease incidence using heteroskedastic
Gaussian processes: a dengue case study" contains the R implementation of the
hetGP library.  In particular, gp.R contains interfaces to the underlying C
library implemented in the src directory.  Compilation instructions can be
found in src, and example uses of this library is in the examples
sub-directory.

To compile the library in the format used by the R i interface functions in the
R sub-directory execute the following.

R CMD SHLIB -o hetGP.so *.c

It is highly recommended to use an R which is linked to an optimized linear
algrbra (BLAS).  Our own executable, used in the experiments for the timings
reported in the paper, linked to Intel's MKL.  With default linear algebra
instead the timings could rize from under ten minutes (for both Iquitos and San
Juan) to more than an hour.
