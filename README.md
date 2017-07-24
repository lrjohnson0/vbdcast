This directory contains supplementary material for Johnson, et al (2016)
"Phenomenological forecasting of disease incidence using heteroskedastic
Gaussian processes: a dengue case study".  The directory structure is explained
as follows.  In each case, see the README therein for more details.

R: source files for hetGP in the R language for statistical computing. 

src: source files for hetGP in C; the README contains compilation instructions for
building shared objects that can be linked into R, as required for the files in
the R directory.  These files are derived from the laGP package on CRAN.

examples: contains R scripts and data files that are used in those scripts to
perform the GP-based analyses in the paper.  The code in this directory
duplicates many of the pdfs in the results directory.

results: pdf-based visualizations for GP and GLM-based figures, a subset of
which are provided in the manuscript.  In particular, iquitos.pdf and
sanjuan.pdf contain the "movies" of each four-weekly prediction under the hetGP
which are genreated by gp_dengue_iquitos.R and gp_dengue_sanjuan.R in the
examples directory, respectively.
