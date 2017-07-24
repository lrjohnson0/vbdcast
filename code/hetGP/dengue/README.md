# vbdcast: Vector-borne disease forecasting

This sub-directory contains the data files and the `R` implementation specific to the hetGP empirical work reported in the Johnson, et al (2016) paper on phenomenological forecasting of dengue incidence.

* The main files are gp_dengue_iquitos.R and gp_dengue_sanjuan.R.  

* As long as the `C` library shared object is compiled as described in the README for the `../src` directory, and the `R` files are sourced from this examples directory, the results of sourcing those files will be to produce the movies (iquitos.pdf and sanjuan.pdf) stored in the `../../results` sub-direcotry.

* Code here also produces the output files which are processed by the results scripts that generate the other figures in the paper (e.g., `dengue_comparison_plots.R`).
