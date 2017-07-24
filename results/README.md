# vbdcast: Vector-borne disease forecasting

This sub-directory contains the pdf output supporting the Johnson, et al (2016) paper on phenomenological forecasting of dengue incidence. 

* To see how these are created, run the `R` code in the `dengue` sub-directories within `../code/*`.

* Of particular interest may be iquitos.pdf and sanjuan.pdf are the "movies" referenced in the main manuscript.  These are produced by wrappying a pdf command around a source-execution of `gp_dengue_iquitos.R` and `gp_dengue_sanjuan.R`, respectively, withing `../code/hetGP/dengue/`, say.

* Filenames with `glm` in them correspond to the GLM-based forecasts.  Those witout correspond to `hetGP`.

* The files in *csvouts* are useful in creating these pdfs.

* The code in `R` contains the aggregation code that reads the CSV files in `csvouts` to combine `hetGP` and GLM outputs.
