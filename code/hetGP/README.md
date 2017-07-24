# vbdcast: Vector-borne disease forecasting

The following sub-directories contain source code and dengue forecasting via heteroskedastic Gaussian processes

* `src`: source files for `hetGP` in `C`; the README therein contains compilation instructions for building shared objects that can be linked into R, as required for the files in the `R` directory.  These files are derived from the laGP package on CRAN.

* `R`: source files for hetGP in the `R` language for statistical computing.

* `dengue`: contains R scripts and data files that are used in those scripts to
perform the GP-based analyses in the dengue forecasting paper.  The code in this directory duplicates many of the pdfs in the `~/results` directory.