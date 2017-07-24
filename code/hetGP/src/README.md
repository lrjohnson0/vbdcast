# vbdcast: Vector-borne disease forecasting

To compile the library in the format used by the R interface functions in the
`../R` sub-directory execute the following.

```
R CMD SHLIB -o hetGP.so *.c
```

* It is highly recommended to use an R which is linked to an optimized linear
algrbra (BLAS).  Our own executable, used in the experiments for the timings
reported in the paper, linked to Intel's MKL.  With default linear algebra
instead the timings could rise from under ten minutes (for both Iquitos and San
Juan) to more than an hour.
