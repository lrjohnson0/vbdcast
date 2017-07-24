# dyn.load("../src/hetGP.so")
# source("../R/prior.R")

## distance:
##
## calculate the distance matrix between the rows of X1
## and X2, or between x1 and itself when X2=NULL

distance <- function(X1, X2=NULL)
  {
    ## coerse arguments and extract dimensions
    X1 <- as.matrix(X1)
    n1 <- nrow(X1)
    m <- ncol(X1)

    if(is.null(X2)) {

      ## calculate D
      outD <- .C("distance_symm_R",
                 X = as.double(t(X1)),
                 n = as.integer(n1),
                 m = as.integer(m),
                 D = double(n1 * n1))

      ## return the distance matrix
      return(matrix(outD$D, ncol=n1, byrow=TRUE))

    } else {

      ## coerse arguments and extract dimensions
      X2 <- as.matrix(X2)
      n2 <- nrow(X2)

      ## check inputs
      if(ncol(X1) != ncol(X2)) stop("col dim mismatch for X1 & X2")

      ## calculate D
      outD <- .C("distance_R",
                 X1 = as.double(t(X1)),
                 n1 = as.integer(n1),
                 X2 = as.double(t(X2)),
                 n2 = as.integer(n2),
                 m = as.integer(m),
                 D = double(n1 * n2))

      ## return the distance matrix
      return(matrix(outD$D, ncol=n2, byrow=TRUE))
    }
  }


## newGP:
##
## build an initial separable.g GP representation on the C-side
## using the X-Z data and d/g paramterization.
## '@export
newGP <- function(X, Z, d, g, gi=NULL, mult=NULL)
  {
    n <- nrow(X)
    m <- ncol(X)
    if(is.null(n)) stop("X must be a matrix")
    if(length(Z) != n) stop("must have nrow(X) = length(Z)")
    if(length(d) == 1) d <- rep(d, m)
    else if(length(d) != m) stop("must have length(d) = ncol(X)")

    if(!is.null(gi) && length(gi) != n)
      stop("length(gi) must be equal to length(Z)")

    if(length(g) == 1) gi <- rep(1, n)
    else {
      ugi <- unique(gi)
      if(length(setdiff(1:length(g), ugi)) > 0)
        stop("unique gi values must span 1:length(g)")
    }

    if(is.null(mult)) mult <- rep(1.0, n)
    else if(length(mult) != n || any(mult <= 0)) stop("mult must be a positive n-vector")

    out <- .C("newGP_R",
              m = as.integer(m),
              n = as.integer(n),
              X = as.double(t(X)),
              Z = as.double(Z),
              d = as.double(d),
              g = as.double(g),
              glen = as.integer(length(g)),
              gi = as.integer(gi-1),
              mult = as.double(mult),
              gpi = integer(1))

    ## return C-side GP index
    return(out$gpi)
  }


## buildkGP:
##
## allocates/calculates the C-side derivative info (only) for particular
## separable.g GP

buildKGP <- function(gpi)
  {
    .C("buildKGP_R",
       gpi = as.integer(gpi),
       PACKAGE = "laGP")
    invisible(NULL)
  }


## deletedkGP:
##
## deletes the C-side derivative info (only) for particular separable.g GP

deletedkGP <- function(gpi)
  {
    .C("deletedKGP_R",
       gpi = as.integer(gpi),
       PACKAGE = "laGP")
    invisible(NULL)
  }

## deleteGP:
##
## deletes the C-side of a particular separable.g GP

deleteGP <- function(gpi)
  {
    .C("deleteGP_R",
       gpi = as.integer(gpi))
    invisible(NULL)
  }


## deleteGPs:
##
## deletes all gps on the C side

deleteGPs <- function()
  {
    .C("deleteGPs_R")
    invisible(NULL)
  }


## llikGP:
##
## calculate the log likelihood of the GP

llikGP <- function(gpi, dab=c(0,0), gab=c(0,0))
  {
    r <- .C("llikGP_R",
            gpi = as.integer(gpi),
            dab = as.double(dab),
            gab = as.double(gab),
            llik = double(1))

    return(r$llik)
  }


## getmGP
##
## acces the input dimension of a separable.g GP
##
## totall new to GP

getmGP <- function(gpi)
  {
    .C("getmGP_R", gpi = as.integer(gpi), m = integer(1))$m
  }


## getphiGP
##
## acces the input dimension of a separable.g GP
##
## totall new to GP

getphiGP <- function(gpi)
  {
    .C("getphiGP_R", gpi = as.integer(gpi), phi = double(1))$phi
  }


## getdGP
##
## acces the separable.g lengthscale of a arable gp
##
## totally new to GP

getdGP <- function(gpi)
  {
    m <- getmGP(gpi)
    .C("getdGP_R", gpi = as.integer(gpi), d = double(m))$d
  }


## getgGP
##
## acces the input dimension of a separable.g GP
##
## totally new to GP

getgGP <- function(gpi)
  {
    glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen
    .C("getgGP_R", gpi = as.integer(gpi), glen = glen, g = double(glen))$g
  }


## dllikGP:
##
## calculate the first and second derivative of the
## log likelihood of the GP with respect to d, the
## lengthscale parameter
##
## SIMILAR to dllikGP except with vector d and gp
## isntead of gp

dllikGP <- function(gpi, ab=c(0,0), param=c("d", "g"), d2nug=FALSE)
  {
    param <- match.arg(param)
    if(param == "d") {
      dim <- getmGP(gpi)
      r <- .C("dllikGP_R",
            gpi = as.integer(gpi),
            ab = as.double(ab),
            d = double(dim))
      return(r$d)
    } else {
      if(d2nug) d2 <- 1
      else d2 <- 0
      r <- .C("dllikGP_nug_R",
            gpi = as.integer(gpi),
            ab = as.double(ab),
            d = double(1),
            d2 = as.double(d2))
      if(d2nug) return(list(d=r$d, d2=r$r2))
      else return(r$d)
    }
  }


## gradllikGP:
##
## calculate the first and second derivative of the
## log likelihood of the GP with respect to d, the
## lengthscale parameter


gradllikGP <- function(gpi, dab=c(0,0), gab=c(0,0))
  {
    m <- getmGP(gpi)
    glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen
    r <- .C("grad_llikGP_R",
            gpi = as.integer(gpi),
            ab = as.double(c(dab, gab)),
            dtheta = double(m+glen))
    return(r$dtheta)
  }


## newparamsGP:
##
## change the separable.g GP lengthscale and nugget parameerization
## (without destroying the object and creating a new one)

newparamsGP <- function(gpi, d, g=-1)
  {
    if(all(d <= 0) & all(g < 0)) stop("one of d or g must be new")
    m <- getmGP(gpi)
    if(length(d) != m) stop("length(d) !=", m)

    glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen
    if(length(g) != glen) stop("length(g) !=", glen)

    ## if(g == -1) {
      ## glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen
      ## g <- rep(g, glen)
    ## }

    r <- .C("newparamsGP_R",
            gpi = as.integer(gpi),
            d = as.double(d),
            g = as.double(g),
            glen = as.integer(length(g)))

    invisible(NULL)
  }


## mleGP.R:
##
## updates the separable.g GP to use its MLE lengthscale
## parameterization using the current data;
##
## differs substantially from mleGP in that L-BFGS-B from
## optim is used to optimize over the separable.g lengthscale;
## an option is also provided to include the nugget in that
## optimization, or do to a mleGP style profile optimization
## for the the nugget instead
##
## in this ".R" version the optim command is used; in mleGP
## an internal C-side call to lbfgsb is used

mleGP.R <- function(gpi, param=c("d", "g", "both"),
                  tmin=sqrt(.Machine$double.eps),
                  tmax=-1, ab=c(0,0), maxit=100, verb=0)
  {
    param <- match.arg(param)

    if(param == "both") { ## lengthscale and multiple nugget jointly with L-BFGS-B

      theta <- c(getdGP(gpi), getgGP(gpi))
      m <- getmGP(gpi)
      if(length(ab) != 4 || any(ab < 0)) stop("ab should be a positive 4-vector")

      ## objective
      f <- function(theta, gpi, ab, m)
        {
          ## cat("theta=(", paste(signif(theta, 3), collapse=","), ")\n", sep="")
          newparamsGP(gpi, d=theta[1:m], g=theta[(m+1):length(theta)])
          -llikGP(gpi, dab=ab[1:2], gab=ab[3:4])
        }
      ## gradient of objective
      g <- function(theta, gpi, ab, m)
        {
          newparamsGP(gpi, d=theta[1:m], g=theta[(m+1):length(theta)])
          -gradllikGP(gpi, dab=ab[1:2], gab=ab[3:4])
        }

      ## for compatibility with mleGP
      tmax[tmax < 0] <- Inf

      ## check tmin and tmax
      if(length(tmin) == 2 && length(tmin) < length(theta))
        tmin <- c(rep(tmin[1], m), rep(tmin[2], length(theta)-m))
      if(length(tmax) == 2 && length(tmax) < length(theta))
        tmax <- c(rep(tmax[1], m), rep(tmax[2], length(theta)-m))

      ## possibly print progress meter
      if(verb > 0) {
        cat("(d=[", paste(signif(theta, 3), collapse=","), "], llik=",
          llikGP(gpi, dab=ab[1:2], gab=ab[3:4]), ") ", sep="")
      }

      ## call R's optim function
      out <- optim(theta, fn=f, gr=g, method="L-BFGS-B",
        control=list(trace=max(verb-1,0), maxit=maxit), lower=tmin, upper=tmax,
        gpi=gpi, ab=ab, m=m)

      ## sanity check completion of scheme
      if(sqrt(mean((out$par - c(getdGP(gpi), getgGP(gpi)))^2)) > sqrt(.Machine$double.eps))
        warning("stored c(d,g) not same as theta-hat")

      ## check that we moved somewhere
      if(sqrt(mean((out$par - theta)^2)) < sqrt(.Machine$double.eps)) {
        out$convergence <- 0
        out$counts <- c(0,0)
        out$message <- "optim initialized at minima"
      }

      ## possibly print progress meter
      if(verb > 0) {
        cat("-> ", out$counts[1], " lbfgs.R its -> (theta=[",
          paste(signif(theta, 3), collapse=","), "], llik=",
          llikGP(gpi, dab=ab), ")\n", sep="")
      }

    } else if(param == "d") { ## lengthscale with L-BFGS-B given nugget

      theta <- getdGP(gpi)
      if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector")

      ## objective
      f <- function(theta, gpi, dab)
        {
          newparamsGP(gpi, d=theta)
          -llikGP(gpi, dab=dab)
        }
      ## gradient of objective
      g <- function(theta, gpi, dab)
        {
          newparamsGP(gpi, d=theta)
          -dllikGP(gpi, param="d", ab=dab)
        }

      ## for compatibility with mleGP
      tmax[tmax < 0] <- Inf

      ## possibly print progress meter
      if(verb > 0) {
        cat("(d=[", paste(signif(theta, 3), collapse=","), "], llik=",
          llikGP(gpi, dab=ab), ") ", sep="")
      }

      ## call R's optim function
      out <- optim(theta, fn=f, gr=g, method="L-BFGS-B",
        control=list(trace=max(verb-1,0), maxit=maxit), lower=tmin, upper=tmax,
        gpi=gpi, dab=ab)

      ## sanity check completion of scheme
      if(sqrt(mean((out$par - getdGP(gpi))^2)) > sqrt(.Machine$double.eps))
        warning("stored d not same as d-hat")

      ## check that we moved somewhere
      if(sqrt(mean((out$par - theta)^2)) < sqrt(.Machine$double.eps)) {
        out$convergence <- 0
        out$counts <- c(0,0)
        out$message <- "optim initialized at minima"
      }

      ## possibly print progress meter
      if(verb > 0) {
        cat("-> ", out$counts[1], " lbfgs.R its -> (d=[",
          paste(signif(theta, 3), collapse=","), "], llik=",
          llikGP(gpi, dab=ab), ")\n", sep="")
      }

    }
    else { ## nugget conditionally on lengthscale

      ## sanity check
      if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector");

      glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen

      r <- .C("mleGP_nug_R",
            gpi = as.integer(gpi),
            verb = as.integer(verb),
            tmin = as.double(tmin),
            tmax = as.double(tmax),
            ab = as.double(ab),
            g = double(glen),
            its = integer(1))
    }

    ## build object for returning
    if(param == "d")
      return(list(d=out$par, its=max(out$counts), msg=out$message, conv=out$convergence))
    else if(param == "both")
      return(list(theta=out$par, its=max(out$counts), msg=out$message, conv=out$convergence))
    else
      return(list(g=r$g, its=r$its))
  }


## mleGP:
##
## updates the separable.g GP to use its MLE lengthscale
## parameterization using the current data;
##
## differs substantially from mleGP in that lbfgsb is used
## to optimize over the separable.g lengthscale;
## an option is also provided to include the nugget in that
## optimization, or do to a mleGP style profile optimization
## for the the nugget instead
##
## this is a mostly C verision
## '@export
mleGP <- function(gpi, param=c("d", "g", "both"),
                  tmin=rep(sqrt(.Machine$double.eps), 2),
                  tmax=c(-1, 1), ab=rep(0, 4), maxit=100, verb=0)
  {
    param <- match.arg(param)

    if(param == "both") { ## lengthscale and multiple nugget jointly with L-BFGS-B

      ## sanity checking, and length checking for tmax and tmin
      m <- getmGP(gpi)
      glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen
      if(length(tmax) == 2) tmax <- c(rep(tmax[1], m), rep(tmax[2], glen))
      else if(length(tmax) != m+glen) stop("length(tmax) should be 2 or m+glen")
      if(length(tmin) == 2) tmin <- c(rep(tmin[1], m), rep(tmin[2], glen))
      else if(length(tmin) != m+glen) stop("length(tmin) should be 2 or m+glen")
      if(length(ab) != 4 || any(ab < 0)) stop("ab should be a positive 4-vector")

      out <- .C("mleGP_both_R",
                gpi = as.integer(gpi),
                maxit = as.integer(maxit),
                verb = as.integer(verb),
                tmin = as.double(tmin),
                tmax = as.double(tmax),
                ab = as.double(ab),
                par = double(m+glen),
                counts = integer(2),
                msg = character(1),
                convergence = integer(1))

      ## sanity check completion of scheme
      if(sqrt(mean((out$par - c(getdGP(gpi), getgGP(gpi)))^2)) > sqrt(.Machine$double.eps))
        warning("stored theta not same as theta-hat")

    } else if(param == "d") { ## lengthscale with L-BFGS-B given nugget

      ## sanity checking
      m <- getmGP(gpi)
      if(length(tmax) == 1 || (length(tmax) == 2 && m != 2)) tmax <- rep(tmax[1], m)
      else if(length(tmax) != m) stop("length(tmax) should be 1 or m")
      if(length(tmin) == 1 || (length(min ) == 2 && m != 2)) tmin <- rep(tmin[1], m)
      else if(length(tmin) != m) stop("length(tmin) should be 1 or m")
      if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector")
        
      out <- .C("mleGP_R",
                gpi = as.integer(gpi),
                maxit = as.integer(maxit),
                verb = as.integer(verb),
                dmin = as.double(tmin),
                dmax = as.double(tmax),
                ab = as.double(ab),
                par = double(m),
                counts = integer(2),
                msg = paste(rep(0,60), collapse=""),
                convergence = integer(1))

      ## sanity check completion of scheme
      if(sqrt(mean((out$par - getdGP(gpi))^2)) > sqrt(.Machine$double.eps))
        warning("stored d not same as theta-hat")
    }
    else { ## nugget conditionally on lengthscale

      ## sanity check
      if(length(ab) != 2 || any(ab < 0)) stop("ab should be a positive 2-vector");

      glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen

      r <- .C("mleGP_nug_R",
            gpi = as.integer(gpi),
            verb = as.integer(verb),
            tmin = as.double(tmin),
            tmax = as.double(tmax),
            ab = as.double(ab),
            g = double(glen),
            its = integer(1))
    }

    ## build object for returning
    if(param == "both") return(list(theta=out$par, its=max(out$counts), msg=out$msg, conv=out$convergence))
    else if(param == "d") return(list(d=out$par, its=max(out$counts), msg=out$msg, conv=out$convergence))
    else return(list(g=r$g, its=r$its))
  }


## jmleGP.R:
##
## joint MLE for lengthscale (d) and nugget (g) parameters;
## updates the internal GP parameterization (since mleGP does);
## R-only version

jmleGP.R <- function(gpi, N=100, drange=c(sqrt(.Machine$double.eps), 10),
  grange=c(sqrt(.Machine$double.eps), 1), dab=c(0,0), gab=c(0,0), maxit=100,
  mleGP=mleGP.R, verb=0)
  {
    ## sanity check N
    if(length(N) != 1 && N > 0)
      stop("N should be a positive scalar integer")
    m <- getmGP(gpi)
    glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen

    dmle <- matrix(NA, nrow=N, ncol=m)
    gmle <- matrix(NA, nrow=N, ncol=glen)
    dits <- dconv <- gits <- rep(NA, N)

    ## sanity check tmin and tmax
    if(length(drange) != 2) stop("drange should be a 2-vector for c(min,max)")
    if(length(grange) != 2) stop("grange should be a 2-vector for c(min,max)")

    ## loop over outer interations
    for(i in 1:N) {
      d <- mleGP(gpi, param="d", tmin=drange[1], tmax=drange[2],
                    ab=dab, maxit=maxit, verb=verb)
      dmle[i,] <- d$d; dits[i] <- d$its; dconv[i] <- d$conv
      g <- mleGP(gpi, param="g", tmin=grange[1], tmax=grange[2],
                    ab=gab, verb=verb)
      gmle[i,] <- g$g; gits[i] <- g$its
      if((gits[i] <= 2*glen && (dits[i] <= m+1 && dconv[i] == 0)) || dconv[i] > 1) break;
    }

    ## check if not converged
    if(i == N) warning("max outer its (N=", N, ") reached", sep="")
    else {
      dmle <- dmle[1:i,]; dits <- dits[1:i]; dconv <- dconv[1:i]
      gmle <- gmle[1:i,]; gits <- gits[1:i]
    }

    ## total iteration count
    totits <- sum(c(dits, gits), na.rm=TRUE)

    ## assemble return objects
    return(list(mle=data.frame(d=dmle[i,,drop=FALSE], g=gmle[i,,drop=FALSE], tot.its=totits,
      conv=dconv[i]), prog=data.frame(dmle=dmle, dits=dits, dconv=dconv, gmle=gmle,
      gits=gits)))
  }


## jmleGP
##
## interface to C-version for jmleGP;
## right now doesn't take an N argument -- the C-side hard-codes
## N=100

jmleGP <- function(gpi, drange=c(sqrt(.Machine$double.eps), 10),
   grange=c(sqrt(.Machine$double.eps), 1), dab=c(0,0), gab=c(0,0), maxit=100,
   verb=0)
  {
    ## sanity check tmin and tmax
    m <- getmGP(gpi)
    if(length(drange) != 2) stop("drange should be a two vector for c(dmin, dmax)")
    dmin <- rep(drange[1], m)
    dmax <- rep(drange[2], m)
    if(length(grange) != 2) stop("grange should be a 2-vector for c(gmin, gmax)")

    ## sanity check ab
    if(length(dab) != 2 || any(dab < 0)) stop("dab should be a positive 2-vector")
    if(length(gab) != 2 || any(gab < 0)) stop("gab should be a positive 2-vector")

    ## get number of nuggets
    glen <- .C("getglenGP_R", gpi=as.integer(gpi), glen=integer(1))$glen

    ## call the C-side function
    r <- .C("jmleGP_R",
            gpi = as.integer(gpi),
            maxit = as.integer(maxit),
            verb = as.integer(verb),
            dmin = as.double(dmin),
            dmax = as.double(dmax),
            grange = as.double(grange),
            dab = as.double(dab),
            gab = as.double(gab),
            d = double(m),
            g = double(glen),
            dits = integer(1),
            gits = integer(1),
            dconv = integer(1))

    return(data.frame(d=t(r$d), g=t(r$g), tot.its=r$dits+r$gits,
                      dits=r$dits, gits=r$gits, dconv=r$dconv))
  }


## predGP
##
## obtain the parameters to a multivariate-t
## distribution describing the predictive surface
## of the fitted GP model
## ' @export
predGP <- function(gpi, XX, ggi=NULL, lite=FALSE)
  {
    nn <- nrow(XX)

    if(is.null(ggi)) ggi <- rep(1, nn)
    else if(length(ggi) == 1) ggi <- rep(ggi, nn)
    if(length(ggi) != nn) stop("length(ggi) should be 1 or equal nrow(XX)")

    if(lite) {  ## lite means does not compute full Sigma, only diag
      out <- .C("predGP_R",
                gpi = as.integer(gpi),
                m = as.integer(ncol(XX)),
                nn = as.integer(nn),
                XX = as.double(t(XX)),
                ggi = as.integer(ggi-1),
                lite = as.integer(TRUE),
                mean = double(nn),
                s2 = double(nn),
                df = double(1),
                llik = double(1))

      ## coerce matrix output
      return(list(mean=out$mean, s2=out$s2, df=out$df, llik=out$llik))

    } else { ## compute full predictive covariance matrix

      out <- .C("predGP_R",
                gpi = as.integer(gpi),
                m = as.integer(ncol(XX)),
                nn = as.integer(nn),
                XX = as.double(t(XX)),
                ggi = as.integer(ggi-1),
                lite = as.integer(FALSE),
                mean = double(nn),
                Sigma = double(nn*nn),
                df = double(1),
                llik = double(1))

      ## coerce matrix output
      Sigma <- matrix(out$Sigma, ncol=nn)

      ## return parameterization
      return(list(mean=out$mean, Sigma=Sigma, df=out$df, llik=out$llik))
    }
  }


