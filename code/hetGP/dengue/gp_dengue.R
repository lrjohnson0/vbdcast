dyn.load("../../src/hetGP.so")
source("../../R/gp.R")
source("../../R/transform.R")
source("../../R/metrics.R")
source("../../R/prior.R")

## create score file names and remove old one
scorefile <- paste("newscores_", location, ".csv", sep="")
unlink(scorefile)

## create quants file name and remove old one
quantsfile <- paste("quants_", location, ".csv", sep="")
unlink(quantsfile)

## create PIT file name and remove old one
pitfile <- paste("pit_", location, ".csv", sep="")
unlink(pitfile)

ratcols <- c("forestgreen", "purple", "red")

## initialize MAP hyperparameters
dmle <- gmle <- NULL

## init output object
out.obj <- list()

## loop over training seasons
for(si in 1:length(tr.seas)) {

	## extract out training and testing data
	nstart <-max(which(data$season==tr.seas[si]))
	seasons <- unique(data$season[1:nstart])

	## constructing the response matrix and predictor variables
	rating <- offset <- lY <- Y <- matrix(NA, nrow=length(seasons), ncol=52)
	for(i in 1:length(seasons)) {
		season <- which(data$season == seasons[i])
		Y[i,] <- data$total_cases[season]
		lY[i,] <- f2(Y[i,] + 1)
		offset[i,] <- mean(lY[i,1:4])
		if(sum(Y[i,] > breaks[1]) >= 2) rating[i,] <- 1
		else if(max(Y[i,]) <= breaks[2]) rating[i,] <- -1
		else rating[i,] <- 0
	}

	## normalize offsets
	o.norm <- min(offset)
	if(o.norm > 0) offset <- offset/o.norm - 1

	## fit offset-severity relationship
	odf <- data.frame(x=offset[,1], y=apply(lY, 1, max))
	ofit <- lm(y ~ x, data=odf)

	pdf("prior.pdf", width=5, height=5)
	plot(odf, xlab="x3: starting level", ylab="peak incidence")
	abline(h=f2(breaks+1), col=2, lty=2, lwd=2)
	text(x=rep(6,3), y=c(3.25,6.5,11), labels=c("x4=-1", "x4=0", "x4=+1"), col=4)
	stop()

	## converting the year-matrices above into flat predictors
	xw <- rep((1:52)/52, nrow(lY))
	xo <- as.numeric(t(offset))
	xrating <- as.numeric(t(rating))
	gi <- xrating+2
	xs <- rep(0.5*(1-cos((0:51)*2*pi/52)), length(seasons))

	## assemble the actual design matrix and respone for the GP regression
	X <- cbind(xw, xo, xrating, xs)
	## X <- cbind(xw, xo, xs)
	Z <- as.numeric(t(lY))

	## default priors
	da <- darg(list(mle=TRUE), X=X, samp.size=nrow(X))
	ga <- garg(list(mle=TRUE), y=Z)
	glen <- length(unique(gi))

	## initial values, taken from mle's found in previous iteration;
	## not sure if this is a good idea (multi-modal)
	if(is.null(dmle)) dstart <- da$start
	else dstart <- dmle
	if(is.null(gmle)) gstart <- rep(10*ga$start, glen) 
	else gstart <- gmle

	## initialize a GP
	gpi <- newGP(X=X, Z=Z, d=dstart, g=gstart, gi=gi) ##, dK=TRUE)

	## progress meter
	cat("performing MAP calculations for season ", te.seas[si], 
		"\n", sep="")

	## maximum a' posteriori calculations for hyperparameters
	if(is.null(dmle) || is.null(gmle)) {
		mle <- mleGP(gpi, ab=c(da$ab, ga$ab), tmin=c(da$min, ga$min), 
					 tmax=c(da$max, ga$max), param="both")
 		mle <- mle$theta
		print(mle)
	}

	## the first call does not converge so we do another one
	mle <- mleGP(gpi, ab=c(da$ab, ga$ab), tmin=c(da$min, ga$min), 
				 tmax=c(da$max, ga$max), param="both")
 	mle <- mle$theta
	print(mle)

	## store for easier access later
	dmle <- mle[1:ncol(X)] 
	gmle <- mle[(ncol(X)+1):(ncol(X)+glen)]

	## initial plot of training data
	skip <- 0 ## prod(dim(Y)) - 52*3
	ptrain <- predGP(gpi, XX=X, ggi=gi)
	alpha <- 0.05
	q1train <- qt(alpha/2, df=ptrain$df)*sqrt(diag(ptrain$Sigma)) + ptrain$mean
	q2train <- qt(1-(alpha/2), df=ptrain$df)*sqrt(diag(ptrain$Sigma)) + ptrain$mean
	plot(as.numeric(t(Y)), xlim=c(skip,nrow(Y)*ncol(Y)+52), col=ratcols[gi], 
		ylim=ylim, main=paste(location, "incidence"), xaxt="n", ylab="incidence")
	at <- seq(1, prod(dim(Y))+52, by=52)
  	axis(1, at=at+26, labels=c(as.character(seasons), te.seas[si]), tick=FALSE)
	lines(f2i(ptrain$mean)-1)
	lines(f2i(q1train)-1, lty=2)
	lines(f2i(q2train)-1, lty=2)
	abline(v=seq(1,(nrow(Y)*ncol(Y)+52), by=52), lty=2, col="gray")

	## look at level of heteroskedasticity in training data via residuals
	resid <- as.numeric(t(lY)) - ptrain$mean
	pdf(paste("resid_", location, ".pdf", sep=""), width=5)
	resid.list <- list(low=resid[gi == 1], mid=resid[gi == 2], high=resid[gi == 3])
	boxplot(resid.list, col=ratcols, xlab="severity", ylab="residuals", main=location)
	dev.off()

	## get predictive equations for each GP
	season <- which(data$season==te.seas[si])
	x.new <- (nrow(Y)*ncol(Y)) + c(1:52)
	y.new <- data$total_cases[season]
	ly.new <- f2(y.new + 1)
	o.new <- mean(Z[(length(Z)-3):length(Z)]) ## ly.new[1]
	if(o.norm > 0) o.new <- o.new/o.norm - 1
	ps <- predict(ofit, newdata=data.frame(x=o.new))
	ps <- f2i(ps)-1
	xo3s <- unique(xo[gi == 3])
	if(length(xo3s) > 0) o.new <- xo3s[which.min(abs(o.new - xo3s))]
	segments(length(Y)-3, f2i(o.new)-1, length(Y), f2i(o.new)-1, lwd=3)
	s.new <- 0.5*(1-cos((0:51)*2*pi/52))

	## a function which evaluates how good a particular "fit" is for
	## a certaint set of the the latent variable pair (latent for the
	## input xrating, and gi for the nugget indicator)
	library(mvtnorm)
	obj <- function(latent, X, Z, back, d, g, gi)
	{
		ii <- (nrow(X)-back+1):nrow(X)
		X[ii,3] <- latent
		gpi <- newGP(X=X, Z=Z, d=d, g=g, gi=gi)
		p <- predGP(gpi, X[ii,], ggi=gi[ii])  ## considers only new data
		ll <- dmvnorm(Z[ii], p$mean, p$Sigma, log=TRUE)
		deleteGP(gpi)
		return(-ll)
	}

	## prime the new year
	s.prime <- 0.5*(1-cos((51:(51-4))*2*pi/52))
	olatents <- c(-1, 0, 1) ## c(-0.5, 0, 0.5)
	nas <- rep(NA, length(olatents))
	GPs <- data.frame(lprob=nas, olatent=nas)
	deleteGP(gpi)

	cat("forecasting for season ", te.seas[si], ", aug ", 0, 
		", o.new=", o.new, "\n", sep="")

	## get truth for benchmarking
	tc <- data$total_cases[data$season == te.seas[si]]

	## for accumulating predictions
	pnew <- list()

	## the priming is applied for each latent setting
	for(oi in 1:length(olatents)) {

		## skip if none at this noise have been observed
		if(!(any(gi == oi))) {
			GPs[oi,] <- c(Inf, olatents[oi])
			pnew[[oi]] <- NULL
			next;
		}

		X.prime <- cbind((0:(-4))/52, o.new, olatents[oi], s.prime)
		Z.prime <- Z[length(Z):(length(Z)-4)]
		gi.prime <- rep(oi, length(Z.prime))

		## build GP and evaluate predictive probability for the prime data
		gpi <- newGP(X=rbind(X, X.prime), Z=c(Z, Z.prime), d=dmle, g=gmle,
					 gi=c(gi, gi.prime)) ##, dK=TRUE)
		out <- obj(olatents[oi], rbind(X, X.prime), c(Z, Z.prime), 
				   back=5, d=dmle, g=gmle, gi=c(gi, gi.prime))
		GPs[oi,] <- c(out, olatents[oi])

		## forecast ahead uusing the current latent(s)
		XX <- cbind((1:52)/52, o.new, olatents[oi], s.new)
		points(x.new, y.new)
		pnew[[oi]] <- predGP(gpi, XX=XX, ggi=oi)

		## clean up
		deleteGP(gpi)
	}

	## calculate the weights with the prior
	if(length(pnew) == 3) {
		if(ps > breaks[1]) w.prior <- c(0.25, 0.25, 0.5)
		else if(ps > breaks[2]) w.prior <- c(0.25, 0.5, 0.25)
		else w.prior <- c(0.5, 0.25, 0.25)
	} else if(length(pnew) == 2) {
		if(ps > breaks[1]) w.prior <- c(0.4, 0.6)
		else wprior <- c(0.6, 0.4)
	} else stop("length(pnew) = 1 not supported")

	w <- exp(-GPs$lprob - max(-GPs$lprob))
	dlen <- 1
	w <- (w/sum(w))*(dlen/52) + w.prior*((52-dlen))/52
	cat("weights (latents) on the three processes: \n")
	for(k in 1:3) cat(signif(w[k], 4), " (", signif(GPs$olatent[k], 4), ") ", sep="")
	cat("\n")	

	## calculate metrics and plot
	mets <- calc.metrics(tc, 0, pnew, w, location=location, alpha=alpha)
	plot.metrics(mets, nrow(Y)*ncol(Y))
	save.logscore(mets, 0, te.seas[si], location=location, file=scorefile)
	save.quants(mets, 0, te.seas[si], file=quantsfile)
	save.pit(mets, 0, te.seas[si], file=pitfile)

	## wait for RETURN?
	if(stepthrough) readline("press RETURN: ")

	o.new.actual <- mean(ly.new[1:4]) ## ly.new[1]
	if(length(xo3s) > 0) o.new.x3 <- xo3s[which.min(abs(o.new.actual - xo3s))]
	else o.new.x3 <- o.new.actual

	## now loop over the forecasting weeks
	augs <- c(0, seq(4,48,by=4))
	for(i in 2:length(augs)) {

		o.new <- o.new.actual*augs[i]/max(augs) + o.new.x3*(1 - augs[i]/max(augs))

		## progress meter
		cat("forecasting for season ", te.seas[si], ", aug ", augs[i], 
			", o.new=", o.new, "\n", sep="")

		## for accumulating predictions
		pnew <- list()

		## each pairing of latent variables
		for(oi in 1:length(olatents)) {

			## skip if none at this noise have been observed
			if(!(any(gi == oi))) {
				GPs[oi,] <- c(Inf, olatents[oi])
				pnew[[oi]] <- NULL
				next;
			}

			## copied from above, priming data using chosen latents
			X.prime <- cbind((0:(-4))/52, o.new, olatents[oi], s.prime)
			Z.prime <- Z[length(Z):(length(Z)-4)]
			gi.prime <- rep(oi, length(Z.prime))

			## optimize the latent variable
			olatent <- optimize(obj, 
				interval=c(GPs$olatent[oi]-0.25,GPs$olatent[oi]+0.25), 
				X=rbind(X, X.prime, XX[1:augs[i],]), 
				Z=c(Z, Z.prime, ly.new[1:augs[i]]), 
				back=5+augs[i], d=dmle, g=gmle, 
				gi=c(gi, gi.prime, rep(oi, augs[i])))$minimum
			## olatent <- olatents[oi]

			## (re-) initialize a gp under the chosen latent setting(s)
			XX[1:52,3] <- olatent
			X.prime[,3] <- olatent
			gpi <- newGP(X=rbind(X, X.prime, XX[1:augs[i],]), 
						 Z=c(Z, Z.prime, ly.new[1:augs[i]]), 
						 d=dmle, g=gmle, 
						 gi=c(gi, gi.prime, rep(oi, augs[i]))) 
			out <- obj(olatent, rbind(X, X.prime, XX[1:augs[i],]), 
					   Z=c(Z, Z.prime, ly.new[1:augs[i]]), 
					   back=5+augs[i], d=dmle, g=gmle, 
					   gi=c(gi, gi.prime, rep(oi, augs[i])))
			GPs[oi,] <- c(out, olatent)

			## (re-) initialize the plot of training and forecasts, but only
			## when we are in the first iteration of the inner latent loop
			if(oi == 1) {
				plot(as.numeric(t(Y)), xlim=c(skip,nrow(Y)*ncol(Y)+52), 
					col=ratcols[xrating+2], ylim=ylim, main=paste(location, "incidence"),
					xaxt="n")
				at <- seq(1, prod(dim(Y))+52, by=52)
  				axis(1, at=at+26, labels=c(as.character(seasons), te.seas[si]), tick=FALSE)
				lines(f2i(ptrain$mean)-1)
				lines(f2i(q1train)-1, lty=2)
				lines(f2i(q2train)-1, lty=2)
				points(x.new, y.new)
				points(x.new[(augs[1]+1):augs[i]], 
					   y.new[(augs[1]+1):augs[i]], col=1, pch=19, cex=0.5)
				abline(v=seq(1,(nrow(Y)*ncol(Y)+52), by=52), lty=2, 
					   col="gray")
				segments(x.new[1], f2i(o.new)-1, x.new[4], f2i(o.new)-1, 
						 lwd=3)
			}

			## calculate forecasts for the current set of latents
			pnew[[oi]] <- predGP(gpi, XX=XX, ggi=oi)
		
			## clean up
			deleteGP(gpi)
		}

		## calculate weights with prior
		w <- exp(-GPs$lprob - max(-GPs$lprob))
		dlen <- augs[i]+1
		w <- (w/sum(w))*(dlen/52) + w.prior*((52-dlen))/52
		cat("weights (latents) on the three processes: \n")
		for(k in 1:3) cat(signif(w[k], 4), " (", signif(GPs$olatent[k], 4), ") ", sep="")
		cat("\n")	
		## print(signif(w, 4))

		## calculate metrics
		mets <- calc.metrics(tc, augs[i], pnew, w, location=location, alpha=alpha)
		plot.metrics(mets, nrow(Y)*ncol(Y))
		save.logscore(mets, augs[i], te.seas[si], location=location, file=scorefile)
		save.quants(mets, augs[i], te.seas[si], file=quantsfile)
		save.pit(mets, augs[i], te.seas[si], file=pitfile)

		## wait for RETURN?
		if(stepthrough) readline("press RETURN: ")
	}
}
