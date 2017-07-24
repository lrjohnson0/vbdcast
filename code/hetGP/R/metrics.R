library(mvtnorm)

bins.pw <- 1:52
bins.pi <- bins.ti <-list()
bins.pi[["iquitos"]] <- c(0, seq(15,150,by=15), Inf)
bins.ti[["iquitos"]] <- c(0, seq(100,1000,by=100), Inf)
bins.pi[["sanjuan"]] <- c(0, seq(50,500,by=50), Inf)
bins.ti[["sanjuan"]] <- c(0, seq(1000,10000,by=1000), Inf)

calc.metrics <- function(tc, week, p, w, N=100000, 
	location=c("sanjuan", "iquitos"), alpha=0.05)
  {
    ## truth
	  pw <- which.max(tc)
    pi <- tc[pw]
    ti <- sum(tc)

    ## true bin
    location <- match.arg(location)
    pwu <- which(pw <= bins.pw)[1]
    piu <- which(pi < bins.pi[[location]])[1]
    bin.pi <- c(bins.pi[[location]][piu-1], bins.pi[[location]][piu])
    tiu <- which(ti < bins.ti[[location]])[1]
    bin.ti <- c(bins.ti[[location]][tiu-1], bins.ti[[location]][tiu])
	
    ## predictive sample with first weight 
    Y <- w[1] * (rmvt(N, p[[1]]$Sigma, df=p[[1]]$df) + matrix(rep(p[[1]]$mean, N), byrow=TRUE, ncol=52))
	
    ## accumulate for other weight
    if(length(w) > 1) {
      for(i in 2:length(w)) {
        if(w[i] <= 0) next;
        Y <- Y + w[i] * (rmvt(N, p[[i]]$Sigma, df=p[[i]]$df) + matrix(rep(p[[i]]$mean, N), byrow=TRUE, ncol=52))
      }
    } else if(w[1] != 1) stop("single weight not 1")

    ## for plotting mean and error bars
    Y <- f2i(Y)-1
    Ymean <- colMeans(Y)
    Yq1 <- apply(Y, 2, quantile, probs=alpha/2) 
    Yq2 <- apply(Y, 2, quantile, probs=1-alpha/2)

    ## empirical CDFs for PIT calculation before incorporating deterministic
    ## past weeks
    ticdf.nopast <- ecdf(rowSums(Y))
    pwcdf.nopast <- ecdf(apply(Y, 1, which.max))
    picdf.nopast <- ecdf(apply(Y, 1, max))

    ## augment with deterministic past weeks
    if(week > 0) {
      Yact <- matrix(rep(tc[1:week], N), byrow=TRUE, ncol=week)
      Y[,1:week] <- Yact
    }

    ## calculate distributions of incidence 
    tiY <- rowSums(Y)
    pwY <- apply(Y, 1, which.max)
    piY <- apply(Y, 1, max)

    ## empirical CDFs for PIT calculation
    ticdf <- ecdf(tiY)
    pwcdf <- ecdf(pwY)
    picdf <- ecdf(piY)

    ## point predictions via medians for MAE metric
    pti <- median(tiY)
    ppw <- median(pwY)
    ppi <- median(piY)

    ## summary of distributions
    tiQ <- quantile(tiY, c(alpha/2, 0.5, 1-alpha/2))
    piQ <- quantile(piY, c(alpha/2, 0.5, 1-alpha/2))
    pwQ <- quantile(pwY, c(alpha/2, 0.5, 1-alpha/2))

    ## augment with prior on each metric
    reps <- round(0.01*N*(52 - week - 4)/52)
    mtcw <- stcw <- 0
    if(week > 0) {
      mtcw <- max(tc[1:week])
      stcw <- sum(tc[1:week])
    }
    ti.reps <- bins.ti[[location]][bins.ti[[location]] > stcw]
    tiY <- c(tiY, rep(ti.reps, reps))
    pw.reps <- bins.pw[bins.pw > week]
    pwY <- c(pwY, rep(pw.reps, reps))
    pi.reps <- bins.pi[[location]][bins.pi[[location]] > mtcw]
    piY <- c(piY, rep(pi.reps, reps))

    ## calculate probabilities of true bins
    pw.p <- mean(pwY == pw)
    pi.p <- mean(piY <= bin.pi[2] & piY > bin.pi[1])
    ti.p <- mean(tiY <= bin.ti[2] & tiY > bin.ti[1])

    ## progress print
    cat("pw=", pw, ", ival=", bins.pw[pwu], ", p=", pw.p, "\n", sep="")
    cat("pi=", pi, ", ival=(", bin.pi[1], ", ", bin.pi[2], "), p=", pi.p, "\n", sep="")
    cat("ti=", ti, ", ival=(", bin.ti[1], ", ", bin.ti[2], "), p=", ti.p, "\n", sep="")

    return(list(tiQ=tiQ, piQ=piQ, pwQ=pwQ, ti=ti, pw=pw, pi=pi, ppi=ppi, ppw=ppw, pti=pti,
      pi.pit=picdf(pi), pw.pit=pwcdf(pw), ti.pit=ticdf(ti), 
      pi.pit.nopast=picdf.nopast(pi), pw.pit.nopast=pwcdf.nopast(pw), ti.pit.nopast=ticdf.nopast(ti), 
      pi.ls=log(pi.p), pw.ls=log(pw.p), ti.ls=log(ti.p), 
      Ymean=Ymean, Yq1=Yq1, Yq2=Yq2))
  }


plot.metrics <- function(mets, offset)
  {

  	## plot the metric medians and error bars
    points(offset+mets$pwQ[2], 0, pch=10, col=4)
    if(mets$pwQ[3] > mets$pwQ[1]) {
      arrows(offset+mets$pwQ[1], 0, offset+mets$pwQ[3], 0, angle=90, col=4, length=0.1)
      arrows(offset+mets$pwQ[3], 0, offset+mets$pwQ[1], 0, angle=90, col=4, length=0.1)
    }
    points(offset+mets$pwQ[2], mets$piQ[2], pch=10, col=4)
    if(mets$piQ[3] > mets$piQ[1]) {
      arrows(offset+mets$pwQ[2], mets$piQ[1], offset+mets$pwQ[2], mets$piQ[3], angle=90, col=4, length=0.1)
      arrows(offset+mets$pwQ[2], mets$piQ[3], offset+mets$pwQ[2], mets$piQ[1], angle=90, col=4, length=0.1)
    }
    points(offset+55, mets$tiQ[2]/52, pch=10, col=4)
    if(mets$tiQ[3] > mets$tiQ[1]) {
      arrows(offset+55, mets$tiQ[1]/52, offset+55, mets$tiQ[3]/52, angle=90, col=4, length=0.1)
      arrows(offset+55, mets$tiQ[3]/52, offset+55, mets$tiQ[1]/52, angle=90, col=4, length=0.1)
    }

    ## truth
    points(offset+55, mets$ti/52, pch=9, col=2)
    points(offset+mets$pw, mets$pi, pch=9, col=2)

    ## plot the full predictive distribution
    lines(offset+(1:52), mets$Ymean, col=4, lwd=2)
    lines(offset+(1:52), mets$Yq1, col=4, lty=2, lwd=2)
    lines(offset+(1:52), mets$Yq2, col=4, lty=2, lwd=2)
  }

save.quants <- function(mets, week, season, file="quants.csv")
  {
  	df <- data.frame(
  		season=season, week=week, 
  		pwq1=mets$pwQ[1], pwm=mets$pwQ[2], pwq2=mets$pwQ[3], pwtrue=mets$pw, pwpred=mets$ppw,
  		piq1=mets$piQ[1], pim=mets$piQ[2], piq2=mets$piQ[3], pitrue=mets$pi, pipred=mets$ppi,
  		tiq1=mets$tiQ[1], tim=mets$tiQ[2], tiq2=mets$tiQ[3], titrue=mets$ti, tipred=mets$pti)

  	append <- FALSE
  	if(file.exists(file)) append <- TRUE

  	 write.table(df, file=file, append=append, col.names=!append, row.names=FALSE, sep=",")
  }

save.pit <- function(mets, week, season, file="pits.csv")
  {
    df <- data.frame(
      season=season, week=week, 
      pw=mets$pw.pit, pi=mets$pi.pit, ti=mets$ti.pit,
      pw.nopast=mets$pw.pit.nopast, pi.nopast=mets$pi.pit.nopast, ti.nopast=mets$ti.pit.nopast)

    append <- FALSE
    if(file.exists(file)) append <- TRUE

     write.table(df, file=file, append=append, col.names=!append, row.names=FALSE, sep=",")
  }

save.logscore <- function(mets, week, season, location, file="newscores.csv", addAE=TRUE)
  {
  	data.type <- "test"
  	if(season %in% c("2005/2006", "2006/2007", "2007/2008", "2008/2009"))
  		data.type <- "train"

  	df <- data.frame(
  		team="new", metric.type="LS", forecast.week=week, season=season,
  		metric=c(mets$pw.ls, mets$pi.ls, mets$ti.ls), location=location, 
  		data.type=data.type, prediction.target=c("peakweek", "peakinc", "seasoninc"))

  	append <- FALSE
  	if(file.exists(file)) append <- TRUE

  	write.table(df, file=file, append=append, col.names=!append, row.names=FALSE, sep=",")

    if(addAE) {
      df <- data.frame(
        team="new", metric.type="AE", forecast.week=week, season=season,
        metric=c(abs(mets$pw - mets$ppw), abs(mets$pi - mets$ppi), abs(mets$ti - mets$pti)), 
        location=location, data.type=data.type, prediction.target=c("peakweek", "peakinc", "seasoninc"))

        write.table(df, file=file, append=TRUE, col.names=FALSE, row.names=FALSE, sep=",")
    }
  }

  plot.peakweek <- function(quants, location="San Juan", ttword=TRUE, legendloc="topright")
  	{
  		plot(quants$pwtrue, col=2, pch=18, ylim=range(c(quants$pwq1, quants$pwq2)),
  			ylab="peak week", xlab="observation weeks by year", xaxt="n",
        main=location)
  		points(quants$pwm, col=4)
  		segments(1:nrow(quants), y0=quants$pwq1, y1=quants$pwq2, col=4, lwd=2)

  		## work on x-axis
  		seasons <- unique(quants$season)
  		at <- seq(1, nrow(quants), by=52/4)
  		axis(1, at=at+6, labels=seasons, tick=FALSE)

  		## put train/test divide
  		w <- which(seasons == "2009/2010")
  		abline(v=at[w]-1/2, col=3)
  		if(ttword) text(c(at[w]-3.5, at[w]+2.25), c(14,14), c("train", "test"), col=3, cex=2)

  		## put "when" on x-axis
  		pwtrue <- quants$pwtrue[at]
  		points(pwtrue/4 + seq(0, (52/4)*(length(at)-1), by=13), pwtrue, col=2, pch=4, cex=2, lwd=2)

  		legend(legendloc, pch=c(1,NA,18,4), col=c(4,4,2,2), pt.lwd=c(1,0,1,2), lty=c(0,1,0,0), lwd=c(0,2,0,0),
  			c("predicted peak week", "peak week 95% interval", "actual peak week", "actual peak week on x-axis"))
  	}


  plot.peakincidence <- function(quants, location="San Juan", ttword=TRUE, legendloc="topright")
  	{
  		plot(quants$pitrue, col=2, pch=18, ylim=range(c(quants$piq1, quants$piq2)),
  			ylab="peak incidence", xlab="observation weeks by year", xaxt="n",
        main=location)
  		points(quants$pim, col=4)
  		segments(1:nrow(quants), y0=quants$piq1, y1=quants$piq2, col=4, lwd=2)

  		## work on x-axis
  		seasons <- unique(quants$season)
  		at <- seq(1, nrow(quants), by=52/4)
  		axis(1, at=at+6, labels=seasons, tick=FALSE)

  		## put train/test divide
  		w <- which(seasons == "2009/2010")
  		abline(v=at[w]-1/2, col=3)
  		if(ttword) text(c(at[w]-3.5, at[w]+2.25), c(325,325), c("train", "test"), col=3, cex=2)

  		## put "when" on x-axis
  		pwtrue <- quants$pwtrue[at]
  		pitrue <- quants$pitrue[at]
  		points(pwtrue/4 + seq(0, (52/4)*(length(at)-1), by=13), pitrue, col=2, pch=4, cex=2, lwd=2)

  		legend(legendloc, pch=c(1,NA,18,4), col=c(4,4,2,2), pt.lwd=c(1,0,1,2), 
  			lty=c(0,1,0,0), lwd=c(0,2,0,0),
  			c("predicted peak incidence", "peak incidence 95% interval", 
  				"actual peak incidence", "actual peak incidence on x-axis"))
  	}

plot.totalincidence <- function(quants, location="San Juan", ttword=TRUE, legendloc="topright")
  	{
  		plot(quants$titrue, col=2, pch=18, ylim=range(c(quants$tiq1, quants$tiq2)),
  			ylab="total incidence", xlab="observation weeks by year", xaxt="n",
        main=location)
  		points(quants$tim, col=4)
  		segments(1:nrow(quants), y0=quants$tiq1, y1=quants$tiq2, col=4, lwd=2)

  		## work on x-axis
  		seasons <- unique(quants$season)
  		at <- seq(1, nrow(quants), by=52/4)
  		axis(1, at=at+6, labels=seasons, tick=FALSE)

  		## put train/test divide
  		w <- which(seasons == "2009/2010")
  		abline(v=at[w]-1/2, col=3)
  		if(ttword) text(c(at[w]-3.5, at[w]+2.25), c(7500,7500), c("train", "test"), col=3, cex=2)

  		legend(legendloc, pch=c(1,NA,18), col=c(4,4,2), pt.lwd=c(1,0,1), 
  			lty=c(0,1,0), lwd=c(0,2,0),
  			c("predicted total incidence", "total incidence 95% interval", 
  				"actual total incidence"))
  	}