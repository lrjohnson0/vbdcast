## The code in this file converts draws from the hybrid GP/GLM
## draws from the predictive distribution of each target, and
## calculates a probability integral transform (PIT) measuring
## goodness of fit to the three (true) targets).

## First for San Juan

## combining training and testing, treatin both as out-of-sample
load("~/Dropbox/dengue_challenge/src/hybrid_dengue_sanjuan.RData")
obj <- out.obj
load("~/Dropbox/dengue_challenge/src/hybrid_dengue_sanjuan_testing.RData")
obj <- c(obj, out.obj)

## build up PIT through empirical cdf
pi.pit <- pw.pit <- ti.pit <- c()
for(i in 1:length(obj)) {

	pi <- max(obj[[i]]$ytrue)
	pw <- (which.max(obj[[i]]$ytrue))[1]
	ti <- sum(obj[[i]]$ytrue)

	picdf <- ecdf(apply(obj[[i]]$tc, 1, max))
	pwcdf <- ecdf(apply(obj[[i]]$tc, 1, which.max))
	ticdf <- ecdf(rowSums(obj[[i]]$tc))

	pi.pit <- c(pi.pit, picdf(pi))
	pw.pit <- c(pw.pit, pwcdf(pw))
	ti.pit <- c(ti.pit, ticdf(ti))

}

## output to sanjan file
pdf("pit_sanjuan_hybrid.pdf", width=5, height=14)
par(mfrow=c(3,1), mar=c(2,5,1,1))
hist(pi.pit, xlab="", ylab="peak incidence", cex.lab=2, main="", ylim=c(0,50))
hist(pw.pit, xlab="", ylab="peak week", cex.lab=2, main="", ylim=c(0,50))
hist(ti.pit, xlab="", ylab="total incidence", cex.lab=2, main="", ylim=c(0,50))
dev.off()

## Then for Iquitos

## combining training and testing, treating both as out-of-sample
load("~/Dropbox/dengue_challenge/src/hybrid_dengue_iquitos.RData")
obj <- out.obj
load("~/Dropbox/dengue_challenge/src/hybrid_dengue_iquitos_testing.RData")
obj <- c(obj, out.obj)

## build up PIT through empirical cdf
pi.pit <- pw.pit <- ti.pit <- c()
for(i in 1:length(obj)) {

	pi <- max(obj[[i]]$ytrue)
	pw <- (which.max(obj[[i]]$ytrue))[1]
	ti <- sum(obj[[i]]$ytrue)

	picdf <- ecdf(apply(obj[[i]]$tc, 1, max))
	pwcdf <- ecdf(apply(obj[[i]]$tc, 1, which.max))
	ticdf <- ecdf(rowSums(obj[[i]]$tc))

	pi.pit <- c(pi.pit, picdf(pi))
	pw.pit <- c(pw.pit, pwcdf(pw))
	ti.pit <- c(ti.pit, ticdf(ti))

}

## output to iquitos file
pdf("pit_iquitos_hybrid.pdf", width=5, height=14)
par(mfrow=c(3,1), mar=c(2,5,1,1))
hist(pi.pit, xlab="", ylab="peak incidence", cex.lab=2, main="", ylim=c(0,60))
hist(pw.pit, xlab="", ylab="peak week", cex.lab=2, main="", ylim=c(0,60))
hist(ti.pit, xlab="", ylab="total incidence", cex.lab=2, main="", ylim=c(0,60))
dev.off()