rm(list=ls())

## data <- read.csv("../../../data/sanjuan_training.csv")
data <- read.csv("../../../data/sanjuan_testing.csv")

tr.seas <- c("2008/2009", "2009/2010", "2010/2011", "2011/2012")
tr.seas <- c("2004/2005", "2005/2006", "2006/2007", "2007/2008", tr.seas)
te.seas <- c("2009/2010", "2010/2011", "2011/2012", "2012/2013")
te.seas <- c("2005/2006", "2006/2007", "2007/2008", "2008/2009", te.seas)

breaks <- c(100, 25)
location <- "sanjuan"
ylim <- c(0,300)
stepthrough <- FALSE

source("gp_dengue.R")

pdf("sanjuan_peak_incidence.pdf", height=6, width=14)
quants <- read.csv("quants_sanjuan.csv")
plot.peakincidence(quants, location="San Juan")
dev.off()
pdf("sanjuan_peak_week.pdf", height=6, width=14)
plot.peakweek(quants, location="San Juan", ttword=FALSE)
dev.off()
pdf("sanjuan_total_incidence.pdf", height=6, width=14)
plot.totalincidence(quants, location="San Juan", ttword=FALSE)
dev.off()

## for visualing the truth (plots in appendix)

## extract out training and testing data
seasons <- unique(data$season)

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

pi <- apply(Y, 1, max)
pw <- apply(Y, 1, which.max)

ratcols <- c("forestgreen", "purple", "red")

pdf("sanjuan_true_offset_pi.pdf", width=5, height=5)
plot(offset[1:15,1], pw[1:15], xlab="starting level", ylab="peak week", 
 	main="San Juan", xlim=c(0,7), col=ratcols[rating[1:15,1]+2])
text(offset[16:23,1], pw[16:23], ((16:23)-15)+4, col=ratcols[rating[16:23,1]+2])
dev.off()
pdf("sanjuan_true_offset_pw.pdf", width=5, height=5)
plot(offset[1:15,1], pi[1:15], xlab="starting level", ylab="peak incidence", 
 	main="San Juan", xlim=c(0,7), col=ratcols[rating[1:15,1]+2])
text(offset[16:23,1], pi[16:23], ((16:23)-15)+4, col=ratcols[rating[16:23,2]+2])
dev.off()
pdf("sanjuan_true_pw_pi.pdf", width=5, height=5)
plot(pw[1:15], pi[1:15], xlim=c(14,45), xlab="peak week", ylab="peak incidence",
 	main="San Juan", col=ratcols[rating[1:15,1]+2])
text(pw[16:23], pi[16:23], ((16:23)-15)+4, col=ratcols[rating[16:23,1]+2])
dev.off()


## histograms of PIT
## pdf("pit_sanjuan.pdf", width=5, height=14)
pit <- read.csv("pit_sanjuan.csv")
pdf("pit_sanjuan.pdf", width=14, height=5)
## par(mfrow=c(3,1), mar=c(2,5,1,1))
pit <- read.csv("pit_sanjuan.csv")
par(mfrow=c(1,3), mar=c(2,5,1,1))
hist(pit$pi.nopast, xlab="", ylab="peak incidence - hetGP", cex.lab=2, main="", ylim=c(0,50))
hist(pit$pw.nopast, xlab="", ylab="peak week - hetGP", cex.lab=2, main="", ylim=c(0,50))
hist(pit$ti.nopast, xlab="", ylab="total incidence - hetGP", cex.lab=2, main="", ylim=c(0,50))
dev.off()
