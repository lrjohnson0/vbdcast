rm(list=ls())

## data <- read.csv("../../../data/dengue/iquitos_training.csv")
data <- read.csv("../../../data/dengue/iquitos_testing.csv")

tr.seas <- c("2008/2009", "2009/2010", "2010/2011", "2011/2012")
tr.seas <- c("2004/2005", "2005/2006", "2006/2007", "2007/2008", tr.seas)
te.seas <- c("2009/2010", "2010/2011", "2011/2012", "2012/2013")
te.seas <- c("2005/2006", "2006/2007", "2007/2008", "2008/2009", te.seas)

breaks <- c(25,10)
location <- "iquitos"
ylim <- c(0,100)
stepthrough <- FALSE

source("gp_dengue.R")

pdf("iquitos_peak_incidence.pdf", height=6, width=14)
quants <- read.csv("quants_iquitos.csv")
plot.peakincidence(quants, location="Iquitos")
dev.off()
pdf("iquitos_peak_week.pdf", height=6, width=14)
plot.peakweek(quants, location="Iquitos", ttword=FALSE, legendloc="bottomright")
dev.off()
pdf("iquitos_total_incidence.pdf", height=6, width=14)
plot.totalincidence(quants, location="Iquitos", ttword=FALSE)
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

pdf("iquitos_true_offset_pi.pdf", width=5, height=5)
plot(offset[1:5,1], pw[1:5], xlab="starting level", ylab="peak week", 
 	main="Iquitos", xlim=range(offset), ylim=range(pw), col=ratcols[rating[1:5,1]+2])
text(offset[6:13,1], pw[6:13], ((6:13)-5)+4, col=ratcols[rating[6:13,1]+2])
dev.off()
pdf("iquitos_true_offset_pw.pdf", width=5, height=5)
plot(offset[1:5,1], pi[1:5], xlab="starting level", ylab="peak incidence", 
 	main="Iquitos", xlim=range(offset), ylim=range(pi), col=ratcols[rating[1:5,1]+2])
text(offset[6:13,1], pi[6:13], ((6:13)-5)+4, col=ratcols[rating[6:13,2]+2])
dev.off()
pdf("iquitos_true_pw_pi.pdf", width=5, height=5)
plot(pw[1:5], pi[1:5], xlab="peak week", ylab="peak incidence",
 	main="Iquitos", xlim=range(pw), ylim=range(pi), col=ratcols[rating[1:5,1]+2])
text(pw[6:13], pi[6:13], ((6:13)-5)+4, col=ratcols[rating[6:13,1]+2])
dev.off()


## histograms of PIT
pit <- read.csv("pit_iquitos.csv")
pdf("pit_iquitos.pdf", width=5, height=14)
par(mfrow=c(3,1), mar=c(2,5,1,1))
hist(pit$pi.nopast, xlab="", ylab="peak incidence", cex.lab=2, main="", ylim=c(0,60))
hist(pit$pw.nopast, xlab="", ylab="peak week", cex.lab=2, main="", ylim=c(0,60))
hist(pit$ti.nopast, xlab="", ylab="total incidence", cex.lab=2, main="", ylim=c(0,60))
dev.off()