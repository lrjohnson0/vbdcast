## choose to collect results for sanjuan or iquitos
loc <- "sanjuan" ## "iquitos"
if(loc == "iquitos") { location <- "Iquitos"
} else if(loc == "sanjuan") { location <- "San Juan"
} else stop("bad loc")

## Read in the data
res <- read.csv("../csvouts/dengue_forecast_scores.csv", header=T)
res.new <- read.csv(paste("../csvouts/newscores_", loc, ".csv", sep=""), header=T)
res <- rbind(res, res.new)
res.new <- read.csv(paste("../csvouts/glm_newscores_", loc, ".csv", sep=""), header=T)
res <- rbind(res, res.new)

## set up utility variables for later
targets <-levels(res$prediction.target)
names <- targets
names[names == "peakinc"] <- "Peak Incidence"
names[names == "peakweek"] <- "Peak Week"
names[names == "seasoninc"] <- "Total Season Incidence"
seasons <- levels(res$season)
weeks <- unique(res$forecast.week)
teams <- levels(res$team)

## Log Score (LS)

### The competition was mainly evaluated on mean log-score for the
### testing data, so I pull that out here:
ls.res <- res[which(res$metric.type=="LS"),]
ls.res <- ls.res[which(ls.res$data.type=="test"),]
### Then parse out the data by location
ls.res <- ls.res[which(ls.res$location==loc),]
ls.dat <- list()

for(k in 1:length(targets)) {

    ##pick out the target
    res.temp<-ls.res[which(ls.res$prediction.target==targets[k]),]

    ## build a data frame to save the weekly results
    temp<-data.frame(matrix(NA, ncol=length(teams)+1, nrow=length(weeks)))
    names(temp)<-c("week", teams)
    temp$week <- weeks
                            
    for(j in 1:length(teams)) {
        
        newdat<-subset(res.temp, team==teams[j], select=c(forecast.week, metric))    
        for(i in 1:length(weeks))
            temp[i,j+1]<-mean(newdat$metric[which(newdat$forecast.week==weeks[i])], na.rm=TRUE)
    }

    ls.dat[[k]]<-list(location=location, metric="LS", data.type="test",
                      target=targets[k], name=names[k], results=temp)
}

## formerly _results rather than _ls_results
save(ls.dat, file=paste(location, "_ls_results.RData", sep="")) 

## make LS plots by week for each comparator
## par(mfrow=c(3,1), bty="n")
pdf(paste(loc, "_results.pdf", sep=""))
for(j in 1:3){
    plot(ls.dat[[j]]$results$week, ls.dat[[j]]$results[,2], ylim=c(-6, 0), type="l",
         xlab="weeks", ylab=paste(location, "log score"), main=ls.dat[[j]]$name)
    for(i in 1:length(teams)) lines(ls.dat[[j]]$results$week, ls.dat[[j]]$results[,i+1])
    lines(ls.dat[[j]]$results$week, ls.dat[[j]]$results$E, col=2, lwd=3)
    lines(ls.dat[[j]]$results$week, ls.dat[[j]]$results$new, col=3, lwd=3)
    lines(ls.dat[[j]]$results$week, ls.dat[[j]]$results$glmnew, col=4, lwd=3)
    lines(ls.dat[[j]]$results$week, ls.dat[[j]]$results$baseline, col=5, lwd=3)

    if(j == 2) 
        legend("topleft", c("comparators", "hybrid", "hetGP", "GLM", "baseline"), col=1:5, lwd=c(1,3,3,3,3),
            bty="n")
}
dev.off()

## table of ranks of LS
ls.rankstab <- cbind(peakinc=rowMeans(apply(ls.dat[[1]]$results[,-1], 1, rank)),
    peakweek=rowMeans(apply(ls.dat[[2]]$results[,-1], 1, rank)),
    seasoninc=rowMeans(apply(ls.dat[[3]]$results[,-1], 1, rank)))
ls.rankstab <- ls.rankstab[c(1,2,4,5,7,8,9,10,11,12,13,14,15,16,17,6,18,19,3),]

library(xtable)
xtable(ls.rankstab)

## table of averages of LS
ls.tab <- cbind(peakinc=colMeans(ls.dat[[1]]$results[,-1]),
    peakweek=colMeans(ls.dat[[2]]$results[,-1]),
    seasoninc=colMeans(ls.dat[[3]]$results[,-1]))
ls.tab <- ls.tab[c(1,2,4,5,7,8,9,10,11,12,13,14,15,16,17,6,18,19,3),]

## Mean Absolute Error

### Similar loops for average error for revised version of the paper,
### again for testing set only
ae.res <- res[which(res$metric.type=="AE"),]
ae.res <- ae.res[which(ae.res$data.type=="test"),]
### Then parse out the data by location
ae.res <- ae.res[which(ae.res$location==loc),]
ae.dat <- list()

for(k in 1:length(targets)) {

    ##pick out the target
    res.temp<-ae.res[which(ae.res$prediction.target==targets[k]),]

    ## build a data frame to save the weekly results
    temp<-data.frame(matrix(NA, ncol=length(teams)+1, nrow=length(weeks)))
    names(temp)<-c("week", teams)
    temp$week<-weeks
                            
    for(j in 1:length(teams)) {
        
        newdat<-subset(res.temp, team==teams[j], select=c(forecast.week, metric))
        for(i in 1:length(weeks))
            temp[i,j+1]<-mean(newdat$metric[which(newdat$forecast.week==weeks[i])], na.rm=TRUE)
    }

    ae.dat[[k]]<-list(location=location, metric="AE", data.type="test",
                      target=targets[k], name=names[k], results=temp)
}

save(ae.dat, file=paste(location, "_ae_results.RData", sep=""))

## make MAE plots by week for each comparator, but these look very similar to
## the LS plots above, so we will not include them in the paper
## par(mfrow=c(3,1), bty="n")
for(j in 1:3){
    plot(ae.dat[[j]]$results$week, log(ae.dat[[j]]$results[,2]+1), ylim=range(log(ae.dat[[j]]$results+1), na.rm=TRUE), type="l",
         xlab="weeks", ylab=paste(location, "mean absolute error"), main=ae.dat[[j]]$name)
    for(i in 1:length(teams)) lines(ae.dat[[j]]$results$week, log(ae.dat[[j]]$results[,i+1]+1))
    lines(ae.dat[[j]]$results$week, log(ae.dat[[j]]$results$E+1), col=2, lwd=3)
    lines(ae.dat[[j]]$results$week, log(ae.dat[[j]]$results$new+1), col=3, lwd=3)
    lines(ae.dat[[j]]$results$week, log(ae.dat[[j]]$results$glmnew+1), col=4, lwd=3)
    lines(ae.dat[[j]]$results$week, log(ae.dat[[j]]$results$baseline+1), col=5, lwd=3)

    if(j == 2) 
        legend("topleft", c("comparators", "hybrid", "hetGP", "GLM", "baseline"), col=1:5, lwd=c(1,3,3,3,3),
            bty="n")
}

## table of averages of AEs (MAE)
ae.tab <- cbind(peakinc=colMeans(ae.dat[[1]]$results[,-1], na.rm=TRUE),
    peakweek=colMeans(ae.dat[[2]]$results[,-1], na.rm=TRUE),
    seasoninc=colMeans(ae.dat[[3]]$results[,-1], na.rm=TRUE))
ae.tab <- ae.tab[c(1,2,4,5,7,8,9,10,11,12,13,14,15,16,17,6,18,19,3),]

## augment for yamana data; not available for Iquitos
if(loc == "sanjuan") {

    ## MAE comparison to yamana, discarding training proportion
    df <- read.csv("../csvouts/quants_sanjuan.csv")[-(1:52),]
    df.y <- read.csv("../../data/yamana/sanjuan.csv")[-(1:52),]

    pi.mae <- c(hetGP=mean(abs(df$pitrue - df$pipred)),
        F1=mean(abs(df$pitrue - df.y$F1.pi)),
        F2=mean(abs(df$pitrue - df.y$F2.pi)),
        F3=mean(abs(df$pitrue - df.y$F3.pi)),
        SEF12=mean(abs(df$pitrue - df.y$SEF12.pi)),
        SEF123=mean(abs(df$pitrue - df.y$SEF123.pi)))

    pw.mae <- c(hetGP=mean(abs(df$pwtrue - df$pwpred)),
        F1=mean(abs(df$pwtrue - df.y$F1.pw)),
        F2=mean(abs(df$pwtrue - df.y$F2.pw)),
        F3=mean(abs(df$pwtrue - df.y$F3.pw)),
        SEF12=mean(abs(df$pwtrue - df.y$SEF12.pw)),
        SEF123=mean(abs(df$pwtrue - df.y$SEF123.pw)))

    ti.mae <- c(hetGP=mean(abs(df$titrue - df$tipred)),
        F1=mean(abs(df$titrue - df.y$F1.ti)),
        F2=mean(abs(df$titrue - df.y$F2.ti)),
        F3=mean(abs(df$titrue - df.y$F3.ti)),
        SEF12=mean(abs(df$titrue - df.y$SEF12.ti)),
        SEF123=mean(abs(df$titrue - df.y$SEF123.ti)))

    mae <- data.frame(pi=pi.mae, pw=pw.mae, ti=ti.mae)
}

xtable(ls.tab)
