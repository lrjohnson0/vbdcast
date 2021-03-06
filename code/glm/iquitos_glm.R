## Set the observation times
augs <- c(0, seq(4, 48, by=4))

## which season/s to make predictions for
## for testing we show one here. 
seasons<- c( "2007/2008")

## all seasons are here

##seasons<-c("2005/2006", "2006/2007", "2007/2008", "2008/2009", "2009/2010","2010/2011" "2011/2012" "2012/2013")

## define the number of replicates for prediction, set to 10 for testing. Paper results based on 1000 replicates
reps <- 10

## read in the command line arguments
args <- commandArgs(TRUE)
if(length(args) > 0)
  for(i in 1:length(args))
    eval(parse(text=args[[i]]))

## echo arguments
cat("augs is ", paste(augs, collapse=", "), "\n", sep="")
cat("seasons is ", paste(seasons, collapse=", "), "\n", sep="")
season.names <- sub("/", "-", seasons)
cat("season.names is ", paste(season.names, collapse=", "), "\n", sep="")
cat("reps is ", reps, "\n", sep="")

source("fits_iquitos.R")
source("fits_env.R")
## the testing data
iquitos <- read.csv("../data/combined_iquitos_testing.csv")

## simple fill in for remaining NAs in the NDVI data for Iquitos:

w<-which(is.na(iquitos$NDVI...3.75..73.24.))
for(i in w) iquitos$NDVI...3.75..73.24.[i]<-mean(iquitos$NDVI...3.75..73.24.[(i-2):(i+2)], na.rm=TRUE)

##training
##NDVIS<-iquitos[,15:18]

##testing
NDVIS<-NULL
for(i in 16:19) NDVIS<-cbind(NDVIS, iquitos[,i])

iquitos$NDVIavg<-rowMeans(NDVIS)
w<-which(is.na(iquitos$NDVIavg))
for(i in w) iquitos$NDVIavg[i]<-mean(iquitos$NDVIavg[(i-1):(i+1)], na.rm=TRUE)

w<-which(is.na(iquitos$adjpop))
iquitos$adjpop[w]<-iquitos$pop[w]

## Add in the rolling grand mean of the log cases
iquitos$lgm<-build.LGM(iquitos)


for(ss in 1:length(seasons)){
    
    ## initialize the output object
    output<-list()
  
    ## set the counter for the output
    cnt<-1
    
    s.lim<-min(which(iquitos$season==seasons[ss]))
    
    for(i in 1:length(augs)){
        
        aug<-augs[i]
        
        print(paste("season: ", seasons[ss], "; week: ", aug, sep=""))
        
        nstart <-min(which(iquitos$season==seasons[ss]))-1 + aug
        fit.p<-fit.t<-fit.ndvi50<-fit.ndvi40<-fit.n12<-fit.s<-fit<-NULL


        ## build the new data from of the appropriate length for the
        ## season and week, and fill in with things that we know
        data<-NULL
        data <- iquitos[1:nstart,]
        
        new.df <- data.frame(matrix(NA, ncol=ncol(data), nrow=52-aug))
        names(new.df) <- names(data)
        data <- rbind(data, new.df)
        data[,1:5] <- iquitos[1:nrow(data),1:5]
        data$lgm<-iquitos$lgm[1:nrow(data)]
        data$pop<-iquitos$pop[1:nrow(data)]
        data$adjpop<-iquitos$adjpop[1:nrow(data)]
        
        
        
        ## define the filter that we want to use to smooth the predictors
        f<-rep(1/10, 10)

        ## fit models for the predictors
        fit.p <- fit.precip(data[1:nstart,], f=f)
        
        fit.t <- fit.temp(data[1:nstart,], f=f, my.phase=7)
        
        fit.ndviavg <- fit.ndvi(data[1:nstart,], w="NDVIavg", f=f)
        
        fit.n12 <- fit.nino12(data[1:nstart,], f=f)
        
        fit.s <- fit.soi(data[1:nstart,], f=f)

        ## fit the model for total cases
        fit <- fit.totalcases.nb(data[1:nstart,], f=f)
        
        
        total_cases <- matrix(NA, nrow=reps, ncol=52-aug)
        nino12 <- soi <-prec <- temp <- ndvi<- total_cases
        
        
        for(r in 1:reps) {
            ## make sure these are reset to NA for each replicate,
            ## since we have the full data already in the dataframe
            data$prec[(nstart+1):nrow(data)]<-data$tavg[(nstart+1):nrow(data)]<-rep(NA, length((nstart+1):nrow(data)))
            data$NDVIavg[(nstart+1):nrow(data)]<-data$nino12[(nstart+1):nrow(data)]<-rep(NA, length((nstart+1):nrow(data)))
            data$soi[(nstart+1):nrow(data)]<-data$total_cases[(nstart+1):nrow(data)]<-rep(NA, length((nstart+1):nrow(data)))
            
            for(n in (nstart+1):nrow(data)) {
                s<-NULL    
                p.precip <- fcast.precip(fit.p, data[1:n,])
                s <- sqrt(p.precip$se.fit^2 + p.precip$residual.scale^2)
                data$prec[n] <- exp(rnorm(1, p.precip$fit, s))-1  ## precip is log-normal
                
                s<-NULL
                p.temp <- fcast.temp(fit.t, data[1:n,])
                s <- sqrt(p.temp$se.fit^2 + p.temp$residual.scale^2)
                data$tavg[n] <-  rnorm(1, p.temp$fit, s) 
                
                s<-NULL
                p.ndvi <- fcast.ndvi(fit.ndviavg, data[1:n,])
                s <- sqrt(p.ndvi$se.fit^2 + p.ndvi$residual.scale^2)
                data$NDVIavg[n] <- max(0, rnorm(1, p.ndvi$fit, p.ndvi$se.fit)) ## very
                                                                               ## small
                                                                               ## NDVIs,
                                                                               ## but
                                                                               ## this
                                                                               ## cannot
                                                                               ## go
                                                                               ## below
                                                                               ## 0
                
                s<-NULL
                p.nino12 <- fcast.nino12(fit.n12, data[1:n,])
                s <- sqrt(p.nino12$se.fit^2 + p.nino12$residual.scale^2)
                data$nino12[n] <-   rnorm(1, p.nino12$fit, s) 
                
                s<-NULL
                p.soi <- fcast.soi(fit.s, data[1:n,])
                s <- sqrt(p.soi$se.fit^2 + p.soi$residual.scale^2)
                data$soi[n] <-   rnorm(1, p.soi$fit, s) 
                
                p<-NULL
                p <- fcast.totalcases(fit, data[1:n,])
                data$total_cases[n] <- rnbinom(1, mu=p$fit, size=fit$bic$theta)
            }
            
            
            total_cases[r,] <- data$total_cases[(nstart+1):nrow(data)]
            prec[r,] <- data$prec[(nstart+1):nrow(data)]
            temp[r,] <- data$tavg[(nstart+1):nrow(data)]
            nino12[r,] <- data$nino12[(nstart+1):nrow(data)]
            soi[r,] <- data$soi[(nstart+1):nrow(data)]
            ndvi[r,] <- data$NDVIavg[(nstart+1):nrow(data)]
            
            if(r%%10==0) print(r)
        }
        
        if(aug!=0){
            dats<-matrix(NA, nrow=reps, ncol=aug)
            for(j in 1:aug) dats[,j] <- rep(data$total_cases[nstart-aug+j], reps)
            total_cases<-cbind(dats, total_cases)
        }
        
        output[[cnt]]<-list(season=seasons[ss], week=aug, total.cases=total_cases)
        f.name<-paste("iquitos_forecast_", season.names[ss], "_wk", aug, "_outs_testing.RData", sep="")
        save(output, file=f.name) 
        cnt<-cnt+1
        
    }
    
    output[[cnt]]<-list(season=seasons[ss], week=aug, total.cases=total_cases)
    f.name<-paste("iquitos_forecast_", season.names[ss], "_wks", augs[1], "-", aug, "_outs_testing.RData", sep="")
    cnt<-cnt+1
}

