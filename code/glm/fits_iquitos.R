library(MASS)
## Here is code to make the R0 predictor

source("fits_env.R")

build.R0 <- function(data)
{
    R0.T <- read.csv("../data/DENV_RelR0_scaled_NEW.csv")
    R0.T$temp <- round(R0.T$T+273,1)
    tavg.round<-round(data$tavg, 1)
    T.start<-which(R0.T$temp==min(tavg.round))
    if(length(T.start)==0) T.start<-which(R0.T$temp==min(R0.T$temp))
    T.end<-which(R0.T$temp==max(tavg.round))
    if(length(T.end)==0) T.end<-which(R0.T$temp==max(R0.T$temp))

    R0<-rep(NA, nrow(data))
    for(i in T.start:T.end){
        ww<-NULL
        tt<-round(R0.T$temp[i],1)
        ww<-which(tavg.round==tt)
        if(length(ww)!= 0) R0[ww]<-R0.T$scaleR0[i]  
    }

    if( T.start==1 | T.end == length(R0.T$temp) ){
        ww<-NULL
        ww<-c( which(tavg.round<R0.T$temp[T.start]), which(tavg.round>R0.T$temp[T.end]) )
        R0[ww]<-0
    }
    
    return(R0)
}

build.GM <- function(data){

  GM<- rep(NA, length(data$total_cases))
  s<-53
  for(i in 2:length(levels(data$season))){
    nstart <-max(which(data$season==levels(data$season)[i-1]))
    GM[(s):(s+51)]<-as.numeric(tapply(data$total_cases[1:nstart], INDEX=data$season_week[1:nstart], FUN=mean))
    s<-s+52             
  }
  return(GM)
}

build.LGM <- function(data){

  LGM<- rep(NA, length(data$total_cases))
  s<-53
  for(i in 2:length(levels(data$season))){
    nstart <-max(which(data$season==levels(data$season)[i-1]))
    LGM[(s):(s+51)]<-as.numeric(tapply(log(data$total_cases[1:nstart]+1), INDEX=data$season_week[1:nstart], FUN=mean))
    s<-s+52             
  }
  return(LGM)
}

preds.totalcases <- function(data, R0, f=rep(1/10, 10), fdy=rep(1/4,4))
  {

    ## build a data frame for smoothing predictors
    Xy.raw <- data.frame(
                         y=data$total_cases,
                         ly=log(data$total_cases+1),
                         lgm=data$lgm, 
                         season=data$season,
                         w=data$season_week,
                         ##epi=as.factor(data$epi),
                         lpop=log(data$pop+1),
                         lp=log(data$prec+1),
                         tavg=data$tavg,
                         ndvi.avg=data$NDVIavg,
                         R0=R0,
                         nino12=data$nino12,
                         ##nino34=data$nino34,
                         soi=data$soi)

    ##w<-which(is.na(Xy.raw$ndvi50))
    ##for(i in w) Xy.raw$ndvi50[i]<-mean(sanjuan$NDVI.18.50..66.14[(i-1):(i+1)], na.rm=TRUE)

    if(any(is.na(Xy.raw$ndvi.avg))) stop("missing value in raw NDVI series")

    ## filter the predictors
    Xy.smooth <- Xy.raw
    for(nm in names(Xy.raw)) {
      if(is.element(nm, c("y", "season", "w", "epi"))) next; ##, "gm", "lgm"
      Xy.smooth[nm] <- filter(Xy.raw[nm], f, sides=1)
    }
    
    ## add in cumulative infected variable over years
    Xy.smooth$ci <- filter(Xy.raw$y, rep(1, 52), sides=1)

    ## add in slightly smoothed changes in observed cases
    if(is.null(data$dy)){
      Xy.smooth$dy <- filter(c(NA,diff(Xy.raw$y, lag=1)), f=fdy, sides=1)
    }else Xy.smooth$dy<-data$dy

    ## trim to remove NAs
    Xy.smooth <- Xy.smooth[-(1:length(f)),]
    n <- nrow(Xy.smooth)

    if(any(is.na(Xy.smooth$ndvi.avg))) stop("missing value in smoothed NDVI series")

    l.lim<-52*3

    phase<-6*(2*pi/52)
    ## only consider the subset
    Xy <- data.frame(
                     t=1:(n-l.lim),
                     ##epi.m1=Xy.smooth$epi[52:(n-1)],
                     ci.m1=Xy.smooth$ci[l.lim:(n-1)],
                     ##ci2.m1=Xy.smooth$ci[l.lim:(n-1)]^2,
                     y=Xy.smooth$y[(l.lim+1):n],
                     ly.m1=Xy.smooth$ly[l.lim:(n-1)],
                     ##ly2.m1=(Xy.smooth$ly[l.lim:(n-1)])^2,
                     ##ly.m52=Xy.smooth$ly[1:(n-l.lim)],
                     lgm=Xy.smooth$lgm[(l.lim+1):n],
                     ##dy.m1=Xy.smooth$dy[l.lim:(n-1)],
                     ##dy.m10=Xy.smooth$dy[(l.lim-9):(n-10)],
                     ## sp.m1=sanjuan$satprecip[52:(n-1)], ## missing data
                     ##lpop=Xy.smooth$lpop[(l.lim+1):n],
                     lp.m1=Xy.smooth$lp[l.lim:(n-1)],
                     lp2.m1=(Xy.smooth$lp[l.lim:(n-1)])^2,
                     tavg.m1=Xy.smooth$tavg[(l.lim):(n-1)],
                     ##tavg.m11=Xy.smooth$tavg[(l.lim-10):(n-1-10)],
                     tavg2.m1=(Xy.smooth$tavg[(l.lim):(n-1)])^2,
                     ndvi.avg.m1=Xy.smooth$ndvi.avg[l.lim:(n-1)], ## missing data
                     R0.m1=Xy.smooth$R0[l.lim:(n-1)],
                     R0.m5=Xy.smooth$R0[(l.lim-4):(n-5)],
                     nino12.m1=Xy.smooth$nino12[l.lim:(n-1)],
                     ##nino12.m6=Xy.smooth$nino12[(l.lim-5):(n-6)],
                     ##nino12.m32=Xy.smooth$nino12[(l.lim-31):(n-32)],
                     ##nino34.m1=Xy.smooth$nino34[52:(n-1)],
                     soi.m1=Xy.smooth$soi[l.lim:(n-1)],
                     ##soi.m24=Xy.smooth$soi[(l.lim-23):(n-24)],
                     s05=sin(((l.lim+1):n)*2*pi/26-phase*2), c05=cos(((l.lim+1):n)*2*pi/26-phase*2),
                     s=sin(((l.lim+1):n)*2*pi/52-phase), c=cos(((l.lim+1):n)*2*pi/52-phase)
                     ##s2=sin(((l.lim+1):n)*2*pi/104-phase/2), c2=cos(((l.lim+1):n)*2*pi/104-phase/2),
                     ##s4=sin(((l.lim+1):n)*2*pi/208-phase/4), c4=cos(((l.lim+1):n)*2*pi/208-phase/4)
          )

    return(Xy)
  }


fit.totalcases <- function(data, f=rep(1/10, 10))
{
  R0 <- build.R0(data)

  Xy <- preds.totalcases(data, R0, f)

  fit <- glm(y~., data=Xy, family="poisson")
  fit.m <- glm(y ~ lgm, data=Xy, family="poisson")

  fit.bic <- step(fit, scope=list(upper=formula(fit), lower=formula(fit.m)), direction="both", 
    k=log(nrow(Xy)), trace=0)
  print(summary(fit.bic))
  ## fit.aic <- step(fit, scope=formula(fit), direction="both")


  ##plot(fit.bic$fitted, type="l", ylim=range(Xy$y))
  ##points(Xy$y, pch=19, cex=0.5)

  return(list(orig=fit, bic=fit.bic, min=fit.m, f=f))
}

fcast.totalcases <- function(fit, data)
{
  R0 <- build.R0(data)

  ## build a data frame for smoothing predictors
  Xy <- preds.totalcases(data, R0, fit$f)
  if(any(is.na(Xy[,-which(names(Xy)=="y")]))) stop("NA values in preds")

  p <- predict(fit$bic, newdata=Xy[nrow(Xy),,drop=FALSE], 
    type="response" , se.fit=TRUE)

  return(p)
}


fit.totalcases.nb <- function(data, f=rep(1/10, 10))
{
  R0 <- build.R0(data)

  Xy <- preds.totalcases(data, R0, f)

  fit <- glm.nb(y~. ,data=Xy)
  fit.m <- glm.nb(y~ lgm, data=Xy)

  fit.bic <- step(fit, scope=list(upper=formula(fit), lower=formula(fit.m)),## formula(fit),
                  direction="both", k=log(nrow(Xy)), trace=0)
  print(summary(fit.bic))
  ## fit.aic <- step(fit, scope=formula(fit), direction="both")


  #plot(fit.bic$fitted, type="l", ylim=range(Xy$y))
  #points(Xy$y, cex=0.5)
  #lines(qnbinom(p=0.025, size=fit.bic$theta, mu=fit.bic$fitted), lty=2, col=2)
  #lines(qnbinom(p=0.925, size=fit.bic$theta, mu=fit.bic$fitted), lty=2, col=2)

  return(list(orig=fit, bic=fit.bic, f=f))
}

fcast.totalcases.nb <- function(fit, data)
{
  R0 <- build.R0(data)

  ## build a data frame for smoothing predictors
  Xy <- preds.totalcases(data, R0, fit$f)
  if(any(is.na(Xy[,-which(names(Xy)=="y")]))) stop("NA values in preds")

  p <- predict(fit$bic, newdata=Xy[nrow(Xy),,drop=FALSE], 
    type="response", se.fit=TRUE)

  return(p)
}

