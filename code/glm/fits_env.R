## This code contains functions to fit models and make predictions for
## the environmental covariates.

preds.precip <- function(data, f=rep(1/10, 10))
{
  n <- nrow(data)

  Xy.raw <- data.frame(
    ly=log(data$prec+1),
    temp=data$tavg
  )

  Xy.smooth <- Xy.raw
  # f <- c(0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1)
  ## f<- c(0.25, 0.15,  rep(0.1, 4), rep(0.05, 4))
  for(nm in names(Xy.raw)) {
    Xy.smooth[nm] <- filter(Xy.raw[nm], f, sides=1)
  }
  Xy.smooth <- Xy.smooth[-(1:length(f)),]
  n <- nrow(Xy.smooth)

  Xy <- data.frame(
      t=1:(n-52),
      ly=log(data$prec+1)[(53+length(f)):nrow(data)],
      ly.m1=Xy.smooth$ly[52:(n-1)],
      ly.m52=Xy.smooth$ly[1:(n-52)],
      temp.m1=Xy.smooth$temp[52:(n-1)],
      temp2.m1=Xy.smooth$temp[52:(n-1)]^2,
      s05=sin((53:n)*2*pi/26), c05=cos((53:n)*2*pi/26),
      s=sin((53:n)*2*pi/52), c=cos((53:n)*2*pi/52),
      s2=sin((53:n)*2*pi/104), c2=cos((53:n)*2*pi/104),
      s4=sin((53:n)*2*pi/208), c4=cos((53:n)*2*pi/208)
    )

  return(Xy)
}


fit.precip <- function(data, f=rep(1/10, 10))
{
  Xy <- preds.precip(data, f)

  fit <- lm(ly~., data=Xy)
  
  fit.bic <- step(fit, scope=formula(fit), direction="both", 
    k=log(nrow(Xy)), trace=0)

  return(list(orig=fit, bic=fit.bic, f=f))
}


fcast.precip <- function(fit, data)
{

  ## build a data frame for smoothing predictors
  Xy <- preds.precip(data, fit$f)
  if(any(is.na(Xy[,-which(names(Xy)=="ly")]))) stop("NA values in preds")

  p <- predict(fit$bic, newdata=Xy[nrow(Xy),,drop=FALSE], se.fit=TRUE)

  return(p)
}


preds.temp <- function(data, f=rep(1/10,10), my.phase=2)
{

  Xy.raw <- data.frame(
    lprecip=log(data$prec+1),
    temp=data$tavg
  )

  Xy.smooth <- Xy.raw
  # f <- c(0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1)
  ## f<- c(0.25, 0.15,  rep(0.1, 4), rep(0.05, 4))
  for(nm in names(Xy.raw)) {
    Xy.smooth[nm] <- filter(Xy.raw[nm], f, sides=1)
  }
  Xy.smooth <- Xy.smooth[-(1:length(f)),]
  n <- nrow(Xy.smooth)

  phase<-my.phase*(2*pi)/52
  Xy <- data.frame(
      t=1:(n-52),
      lprecip.m1=Xy.smooth$lprecip[52:(n-1)],
      temp=data$tavg[(53+length(f)):nrow(data)],
      temp.m1=Xy.smooth$temp[52:(n-1)],
      temp.m52=Xy.smooth$temp[1:(n-52)],
      lprecip2.m1=Xy.smooth$lprecip[52:(n-1)]^2,
      s05=sin((53:n)*2*pi/26-phase*2), c05=cos((53:n)*2*pi/26-phase*2),
      s=sin((53:n)*2*pi/52-phase), c=cos((53:n)*2*pi/52-phase),
      s2=sin((53:n)*2*pi/104-phase/2), c2=cos((53:n)*2*pi/104-phase/2),
      s4=sin((53:n)*2*pi/208-phase/4), c4=cos((53:n)*2*pi/208-phase/4)
    )

  return(Xy)
}


fit.temp <- function(data, f=rep(1/10, 10), my.phase=2)
{

  Xy <- preds.temp(data, f, my.phase=my.phase)

  fit <- lm(temp~., data=Xy)
  fit.m<- lm(temp ~ t + s + c , data=Xy)
  
  fit.bic <- step(fit, scope=list(upper=formula(fit), lower=formula(fit.m)), direction="both", 
    k=log(nrow(Xy)), trace=0)

  ##plot(fit.bic$fitted, type="l")
  ##points(Xy$temp, col=2)

  return(list(orig=fit, bic=fit.bic, f=f, my.phase=my.phase))
}


fcast.temp <- function(fit, data)
{

  ## build a data frame for smoothing predictors
  Xy <- preds.temp(data, fit$f, my.phase=fit$my.phase)
  if(any(is.na(Xy[,-which(names(Xy)=="temp")]))) stop("NA values in preds")

  p <- predict(fit$bic, newdata=Xy[nrow(Xy),,drop=FALSE], se.fit=TRUE)

  return(p)
}


preds.ndvi <- function(data, w="NDVI.18.45.66.14.", f=rep(1/10, 10))
{
  data$ndvi <- data[[w]]

  Xy.raw <- data.frame(
    ndvi=data$ndvi,
    ly=log(data$prec+1),
    temp=data$tavg
  )

  Xy.smooth <- Xy.raw
  for(nm in names(Xy.raw)) {
    Xy.smooth[nm] <- filter(Xy.raw[nm], f, sides=1)
  }
  Xy.smooth <- Xy.smooth[-(1:length(f)),]
  n <- nrow(Xy.smooth)

  Xy <- data.frame(
      t=1:(n-52),
      ndvi=data$ndvi[(53+length(f)):nrow(data)],
      ndvi.m1=Xy.smooth$ndvi[52:(n-1)],
      ## ndvi.m1=ndvi[52:(n-1)],
      ndvi.m52=Xy.smooth$ndvi[1:(n-52)],
      ly.m1=Xy.smooth$ly[52:(n-1)],
      temp.m1=Xy.smooth$temp[52:(n-1)],
      temp2.m1=Xy.smooth$temp[52:(n-1)]^2,
      ly2.m1=Xy.smooth$ly[52:(n-1)]^2,
      s05=sin((53:n)*2*pi/26), c05=cos((53:n)*2*pi/26),
      s=sin((53:n)*2*pi/52), c=cos((53:n)*2*pi/52),
      s2=sin((53:n)*2*pi/104), c2=cos((53:n)*2*pi/104),
      s4=sin((53:n)*2*pi/208), c4=cos((53:n)*2*pi/208)
    )

  return(Xy)
}


fit.ndvi <- function(data, w="NDVI.18.45.66.14.", f=rep(1/10, 10))
{
  Xy <- preds.ndvi(data, w, f)

  fit <- lm(ndvi~., data=Xy)
  
  fit.bic <- step(fit, scope=formula(fit), direction="both", 
    k=log(nrow(Xy)), trace=0)

  # plot(fit.bic$fitted, type="l")
  # points(Xy$ndvi, col=2)

  return(list(orig=fit, bic=fit.bic, f=f, w=w))
}

fcast.ndvi <- function(fit, data)
{
  ## build a data frame for smoothing predictors
  Xy <- preds.ndvi(data, w=fit$w, fit$f)
  if(any(is.na(Xy[,-which(names(Xy)=="ndvi")]))) stop("NA values in preds")

  p <- predict(fit$bic, newdata=Xy[nrow(Xy),,drop=FALSE], se.fit=TRUE)

  return(p)
}




preds.nino12 <- function(data, f=rep(1/10, 10))
{
  Xy.raw <- data.frame(
                       nino12=data$nino12,
                       ##nino34=data$nino34,
                       soi=data$soi
                       )

  Xy.smooth <- Xy.raw
  for(nm in names(Xy.raw)) {
    Xy.smooth[nm] <- filter(Xy.raw[nm], f, sides=1)
  }
  Xy.smooth <- Xy.smooth[-(1:length(f)),]
  n <- nrow(Xy.smooth)

  phase<-2
  Xy <- data.frame(
                   t=1:(n-52),
                   nino12=data$nino12[(53+length(f)):nrow(data)],
                   nino12.m1=Xy.smooth$nino12[52:(n-1)],
                   ##nino34.m1=Xy.smooth$nino34[52:(n-1)],
                   soi.m1=Xy.smooth$soi[52:(n-1)],
                   s05=sin((53:n-phase)*2*pi/26), c05=cos((53:n-phase)*2*pi/26),
                   s=sin((53:n-phase)*2*pi/52), c=cos((53:n-phase)*2*pi/52),
                   s2=sin((53:n-phase)*2*pi/104), c2=cos((53:n-phase)*2*pi/104),
                   s4=sin((53:n-phase)*2*pi/208), c4=cos((53:n-phase)*2*pi/208)
    )

  return(Xy)
}


fit.nino12 <- function(data, f=rep(1/10, 10))
{
  Xy <- preds.nino12(data, f)

  fit <- lm(nino12 ~., data=Xy)
  fit.m<-lm(nino12 ~ s + c, data=Xy)
  
  fit.bic <- step(fit, scope=list(upper=formula(fit), lower=formula(fit.m)), direction="both", 
    k=log(nrow(Xy)), trace=0)

  ##plot(fit.bic$fitted, type="l")
  ##points(Xy$nino12, col=2)

  return(list(orig=fit, bic=fit.bic, f=f))
}


fcast.nino12 <- function(fit, data)
{
  ## build a data frame for smoothing predictors
  Xy <- preds.nino12(data, fit$f)
  if(any(is.na(Xy[,-which(names(Xy)=="nino12")]))) stop("NA values in preds")

  p <- predict(fit$bic, newdata=Xy[nrow(Xy),,drop=FALSE], se.fit=TRUE)

  return(p)
}





preds.soi <- function(data, f=rep(1/10, 10))
{
  Xy.raw <- data.frame(
                       nino12=data$nino12,
                       ##nino34=data$nino34,
                       soi=data$soi
                       )

  Xy.smooth <- Xy.raw
  for(nm in names(Xy.raw)) {
    Xy.smooth[nm] <- filter(Xy.raw[nm], f, sides=1)
  }
  Xy.smooth <- Xy.smooth[-(1:length(f)),]
  n <- nrow(Xy.smooth)

  phase<-2*(2*pi)
  Xy <- data.frame(
                   t=1:(n-52),
                   soi=data$soi[(53+length(f)):nrow(data)],
                   nino12.m1=Xy.smooth$nino12[52:(n-1)],
                   ##nino34.m1=Xy.smooth$nino34[52:(n-1)],
                   soi.m1=Xy.smooth$soi[52:(n-1)],
                   soi.m2=Xy.smooth$soi[51:(n-2)],
                   soi.m3=Xy.smooth$soi[50:(n-3)],
                   soi.m10=Xy.smooth$soi[(53-10):(n-10)],
                   s05=sin((53:n)*2*pi/26 +phase/26), c05=cos((53:n)*2*pi/26+phase/26),
                   s=sin((53:n)*2*pi/52+phase/52), c=cos((53:n)*2*pi/52+phase/52),
                   s2=sin((53:n)*2*pi/104+phase/104), c2=cos((53:n)*2*pi/104+phase/104),
                   s4=sin((53:n)*2*pi/208+phase/208), c4=cos((53:n)*2*pi/208+phase/208)
    )

  return(Xy)
}


fit.soi <- function(data, f=rep(1/10, 10))
{
  Xy <- preds.soi(data, f)

  fit <- lm(soi ~., data=Xy)

  fit.bic <- step(fit, scope=formula(fit), direction="both", 
    k=log(nrow(Xy)), trace=0)

  ##plot(fit.bic$fitted, type="l")
  ##points(Xy$soi, col=2)

  return(list(orig=fit, bic=fit.bic, f=f))
}


fcast.soi <- function(fit, data)
{
  ## build a data frame for smoothing predictors
  Xy <- preds.soi(data, fit$f)
  if(any(is.na(Xy[,-which(names(Xy)=="soi")]))) stop("NA values in preds")

  p <- predict(fit$bic, newdata=Xy[nrow(Xy),,drop=FALSE], se.fit=TRUE)

  return(p)
}




