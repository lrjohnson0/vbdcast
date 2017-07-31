augs<-c(0, seq(4, 48, by=4))

loc<-"Iquitos"
if(loc=="Iquitos") dat <- read.csv("../data/combined_iquitos_testing.csv")
if(loc=="San Juan") dat <- read.csv("../data/combined_sanjuan_testing.csv")

## which seasons to combine. Just one for testing.
seasons<- c("2007-2008")

## all the seasons
## seasons<-c("2005-2006", "2006-2007", "2007-2008", "2008-2009",
##           "2009-2010", "2010-2011", "2011-2012", "2012-2013") 


out.obj<-list()

cnt<-1
for(j in 1:length(seasons)){
    if(loc=="San Juan"){
        fname<-paste("sanjuan_forecast_", seasons[j], "_wk", 48,
                     "_outs_testing.RData", sep="")
    }
    if(loc=="Iquitos"){
        fname<-paste("iquitos_forecast_", seasons[j], "_wk", 48,
                     "_outs_testing.RData", sep="")
    }
    
    load(fname)
    for(i in 1:length(augs)){
        w<-which(dat$season==output[[i]]$season)
        out.obj[[cnt]]<-list(season=output[[i]]$season, week=output[[i]]$week,
                     tc=output[[i]]$total.cases,
                     ytrue=dat$total_cases[w])
    cnt<-cnt+1
  }
}

if(loc=="San Juan") save(out.obj, file="glm_dengue_sanjuan_testing_final.RData")
if(loc=="Iquitos") save(out.obj, file="glm_dengue_iquitos_testing_final.RData")

