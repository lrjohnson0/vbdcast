wk <- rep(0:51, 9)
seasons <- rep("2004/2005", 52)
seasons <- c(seasons, rep("2005/2006", 52))
seasons <- c(seasons, rep("2006/2007", 52))
seasons <- c(seasons, rep("2007/2008", 52))
seasons <- c(seasons, rep("2008/2009", 52))
seasons <- c(seasons, rep("2009/2010", 52))
seasons <- c(seasons, rep("2010/2011", 52))
seasons <- c(seasons, rep("2011/2012", 52))
seasons <- c(seasons, rep("2012/2013", 52))

## Peak Incidence

## read the data
pif <- list.files(pattern="-pi")
df <- NULL
nms <- NULL
for(f in pif) {
	nms <- c(nms, strsplit(f, "-")[[1]][1])
	tab <- read.csv(f, header=FALSE)
	df <- cbind(df, as.numeric(as.matrix(tab)))
}
df <- data.frame(cbind(seasons, wk, df))
names(df) <- c("season", "week", paste(nms, "-pi", sep=""))

## remove 2004/5
df <- df[-c(1:52),]
## only save every 4
df <- df[seq(1,nrow(df), by=4),]

## save result
## write.csv(df, file="pi_sanjuan.csv", row.names=FALSE, quote=FALSE)
df.all <- df

## Peak Week

pif <- list.files(pattern="-pw")
df <- NULL
nms <- NULL
for(f in pif) {
	nms <- c(nms, strsplit(f, "-")[[1]][1])
	tab <- read.csv(f, header=FALSE)
	df <- cbind(df, as.numeric(as.matrix(tab)))
}

df <- data.frame(cbind(seasons, wk, df))
names(df) <- c("season", "week", paste(nms, "-pw", sep=""))

## remove 2004/5
df <- df[-c(1:52),]
## only save every 4
df <- df[seq(1,nrow(df), by=4),]

## save result
## write.csv(df, file="pw_sanjuan.csv", row.names=FALSE, quote=FALSE)
df.all <- cbind(df.all, df)

## Total Incidence

pif <- list.files(pattern="-ti")
df <- NULL
nms <- NULL
for(f in pif) {
	nms <- c(nms, strsplit(f, "-")[[1]][1])
	tab <- read.csv(f, header=FALSE)
	df <- cbind(df, as.numeric(as.matrix(tab)))
}

df <- data.frame(cbind(seasons, wk, df))
names(df) <- c("season", "week", paste(nms, "-ti", sep=""))

## remove 2004/5
df <- df[-c(1:52),]
## only save every 4
df <- df[seq(1,nrow(df), by=4),]

## save result
## write.csv(df, file="ti_sanjuan.csv", row.names=FALSE, quote=FALSE)
df.all <- cbind(df.all, df)

## write out
write.csv(df.all, file="sanjuan.csv", row.names=FALSE, quote=FALSE)
