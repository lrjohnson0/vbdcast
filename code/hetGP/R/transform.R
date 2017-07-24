f <- function(x) {
	y <- x-1
	y[x < 1] <- log(x[x < 1])
	return(y)
}

fi <- function(y) {
	x <- y+1
	x[y < 0] <- exp(y[y < 0])
	return(x)
}

f2 <- function(x) {
	y <- sqrt(x)-1
	y[x < 1] <- log(x[x < 1])
	return(y)
}

f2i <- function(y) {
	x <- (y+1)^2
	x[y < 0] <- exp(y[y < 0])
	return(x)
}

## for visualization
if(0) {
x <- seq(0,100,length=1000)
y <- f2(x+1)
pdf("transform.pdf", height=5, width=10.5)
par(mfrow=c(1,2))
plot(x, y, type="l", lwd=2, xlab="x", ylab="f(x)", main="forward")
yprime <- seq(-2,8.5, length=100)
xprime <- f2i(yprime)-1
plot(yprime, xprime, type="l", lwd=2, xlab="y", ylab="fi(y)", main="inverse")
dev.off()
}