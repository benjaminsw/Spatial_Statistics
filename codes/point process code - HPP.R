## ======================================================================
## A simulated homogeneous Poisson process
## ======================================================================
library(spatstat)


## Generate and plot the point pattern
hpp.50 <- rpoispp(lambda=50, win=owin(xrange=c(0,1), yrange=c(0,1)))
hpp.100 <- rpoispp(lambda=100, win=owin(xrange=c(0,1), yrange=c(0,1)))
hpp.200 <- rpoispp(lambda=200, win=owin(xrange=c(0,1), yrange=c(0,1)))
hpp.500 <- rpoispp(lambda=500, win=owin(xrange=c(0,1), yrange=c(0,1)))

par(mfrow=c(2,2))
plot(hpp.50, pch=19, main="lambda=50")
plot(hpp.100, pch=19, main="lambda=100")
plot(hpp.200, pch=19, main="lambda=200")
plot(hpp.500, pch=19, main="lambda=500")




## Quadrat test
Q <- quadratcount(hpp.100, nx=3, ny=3)
par(mfrow=c(1,1))
plot(hpp.100, pch="+", main="")
plot(Q, add=T, col="red")

quadrat.test(hpp.100, nx=3, ny=3, method="MonteCarlo", nsim=5000)



## Estimate the K function
par(mfrow=c(2,1))
Kc <- Kest(hpp.100, correction="Ripley")
plot(envelope(hpp.100, Kest,nsim=39))
plot(Kc,  iso - theo ~ r, main="Deviation from CSR")
abline(h=0, lty=2, col="gray40")



plot(density(hpp.sim, bw.diggle))


#### IPP
ipp <- rLGCP(model="exp", mu=7, param=list(var=1, scale=0.1),
             win=owin(xrange=c(0,1), yrange=c(0,1)))
plot(ipp, pch=19, main="Realisation of a IPP")
