#############################
#### Analyse the Gorilla data
#############################
library(spatstat)
data(gorillas)


#### Plot the data
plot(gorillas, main="", cols=c("black"), use.marks=F, pch=4)


#### Do the quadrat test
Q <- quadratcount(gorillas, nx=5, ny=5)
par(mfrow=c(1,1))
plot(gorillas, pch="+", cols=c("black"), use.marks=F, main="")
plot(Q, add=T, col="red")
quadrat.test(gorillas, nx=5, ny=5, method="MonteCarlo", nsim=5000)


#### Estimate the K function
Kc <- Kest(gorillas, correction="Ripley")
plot(envelope(gorillas, Kest,nsim=39))
plot(Kc,  iso - theo ~ r, main="Deviation from CSR")
abline(h=0, lty=2, col="gray40")


#### Plot the intensity function
plot(density(gorillas, bw.diggle), main="Estimated intensity")



