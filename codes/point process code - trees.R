##################################
#### Analyse the Finland tree data
##################################
library(spatstat)
data(hyytiala)

#### Plot the data
plot(hyytiala)


#### Estimate the K function
par(mfrow=c(2,1))
Kc <- Kest(hyytiala, correction="Ripley")
plot(envelope(hyytiala, Kest,nsim=39))
plot(Kc,  iso - theo ~ r, main="Deviation from CSR")
abline(h=0, lty=2, col="gray40")


#### Plot the estimated intensity surface
plot(density(hyytiala, bw.diggle), main="Estimated intensity")
