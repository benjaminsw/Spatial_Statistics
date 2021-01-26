#############################################
#### Analyse the South Lancashire cancer data
#############################################
library(splancs)
data(southlancs)


#### Plot the data
par(mfrow=c(1,1))
plot(southlancs.bdy, type="l", xlab="Easting", ylab="Northing")
points(southlancs[ ,-3], pch=19)


#### Plot the intensity function
window <- owin(poly=southlancs.bdy[seq(345,1,-1), ])
data <- ppp(x=southlancs$x, y=southlancs$y, window=window)
plot(density(data, bw.diggle), main="Estimated intensity")


#### Estimate the K function
Kc <- Kest(data, correction="Ripley")
plot(Kc, main="Estimated K function")
plot(Kc,  iso - theo ~ r, main="Deviation from CSR")
abline(h=0, lty=2, col="gray40")
