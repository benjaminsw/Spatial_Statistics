##############################
#### Beilschmiedia forest data
##############################

#### Load the data and library
library(spatstat)
data(bei)
plot(bei, pch="+",)
plot(bei.extra, main="")



#### Fit the model - covariates
fit <- ppm(bei~elev + grad, data=bei.extra)
summary(fit)
plot(fit)


fit2 <- ppm(bei~elev + grad + I(elev^2) + I(grad^2), data=bei.extra)
print(fit2)
plot(fit2)


#### Fit the model - spatial process
fit2 <- kppm(bei, trend=~elev + grad, data=bei.extra, clusters="LGCP", method="palm")
print(fit2)
plot(fit2)
summary(fit2)


par(mfrow=c(2,2))
plot(fit2)
plot(simulate(fit2))
plot(simulate(fit2))
plot(simulate(fit2), main="")



