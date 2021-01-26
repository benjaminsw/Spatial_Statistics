###############################
#### Soil contamination example
###############################
library(geoR)

## Read in the data and create a geoR object
soil.temp1 <- read.csv(file="geostat data 2 - soil.csv")
head(soil.temp1)
soil.temp2 <- read.csv(file="geostat data 2 - soil borders.csv")
head(soil.temp2)
soil <- as.geodata(soil.temp1, coords.col=1:2, data.col=3, covar.col=4:5, borders=TRUE)
soil$borders <- soil.temp2


## Plot the data
plot(soil, lowess=TRUE)


## Estimate covariate effects and re-plot
model <- lm(calcium~factor(region)+elevation, data=soil.temp1)
summary(model)

model <- lm(calcium~factor(region), data=soil.temp1)
summary(model)


#### Construct a geodata object of the residuals
resid <- residuals(model)
soil.temp1$resid <-resid 
soil2 <- as.geodata(soil.temp1, coords.col=1:2, data.col=6, covar.col=4:5, borders=TRUE)
soil2$borders <- soil.temp2


#### Plot the variogram of the residuals
vari <- variog(soil2)
vari.mc <- variog.mc.env(soil2, obj.variog=vari)
plot(vari, envelope.obj=vari.mc)


## Directional variogram
vari <- variog4(soil2)
plot(vari, envelope.obj=vari.mc)


## Estimate the parameters using maximum likelihood
model <- likfit(geodata=soil, trend=~ factor(soil$covariate$region), ini.cov.pars=variofit(variog(soil2), cov.model="exponential"), fix.nugget=FALSE, cov.model="exponential")
summary(model)


## Construct 95% confidence intervals
SE <- sqrt(diag(model$beta.var))
round(cbind(model$beta - qt(0.975, 174)*SE, model$beta + qt(0.975, 174)*SE),3)


## Create the uncorrelated innovations and assess model suitability.
fitted <- as.numeric(trend.spatial(trend=~factor(soil$covariate$region), geodata=soil) %*% model$beta)
resid <- soil$data - fitted
Sigma.est <- varcov.spatial(soil$coords, cov.model="exponential", cov.pars=model$cov.pars, nugget=model$nugget)$varcov
innovation <- as.numeric(solve(chol(Sigma.est)) %*% resid)
innovation.geo <- soil
innovation.geo$data <- innovation


## Conduct residual checks
par(mfrow=c(1,2))
qqnorm(innovation)
vari <- variog(innovation.geo)
vari.mc <- variog.mc.env(innovation.geo, obj.variog=vari)
plot(vari, envelope.obj=vari.mc)



## Comparing different covariance models
model1 <- likfit(geodata=soil, trend=~factor(soil$covariate$region), ini.cov.pars=variofit(variog(soil2), cov.model="exponential"), fix.nugget=FALSE, cov.model="exponential")
model2 <- likfit(geodata=soil, trend=~factor(soil$covariate$region), ini.cov.pars=variofit(variog(soil2), cov.model="gaussian"), fix.nugget=FALSE, cov.model="gaussian")
model3 <- likfit(geodata=soil, trend=~factor(soil$covariate$region), ini.cov.pars=variofit(variog(soil2), cov.model="spherical"), fix.nugget=FALSE, cov.model="spherical")
AIC(model1)
AIC(model2)
AIC(model3)




## Predict the surface
x.grid <- seq(4900, 6000, 10)
y.grid <- seq(4800,6000, 10)
predgrid <- expand.grid(x=x.grid,  y=y.grid)
par(mfrow=c(1,1))
plot(predgrid)

## Subset the surface to only the locations inside the boundary
library(mgcv)
which.in <- in.out(bnd=as.matrix(soil$borders), x=as.matrix(predgrid))
predgrid2 <- predgrid[which.in, ]


## Plot the predicton locations
plot(predgrid2, pch=19, col="red", main="Prediction and data locations", xlim=c(4900, 6500))
points(soil$coords, col="blue", pch=19)
lines(soil$borders)


## Do the Kriging
model <- likfit(geodata=soil, trend="cte", ini.cov.pars=variofit(variog(soil), cov.model="exponential"), fix.nugget=FALSE, cov.model="exponential")
control <- krige.control(type.krige="OK", trend.d="cte", trend.l="cte", obj.model=model)
kriging <- krige.conv(geodata=soil, locations=predgrid2, krige=control)
summary(kriging)


## Create a data frame of predictions
pred.data <- data.frame(x=predgrid2[ ,1], y=predgrid2[ ,2], predictions=kriging$predict, sd=sqrt(kriging$krige.var))


## Map the predicted values
library(ggplot2)
library(RColorBrewer)
ggplot(aes(x = x, y = y), data = pred.data) + 
    geom_tile(aes(fill = predictions)) + 
    coord_equal() + 
    labs(title = "Predicted nitrogen levels", fill = "Nitrogen") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="YlGnBu"))


## Map the prediction standard deviations
ggplot(aes(x = x, y = y), data = pred.data) + 
    geom_tile(aes(fill = sd)) + 
    coord_equal() + 
    labs(title = "Prediction standard deviation", fill = "Nitrogen") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="YlGnBu"))





#### Model comparison via cross validation - Exponential vs Gaussian vs Spherical
## Set up the results matrix
m <- 30
results <- array(NA, c(m,3))
colnames(results) <- c("Data", "Exponential", "Spherical")
results[ ,1] <- soil$data[1:m]


## Estimate initial covariance parameters
inits.exp <- variofit(variog(soil), cov.model="exponential")
inits.sph <- variofit(variog(soil), cov.model="spherical")


## Run the LOOCV
     for(i in 1:m)
     {
     ## Create the data
     soil.data.temp <- as.geodata(soil.temp1[-i, ], coords.col=1:2, data.col=3, covar.col=4:5)
     
     ## Exponential
     model.exp <- likfit(geodata=soil.data.temp, trend="cte", ini.cov.pars=inits.exp, fix.nugget=FALSE, cov.model="exponential")
     control <- krige.control(type.krige="OK", trend.d="cte", trend.l="cte", obj.model=model.exp)
     kriging <- krige.conv(geodata=soil.data.temp, locations=c(soil.temp1[i, 1:2]), krige=control)
     results[i, 2] <- kriging$predict

     ## Spherical
     model.sph <- likfit(geodata=soil.data.temp, trend="cte", ini.cov.pars=inits.sph, fix.nugget=FALSE, cov.model="spherical")
     control <- krige.control(type.krige="OK", trend.d="cte", trend.l="cte", obj.model=model.sph)
     kriging <- krige.conv(geodata=soil.data.temp, locations=c(soil.temp1[i, 1:2]), krige=control)
     results[i, 3] <- kriging$predict 
     }

## Compute the RMSE
sqrt(mean((results[ ,1] - results[ ,2])^2))
sqrt(mean((results[ ,1] - results[ ,3])^2))

