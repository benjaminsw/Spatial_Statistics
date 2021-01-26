###########################
#### Chesapeake Bay example
###########################


#### Read in the
data1 <- read.csv(file="chesbay.csv")
head(data1)
data2 <- read.csv(file="chesbayboundary.csv")
head(data2)


#### Create a geodata object
library(geoR)
data.geodata <- as.geodata(data1, coords.col=2:3, data.col=4, borders=TRUE)
data.geodata$borders <- data2


#### Plot the geodata object
plot(data.geodata)
## The data have a prominant north/south trend which is non-linear.
## No east/west trend though, although this direction is very narrow.
## The data also appear to be skewed from the histogram.



#### Leaflet mapping
library(leaflet)
library(RColorBrewer)
library(sp)

sp.dat <- SpatialPointsDataFrame(coords=data1[ ,2:3], data=data1)
proj4string(sp.dat) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

colours <- colorNumeric(palette = "YlOrRd", domain = sp.dat@data$nitrogen, reverse=FALSE)
leaflet(data=sp.dat) %>% 
    addTiles() %>% 
    addCircles(~longitude, ~latitude, color = ~colours(nitrogen), opacity = 1, radius=80) %>%
    addLegend(pal = colours, values = sp.dat@data$nitrogen, 
              opacity = 1, title="Nitrogen")


#### Assess a suitable linear trend model
model1 <- lm(nitrogen~longitude+latitude, data=data1)
summary(model1)
qqnorm(residuals(model1))
## High R2 but the residuals look non-normal. Given the first histogram


#### try a log transform of nitrogen as the response.
model2 <- lm(log(nitrogen)~longitude+latitude, data=data1)
summary(model2)
qqnorm(residuals(model2))
## This looks better, clear latitude effect but no longitude effect.
## The latter will be left in to provide a rotationally invariant surface.


#### Create a new residual geodata set
resid.dat <- data.frame(residuals=residuals(model2), longitude=data1$longitude, latitude=data1$latitude)
resid.geodata <- as.geodata(resid.dat, coords.col=2:3, data.col=1, borders=TRUE)
resid.geodata$borders <- data2


#### Examine the presence of residual spatial correlation
vari <- variog(resid.geodata)
vari.mc <- variog.mc.env(resid.geodata, obj.variog=vari)
plot(vari, envelope.obj=vari.mc)
## The variogram as always is messy, and shows some evidence of correlation
## as some of the points are outside the monte carlo envelopes.


#### Fit an exponential model to the data
model3 <- likfit(geodata=data.geodata, trend="1st", ini.cov.pars=c(0.04, 0.5), fix.nugget=FALSE, lambda=0, cov.model="exponential")
summary(model3)


#### Set up a regular prediction grid
predgrid <- expand.grid(x=seq(-77.5, -75, 0.01), y=seq(36.8, 40.2, 0.01))


#### Subset it to only locations within the boundary
library(mgcv)
which.in <- in.out(bnd=as.matrix(data.geodata$borders),x=as.matrix(predgrid))
predgrid2 <- predgrid[which.in, ]


#### Plot the prediction locations
plot(data.geodata$borders, type="l")
points(predgrid2, pch=19, col="red", cex=0.1)


#### Make the predictions
control <- krige.control(type.krige="OK", trend.d="1st", trend.l="1st", obj.model=model3, lambda=0)
kriging <- krige.conv(geodata=data.geodata, locations=predgrid2, krige=control)


#### Make a prediction data frame
pred.data <- data.frame(x=predgrid2[ ,1], y=predgrid2[ ,2], predictions=kriging$predict, sd=sqrt(kriging$krige.var))


#### Simple graph
library(ggplot2)
ggplot(aes(x = x, y = y), data = pred.data) + 
    geom_tile(aes(fill = predictions)) 


#### Change the colour scale
library(RColorBrewer)
ggplot(aes(x = x, y = y), data = pred.data) + 
    geom_tile(aes(fill = predictions)) + 
    coord_equal() + 
    xlab("Longitude (degrees)") + 
    ylab("Latitude (degrees)") + 
    labs(title = "Predicted nitrogen levels", fill = "Nitrogen") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="YlGnBu"))


#### Map the standard deviation
ggplot(aes(x = x, y = y), data = pred.data) + 
    geom_tile(aes(fill = sd)) + 
    coord_equal() + 
    xlab("Longitude (degrees)") + 
    ylab("Latitude (degrees)") + 
    labs(title = "Predictive standard deviations", fill = "Uncertainty") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="RdPu"))




#### Cross validation to compare exponential and spherical models
## Set up the results matrix
m <- nrow(data1)
results <- array(NA, c(m,3))
colnames(results) <- c("data", "exponential", "spherical")
results[ ,1] <- data1$nitrogen

## Run the predictions
for(i in 1:m)
{
    ## Set up the test data
    test <- as.geodata(data1[-i, ], coords.col=2:3, data.col=4, borders=TRUE)
    
    ## Run the exponential model
    model.exp <- likfit(geodata=test, trend="1st", ini.cov.pars=c(0.04, 0.5), fix.nugget=FALSE, lambda=0, cov.model="exponential")
    control <- krige.control(type.krige="OK", trend.d="1st", trend.l="1st", obj.model=model.exp, lambda=0)
    kriging <- krige.conv(geodata=test, locations=data1[i, 2:3], krige=control)
    results[i, 2] <- kriging$predict 
    
    ## Run the spherical model
    model.sph <- likfit(geodata=test, trend="1st", ini.cov.pars=c(0.04, 0.5), fix.nugget=FALSE, lambda=0, cov.model="spherical")
    control <- krige.control(type.krige="OK", trend.d="1st", trend.l="1st", obj.model=model.sph, lambda=0)
    kriging <- krige.conv(geodata=test, locations=data1[i, 2:3], krige=control)
    results[i, 3] <- kriging$predict 
}

## Compute the RMSE
sqrt(mean((results[ ,1] - results[ ,2])^2))
sqrt(mean((results[ ,1] - results[ ,3])^2))


