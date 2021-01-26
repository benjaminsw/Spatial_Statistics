####################
#### Glasgow example
####################

#### Read in the data
dat <- read.csv(file="glasgowdata.csv")
head(dat)


#### Compute the SIR
dat$sir <- dat$Y / dat$E


#### Read in the shapefiles
library(sp)
library(rgdal)
shape <- readOGR(dsn = "SG.shp")
class(shape)
head(shape@data)

#### Change the coordinates
proj4string(shape)
shape2 <- spTransform(shape, CRS("+proj=longlat +datum=WGS84 +no_defs"))
proj4string(shape2)


#### Merge the data and shapefile
sp.dat <- merge(shape2, dat, all.x=FALSE, by.x="IZ_CODE", by.y="IG")
class(sp.dat)
head(sp.dat@data)
plot(sp.dat)

#### Map the SIR
spplot(sp.dat, "sir")


#### Load mapping libraries
library(ggplot2)
library(rgeos)
library(maptools)

#### Turn it back into a data.frame for mapping
sp.dat@data$id <- rownames(sp.dat@data)
temp1 <- fortify(sp.dat, region = "id")
sp.dat2 <- merge(temp1, sp.dat@data, by = "id")


#### A simple map
ggplot(data = sp.dat2, aes(x=long, y=lat, goup=group, fill = sir)) + 
    geom_polygon()



#### A nicer colour scale
library(RColorBrewer)
ggplot(data = sp.dat2, aes(x=long, y=lat, goup=group, fill = sir)) + 
    geom_polygon() + 
    xlab("Longitude") + 
    ylab("Latitude") + 
    labs(title = "SIR for cancer risk", fill = "SIR") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="Reds"))


#### Interactive mapping
library(leaflet)
colours <- colorNumeric(palette = "OrRd", domain = sp.dat@data$sir, reverse=FALSE)
leaflet(data=sp.dat) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(sir), 
                color="black", weight=1, 
                fillOpacity = 0.7) %>%
    addLegend(pal = colours, values = sp.dat@data$sir, 
              opacity = 1, title="SIR") %>%
    addScaleBar(position="bottomleft")


#### Create the spatial list information for constructing Moran's I
library(spdep)
W.nb <- poly2nb(sp.dat, row.names = rownames(sp.dat@data))
W.list <- nb2listw(W.nb, style = "B")


#### Compute Moran's I
moran.mc(x = sp.dat@data$sir, listw = W.list, nsim = 10000)


#### Fit a GLM model
form <- Y~offset(log(E))+pm10+smoke+ethnic
model1 <- glm(formula=form, family=poisson, data=dat)
summary(model1)


#### Check the residuals for correlation
moran.mc(x = residuals(model1), listw = W.list, nsim = 10000)


#### Compute the neighbourhood matrix
W <- nb2mat(W.nb, style = "B")


#### Fit a spatial correlation model
library(CARBayes)
model2 <- S.CARleroux(formula=form, family="poisson", data=sp.dat@data, W=W, 
                      burnin=20000, n.sample=120000, thin=10,  verbose=TRUE)
print(model2)


#### View model components
summary(model2)
summary(model2$samples)


#### Check convergence
plot(model2$samples$rho)


#### Covariate effects
exp(model2$summary.results[2:4 , 1:3])

summary(dat$pm10)
summary(dat$smoke)

exp(sd(dat$pm10) * model2$summary.results[2 , 1:3])
exp(sd(dat$smoke) * model2$summary.results[3 , 1:3])



#### What happens if you fit the intrinsic model?
model3 <- S.CARleroux(formula=form, family="poisson", data=sp.dat@data, W=W, 
                      burnin=20000, n.sample=120000, thin=10,  rho=1, verbose=TRUE)
print(model3)

