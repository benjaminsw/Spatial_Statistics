#########################################
#### Analyse the Glasgow respiratory data
#########################################


################################
#### Read in and format the data
################################
#### Read the data in
dat <- read.csv(file="areal data - glasgow.csv")
head(dat)


#### Read in the shapefiles
library(sp)
library(rgdal)
shape <- readOGR(dsn = "SG_IntermediateZone_Bdry_2011.shp")
proj4string(shape)
shape2 <- spTransform(shape, CRS("+proj=longlat +datum=WGS84 +no_defs"))
proj4string(shape2)
head(shape2@data)


#### Merge the data and shapefile 
sp.dat <- merge(shape2, dat, all.x=FALSE, by.x="InterZone", by.y="IZ")
head(sp.dat@data)



##########################
#### Mapping using ggplot2
##########################
#### Turn the SpatialPolygonsDataFrame object into a data.frame object that contains the spaital information
library(ggplot2)
sp.dat@data$id <- rownames(sp.dat@data)
temp1 <- fortify(sp.dat, region = "id")
sp.dat2 <- merge(temp1, sp.dat@data, by = "id")



#### Change the colour scale
library(RColorBrewer)
ggplot(data = sp.dat2, aes(x=long, y=lat, goup=group, fill = smr)) + 
    geom_polygon() + 
    xlab("Longitude (degrees)") + 
    ylab("Latitude (degrees)") + 
    labs(title = "SMR for respiratory disease in 2015-16", fill = "SMR") +  
    theme(title = element_text(size=14)) + 
    scale_fill_gradientn(colors=brewer.pal(n=9, name="YlOrRd"))



#################################
#### Interactive map with leaflet
#################################
library(leaflet)
colours <- colorNumeric(palette = "PuRd", domain = sp.dat@data$smr, reverse=FALSE)
leaflet(data=sp.dat) %>% 
    addTiles() %>% 
    addPolygons(fillColor = ~colours(smr), 
                color="", weight=1, 
                fillOpacity = 0.6) %>%
    addLegend(pal = colours, values = sp.dat@data$smr, 
              opacity = 1, title="SMR") %>%
    addScaleBar(position="bottomleft")



####################################
#### Create the neighbourhood matrix
####################################
library(spdep)
W.nb <- poly2nb(sp.dat, row.names = rownames(sp.dat@data))
summary(W.nb)
W <- nb2mat(W.nb, style = "B")
class(W)
dim(W)
W[1:5, 1:5]

#### Plot the neighbourhood structure
plot(sp.dat)
plot.nb(W.nb, coords=coordinates(sp.dat), col="blue", add=TRUE)


######################
#### Compute Moran's I
######################
W.list <- nb2listw(W.nb, style = "B")
moran.mc(x = sp.dat@data$smr, listw = W.list, nsim = 10000)



####################################################
#### Fit a simple glm to assess residual correlation
####################################################
form <- Z~offset(log(E))+income + pm25 
model.glm <- glm(form, family="poisson", data=sp.dat@data)
summary(model.glm)



######################
#### Compute Moran's I
######################
moran.mc(x = residuals(model.glm), listw = W.list, nsim = 10000)



################################
#### Fit the Intrinsic CAR prior
################################
library(CARBayes)
model.car <- S.CARleroux(formula=form, family="poisson", data=sp.dat@data, W=W, rho=1, burnin=100000, n.sample=500000, thin=40)
print(model.car)
plot(model.car$samples$beta)
plot(model.car$samples$tau2)

## Relative risks - 10 unit increase for income
round(exp(10*model.car$summary.results[2, 1:3]),3)

## Relative risks - 1 unit increase for pm25
round(exp(1*model.car$summary.results[3, 1:3]),3)

################################
#### Fit the Leroux CAR prior
################################
model.car2 <- S.CARleroux(formula=form, family="poisson", data=sp.dat@data, W=W, burnin=100000, n.sample=500000, thin=40)
print(model.car2)
plot(model.car2$samples$beta)
plot(model.car2$samples$tau2)
plot(model.car2$samples$rho)

## Relative risks - 10 unit increase for income
round(exp(10*model.car2$summary.results[2, 1:3]),3)

## Relative risks - 1 unit increase for pm25
round(exp(1*model.car2$summary.results[3, 1:3]),3)