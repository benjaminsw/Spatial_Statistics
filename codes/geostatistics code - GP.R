#######################################################
#### Read in and visualise the GP data on a leaflet map
#######################################################


setwd("~/OneDrive - University of Glasgow/teaching/Spatial statistics 2021/Code and data")


library(sp)
library(leaflet)
library(RColorBrewer)


################################
#### Read in and format the data
################################
#### Read the data in
dat <- read.csv(file="geostat data 3 GP.csv")
head(dat)


#### Format as a spatialpointsdataframe
sp.dat <- SpatialPointsDataFrame(coords=dat[ ,2:3], data=dat)
proj4string(sp.dat) <- CRS("+proj=longlat +datum=WGS84 +no_defs")


#### Map the data
colours <- colorNumeric(palette = "YlOrRd", domain = sp.dat@data$SPR)
leaflet(data=sp.dat) %>% 
    addTiles() %>% 
    addCircles(~lon, ~lat, color = ~colours(SPR), opacity = 1, radius=90) %>%
    addLegend(pal = colours, values = sp.dat@data$SPR, 
              opacity = 1, title="SPR")

