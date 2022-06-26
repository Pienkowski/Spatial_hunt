################################################
###### Hunting project: Data exploration #######
################################################

# This script explores hunting data for the project XXX.

### It proceeds through the following steps:
# 1) Set environment 
# 2) Prepare the spatial data
#   2.1) Convert all data to UTM 
#   2.2) Preparing the land cover data 
#   2.3) Preparing the ridge and valley distance data 
#   2.4) Preparing the roads distance data
#   2.5) Preparing the travel time data 
# 3) Save all data for analysis


######### 1) Set environment ######### 

# Set working directory 
main.dir <- "C:/Users/wolf5246/Dropbox/Postdocs/Stirling CIFOR/Hunting_project/Hunting_analysis/Data_preperation"
setwd(main.dir)
### !!! Change file directory to wherever the data is stored !!! ###

### Load packages
library(dplyr)
library(sp)
library(ggplot2)
library(rgdal)
library(INLA)
library(inlabru)
library(RColorBrewer)
library(rgeos)
library(raster)
library(egg)
library(sf)

bru_options_set(inla.mode = "experimental")


#########  2) Prepare the spatial data ######### 

###### 2.1) Convert all data to UTM ######
# Trap data
Trap_Hunt_spatial <- readOGR("XXX/Input_layers/Trap_Hunt_spatial_rep.shp")
### !!! Change file directory to wherever the data is stored !!! ###

# Boundary area
Boundary <- readOGR("XXX/Input_layers/Boundary_rep.shp")
### !!! Change file directory to wherever the data is stored !!! ###

# Landcover
landcover <- raster("XXX/Input_layers/dense_dis.tif")
### !!! Change file directory to wherever the data is stored !!! ###

# Valley and ridge distance
valley_ridge_dis <- raster("XXX/Input_layers/vall_rid_dis.tif")

# Road distance
Road_dis <- raster("XXX/Input_layers/road_dis.tif")
### !!! Change file directory to wherever the data is stored !!! ###

# Resistance layer
landfeatures <- raster("XXX/Input_layers/land_mosaic.tif")
### !!! Change file directory to wherever the data is stored !!! ###

# Source points - communities and hunting camps
Sources <- readOGR("XXX/Input_layers/Source.shp")
### !!! Change file directory to wherever the data is stored !!! ###

# DEM
elevation <- raster("XXX/Input_layers/elev.tif")
### !!! Change file directory to wherever the data is stored !!! ###

# Save original CRS
orig_CRS <- crs(Trap_Hunt_spatial)

### Change map units from m to km 
Trap_Hunt_spatial <- spTransform(Trap_Hunt_spatial,"+proj=utm +zone=33 +south +datum=WGS84 +units=km +no_defs")

### Convert Boundary to UTM
Boundary <- spTransform(Boundary, CRSobj = crs(Trap_Hunt_spatial))

# Plot 
ggplot() +
  gg(Boundary) +
  gg(Trap_Hunt_spatial, color = "black", cex = 0.01) +
  coord_equal()

### Convert landcover to UTM
sr <- crs(Boundary)
landcover <- projectRaster(landcover, crs = sr,  method = "ngb")

# Plot 
plot(landcover)

### Convert valley/ridge distance to UTM
valley_ridge_dis <- projectRaster(valley_ridge_dis, crs = sr,  method = "ngb")

# Plot 
plot(valley_ridge_dis)

### Convert valley/ridge distance to UTM
Road_dis <- projectRaster(Road_dis, crs = sr,  method = "ngb")

# Plot 
plot(Road_dis)

### Convert resistance layer to UTM
landfeatures <- projectRaster(landfeatures, crs = sr,  method = "ngb")

# Plot 
plot(landfeatures)

### Convert Source points to UTM
Sources <- spTransform(Sources, CRSobj = crs(Trap_Hunt_spatial))

# Plot 
ggplot() +
  gg(Sources) +
  coord_equal()

### Convert DEM to UTM
elevation <- projectRaster(elevation, crs = orig_CRS,  method = "ngb")

# Plot 
plot(elevation)

######  2.2) Preparing the land cover data ###### 
# Convert to SpatialPixelsDataFrame
landcover <- as(landcover, "SpatialPixelsDataFrame")

# Plot 
ggplot() +
  gg(crop(landcover, Boundary)) +  
  gg(Trap_Hunt_spatial, color = "black", cex = 0.01) +
  coord_equal()


######  2.3) Preparing the ridge and valley distance data ###### 
# Convert to SpatialPixelsDataFrame
valley_ridge_dis <- as(valley_ridge_dis, "SpatialPixelsDataFrame")

# Plot 
ggplot() +
  gg(crop(valley_ridge_dis, Boundary)) +  scale_fill_distiller(palette = "clarity", direction=1) +
  gg(Trap_Hunt_spatial, color = "black", cex = 0.01) +
  coord_equal()


######  2.4) Preparing the roads distance data ###### 
# Convert to SpatialPixelsDataFrame
Road_dis <- as(Road_dis, "SpatialPixelsDataFrame")

# Plot 
ggplot() +
  gg(crop(Road_dis, Boundary))+ scale_fill_distiller(palette = "clarity", direction=1) +
  gg(Trap_Hunt_spatial, color = "black", cex = 0.01) +
  coord_equal()


######  2.5) Preparing the travel time data ###### 
### Time it takes to walk across one meter in minutes
road_sp <- 1/(4/60) # 4 kmph 
ridge_sp <- 1/(2/60)# 2 kmph 
valley_sp <- 1/(2/60) # 2 kmph 
bare_sp <- 1/(3/60) # 3 kmph 
savannah_sp <- 1/(3/60)# 3 kmph 
open_forest_sp <- 1/(2/60) # 2 kmph 
dense_forest_sp <- 1/(1/60) # 1 kmph 

### Create the slope layer### 
slope <- terrain(elevation, "slope", unit="radians", neighbors=8)
plot(slope, main="Slope")

# Re-project slope 
### Convert Source points to UTM
slope <- projectRaster(slope, crs = sr,  method = "ngb")
plot(slope, main="Slope")

###### Speed accross each type of terrain ###  
### Road resistence ### 
landfeatures <- reclassify(landfeatures, cbind(1, road_sp ))

### Ridge resistance ### 
landfeatures <- reclassify(landfeatures, cbind(2, ridge_sp ))

### Valley resistance ### 
landfeatures <- reclassify(landfeatures, cbind(3, valley_sp ))

### Bare land resistance ### 
landfeatures <- reclassify(landfeatures, cbind(4, bare_sp ))

### Savannah resistance ### 
landfeatures <- reclassify(landfeatures, cbind(6, savannah_sp ))

### Dense forest ### 
landfeatures <- reclassify(landfeatures, cbind(9, dense_forest_sp ))

### Open mixed forest ### 
landfeatures <- reclassify(landfeatures, cbind(10, open_forest_sp ))

### Plot
plot(landfeatures)

### Travel times 
table(landfeatures@data@values)

### Weighting each layer by slope ### 
# The  problem with this approach is that it treats walking uphill and downhill the same as we do not know the direction of approach.
slope_cost <- exp(-3.5*(abs(slope+0.05)))
names(slope_cost) <- "slope_cost"

# Land use multiplied by slope cost
cost <- landfeatures * slope_cost
plot(cost)

# For sense checking - crop the raster by the boundary of the study area
cost_boundary <- crop(cost, Boundary)

# Histogram of how long each cell takes to cross, weighted by slope
hist(cost_boundary@data@values)

# Plot resistance per pixel 
ggplot() +
  gg(cost_boundary) + scale_fill_distiller(palette = "clarity", direction=1) +
  gg(Trap_Hunt_spatial, color = "black", cex = 0.01) +
  coord_equal()


###### Create the travel time map ###### 
# Following: https://medium.com/@abertozz/mapping-travel-times-with-malariaatlas-and-friction-surfaces-f4960f584f08
### Create the transition layer using the Moore neighbourhood (eight orthogonal and diagonal nearest neighbours)
# The 1/mean in the below automatically converts resistance to concordance 
cost_trans <- gdistance::transition(cost, function(x) 1/mean(x), 8) 
cost_trans.GC <- gdistance::geoCorrection(cost_trans)   

# The accessibility algorithm needs a matrix of lat-longs
points <- as.matrix(Sources@coords)

# Calculate the accumulated cost surface
travel_time <- gdistance::accCost(cost_trans.GC, points)

# Travel times for each trap 
rasValue <- extract(travel_time, Trap_Hunt_spatial)

# Histogram of travel times to traps
hist(rasValue)

# Mean travel time to traps
mean(rasValue)

### Set the CRS for the layer
travel_time <- projectRaster(travel_time, crs = sr,  method = "ngb")

## Convert Inf values (on the border) to the mean of the surrounding cells
# Inf to arbitrary (but plausible travel times)
travel_time <- reclassify(travel_time, cbind(Inf, NA), right=FALSE)

# Save
writeRaster(travel_time, 'C:/Users/wolf5246/Dropbox/Postdocs/Stirling CIFOR/Hunting_project/Hunting_analysis/Spatial/Travel_time/travel_time.tif',options=c('TFW=YES'),  overwrite=TRUE)
### !!! Change file directory to wherever the data is stored !!! ###

# As SpatialPixelsDataFrame
travel_time <- as(travel_time, "SpatialPixelsDataFrame")

# First, scale travel time to hours 
travel_time$layer <- travel_time$layer/60


# Plot travel time 
ggplot() +
  gg(crop(travel_time, Boundary)) + scale_fill_distiller(palette = "clarity", direction=1) +
  gg(Trap_Hunt_spatial, color = "black", cex = 0.01) +
  coord_equal() 

### convert dense forest distance from meters to 100 meters 
landcover$dense_dis <- landcover$dense_dis/100
summary(landcover$dense_dis)

### convert valley distance from meters to 100 meters 
valley_ridge_dis$vall_rid_dis <- valley_ridge_dis$vall_rid_dis/100
summary(valley_ridge_dis$vall_rid_dis)

### convert valley distance from meters to 1 km 
Road_dis$road_dis <- Road_dis$road_dis/1000
summary(Road_dis$road_dis)


######### 3) Save all data for analysis #########
# Create a list containing all the spatial objects
trapping_data <- list(Boundary, Trap_Hunt_spatial, landcover, valley_ridge_dis, Road_dis, travel_time, Sources)
names(trapping_data) <- c("Boundary", "Trap_Hunt_spatial", "landcover", "valley_ridge_dis", "Road_dis", "travel_time", "Sources")

# Convert to Raster
landcover.R <- raster(trapping_data$landcover)
valley_ridge_dis.R <- raster(trapping_data$valley_ridge_dis)
Road_dis.R <- raster(trapping_data$Road_dis)
travel_time.R <- raster(trapping_data$travel_time)

# Resample travel time 
travel_time.R <- resample(travel_time.R, landcover.R, method = "bilinear")

# Make sure they are all the same extent 
ex = extent(landcover.R)

# Crop everything to travel time
travel_time.R <- crop(travel_time.R, ex)

# Mask 
travel_time.R = mask(travel_time.R, landcover.R)

# Converting directly from raster to spatial pixed drops NA values
# So, first converting to data frame then to 
travel_time.df = as.data.frame(travel_time.R,xy=TRUE); coordinates(travel_time.df)=~x+y
travel_time.df@proj4string <- sr
travel_time.t <- as(travel_time.df, "SpatialPixelsDataFrame")
trapping_data$travel_time <- travel_time.t

# Save the data
save(trapping_data, file = "trapping_data.Rda")


