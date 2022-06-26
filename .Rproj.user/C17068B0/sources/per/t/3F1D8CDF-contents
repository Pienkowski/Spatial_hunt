################################################
###### Hunting project: Data analysis #######
################################################

# This script explores hunting data for the project XXX.

### It proceeds through the following steps:
# 0) Set environment 
# 1) Check correlation between rasters 
# 2) Main model
# 3) Age model 

######### 0) Set environment ######### 

# Set working directory 
main.dir <- "XXX/Data_preperation"
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
library(usdm)

# # Specify experimental BRU options 
# bru_options_set(inla.mode = "experimental")

# Load the data
load("trapping_data.Rda")

# Create the individual objects 
mesh <- trapping_data$mesh
Boundary <- trapping_data$Boundary
Trap_Hunt_spatial <- trapping_data$Trap_Hunt_spatial
landcover <- trapping_data$landcover
valley_ridge_dis <- trapping_data$valley_ridge_dis
Road_dis <- trapping_data$Road_dis
travel_time <- trapping_data$travel_time
Sources <- trapping_data$Sources

######### 1) Check correlation between rasters ######### 

### Convert to raster object 
landcover_ras <- raster(landcover)
valley_ridge_dis_ras <- raster(valley_ridge_dis)
Road_dis_ras <- raster(Road_dis)
travel_time_ras <- raster(travel_time)

# Create a raster stack 
rast_stack <- stack(list(landcover_ras, valley_ridge_dis_ras, Road_dis_ras, travel_time_ras))

# Pearson correlation between rasters
jnk=layerStats(rast_stack, 'pearson', na.rm=T)
corr_matrix=jnk$'pearson correlation coefficient'

# Examine the correlation - all seem pretty modest
corr_matrix

# Collinearity - all seem pretty modest
usdm::vif(rast_stack)

######### 2) Main model######### 
# Define the matrix
mesh <- inla.mesh.2d(max.edge = c(0.5, 1),
                     cutoff= 0.2,
                     boundary = Boundary,
                     min.angle = c(33, 25),
                     crs = wkt(obj = Boundary))

# Plot travel time 
ggplot() +
  gg(mesh) +
  coord_equal() 

# Save the mesh 
save(mesh, file = "mesh.Rda")

# Matern 
matern <- inla.spde2.matern(mesh)

### Run each model
# Model 1
start_time <- Sys.time()
mod_1 <- coordinates ~ dense_dis(landcover, model = "linear")  + mySmooth(coordinates, model = matern)  
mod_1_fit <- lgcp(mod_1, Trap_Hunt_spatial, samplers = Boundary, domain = list(coordinates = mesh))
end_time <- Sys.time()
end_time - start_time
save(mod_1_fit, file = "mod_1_fit.Rda")

# Model 2
start_time <- Sys.time()
mod_2 <- coordinates ~ layer(travel_time, model = "linear") + mySmooth(coordinates, model = matern)  
mod_2_fit <- lgcp(mod_2, Trap_Hunt_spatial, samplers = Boundary, domain = list(coordinates = mesh))
end_time <- Sys.time()
end_time - start_time
save(mod_2_fit, file = "mod_2_fit.Rda")

# Model 3
start_time <- Sys.time()
mod_3 <- coordinates ~ dense_dis(landcover, model = "linear") + layer(travel_time, model = "linear") + mySmooth(coordinates, model = matern)  
mod_3_fit <- lgcp(mod_3, Trap_Hunt_spatial, samplers = Boundary, domain = list(coordinates = mesh))
end_time <- Sys.time()
end_time - start_time
save(mod_3_fit, file = "mod_3_fit.Rda")

# Model 4
start_time <- Sys.time()
mod_4 <- coordinates ~ dense_dis(landcover, model = "linear") + layer(travel_time, model = "linear") + vall_rid_dis(valley_ridge_dis, model = "linear") +  mySmooth(coordinates, model = matern)  
mod_4_fit <- lgcp(mod_4, Trap_Hunt_spatial, samplers = Boundary, domain = list(coordinates = mesh))
end_time <- Sys.time()
end_time - start_time
save(mod_4_fit, file = "mod_4_fit.Rda")

# Model 5
start_time <- Sys.time()
mod_5 <- coordinates ~ dense_dis(landcover, model = "linear") + layer(travel_time, model = "linear") + vall_rid_dis(valley_ridge_dis, model = "linear") + road_dis(Road_dis, model = "linear") + mySmooth(coordinates, model = matern)  
mod_5_fit <- lgcp(mod_5, Trap_Hunt_spatial, samplers = Boundary, domain = list(coordinates = mesh))
end_time <- Sys.time()
end_time - start_time
save(mod_5_fit, file = "mod_5_fit.Rda")

# Create list of fitted models 
fit_list <- list(mod_1_fit, mod_2_fit, mod_3_fit, mod_4_fit, mod_5_fit)

# # Save 
save(fit_list, file = "fit_list.Rda")

# Model AIC
model_aic <- lapply(fit_list, deltaIC)
model_aic_df <- do.call("rbind", model_aic)

# Save optimal model 
optimal_mod <- fit_list[[4]]
save(optimal_mod, file = "optimal_mod.Rda")


######### 3) Age model #########
### Calculate the mean age of hunters

# Subset to hunters
Trap_Hunt_data_df <- data.frame(Trap_Hunt_spatial)
Trap_Hunt_data_df <- Trap_Hunt_data_df %>% 
  group_by(Huntrnm) %>%
  filter(row_number()==1)

# mean age 
age_mean <- mean(Trap_Hunt_data_df$Age)

# Split the sample into younger and older respondents, based on the median
Trap_Hunt_spatial_older <- Trap_Hunt_spatial[Trap_Hunt_spatial$Age > age_mean ,]
Trap_Hunt_spatial_younger <- Trap_Hunt_spatial[Trap_Hunt_spatial$Age <= age_mean,]


### Run 

# Model 4 - older 
start_time <- Sys.time()
mod_4 <- coordinates ~ dense_dis(landcover, model = "linear") + layer(travel_time, model = "linear") + vall_rid_dis(valley_ridge_dis, model = "linear") +  mySmooth(coordinates, model = matern)  
mod_4_older_fit <- lgcp(mod_4, Trap_Hunt_spatial, samplers = Boundary, domain = list(coordinates = mesh))
end_time <- Sys.time()
end_time - start_time
save(mod_4_older_fit, file = "mod_4_older_fit.Rda")

# Model 4 - younger 
start_time <- Sys.time()
mod_4 <- coordinates ~ dense_dis(landcover, model = "linear") + layer(travel_time, model = "linear") + vall_rid_dis(valley_ridge_dis, model = "linear") +  mySmooth(coordinates, model = matern)  
mod_4_younger_fit <- lgcp(mod_4, Trap_Hunt_spatial, samplers = Boundary, domain = list(coordinates = mesh))
end_time <- Sys.time()
end_time - start_time
save(mod_4_younger_fit, file = "mod_4_younger_fit.Rda")

