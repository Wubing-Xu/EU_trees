rm(list= ls())

library(tidyverse)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
library(maps)
library(sf)


load("data/EU_tree_distributions.RDATA")

# read the raster of human footprint 2009
hfp <- stack("data/HFP2009/HFP2009.tif")

# read the raster of bioclimatic variables from worldclim
# as the size of the files is very big, the rasters were firstly croped into europe and the keep the cropped in folders
# climate_dir <-paste("F:/Study/DATA/Environments/WorldClim/wc2.1_30s_bio/wc2.1_30s_bio_", 1:19, ".tif", sep="")
# climate <- stack(climate_dir)
# crop the global map into the Europe
# climate_crop <- crop(x = climate, y = extent(st_transform(tree_points, projection(climate)))) 
# writeRaster(climate_crop, file = "data/wc2.1_30s_bio/climate_crop.tif")

climate_crop <- stack("data/wc2.1_30s_bio/climate_crop.tif")
# climate_crop <- readAll(climate_crop)


# project the raster of HPF into the raster in different resolutions that used to determine tree distributions 
# crop the global map into the Europe
hfp_crop <- crop(x = hfp, y = extent(st_transform(tree_points, projection(hfp)))) 
# 1 km
hfp_1km <- projectRaster(from = hfp_crop, to = eu_raster_1km)
hfp_1km[hfp_1km<0] <- 0
# 5 km
hfp_5km <- raster::aggregate(hfp_crop, 5)
hfp_5km <- projectRaster(from = hfp_5km, to = eu_raster_5km)
hfp_5km[hfp_5km<0] <- 0
# 10 km
hfp_10km <- raster::aggregate(hfp_crop, 10)
hfp_10km <- projectRaster(from = hfp_10km, to = eu_raster_10km)
hfp_10km[hfp_10km<0] <- 0
# 20 km
hfp_20km <- raster::aggregate(hfp_crop, 20)
hfp_20km <- projectRaster(from = hfp_20km, to = eu_raster_20km)
hfp_20km[hfp_20km<0] <- 0

# get the HFP for each occupied gird-cells
hfp_cells <- c(hfp_1km[tree_occs_cells %>% filter(scale == 1) %>% pull(cellid)],
               hfp_5km[tree_occs_cells %>% filter(scale == 5) %>% pull(cellid)],
               hfp_10km[tree_occs_cells %>% filter(scale == 10) %>% pull(cellid)],
               hfp_20km[tree_occs_cells %>% filter(scale == 20) %>% pull(cellid)])


# project the raster of bioclimatic variables into the raster in different resolutions that used to determine tree distributions 
# 1 km
climate_1km <- projectRaster(from = climate_crop, to = eu_raster_1km)
climate_1km[climate_1km<0] <- 0
# 5 km
climate_5km <- raster::aggregate(climate_crop, 5)
climate_5km <- projectRaster(from = climate_5km, to = eu_raster_5km)
climate_5km[climate_5km<0] <- 0
# 10 km
climate_10km <- raster::aggregate(climate_crop, 10)
climate_10km <- projectRaster(from = climate_10km, to = eu_raster_10km)
climate_10km[climate_10km<0] <- 0
# 20 km
climate_20km <- raster::aggregate(climate_crop, 20)
climate_20km <- projectRaster(from = climate_20km, to = eu_raster_20km)
climate_20km[climate_20km<0] <- 0

# get the climate for each occupied gird-cells
climate_cells <- rbind(climate_1km[tree_occs_cells %>% filter(scale == 1) %>% pull(cellid)],
               climate_5km[tree_occs_cells %>% filter(scale == 5) %>% pull(cellid)],
               climate_10km[tree_occs_cells %>% filter(scale == 10) %>% pull(cellid)],
               climate_20km[tree_occs_cells %>% filter(scale == 20) %>% pull(cellid)])

colnames(climate_cells)<- paste("bio_", 1:19, sep = "")


# save(hfp_cells, climate_cells, file = "data/environment_cells.RDATA")


# combine environments with occurrences
tree_occs_cells <- tree_occs_cells %>% 
  mutate(hfp = hfp_cells) %>%
  bind_cols(climate_cells %>% as_tibble())


save(tree_occs, tree_occs_cells, tree_points, eu, 
     eu_raster_1km, eu_raster_5km, eu_raster_10km, eu_raster_20km,
     file = "data/EU_tree_distributions_environments.RDATA")
