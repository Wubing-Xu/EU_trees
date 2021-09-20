rm(list= ls())

library(tidyverse)
library(raster)
library(rgeos)
library(maptools)
library(rgdal)
library(maps)
library(sf)


# tree species list
species <- dir("data/Occurrences",pattern = ".shp",full.names = FALSE)
species <- gsub(".shp", "", species)

# read the spatial points of all tree species
tree_points <-lapply(species, function(x) sf::st_read("data/Occurrences", x))
tree_points <- do.call(rbind, tree_points)

# the Lambert Azimuthal Equal Area projection in Europe
# https://gis.stackexchange.com/questions/182417/what-is-the-proj4string-of-the-lambert-azimuthal-equal-area-projection
laea <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs")

# the spatial polygons of European countries
data(wrld_simpl)
eu <- wrld_simpl[wrld_simpl@data$REGION == 150 & wrld_simpl@data$NAME != "Russia", ] #0   2   9  19 142 150
eu <- eu[eu$LAT < 75, ]
eu <- spTransform(eu, laea)
plot(eu, axes =TRUE)


# generate European raster in different resolutions
eu_raster <- raster(ext = extent(tree_points))
# 1 km
eu_raster_1km <- projectRaster(eu_raster, crs = laea, res = 1)
eu_raster_1km[] <- 1
eu_raster_1km <- mask(eu_raster_1km, eu)
# 5 km
eu_raster_5km <- projectRaster(eu_raster, crs = laea, res = 5)
eu_raster_5km[] <- 1
eu_raster_5km <- mask(eu_raster_5km, eu)
# 10 km
eu_raster_10km <- projectRaster(eu_raster, crs = laea, res = 10)
eu_raster_10km[] <- 1
eu_raster_10km <- mask(eu_raster_10km, eu)
# 20 km
eu_raster_20km <- projectRaster(eu_raster, crs = laea, res = 20)
eu_raster_20km[] <- 1
eu_raster_20km <- mask(eu_raster_20km, eu)

plot(eu_raster_20km)


# transform projection of tree distributions of points using equal-area projection
tree_points_latea <- st_transform(tree_points, laea)

# get the longitude and latitude of each occurrences and the transfomred coordinates X and Y
spid <- tibble(species = tree_points$Name, spid = 1:length(tree_points$Name))
tree_occs <- st_coordinates(tree_points) %>% 
  as_tibble() %>% 
  dplyr::select(longitude = X, latitude = Y) %>%
  bind_cols(st_coordinates(tree_points_latea) %>% as_tibble()) %>%
  left_join(spid, by = c("L1" = "spid")) %>%
  dplyr::select(species, longitude, latitude, X, Y) %>%
  mutate(species = gsub(" ", "_", species))

# get the ID of cells in different resolutions that occurrences located
tree_occs <- tree_occs %>%
  mutate(cell_1km = cellFromXY(eu_raster_1km, xy = as.matrix(tree_occs[, c("X", "Y")])),
         cell_5km = cellFromXY(eu_raster_5km, xy = as.matrix(tree_occs[, c("X", "Y")])),
         cell_10km = cellFromXY(eu_raster_10km, xy = as.matrix(tree_occs[, c("X", "Y")])),
         cell_20km = cellFromXY(eu_raster_20km, xy = as.matrix(tree_occs[, c("X", "Y")])))

# try to get the presence/absence of each species across grid cells. There are too many absences, requiring huge memory.
# decide to probide only presence
tree_occs %>%
  dplyr::select(species, cell_10km) %>%
  count(species, cell_10km) %>%
  pivot_wider(names_from = species, values_from = n, values_fill = 0) %>%
  pivot_longer(cols =-1, names_to = "species", values_to = "nocc") %>%
  mutate(presence = ifelse(nocc > 0, 1, 0))

# get the coordinates and number of occurrences of each species in grid-cells in different resolutions
# 1 km
tree_occs_1km <- tree_occs %>%
  dplyr::select(species, cell_1km) %>%
  count(species, cell_1km)
tree_occs_1km <- tree_occs_1km %>%
  bind_cols(xyFromCell(eu_raster_1km, tree_occs_1km$cell_1km) %>%
              as_tibble()) %>% 
  mutate(scale = 1) %>%
  dplyr::select(species, cellid = cell_1km, x, y, scale, nocc = n)

# 5 km
tree_occs_5km <- tree_occs %>%
  dplyr::select(species, cell_5km) %>%
  count(species, cell_5km)
tree_occs_5km <- tree_occs_5km %>%
  bind_cols(xyFromCell(eu_raster_5km, tree_occs_5km$cell_5km) %>%
              as_tibble()) %>% 
  mutate(scale = 5) %>%
  dplyr::select(species, cellid = cell_5km, x, y, scale, nocc = n)

# 10 km
tree_occs_10km <- tree_occs %>%
  dplyr::select(species, cell_10km) %>%
  count(species, cell_10km)
tree_occs_10km <- tree_occs_10km %>%
  bind_cols(xyFromCell(eu_raster_10km, tree_occs_10km$cell_10km) %>%
              as_tibble()) %>% 
  mutate(scale = 10) %>%
  dplyr::select(species, cellid = cell_10km, x, y, scale, nocc = n)

# 20 km
tree_occs_20km <- tree_occs %>%
  dplyr::select(species, cell_20km) %>%
  count(species, cell_20km)
tree_occs_20km <- tree_occs_20km %>%
  bind_cols(xyFromCell(eu_raster_20km, tree_occs_20km$cell_20km) %>%
              as_tibble()) %>% 
  mutate(scale = 20) %>%
  dplyr::select(species, cellid = cell_20km, x, y, scale, nocc = n)


# combine tree occurrences in different resolution
tree_occs_cells <- bind_rows(tree_occs_1km, tree_occs_5km, tree_occs_10km, tree_occs_20km)


save(tree_occs, tree_occs_cells, tree_points, eu, 
     eu_raster_1km, eu_raster_5km, eu_raster_10km, eu_raster_20km,
     file = "data/EU_tree_distributions.RDATA")


##############
# if require the spatial polygons, convert the raster of European coverage in different resolutions into spatial polygons
# the size of spatial polygons is quite big
# 1 km
eu_shape_1km <- rasterToPolygons(eu_raster_1km) 
colnames(eu_shape_1km@data)[1] <- "cellid"
eu_shape_1km$cellid <- which(as.vector(!is.na(eu_raster_1km)))

# 5 km
eu_shape_5km <- rasterToPolygons(eu_raster_5km) 
colnames(eu_shape_5km@data)[1] <- "cellid"
eu_shape_5km$cellid <- which(as.vector(!is.na(eu_raster_5km)))

# 10 km
eu_shape_10km <- rasterToPolygons(eu_raster_10km) 
colnames(eu_shape_10km@data)[1] <- "cellid"
eu_shape_10km$cellid <- which(as.vector(!is.na(eu_raster_10km)))

# 20 km
eu_shape_20km <- rasterToPolygons(eu_raster_20km) 
colnames(eu_shape_20km@data)[1] <- "cellid"
eu_shape_20km$cellid <- which(as.vector(!is.na(eu_raster_20km)))
