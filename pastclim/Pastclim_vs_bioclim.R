# 1. download the pastclim data (The Aohan climate layers can be found in the here("pastclim") folder.)
# Install and load required packages
install.packages("devtools")
devtools::install_github("EvolEcolGroup/pastclim")
devtools::install_github("EvolEcolGroup/pastclim", build_vignettes = TRUE)
library(terra)
library(pastclim)
library(here)
# Load vignettes
vignette("pastclim_overview", package = "pastclim")
vignette("custom_datasets", package = "pastclim")
vignette("available_datasets", package = "pastclim")
citation("pastclim")
# Set data path
set_data_path(path_to_nc = here("pastclim", "my_reconstructions"))
# Get dataset variables
get_vars_for_dataset(dataset = "Beyer2020")
# Read vector file
Aohan_vect <- terra::vect(read_sf(here('pastclim', 'aohanbianjie.shp')))

# Beyer 2020
# Annual precipitation
download_dataset(dataset = "Beyer2020", bio_variables = "bio12")
# BP4k
climate_BP4k <- region_slice(time_bp = -4000, bio_variables = "bio12", dataset = "Beyer2020", crop = Aohan_vect)
terra::writeRaster(climate_BP4k, here("pastclim", "climate_BP4k.tif"), overwrite=TRUE)
# BP3k
climate_BP3k <- region_slice(time_bp = -3000, bio_variables = "bio12", dataset = "Beyer2020", crop = Aohan_vect)
terra::writeRaster(climate_BP3k, here("pastclim", "climate_BP3k.tif"), overwrite=TRUE)
# Annual mean temperature
download_dataset(dataset = "Beyer2020", bio_variables = "bio01")
# BP4k
climate_BP4k_meantem <- region_slice(time_bp = -4000, bio_variables = "bio01", dataset = "Beyer2020", crop = Aohan_vect)
terra::writeRaster(climate_BP4k_meantem, here("pastclim", "climate_BP4k_meantem.tif"), overwrite=TRUE)
# BP3k
climate_BP3k_meantem <- region_slice(time_bp = -3000, bio_variables = "bio01", dataset = "Beyer2020", crop = Aohan_vect)
terra::writeRaster(climate_BP3k_meantem, here("pastclim", "climate_BP3k_meantem.tif"), overwrite=TRUE)

# Krapp2021
get_vars_for_dataset(dataset = "Krapp2021")
# BC2000-1400/BP3950-3350
Aohan_vect <- terra::vect(read_sf(here('pastclim', 'aohanbianjie.shp')))
# Annual precipitation
download_dataset(dataset = "Krapp2021", bio_variables = "bio12")
# BP4k
climate_BP4k_Krapp <- region_slice(time_bp = -4000, bio_variables = "bio12", dataset = "Krapp2021", crop = Aohan_vect)
terra::writeRaster(climate_BP4k_Krapp, here("pastclim", "climate_BP4k_Krapp.tif"), overwrite=TRUE)
# BP3k
climate_BP3k_Krapp <- region_slice(time_bp = -3000, bio_variables = "bio12", dataset = "Krapp2021", crop = Aohan_vect)
terra::writeRaster(climate_BP3k_Krapp, here("pastclim", "climate_BP3k_Krapp.tif"), overwrite=TRUE)
# Annual mean temperature
download_dataset(dataset = "Krapp2021", bio_variables = "bio01")
# BP4k
climate_BP4k_meantem_Krapp <- region_slice(time_bp = -4000, bio_variables = "bio01", dataset = "Krapp2021", crop = Aohan_vect)
terra::writeRaster(climate_BP4k_meantem_Krapp, here("pastclim", "climate_BP4k_meantem_Krapp.tif"), overwrite=TRUE)
# BP3k
climate_BP3k_meantem_Krapp <- region_slice(time_bp = -3000, bio_variables = "bio01", dataset = "Krapp2021", crop = Aohan_vect)
terra::writeRaster(climate_BP3k_meantem_Krapp, here("pastclim", "climate_BP3k_meantem_Krapp.tif"), overwrite=TRUE)

# 2. read and prepare the Aohan climate layers
# Beyer 2020
# Annual precipitation
# BP4k
climate_BP4k_raster <- raster(here("pastclim", 'climate_BP4k.tif'), window=Aohan_owin)
climate_BP4k_ras_proj <- projectRaster(climate_BP4k_raster, crs=crs(DEM))
summary(climate_BP4k_ras_proj)
#            bio12
# Min.    406.1114
# 1st Qu. 418.2406
# Median  419.3803
# 3rd Qu. 427.2933
# Max.    434.4400
# NA's      0.0000
res(climate_BP4k_ras_proj)
# [1] 41200 55500
jpeg(here("pastclim", "plot", "Past_precipitation_Beyer_BP4k.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(climate_BP4k_ras_proj, main = NULL)
dev.off()

# BP3k
climate_BP3k_raster <- raster(here("pastclim", 'climate_BP3k.tif'), window=Aohan_owin)
climate_BP3k_ras_proj <- projectRaster(climate_BP3k_raster, crs=crs(DEM))
summary(climate_BP3k_ras_proj)
#            bio12
# Min.    407.6532
# 1st Qu. 422.2110
# Median  428.0567
# 3rd Qu. 440.0047
# Max.    447.4967
# NA's      0.0000
res(climate_BP3k_ras_proj)
# [1] 41200 55500
jpeg(here("pastclim", "plot", "Past_precipitation_Beyer_BP3k.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(climate_BP3k_ras_proj, main = NULL)
dev.off()


# Annual mean temperature
# BP4k
climate_BP4k_meantem_raster <- raster(here("pastclim", 'climate_BP4k_meantem.tif'), window=Aohan_owin)
climate_BP4k_meantem_ras_proj <- projectRaster(climate_BP4k_meantem_raster, crs=crs(DEM))
summary(climate_BP4k_meantem_ras_proj)
#            bio01
# Min.    3.727445
# 1st Qu. 4.443860
# Median  5.176706
# 3rd Qu. 5.436998
# Max.    5.706034
# NA's    0.000000
res(climate_BP4k_meantem_ras_proj)
# [1] 41200 55500
jpeg(here("pastclim", "plot", "Past_temperature_Beyer_BP4k.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(climate_BP4k_meantem_ras_proj, main = NULL)
dev.off()

# BP3k
climate_BP3k_meantem_raster <- raster(here("pastclim", 'climate_BP3k_meantem.tif'), window=Aohan_owin)
climate_BP3k_meantem_ras_proj <- projectRaster(climate_BP3k_meantem_raster, crs=crs(DEM))
summary(climate_BP3k_meantem_ras_proj)
#            bio01
# Min.    3.985152
# 1st Qu. 4.711642
# Median  5.351734
# 3rd Qu. 5.654961
# Max.    5.925423
# NA's    0.000000
res(climate_BP3k_meantem_ras_proj)
# [1] 41200 55500
jpeg(here("pastclim", "plot", "Past_temperature_Beyer_BP3k.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(climate_BP3k_meantem_ras_proj, main = NULL)
dev.off()

# Krapp2021
# Annual precipitation
# BP4k
climate_BP4k_raster_Krapp <- raster(here("pastclim", 'climate_BP4k_Krapp.tif'), window=Aohan_owin)
climate_BP4k_ras_proj_Krapp <- projectRaster(climate_BP4k_raster_Krapp, crs=crs(DEM))
summary(climate_BP4k_ras_proj_Krapp)
#            bio12
# Min.    401.8155
# 1st Qu. 440.9114
# Median  484.9218
# 3rd Qu. 490.9913
# Max.    565.0544
# NA's      0.0000
res(climate_BP4k_ras_proj_Krapp)
# [1] 41200 55500
jpeg(here("pastclim", "plot", "Past_precipitation_Krapp_BP4k.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(climate_BP4k_ras_proj_Krapp, main = NULL)
dev.off()
# BP3k
climate_BP3k_raster_Krapp <- raster(here("pastclim", 'climate_BP3k_Krapp.tif'), window=Aohan_owin)
climate_BP3k_ras_proj_Krapp <- projectRaster(climate_BP3k_raster_Krapp, crs=crs(DEM))
summary(climate_BP3k_ras_proj_Krapp)
#            bio12
# Min.    399.0624
# 1st Qu. 439.4001
# Median  481.7256
# 3rd Qu. 487.8648
# Max.    558.4910
# NA's      0.0000
res(climate_BP3k_ras_proj_Krapp)
# [1] 41200 55500
jpeg(here("pastclim", "plot", "Past_precipitation_Krapp_BP3k.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(climate_BP3k_ras_proj_Krapp, main = NULL)
dev.off()

# Annual mean temperature
# BP4k
climate_BP4k_meantem_raster_Krapp <- raster(here("pastclim", 'climate_BP4k_meantem_Krapp.tif'), window=Aohan_owin)
climate_BP4k_meantem_ras_proj_Krapp <- projectRaster(climate_BP4k_meantem_raster_Krapp, crs=crs(DEM))
summary(climate_BP4k_meantem_ras_proj_Krapp)
#            bio01
# Min.    4.680763
# 1st Qu. 5.436866
# Median  6.124183
# 3rd Qu. 6.323543
# Max.    6.642521
# NA's    0.000000
res(climate_BP4k_meantem_ras_proj_Krapp)
# [1] 41200 55500
jpeg(here("pastclim", "plot", "Past_temperature_Krapp_BP4k.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(climate_BP4k_meantem_ras_proj_Krapp, main = NULL)
dev.off()


# BP3k
climate_BP3k_meantem_raster_Krapp <- raster(here("pastclim", 'climate_BP3k_meantem_Krapp.tif'), window=Aohan_owin)
climate_BP3k_meantem_ras_proj_Krapp <- projectRaster(climate_BP3k_meantem_raster_Krapp, crs=crs(DEM))
summary(climate_BP3k_meantem_ras_proj_Krapp)
#            bio01
# Min.    4.852989
# 1st Qu. 5.600646
# Median  6.297112
# 3rd Qu. 6.508182
# Max.    6.806907
# NA's    0.000000
res(climate_BP3k_meantem_ras_proj_Krapp)
# [1] 41200 55500
jpeg(here("pastclim", "plot", "Past_temperature_Krapp_BP3k.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(climate_BP3k_meantem_ras_proj_Krapp, main = NULL)
dev.off()


# 3.Comparison with Bioclim data
# read Aohan Bioclim layers
precipitation_raw <- mask(raster(here("GIS_layers","1_10_precipitation_clipped_qgisproj.tif")), Aohan)
temperature_raw <- mask(raster(here("GIS_layers","1_11_mean_tem_clipped_qgisproj.tif")), Aohan)
# 3.1 resolution
res(precipitation_raw)
# [1] 809.456 809.456
res(temperature)
# [1] 809.456 809.456

# 3.2 the values
summary(precipitation_raw)
# Min.                                     377
# 1st Qu.                                  425
# Median                                   440
# 3rd Qu.                                  460
# Max.                                     503
# NA's                                   13360
summary(temperature_raw)
# Min.                           4.625000
# 1st Qu.                        6.762500
# Median                         7.095833
# 3rd Qu.                        7.245833
# Max.                           7.937500
# NA's                       13360.000000

jpeg(here("pastclim", "plot", "Bio_precipitation_raw.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(precipitation_raw, main = NULL)
dev.off()

jpeg(here("pastclim", "plot", "Bio_temperature_raw.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(temperature_raw, main = NULL)
dev.off()





