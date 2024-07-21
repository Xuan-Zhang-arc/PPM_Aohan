# R script
## article title: A comparative analyses of stone- and earth-wall settlement locations of the Lower Xiajiadian Culture in Aohan Banner, China
## journal name: Archaeological and Anthropological Sciences
## author names: Xuan Zhang; Yukun Zhang; Lifeng Tan; Enrico R. Crema; Yanguo Tian; Ze Wang
## corresponding authors: corresponding author at School of Architecture, Tianjin University, China. E-mail address: 961295@tju.edu.cn (Yukun Zhang); corresponding author at McDonald Institute for Archaeological Research, Department of Archaeology,  University of Cambridge, UK. E-mail address: erc62@cam.ac.uk (Enrico R. Crema).


# 1. Prepare the packages and file paths
library(here)
library(sf)
library(spatstat)
library(spatstat.model)
library(raster)
library(terra)
library(sp)
library(MuMIn)
library(MASS)

# Versions of pacakges
packageVersion("here")
# [1] '1.0.1'
packageVersion("sf")
# [1] '1.0.9'
packageVersion("spatstat")
# [1] '3.0.2'
packageVersion("spatstat.model")
# [1] '3.2.1.7'
packageVersion("raster")
# [1] '3.6.3'
packageVersion("terra")
# [1] '1.7.18'
packageVersion("sp")
# [1] '1.5.1'
packageVersion("MuMIn")
# [1] '1.47.5'
packageVersion("MASS")
# [1] '7.3.58.1'






# 2. Read the objects and convert them to spatstat format ----
# The observed window
Aohan <-read_sf(here("GIS_layers","1_1_Aohan_qgisproj.shp"))
Aohan_owin=as.owin(Aohan)
# The CRS
crs_obj <- crs(Aohan)
crs_obj
# [1] "PROJCRS[\"WGS_84_UTM_zone_50N_and_a_bit\",\n    BASEGEOGCRS[\"WGS 84\",\n  
#       DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",63
# 78137,298.257223563,\n                LENGTHUNIT[\"metre\",1]],\n            ID[
# \"EPSG\",6326]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"Degre
# e\",0.0174532925199433]]],\n    CONVERSION[\"unnamed\",\n        METHOD[\"Transv
# erse Mercator\",\n            ID[\"EPSG\",9807]],\n        PARAMETER[\"Latitude 
# of natural origin\",0,\n            ANGLEUNIT[\"Degree\",0.0174532925199433],\n 
#            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\
# ",120,\n            ANGLEUNIT[\"Degree\",0.0174532925199433],\n            ID[\"
# EPSG\",8802]],\n        PARAMETER[\"Scale factor at natural origin\",0.9996,\n  
#           SCALEUNIT[\"unity\",1],\n            ID[\"EPSG\",8805]],\n        PARA
# METER[\"False easting\",500000,\n            LENGTHUNIT[\"metre\",1],\n
#    ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LEN
# GTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n 
#        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre
# \",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"(N)\",north,\n      
#       ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\
# ",9001]]]]"
 

# Stone-wall sites
stone_wall <- read.csv(here("Data", "1_xiaxia_stone_wall.csv"))
stone.sf <- st_as_sf(stone_wall, coords = c("Longitude", "Latitude"), crs = 4326)
stone.sf <- st_transform(stone.sf, crs = crs_obj)
stone.xy <- st_coordinates(stone.sf)
TypeA_ppp <- ppp(x=stone.xy[,1],y=stone.xy[,2], window=Aohan_owin)
# Warning message:
# 5 points were rejected as lying outside the specified window
summary(TypeA_ppp)
# Planar point pattern:  497 points
# Average intensity 6.037471e-08 points per square unit
# 
# Coordinates are given to 1 decimal place
# i.e. rounded to the nearest multiple of 0.1 units
# 
# Window: polygonal boundary
# single connected closed polygon with 2615 vertices
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
#                      (114600 x 148300 units)
# Window area = 8231920000 square units
# Fraction of frame area: 0.485
# 
# *** 5 illegal points stored in attr(,"rejects") ***
rejects <- attr(TypeA_ppp, "rejects")
jpeg(here("Output", "Type_A_ppp_rejects.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(TypeA_ppp$window, main=NULL)
plot(rejects, col='red',add=TRUE)
dev.off()
# From the plot() results, it can be seen that the five rejected points are located on the edge of the observed window, or rather, on the Aohanqi administrative boundary. Since these five points constitute less than 1% of all the 502 stone-wall sites, we used Type_A_ppp for subsequent analysis. 
# Earth-wall sites
earth_wall <- read.csv(here("Data", "1_xiaxia_earth_wall.csv"))
earth.sf <- st_as_sf(earth_wall, coords = c("Longitude", "Latitude"), crs = 4326)
earth.sf <- st_transform(earth.sf, crs = crs_obj)
earth.xy <- st_coordinates(earth.sf)
TypeB_ppp <- ppp(x=earth.xy[,1],y=earth.xy[,2], window=Aohan_owin)
# Stone- and earth-wall sites
AB.multitype<-superimpose(stonewall=TypeA_ppp,earthwall=TypeB_ppp)

# 1.1 elevation
DEM <- mask(raster(here("GIS_layers","1_1_DEM_qgisproj.tif")), Aohan)
res(DEM)
# [1] 26.99211 26.99211
DEM_matrix <- apply(as.matrix(DEM), 2, rev)
elevation_im <- as.im(DEM_matrix,Aohan_owin)
elevation_im$v[is.na(elevation_im$v)] <- 0
elevation_im <- as.im(elevation_im,Aohan_owin)
class(elevation_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_1_elevation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(DEM, main = NULL, col=terrain.colors(12), axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.2 slope
slope <- mask(raster(here("GIS_layers","1_2_slope_qgisproj2.tif")), Aohan)
res(slope)
# [1] 26.99211 26.99211
slope_matrix <- apply(as.matrix(slope), 2, rev)
slope_im<-as.im(slope_matrix, Aohan_owin)
slope_im$v[is.na(slope_im$v)] <- 0
slope_im <- as.im(slope_im,Aohan_owin)
class(slope_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_2_slope.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(slope, main = NULL, col=terrain.colors(6), axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.3 aspect (south and north)
aspect <- mask(raster(here("GIS_layers","1_3_aspect_qgisproj.tif")), Aohan)
res(aspect)
# [1] 26.99211 26.99211
aspect_rad <- aspect * (pi / 180)
north_south <- cos(aspect_rad)*(-1)
north_south_matrix <- apply(as.matrix(north_south), 2, rev)
aspect_im <- as.im(north_south_matrix, Aohan_owin)
aspect_im$v[is.na(aspect_im$v)] <- 0
aspect_im <- as.im(aspect_im,Aohan_owin)
class(aspect_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_3_aspect.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(north_south, main = NULL, col=terrain.colors(4), axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.4 incoming visibility index
visibility_in<-mask(raster(here("GIS_layers","1_4_in_vis_index_qgisproj.tif"),window=Aohan_owin), Aohan)
res(visibility_in)
# [1] 26.99211 26.99211
vis_in_matrix <- apply(as.matrix(visibility_in), 2, rev)
vis_in_im <- as.im(vis_in_matrix, Aohan_owin)
vis_in_im$v[is.na(vis_in_im$v)] <- 0
vis_in_im <- as.im(vis_in_im,Aohan_owin)
class(vis_in_im$v)
# [1] "matrix" "array"
# code for Figure 4
custom_colors <- colorRampPalette(c("white", "green","red"))(7)
jpeg(here("Output", "1_covariates", "1_4_vis_in.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(visibility_in, main = NULL, col=custom_colors, axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.5 outgoing visbility index
visibility_out<-mask(raster(here("GIS_layers","1_5_out_vis_index_qgisproj.tif"),window=Aohan_owin), Aohan)
res(visibility_out)
# [1] 26.99211 26.99211
vis_out_matrix <- apply(as.matrix(visibility_out), 2, rev)
vis_out_im <- as.im(vis_out_matrix, Aohan_owin)
vis_out_im$v[is.na(vis_out_im$v)] <- 0
vis_out_im <- as.im(vis_out_im,Aohan_owin)
class(vis_out_im$v)
# [1] "matrix" "array"
# code for Figure 4
custom_colors <- colorRampPalette(c("grey", "yellow","red"))(7)
jpeg(here("Output", "1_covariates", "1_5_vis_out.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(visibility_out, main = NULL, col=custom_colors, axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.6 landforms
landforms <- mask(raster(here("GIS_layers","1_6_landforms_qgisproj_reclass.tif")), Aohan)
res(landforms)
# [1] 26.99211 26.99211
landforms_matrix <- apply(as.matrix(landforms), 2, rev)
landforms_im<-as.im(landforms_matrix, Aohan_owin)
landforms_im$v[is.na(landforms_im$v)] <- 0
landforms_im <- as.im(landforms_im,Aohan_owin)
landforms_recla_im <- cut(landforms_im$v, 
                          breaks = c(-Inf, 2, 3, Inf), 
                          labels = c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder"), 
                          include.lowest = TRUE)

landforms_recla_im <- im(factor(landforms_recla_im, levels = c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), 
                         xcol = landforms_im$xcol, 
                         yrow = landforms_im$yrow)

class(landforms_recla_im$v)
# [1] "factor"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_6_landforms.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(landforms_recla_im, main = NULL, col = terrain.colors(3), axes = TRUE,
     ribside = "right", 
     ribsep = 0.05, 
     ribargs = list(at = 1:3, 
                    labels = c("flat/pit/\nvalley/footslope", "spur/slope/\nhollow", "peak/ridge/\nshoulder"), 
                    cex.axis = 1, 
                    width = 6, 
                    mgp = c(3, 1.5, 0)))  
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.7 river (distance to rivers)
river <- mask(raster(here("GIS_layers","1_7_river_level_1_to_5_clipped_qgisproj2_distance.tif")), Aohan)
res(river)
# [1] 26.99211 26.99211
river_matrix <- apply(as.matrix(river), 2, rev)
river_im<-as.im(river_matrix, Aohan_owin)
river_im$v[is.na(river_im$v)] <- 0
river_im <- as.im(river_im,Aohan_owin)
class(river_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_7_river.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(river, main = NULL, axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.8 soil (presence and absence of soil possible for growing millet)
soil <- mask(raster(here("GIS_layers","1_8_soil_reclass_qgisproj.tif")), Aohan)
soil_resample <- resample(soil, DEM, method='bilinear')
res(soil_resample)
# [1] 26.99211 26.99211
class(values(soil_resample))
# [1] "numeric"
soil_resample_matrix <- apply(as.matrix(soil_resample), 2, rev)
soil_im<-as.im(soil_resample_matrix, Aohan_owin)
soil_im$v[is.na(soil_im$v)] <- 0
soil_im <- as.im(soil_im,Aohan_owin)
class(soil_im$v)
# [1] "matrix" "array"
soil_categorical <- function(x) {
   as.factor(ifelse(x > 0, "presence", "absence"))
 }
soil_recla_im <- eval.im(soil_categorical(soil_im))
class(soil_recla_im$v)
# [1] "factor"
# code for Figure 4
custom_colors <- colorRampPalette(c("light grey","#D2B48C" ))(2)
jpeg(here("Output", "1_covariates", "1_8_soil.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(soil_recla_im, main = NULL, col=custom_colors, axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.9 geology (distance to source of building materials)
# 1.9.1 distance to materials for constructing stone-wall sites
geo_rock <- mask(raster(here("GIS_layers","1_9_building_rock_qgisproj_distance.tif")), Aohan)
geo_rock_resample <- resample(geo_rock, DEM, method='bilinear')
res(geo_rock_resample)
# [1] 26.99211 26.99211
geo_rock_resample_matrix <- apply(as.matrix(geo_rock_resample), 2, rev)
geo_rock_im<-as.im(geo_rock_resample_matrix, Aohan_owin)
geo_rock_im$v[is.na(geo_rock_im$v)] <- 0
geo_rock_im <- as.im(geo_rock_im,Aohan_owin)
class(geo_rock_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_9_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(geo_rock_resample, main = NULL, axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()
# 1.9.2 distance to materials for constructing earth-wall sites
geo_earth <- mask(raster(here("GIS_layers","1_9_building_earth_qgisproj_distance.tif")), Aohan)
geo_earth_resample <- resample(geo_earth, DEM, method='bilinear')
res(geo_earth_resample)
# [1] 26.99211 26.99211
geo_earth_resample_matrix <- apply(as.matrix(geo_earth_resample), 2, rev)
geo_earth_im<-as.im(geo_earth_resample_matrix, Aohan_owin)
geo_earth_im$v[is.na(geo_earth_im$v)] <- 0
geo_earth_im <- as.im(geo_earth_im,Aohan_owin)
class(geo_earth_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_9_geo_earth.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(geo_earth_resample, main = NULL, axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.10 modern precipitation
precipitation <- mask(raster(here("GIS_layers","1_10_precipitation_clipped_qgisproj.tif")), Aohan)
precipitation_resample <- resample(precipitation, DEM, method='bilinear')
res(precipitation_resample)
# [1] 26.99211 26.99211
precipitation_resample_matrix <- apply(as.matrix(precipitation_resample), 2, rev)
precipitation_im <- as.im(precipitation_resample_matrix, Aohan_owin)
precipitation_im$v[is.na(precipitation_im$v)] <- 0
precipitation_im <- as.im(precipitation_im,Aohan_owin)
class(precipitation_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_10_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(precipitation_resample, main = NULL, col=rev(terrain.colors(10)),axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 1.11 temperature 
temperature <- mask(raster(here("GIS_layers","1_11_mean_tem_clipped_qgisproj.tif")), Aohan)
temperature_resample <- resample(temperature, DEM, method='bilinear')
res(temperature_resample)
# [1] 26.99211 26.99211
temperature_resample_matrix <- apply(as.matrix(temperature_resample), 2, rev)
temperature_im <- as.im(temperature_resample_matrix, Aohan_owin)
temperature_im$v[is.na(temperature_im$v)] <- 0
temperature_im <- as.im(temperature_im,Aohan_owin)
class(temperature_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_11_temperature.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(temperature_resample, main = NULL, col=terrain.colors(10), axes = TRUE)
plot(Aohan_owin, add = TRUE)
dev.off()

# 3. Function form selection: LM or GAM
# 3.1 stone-wall sites
# 1) elevation: GAM  
ppm_A_elevation<-ppm(TypeA_ppp~elevation_im)
AICc(ppm_A_elevation)
# [1] 17158.08
residuals.ppm(ppm_A_elevation)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -8.5181e-09
gam_A_elevation<-ppm(TypeA_ppp~s(elevation_im),use.gam=TRUE)
AICc(gam_A_elevation)
# [1] 17024.77
residuals.ppm(gam_A_elevation)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = 1.8169e-12

# 2) slope: LM 
ppm_A_slope<-ppm(TypeA_ppp~slope_im)
AICc(ppm_A_slope)
# [1] 17316.38
residuals.ppm(ppm_A_slope)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -5.2449e-11
gam_A_slope<-ppm(TypeA_ppp~s(slope_im),use.gam=TRUE)
AICc(gam_A_slope)
# [1] 17296.42
residuals.ppm(gam_A_slope)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -6.5166e-11

# 3) aspect: LM 
ppm_A_aspect<-ppm(TypeA_ppp~aspect_im)
AICc(ppm_A_aspect)
# [1] 17509.44
residuals.ppm(ppm_A_aspect)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = 1.5069e-11
gam_A_aspect<-ppm(TypeA_ppp~s(aspect_im),use.gam=TRUE)
AICc(gam_A_aspect)
# [1] 17525.87
residuals.ppm(gam_A_aspect)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = 1.841e-11

# 4) incoming visibility index: LM 
ppm_A_vis_in<-ppm(TypeA_ppp~vis_in_im)
AICc(ppm_A_vis_in)
# [1] 17432.98
residuals.ppm(ppm_A_vis_in)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = 1.6047e-11
gam_A_vis_in<-ppm(TypeA_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gam_A_vis_in)
# [1] 17369.39
residuals.ppm(gam_A_vis_in)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = 4.0643e-11

# 5) outgoing visibility index: LM
ppm_A_vis_out<-ppm(TypeA_ppp~vis_out_im)
AICc(ppm_A_vis_out)
# [1] 17434.45
residuals.ppm(ppm_A_vis_out)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -2.735e-12
gam_A_vis_out<-ppm(TypeA_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gam_A_vis_out)
# [1] 17385.04
residuals.ppm(gam_A_vis_out)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = 1.2347e-11

# 6) landforms: LM
ppm_A_landforms<-ppm(TypeA_ppp~landforms_recla_im)
AICc(ppm_A_landforms)
# [1] 17482.02

# 7) river:  GAM
ppm_A_river<-ppm(TypeA_ppp~river_im)
AICc(ppm_A_river)
# [1] 17503.9
residuals.ppm(ppm_A_river)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -1.2148e-11

gam_A_river<-ppm(TypeA_ppp~s(river_im),use.gam=TRUE)
AICc(gam_A_river)
# [1] 17495.25
residuals.ppm(gam_A_river)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -3.6381e-12

# 8) soil: LM
ppm_A_soil<-ppm(TypeA_ppp~soil_recla_im)
AICc(ppm_A_soil)
# [1] 17309.31

# 9) geology
# building rock: GAM
ppm_A_geo_rock <- ppm(TypeA_ppp~geo_rock_im)
AICc(ppm_A_geo_rock)
# [1] 17191.7
residuals.ppm(ppm_A_geo_rock)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -2.4607e-06
gam_A_geo_rock <- ppm(TypeA_ppp~s(geo_rock_im), use.gam=TRUE)
AICc(gam_A_geo_rock)
# [1] 17208.13
residuals.ppm(gam_A_geo_rock)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -6.1726e-11

# 10) modern annual precipitation: GAM 
ppm_A_precipitation<-ppm(TypeA_ppp~precipitation_im)
AICc(ppm_A_precipitation)
# [1] 17067.29
residuals.ppm(ppm_A_precipitation)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = 1.3055e-11

gam_A_precipitation<-ppm(TypeA_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gam_A_precipitation)
# [1] 17083.72
residuals.ppm(gam_A_precipitation)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -1.017e-11

# 11) modern annual mean temperature: GAM 
ppm_A_temperature<-ppm(TypeA_ppp~temperature_im)
AICc(ppm_A_temperature)
# [1] 17436.23
residuals.ppm(ppm_A_temperature)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -3.7529e-11

gam_A_temperature<-ppm(TypeA_ppp~s(temperature_im),use.gam=TRUE)
AICc(gam_A_temperature)
# [1] 17040.92
residuals.ppm(gam_A_temperature)
# Scalar-valued measure
# Approximated by 2566 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 497 atoms
# Total mass:
# discrete = 497   continuous = -497   total = -5.3767e-12

# 3.2 earth-wall sites
# 1) elevation: GAM 
ppm_B_elevation<-ppm(TypeB_ppp~elevation_im)
AICc(ppm_B_elevation)
# [1] 2626.254
residuals.ppm(ppm_B_elevation)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -2.8533e-12

gam_B_elevation<-ppm(TypeB_ppp~s(elevation_im),use.gam=TRUE)
AICc(gam_B_elevation)
# [1] 2645.995
residuals.ppm(gam_B_elevation)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 1.296e-12

# 2) slope: GAM 
ppm_B_slope<-ppm(TypeB_ppp~slope_im)
AICc(ppm_B_slope)
# [1] 2633.993
residuals.ppm(ppm_B_slope)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -1.0237e-06

gam_B_slope<-ppm(TypeB_ppp~s(slope_im),use.gam=TRUE)
AICc(gam_B_slope)
# [1] 2653.734
residuals.ppm(gam_B_slope)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 2.3075e-13

# 3) aspect:  LM
ppm_B_aspect<-ppm(TypeB_ppp~aspect_im)
AICc(ppm_B_aspect)
# [1] 2633.21
residuals.ppm(ppm_B_aspect)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -1.1309e-12

gam_B_aspect<-ppm(TypeB_ppp~s(aspect_im),use.gam=TRUE)
AICc(gam_B_aspect)
# [1] 2652.951
residuals.ppm(gam_B_aspect)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 2.2679e-12

# 4) incoming visibility index: GAM
ppm_B_vis_in<-ppm(TypeB_ppp~vis_in_im)
AICc(ppm_B_vis_in)
# [1] 2634.022
residuals.ppm(ppm_B_vis_in)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -9.9415e-07

gam_B_vis_in<-ppm(TypeB_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gam_B_vis_in)
# [1] 2638.175
residuals.ppm(gam_B_vis_in)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -3.2636e-12

# 5) outgoing visibility index:GAM 
ppm_B_vis_out<-ppm(TypeB_ppp~vis_out_im)
AICc(ppm_B_vis_out)
# [1] 2633.61
residuals.ppm(ppm_B_vis_out)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -1.219e-06

gam_B_vis_out<-ppm(TypeB_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gam_B_vis_out)
# [1] 2639.044
residuals.ppm(gam_B_vis_out)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 6.4198e-13

# 6) landforms: LM
ppm_B_landforms<-ppm(TypeB_ppp~landforms_recla_im)
AICc(ppm_B_landforms)
# [1] 2634.567

# 7) river:  GAM
ppm_B_river<-ppm(TypeB_ppp~river_im)
AICc(ppm_B_river)
# [1] 2624.37
residuals.ppm(ppm_B_river)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 5.7066e-13

gam_B_river<-ppm(TypeB_ppp~s(river_im),use.gam=TRUE)
AICc(gam_B_river)
# [1] 2644.111
residuals.ppm(gam_B_river)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 3.6776e-13

# 8) soil: LM
ppm_B_soil<-ppm(TypeB_ppp~soil_recla_im)
AICc(ppm_B_soil)
# [1] 2617.957

# 9) geology
# building earth: LM
ppm_B_geo_earth <- ppm(TypeB_ppp~geo_earth_im)
AICc(ppm_B_geo_earth)
# [1] 2632.468
residuals.ppm(ppm_B_geo_earth)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 1.5052e-13
gam_B_geo_earth <- ppm(TypeB_ppp~s(geo_earth_im), use.gam=TRUE)
AICc(gam_B_geo_earth)
# [1] 2652.209
residuals.ppm(gam_B_geo_earth)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 2.1673e-12

# 10) modern annual precipitation:  GAM
ppm_B_precipitation<-ppm(TypeB_ppp~precipitation_im)
AICc(ppm_B_precipitation)
# [1] 2581.603
residuals.ppm(ppm_B_precipitation)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -2.2591e-12

gam_B_precipitation<-ppm(TypeB_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gam_B_precipitation)
# [1] 2601.344
residuals.ppm(gam_B_precipitation)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = 4.743e-13

# 11) modern annual mean temperature: GAM 
ppm_B_temperature<-ppm(TypeB_ppp~temperature_im)
AICc(ppm_B_temperature)
# [1] 2632.182
residuals.ppm(ppm_B_temperature)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -1.2032e-06

gam_B_temperature<-ppm(TypeB_ppp~s(temperature_im, k=3),use.gam=TRUE)
AICc(gam_B_temperature)
# [1] 2634.375
residuals.ppm(gam_B_temperature)
# Scalar-valued measure
# Approximated by 620 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 67 atoms
# Total mass:
# discrete = 67   continuous = -67   total = -1.6247e-12

# 4. Relative Variable Importance (This section may take 10 or more hours.)   
# 4.1 stone-wall sites
global.model.A<-ppm(TypeA_ppp~s(elevation_im, k=5)+slope_im+aspect_im+vis_in_im+vis_out_im+landforms_recla_im+s(river_im, k=5)+soil_recla_im+s(geo_rock_im, k=5)+s(precipitation_im, k=5)+s(temperature_im, k=5), use.gam=TRUE) # The k value was defined to ensure that all models were successfully fitted. This k value was the minimum value selected after trying different values.

all.model.A <- dredge(global.model.A,evaluate=FALSE)
all.model2.A <- lapply(all.model.A,function(x){as.character(x)[4]})
candidates.aic.A  <- numeric(length=length(all.model2.A))
for (i in 1:length(candidates.aic.A))
{
	candidates.aic.A[i] <- AICc(ppm(Q=as.formula(paste0('TypeA_ppp',all.model2.A[[i]])),use.gam=TRUE))
}

params.A <- getAllTerms(global.model.A) |> as.character()
permuts.A <- vector('list', length=length(params.A))
for (i in 1:length(permuts.A)) {
  permuts.A[[i]] <- c(0, 1)
}
param.comb.A <- expand.grid(permuts.A)
colnames(param.comb.A) <- params.A
raw_weights <- Weights(candidates.aic.A)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.A$weights <- rounded_weights / sum(rounded_weights)


param.importance.A <- data.frame(params.A=params.A, aic.imp=NA)

for (i in 1:length(params.A)) {
  param.importance.A$aic.imp[i] <- sum(param.comb.A[, i] * param.comb.A$weights)
}

param.importance.A 
#                      params.A   aic.imp
# 1                   aspect_im 0.9666867
# 2          landforms_recla_im 0.9361745
# 3      s(elevation_im, k = 5) 1.0000000
# 4       s(geo_rock_im, k = 5) 0.9895958
# 5  s(precipitation_im, k = 5) 1.0000000
# 6          s(river_im, k = 5) 0.5671269
# 7    s(temperature_im, k = 5) 0.3890556
# 8                    slope_im 0.9933974
# 9               soil_recla_im 0.3221289
# 10                  vis_in_im 0.8163265
# 11                 vis_out_im 0.3695478



sum(param.importance.A$aic.imp)
sum(param.comb.A$weights)
head(param.comb.A)
elevation.results=subset(param.comb.A, param.comb.A[,3]==1)
sum(elevation.results$weights) # 1
nrow(elevation.results) # [1] 1024
slope.results=subset(param.comb.A, param.comb.A[,8]==1)
sum(slope.results$weights) # [1] 0.9933974
nrow(elevation.results) # [1] 1024
river.results=subset(param.comb.A, param.comb.A[,6]==1)
sum(river.results$weights) # [1] 0.5671269
nrow(river.results) # [1] 1024

# 4.2 earth-wall sites
global.model.B<-ppm(TypeB_ppp~s(elevation_im, k=3)+s(slope_im, k=3)+aspect_im+s(vis_in_im, k=3)+s(vis_out_im, k=3)+landforms_recla_im+s(river_im, k=3)+soil_recla_im+geo_earth_im+s(precipitation_im, k=3)+s(temperature_im, k=3), use.gam=TRUE) # The k value was defined to ensure that all models were successfully fitted. This k value was the minimum value selected after trying different values. 

all.model.B <- dredge(global.model.B,evaluate=FALSE)
all.model2.B <- lapply(all.model.B,function(x){as.character(x)[4]})

candidates.aic.B  <- numeric(length=length(all.model2.B))
for (i in 1:length(candidates.aic.B))
{
	candidates.aic.B[i] <- AICc(ppm(Q=as.formula(paste0('TypeB_ppp',all.model2.B[[i]])),use.gam=TRUE))
}

params.B <- getAllTerms(global.model.B) |> as.character()
permuts.B <- vector('list', length=length(params.B))
for (i in 1:length(permuts.B)) {
  permuts.B[[i]] <- c(0, 1)
}

param.comb.B <- expand.grid(permuts.B)
colnames(param.comb.B) <- params.B
raw_weights <- Weights(candidates.aic.B)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.B$weights <- rounded_weights / sum(rounded_weights)

param.importance.B <- data.frame(params.B=params.B, aic.imp=NA)

for (i in 1:length(params.B)) {
  param.importance.B$aic.imp[i] <- sum(param.comb.B[, i] * param.comb.B$weights)
}

param.importance.B
#                      params.B    aic.imp
# 1                   aspect_im 0.30460520
# 2                geo_earth_im 0.84067422
# 3          landforms_recla_im 0.39791311
# 4      s(elevation_im, k = 3) 0.12180195
# 5  s(precipitation_im, k = 3) 1.00000000
# 6          s(river_im, k = 3) 0.07524832
# 7          s(slope_im, k = 3) 0.21400622
# 8    s(temperature_im, k = 3) 0.19042841
# 9         s(vis_in_im, k = 3) 0.12732016
# 10       s(vis_out_im, k = 3) 0.09270593
# 11              soil_recla_im 0.22785191

# 4.3 comparison between stone-wall sites and earth-wall sites
# code for Figure 5
importance.A.df <- data.frame(covariates=c("elevation", "slope", "aspect", "visibility_in", "visibility_out", "landforms", "river", "soil", "geology", "precipitation", "temperature"),importance.A=c(1.0000000, 0.9933974, 0.9666867, 0.8163265, 0.3695478, 0.9361745, 0.5671269, 0.3221289, 0.9895958, 1.0000000, 0.3890556))
importance.B.df <- data.frame(covariates=c("elevation", "slope", "aspect", "visibility_in", "visibility_out", "landforms", "river", "soil", "geology", "precipitation", "temperature"),importance.B=c(0.12180195, 0.21400622, 0.30460520, 0.12732016, 0.09270593, 0.39791311, 0.07524832, 0.22785191, 0.84067422, 1.00000000, 0.19042841))
importance.A.B.df <- data.frame(covariates=factor(importance.A.df$covariates), importance.A=as.numeric(importance.A.df$importance.A), importance.B=as.numeric(importance.B.df$importance.B))
jpeg(here("Output", "2_first_order_interactions", "AICc_importance.jpg"),width = 12, height = 7, units = "in", res = 1200)
par(mar = c(5, 10, 4, 3) + 0.1) 
barplot(
	beside = TRUE,
	horiz = TRUE,
	height = rbind(importance.A.B.df$importance.A, importance.A.B.df$importance.B), 
	names.arg = importance.A.B.df$covariates,
	col = c("white", "gray"),
	xlab = "AICc importance",
	ylab = mtext("Covariates", side = 2, line = 8, las = 0),
	main = "AICc importance comparison",
	legend.text = TRUE,las = 1, 
	xlim = c(0, 1)
)
abline(v = 0.8, col = "red", lty = 2)
legend("bottomright", legend = c("Stone-wall sites", "earth-wall sites"), fill = c("white", "gray"))
dev.off()


# 5. Final models and goodness-of-fit
# 5.1 stone-wall sites
model.A<-ppm(TypeA_ppp~aspect_im+landforms_recla_im+s(elevation_im, k=5)+s(geo_rock_im, k=5)+s(precipitation_im, k=5)+slope_im+vis_in_im,use.gam=TRUE) # Keep the k value in the global model.
pseudoR2(model.A)
# [1] 0.2723642
model.A.Kres.env <-envelope(model.A,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "2_first_order_interactions", "model_A_Kres.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(model.A.Kres.env, main = NULL)
dev.off()

model.A2<-ppm(TypeA_ppp~aspect_im+landforms_recla_im+s(elevation_im)+s(geo_rock_im)+s(precipitation_im)+slope_im+vis_in_im,use.gam=TRUE) # Remove or reduce the k value to get a better fitted final model.
pseudoR2(model.A2)
# [1] 0.3077873
# code for Figure 6
model.A2.Kres.env <-envelope(model.A2,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "2_first_order_interactions", "model_A2_Kres.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(model.A2.Kres.env, main = NULL)
dev.off()

global.model.A.Kres.env <-envelope(global.model.A,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "2_first_order_interactions", "global_model_A_Kres.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(global.model.A.Kres.env,main=NULL)
dev.off()
# code for Figure 7
jpeg(here("Output", "2_first_order_interactions", "effectfun_A2_aspect.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(model.A2,covname='aspect_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), elevation_im=mean(elevation_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A2_landforms.jpg"),width = 7, height = 7, units = "in", res = 1200)
eff <- effectfun(model.A2,covname='landforms_recla_im', aspect_im=mean(aspect_im), elevation_im=mean(elevation_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE)
barplot(eff$lambda,names.arg=eff$landforms_recla_im)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A2_elevation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(model.A2,covname='elevation_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A2_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(model.A2,covname='geo_rock_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A2_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(model.A2,covname='precipitation_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A2_slope.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(model.A2,covname='slope_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), precipitation_im=mean(precipitation_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A2_vis_in.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(model.A2,covname='vis_in_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), precipitation_im=mean(precipitation_im), slope_im=mean(slope_im),se.fit=TRUE),main=NULL)
dev.off()

# 5.2 earth-wall sites
model.B<-ppm(TypeB_ppp~geo_earth_im+s(precipitation_im, k = 3),use.gam=TRUE) # Keep the k value in the global model.
pseudoR2(model.B)
# [1] 0.1554673
model.B.Kres.env <-envelope(model.B,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "2_first_order_interactions", "model_B_Kres.jpg"),width = 7, height = 7, units = "in",res = 1200)
plot(model.B.Kres.env, main = NULL)
dev.off()

model.B2<-ppm(TypeB_ppp~geo_earth_im+s(precipitation_im),use.gam=TRUE) # Remove or reduce the k value to get a better fitted final model. 
pseudoR2(model.B2)
# [1] 0.1554673
# code for Figure 6
model.B2.Kres.env <-envelope(model.B2,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "2_first_order_interactions", "model_B2_Kres.jpg"),width = 7, height = 7, units = "in",res = 1200)
plot(model.B2.Kres.env, main = NULL)
dev.off()

global.model.B.Kres.env <-envelope(global.model.B,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "2_first_order_interactions", "global_model_B_Kres.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(global.model.B.Kres.env,main=NULL)
dev.off()

# code for Figure 7
jpeg(here("Output", "2_first_order_interactions", "effectfun_B2_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(model.B2,covname='precipitation_im',geo_earth_im=mean(geo_earth_im), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_B2_geo_earth.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(model.B2,covname='geo_earth_im',precipitation_im=mean(precipitation_im), se.fit=TRUE),main=NULL)
dev.off()

# 6. second-order interaction
# 6.1 stone-wall sites
# 6.1.1 Strauss model
rvals <- data.frame(r=seq(100, 15000, by=100))
distance_A <-profilepl(rvals, Strauss, TypeA_ppp~aspect_im+landforms_recla_im+s(elevation_im, k=3)+s(geo_rock_im, k=3)+s(precipitation_im, k=3)+slope_im+vis_in_im, use.gam=TRUE) # The k value was defined to ensure that all models were successfully fitted. This k value was the minimum value selected after trying different values. 
print(distance_A)
# profile log pseudolikelihood for model:
# ppm(TypeA_ppp ~ aspect_im + landforms_recla_im + s(elevation_im,       k = 3) + 
# s(geo_rock_im,  k = 3) + s(precipitation_im,  k = 3) +      slope_im + 
# vis_in_im,  use.gam = TRUE,  interaction = Strauss)
# fitted with rbord = 15000
# interaction: Strauss process
# irregular parameter: r in [100, 15000]
# optimum value of irregular parameter:  r = 1600
jpeg(here("Output", "4_second_order_interactions", "Gibbs_A_distance.jpg"), width = 7, height = 7, units = "in", res = 1200)
plot(distance_A, main="")
text(x = 1600, y = par("usr")[3], labels = 1600, pos = 3, col = "black")
dev.off()
Gibbs_A <- as.ppm(distance_A)
pseudoR2(Gibbs_A)
# [1] 0.2285206
# kres
Gibbs_A.Kres.env <-envelope(Gibbs_A,fun=Kres,nsim=50,correction='best')
# Error: This calculation is not supported for GAM fits
Gibbs_A.Kres <- Kres(Gibbs_A,correction='best')
jpeg(here("Output", "4_second_order_interactions", "Gibbs_A_Kres_noenv.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(Gibbs_A.Kres, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_aspect.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A,covname='aspect_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), elevation_im=mean(elevation_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_landforms.jpg"),width = 7, height = 7, units = "in", res = 1200)
eff <- effectfun(Gibbs_A,covname='landforms_recla_im', aspect_im=mean(aspect_im), elevation_im=mean(elevation_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE)
barplot(eff$lambda,names.arg=eff$landforms_recla_im)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_elevation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A,covname='elevation_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A,covname='geo_rock_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A,covname='precipitation_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_slope.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A,covname='slope_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), precipitation_im=mean(precipitation_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_vis_in.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A,covname='vis_in_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), precipitation_im=mean(precipitation_im), slope_im=mean(slope_im),se.fit=TRUE),main=NULL)
dev.off()



# 6.1.2 area interaction model
rvals <- data.frame(r=seq(100, 15000, by=100))
distance_A2 <-profilepl(rvals, AreaInter, TypeA_ppp~aspect_im+landforms_recla_im+s(elevation_im, k=3)+s(geo_rock_im, k=3)+s(precipitation_im, k=3)+slope_im+vis_in_im, use.gam=TRUE) # The k value was defined to ensure that all models were successfully fitted. This k value was the minimum value selected after trying different values. 
# (computing rbord)
# comparing 150 models...
# 1, Creating a 33.5-megapixel window mask
# Creating a 33.5-megapixel window mask
# Creating a 35.6-megapixel window mask
# Creating a 35.6-megapixel window mask
# Error: cannot allocate vector of size 2.2 Gb
# In addition: Warning messages:
# 1: Creating a 33.5-megapixel window mask
# 2: Creating a 33.5-megapixel window mask
# 3: Creating a 35.6-megapixel window mask
# 4: Creating a 35.6-megapixel window mask

rvals <- data.frame(r=seq(100, 15000, by=200))
distance_A2 <-profilepl(rvals, AreaInter, TypeA_ppp~aspect_im+landforms_recla_im+s(elevation_im)+s(geo_rock_im)+s(precipitation_im)+slope_im+vis_in_im, use.gam=TRUE)

print(distance_A2)

jpeg(here("Output", "4_second_order_interactions", "Gibbs_A2_distance.jpg"), width = 7, height = 7, units = "in", res = 1200)
plot(distance_A2, main=NULL)
dev.off()
Gibbs_A2<- as.ppm(distance_A2)
pseudoR2(Gibbs_A2)
#
# kres
Gibbs_AGibbs_A2.Kres.env <-envelope(Gibbs_A2,fun=Kres,nsim=50,correction='best')
# 
Gibbs_A2.Kres <- Kres(Gibbs_A2,correction='best')
jpeg(here("Output", "4_second_order_interactions", "Gibbs_A2_Kres_noenv.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(Gibbs_A2.Kres, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A2_aspect.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A2,covname='aspect_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), elevation_im=mean(elevation_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A2_landforms.jpg"),width = 7, height = 7, units = "in", res = 1200)
eff <- effectfun(Gibbs_A2,covname='landforms_recla_im', aspect_im=mean(aspect_im), elevation_im=mean(elevation_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE)
barplot(eff$lambda,names.arg=eff$landforms_recla_im)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A2_elevation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A2,covname='elevation_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A2_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A2,covname='geo_rock_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A2_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A2,covname='precipitation_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), slope_im =mean(slope_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A2_slope.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A2,covname='slope_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), precipitation_im=mean(precipitation_im), vis_in_im=mean(vis_in_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A2_vis_in.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(effectfun(Gibbs_A2,covname='vis_in_im',landforms_recla_im =as.factor(c("flat/pit/valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), aspect_im=mean(aspect_im),elevation_im=mean(elevation_im),geo_rock_im = mean(geo_rock_im), precipitation_im=mean(precipitation_im), slope_im=mean(slope_im),se.fit=TRUE),main=NULL)
dev.off()





# 6.2 interaction between types
AB_kcross <- envelope(AB.multitype,fun=Kcross,nsim=50,correction='best')
jpeg(here("Output", "4_second_order_interactions", "AB_kcross.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(AB_kcross,main=NULL)
dev.off()


# 7. site size and model
# 7.1 divide the sites based on size
# Stone-wall sites
stone_wall_filtered <- subset(stone_wall, Area != 0)
stone_wall_sorted <- stone_wall_filtered[order(-stone_wall_filtered$Area), ]
n <- nrow(stone_wall_sorted)
n1 <- ceiling(n / 3)  
n2 <- ceiling((n - n1) / 2)  
A_large <- stone_wall_sorted[1:n1, ]
A_middle <- stone_wall_sorted[(n1 + 1):(n1 + n2), ]
A_small <- stone_wall_sorted[(n1 + n2 + 1):n, ]

A_large.sf <- st_as_sf(A_large, coords = c("Longitude", "Latitude"), crs = 4326)
A_large.sf <- st_transform(A_large.sf, crs = crs_obj)
A_large.xy <- st_coordinates(A_large.sf)
A_large_ppp <- ppp(x=A_large.xy[,1],y=A_large.xy[,2], window=Aohan_owin)
plot(A_large_ppp)

A_middle.sf <- st_as_sf(A_middle, coords = c("Longitude", "Latitude"), crs = 4326)
A_middle.sf <- st_transform(A_middle.sf, crs = crs_obj)
A_middle.xy <- st_coordinates(A_middle.sf)
A_middle_ppp <- ppp(x=A_middle.xy[,1],y=A_middle.xy[,2], window=Aohan_owin)
plot(A_middle_ppp)

A_small.sf <- st_as_sf(A_small, coords = c("Longitude", "Latitude"), crs = 4326)
A_small.sf <- st_transform(A_small.sf, crs = crs_obj)
A_small.xy <- st_coordinates(A_small.sf)
A_small_ppp <- ppp(x=A_small.xy[,1],y=A_small.xy[,2], window=Aohan_owin)
plot(A_small_ppp)

# Earth-wall sites
earth_wall_filtered <- subset(earth_wall, Area != 0)
earth_wall_sorted <- earth_wall_filtered[order(-earth_wall_filtered$Area), ]
n <- nrow(earth_wall_sorted)
n1 <- ceiling(n / 3)  
n2 <- ceiling((n - n1) / 2)  
B_large <- earth_wall_sorted[1:n1, ]
B_middle <- earth_wall_sorted[(n1 + 1):(n1 + n2), ]
B_small <- earth_wall_sorted[(n1 + n2 + 1):n, ]

B_large.sf <- st_as_sf(B_large, coords = c("Longitude", "Latitude"), crs = 4326)
B_large.sf <- st_transform(B_large.sf, crs = crs_obj)
B_large.xy <- st_coordinates(B_large.sf)
B_large_ppp <- ppp(x=B_large.xy[,1],y=B_large.xy[,2], window=Aohan_owin)
plot(B_large_ppp)

B_middle.sf <- st_as_sf(B_middle, coords = c("Longitude", "Latitude"), crs = 4326)
B_middle.sf <- st_transform(B_middle.sf, crs = crs_obj)
B_middle.xy <- st_coordinates(B_middle.sf)
B_middle_ppp <- ppp(x=B_middle.xy[,1],y=B_middle.xy[,2], window=Aohan_owin)
plot(B_middle_ppp)

B_small.sf <- st_as_sf(B_small, coords = c("Longitude", "Latitude"), crs = 4326)
B_small.sf <- st_transform(B_small.sf, crs = crs_obj)
B_small.xy <- st_coordinates(B_small.sf)
B_small_ppp <- ppp(x=B_small.xy[,1],y=B_small.xy[,2], window=Aohan_owin)
plot(B_small_ppp)

# 7.2 test the model performance using the selected covariates above through relative variable improtance analysis as suggested by reviewer 1
# 7.2.1 stone-wall sites
# 1) large
model.A_large<-ppm(A_large_ppp~aspect_im+landforms_recla_im+s(elevation_im)+s(geo_rock_im)+s(precipitation_im)+slope_im+vis_in_im,use.gam=TRUE)
pseudoR2(model.A_large)
# [1] 0.2940389
# kres
model.A_large.Kres.env <-envelope(model.A_large,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A_large_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A_large.Kres.env, main = NULL)
dev.off()

# 2) middle
model.A_middle<-ppm(A_middle_ppp~aspect_im+landforms_recla_im+s(elevation_im, k=3)+s(geo_rock_im, k=3)+s(precipitation_im, k=3)+slope_im+vis_in_im,use.gam=TRUE)
pseudoR2(model.A_middle)
# [1] 0.3933742
# kres
model.A_middle.Kres.env <-envelope(model.A_middle,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A_middle_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A_middle.Kres.env, main = NULL)
dev.off()

# 3) small
model.A_small<-ppm(A_small_ppp~aspect_im+landforms_recla_im+s(elevation_im)+s(geo_rock_im)+s(precipitation_im)+slope_im+vis_in_im,use.gam=TRUE)
pseudoR2(model.A_small)
# [1] 0.3260811
# kres
model.A_small.Kres.env <-envelope(model.A_small,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A_small_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A_small.Kres.env, main = NULL)
dev.off()


# 7.2.2 earth-wall sites
# 1) large
model.B_large<-ppm(B_large_ppp~geo_earth_im+s(precipitation_im),use.gam=TRUE)
pseudoR2(model.B_large)
# [1] 0.07658564
# kres
model.B_large.Kres.env <-envelope(model.B_large,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "model_B_large_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B_large.Kres.env, main = NULL)
dev.off()



# 2) middle
model.B_middle<-ppm(B_middle_ppp~geo_earth_im+s(precipitation_im),use.gam=TRUE)
pseudoR2(model.B_middle)
# [1] 0.137227
# kres
model.B_middle.Kres.env <-envelope(model.B_middle,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "model_B_middle_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B_middle.Kres.env, main = NULL)
dev.off()



# 3) small
model.B_small<-ppm(B_small_ppp~geo_earth_im+s(precipitation_im),use.gam=TRUE)
pseudoR2(model.B_small)
# [1] 0.1293161
# kres
model.B_small.Kres.env <-envelope(model.B_small,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "model_B_small_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B_small.Kres.env, main = NULL)
dev.off()

























