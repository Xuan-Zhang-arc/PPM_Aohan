# R script
## article title: A comparative analyses of stone- and earth-wall settlement locations of the Lower Xiajiadian Culture in Aohan Banner, China
## journal name: Journal of Archaeological Method and Theory
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

# 2. Read the objects and convert them to spatstat format ----
# The observed window
Aohan <-read_sf(here("GIS_layers","1_1_Aohan_qgisproj.shp"))
Aohan_owin=as.owin(Aohan)
# The CRS
crs_obj <- crs(Aohan)
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
plot(elevation_im, main = NULL, col=terrain.colors(10))
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
plot(slope_im, main = NULL, col=terrain.colors(8))
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
plot(aspect_im, main = NULL, col=terrain.colors(4))
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
custom_colors <- colorRampPalette(c("grey", "green","yellow", "red"))(10)
jpeg(here("Output", "1_covariates", "1_4_vis_in.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(vis_in_im, main = NULL, col=custom_colors)
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
custom_colors <- colorRampPalette(c("grey", "yellow","orange", "red"))(10)
jpeg(here("Output", "1_covariates", "1_5_vis_out.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(vis_out_im, main = NULL, col=custom_colors)
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
custom_colors <- colorRampPalette(c("light green","light yellow", "red"))(3)
jpeg(here("Output", "1_covariates", "1_6_landforms.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(landforms_recla_im, main = NULL, col = custom_colors, 
     ribside = "right", 
     ribsep = 0.05, 
     ribargs = list(at = 1:3, 
                    labels = c("flat/pit/\nvalley/footslope", "spur/slope/\nhollow", "peak/ridge/\nshoulder"), 
                    cex.axis = 1, 
                    width = 6, 
                    mgp = c(3, 1.5, 0)))  
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
plot(river_im, main = NULL)
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
custom_colors <- colorRampPalette(c("grey", "light green"))(2)
jpeg(here("Output", "1_covariates", "1_8_soil.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(soil_recla_im, main = NULL, col=custom_colors)
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
plot(geo_rock_im, main = NULL)
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
plot(geo_earth_im, main = NULL)
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
custom_colors <- colorRampPalette(c("white","yellow", "light blue"))(10)
jpeg(here("Output", "1_covariates", "1_10_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(precipitation_im, main = NULL, col=custom_colors)
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
custom_colors <- colorRampPalette(c("white","light green","light yellow","orange"))(16)
jpeg(here("Output", "1_covariates", "1_11_temperature.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(temperature_im, main = NULL, col=custom_colors)
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
plot(distance_A, main=NULL)
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

# 6.2 interaction between types
AB_kcross <- envelope(AB.multitype,fun=Kcross,nsim=50,correction='best')
jpeg(here("Output", "4_second_order_interactions", "AB_kcross.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(AB_kcross,main=NULL)
dev.off()
# 6.2.1 function form decision 
# 1) elevation: GAM  
ppmAB_elevation<-ppm(AB.multitype~elevation_im)
residuals.ppm(ppmAB_elevation)
# Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -5.6731e-11
gamAB_elevation<-ppm(AB.multitype~s(elevation_im),use.gam=TRUE)
residuals.ppm(gamAB_elevation)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = 1.9029e-11

# 2) slope: GAM    
ppmAB_slope<-ppm(AB.multitype~slope_im)
residuals.ppm(ppmAB_slope)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -5.2827e-07
gamAB_slope<-ppm(AB.multitype~s(slope_im),use.gam=TRUE)
residuals.ppm(gamAB_slope)
# Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -5.9508e-11

# 3) aspect: GAM   
ppmAB_aspect<-ppm(AB.multitype~aspect_im)
residuals.ppm(ppmAB_aspect)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -4.6137e-08
gamAB_aspect<-ppm(AB.multitype~s(aspect_im),use.gam=TRUE)
residuals.ppm(gamAB_aspect)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = 3.2744e-11

# 4) incoming visibility index: GAM    
ppmAB_vis_in<-ppm(AB.multitype~vis_in_im)
residuals.ppm(ppmAB_vis_in)
# Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -5.5297e-08
gamAB_vis_in<-ppm(AB.multitype~s(vis_in_im),use.gam=TRUE)
residuals.ppm(gamAB_vis_in)
# Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = 2.4083e-12

# 5) outgoing visibility index: GAM  
ppmAB_vis_out<-ppm(AB.multitype~vis_out_im)
residuals.ppm(ppmAB_vis_out)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall 
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -4.4554e-08
gamAB_vis_out<-ppm(AB.multitype~s(vis_out_im),use.gam=TRUE)
residuals.ppm(gamAB_vis_out)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -7.412e-11

# 6) landforms:LM  
ppmAB_landforms<-ppm(AB.multitype~landforms_recla_im)
residuals.ppm(ppmAB_landforms)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall 
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -5.8304e-07

# 7) river: GAM    
ppmAB_river<-ppm(AB.multitype~river_im)
residuals.ppm(ppmAB_river)
#Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall 
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -3.4299e-08 
gamAB_river<-ppm(AB.multitype~s(river_im),use.gam=TRUE)
residuals.ppm(gamAB_river)
#Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -8.2234e-12 

# 8) soil: LM   
ppmAB_soil<-ppm(AB.multitype~soil_recla_im)
residuals.ppm(ppmAB_soil)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall 
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -1.3662e-07

# 9) geology
# building earth:  GAM   
ppmAB_geo_rock <- ppm(AB.multitype~geo_rock_im)
residuals.ppm(ppmAB_geo_rock)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -3.1687e-09
gamAB_geo_rock <- ppm(AB.multitype~s(geo_rock_im), use.gam=TRUE)
residuals.ppm(gamAB_geo_rock)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = 2.6622e-10

# building earth:  GAM   
ppmAB_geo_earth <- ppm(AB.multitype~geo_earth_im)
residuals.ppm(ppmAB_geo_earth)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -2.6849e-07
gamAB_geo_earth <- ppm(AB.multitype~s(geo_earth_im), use.gam=TRUE)
residuals.ppm(gamAB_geo_earth)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = 8.7619e-11

# 10) modern annual precipitation:  GAM
ppmAB_precipitation<-ppm(AB.multitype~precipitation_im)
residuals.ppm(ppmAB_precipitation)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -1.4609e-11
gamAB_precipitation<-ppm(AB.multitype~s(precipitation_im),use.gam=TRUE)
residuals.ppm(gamAB_precipitation)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -1.0785e-11

# 11) modern annual mean temperature: GAM    
ppmAB_temperature<-ppm(AB.multitype~temperature_im)
residuals.ppm(ppmAB_temperature)
#  Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = -2.5045e-08
gamAB_temperature<-ppm(AB.multitype~s(temperature_im),use.gam=TRUE)
residuals.ppm(gamAB_temperature)
# Scalar-valued measure
#         Defined on 2-dimensional space x marks
#         Possible marks:  stonewall, earthwall
# Approximated by 5266 quadrature points
# window: polygonal boundary
# enclosing rectangle: [458507.4, 573085.7] x [4615669, 4763928] units
# 564 atoms
# Total mass:
# discrete = 564   continuous = -564   total = 2.8442e-11

# 6.2.2 relative variable importance
global.model.AB<-ppm(AB.multitype~s(elevation_im)+s(slope_im)+aspect_im+s(vis_in_im)+s(vis_out_im)+landforms_recla_im+s(river_im)+soil_recla_im+geo_rock_im+geo_earth_im+precipitation_im+s(temperature_im), use.gam=TRUE)
all.model.AB <- dredge(global.model.AB,evaluate=FALSE)
all.model2.AB <- lapply(all.model.AB,function(x){as.character(x)[4]})
candidates.aic.AB  <- numeric(length=length(all.model2.AB))
for (i in 1:length(candidates.aic.AB))
{
	candidates.aic.AB[i] <- AICc(ppm(Q=as.formula(paste0('AB.multitype',all.model2.AB[[i]])),use.gam=TRUE))
}

params.AB <- getAllTerms(global.model.AB) |> as.character()
permuts.AB <- vector('list', length=length(params.AB))
for (i in 1:length(permuts.AB)) {
  permuts.AB[[i]] <- c(0, 1)
}
param.comb.AB <- expand.grid(permuts.AB)
colnames(param.comb.AB) <- params.AB
raw_weights <- Weights(candidates.aic.AB)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.AB$weights <- rounded_weights / sum(rounded_weights)

param.importance.AB <- data.frame(params.AB=params.AB, aic.imp=NA)

for (i in 1:length(params.AB)) {
  param.importance.AB$aic.imp[i] <- sum(param.comb.AB[, i] * param.comb.AB$weights)
}
param.importance.AB
#             params.AB    aic.imp
# 1           aspect_im 0.26527958
# 2        geo_earth_im 0.29368811
# 3         geo_rock_im 0.60998299
# 4  landforms_recla_im 0.04031209
# 5    precipitation_im 1.00000000
# 6     s(elevation_im) 1.00000000
# 7         s(river_im) 0.99989997
# 8         s(slope_im) 0.60058017
# 9   s(temperature_im) 0.00020006
# 10       s(vis_in_im) 0.41202361
# 11      s(vis_out_im) 0.97799340
# 12      soil_recla_im 1.00000000

# 6.2.3 models and goodness of fit
model.AB <- ppm(AB.multitype~precipitation_im+s(elevation_im)+s(river_im)+s(vis_out_im)+soil_recla_im, use.gam=TRUE) 

model.AB.Kres.env <-envelope(model.AB,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "4_second_order_interactions", "model_AB_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.AB.Kres.env, main = NULL)
dev.off()

dist.AB <- pairdist(AB.multitype)
mean(dist.AB) # [1] 34812.71
max(dist.AB) # [1] 107376
AB.RR <- data.frame(R=seq(200,34812.71, by=200))
AB_MS <- function(R) { MultiStrauss(radii=diag(c(R,R))) }
AB.pm <- profilepl(AB.RR, AB_MS,AB.multitype~marks*(precipitation_im+s(elevation_im, k=3)+s(river_im, k=3)+s(vis_out_im, k=3)+soil_recla_im), use.gam=TRUE,rbord=NULL,correction="Ripley")
#comparing 174 models...
# 1, Error in model.frame.default(formula = .mpl.Y ~ marks + precipitation_im + 
#  :
#   invalid type (list) for variable 's(elevation_im, k = 3)'
plot(AB.pm)
as.ppm(AB.pm)


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
model.A_large<-ppm(A_large_ppp~geo_rock_im+s(elevation_im) +s(precipitation_im) +s(river_im) +s(vis_out_im),use.gam=TRUE)
pseudoR2(model.A_large)
# [1] 0.3161105
# kres
model.A_large.Kres.env <-envelope(model.A_large,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A_large_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A_large.Kres.env, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_A_large_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_large,covname='geo_rock_im',elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_large_elevation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_large,covname='elevation_im',geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_large_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_large,covname='precipitation_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()



# 2) middle
model.A_middle<-ppm(A_middle_ppp~geo_rock_im+s(elevation_im) +s(precipitation_im) +s(river_im) +s(vis_out_im),use.gam=TRUE)
pseudoR2(model.A_middle)
# [1] 0.4836151
# kres
model.A_middle.Kres.env <-envelope(model.A_middle,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A_middle_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A_middle.Kres.env, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_A_middle_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_middle,covname='geo_rock_im',elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_middle_elevation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_middle,covname='elevation_im',geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_middle_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_middle,covname='precipitation_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()


# 3) small
model.A_small<-ppm(A_small_ppp~geo_rock_im+s(elevation_im) +s(precipitation_im) +s(river_im) +s(vis_out_im),use.gam=TRUE)
pseudoR2(model.A_small)
# [1] 0.4550205
# kres
model.A_small.Kres.env <-envelope(model.A_small,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A_small_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A_small.Kres.env, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_A_small_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_small,covname='geo_rock_im',elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_small_elevation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_small,covname='elevation_im',geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_small_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_small,covname='precipitation_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

# 7.2.2 earth-wall sites
# 1) large
model.B_large<-ppm(B_large_ppp~geo_earth_im + precipitation_im+ slope_im+ soil_recla_im,use.gam=TRUE)
pseudoR2(model.B_large)
# [1] 0.1382876
# kres
model.B_large.Kres.env <-envelope(model.B_large,fun=Kres,nsim=500,correction='best')
jpeg("model_B_large_Kres.jpg",width = 7, height = 7, units = "in", res = 1000)
plot(model.B_large.Kres.env, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_B_large_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_large,covname='precipitation_im',geo_earth_im=mean(geo_earth_im),slope_im=mean(slope_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_B_large_geo_earth.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_large,covname='geo_earth_im',precipitation_im=mean(precipitation_im),slope_im=mean(slope_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()


# 2) middle
model.B_middle<-ppm(B_middle_ppp~geo_earth_im + precipitation_im+ slope_im+ soil_recla_im,use.gam=TRUE)
pseudoR2(model.B_middle)
# [1] 0.1427192
# kres
model.B_middle.Kres.env <-envelope(model.B_middle,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "model_B_middle_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B_middle.Kres.env, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_B_middle_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_middle,covname='precipitation_im',geo_earth_im=mean(geo_earth_im),slope_im=mean(slope_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_B_middle_geo_earth.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_middle,covname='geo_earth_im',precipitation_im=mean(precipitation_im),slope_im=mean(slope_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()


# 3) small
model.B_small<-ppm(B_small_ppp~geo_earth_im + precipitation_im+ slope_im+ soil_recla_im,use.gam=TRUE)
pseudoR2(model.B_small)
#[1] 0.1510483 
# kres
model.B_small.Kres.env <-envelope(model.B_small,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "model_B_small_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B_small.Kres.env, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_B_small_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_small,covname='precipitation_im',geo_earth_im=mean(geo_earth_im),slope_im=mean(slope_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_B_small_geo_earth.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_small,covname='geo_earth_im',precipitation_im=mean(precipitation_im),slope_im=mean(slope_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()



# The following code will be rerunning the whole workflow to compare among different size categories as suggested by reviewer 2.
# 7.3 function form
# 7.3.1 stone-wall sites
# 7.3.1.1 large
# 1) elevation:  
ppmA_large_elevation<-ppm(A_large_ppp~elevation_im)
residuals.ppm(ppmA_large_elevation)
# 
gamA_large_elevation<-ppm(A_large_ppp~s(elevation_im),use.gam=TRUE)
residuals.ppm(gamA_large_elevation)
# 

# 2) slope:   
ppmA_large_slope<-ppm(A_large_ppp~slope_im)
residuals.ppm(ppmA_large_slope)
# 
gamA_large_slope<-ppm(A_large_ppp~s(slope_im),use.gam=TRUE)
residuals.ppm(gamA_large_slope)
# 

# 3) aspect:   
ppmA_large_aspect<-ppm(A_large_ppp~aspect_im)
residuals.ppm(ppmA_large_aspect)
# 
gamA_large_aspect<-ppm(A_large_ppp~s(aspect_im),use.gam=TRUE)
residuals.ppm(gamA_large_aspect)
# 

# 4) incoming visibility index:  
ppmA_large_vis_in<-ppm(A_large_ppp~vis_in_im)
residuals.ppm(ppmA_large_vis_in)
# 
gamA_large_vis_in<-ppm(A_large_ppp~s(vis_in_im),use.gam=TRUE)
residuals.ppm(gamA_large_vis_in)
# 

# 5) outgoing visibility index:  
ppmA_large_vis_out<-ppm(A_large_ppp~vis_out_im)
residuals.ppm(ppmA_large_vis_out)
# 
gamA_large_vis_out<-ppm(A_large_ppp~s(vis_out_im),use.gam=TRUE)
residuals.ppm(gamA_large_vis_out)
# 

# 6) landforms:  
ppmA_large_landforms<-ppm(A_large_ppp~landforms_recla_im)
residuals.ppm(ppmA_large_landforms)
# 

# 7) river:   
ppmA_large_river<-ppm(A_large_ppp~river_im)
residuals.ppm(ppmA_large_river)
# 
gamA_large_river<-ppm(A_large_ppp~s(river_im),use.gam=TRUE)
residuals.ppm(gamA_large_river)
# 

# 8) soil:  
ppmA_large_soil<-ppm(A_large_ppp~soil_recla_im)
residuals.ppm(ppmA_large_soil)
# 

# 9) geology
# building rock:  
ppmA_large_geo_rock <- ppm(A_large_ppp~geo_rock_im)
residuals.ppm(ppmA_large_geo_rock)
# 
gamA_large_geo_rock <- ppm(A_large_ppp~s(geo_rock_im), use.gam=TRUE)
residuals.ppm(gamA_large_geo_rock)
# 

# 10) modern annual precipitation:   
ppmA_large_precipitation<-ppm(A_large_ppp~precipitation_im)
residuals.ppm(ppmA_large_precipitation)
# 
gamA_large_precipitation<-ppm(A_large_ppp~s(precipitation_im),use.gam=TRUE)
residuals.ppm(gamA_large_precipitation)
# 

# 11) modern annual mean temperature:   
ppmA_large_temperature<-ppm(A_large_ppp~temperature_im)
residuals.ppm(ppmA_large_temperature)
# 
gamA_large_temperature<-ppm(A_large_ppp~s(temperature_im),use.gam=TRUE)
residuals.ppm(gamA_large_temperature)
# 

# 7.3.1.2 middle
# 1) elevation: 
ppmA_middle_elevation<-ppm(A_middle_ppp~elevation_im)
residuals.ppm(ppmA_middle_elevation)
# 
gamA_middle_elevation<-ppm(A_middle_ppp~s(elevation_im),use.gam=TRUE)
residuals.ppm(gamA_middle_elevation)
# 

# 2) slope:   
ppmA_middle_slope<-ppm(A_middle_ppp~slope_im)
residuals.ppm(ppmA_middle_slope)
# 
gamA_middle_slope<-ppm(A_middle_ppp~s(slope_im),use.gam=TRUE)
residuals.ppm(gamA_middle_slope)
# 

# 3) aspect:   
ppmA_middle_aspect<-ppm(A_middle_ppp~aspect_im)
residuals.ppm(ppmA_middle_aspect)
# 
gamA_middle_aspect<-ppm(A_middle_ppp~s(aspect_im),use.gam=TRUE)
residuals.ppm(gamA_middle_aspect)
# 

# 4) incoming visibility index:  
ppmA_middle_vis_in<-ppm(A_middle_ppp~vis_in_im)
residuals.ppm(ppmA_middle_vis_in)
# 
gamA_middle_vis_in<-ppm(A_middle_ppp~s(vis_in_im),use.gam=TRUE)
residuals.ppm(gamA_middle_vis_in)
# 

# 5) outgoing visibility index:  
ppmA_middle_vis_out<-ppm(A_middle_ppp~vis_out_im)
residuals.ppm(ppmA_middle_vis_out)
# 
gamA_middle_vis_out<-ppm(A_middle_ppp~s(vis_out_im),use.gam=TRUE)
residuals.ppm(gamA_middle_vis_out)
# 

# 6) landforms:  
ppmA_middle_landforms<-ppm(A_middle_ppp~landforms_recla_im)
residuals.ppm(ppmA_middle_landforms)
# 

# 7) river:   
ppmA_middle_river<-ppm(A_middle_ppp~river_im)
residuals.ppm(ppmA_middle_river)
# 
gamA_middle_river<-ppm(A_middle_ppp~s(river_im),use.gam=TRUE)
residuals.ppm(gamA_middle_river)
# 

# 8) soil:  
ppmA_middle_soil<-ppm(A_middle_ppp~soil_recla_im)
residuals.ppm(ppmA_middle_soil)
# 

# 9) geology
# building rock:  
ppmA_large_geo_rock <- ppm(A_large_ppp~geo_rock_im)
residuals.ppm(ppmA_large_geo_rock)
# 
gamA_large_geo_rock <- ppm(A_large_ppp~s(geo_rock_im), use.gam=TRUE)
residuals.ppm(gamA_large_geo_rock)
# 

# 10) modern annual precipitation:   
ppmA_middle_precipitation<-ppm(A_middle_ppp~precipitation_im)
residuals.ppm(ppmA_middle_precipitation)
# 
gamA_middle_precipitation<-ppm(A_middle_ppp~s(precipitation_im),use.gam=TRUE)
residuals.ppm(gamA_middle_precipitation)
# 

# 11) modern annual mean temperature:    
ppmA_middle_temperature<-ppm(A_middle_ppp~temperature_im)
residuals.ppm(ppmA_middle_temperature)
# 
gamA_middle_temperature<-ppm(A_middle_ppp~s(temperature_im),use.gam=TRUE)
residuals.ppm(gamA_middle_temperature)
#

# 7.3.1.3 small
# 1) elevation:  
ppmA_small_elevation<-ppm(A_small_ppp~elevation_im)
residuals.ppm(ppmA_small_elevation)
# 
gamA_small_elevation<-ppm(A_small_ppp~s(elevation_im),use.gam=TRUE)
residuals.ppm(gamA_small_elevation)
# 

# 2) slope:   
ppmA_small_slope<-ppm(A_small_ppp~slope_im)
residuals.ppm(ppmA_small_slope)
# 
gamA_small_slope<-ppm(A_small_ppp~s(slope_im),use.gam=TRUE)
residuals.ppm(gamA_small_slope)
# 

# 3) aspect:   
ppmA_small_aspect<-ppm(A_small_ppp~aspect_im)
residuals.ppm(ppmA_small_aspect)
# 
gamA_small_aspect<-ppm(A_small_ppp~s(aspect_im),use.gam=TRUE)
residuals.ppm(gamA_small_aspect)
# 

# 4) incoming visibility index:  
ppmA_small_vis_in<-ppm(A_small_ppp~vis_in_im)
residuals.ppm(ppmA_small_vis_in)
# 
gamA_small_vis_in<-ppm(A_small_ppp~s(vis_in_im),use.gam=TRUE)
residuals.ppm(gamA_small_vis_in)
# 

# 5) outgoing visibility index:  
ppmA_small_vis_out<-ppm(A_small_ppp~vis_out_im)
residuals.ppm(ppmA_small_vis_out)
# 
gamA_small_vis_out<-ppm(A_small_ppp~s(vis_out_im),use.gam=TRUE)
residuals.ppm(gamA_small_vis_out)
# 

# 6) landforms:  
ppmA_small_landforms<-ppm(A_small_ppp~landforms_recla_im)
residuals.ppm(ppmA_small_landforms)
# 

# 7) river:   
ppmA_small_river<-ppm(A_small_ppp~river_im)
residuals.ppm(ppmA_small_river)
# 
gamA_small_river<-ppm(A_small_ppp~s(river_im),use.gam=TRUE)
residuals.ppm(gamA_small_river)
# 

# 8) soil:  
ppmA_small_soil<-ppm(A_small_ppp~soil_recla_im)
residuals.ppm(ppmA_small_soil)
# 

# 9) geology
# building rock:  
ppmA_small_geo_rock <- ppm(A_small_ppp~geo_rock_im)
residuals.ppm(ppmA_small_geo_rock)
# 
gamA_small_geo_rock <- ppm(A_small_ppp~s(geo_rock_im), use.gam=TRUE)
residuals.ppm(gamA_small_geo_rock)
# 

# 10) modern annual precipitation:   
ppmA_small_precipitation<-ppm(A_small_ppp~precipitation_im)
residuals.ppm(ppmA_small_precipitation)
# 
gamA_small_precipitation<-ppm(A_small_ppp~s(precipitation_im),use.gam=TRUE)
residuals.ppm(gamA_small_precipitation)
# 

# 11) modern annual mean temperature:   
ppmA_small_temperature<-ppm(A_small_ppp~temperature_im)
residuals.ppm(ppmA_small_temperature)
# 

gamA_small_temperature<-ppm(A_small_ppp~s(temperature_im),use.gam=TRUE)
residuals.ppm(gamA_small_temperature)
#

# 7.3.2 earth-wall sites
# 7.3.2.1 large
# 1) elevation:   
ppmB_large_elevation<-ppm(B_large_ppp~elevation_im)
residuals.ppm(ppmB_large_elevation)
# 
gamB_large_elevation<-ppm(B_large_ppp~s(elevation_im),use.gam=TRUE)
residuals.ppm(gamB_large_elevation)
# 

# 2) slope:   
ppmB_large_slope<-ppm(B_large_ppp~slope_im)
residuals.ppm(ppmB_large_slope)
# 
gamB_large_slope<-ppm(B_large_ppp~s(slope_im),use.gam=TRUE)
residuals.ppm(gamB_large_slope)
# 

# 3) aspect:   
ppmB_large_aspect<-ppm(B_large_ppp~aspect_im)
residuals.ppm(ppmB_large_aspect)
# 
gamB_large_aspect<-ppm(B_large_ppp~s(aspect_im),use.gam=TRUE)
residuals.ppm(gamB_large_aspect)
# 

# 4) incoming visibility index:   
ppmB_large_vis_in<-ppm(B_large_ppp~vis_in_im)
residuals.ppm(ppmB_large_vis_in)
# 
gamB_large_vis_in<-ppm(B_large_ppp~s(vis_in_im),use.gam=TRUE)
residuals.ppm(gamB_large_vis_in)
# 

# 5) outgoing visibility index:  
ppmB_large_vis_out<-ppm(B_large_ppp~vis_out_im)
residuals.ppm(ppmB_large_vis_out)
# 
gamB_large_vis_out<-ppm(B_large_ppp~s(vis_out_im),use.gam=TRUE)
residuals.ppm(gamB_large_vis_out)
# 

# 6) landforms:  
ppmB_large_landforms<-ppm(B_large_ppp~landforms_recla_im)
residuals.ppm(ppmB_large_landforms)
# 

# 7) river:   
ppmB_large_river<-ppm(B_large_ppp~river_im)
residuals.ppm(ppmB_large_river)
# 
gamB_large_river<-ppm(B_large_ppp~s(river_im),use.gam=TRUE)
residuals.ppm(gamB_large_river)
# 

# 8) soil:  
ppmB_large_soil<-ppm(B_large_ppp~soil_recla_im)
residuals.ppm(ppmB_large_soil)
# 

# 9) geology
# building earth:  
ppmB_large_geo_earth <- ppm(B_large_ppp~geo_earth_im)
residuals.ppm(ppmB_large_geo_earth)
# 
gamB_large_geo_earth <- ppm(B_large_ppp~s(geo_earth_im), use.gam=TRUE)
residuals.ppm(gamB_large_geo_earth)
# 

# 10) modern annual precipitation:    
ppmB_large_precipitation<-ppm(B_large_ppp~precipitation_im)
residuals.ppm(ppmB_large_precipitation)
# 
gamB_large_precipitation<-ppm(B_large_ppp~s(precipitation_im),use.gam=TRUE)
residuals.ppm(gamB_large_precipitation)
# 

# 11) modern annual mean temperature:   
ppmB_large_temperature<-ppm(B_large_ppp~temperature_im)
residuals.ppm(ppmB_large_temperature)
# 
gamB_large_temperature<-ppm(B_large_ppp~s(temperature_im),use.gam=TRUE)
residuals.ppm(gamB_large_temperature)
# 
# 7.3.2.2 middle
# 1) elevation: 
ppmB_middle_elevation<-ppm(B_middle_ppp~elevation_im)
residuals.ppm(ppmB_middle_elevation)
# 
gamB_middle_elevation<-ppm(B_middle_ppp~s(elevation_im),use.gam=TRUE)
residuals.ppm(gamB_middle_elevation)
# 

# 2) slope:   
ppmB_middle_slope<-ppm(B_middle_ppp~slope_im)
residuals.ppm(ppmB_middle_slope)
# 
gamB_middle_slope<-ppm(B_middle_ppp~s(slope_im),use.gam=TRUE)
residuals.ppm(gamB_middle_slope)
# 

# 3) aspect:   
ppmB_middle_aspect<-ppm(B_middle_ppp~aspect_im)
residuals.ppm(ppmB_middle_aspect)
# 
gamB_middle_aspect<-ppm(B_middle_ppp~s(aspect_im),use.gam=TRUE)
residuals.ppm(gamB_middle_aspect)
# 

# 4) incoming visibility index:   
ppmB_middle_vis_in<-ppm(B_middle_ppp~vis_in_im)
residuals.ppm(ppmB_middle_vis_in)
# 
gamB_middle_vis_in<-ppm(B_middle_ppp~s(vis_in_im),use.gam=TRUE)
residuals.ppm(gamB_middle_vis_in)
# 

# 5) outgoing visibility index:   
ppmB_middle_vis_out<-ppm(B_middle_ppp~vis_out_im)
residuals.ppm(ppmB_middle_vis_out)
# 
gamB_middle_vis_out<-ppm(B_middle_ppp~s(vis_out_im),use.gam=TRUE)
residuals.ppm(gamB_middle_vis_out)
# 

# 6) landforms: 
ppmB_middle_landforms<-ppm(B_middle_ppp~landforms_recla_im)
residuals.ppm(ppmB_middle_landforms)
# 

# 7) river:   
ppmB_middle_river<-ppm(B_middle_ppp~river_im)
residuals.ppm(ppmB_middle_river)
# 
gamB_middle_river<-ppm(B_middle_ppp~s(river_im),use.gam=TRUE)
residuals.ppm(gamB_middle_river)
# 

# 8) soil: 
ppmB_middle_soil<-ppm(B_middle_ppp~soil_recla_im)
residuals.ppm(ppmB_middle_soil)
# 

# 9) geology
# building earth:  
ppmB_middle_geo_earth <- ppm(B_middle_ppp~geo_earth_im)
residuals.ppm(ppmB_middle_geo_earth)
# 
gamB_middle_geo_earth <- ppm(B_middle_ppp~s(geo_earth_im), use.gam=TRUE)
residuals.ppm(gamB_middle_geo_earth)
# 

# 10) modern annual precipitation:   
ppmB_middle_precipitation<-ppm(B_middle_ppp~precipitation_im)
residuals.ppm(ppmB_middle_precipitation)
# 
gamB_middle_precipitation<-ppm(B_middle_ppp~s(precipitation_im),use.gam=TRUE)
residuals.ppm(gamB_middle_precipitation)
# 

# 11) modern annual mean temperature:   
ppmB_middle_temperature<-ppm(B_middle_ppp~temperature_im)
residuals.ppm(ppmB_middle_temperature)
# 
gamB_middle_temperature<-ppm(B_middle_ppp~s(temperature_im),use.gam=TRUE)
residuals.ppm(gamB_middle_temperature)
#

# 7.3.2.3 small
# 1) elevation: 
ppmB_small_elevation<-ppm(B_small_ppp~elevation_im)
residuals.ppm(ppmB_small_elevation)
# 
gamB_small_elevation<-ppm(B_small_ppp~s(elevation_im),use.gam=TRUE)
residuals.ppm(gamB_small_elevation)
# 

# 2) slope:   
ppmB_small_slope<-ppm(B_small_ppp~slope_im)
residuals.ppm(ppmB_small_slope)
# 
gamB_small_slope<-ppm(B_small_ppp~s(slope_im),use.gam=TRUE)
residuals.ppm(gamB_small_slope)
# 

# 3) aspect:   
ppmB_small_aspect<-ppm(B_small_ppp~aspect_im)
residuals.ppm(ppmB_small_aspect)
# 
gamB_small_aspect<-ppm(B_small_ppp~s(aspect_im),use.gam=TRUE)
residuals.ppm(gamB_small_aspect)
# 

# 4) incoming visibility index:   
ppmB_small_vis_in<-ppm(B_small_ppp~vis_in_im)
residuals.ppm(ppmB_small_vis_in)
# 
gamB_small_vis_in<-ppm(B_small_ppp~s(vis_in_im),use.gam=TRUE)
residuals.ppm(gamB_small_vis_in)
# 

# 5) outgoing visibility index:  
ppmB_small_vis_out<-ppm(B_small_ppp~vis_out_im)
residuals.ppm(ppmB_small_vis_out)
# 
gamB_small_vis_out<-ppm(B_small_ppp~s(vis_out_im),use.gam=TRUE)
residuals.ppm(gamB_small_vis_out)
# 

# 6) landforms:  
ppmB_small_landforms<-ppm(B_small_ppp~landforms_recla_im)
residuals.ppm(ppmB_small_landforms)
# 

# 7) river:   
ppmB_small_river<-ppm(B_small_ppp~river_im)
residuals.ppm(ppmB_small_river)
# 
gamB_small_river<-ppm(B_small_ppp~s(river_im),use.gam=TRUE)
residuals.ppm(gamB_small_river)
# 

# 8) soil:  
ppmB_small_soil<-ppm(B_small_ppp~soil_recla_im)
residuals.ppm(ppmB_small_soil)
# 

# 9) geology
# building earth:   
ppmB_small_geo_earth <- ppm(B_small_ppp~geo_earth_im)
residuals.ppm(ppmB_small_geo_earth)
# 
gamB_small_geo_earth <- ppm(B_small_ppp~s(geo_earth_im), use.gam=TRUE)
residuals.ppm(gamB_small_geo_earth)
# 

# 10) modern annual precipitation:   
ppmB_small_precipitation<-ppm(B_small_ppp~precipitation_im)
residuals.ppm(ppmB_small_precipitation)
# 
gamB_small_precipitation<-ppm(B_small_ppp~s(precipitation_im),use.gam=TRUE)
residuals.ppm(gamB_small_precipitation)
# 

# 11) modern annual mean temperature:   
ppmB_small_temperature<-ppm(B_small_ppp~temperature_im)
residuals.ppm(ppmB_small_temperature)
# 
gamB_small_temperature<-ppm(B_small_ppp~s(temperature_im),use.gam=TRUE)
residuals.ppm(gamB_small_temperature)
#

# 7.4 relative variable importance
# 7.4.1 stone-wall sites
# 1) large
global.model.A_large<-ppm(A_large_ppp~s(elevation_im)+s(slope_im)+aspect_im+s(vis_in_im)+s(vis_out_im)+landforms_recla_im+river_im+soil_recla_im+geo_rock_im+precipitation_im+s(temperature_im), use.gam=TRUE)
all.model.A_large <- dredge(global.model.A_large,evaluate=FALSE)
all.model2.A_large <- lapply(all.model.A_large,function(x){as.character(x)[4]})
candidates.aic.A_large  <- numeric(length=length(all.model2.A_large))
for (i in 1:length(candidates.aic.A_large))
{
	candidates.aic.A_large[i] <- AICc(ppm(Q=as.formula(paste0('A_large_ppp',all.model2.A_large[[i]])),use.gam=TRUE))
}

params.A_large <- getAllTerms(global.model.A_large) |> as.character()
permuts.A_large <- vector('list', length=length(params.A_large))
for (i in 1:length(permuts.A_large)) {
  permuts.A_large[[i]] <- c(0, 1)
}
param.comb.A_large <- expand.grid(permuts.A_large)
colnames(param.comb.A_large) <- params.A_large
raw_weights <- Weights(candidates.aic.A_large)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.A_large$weights <- rounded_weights / sum(rounded_weights)

param.importance.A_large <- data.frame(params.A_large=params.A_large, aic.imp=NA)

for (i in 1:length(params.A_large)) {
  param.importance.A_large$aic.imp[i] <- sum(param.comb.A_large[, i] * param.comb.A_large$weights)
}
param.importance.A_large


# 2) middle
global.model.A_middle<-ppm(A_middle_ppp~s(elevation_im)+s(slope_im)+aspect_im+s(vis_in_im)+vis_out_im+landforms_recla_im+river_im+soil_recla_im+geo_rock_im+s(precipitation_im)+s(temperature_im), use.gam=TRUE)
all.model.A_middle <- dredge(global.model.A_middle,evaluate=FALSE)
all.model2.A_middle <- lapply(all.model.A_middle,function(x){as.character(x)[4]})
candidates.aic.A_middle  <- numeric(length=length(all.model2.A_middle))
for (i in 1:length(candidates.aic.A_middle))
{
	candidates.aic.A_middle[i] <- AICc(ppm(Q=as.formula(paste0('A_middle_ppp',all.model2.A_middle[[i]])),use.gam=TRUE))
}

params.A_middle <- getAllTerms(global.model.A_middle) |> as.character()
permuts.A_middle <- vector('list', length=length(params.A_middle))
for (i in 1:length(permuts.A_middle)) {
  permuts.A_middle[[i]] <- c(0, 1)
}
param.comb.A_middle <- expand.grid(permuts.A_middle)
colnames(param.comb.A_middle) <- params.A_middle
raw_weights <- Weights(candidates.aic.A_middle)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.A_middle$weights <- rounded_weights / sum(rounded_weights)

param.importance.A_middle <- data.frame(params.A_middle=params.A_middle, aic.imp=NA)

for (i in 1:length(params.A_middle)) {
  param.importance.A_middle$aic.imp[i] <- sum(param.comb.A_middle[, i] * param.comb.A_middle$weights)
}
param.importance.A_middle


# 3) small
global.model.A_small<-ppm(A_small_ppp~s(elevation_im)+s(slope_im)+aspect_im+s(vis_in_im)+s(vis_out_im)+landforms_recla_im+s(river_im)+soil_recla_im+geo_rock_im+precipitation_im+s(temperature_im), use.gam=TRUE)
all.model.A_small <- dredge(global.model.A_small,evaluate=FALSE)
all.model2.A_small <- lapply(all.model.A_small,function(x){as.character(x)[4]})
candidates.aic.A_small  <- numeric(length=length(all.model2.A_small))
for (i in 1:length(candidates.aic.A_small))
{
	candidates.aic.A_small[i] <- AICc(ppm(Q=as.formula(paste0('A_small_ppp',all.model2.A_small[[i]])),use.gam=TRUE))
}

params.A_small <- getAllTerms(global.model.A_small) |> as.character()
permuts.A_small <- vector('list', length=length(params.A_small))
for (i in 1:length(permuts.A_small)) {
  permuts.A_small[[i]] <- c(0, 1)
}
param.comb.A_small <- expand.grid(permuts.A_small)
colnames(param.comb.A_small) <- params.A_small
raw_weights <- Weights(candidates.aic.A_small)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.A_small$weights <- rounded_weights / sum(rounded_weights)

param.importance.A_small <- data.frame(params.A_small=params.A_small, aic.imp=NA)

for (i in 1:length(params.A_small)) {
  param.importance.A_small$aic.imp[i] <- sum(param.comb.A_small[, i] * param.comb.A_small$weights)
}
param.importance.A_small

# 
# 7.4.2 earth-wall sites
# 1) large
global.model.B_large<-ppm(B_large_ppp~elevation_im+slope_im+aspect_im+vis_in_im+vis_out_im+landforms_recla_im+river_im+soil_recla_im+geo_earth_im+precipitation_im+temperature_im,use.gam=TRUE)
all.model.B_large <- dredge(global.model.B_large,evaluate=FALSE)
all.model2.B_large <- lapply(all.model.B_large,function(x){as.character(x)[4]})
candidates.aic.B_large  <- numeric(length=length(all.model2.B_large))
for (i in 1:length(candidates.aic.B_large))
{
	candidates.aic.B_large[i] <- AICc(ppm(Q=as.formula(paste0('B_large_ppp',all.model2.B_large[[i]])),use.gam=TRUE))
}

params.B_large <- getAllTerms(global.model.B_large) |> as.character()
permuts.B_large <- vector('list', length=length(params.B_large))
for (i in 1:length(permuts.B_large)) {
  permuts.B_large[[i]] <- c(0, 1)
}
param.comb.B_large <- expand.grid(permuts.B_large)
colnames(param.comb.B_large) <- params.B_large
raw_weights <- Weights(candidates.aic.B_large)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.B_large$weights <- rounded_weights / sum(rounded_weights)

param.importance.B_large <- data.frame(params.B_large=params.B_large, aic.imp=NA)

for (i in 1:length(params.B_large)) {
  param.importance.B_large$aic.imp[i] <- sum(param.comb.B_large[, i] * param.comb.B_large$weights)
}
param.importance.B_large


# 2) middle
global.model.B_middle<-ppm(B_middle_ppp~elevation_im+slope_im+aspect_im+vis_in_im+vis_out_im+landforms_recla_im+river_im+soil_recla_im+geo_earth_im+precipitation_im+temperature_im,use.gam=TRUE)
all.model.B_middle <- dredge(global.model.B_middle,evaluate=FALSE)
all.model2.B_middle <- lapply(all.model.B_middle,function(x){as.character(x)[4]})
candidates.aic.B_middle  <- numeric(length=length(all.model2.B_middle))
for (i in 1:length(candidates.aic.B_middle))
{
	candidates.aic.B_middle[i] <- AICc(ppm(Q=as.formula(paste0('B_middle_ppp',all.model2.B_middle[[i]])),use.gam=TRUE))
}

params.B_middle <- getAllTerms(global.model.B_middle) |> as.character()
permuts.B_middle <- vector('list', length=length(params.B_middle))
for (i in 1:length(permuts.B_middle)) {
  permuts.B_middle[[i]] <- c(0, 1)
}
param.comb.B_middle <- expand.grid(permuts.B_middle)
colnames(param.comb.B_middle) <- params.B_middle
raw_weights <- Weights(candidates.aic.B_middle)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.B_middle$weights <- rounded_weights / sum(rounded_weights)

param.importance.B_middle <- data.frame(params.B_middle=params.B_middle, aic.imp=NA)

for (i in 1:length(params.B_middle)) {
  param.importance.B_middle$aic.imp[i] <- sum(param.comb.B_middle[, i] * param.comb.B_middle$weights)
}
param.importance.B_middle


# 3) small
global.model.B_small<-ppm(B_small_ppp~elevation_im+slope_im+aspect_im+vis_in_im+vis_out_im+landforms_recla_im+river_im+soil_recla_im+geo_earth_im+precipitation_im+temperature_im,use.gam=TRUE)
all.model.B_small <- dredge(global.model.B_small,evaluate=FALSE)
all.model2.B_small <- lapply(all.model.B_small,function(x){as.character(x)[4]})
candidates.aic.B_small  <- numeric(length=length(all.model2.B_small))
for (i in 1:length(candidates.aic.B_small))
{
	candidates.aic.B_small[i] <- AICc(ppm(Q=as.formula(paste0('B_small_ppp',all.model2.B_small[[i]])),use.gam=TRUE))
}

params.B_small <- getAllTerms(global.model.B_small) |> as.character()
permuts.B_small <- vector('list', length=length(params.B_small))
for (i in 1:length(permuts.B_small)) {
  permuts.B_small[[i]] <- c(0, 1)
}
param.comb.B_small <- expand.grid(permuts.B_small)
colnames(param.comb.B_small) <- params.B_small
raw_weights <- Weights(candidates.aic.B_small)
normalized_weights <- raw_weights / sum(raw_weights)
rounded_weights <- round(normalized_weights, 4)
param.comb.B_small$weights <- rounded_weights / sum(rounded_weights)

param.importance.B_small <- data.frame(params.B_small=params.B_small, aic.imp=NA)

for (i in 1:length(params.B_small)) {
  param.importance.B_small$aic.imp[i] <- sum(param.comb.B_small[, i] * param.comb.B_small$weights)
}
param.importance.B_small


# 7.5 final model and goodness of fit
# 7.5.1 stone-wall sites
# 1) large
model.A2_large<-ppm(A_large_ppp~geo_rock_im+precipitation_im,use.gam=TRUE)
pseudoR2(model.A2_large)
# 

# kres
model.A2_large.Kres.env <-envelope(model.A2_large,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A2_large_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A2_large.Kres.env, main = NULL)
dev.off()
global.model.A_large.Kres.env <-envelope(global.model.A_large,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "global_model_A_large_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(global.model.A_large.Kres.env,main=NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_A2_large_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A2_large,covname='geo_rock_im',precipitation_im = mean(precipitation_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A2_large_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A2_large,covname='precipitation_im',geo_rock_im = mean(geo_rock_im), se.fit=TRUE),main=NULL)
dev.off()

# 2) middle
model.A2_middle<-ppm(A_middle_ppp~s(elevation_im)+s(precipitation_im)+vis_out_im,use.gam=TRUE)
pseudoR2(model.A2_middle)
# 

# kres
model.A2_middle.Kres.env <-envelope(model.A2_middle,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A2_middle_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A2_middle.Kres.env, main = NULL)
dev.off()
global.model.A_middle.Kres.env <-envelope(global.model.A_middle,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "global_model_A_middle_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(global.model.A_middle.Kres.env,main=NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_A2_middle_elevation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A2_middle,covname='elevation_im',precipitation_im = mean(precipitation_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A2_middle_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A2_middle,covname='precipitation_im',elevation_im=mean(elevation_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A2_middle_vis_out.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A2_middle,covname='vis_out_im',elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im),se.fit=TRUE),main=NULL)
dev.off()

# 3) small
model.A2_small<-ppm(A_small_ppp~landforms_recla_im+precipitation_im+s(elevation_im)+s(river_im),use.gam=TRUE)
pseudoR2(model.A2_small)
# 
print(model.A2_small)


# kres
model.A2_small.Kres.env <-envelope(model.A2_small,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "model_A2_small_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A2_small.Kres.env, main = NULL)
dev.off()
global.model.A_small.Kres.env <-envelope(global.model.A_small,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "3_size_categories", "global_model_A_small_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(global.model.A_small.Kres.env,main=NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_A2_small_landforms.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff=effectfun(
    model.A2_small,
    covname="landforms_recla_im",
    landforms_recla_im=landforms_recla_im_new, 
    precipitation_im=mean(precipitation_im), 
    elevation_im=mean(elevation_im), 
    river_im=mean(river_im), 
    se.fit=TRUE
)
barplot(eff$lambda,names.arg=eff$landforms_recla_im,main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A2_small_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff = effectfun(
  model.A2_small,
  covname = "precipitation_im",
  landforms_recla_im = factor(c("valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")),
  elevation_im = mean(elevation_im),
  river_im = mean(river_im),
  se.fit = TRUE
)
plot(eff,main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A2_small_elevation.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff = effectfun(
  model.A2_small,
  covname = "elevation_im",
  landforms_recla_im = factor(c("valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")),
  precipitation_im = mean(precipitation_im),
  river_im = mean(river_im),
  se.fit = TRUE
)
plot(eff,main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A2_small_river.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff = effectfun(
  model.A2_small,
  covname = "river_im",
  landforms_recla_im = factor(c("valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")),
  precipitation_im = mean(precipitation_im),
  elevation_im = mean(elevation_im),
  se.fit = TRUE
)
plot(eff,main=NULL)
dev.off()

# 7.5.2 earth-wall sites
# 1) large
model.B2_large<-ppm(B_large_ppp~precipitation_im+soil_recla_im,use.gam=TRUE)
pseudoR2(model.B2_large)
# 

# kres
model.B2_large.Kres.env <-envelope(model.B2_large,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "model_B2_large_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B2_large.Kres.env, main = NULL)
dev.off()
global.model.B_large.Kres.env <-envelope(global.model.B_large,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "global_model_B_large_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(global.model.B_large.Kres.env,main=NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_B2_large_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B2_large, covname='precipitation_im', soil_recla_im=as.factor(c("absence", "presence")),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_B2_large_soil.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff_soil <- effectfun(model.B2_large,covname='soil_recla_im',precipitation_im=mean(precipitation_im), se.fit=TRUE)
barplot(eff_soil$lambda,names.arg=eff_soil$soil_recla_im,main=NULL)
dev.off()

# 2) middle
model.B2_middle<-ppm(B_middle_ppp~precipitation_im,use.gam=TRUE)
pseudoR2(model.B2_middle)
# 

# kres
model.B2_middle.Kres.env <-envelope(model.B2_middle,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "model_B2_middle_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B2_middle.Kres.env, main = NULL)
dev.off()
global.model.B_middle.Kres.env <-envelope(global.model.B_middle,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "global_model_B_middle_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(global.model.B_middle.Kres.env,main=NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_B2_middle_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B2_middle,se.fit=TRUE),main=NULL)
dev.off()

# 3) small
model.B2_small<-ppm(B_small_ppp~precipitation_im,use.gam=TRUE)
pseudoR2(model.B2_small)
# 

# kres
model.B2_small.Kres.env <-envelope(model.B2_small,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "model_B2_small_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B2_small.Kres.env, main = NULL)
dev.off()
global.model.B_small.Kres.env <-envelope(global.model.B_small,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "3_size_categories", "global_model_B_small_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(global.model.B_small.Kres.env,main=NULL)
dev.off()

# effectfun
jpeg(here("Output", "3_size_categories", "effectfun_B2_small_precipitation_im.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B2_small,se.fit=TRUE),main=NULL)
dev.off()


























