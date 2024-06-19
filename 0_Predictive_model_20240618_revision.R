# R script
## article title: A comparative analyses of stone- and earth-wall settlement locations of the Lower Xiajiadian Culture in Aohan Banner, China
## journal name: Journal of Archaeological Method and Theory
## author names: Xuan Zhang; Yukun Zhang; Lifeng Tan; Enrico R. Crema; Yanguo Tian; Ze Wang
## corresponding authors: corresponding author at School of Architecture, Tianjin University, China. E-mail address: 961295@tju.edu.cn (Yukun Zhang); corresponding author at McDonald Institute for Archaeological Research, Department of Archaeology,  University of Cambridge, UK. E-mail address: erc62@cam.ac.uk (Enrico R. Crema).


# 0. Prepare the packages
install.packages("usethis")
library(usethis)
usethis::create_project(".")
install.packages("here")
library(here)
install.packages("sf")
install.packages("maptools")
install.packages("spatstat")
devtools::install_github("spatstat/spatstat.model")
install.packages("raster")
install.packages("terra")
install.packages("sp")
install.packages("MuMIn")
install.packages("MASS")
library(sf)
library(maptools)
library(spatstat)
library(spatstat.model)
library(raster)
library(terra)
library(sp)
library(MuMIn)
library(MASS)
# 1. Read the objects and convert them to spatstat format ----
# The observed window
Aohan <-read_sf(here("GIS_layers","1_1_Aohan_qgisproj.shp"))
Aohan_owin=as.owin(Aohan$geometry)
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
elevation_im <- as.im(DEM, Aohan_owin)
class(elevation_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_1_elevation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(elevation_im, main = NULL)
dev.off()

# 1.2 slope
slope <- mask(raster(here("GIS_layers","1_2_slope_qgisproj.tif")), Aohan)
res(slope)
# [1] 26.99211 26.99211
slope_im<-as.im(slope, Aohan_owin)
class(slope_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_2_slope.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(slope_im, main = NULL)
dev.off()

# 1.3 aspect
aspect <- mask(raster(here("GIS_layers","1_3_aspect_qgisproj.tif")), Aohan)
res(aspect)
# [1] 26.99211 26.99211
north_south <- cos(aspect)*(-1)
aspect_im <- as.im(north_south, Aohan_owin)
class(aspect_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_3_aspect.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(aspect_im, main = NULL)
dev.off()

# 1.4 incoming visibility index
visibility_in<-mask(raster(here("GIS_layers","1_4_in_vis_index_qgisproj.tif"),window=Aohan_owin), Aohan)
res(visibility_in)
# [1] 26.99211 26.99211
vis_in_im <- as.im(visibility_in, Aohan_owin)
plot(vis_in_im)
class(vis_in_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_4_vis_in.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(vis_in_im, main = NULL)
dev.off()


# 1.5 outgoing visbility index
visibility_out<-mask(raster(here("GIS_layers","1_5_out_vis_index_qgisproj.tif"),window=Aohan_owin), Aohan)
res(visibility_out)
# [1] 26.99211 26.99211
vis_out_im <- as.im(visibility_out, Aohan_owin)
plot(vis_out_im)
class(vis_out_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_5_vis_out.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(vis_out_im, main = NULL)
dev.off()




# 1.6 landforms
landforms <- mask(raster(here("GIS_layers","1_6_landforms_qgisproj_reclass.tif")), Aohan)
res(landforms)
# [1] 26.99211 26.99211
landforms_im<-as.im(landforms, Aohan_owin)
landforms_recla_im <- cut(landforms_im$v, 
                          breaks = c(-Inf, 0, 1, 2, 3, Inf), 
                          labels = c("flat", "pit", "valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder"), 
                          include.lowest = TRUE)

landforms_recla_im <- im(factor(landforms_recla_im, levels = c("flat", "pit", "valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder")), 
                         xcol = landforms_im$xcol, 
                         yrow = landforms_im$yrow)

class(landforms_recla_im$v)
# [1] "factor"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_6_landforms.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(landforms_recla_im, main = "landforms_recla_im", 
     ribside = "bottom", 
     ribsep = 0.05, 
     ribargs = list(at = 1:5, 
		    labels = c("flat", "pit", "valley/footslope", "spur/slope/hollow", "peak/ridge/shoulder"), 
		    cex.axis = 0.8, 
		    width = 1.5)) 

dev.off()

# 1.7 river
river <- mask(raster(here("GIS_layers","1_7_river_level_1_to_5_clipped_qgisproj2_distance.tif")), Aohan)
res(river)
# [1] 26.99211 26.99211
river_im<-as.im(river, Aohan_owin)
class(river_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_7_river.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(river_im, main = NULL)
dev.off()

# 1.8 soil
soil <- mask(raster(here("GIS_layers","1_8_soil_reclass_qgisproj.tif")), Aohan)
soil_resample <- resample(soil, DEM, method='bilinear')
res(soil_resample)
# [1] 26.99211 26.99211
class(values(soil_resample))
# [1] "numeric"
soil_im<-as.im(soil_resample, Aohan_owin)
class(soil_im$v)
# [1] "matrix" "array"
soil_categorical <- function(x) {
   as.factor(ifelse(x > 0, "presence", "absence"))
 }
soil_recla_im <- eval.im(soil_categorical(soil_im))
class(soil_recla_im$v)
# [1] "factor"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_8_soil.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(soil_recla_im, main = NULL)
dev.off()

# 1.9 geology
geo_rock <- mask(raster(here("GIS_layers","1_9_building_rock_qgisproj_distance.tif")), Aohan)
geo_rock_resample <- resample(geo_rock, DEM, method='bilinear')
res(geo_rock_resample)
# [1] 26.99211 26.99211
geo_rock_im<-as.im(geo_rock_resample, Aohan_owin)
class(geo_rock_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_9_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(geo_rock_im, main = NULL)
dev.off()

geo_earth <- mask(raster(here("GIS_layers","1_9_building_earth_qgisproj_distance.tif")), Aohan)
geo_earth_resample <- resample(geo_earth, DEM, method='bilinear')
res(geo_earth_resample)
# [1] 26.99211 26.99211
geo_earth_im<-as.im(geo_earth_resample, Aohan_owin)
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
precipitation_im <- as.im(precipitation_resample, Aohan_owin)
class(precipitation_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_10_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(precipitation_im, main = NULL)
dev.off()

# 1.11 temperature 
temperature <- mask(raster(here("GIS_layers","1_11_mean_tem_clipped_qgisproj.tif")), Aohan)
temperature_resample <- resample(temperature, DEM, method='bilinear')
res(temperature_resample)
# [1] 26.99211 26.99211
temperature_im <- as.im(temperature_resample, Aohan_owin)
class(temperature_im$v)
# [1] "matrix" "array"
# code for Figure 4
jpeg(here("Output", "1_covariates", "1_11_temperature.jpg"),width = 7, height = 7, units = "in", res = 1200)
plot(temperature_im, main = NULL)
dev.off()


# 3.Function form selection: LM or GAM
# 3.1 stone-wall sites
# 1) elevation: GAM 
ppm_A_elevation<-ppm(TypeA_ppp~elevation_im)
AICc(ppm_A_elevation)
# [1] 17088.05
gam_A_elevation<-ppm(TypeA_ppp~s(elevation_im),use.gam=TRUE)
AICc(gam_A_elevation)
# [1] 16931.66

# 2) slope: GAM 
ppm_A_slope<-ppm(TypeA_ppp~slope_im)
AICc(ppm_A_slope)
# [1] 17406.61
gam_A_slope<-ppm(TypeA_ppp~s(slope_im),use.gam=TRUE)
AICc(gam_A_slope)
# [1] 17347.26


# 3) aspect: LM 
ppm_A_aspect<-ppm(TypeA_ppp~aspect_im)
AICc(ppm_A_aspect)
# [1] 17519.14
gam_A_aspect<-ppm(TypeA_ppp~s(aspect_im),use.gam=TRUE)
AICc(gam_A_aspect)
# [1] 17535.57

# 4) incoming visibility index: GAM
ppm_A_vis_in<-ppm(TypeA_ppp~vis_in_im)
AICc(ppm_A_vis_in)
# [1] 17140.85
gam_A_vis_in<-ppm(TypeA_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gam_A_vis_in)
# [1] 17032.98

# 5) outgoing visibility index: GAM
ppm_A_vis_out<-ppm(TypeA_ppp~vis_out_im)
AICc(ppm_A_vis_out)
# [1] 17097.26
gam_A_vis_out<-ppm(TypeA_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gam_A_vis_out)
# [1] 17016.48

# 6) landforms: LM
ppm_A_landforms<-ppm(TypeA_ppp~landforms_recla_im)
AICc(ppm_A_landforms)
# [1] 17415.43


# 7) river: GAM 
ppm_A_river<-ppm(TypeA_ppp~river_im)
AICc(ppm_A_river)
# [1] 17504.86
gam_A_river<-ppm(TypeA_ppp~s(river_im),use.gam=TRUE)
AICc(gam_A_river)
# [1] 17493.21


# 8) soil: LM
ppm_A_soil<-ppm(TypeA_ppp~soil_recla_im)
AICc(ppm_A_soil)
# [1] 17310.82


# 9) geology
# building rock: LM
ppm_A_geo_rock <- ppm(TypeA_ppp~geo_rock_im)
AICc(ppm_A_geo_rock)
# [1] 17186.43
gam_A_geo_rock <- ppm(TypeA_ppp~s(geo_rock_im), use.gam=TRUE)
AICc(gam_A_geo_rock)
# [1] 17202.86


# 10) modern annual precipitation: GAM 
ppm_A_precipitation<-ppm(TypeA_ppp~precipitation_im)
AICc(ppm_A_precipitation)
# [1] 17026.54


gam_A_precipitation<-ppm(TypeA_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gam_A_precipitation)
# [1] 16938.72




# 11) modern annual mean temperature: GAM 
ppm_A_temperature<-ppm(TypeA_ppp~temperature_im)
AICc(ppm_A_temperature)
# [1] 17227.72


gam_A_temperature<-ppm(TypeA_ppp~s(temperature_im),use.gam=TRUE)
AICc(gam_A_temperature)
# [1] 17024.94



# 3.2 earth-wall sites
# 1) elevation: LM 
ppm_B_elevation<-ppm(TypeB_ppp~elevation_im)
AICc(ppm_B_elevation)
# [1] 2624.926
gam_B_elevation<-ppm(TypeB_ppp~s(elevation_im),use.gam=TRUE)
AICc(gam_B_elevation)
# [1] 2644.667

# 2) slope: LM 
ppm_B_slope<-ppm(TypeB_ppp~slope_im)
AICc(ppm_B_slope)
# [1] 2633.698
gam_B_slope<-ppm(TypeB_ppp~s(slope_im),use.gam=TRUE)
AICc(gam_B_slope)
# [1] 2653.439



# 3) aspect: LM 
ppm_B_aspect<-ppm(TypeB_ppp~aspect_im)
AICc(ppm_B_aspect)
# [1] 2634.049
gam_B_aspect<-ppm(TypeB_ppp~s(aspect_im),use.gam=TRUE)
AICc(gam_B_aspect)
# [1] 2653.79

# 4) incoming visibility index: LM
ppm_B_vis_in<-ppm(TypeB_ppp~vis_in_im)
AICc(ppm_B_vis_in)
# [1] 2618.28
gam_B_vis_in<-ppm(TypeB_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gam_B_vis_in)
# [1] 2621.042

# 5) outgoing visibility index: LM
ppm_B_vis_out<-ppm(TypeB_ppp~vis_out_im)
AICc(ppm_B_vis_out)
# [1] 2617.355
gam_B_vis_out<-ppm(TypeB_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gam_B_vis_out)
# [1] 2637.096

# 6) landforms: LM
ppm_B_landforms<-ppm(TypeB_ppp~landforms_recla_im)
AICc(ppm_B_landforms)
# [1] 2626.582


# 7) river: LM 
ppm_B_river<-ppm(TypeB_ppp~river_im)
AICc(ppm_B_river)
# [1] 2624.92
gam_B_river<-ppm(TypeB_ppp~s(river_im),use.gam=TRUE)
AICc(gam_B_river)
# [1] 2644.661


# 8) soil: LM
ppm_B_soil<-ppm(TypeB_ppp~soil_recla_im)
AICc(ppm_B_soil)
# [1] 2577.265

# 9) geology
# building earth: LM
ppm_B_geo_earth <- ppm(TypeB_ppp~geo_earth_im)
AICc(ppm_B_geo_earth)
# [1] 2632.412


gam_B_geo_earth <- ppm(TypeB_ppp~s(geo_earth_im), use.gam=TRUE)
AICc(gam_B_geo_earth)
# [1] 2652.153


# 10) modern annual precipitation: LM 
ppm_B_precipitation<-ppm(TypeB_ppp~precipitation_im)
AICc(ppm_B_precipitation)
# [1] 2581.457


gam_B_precipitation<-ppm(TypeB_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gam_B_precipitation)
# [1] 2601.198




# 11) modern annual mean temperature: LM 
ppm_B_temperature<-ppm(TypeB_ppp~temperature_im)
AICc(ppm_B_temperature)
# [1] 2627.685


gam_B_temperature<-ppm(TypeB_ppp~s(temperature_im),use.gam=TRUE)
AICc(gam_B_temperature)
# [1] 2647.426







# 4.Relative Variable Importance  
# 4.1 stone-wall sites
global.model.A<-ppm(TypeA_ppp~s(elevation_im)+s(slope_im)+aspect_im+s(vis_in_im)+s(vis_out_im)+landforms_recla_im+s(river_im)+soil_recla_im+geo_rock_im+s(precipitation_im)+s(temperature_im), use.gam=TRUE)
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
# params.A    aic.imp
# 1            aspect_im 0.21417858
# 2          geo_rock_im 0.93530647
# 3   landforms_recla_im 0.12108789
# 4      s(elevation_im) 1.00000000
# 5  s(precipitation_im) 1.00000000
# 6          s(river_im) 0.96390361
# 7          s(slope_im) 0.01009899
# 8    s(temperature_im) 0.51114889
# 9         s(vis_in_im) 0.16688331
# 10       s(vis_out_im) 0.87581242
# 11       soil_recla_im 0.27897210

sum(param.importance.A$aic.imp)

sum(param.comb.A$weights)

head(param.comb.A)
elevation.results=subset(param.comb.A, param.comb.A[,4]==1)
sum(elevation.results$weights) # 1
nrow(elevation.results)



# 4.2 earth-wall sites
global.model.B<-ppm(TypeB_ppp~elevation_im+slope_im+aspect_im+vis_in_im+vis_out_im+landforms_recla_im+river_im+soil_recla_im+geo_earth_im+precipitation_im+temperature_im, use.gam=TRUE)
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
#              params.B   aic.imp
# 1           aspect_im 0.2218102
# 2        elevation_im 0.2314323
# 3        geo_earth_im 0.9433698
# 4  landforms_recla_im 0.1582640
# 5    precipitation_im 1.0000000
# 6            river_im 0.2405533
# 7            slope_im 0.9225218
# 8       soil_recla_im 1.0000000
# 9      temperature_im 0.2272226
# 10          vis_in_im 0.6640273
# 11         vis_out_im 0.4541445

# 4.3 comparison between stone-wall sites and earth-wall sites
# code for Figure 5
importance.A.df <- data.frame(covariates=c("Aspect","Elevation", "Geology", "Landforms","Precipitation","River","Slope","Soil","Temperature","Visibility_in","Visibility_out"),importance.A=c(0.21417858, 1.00000000, 0.93530647, 0.12108789, 1.00000000, 0.96390361, 0.01009899, 0.27897210, 0.51114889, 0.16688331, 0.87581242))
importance.B.df <- data.frame(covariates=c("Aspect","Elevation", "Geology", "Landforms","Precipitation","River","Slope","Soil","Temperature","Visibility_in","Visibility_out"),importance.B=c(0.2218102, 0.2314323, 0.9433698, 0.1582640, 1.0000000, 0.2405533, 0.9225218, 1.0000000, 0.2272226, 0.6640273, 0.4541445))
importance.A.B.df <- data.frame(covariates=factor(importance.A.df$covariates), importance.A=as.numeric(importance.A.df$importance.A), importance.B=as.numeric(importance.B.df$importance.B))
tiff(here("Output", "2_first_order_interactions", "AICc_importance.tiff"),width = 12, height = 7, units = "in", res = 1000)
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
model.A<-ppm(TypeA_ppp~geo_rock_im+s(elevation_im) +s(precipitation_im) +s(river_im) +s(vis_out_im),use.gam=TRUE)
AICc(model.A)
# [1] 16590.23
pseudoR2(model.A)
# [1] 0.3967027

# code for Figure 6
model.A.Kres.env <-envelope(model.A,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "2_first_order_interactions", "model_A_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.A.Kres.env, main = NULL)
dev.off()
global.model.A.Kres.env <-envelope(global.model.A,fun=Kres,nsim=50,correction='best')
jpeg(here("Output", "2_first_order_interactions", "global_model_A_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(global.model.A.Kres.env,main=NULL)
dev.off()


# code for Figure 7
jpeg(here("Output", "2_first_order_interactions", "effectfun_A_geo_rock.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A,covname='geo_rock_im',elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A_elevation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A,covname='elevation_im',geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A,covname='precipitation_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()


jpeg(here("Output", "2_first_order_interactions", "effectfun_A_river.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A,covname='river_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_A_vis_out.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A,covname='vis_out_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), river_im=mean(river_im),se.fit=TRUE),main=NULL)
dev.off()

# 5.2 earth-wall sites
model.B<-ppm(TypeB_ppp~geo_earth_im + precipitation_im+ slope_im+ soil_recla_im)
pseudoR2(model.B)
# [1] 0.1836129

install.packages("ppcor")
library(ppcor)
B_geo_earth <- geo_earth_im[TypeB_ppp]
B_precipitation <- precipitation_im[TypeB_ppp]
B_slope <-slope_im[TypeB_ppp] 
B_soil <- soil_recla_im[TypeB_ppp]
B_soil <- extract(soil_resample,earth.sf)
B_soil_transformed <- ifelse(is.na(B_soil), 0, ifelse(B_soil > 0, 1, 0))
B_cor.df <- data.frame(matrix(nrow=67, ncol=4))
B_cor.df$geo_earth <- B_geo_earth
B_cor.df$precipitation <-B_precipitation 
B_cor.df$slope <- B_slope
B_cor.df$soil <- B_soil_transformed
# geology & precipitation: significant
pcor.test(B_cor.df$geo_earth, B_cor.df$precipitation, c(B_cor.df$slope,B_cor.df$soil), method="kendall")
#    estimate     p.value statistic   n gp  Method
# 1 0.4788054 2.95602e-16  8.175108 134  1 kendall

# geology & slope: significant
pcor.test(B_cor.df$geo_earth, B_cor.df$slope, c(B_cor.df$precipitation,B_cor.df$soil), method="kendall")
#    estimate     p.value statistic   n gp  Method
# 1 0.1699798 0.003705177  2.902229 134  1 kendall

# geology & soil: not significant
pcor.test(B_cor.df$geo_earth, B_cor.df$soil, c(B_cor.df$slope,B_cor.df$precipitation), method="kendall")
#     estimate   p.value statistic   n gp  Method
# 1 0.03823126 0.5139115 0.6527592 134  1 kendall

# precipitation & slope: significant
pcor.test(B_cor.df$precipitation, B_cor.df$slope,c(B_cor.df$geo_earth, B_cor.df$soil), method="kendall")
#    estimate    p.value statistic   n gp  Method
# 1 0.1193576 0.04155912  2.037908 134  1 kendall

# precipitation & soil: significant
pcor.test( B_cor.df$precipitation, B_cor.df$soil, c(B_cor.df$geo_earth,B_cor.df$slope), method="kendall")
#    estimate      p.value statistic   n gp  Method
# 1 0.2424726 3.473541e-05  4.139968 134  1 kendall

# slope & soil: not significant
pcor.test(B_cor.df$slope,B_cor.df$soil, c(B_cor.df$geo_earth, B_cor.df$precipitation), method="kendall")
#    estimate   p.value statistic   n gp  Method
# 1 0.0149666 0.7983067 0.2555391 134  1 kendall

# code for Figure 6
model.B.Kres.env <-envelope(model.B,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "2_first_order_interactions", "model_B_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(model.B.Kres.env, main = NULL)
dev.off()
global.model.B.Kres.env <-envelope(global.model.B,fun=Kres,nsim=500,correction='best')
jpeg(here("Output", "2_first_order_interactions", "global_model_B_Kres.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(global.model.B.Kres.env,main=NULL)
dev.off()




# code for Figure 7
jpeg(here("Output", "2_first_order_interactions", "effectfun_B_precipitation.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B,covname='precipitation_im',geo_earth_im=mean(geo_earth_im),slope_im=mean(slope_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_B_geo_earth.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B,covname='geo_earth_im',precipitation_im=mean(precipitation_im),slope_im=mean(slope_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_B_slope.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B,covname='slope_im',precipitation_im=mean(precipitation_im),geo_earth_im=mean(geo_earth_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "2_first_order_interactions", "effectfun_B_soil.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff_soil <- effectfun(model.B,covname='soil_recla_im',geo_earth_im=mean(geo_earth_im),slope_im=mean(slope_im),precipitation_im=mean(precipitation_im), se.fit=TRUE)
barplot(eff_soil$lambda,names.arg=eff_soil$soil_recla_im,main=NULL)
dev.off()


# 6.second-order interaction
# 6.1 stone-wall sites
distance_A_mean <- mean(st_distance(stone.sf))
distance_A_mean
# 33942.92 [m]
rvals <- data.frame(r=seq(50, 10000, by=50))
distance_A <-profilepl(rvals, Strauss, TypeA_ppp~geo_rock_im+s(elevation_im, k=3) +s(precipitation_im, k=3) +s(river_im, k=3) +s(vis_out_im, k=3), use.gam=TRUE) 
print(distance_A)
# profile log pseudolikelihood for model:
# ppm(TypeA_ppp ~ geo_rock_im + s(elevation_im,  k = 3) + s(precipitation_im,     
#   k = 3) + s(river_im,  k = 3) + s(vis_out_im,  k = 3),  use.gam = TRUE,       
# interaction = Strauss)
# fitted with rbord = 10000
# interaction: Strauss process
# irregular parameter: r in [50, 10000]
# optimum value of irregular parameter:  r = 850
jpeg(here("Output", "4_second_order_interactions", "Gibbs_A_distance.jpg"), width = 7, height = 7, units = "in", res = 1000)
plot(distance_A)
dev.off()

Gibbs_A <- as.ppm(distance_A)
AICc(Gibbs_A)
# [1] 7839.674
pseudoR2(Gibbs_A)
# [1] 0.3526079

# kres
Gibbs_A.Kres.env <-envelope(Gibbs_A,fun=Kres,nsim=50,correction='best')
# Error: This calculation is not supported for GAM fits
Gibbs_A.Kres <- Kres(Gibbs_A,correction='best')
plot(Gibbs_A.Kres)

jpeg(here("Output", "4_second_order_interactions", "Gibbs_A_Kres_noenv.jpg",width = 7, height = 7, units = "in", res = 1000)
plot(Gibbs_A.Kres, main = NULL)
dev.off()

# effectfun
jpeg(here("Output", "4_second_order_interactions", "Deffectfun_Gibbs_A_geo_rock.jpg"), width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(Gibbs_A,covname='geo_rock_im',elevation_im=mean(elevation_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_elevation.jpg"), width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(Gibbs_A,covname='elevation_im',geo_rock_im=mean(geo_rock_im),precipitation_im = mean(precipitation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_precipitation.jpg"), width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(Gibbs_A,covname='precipitation_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), river_im  =mean(river_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()


jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_river.jpg"), width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(Gibbs_A,covname='river_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "4_second_order_interactions", "effectfun_Gibbs_A_vis_out.jpg"), width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(Gibbs_A,covname='vis_out_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), river_im=mean(river_im),se.fit=TRUE),main=NULL)
dev.off()





# 6.2 interaction between types
AB_kcross <- envelope(AB.multitype,fun=Kcross,nsim=50,correction='best')
jpeg(here("Output", "4_second_order_interactions", "AB_kcross.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(AB_kcross,main=NULL)
dev.off()
# 6.2.1 function form decision 
# 1) elevation: GAM 
ppmAB_elevation<-ppm(AB.multitype~elevation_im)
AICc(ppmAB_elevation)
# [1] 20096.94
gamAB_elevation<-ppm(AB.multitype~s(elevation_im),use.gam=TRUE)
AICc(gamAB_elevation)
# [1] 19945


# 2) slope: GAM   
ppmAB_slope<-ppm(AB.multitype~slope_im)
AICc(ppmAB_slope)
# [1] 20409.18
gamAB_slope<-ppm(AB.multitype~s(slope_im),use.gam=TRUE)
AICc(gamAB_slope)
#[1] 20354.17

# 3) aspect:  LM 
ppmAB_aspect<-ppm(AB.multitype~aspect_im)
AICc(ppmAB_aspect)
# [1] 20519.56
gamAB_aspect<-ppm(AB.multitype~s(aspect_im),use.gam=TRUE)
AICc(gamAB_aspect)
# [1] 20535.94

# 4) incoming visibility index: GAM   
ppmAB_vis_in<-ppm(AB.multitype~vis_in_im)
AICc(ppmAB_vis_in)
#[1] 20144.32 
gamAB_vis_in<-ppm(AB.multitype~s(vis_in_im),use.gam=TRUE)
AICc(gamAB_vis_in)
#[1] 20022.71 

# 5) outgoing visibility index: GAM 
ppmAB_vis_out<-ppm(AB.multitype~vis_out_im)
AICc(ppmAB_vis_out)
# [1] 20102.46
gamAB_vis_out<-ppm(AB.multitype~s(vis_out_im),use.gam=TRUE)
AICc(gamAB_vis_out)
# [1] 20012.32

# 6) landforms:  LM
ppmAB_landforms<-ppm(AB.multitype~landforms_recla_im)
AICc(ppmAB_landforms)
# [1] 20411.72


# 7) river: GAM   
ppmAB_river<-ppm(AB.multitype~river_im)
AICc(ppmAB_river)
# [1] 20498.14
gamAB_river<-ppm(AB.multitype~s(river_im),use.gam=TRUE)
AICc(gamAB_river)
# [1] 20487.04


# 8) soil:  LM 
ppmAB_soil<-ppm(AB.multitype~soil_recla_im)
AICc(ppmAB_soil)
#[1] 20261.6 


# 9) geology
# building earth: LM   
ppmAB_geo_rock <- ppm(AB.multitype~geo_rock_im)
AICc(ppmAB_geo_rock)
# [1] 20152.58
gamAB_geo_rock <- ppm(AB.multitype~s(geo_rock_im), use.gam=TRUE)
AICc(gamAB_geo_rock)
# [1] 20168.96


# building earth: LM   
ppmAB_geo_earth <- ppm(AB.multitype~geo_earth_im)
AICc(ppmAB_geo_earth)
# [1] 20480.31
gamAB_geo_earth <- ppm(AB.multitype~s(geo_earth_im), use.gam=TRUE)
AICc(gamAB_geo_earth)
# [1] 20496.68


# 10) modern annual precipitation: LM
ppmAB_precipitation<-ppm(AB.multitype~precipitation_im)
AICc(ppmAB_precipitation)
# [1] 19975.48


gamAB_precipitation<-ppm(AB.multitype~s(precipitation_im),use.gam=TRUE)
AICc(gamAB_precipitation)
# [1] 19991.86




# 11) modern annual mean temperature: GAM    
ppmAB_temperature<-ppm(AB.multitype~temperature_im)
AICc(ppmAB_temperature)
# [1] 20231.48


gamAB_temperature<-ppm(AB.multitype~s(temperature_im),use.gam=TRUE)
AICc(gamAB_temperature)
#[1] 20007.1




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



# 7.site size and model
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

# 7.2 test the model performance using the selected covariates above through relative variable improtance analysis
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


jpeg(here("Output", "3_size_categories", "effectfun_A_large_river.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_large,covname='river_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_large_vis_out.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_large,covname='vis_out_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), river_im=mean(river_im),se.fit=TRUE),main=NULL)
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


jpeg(here("Output", "3_size_categories", "effectfun_A_middle_river.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_middle,covname='river_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_middle_vis_out.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_middle,covname='vis_out_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), river_im=mean(river_im),se.fit=TRUE),main=NULL)
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


jpeg(here("Output", "3_size_categories", "effectfun_A_small_river.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_small,covname='river_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), vis_out_im=mean(vis_out_im),se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_A_small_vis_out.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.A_small,covname='vis_out_im',geo_rock_im=mean(geo_rock_im),elevation_im = mean(elevation_im), precipitation_im =mean(precipitation_im), river_im=mean(river_im),se.fit=TRUE),main=NULL)
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

jpeg(here("Output", "3_size_categories", "effectfun_B_large_slope.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_large,covname='slope_im',precipitation_im=mean(precipitation_im),geo_earth_im=mean(geo_earth_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_B_large_soil.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff_soil <- effectfun(model.B_large,covname='soil_recla_im',geo_earth_im=mean(geo_earth_im),slope_im=mean(slope_im),precipitation_im=mean(precipitation_im), se.fit=TRUE)
barplot(eff_soil$lambda,names.arg=eff_soil$soil_recla_im,main=NULL)
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

jpeg(here("Output", "3_size_categories", "effectfun_B_middle_slope.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_middle,covname='slope_im',precipitation_im=mean(precipitation_im),geo_earth_im=mean(geo_earth_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_B_middle_soil.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff_soil <- effectfun(model.B_middle,covname='soil_recla_im',geo_earth_im=mean(geo_earth_im),slope_im=mean(slope_im),precipitation_im=mean(precipitation_im), se.fit=TRUE)
barplot(eff_soil$lambda,names.arg=eff_soil$soil_recla_im,main=NULL)
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

jpeg(here("Output", "3_size_categories", "effectfun_B_small_slope.jpg"),width = 7, height = 7, units = "in", res = 1000)
plot(effectfun(model.B_small,covname='slope_im',precipitation_im=mean(precipitation_im),geo_earth_im=mean(geo_earth_im),soil_recla_im=as.factor(c("absence", "presence")), se.fit=TRUE),main=NULL)
dev.off()

jpeg(here("Output", "3_size_categories", "effectfun_B_small_soil.jpg"),width = 7, height = 7, units = "in", res = 1000)
eff_soil <- effectfun(model.B_small,covname='soil_recla_im',geo_earth_im=mean(geo_earth_im),slope_im=mean(slope_im),precipitation_im=mean(precipitation_im), se.fit=TRUE)
barplot(eff_soil$lambda,names.arg=eff_soil$soil_recla_im,main=NULL)
dev.off()



# The following code will be rerunning the whole workflow to compare among different size categories as suggested by reviewer 2.
# 7.3 function form
# 7.3.1 stone-wall sites
# 7.3.1.1 large
# 1) elevation: GAM 
ppmA_large_elevation<-ppm(A_large_ppp~elevation_im)
AICc(ppmA_large_elevation)
#[1] 5837.845 
gamA_large_elevation<-ppm(A_large_ppp~s(elevation_im),use.gam=TRUE)
AICc(gamA_large_elevation)
# [1] 5803.299

# 2) slope: GAM  
ppmA_large_slope<-ppm(A_large_ppp~slope_im)
AICc(ppmA_large_slope)
#[1] 5914.925 
gamA_large_slope<-ppm(A_large_ppp~s(slope_im),use.gam=TRUE)
AICc(gamA_large_slope)
# [1] 5913.509


# 3) aspect: LM  
ppmA_large_aspect<-ppm(A_large_ppp~aspect_im)
AICc(ppmA_large_aspect)
#[1] 5934.773 
gamA_large_aspect<-ppm(A_large_ppp~s(aspect_im),use.gam=TRUE)
AICc(gamA_large_aspect)
# [1] 5952.192

# 4) incoming visibility index: GAM 
ppmA_large_vis_in<-ppm(A_large_ppp~vis_in_im)
AICc(ppmA_large_vis_in)
# [1] 5828.997
gamA_large_vis_in<-ppm(A_large_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gamA_large_vis_in)
# 5806.729

# 5) outgoing visibility index: GAM 
ppmA_large_vis_out<-ppm(A_large_ppp~vis_out_im)
AICc(ppmA_large_vis_out)
# [1] 5820.753
gamA_large_vis_out<-ppm(A_large_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gamA_large_vis_out)
# [1] 5810.293

# 6) landforms: LM 
ppmA_large_landforms<-ppm(A_large_ppp~landforms_recla_im)
AICc(ppmA_large_landforms)
# [1] 5897.547


# 7) river: LM  
ppmA_large_river<-ppm(A_large_ppp~river_im)
AICc(ppmA_large_river)
# [1] 5933.795
gamA_large_river<-ppm(A_large_ppp~s(river_im),use.gam=TRUE)
AICc(gamA_large_river)
# [1] 5945.167


# 8) soil: LM 
ppmA_large_soil<-ppm(A_large_ppp~soil_recla_im)
AICc(ppmA_large_soil)
# [1] 5870.646


# 9) geology
# building rock:  LM
ppmA_large_geo_rock <- ppm(A_large_ppp~geo_rock_im)
AICc(ppmA_large_geo_rock)
# [1] 5830.514
gamA_large_geo_rock <- ppm(A_large_ppp~s(geo_rock_im), use.gam=TRUE)
AICc(gamA_large_geo_rock)
# [1] 5847.933



# 10) modern annual precipitation: LM  
ppmA_large_precipitation<-ppm(A_large_ppp~precipitation_im)
AICc(ppmA_large_precipitation)
# [1] 5840.826


gamA_large_precipitation<-ppm(A_large_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gamA_large_precipitation)
# [1] 5858.246




# 11) modern annual mean temperature: GAM  
ppmA_large_temperature<-ppm(A_large_ppp~temperature_im)
AICc(ppmA_large_temperature)
#[1] 5884.647 


gamA_large_temperature<-ppm(A_large_ppp~s(temperature_im),use.gam=TRUE)
AICc(gamA_large_temperature)
#[1] 5828.167 



# 7.3.1.2 middle
# 1) elevation: GAM
ppmA_middle_elevation<-ppm(A_middle_ppp~elevation_im)
AICc(ppmA_middle_elevation)
#[1] 5764.782 
gamA_middle_elevation<-ppm(A_middle_ppp~s(elevation_im),use.gam=TRUE)
AICc(gamA_middle_elevation)
#[1] 5731.321 

# 2) slope:  GAM 
ppmA_middle_slope<-ppm(A_middle_ppp~slope_im)
AICc(ppmA_middle_slope)
#[1] 5897.752 
gamA_middle_slope<-ppm(A_middle_ppp~s(slope_im),use.gam=TRUE)
AICc(gamA_middle_slope)
# [1] 5889.973


# 3) aspect: LM  
ppmA_middle_aspect<-ppm(A_middle_ppp~aspect_im)
AICc(ppmA_middle_aspect)
# [1] 5934.848
gamA_middle_aspect<-ppm(A_middle_ppp~s(aspect_im),use.gam=TRUE)
AICc(gamA_middle_aspect)
# [1] 5952.267

# 4) incoming visibility index: GAM 
ppmA_middle_vis_in<-ppm(A_middle_ppp~vis_in_im)
AICc(ppmA_middle_vis_in)
# [1] 5824.141
gamA_middle_vis_in<-ppm(A_middle_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gamA_middle_vis_in)
# [1] 5811.493

# 5) outgoing visibility index: LM 
ppmA_middle_vis_out<-ppm(A_middle_ppp~vis_out_im)
AICc(ppmA_middle_vis_out)
# [1] 5819.025
gamA_middle_vis_out<-ppm(A_middle_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gamA_middle_vis_out)
# [1] 5819.187

# 6) landforms: LM 
ppmA_middle_landforms<-ppm(A_middle_ppp~landforms_recla_im)
AICc(ppmA_middle_landforms)
#[1] 5919.562 


# 7) river: LM  
ppmA_middle_river<-ppm(A_middle_ppp~river_im)
AICc(ppmA_middle_river)
# [1] 5933.332
gamA_middle_river<-ppm(A_middle_ppp~s(river_im),use.gam=TRUE)
AICc(gamA_middle_river)
# [1] 5942.639


# 8) soil:  LM
ppmA_middle_soil<-ppm(A_middle_ppp~soil_recla_im)
AICc(ppmA_middle_soil)
# [1] 5856.275


# 9) geology
# building rock:  LM
ppmA_large_geo_rock <- ppm(A_large_ppp~geo_rock_im)
AICc(ppmA_large_geo_rock)
# [1] 5830.514
gamA_large_geo_rock <- ppm(A_large_ppp~s(geo_rock_im), use.gam=TRUE)
AICc(gamA_large_geo_rock)
# [1] 5847.933


# 10) modern annual precipitation:  GAM 
ppmA_middle_precipitation<-ppm(A_middle_ppp~precipitation_im)
AICc(ppmA_middle_precipitation)
# [1] 5690.922


gamA_middle_precipitation<-ppm(A_middle_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gamA_middle_precipitation)
# [1] 5666.542




# 11) modern annual mean temperature:  GAM  
ppmA_middle_temperature<-ppm(A_middle_ppp~temperature_im)
AICc(ppmA_middle_temperature)
# [1] 5803.593



gamA_middle_temperature<-ppm(A_middle_ppp~s(temperature_im),use.gam=TRUE)
AICc(gamA_middle_temperature)
#[1] 5757.494


# 7.3.1.3 small
# 1) elevation:  GAM
ppmA_small_elevation<-ppm(A_small_ppp~elevation_im)
AICc(ppmA_small_elevation)
# [1] 5799.311
gamA_small_elevation<-ppm(A_small_ppp~s(elevation_im),use.gam=TRUE)
AICc(gamA_small_elevation)
# [1] 5778.092

# 2) slope: GAM  
ppmA_small_slope<-ppm(A_small_ppp~slope_im)
AICc(ppmA_small_slope)
# [1] 5926.461
gamA_small_slope<-ppm(A_small_ppp~s(slope_im),use.gam=TRUE)
AICc(gamA_small_slope)
# [1] 5917.161


# 3) aspect: LM  
ppmA_small_aspect<-ppm(A_small_ppp~aspect_im)
AICc(ppmA_small_aspect)
# [1] 5969.767
gamA_small_aspect<-ppm(A_small_ppp~s(aspect_im),use.gam=TRUE)
AICc(gamA_small_aspect)
# [1] 5987.177

# 4) incoming visibility index: GAM 
ppmA_small_vis_in<-ppm(A_small_ppp~vis_in_im)
AICc(ppmA_small_vis_in)
# [1] 5854.403
gamA_small_vis_in<-ppm(A_small_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gamA_small_vis_in)
# [1] 5843.908

# 5) outgoing visibility index: GAM 
ppmA_small_vis_out<-ppm(A_small_ppp~vis_out_im)
AICc(ppmA_small_vis_out)
# [1] 5842.809
gamA_small_vis_out<-ppm(A_small_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gamA_small_vis_out)
# [1] 5836.268

# 6) landforms: LM 
ppmA_small_landforms<-ppm(A_small_ppp~landforms_recla_im)
AICc(ppmA_small_landforms)
# [1] 5933.57


# 7) river:  GAM 
ppmA_small_river<-ppm(A_small_ppp~river_im)
AICc(ppmA_small_river)
# [1] 5954.971
gamA_small_river<-ppm(A_small_ppp~s(river_im),use.gam=TRUE)
AICc(gamA_small_river)
# [1] 5952.314


# 8) soil: LM 
ppmA_small_soil<-ppm(A_small_ppp~soil_recla_im)
AICc(ppmA_small_soil)
# [1] 5898.741


# 9) geology
# building rock: LM 
ppmA_small_geo_rock <- ppm(A_small_ppp~geo_rock_im)
AICc(ppmA_small_geo_rock)
# [1] 5873.034
gamA_small_geo_rock <- ppm(A_small_ppp~s(geo_rock_im), use.gam=TRUE)
AICc(gamA_small_geo_rock)
# [1] 5890.444



# 10) modern annual precipitation: LM  
ppmA_small_precipitation<-ppm(A_small_ppp~precipitation_im)
AICc(ppmA_small_precipitation)
# [1] 5791.129


gamA_small_precipitation<-ppm(A_small_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gamA_small_precipitation)
# [1] 5808.539




# 11) modern annual mean temperature: GAM  
ppmA_small_temperature<-ppm(A_small_ppp~temperature_im)
AICc(ppmA_small_temperature)
# [1] 5848.052


gamA_small_temperature<-ppm(A_small_ppp~s(temperature_im),use.gam=TRUE)
AICc(gamA_small_temperature)
#[1] 5801.619




# 7.3.2 earth-wall sites
# 7.3.2.1 large
# 1) elevation: LM  
ppmB_large_elevation<-ppm(B_large_ppp~elevation_im)
AICc(ppmB_large_elevation)
# [1] 915.5857
gamB_large_elevation<-ppm(B_large_ppp~s(elevation_im),use.gam=TRUE)
AICc(gamB_large_elevation)
# [1] 950.9541

# 2) slope:  LM 
ppmB_large_slope<-ppm(B_large_ppp~slope_im)
AICc(ppmB_large_slope)
# [1] 917.1239
gamB_large_slope<-ppm(B_large_ppp~s(slope_im),use.gam=TRUE)
AICc(gamB_large_slope)
# [1] 952.4923


# 3) aspect:  LM 
ppmB_large_aspect<-ppm(B_large_ppp~aspect_im)
AICc(ppmB_large_aspect)
# [1] 916.8271
gamB_large_aspect<-ppm(B_large_ppp~s(aspect_im),use.gam=TRUE)
AICc(gamB_large_aspect)
# [1] 952.1955

# 4) incoming visibility index:  LM 
ppmB_large_vis_in<-ppm(B_large_ppp~vis_in_im)
AICc(ppmB_large_vis_in)
# [1] 915.1926
gamB_large_vis_in<-ppm(B_large_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gamB_large_vis_in)
# [1] 950.561

# 5) outgoing visibility index: LM 
ppmB_large_vis_out<-ppm(B_large_ppp~vis_out_im)
AICc(ppmB_large_vis_out)
# [1] 914.022
gamB_large_vis_out<-ppm(B_large_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gamB_large_vis_out)
# [1] 949.3904

# 6) landforms:  LM
ppmB_large_landforms<-ppm(B_large_ppp~landforms_recla_im)
AICc(ppmB_large_landforms)
# [1] 921.4442


# 7) river:  LM 
ppmB_large_river<-ppm(B_large_ppp~river_im)
AICc(ppmB_large_river)
# [1] 916.3852
gamB_large_river<-ppm(B_large_ppp~s(river_im),use.gam=TRUE)
AICc(gamB_large_river)
# [1] 951.7536


# 8) soil:  LM
ppmB_large_soil<-ppm(B_large_ppp~soil_recla_im)
AICc(ppmB_large_soil)
# [1] 874.5221


# 9) geology
# building earth:  LM
ppmB_large_geo_earth <- ppm(B_large_ppp~geo_earth_im)
AICc(ppmB_large_geo_earth)
# [1] 916.5294


gamB_large_geo_earth <- ppm(B_large_ppp~s(geo_earth_im), use.gam=TRUE)
AICc(gamB_large_geo_earth)
# [1] 951.8978


# 10) modern annual precipitation:  LM  
ppmB_large_precipitation<-ppm(B_large_ppp~precipitation_im)
AICc(ppmB_large_precipitation)
# [1] 905.3853


gamB_large_precipitation<-ppm(B_large_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gamB_large_precipitation)
# [1] 940.7537




# 11) modern annual mean temperature: LM  
ppmB_large_temperature<-ppm(B_large_ppp~temperature_im)
AICc(ppmB_large_temperature)
# [1] 916.8718


gamB_large_temperature<-ppm(B_large_ppp~s(temperature_im),use.gam=TRUE)
AICc(gamB_large_temperature)
# [1] 952.2402




# 7.3.2.2 middle
# 1) elevation: LM
ppmB_middle_elevation<-ppm(B_middle_ppp~elevation_im)
AICc(ppmB_middle_elevation)
# [1] 914.176
gamB_middle_elevation<-ppm(B_middle_ppp~s(elevation_im),use.gam=TRUE)
AICc(gamB_middle_elevation)
# [1] 949.5444

# 2) slope: LM  
ppmB_middle_slope<-ppm(B_middle_ppp~slope_im)
AICc(ppmB_middle_slope)
# [1] 916.872
gamB_middle_slope<-ppm(B_middle_ppp~s(slope_im),use.gam=TRUE)
AICc(gamB_middle_slope)
# [1] 952.2405


# 3) aspect:  LM 
ppmB_middle_aspect<-ppm(B_middle_ppp~aspect_im)
AICc(ppmB_middle_aspect)
# [1] 914.7698
gamB_middle_aspect<-ppm(B_middle_ppp~s(aspect_im),use.gam=TRUE)
AICc(gamB_middle_aspect)
# [1] 950.1382

# 4) incoming visibility index:  LM 
ppmB_middle_vis_in<-ppm(B_middle_ppp~vis_in_im)
AICc(ppmB_middle_vis_in)
# [1] 909.5991
gamB_middle_vis_in<-ppm(B_middle_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gamB_middle_vis_in)
# [1] 944.9676

# 5) outgoing visibility index:  LM 
ppmB_middle_vis_out<-ppm(B_middle_ppp~vis_out_im)
AICc(ppmB_middle_vis_out)
# [1] 909.8802
gamB_middle_vis_out<-ppm(B_middle_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gamB_middle_vis_out)
# [1] 945.2486

# 6) landforms: LM 
ppmB_middle_landforms<-ppm(B_middle_ppp~landforms_recla_im)
AICc(ppmB_middle_landforms)
# [1] 914.434


# 7) river:  LM 
ppmB_middle_river<-ppm(B_middle_ppp~river_im)
AICc(ppmB_middle_river)
# [1] 916.5239
gamB_middle_river<-ppm(B_middle_ppp~s(river_im),use.gam=TRUE)
AICc(gamB_middle_river)
# [1] 951.8923


# 8) soil: LM 
ppmB_middle_soil<-ppm(B_middle_ppp~soil_recla_im)
AICc(ppmB_middle_soil)
# [1] 910.7425


# 9) geology
# building earth: LM 
ppmB_middle_geo_earth <- ppm(B_middle_ppp~geo_earth_im)
AICc(ppmB_middle_geo_earth)
# [1] 915.5148


gamB_middle_geo_earth <- ppm(B_middle_ppp~s(geo_earth_im), use.gam=TRUE)
AICc(gamB_middle_geo_earth)
# [1] 950.8832


# 10) modern annual precipitation: LM  
ppmB_middle_precipitation<-ppm(B_middle_ppp~precipitation_im)
AICc(ppmB_middle_precipitation)
# [1] 898.5571


gamB_middle_precipitation<-ppm(B_middle_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gamB_middle_precipitation)
# [1] 933.9255




# 11) modern annual mean temperature: LM  
ppmB_middle_temperature<-ppm(B_middle_ppp~temperature_im)
AICc(ppmB_middle_temperature)
# [1] 913.881


gamB_middle_temperature<-ppm(B_middle_ppp~s(temperature_im),use.gam=TRUE)
AICc(gamB_middle_temperature)
#[1] 949.2494


# 7.3.2.3 small
# 1) elevation: LM
ppmB_small_elevation<-ppm(B_small_ppp~elevation_im)
AICc(ppmB_small_elevation)
# [1] 913.0893
gamB_small_elevation<-ppm(B_small_ppp~s(elevation_im),use.gam=TRUE)
AICc(gamB_small_elevation)
# [1] 948.4577

# 2) slope: LM  
ppmB_small_slope<-ppm(B_small_ppp~slope_im)
AICc(ppmB_small_slope)
# [1] 916.446
gamB_small_slope<-ppm(B_small_ppp~s(slope_im),use.gam=TRUE)
AICc(gamB_small_slope)
# [1] 951.8144


# 3) aspect:  LM 
ppmB_small_aspect<-ppm(B_small_ppp~aspect_im)
AICc(ppmB_small_aspect)
# [1] 916.7937
gamB_small_aspect<-ppm(B_small_ppp~s(aspect_im),use.gam=TRUE)
AICc(gamB_small_aspect)
# [1] 952.1621

# 4) incoming visibility index:  LM 
ppmB_small_vis_in<-ppm(B_small_ppp~vis_in_im)
AICc(ppmB_small_vis_in)
# [1] 909.4742
gamB_small_vis_in<-ppm(B_small_ppp~s(vis_in_im),use.gam=TRUE)
AICc(gamB_small_vis_in)
# [1] 944.8426

# 5) outgoing visibility index: LM 
ppmB_small_vis_out<-ppm(B_small_ppp~vis_out_im)
AICc(ppmB_small_vis_out)
# [1] 910.5368
gamB_small_vis_out<-ppm(B_small_ppp~s(vis_out_im),use.gam=TRUE)
AICc(gamB_small_vis_out)
# [1] 945.9053

# 6) landforms:  LM
ppmB_small_landforms<-ppm(B_small_ppp~landforms_recla_im)
AICc(ppmB_small_landforms)
# [1] 911.9101


# 7) river:  LM 
ppmB_small_river<-ppm(B_small_ppp~river_im)
AICc(ppmB_small_river)
# [1] 906.5912
gamB_small_river<-ppm(B_small_ppp~s(river_im),use.gam=TRUE)
AICc(gamB_small_river)
# [1] 941.9597


# 8) soil:  LM
ppmB_small_soil<-ppm(B_small_ppp~soil_recla_im)
AICc(ppmB_small_soil)
# [1] 904.8762


# 9) geology
# building earth:  LM 
ppmB_small_geo_earth <- ppm(B_small_ppp~geo_earth_im)
AICc(ppmB_small_geo_earth)
# [1] 916.8589


gamB_small_geo_earth <- ppm(B_small_ppp~s(geo_earth_im), use.gam=TRUE)
AICc(gamB_small_geo_earth)
# [1] 952.2273


# 10) modern annual precipitation: LM  
ppmB_small_precipitation<-ppm(B_small_ppp~precipitation_im)
AICc(ppmB_small_precipitation)
# [1] 897.1254


gamB_small_precipitation<-ppm(B_small_ppp~s(precipitation_im),use.gam=TRUE)
AICc(gamB_small_precipitation)
# [1] 932.4938




# 11) modern annual mean temperature: LM  
ppmB_small_temperature<-ppm(B_small_ppp~temperature_im)
AICc(ppmB_small_temperature)
# [1] 914.3507


gamB_small_temperature<-ppm(B_small_ppp~s(temperature_im),use.gam=TRUE)
AICc(gamB_small_temperature)
#[1] 949.7191




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
#        params.A_large    aic.imp
# 1           aspect_im 0.23680785
# 2         geo_rock_im 0.95494142
# 3  landforms_recla_im 0.09221989
# 4    precipitation_im 0.99539401
# 5            river_im 0.24461800
# 6     s(elevation_im) 0.15690398
# 7         s(slope_im) 0.00580755
# 8   s(temperature_im) 0.12075698
# 9        s(vis_in_im) 0.79633524
# 10      s(vis_out_im) 0.12015620
# 11      soil_recla_im 0.70681886
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
#        params.A_middle      aic.imp
# 1            aspect_im 0.3245973792
# 2          geo_rock_im 0.4261278384
# 3   landforms_recla_im 0.1258377513
# 4             river_im 0.4295288587
# 5      s(elevation_im) 0.9519855957
# 6  s(precipitation_im) 1.0000000000
# 7          s(slope_im) 0.0111033310
# 8    s(temperature_im) 0.0006001801
# 9         s(vis_in_im) 0.0000000000
# 10       soil_recla_im 0.1033309993
# 11          vis_out_im 0.9981994598
# 

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
#        params.A_small    aic.imp
# 1           aspect_im 0.21518607
# 2         geo_rock_im 0.24249700
# 3  landforms_recla_im 0.94937975
# 4    precipitation_im 1.00000000
# 5     s(elevation_im) 1.00000000
# 6         s(river_im) 1.00000000
# 7         s(slope_im) 0.00270108
# 8   s(temperature_im) 0.00010004
# 9        s(vis_in_im) 0.02200880
# 10      s(vis_out_im) 0.04631853
# 11      soil_recla_im 0.26830732
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
#        params.B_large     aic.imp
# 1           aspect_im 0.210299569
# 2        elevation_im 0.163911432
# 3        geo_earth_im 0.464282136
# 4  landforms_recla_im 0.003306282
# 5    precipitation_im 0.946398156
# 6            river_im 0.147680593
# 7            slope_im 0.455365194
# 8       soil_recla_im 1.000000000
# 9      temperature_im 0.181544935
# 10          vis_in_im 0.200080152
# 11         vis_out_im 0.243763150





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
#       params.B_middle    aic.imp
# 1           aspect_im 0.20485554
# 2        elevation_im 0.17275281
# 3        geo_earth_im 0.78611557
# 4  landforms_recla_im 0.03190209
# 5    precipitation_im 0.99979936
# 6            river_im 0.27608347
# 7            slope_im 0.38864366
# 8       soil_recla_im 0.13563403
# 9      temperature_im 0.14977929
# 10          vis_in_im 0.54263644
# 11         vis_out_im 0.36807785


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
#        params.B_small   aic.imp
# 1           aspect_im 0.1633803
# 2        elevation_im 0.1543260
# 3        geo_earth_im 0.1955734
# 4  landforms_recla_im 0.3110664
# 5    precipitation_im 0.9363179
# 6            river_im 0.5139839
# 7            slope_im 0.1985915
# 8       soil_recla_im 0.4316901
# 9      temperature_im 0.1632797
# 10          vis_in_im 0.3988934
# 11         vis_out_im 0.2529175
# 


# 7.5 final model and goodness of fit
# 7.5.1 stone-wall sites
# 1) large
model.A2_large<-ppm(A_large_ppp~geo_rock_im+precipitation_im,use.gam=TRUE)
pseudoR2(model.A2_large)
# [1] 0.188095

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
# [1] 0.4474158

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
# [1] 0.4400576
print(model.A2_small)
# Nonstationary Poisson process
# Fitted to point pattern dataset 'A_small_ppp'
# 
# Log intensity:  ~landforms_recla_im + precipitation_im +
# s(elevation_im) + s(river_im)
# 
# Fitted trend coefficients:
#                           (Intercept)
#                         -1.501280e+02
#    landforms_recla_imvalley/footslope 
#                          1.177073e+02
#   landforms_recla_imspur/slope/hollow
#                          1.183553e+02
# landforms_recla_impeak/ridge/shoulder
#                          1.190215e+02
#                      precipitation_im
#                          2.817669e-02
#                     s(elevation_im).1 
#                         -3.906354e+00
#                     s(elevation_im).2
#                          9.491719e+00
#                     s(elevation_im).3
#                         -1.644905e+00
#                     s(elevation_im).4
#                         -5.985818e+00
#                     s(elevation_im).5 
#                         -1.209152e+00
#                     s(elevation_im).6
#                         -3.522385e+00
#                     s(elevation_im).7
#                          4.176305e-01
#                     s(elevation_im).8
#                          1.710496e+01 
#                     s(elevation_im).9
#                          1.043936e+01
#                         s(river_im).1
#                          7.134894e-01
#                         s(river_im).2
#                          7.705912e-01 
#                         s(river_im).3
#                          4.713532e-02
#                         s(river_im).4
#                         -1.461531e-01
#                         s(river_im).5
#                          5.544825e-02
#                         s(river_im).6 
#                          3.953307e-01
#                         s(river_im).7
#                         -6.926657e-02
#                         s(river_im).8
#                         -1.125573e+00
#                         s(river_im).9
#                          9.525765e-03

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
# [1] 0.1018954

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
# [1] 0.1075001

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
# [1] 0.1154701

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











































