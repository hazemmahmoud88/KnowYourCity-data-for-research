##############################################################
### Comparison of estimates populations in Nigerian slums  ###
### By: Dana R Thomson (dana.r.thomson@gmail.com)          ###
### 07 Oct 2020                                            ###
##############################################################

# This project is funded by CIESIN, Columbia University
# Objectives:
#   1. Compare gridded population estimates to field population estiamtes in Lagos, Port 
#      Harcourt, and Nairobi where SDI enumerated slums (https://knowyourcity.info/explore-our-data/)
#   2. Compare gridded population estimates for all "slum" areas in Lagos


####################################    SET-UP    #######################################


######################
## Step 0. Initalise
######################

#rm(list=ls())

setwd("D:/Dropbox/Work_main/p_ciesin_pop_slum_compare/data")

proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
proj.utm31 <- "+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" #WGS84/UTM 31N # Lagos
proj.utm32 <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" #WGS84/UTM 32N # Port Harcourt
proj.utm37 <- "+proj=utm +zone=37 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" #WGS84/UTM 37S # Nairobi

library(raster)
library(rgeos)
library(rgdal)
library(dplyr)
library(sp)
library(exactextractr)
library(tmaptools)
library(sf)
library(ggplot2)
library(tmap)
library(data.table)
library(plyr)
library(haven)


####################################    ANALYSIS 1    #######################################

######################
## Step 1. Read & join pop data to boundaries
######################

##### LAGOS ######

# Lagos fishnet 200x200m 
#lag_grid <- readOGR("working/lag_fish2.shp")
#lag_grid <- st_as_sf(spTransform(lag_grid, CRS = proj.utm31))

# KnowYourCity Campaign slum boundaries
slums <- readOGR("kyc_campaign/kyc_nga.shp")
slums <- slums[slums$Id <= 35,] 

#create unprojected bounding box to crop rasters
bbox <- bb(slums, ext=1.2, output=c("extent"))

#create sf object for exact_extract calcs
slums <- st_as_sf(spTransform(slums, CRS = proj.utm31))

#GPWv4.11 (2020)
gpw20 <- raster("gpw/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")
gpw20 <- crop(gpw20, bbox)
gpw20 <- projectRaster(gpw20, crs=proj.utm31)
slums$gpw20 <- exact_extract(gpw20, slums, 'sum')
#lag_grid$gpw20 <- exact_extract(gpw20, lag_grid, 'sum')
rm(gpw20)

#GPWv4.11 (2015)
gpw15 <- raster("gpw/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec.tif")
gpw15 <- crop(gpw15, bbox)
gpw15 <- projectRaster(gpw15, crs=proj.utm31)
slums$gpw15 <- exact_extract(gpw15, slums, 'sum')
#lag_grid$gpw15 <- exact_extract(gpw15, lag_grid, 'sum')
rm(gpw15)

#GHS-POP (2015)
ghs15 <- raster("ghspop/GHS_POP_E2015_GLOBE_R2019A_54009_250_V1_0_18_8.tif")
ghs15 <- projectRaster(ghs15, crs=proj.utm31)
slums$ghs15 <- exact_extract(ghs15, slums, 'sum')
#lag_grid$ghs15 <- exact_extract(ghs15, lag_grid, 'sum')
rm(ghs15)

#HRSL-Facebook (2018)
hrsl18 <- raster("hrsl/population_nga_2018-10-01.tif")
hrsl18 <- crop(hrsl18, bbox)
hrsl18 <- projectRaster(hrsl18, crs=proj.utm31)
slums$hrsl18 <- exact_extract(hrsl18, slums, 'sum')
#lag_grid$hrsl18 <- exact_extract(hrsl18, lag_grid, 'sum')
rm(hrsl18)

#WP peanutButter (2020)
wp_pb20 <- raster("wp_pb/NGA_population_202010271652.tif")
wp_pb20 <- crop(wp_pb20, bbox)
wp_pb20 <- projectRaster(wp_pb20, crs=proj.utm31)
slums$wp_pb20 <- exact_extract(wp_pb20, slums, 'sum')
#lag_grid$wp_pb20 <- exact_extract(wp_pb20, lag_grid, 'sum')
rm(wp_pb20)

#WP UN-adj unconstrained (2018)
wp_u18 <- raster("wp_u/nga_ppp_2018_UNadj.tif")
wp_u18 <- crop(wp_u18, bbox)
wp_u18 <- projectRaster(wp_u18, crs=proj.utm31)
slums$wp_u18 <- exact_extract(wp_u18, slums, 'sum')
#lag_grid$wp_u18 <- exact_extract(wp_u18, lag_grid, 'sum')
rm(wp_u18)

#WP UN-adj unconstrained (2015)
wp_u15 <- raster("wp_u/nga_ppp_2015_UNadj.tif")
wp_u15 <- crop(wp_u15, bbox)
wp_u15 <- projectRaster(wp_u15, crs=proj.utm31)
slums$wp_u15 <- exact_extract(wp_u15, slums, 'sum')
#lag_grid$wp_u15 <- exact_extract(wp_u15, lag_grid, 'sum')
rm(wp_u15)

#WP UN-adj constrained (2020)
wp_c20 <- raster("wp_c/nga_ppp_2020_UNadj_constrained.tif")
wp_c20 <- crop(wp_c20, bbox)
wp_c20 <- projectRaster(wp_c20, crs=proj.utm31)
slums$wp_c20 <- exact_extract(wp_c20, slums, 'sum')
#lag_grid$wp_c20 <- exact_extract(wp_c20, lag_grid, 'sum')
rm(wp_c20)

#WPE (2016)
wpe16 <- raster("wpe/wpe_31n_2.tif")
wpe16 <- projectRaster(wpe16, crs=proj.utm31)
slums$wpe16 <- exact_extract(wpe16, slums, 'sum')
#lag_grid$wpe16 <- exact_extract(wpe16, lag_grid, 'sum')
rm(wpe16)

#LandScan (2018)
ls18 <- raster("landscan/LandScan Global 2018/lspop2018/w001001.adf")
ls18 <- crop(ls18, bbox)
#writeRaster(ls18, "landscan/LandScan Global 2018/ls18.tif") # for visuals in ArcGIS
ls18 <- projectRaster(ls18, crs=proj.utm31)
slums$ls18 <- exact_extract(ls18, slums, 'sum')
#lag_grid$ls18 <- exact_extract(ls18, lag_grid, 'sum')
rm(ls18)

#LandScan (2015)
ls15 <- raster("landscan/LandScan Global 2015/lspop2015/w001001.adf")
ls15 <- crop(ls15, bbox)
ls15 <- projectRaster(ls15, crs=proj.utm31)
slums$ls15 <- exact_extract(ls15, slums, 'sum')
#lag_grid$ls15 <- exact_extract(ls15, lag_grid, 'sum')
rm(ls15)

#GRID3 Nigeria v1.2 (2016)
grid16 <- raster("grid3/NGA - population - v1.2 - gridded/NGA_population_v1_2_gridded.tif")
grid16 <- crop(grid16, bbox)
grid16 <- projectRaster(grid16, crs=proj.utm31)
slums$grid16 <- exact_extract(grid16, slums, 'sum')
#lag_grid$grid16 <- exact_extract(grid16, lag_grid, 'sum')
rm(grid16)

slums_lag <- slums
slums_lag <- slums_lag[,-9]

#summary(lag_grid)
#rm(lag_grid)


##### PORT HARCOURT #####

# Port Harcourt fishnet 200x200m 
#phc_grid <- readOGR("working/phc_fish2.shp")
#phc_grid <- st_as_sf(spTransform(phc_grid, CRS = proj.utm32))

# KnowYourCity Campaign slum boundaries
slums <- readOGR("kyc_campaign/kyc_nga.shp")
slums <- slums[slums$Id > 35,] 

#create unprojected bounding box to crop rasters
bbox <- bb(slums, ext=1.2, output=c("extent"))

#create sf object for exact_extract calcs
slums <- st_as_sf(spTransform(slums, CRS = proj.utm32))

#GPWv4.11 (2020)
gpw20 <- raster("gpw/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")
gpw20 <- crop(gpw20, bbox)
gpw20 <- projectRaster(gpw20, crs=proj.utm32)
slums$gpw20 <- exact_extract(gpw20, slums, 'sum')
#phc_grid$gpw20 <- exact_extract(gpw20, phc_grid, 'sum')
rm(gpw20)

#GPWv4.11 (2015)
gpw15 <- raster("gpw/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec.tif")
gpw15 <- crop(gpw15, bbox)
gpw15 <- projectRaster(gpw15, crs=proj.utm32)
slums$gpw15 <- exact_extract(gpw15, slums, 'sum')
#phc_grid$gpw15 <- exact_extract(gpw15, phc_grid, 'sum')
rm(gpw15)

#GHS-POP (2015)
ghs15 <- raster("ghspop/GHS_POP_E2015_GLOBE_R2019A_54009_250_V1_0_18_8.tif")
ghs15 <- projectRaster(ghs15, crs=proj.utm32)
slums$ghs15 <- exact_extract(ghs15, slums, 'sum')
#phc_grid$ghs15 <- exact_extract(ghs15, phc_grid, 'sum')
rm(ghs15)

#HRSL-Facebook (2018)
hrsl18 <- raster("hrsl/population_nga_2018-10-01.tif")
hrsl18 <- crop(hrsl18, bbox)
hrsl18 <- projectRaster(hrsl18, crs=proj.utm32)
slums$hrsl18 <- exact_extract(hrsl18, slums, 'sum')
#phc_grid$hrsl18 <- exact_extract(hrsl18, phc_grid, 'sum')
rm(hrsl18)

#WP peanutButter (2020)
wp_pb20 <- raster("wp_pb/NGA_population_202010271652.tif")
wp_pb20 <- crop(wp_pb20, bbox)
wp_pb20 <- projectRaster(wp_pb20, crs=proj.utm32)
slums$wp_pb20 <- exact_extract(wp_pb20, slums, 'sum')
#phc_grid$wp_pb20 <- exact_extract(wp_pb20, phc_grid, 'sum')
rm(wp_pb20)

#WP UN-adj unconstrained (2018)
wp_u18 <- raster("wp_u/nga_ppp_2018_UNadj.tif")
wp_u18 <- crop(wp_u18, bbox)
wp_u18 <- projectRaster(wp_u18, crs=proj.utm32)
slums$wp_u18 <- exact_extract(wp_u18, slums, 'sum')
#phc_grid$wp_u18 <- exact_extract(wp_u18, phc_grid, 'sum')
rm(wp_u18)

#WP UN-adj unconstrained (2015)
wp_u15 <- raster("wp_u/nga_ppp_2015_UNadj.tif")
wp_u15 <- crop(wp_u15, bbox)
wp_u15 <- projectRaster(wp_u15, crs=proj.utm32)
slums$wp_u15 <- exact_extract(wp_u15, slums, 'sum')
#phc_grid$wp_u15 <- exact_extract(wp_u15, phc_grid, 'sum')
rm(wp_u15)

#WP UN-adj constrained (2020)
wp_c20 <- raster("wp_c/nga_ppp_2020_UNadj_constrained.tif")
wp_c20 <- crop(wp_c20, bbox)
wp_c20 <- projectRaster(wp_c20, crs=proj.utm32)
slums$wp_c20 <- exact_extract(wp_c20, slums, 'sum')
#phc_grid$wp_c20 <- exact_extract(wp_c20, phc_grid, 'sum')
rm(wp_c20)

#WPE (2016)
wpe16 <- raster("wpe/wpe_32n_2.tif")
wpe16 <- projectRaster(wpe16, crs=proj.utm32)
slums$wpe16 <- exact_extract(wpe16, slums, 'sum')
#phc_grid$wpe16 <- exact_extract(wpe16, phc_grid, 'sum')
rm(wpe16)

#LandScan (2018)
ls18 <- raster("landscan/LandScan Global 2018/lspop2018/w001001.adf")
ls18 <- crop(ls18, bbox)
ls18 <- projectRaster(ls18, crs=proj.utm32)
slums$ls18 <- exact_extract(ls18, slums, 'sum')
#phc_grid$ls18 <- exact_extract(ls18, phc_grid, 'sum')
rm(ls18)

#LandScan (2015)
ls15 <- raster("landscan/LandScan Global 2015/lspop2015/w001001.adf")
ls15 <- crop(ls15, bbox)
ls15 <- projectRaster(ls15, crs=proj.utm32)
slums$ls15 <- exact_extract(ls15, slums, 'sum')
#phc_grid$ls15 <- exact_extract(ls15, phc_grid, 'sum')
rm(ls15)

#GRID3 Nigeria v1.2 (2016)
grid16 <- raster("grid3/NGA - population - v1.2 - gridded/NGA_population_v1_2_gridded.tif")
grid16 <- crop(grid16, bbox)
grid16 <- projectRaster(grid16, crs=proj.utm32)
slums$grid16 <- exact_extract(grid16, slums, 'sum')
#phc_grid$grid16 <- exact_extract(grid16, phc_grid, 'sum')
rm(grid16)

slums_phc <- slums
slums_phc <- slums_phc[,-9]

#summary(phc_grid)
#rm(phc_grid)


##### NAIROBI #####

# Nairobi fishnet 200x200m 
#nai_grid <- readOGR("working/nai_fish2.shp")
#nai_grid <- st_as_sf(spTransform(nai_grid, CRS = proj.utm37))

# KnowYourCity Campaign slum boundaries
slums <- readOGR("kyc_campaign/kyc_ken.shp")

#create unprojected bounding box to crop rasters
bbox <- bb(slums, ext=1.2, output=c("extent"))

#create sf object for exact_extract calcs
slums <- st_as_sf(spTransform(slums, CRS = proj.utm37))

#GPWv4.11 (2020)
gpw20 <- raster("gpw/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")
gpw20 <- crop(gpw20, bbox)
gpw20 <- projectRaster(gpw20, crs=proj.utm37)
slums$gpw20 <- exact_extract(gpw20, slums, 'sum')
#nai_grid$gpw20 <- exact_extract(gpw20, nai_grid, 'sum')
rm(gpw20)

#GPWv4.11 (2015)
gpw15 <- raster("gpw/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec.tif")
gpw15 <- crop(gpw15, bbox)
gpw15 <- projectRaster(gpw15, crs=proj.utm37)
slums$gpw15 <- exact_extract(gpw15, slums, 'sum')
#nai_grid$gpw15 <- exact_extract(gpw15, nai_grid, 'sum')
rm(gpw15)

#GHS-POP (2015)
ghs15 <- raster("ghspop/GHS_POP_E2015_GLOBE_R2019A_54009_250_V1_0_21_9.tif")
ghs15 <- projectRaster(ghs15, crs=proj.utm37)
slums$ghs15 <- exact_extract(ghs15, slums, 'sum')
#nai_grid$ghs15 <- exact_extract(ghs15, nai_grid, 'sum')
rm(ghs15)

#HRSL-Facebook (2018)
hrsl18 <- raster("hrsl/population_AF20_2018-10-01.tif")
hrsl18 <- crop(hrsl18, bbox)
hrsl18 <- projectRaster(hrsl18, crs=proj.utm37)
slums$hrsl18 <- exact_extract(hrsl18, slums, 'sum')
#nai_grid$hrsl18 <- exact_extract(hrsl18, nai_grid, 'sum')
rm(hrsl18)

#WP peanutButter (2020)
wp_pb20 <- raster("wp_pb/KEN_population_202012031254.tif")
wp_pb20 <- crop(wp_pb20, bbox)
wp_pb20 <- projectRaster(wp_pb20, crs=proj.utm37)
slums$wp_pb20 <- exact_extract(wp_pb20, slums, 'sum')
#nai_grid$wp_pb20 <- exact_extract(wp_pb20, nai_grid, 'sum')
rm(wp_pb20)

#WP UN-adj unconstrained (2018)
wp_u18 <- raster("wp_u/ken_ppp_2018_UNadj.tif")
wp_u18 <- crop(wp_u18, bbox)
wp_u18 <- projectRaster(wp_u18, crs=proj.utm37)
slums$wp_u18 <- exact_extract(wp_u18, slums, 'sum')
#nai_grid$wp_u18 <- exact_extract(wp_u18, nai_grid, 'sum')
rm(wp_u18)

#WP UN-adj unconstrained (2015)
wp_u15 <- raster("wp_u/ken_ppp_2015_UNadj.tif")
wp_u15 <- crop(wp_u15, bbox)
wp_u15 <- projectRaster(wp_u15, crs=proj.utm37)
slums$wp_u15 <- exact_extract(wp_u15, slums, 'sum')
#nai_grid$wp_u15 <- exact_extract(wp_u15, nai_grid, 'sum')
rm(wp_u15)

#WP UN-adj constrained (2020)
wp_c20 <- raster("wp_c/ken_ppp_2020_UNadj_constrained.tif")
wp_c20 <- crop(wp_c20, bbox)
wp_c20 <- projectRaster(wp_c20, crs=proj.utm37)
slums$wp_c20 <- exact_extract(wp_c20, slums, 'sum')
#nai_grid$wp_c20 <- exact_extract(wp_c20, nai_grid, 'sum')
rm(wp_c20)

#WPE (2016)
wpe16 <- raster("wpe/wpe_37s_2.tif")
wpe16 <- projectRaster(wpe16, crs=proj.utm37)
slums$wpe16 <- exact_extract(wpe16, slums, 'sum')
#nai_grid$wpe16 <- exact_extract(wpe16, nai_grid, 'sum')
rm(wpe16)

#LandScan (2018)
ls18 <- raster("landscan/LandScan Global 2018/lspop2018/w001001.adf")
ls18 <- crop(ls18, bbox)
#writeRaster(ls18, "landscan/LandScan Global 2018/ls18.tif") # for visuals in ArcGIS
ls18 <- projectRaster(ls18, crs=proj.utm37)
slums$ls18 <- exact_extract(ls18, slums, 'sum')
#nai_grid$ls18 <- exact_extract(ls18, nai_grid, 'sum')
rm(ls18)

#LandScan (2015)
ls15 <- raster("landscan/LandScan Global 2015/lspop2015/w001001.adf")
ls15 <- crop(ls15, bbox)
ls15 <- projectRaster(ls15, crs=proj.utm37)
slums$ls15 <- exact_extract(ls15, slums, 'sum')
#nai_grid$ls15 <- exact_extract(ls15, nai_grid, 'sum')
rm(ls15)


slums_nai <- slums
slums_nai <- slums_nai[,-13]

#summary(nai_grid)
#rm(nai_grid)


######################
## Step 2. Calculations 
######################

#Calculate area of digitised boundary
#LAGOS
slums_lag$area_m2 <- st_area(slums_lag)
slums_lag <- as.data.frame(slums_lag)   
slums_lag <- slums_lag[,c(-13)]                 # remove geometry
slums_lag$area_m2 <- as.numeric(slums_lag$area_m2) # remove units for comparisons
colnames(slums_lag) <- c("Id", "country", "city", "settlement", "kyc_date", "kyc_day", "kyc_month", "kyc_year", "kyc_pop", "kyc_area_acre", "kyc_density", "kyc_structures", "gpw20", "gpw15", "ghs15", "hrsl18", "wp_pb20", "wp_u18", "wp_u15", "wp_c20", "wpe16", "ls18", "ls15", "grid16", "area_m2")

#PORT HARCOURT
slums_phc$area_m2 <- st_area(slums_phc)
slums_phc <- as.data.frame(slums_phc)   
slums_phc <- slums_phc[,c(-13)]                 # remove geometry
slums_phc$area_m2 <- as.numeric(slums_phc$area_m2) # remove units for comparisons
colnames(slums_phc) <- c("Id", "country", "city", "settlement", "kyc_date", "kyc_day", "kyc_month", "kyc_year", "kyc_pop", "kyc_area_acre", "kyc_density", "kyc_structures", "gpw20", "gpw15", "ghs15", "hrsl18", "wp_pb20", "wp_u18", "wp_u15", "wp_c20", "wpe16", "ls18", "ls15", "grid16", "area_m2")

#NAIROBI
slums_nai$grid16 <- NA
slums_nai$area_m2 <- st_area(slums_nai)
slums_nai <- as.data.frame(slums_nai)  
slums_nai <- slums_nai[,c(-13)]                 # remove geometry
slums_nai$area_m2 <- as.numeric(slums_nai$area_m2) # remove units for comparisons
colnames(slums_nai) <- c("Id", "country", "city", "settlement", "kyc_date", "kyc_day", "kyc_month", "kyc_year", "kyc_pop", "kyc_area_acre", "kyc_density", "kyc_structures", "gpw20", "gpw15", "ghs15", "hrsl18", "wp_pb20", "wp_u18", "wp_u15", "wp_c20", "wpe16", "ls18", "ls15", "grid16", "area_m2")

#combine data from 3 cities and clean up
slums <- rbind(slums_lag, slums_phc, slums_nai)
rm(bbox, slums_lag, slums_phc, slums_nai)

#convert KYC acres to m^2
slums$kyc_area_acre[is.na(slums$kyc_area_acre)] <- 0
slums$kyc_area_m2 <- as.numeric(slums$kyc_area_acre) / 0.00024711

#percent difference in KYC area vs digitised area
slums$area_perdiff <- (slums$area_m2 - slums$kyc_area_m2) / slums$kyc_area_m2 * 100
#write.csv(slums, file="working/slums_v4.csv")   # for data quality check of reported vs digitised area

#mean area
med_own_area_m <- sqrt(median(slums$area_m2))
med_own_area_m # 232 x 232 m
med_kyc_area_m <- sqrt(median(slums$kyc_area_m2))
med_kyc_area_m # 129 x 129 m
#use a fishnet of 200x200m to evaluate gridded pop densities over whole cities


#compare same gridded population dataset for different years
slums$comp_gwp <- slums$gpw20 - slums$gpw15
slums$comp_wp_u <- slums$wp_u18 - slums$wp_u15
slums$comp_ls <- slums$ls18 - slums$ls15

#align gridded estimates by year of KYC data collection
slums$gpw_combo <- ifelse(slums$kyc_year>=2017, slums$gpw20, slums$gpw15)
slums$wp_u_combo <- ifelse(slums$kyc_year>=2017, slums$wp_u18, slums$wp_u15)
slums$ls_combo <- ifelse(slums$kyc_year>=2017, slums$ls18, slums$ls15)
slums$ghs_combo <- as.numeric(ifelse(slums$kyc_year>=2017, "" , slums$ghs15))
slums$wpe_combo <- as.numeric(ifelse(slums$kyc_year>=2017, "" , slums$wpe16))
slums$grid_combo <- as.numeric(ifelse(slums$kyc_year>=2017, "" , slums$grid16))
slums$hrsl_combo <- as.numeric(ifelse(slums$kyc_year>=2017, slums$hrsl18 , ""))
slums$wp_pb_combo <- as.numeric(ifelse(slums$kyc_year>=2017, slums$wp_pb20 , ""))
slums$wp_c_combo <- as.numeric(ifelse(slums$kyc_year>=2017, slums$wp_c20 , ""))


#differences btw KYC population & gridded population estimates
slums$diff_gpw <- slums$gpw_combo - slums$kyc_pop
slums$diff_wp_u <- slums$wp_u_combo - slums$kyc_pop
slums$diff_ls <- slums$ls_combo - slums$kyc_pop

slums$diff_ghs <- slums$ghs_combo - slums$kyc_pop
slums$diff_wpe <- slums$wpe_combo - slums$kyc_pop
slums$diff_grid <- slums$grid_comb - slums$kyc_pop

slums$diff_hrsl <- slums$hrsl_combo - slums$kyc_pop
slums$diff_wp_pb <- slums$wp_pb_combo - slums$kyc_pop
slums$diff_wp_c <- slums$wp_c_comb - slums$kyc_pop

#absolute differences btw KYC population & gridded population estimates
slums$adiff_gpw <- abs(slums$diff_gpw)
slums$adiff_wp_u <- abs(slums$diff_wp_u)
slums$adiff_ls <- abs(slums$diff_ls)

slums$adiff_ghs <- abs(slums$diff_ghs)
slums$adiff_wpe <- abs(slums$diff_wpe)
slums$adiff_grid <- abs(slums$diff_grid)

slums$adiff_hrsl <- abs(slums$diff_hrsl)
slums$adiff_wp_pb <- abs(slums$diff_wp_pb)
slums$adiff_wp_c <- abs(slums$diff_wp_c)


#differences btw KYC population & gridded population estimates - squared
slums$diff2_gpw <- slums$diff_gpw^2
slums$diff2_wp_u <- slums$diff_wp_u^2
slums$diff2_ls <- slums$diff_ls^2

slums$diff2_ghs <- slums$diff_ghs^2
slums$diff2_wpe <- slums$diff_wpe^2
slums$diff2_grid <- slums$diff_grid^2

slums$diff2_hrsl <- slums$diff_hrsl^2
slums$diff2_wp_pb <- slums$diff_wp_pb^2
slums$diff2_wp_c <- slums$diff_wp_c^2

#percent differences btw KYC population & gridded population estimates 
slums$pdiff_gpw <- (slums$gpw_combo - slums$kyc_pop) / slums$kyc_pop * 100
slums$pdiff_wp_u <- (slums$wp_u_combo - slums$kyc_pop) / slums$kyc_pop * 100
slums$pdiff_ls <- (slums$ls_combo - slums$kyc_pop) / slums$kyc_pop * 100

slums$pdiff_ghs <- (slums$ghs_combo - slums$kyc_pop) / slums$kyc_pop * 100
slums$pdiff_wpe <- (slums$wpe_combo - slums$kyc_pop) / slums$kyc_pop * 100
slums$pdiff_grid <- (slums$grid_combo - slums$kyc_pop) / slums$kyc_pop * 100

slums$pdiff_hrsl <- (slums$hrsl_combo - slums$kyc_pop) / slums$kyc_pop * 100
slums$pdiff_wp_pb <- (slums$wp_pb_combo - slums$kyc_pop) / slums$kyc_pop * 100
slums$pdiff_wp_c <- (slums$wp_c_combo - slums$kyc_pop) / slums$kyc_pop * 100


#fraction of KYC population estimated by gridded population 
slums$frac_gpw <- (slums$gpw_combo / slums$kyc_pop) 
slums$frac_wp_u <- (slums$wp_u_combo / slums$kyc_pop) 
slums$frac_ls <- (slums$ls_combo / slums$kyc_pop)

slums$frac_ghs <- (slums$ghs_combo / slums$kyc_pop) 
slums$frac_wpe <- (slums$wpe_combo / slums$kyc_pop) 
slums$frac_grid <- (slums$grid_combo / slums$kyc_pop) 

slums$frac_hrsl <- (slums$hrsl_combo / slums$kyc_pop) 
slums$frac_wp_pb <- (slums$wp_pb_combo / slums$kyc_pop) 
slums$frac_wp_c <- (slums$wp_c_combo / slums$kyc_pop)


# Densities per 200x200m

slums$dens_kyc <- slums$kyc_pop / slums$area_m2 * 40000

slums$dens_gpw <- slums$gpw_combo / slums$area_m2 * 40000
slums$dens_wp_u <- slums$wp_u_combo / slums$area_m2 * 40000
slums$dens_ls <- slums$ls_combo / slums$area_m2 * 40000

slums$dens_ghs <- slums$ghs_combo / slums$area_m2 * 40000
slums$dens_wpe <- slums$wpe_combo / slums$area_m2 * 40000
slums$dens_grid <- slums$grid_combo / slums$area_m2 * 40000

slums$dens_hrsl <- slums$hrsl_combo / slums$area_m2 * 40000
slums$dens_wp_pb <- slums$wp_pb_combo / slums$area_m2 * 40000
slums$dens_wp_c <- slums$wp_c_combo / slums$area_m2 * 40000



######################
## Step 3. Data checks
######################

#graphs & t-tests
graph_gpw <- ggplot(data=slums, aes(x=reorder(Id,comp_gwp), y=comp_gwp)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
graph_gpw  + labs(x="Settlement", y="Difference in population estimate (2015-2020)")
ggsave("comp_gpw_2.tiff", units="in", width=7, height=5, dpi=300, compression = 'lzw')

t.test(x=slums$gpw20, y=slums$gpw15, paired = TRUE, alternative = "two.sided")

graph_wp <- ggplot(data=slums, aes(x=reorder(Id,comp_wp_u), y=comp_wp_u)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
graph_wp  + labs(x="Settlement", y="Difference in population estimate (2015-2018)")
ggsave("comp_wp_2.tiff", units="in", width=7, height=5, dpi=300, compression = 'lzw')

t.test(x=slums$wp_u18, y=slums$wp_u15, paired = TRUE, alternative = "two.sided")

graph_ls <- ggplot(data=slums, aes(x=reorder(Id,comp_ls), y=comp_ls)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
graph_ls  + labs(x="Settlement", y="Difference in population estimate (2015-2018)")
ggsave("comp_ls_2.tiff", units="in", width=7, height=5, dpi=300, compression = 'lzw')

t.test(x=slums$ls18, y=slums$ls15, paired = TRUE, alternative = "two.sided")


######################
## Step 4. Summarise "true" population
######################

# Exclude settlements after visual inspection suggests errors
slums <- slums[!(slums$Id==1 | slums$Id==3 | slums$Id==14 | slums$Id==17 | slums$Id==19 | slums$Id==24 | slums$Id==86 | slums$Id==90 | slums$Id==103 | slums$Id==109 | slums$Id==110 | slums$Id==111 | slums$Id==114 | slums$Id==115 | slums$Id==131 | slums$Id==136),]

#summarise density
tapply(slums$kyc_pop, slums$city, summary)
tapply(slums$area_m2, slums$city, summary)
tapply(slums$dens_kyc, slums$city, summary)

# Create orderid for graphs
slums <- slums[order(slums$city, slums$kyc_pop), ]
slums$orderid <- data.table::rowid(slums$city)

slums <- slums[order(slums$city, slums$dens_kyc), ]
slums$orderdens <- data.table::rowid(slums$city)


# quick visual check of population and density calcs
graph_kyc_pop <- ggplot(slums, aes(x=orderid, y=kyc_pop)) + 
  geom_line(aes(color=city), size=1) +
  geom_point(size=1, alpha=0.3)
graph_kyc_pop

graph_kyc_dens <- ggplot(slums, aes(x=orderdens, y=dens_kyc)) + 
  geom_line(aes(color=city), size=1) +
  geom_point(size=1, alpha=0.3)
graph_kyc_dens


######################
## Step 5. Summarise population estimation errors
######################

#### Estimates 2013 - 2020 by city
slums_lag_2015 <- slums[which(slums$city=="Lagos" & slums$kyc_year<2017),]
slums_lag_2018 <- slums[which(slums$city=="Lagos" & slums$kyc_year>=2017),]

slums_phc_2015 <- slums[which(slums$city=="Port Harcourt" & slums$kyc_year<2017),]
slums_phc_2018 <- slums[which(slums$city=="Port Harcourt" & slums$kyc_year>=2017),]

slums_nai_2015 <- slums[which(slums$city=="Nairobi" & slums$kyc_year<2017),]
slums_nai_2018 <- slums[which(slums$city=="Nairobi" & slums$kyc_year>=2017),]

colors2015 <- c("KYC Reported" = "black", "GPWv4.11"="skyblue3", "WP-Unconstr"="darkred", "LandScan"="seagreen4", "GHS-POP"="hotpink2", "WPE"="darkorange2", "GRID3"="grey45")
colors2018 <- c("KYC Reported" = "black", "GPWv4.11"="skyblue3", "WP-Unconstr"="darkred", "LandScan"="seagreen4", "HRSL"="gold", "WP-PButter"="maroon3", "WP-Constr"="dodgerblue4")


###### LAGOS population comparisons

#test <- ggplot(slums_lag_2015, aes(x=reorder(Id,orderid),y=kyc_pop), size=1.25) + geom_line() + geom_point()
#test

# ------>  FORREST  <-------  In this section I would like to show Settlement ID on the X axis rather than an arbitrary order #


slums_lag_2015 <- slums_lag_2015[order(slums_lag_2015$kyc_pop), ]
slums_lag_2015$orderid <- data.table::rowid(slums_lag_2015$city)
graph_lag_2015 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_lag_2015, aes(x=orderid, y=kyc_pop, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_lag_2015, aes(x=orderid, y=kyc_pop), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_lag_2015, aes(x=orderid, y=gpw_combo, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderid, y=gpw_combo), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_lag_2015, aes(x=orderid, y=wp_u_combo, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderid, y=wp_u_combo), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_lag_2015, aes(x=orderid, y=ls_combo, color="LandScan"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderid, y=ls_combo), size=1, alpha=0.3) +
  # GHS-POP
  geom_line(data=slums_lag_2015, aes(x=orderid, y=ghs_combo, color="GHS-POP"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderid, y=ghs_combo), size=1, alpha=0.3)+
  # WPE
  geom_line(data=slums_lag_2015, aes(x=orderid, y=wpe_combo, color="WPE"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderid, y=wpe_combo), size=1, alpha=0.3) +
  # GRID3
  geom_line(data=slums_lag_2015, aes(x=orderid, y=grid_combo, color="GRID3"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderid, y=grid_combo), size=1, alpha=0.3) +
  ylim(0,30000) +
  labs(x = "Settlement #", y = "Population", color = "Legend") +
  scale_color_manual(values = colors2015) +
  ggtitle("Lagos 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_lag_2015
ggsave("graph_lag_2015.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')


slums_lag_2018 <- slums_lag_2018[order(slums_lag_2018$kyc_pop), ]
slums_lag_2018$orderid <- data.table::rowid(slums_lag_2018$city)
graph_lag_2018 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_lag_2018, aes(x=orderid, y=kyc_pop, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_lag_2018, aes(x=orderid, y=kyc_pop), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_lag_2018, aes(x=orderid, y=gpw_combo, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderid, y=gpw_combo), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_lag_2018, aes(x=orderid, y=wp_u_combo, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderid, y=wp_u_combo), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_lag_2018, aes(x=orderid, y=ls_combo, color="LandScan"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderid, y=ls_combo), size=1, alpha=0.3) +
  # HRSL
  geom_line(data=slums_lag_2018, aes(x=orderid, y=hrsl_combo, color="HRSL"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderid, y=hrsl_combo), size=1, alpha=0.3) +
  # WorldPop-PeanutButter
  geom_line(data=slums_lag_2018, aes(x=orderid, y=wp_pb_combo, color="WP-PButter"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderid, y=wp_pb_combo), size=1, alpha=0.3) +
  # WorldPop-Constrained
  geom_line(data=slums_lag_2018, aes(x=orderid, y=wp_c_combo, color="WP-Constr"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderid, y=wp_c_combo), size=1, alpha=0.3) +
  ylim(0,30000) +
  labs(x = "Settlement #", y = "Population", color = "Legend") +
  scale_color_manual(values = colors2018) +
  ggtitle("Lagos 2017-2020") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_lag_2018
ggsave("graph_lag_2018.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')


###### PORT HARCOURT population comparisons

slums_phc_2015 <- slums_phc_2015[order(slums_phc_2015$kyc_pop), ]
slums_phc_2015$orderid <- data.table::rowid(slums_phc_2015$city)
graph_phc_2015 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_phc_2015, aes(x=orderid, y=kyc_pop, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_phc_2015, aes(x=orderid, y=kyc_pop), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_phc_2015, aes(x=orderid, y=gpw_combo, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderid, y=gpw_combo), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_phc_2015, aes(x=orderid, y=wp_u_combo, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderid, y=wp_u_combo), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_phc_2015, aes(x=orderid, y=ls_combo, color="LandScan"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderid, y=ls_combo), size=1, alpha=0.3) +
  # GHS-POP
  geom_line(data=slums_phc_2015, aes(x=orderid, y=ghs_combo, color="GHS-POP"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderid, y=ghs_combo), size=1, alpha=0.3)+
  # WPE
  geom_line(data=slums_phc_2015, aes(x=orderid, y=wpe_combo, color="WPE"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderid, y=wpe_combo), size=1, alpha=0.3) +
  # GRID3
  geom_line(data=slums_phc_2015, aes(x=orderid, y=grid_combo, color="GRID3"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderid, y=grid_combo), size=1, alpha=0.3) +
  ylim(0,30000) +
  labs(x = "Settlement #", y = "Population", color = "Legend") +
  scale_color_manual(values = colors2015) +
  ggtitle("Port Harcourt 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_phc_2015
ggsave("graph_phc_2015.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')


slums_phc_2018 <- slums_phc_2018[order(slums_phc_2018$kyc_pop), ]
slums_phc_2018$orderid <- data.table::rowid(slums_phc_2018$city)
graph_phc_2018 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_phc_2018, aes(x=orderid, y=kyc_pop, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_phc_2018, aes(x=orderid, y=kyc_pop), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_phc_2018, aes(x=orderid, y=gpw_combo, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderid, y=gpw_combo), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_phc_2018, aes(x=orderid, y=wp_u_combo, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderid, y=wp_u_combo), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_phc_2018, aes(x=orderid, y=ls_combo, color="LandScan"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderid, y=ls_combo), size=1, alpha=0.3) +
  # HRSL
  geom_line(data=slums_phc_2018, aes(x=orderid, y=hrsl_combo, color="HRSL"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderid, y=hrsl_combo), size=1, alpha=0.3) +
  # WorldPop-PeanutButter
  geom_line(data=slums_phc_2018, aes(x=orderid, y=wp_pb_combo, color="WP-PButter"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderid, y=wp_pb_combo), size=1, alpha=0.3) +
  # WorldPop-Constrained
  geom_line(data=slums_phc_2018, aes(x=orderid, y=wp_c_combo, color="WP-Constr"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderid, y=wp_c_combo), size=1, alpha=0.3) +
  ylim(0,30000) +
  labs(x = "Settlement #", y = "Population", color = "Legend") +
  scale_color_manual(values = colors2018) +
  ggtitle("Port Harcourt 2017-2020") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_phc_2018
ggsave("graph_phc_2018.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')


###### NAIROBI pouplation comparisons

slums_nai_2015 <- slums_nai_2015[order(slums_nai_2015$kyc_pop), ]
slums_nai_2015$orderid <- data.table::rowid(slums_nai_2015$city)
graph_nai_2015 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_nai_2015, aes(x=orderid, y=kyc_pop, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_nai_2015, aes(x=orderid, y=kyc_pop), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_nai_2015, aes(x=orderid, y=gpw_combo, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderid, y=gpw_combo), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_nai_2015, aes(x=orderid, y=wp_u_combo, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderid, y=wp_u_combo), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_nai_2015, aes(x=orderid, y=ls_combo, color="LandScan"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderid, y=ls_combo), size=1, alpha=0.3) +
  # GHS-POP
  geom_line(data=slums_nai_2015, aes(x=orderid, y=ghs_combo, color="GHS-POP"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderid, y=ghs_combo), size=1, alpha=0.3)+
  # WPE
  geom_line(data=slums_nai_2015, aes(x=orderid, y=wpe_combo, color="WPE"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderid, y=wpe_combo), size=1, alpha=0.3) +
  ylim(0,30000) +
  labs(x = "Settlement #", y = "Population", color = "Legend") +
  scale_color_manual(values = colors2015) +
  ggtitle("Nairobi 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_nai_2015
ggsave("graph_nai_2015.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')
# run a version without ylim for appendix

slums_nai_2018 <- slums_nai_2018[order(slums_nai_2018$kyc_pop), ]
slums_nai_2018$orderid <- data.table::rowid(slums_nai_2018$city)
graph_nai_2018 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_nai_2018, aes(x=orderid, y=kyc_pop, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_nai_2018, aes(x=orderid, y=kyc_pop), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_nai_2018, aes(x=orderid, y=gpw_combo, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderid, y=gpw_combo), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_nai_2018, aes(x=orderid, y=wp_u_combo, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderid, y=wp_u_combo), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_nai_2018, aes(x=orderid, y=ls_combo, color="LandScan"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderid, y=ls_combo), size=1, alpha=0.3) +
  # HRSL
  geom_line(data=slums_nai_2018, aes(x=orderid, y=hrsl_combo, color="HRSL"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderid, y=hrsl_combo), size=1, alpha=0.3) +
  # WorldPop-PeanutButter
  geom_line(data=slums_nai_2018, aes(x=orderid, y=wp_pb_combo, color="WP-PButter"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderid, y=wp_pb_combo), size=1, alpha=0.3) +
  # WorldPop-Constrained
  geom_line(data=slums_nai_2018, aes(x=orderid, y=wp_c_combo, color="WP-Constr"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderid, y=wp_c_combo), size=1, alpha=0.3) +
  ylim(0,30000) +
  labs(x = "Settlement #", y = "Population", color = "Legend") +
  scale_color_manual(values = colors2018) +
  ggtitle("Nairobi 2017-2020") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_nai_2018
ggsave("graph_nai_2018.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')


####### LAGOS density comparisons

slums_lag_2015 <- slums_lag_2015[order(slums_lag_2015$city, slums_lag_2015$dens_kyc), ]
slums_lag_2015$orderdens <- data.table::rowid(slums_lag_2015$city)
graph_lag_2015 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_lag_2015, aes(x=orderdens, y=dens_kyc, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_lag_2015, aes(x=orderdens, y=dens_kyc), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_lag_2015, aes(x=orderdens, y=dens_gpw, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderdens, y=dens_gpw), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_lag_2015, aes(x=orderdens, y=dens_wp_u, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderdens, y=dens_wp_u), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_lag_2015, aes(x=orderdens, y=dens_ls, color="LandScan"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderdens, y=dens_ls), size=1, alpha=0.3) +
  # GHS-POP
  geom_line(data=slums_lag_2015, aes(x=orderdens, y=dens_ghs, color="GHS-POP"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderdens, y=dens_ghs), size=1, alpha=0.3)+
  # WPE
  geom_line(data=slums_lag_2015, aes(x=orderdens, y=dens_wpe, color="WPE"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderdens, y=dens_wpe), size=1, alpha=0.3) +
  # GRID3
  geom_line(data=slums_lag_2015, aes(x=orderdens, y=dens_grid, color="GRID3"), size=0.75) +
  geom_point(data=slums_lag_2015, aes(x=orderdens, y=dens_grid), size=1, alpha=0.3) +
  ylim(0,.4) +
  labs(x = "Settlement #", y = "Density (pop/m2)", color = "Legend") +
  scale_color_manual(values = colors2015) +
  ggtitle("Lagos 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_lag_2015
ggsave("graph_lag_dens_2015.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')


slums_lag_2018 <- slums_lag_2018[order(slums_lag_2018$city, slums_lag_2018$dens_kyc), ]
slums_lag_2018$orderdens <- data.table::rowid(slums_lag_2018$city)
graph_lag_2018 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_lag_2018, aes(x=orderdens, y=dens_kyc, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_lag_2018, aes(x=orderdens, y=dens_kyc), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_lag_2018, aes(x=orderdens, y=dens_gpw, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderdens, y=dens_gpw), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_lag_2018, aes(x=orderdens, y=dens_wp_u, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderdens, y=dens_wp_u), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_lag_2018, aes(x=orderdens, y=dens_ls, color="LandScan"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderdens, y=dens_ls), size=1, alpha=0.3) +
  # HRSL
  geom_line(data=slums_lag_2018, aes(x=orderdens, y=dens_hrsl, color="HRSL"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderdens, y=dens_hrsl), size=1, alpha=0.3) +
  # WorldPop-PeanutButter
  geom_line(data=slums_lag_2018, aes(x=orderdens, y=dens_wp_pb, color="WP-PButter"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderdens, y=dens_wp_pb), size=1, alpha=0.3) +
  # WorldPop-Constrained
  geom_line(data=slums_lag_2018, aes(x=orderdens, y=dens_wp_c, color="WP-Constr"), size=0.75) +
  geom_point(data=slums_lag_2018, aes(x=orderdens, y=dens_wp_c), size=1, alpha=0.3) +
  ylim(0,.4) +
  labs(x = "Settlement #", y = "Density (pop/m2)", color = "Legend") +
  scale_color_manual(values = colors2018) +
  ggtitle("Lagos 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_lag_2018
ggsave("graph_lag_dens_2018.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')



####### PORT HARCOURT density comparisons

slums_phc_2015 <- slums_phc_2015[order(slums_phc_2015$city, slums_phc_2015$dens_kyc), ]
slums_phc_2015$orderdens <- data.table::rowid(slums_phc_2015$city)
graph_phc_2015 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_phc_2015, aes(x=orderdens, y=dens_kyc, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_phc_2015, aes(x=orderdens, y=dens_kyc), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_phc_2015, aes(x=orderdens, y=dens_gpw, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderdens, y=dens_gpw), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_phc_2015, aes(x=orderdens, y=dens_wp_u, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderdens, y=dens_wp_u), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_phc_2015, aes(x=orderdens, y=dens_ls, color="LandScan"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderdens, y=dens_ls), size=1, alpha=0.3) +
  # GHS-POP
  geom_line(data=slums_phc_2015, aes(x=orderdens, y=dens_ghs, color="GHS-POP"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderdens, y=dens_ghs), size=1, alpha=0.3)+
  # WPE
  geom_line(data=slums_phc_2015, aes(x=orderdens, y=dens_wpe, color="WPE"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderdens, y=dens_wpe), size=1, alpha=0.3) +
  # GRID3
  geom_line(data=slums_phc_2015, aes(x=orderdens, y=dens_grid, color="GRID3"), size=0.75) +
  geom_point(data=slums_phc_2015, aes(x=orderdens, y=dens_grid), size=1, alpha=0.3) +
  ylim(0,.4) +
  labs(x = "Settlement #", y = "Density (pop/m2)", color = "Legend") +
  scale_color_manual(values = colors2015) +
  ggtitle("Port Harcourt 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_phc_2015
ggsave("graph_phc_dens_2015.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')


slums_phc_2018 <- slums_phc_2018[order(slums_phc_2018$city, slums_phc_2018$dens_kyc), ]
slums_phc_2018$orderdens <- data.table::rowid(slums_phc_2018$city)
graph_phc_2018 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_phc_2018, aes(x=orderdens, y=dens_kyc, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_phc_2018, aes(x=orderdens, y=dens_kyc), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_phc_2018, aes(x=orderdens, y=dens_gpw, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderdens, y=dens_gpw), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_phc_2018, aes(x=orderdens, y=dens_wp_u, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderdens, y=dens_wp_u), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_phc_2018, aes(x=orderdens, y=dens_ls, color="LandScan"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderdens, y=dens_ls), size=1, alpha=0.3) +
  # HRSL
  geom_line(data=slums_phc_2018, aes(x=orderdens, y=dens_hrsl, color="HRSL"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderdens, y=dens_hrsl), size=1, alpha=0.3) +
  # WorldPop-PeanutButter
  geom_line(data=slums_phc_2018, aes(x=orderdens, y=dens_wp_pb, color="WP-PButter"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderdens, y=dens_wp_pb), size=1, alpha=0.3) +
  # WorldPop-Constrained
  geom_line(data=slums_phc_2018, aes(x=orderdens, y=dens_wp_c, color="WP-Constr"), size=0.75) +
  geom_point(data=slums_phc_2018, aes(x=orderdens, y=dens_wp_c), size=1, alpha=0.3) +
  ylim(0,.4) +
  labs(x = "Settlement #", y = "Density (pop/m2)", color = "Legend") +
  scale_color_manual(values = colors2018) +
  ggtitle("Port Harcourt 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_phc_2018
ggsave("graph_phc_dens_2018.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')




####### NAIROBI density comparisons

slums_nai_2015 <- slums_nai_2015[order(slums_nai_2015$city, slums_nai_2015$dens_kyc), ]
slums_nai_2015$orderdens <- data.table::rowid(slums_nai_2015$city)
graph_nai_2015 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_nai_2015, aes(x=orderdens, y=dens_kyc, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_nai_2015, aes(x=orderdens, y=dens_kyc), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_nai_2015, aes(x=orderdens, y=dens_gpw, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderdens, y=dens_gpw), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_nai_2015, aes(x=orderdens, y=dens_wp_u, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderdens, y=dens_wp_u), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_nai_2015, aes(x=orderdens, y=dens_ls, color="LandScan"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderdens, y=dens_ls), size=1, alpha=0.3) +
  # GHS-POP
  geom_line(data=slums_nai_2015, aes(x=orderdens, y=dens_ghs, color="GHS-POP"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderdens, y=dens_ghs), size=1, alpha=0.3)+
  # WPE
  geom_line(data=slums_nai_2015, aes(x=orderdens, y=dens_wpe, color="WPE"), size=0.75) +
  geom_point(data=slums_nai_2015, aes(x=orderdens, y=dens_wpe), size=1, alpha=0.3) +
  ylim(0,.9) +
  labs(x = "Settlement #", y = "Density (pop/m2)", color = "Legend") +
  scale_color_manual(values = colors2015) +
  ggtitle("Nairobi 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_nai_2015
ggsave("graph_nai_dens_2015.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')


slums_nai_2018 <- slums_nai_2018[order(slums_nai_2018$city, slums_nai_2018$dens_kyc), ]
slums_nai_2018$orderdens <- data.table::rowid(slums_nai_2018$city)
graph_nai_2018 <- ggplot() + 
  # KYC Campaign
  geom_line(data=slums_nai_2018, aes(x=orderdens, y=dens_kyc, color="KYC Reported"), size=1.25) +
  geom_point(data=slums_nai_2018, aes(x=orderdens, y=dens_kyc), size=1, alpha=0.3) +
  # GWPv4
  geom_line(data=slums_nai_2018, aes(x=orderdens, y=dens_gpw, color="GPWv4.11"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderdens, y=dens_gpw), size=1, alpha=0.3) +
  # WorldPop-Unconstrained
  geom_line(data=slums_nai_2018, aes(x=orderdens, y=dens_wp_u, color="WP-Unconstr"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderdens, y=dens_wp_u), size=1, alpha=0.3) +
  # LandScan
  geom_line(data=slums_nai_2018, aes(x=orderdens, y=dens_ls, color="LandScan"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderdens, y=dens_ls), size=1, alpha=0.3) +
  # HRSL
  geom_line(data=slums_nai_2018, aes(x=orderdens, y=dens_hrsl, color="HRSL"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderdens, y=dens_hrsl), size=1, alpha=0.3) +
  # WorldPop-PeanutButter
  geom_line(data=slums_nai_2018, aes(x=orderdens, y=dens_wp_pb, color="WP-PButter"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderdens, y=dens_wp_pb), size=1, alpha=0.3) +
  # WorldPop-Constrained
  geom_line(data=slums_nai_2018, aes(x=orderdens, y=dens_wp_c, color="WP-Constr"), size=0.75) +
  geom_point(data=slums_nai_2018, aes(x=orderdens, y=dens_wp_c), size=1, alpha=0.3) +
  ylim(0,.9) +
  labs(x = "Settlement #", y = "Density (pop/m2)", color = "Legend") +
  scale_color_manual(values = colors2018) +
  ggtitle("Nairobi 2013-2016") +
  theme(plot.title = element_text(vjust = -10, hjust = 0.5), legend.title=element_blank(), legend.position = c(.13, .69))
graph_nai_2018
ggsave("graph_nai_dens_2018.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw')






# Slum pop MAE by dataset

gpw <- sum(slums$adiff_gpw, na.rm = TRUE) / NROW(which(!is.na(slums$adiff_gpw))) # Note that I previously used "length" instead of "NROW", but "length gave me # of columns
wp_u <- sum(slums$adiff_wp_u, na.rm = TRUE) / NROW(which(!is.na(slums$adiff_wp_u)))
ls <- sum(slums$adiff_ls, na.rm = TRUE) / NROW(which(!is.na(slums$adiff_ls)))

ghs <- sum(slums$adiff_ghs, na.rm = TRUE)  / NROW(which(!is.na(slums$adiff_ghs)))
wpe <- sum(slums$adiff_wpe, na.rm = TRUE) / NROW(which(!is.na(slums$adiff_wpe)))
grid <- sum(slums$adiff_grid, na.rm = TRUE) / NROW(which(!is.na(slums$adiff_grid)))

wp_pb <- sum(slums$adiff_wp_pb, na.rm = TRUE) / NROW(which(!is.na(slums$adiff_wp_pb)))
wp_c <- sum(slums$adiff_wp_c, na.rm = TRUE) / NROW(which(!is.na(slums$adiff_wp_c)))
hrsl <- sum(slums$adiff_hrsl, na.rm = TRUE) / NROW(which(!is.na(slums$adiff_hrsl)))

mae <- as.data.frame(rbind(gpw, wp_u, ls, ghs, wpe, grid, wp_pb, wp_c, hrsl))
names(mae)[1] <- "mae"

# Slum pop RMSE by dataset

gpw <- sqrt( sum(slums$diff2_gpw, na.rm = TRUE) / NROW(which(!is.na(slums$diff2_gpw))) )
wp_u <- sqrt( sum(slums$diff2_wp_u, na.rm = TRUE) / NROW(which(!is.na(slums$diff2_wp_u))) )
ls <- sqrt( sum(slums$diff2_ls, na.rm = TRUE) / NROW(which(!is.na(slums$diff2_ls))) )

ghs <- sqrt( sum(slums$diff2_ghs, na.rm = TRUE)  / NROW(which(!is.na(slums$diff2_ghs))) )
wpe <- sqrt( sum(slums$diff2_wpe, na.rm = TRUE) / NROW(which(!is.na(slums$diff2_wpe))) )
grid <- sqrt( sum(slums$diff2_grid, na.rm = TRUE) / NROW(which(!is.na(slums$diff2_grid))) )

wp_pb <- sqrt( sum(slums$diff2_wp_pb, na.rm = TRUE) / NROW(which(!is.na(slums$diff2_wp_pb))) )
wp_c <- sqrt( sum(slums$diff2_wp_c, na.rm = TRUE) / NROW(which(!is.na(slums$diff2_wp_c))) )
hrsl <- sqrt( sum(slums$diff2_hrsl, na.rm = TRUE) / NROW(which(!is.na(slums$diff2_hrsl))) )

rmse <- as.data.frame(rbind(gpw, wp_u, ls, ghs, wpe, grid, wp_pb, wp_c, hrsl))
names(rmse)[1] <- "rmse"

# Slum pop Bias by dataset

gpw <- mean(slums$diff_gpw, na.rm = TRUE)
wp_u <- mean(slums$diff_wp_u, na.rm = TRUE)
ls <- mean(slums$diff_ls, na.rm = TRUE) 

ghs <- mean(slums$diff_ghs, na.rm = TRUE) 
wpe <- mean(slums$diff_wpe, na.rm = TRUE) 
grid <- mean(slums$diff_grid, na.rm = TRUE)

wp_pb <- mean(slums$diff_wp_pb, na.rm = TRUE)
wp_c <- mean(slums$diff_wp_c, na.rm = TRUE)
hrsl <- mean(slums$diff_hrsl, na.rm = TRUE) 

bias <- as.data.frame(rbind(gpw, wp_u, ls, ghs, wpe, grid, wp_pb, wp_c, hrsl))
names(bias)[1] <- "bias"

# Median fraction of KYC pop estimated by gridded pop

gpw <- median(slums$frac_gpw, na.rm = TRUE)
wp_u <- median(slums$frac_wp_u, na.rm = TRUE)
ls <- median(slums$frac_ls, na.rm = TRUE) 

ghs <- median(slums$frac_ghs, na.rm = TRUE) 
wpe <- median(slums$frac_wpe, na.rm = TRUE) 
grid <- median(slums$frac_grid, na.rm = TRUE)

wp_pb <- median(slums$frac_wp_pb, na.rm = TRUE)
wp_c <- median(slums$frac_wp_c, na.rm = TRUE)
hrsl <- median(slums$frac_hrsl, na.rm = TRUE) 

frac <- as.data.frame(rbind(gpw, wp_u, ls, ghs, wpe, grid, wp_pb, wp_c, hrsl))
names(frac)[1] <- "Median Fraction Pop Est'd"

errorstats_pop <- cbind(mae, rmse, bias, frac)

#write.csv(slums, file="working/slums_v5_long.csv")   # to identify case study areas


####################################    ANALYSIS 2    #######################################



######################
## Step 6. Summarise population estimates
######################


##### LAGOS ######

# KnowYourCity Campaign slum boundaries
slums <- readOGR("Badmos_etal/slum_50m_poly3.shp")

#create unprojected bounding box to crop rasters
bbox <- bb(slums, ext=1.2, output=c("extent"))

#create sf object for exact_extract calcs
slums <- st_as_sf(spTransform(slums, CRS = proj.utm31))

#GPWv4.11 (2020)
gpw20 <- raster("gpw/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_30_sec.tif")
gpw20 <- crop(gpw20, bbox)
gpw20 <- projectRaster(gpw20, crs=proj.utm31)
slums$gpw20 <- exact_extract(gpw20, slums, 'sum')
rm(gpw20)

#GPWv4.11 (2015)
gpw15 <- raster("gpw/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2015_30_sec.tif")
gpw15 <- crop(gpw15, bbox)
gpw15 <- projectRaster(gpw15, crs=proj.utm31)
slums$gpw15 <- exact_extract(gpw15, slums, 'sum')
rm(gpw15)

#GHS-POP (2015)
ghs15 <- raster("ghspop/GHS_POP_E2015_GLOBE_R2019A_54009_250_V1_0_18_8.tif")
ghs15 <- projectRaster(ghs15, crs=proj.utm31)
slums$ghs15 <- exact_extract(ghs15, slums, 'sum')
rm(ghs15)

#HRSL-Facebook (2018)
hrsl18 <- raster("hrsl/population_nga_2018-10-01.tif")
hrsl18 <- crop(hrsl18, bbox)
hrsl18 <- projectRaster(hrsl18, crs=proj.utm31)
slums$hrsl18 <- exact_extract(hrsl18, slums, 'sum')
rm(hrsl18)

#WP peanutButter (2020)
wp_pb20 <- raster("wp_pb/NGA_population_202010271652.tif")
wp_pb20 <- crop(wp_pb20, bbox)
wp_pb20 <- projectRaster(wp_pb20, crs=proj.utm31)
slums$wp_pb20 <- exact_extract(wp_pb20, slums, 'sum')
rm(wp_pb20)

#WP UN-adj unconstrained (2018)
wp_u18 <- raster("wp_u/nga_ppp_2018_UNadj.tif")
wp_u18 <- crop(wp_u18, bbox)
wp_u18 <- projectRaster(wp_u18, crs=proj.utm31)
slums$wp_u18 <- exact_extract(wp_u18, slums, 'sum')
rm(wp_u18)

#WP UN-adj unconstrained (2015)
wp_u15 <- raster("wp_u/nga_ppp_2015_UNadj.tif")
wp_u15 <- crop(wp_u15, bbox)
wp_u15 <- projectRaster(wp_u15, crs=proj.utm31)
slums$wp_u15 <- exact_extract(wp_u15, slums, 'sum')
rm(wp_u15)

#WP UN-adj constrained (2020)
wp_c20 <- raster("wp_c/nga_ppp_2020_UNadj_constrained.tif")
wp_c20 <- crop(wp_c20, bbox)
wp_c20 <- projectRaster(wp_c20, crs=proj.utm31)
slums$wp_c20 <- exact_extract(wp_c20, slums, 'sum')
rm(wp_c20)

#WPE (2016)
wpe16 <- raster("wpe/wpe_31n_2.tif")
wpe16 <- projectRaster(wpe16, crs=proj.utm31)
slums$wpe16 <- exact_extract(wpe16, slums, 'sum')
rm(wpe16)

#LandScan (2018)
ls18 <- raster("landscan/LandScan Global 2018/lspop2018/w001001.adf")
ls18 <- crop(ls18, bbox)
#writeRaster(ls18, "landscan/LandScan Global 2018/ls18.tif") # for visuals in ArcGIS
ls18 <- projectRaster(ls18, crs=proj.utm31)
slums$ls18 <- exact_extract(ls18, slums, 'sum')
rm(ls18)

#LandScan (2015)
ls15 <- raster("landscan/LandScan Global 2015/lspop2015/w001001.adf")
ls15 <- crop(ls15, bbox)
ls15 <- projectRaster(ls15, crs=proj.utm31)
slums$ls15 <- exact_extract(ls15, slums, 'sum')
rm(ls15)

#GRID3 Nigeria v1.2 (2016)
grid16 <- raster("grid3/NGA - population - v1.2 - gridded/NGA_population_v1_2_gridded.tif")
grid16 <- crop(grid16, bbox)
grid16 <- projectRaster(grid16, crs=proj.utm31)
slums$grid16 <- exact_extract(grid16, slums, 'sum')
rm(grid16)


######################
## Step 7. Calculate percent slums
######################

t_slums <- as.data.frame(slums)
t_slums <- transpose(t_slums[,c(-1,-2,-3)])
c_names <- c("Non_Slum", "Slum")
r_names <- colnames(slums[,c(-1,-2,-3)])
r_names <- r_names[1:12]
colnames(t_slums) <- c_names
t_slums <- cbind(r_names, t_slums)
t_slums$Total <- t_slums$Non_Slum + t_slums$Slum
t_slums$Percent <- t_slums$Slum / t_slums$Total * 100


######################
## Step 8. Calculate DHS slum households for reference
######################

# 2018 Nigeria DHS
dhs <- read_dta("D:/Dropbox/Work_main/p_unhab_slum_indices/data/NGA/DHS/NGHR7ADT/NGHR7AFL.DTA")
keep <- c("hv001", "hv002", "hv005", "hv009", "hv023", "hv024", "hv025", "hv201", "hv205", "hv213", "hv216")
dhs <- dhs[keep]
dhs <- dhs[dhs$hv023==65 | dhs$hv023==66,] # these are "Lagos urban" and "Lagos rural" households

#flag slum households
dhs$water <- ifelse(dhs$hv201 %in% c(32, 42, 43, 96), 1, 0)
dhs$toilet <- ifelse(dhs$hv205  %in% c(14, 23, 31, 42, 43, 96), 1, 0)
dhs$floor <- ifelse(dhs$hv009 %in% c(10:29, 41, 96, 99), 1, 0)
dhs$crowd <- ifelse(dhs$hv009 / dhs$hv216 >=3, 1, 0)
dhs$slumhh <- ifelse(dhs$water + dhs$toilet + dhs$floor + dhs$crowd >=1, 1, 0)

#calculate weighted percent slumhh

dhs$aveslumhh <- sum(dhs$slumhh * (dhs$hv005/1000000)) / sum(dhs$hv005/1000000) *100





