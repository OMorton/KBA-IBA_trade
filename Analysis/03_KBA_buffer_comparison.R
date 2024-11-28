
library(sf)
library(terra)
library(tidyverse)
library(data.table)
library(dqrng)


options(scipen = 999)
terraOptions(tempdir = "I:/TERRA", memfrac = .8)
source("Functions.r")

#### Data ####
## metrics
SR_rast <- rast("F:/Data/Spatial/AOH_steps/All_SR/All_SR_ResBrNbr_NOVupd_zero_to_NA.tif")
Thr_rast <- rast("F:/Data/Spatial/AOH_steps/All_ThrScore/Total_ResBrNbr_ThreatScore_NOVupd_zero_to_NA.tif")
RW_rast <- rast("F:/Data/Spatial/AOH_RW_steps/All_RW/All_RW_ResBrNbr_Ave_NOVupd.tif")
ref_crs <- crs(SR_rast)

## Countries
world <- sf::st_read("Data/GADMworld/gadm_410-levels.gpkg", layer = "ADM_0")
# Europe <- rnaturalearth::ne_countries(continent = "Europe", returnclass = "sf") %>% st_transform(crs = ref_crs)
# Asia <- rnaturalearth::ne_countries(continent = "Asia", returnclass = "sf") %>% st_transform(crs = ref_crs)
# Africa <- rnaturalearth::ne_countries(continent = "Africa", returnclass = "sf") %>% st_transform(crs = ref_crs)
# Oceania <- rnaturalearth::ne_countries(continent = "Oceania", returnclass = "sf") %>% st_transform(crs = ref_crs)
# North_America <- rnaturalearth::ne_countries(continent = "North America", returnclass = "sf") %>% st_transform(crs = ref_crs)
# South_America <- rnaturalearth::ne_countries(continent = "South America", returnclass = "sf") %>% st_transform(crs = ref_crs)
# 
# plot(select(Europe, geometry))
# Europe_SR <- crop(SRNBr_rast, ext(Europe), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_SR_nafix_Europe.tif")
# Asia_SR <- crop(SRNBr_rast, ext(Asia), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_SR_nafix_Asia.tif")
# Africa_SR <- crop(SRNBr_rast, ext(Africa), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_SR_nafix_Africa.tif")
# Oceania_SR <- crop(SRNBr_rast, ext(Oceania), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_SR_nafix_Oceania.tif")
# NAm_SR <- crop(SRNBr_rast, ext(North_America), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_SR_nafix_NorthAmerica.tif")
# SAm_SR <- crop(SRNBr_rast, ext(South_America), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_SR_nafix_SouthAmerica.tif")
# 
# Europe_Thr <- crop(ThrNBr_rast, ext(Europe), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_Thr_nafix_Europe.tif")
# Asia_Thr <- crop(ThrNBr_rast, ext(Asia), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_Thr_nafix_Asia.tif")
# Africa_Thr <- crop(ThrNBr_rast, ext(Africa), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_Thr_nafix_Africa.tif")
# Oceania_Thr <- crop(ThrNBr_rast, ext(Oceania), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_Thr_nafix_Oceania.tif")
# NAm_Thr <- crop(ThrNBr_rast, ext(North_America), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_Thr_nafix_NorthAmerica.tif")
# SAm_Thr <- crop(ThrNBr_rast, ext(South_America), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_Thr_nafix_SouthAmerica.tif")
# 
# Europe_RW <- crop(RWNBr_rast, ext(Europe), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_RW_nafix_Europe.tif")
# Asia_RW <- crop(RWNBr_rast, ext(Asia), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_RW_nafix_Asia.tif")
# Africa_RW <- crop(RWNBr_rast, ext(Africa), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_RW_nafix_Africa.tif")
# Oceania_RW <- crop(RWNBr_rast, ext(Oceania), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_RW_nafix_Oceania.tif")
# NAm_RW <- crop(RWNBr_rast, ext(North_America), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_RW_nafix_NorthAmerica.tif")
# SAm_RW <- crop(RWNBr_rast, ext(South_America), filename = "F:/Data/KBA_paper/Continent_Metrics/Total_Nonbreeding_RW_nafix_SouthAmerica.tif")
# 
# w2 <- rnaturalearth::ne_countries(returnclass = "sf", scale = 10)
# wc <- w2 %>% select(iso_a3, adm0_a3, continent, - geometry) %>% as.data.frame()
# check1 <- k %>% left_join(wc, by = c("ISO3" = "adm0_a3")) %>% filter(is.na(continent))

## KBAs and buffers
KBA_only <- sf::st_read("Data/KBA/kbas_minus_nine_invalid_TerrOnly.shp")
KBA_forbuff <- sf::st_read("Data/KBA/kbas_minus_nine_invalid_TerrOnly_10km_s500PT_v2.gpkg")

#### Preparation ####
## check GID0 match and coverage
KBA_only <- KBA_only %>% mutate(ISO3 = case_when(ISO3 == "MAC" ~ "CHN",
                                                 ISO3 == "HKG" ~ "CHN",
                                                 .default = ISO3))
KBA_10kmbuff <- KBA_forbuff %>% mutate(ISO3 = case_when(ISO3 == "MAC" ~ "CHN",
                                                 ISO3 == "HKG" ~ "CHN",
                                                 .default = ISO3))
kba_df <- KBA_only %>% as.data.frame()

k <- data.frame(ISO3 = unique(kba_df$ISO3)) # 237 GID0 have KBAs
world_df <- world %>% as.data.frame()
w <- data.frame(GID0 = unique(world_df$GID_0),check = 1) # 263 GID0 in world
check <- left_join(k, w, by = c("ISO3" = "GID0")) ## all 237 match (HKG MAC = CHN)

## Standardise CRS 

world_reproj <- st_transform(world, crs = ref_crs) %>% vect()
kba_reproj <- st_transform(KBA_only, crs = ref_crs) %>% vect()
KBA_10kmbuff_reproj <- st_transform(KBA_10kmbuff, crs = ref_crs) %>% vect()


#### Country KBA + buffs metric files ####
country.list <- unique(kba_df$ISO3)

start.time <- Sys.time()
#i <- 1
#m <- "RW"
for(m in c("SR", "RW", "Thr")) {
  if (m=="SR") {metric.lyr <- SR_rast}
  if (m=="RW") {metric.lyr <- RW_rast}
  if (m=="Thr") {metric.lyr <- Thr_rast}
  
  for (i in 1:length(country.list)) {
  
  start.time <- Sys.time()
  targ.c <- country.list[i] #get country
  cat("Working on ", targ.c, ": ", i, "out of", length(country.list), '\n')
  
  c.vect <- subset(world_reproj, world_reproj$GID_0 == targ.c) %>% st_as_sf()# country shape
  kba.vect <- subset(kba_reproj, kba_reproj$ISO3 == targ.c) # countries KBAs
  kba.10kmbuff.vect <- subset(KBA_10kmbuff_reproj, KBA_10kmbuff_reproj$ISO3 == targ.c)# countries KBAs+10kmbuff
  
  # cat("Buffering", '\n')
  # get buffer in loop
  # kba10km.vect <- st_buffer(st_as_sf(kba.forbuff.vect), dist = 10000)
  # buff.time <- Sys.time()
  
  # cat("Intersecting", '\n')
  # 
  # ## crop to just KBAs area in the country
  # overlap.c.kba <-  kba.vect %>% st_as_sf() %>%
  #   st_filter(y=c.vect, .predicate = st_intersects)
  # kba.c <- st_intersection(overlap.c.kba, c.vect) %>%
  #   st_make_valid() %>% vect()
  # 
  # ## crop to just KBA+10kmbiff area in the country
  # overlap.c.kba10km <-  kba10km.vect %>% st_as_sf() %>%
  #   st_filter(y=c.vect, .predicate = st_intersects)
  # kba10km.c <- st_intersection(overlap.c.kba10km, c.vect) %>%
  #   st_make_valid() %>% vect()
  # 
  # intersect.time <- Sys.time()
  cat("Masking", '\n')
  
  ## get metric values in kba and kba+10kmbuff
  metric.kba.c <- metric.lyr %>% crop(kba.vect) %>% mask(kba.vect, overwrite = TRUE,
                                                         filename = paste0("F:/Data/KBA_paper/Country_masks/ResBrNbr/", m,"/",targ.c, ".KBAonly.", m, ".tif" ))
  metric.kba10km.c <- metric.lyr %>% crop(kba.10kmbuff.vect) %>% mask(kba.10kmbuff.vect, overwrite = TRUE,
                                                                 filename = paste0("F:/Data/KBA_paper/Country_masks/ResBrNbr/",m,"/",targ.c, ".KBA10kmbuff.", m, ".tif" ))
  
  }
}
end.time <- Sys.time()


#### GET OBS and NULL DIFF - Per country ####
  
## empty data frames for storage. Pre-allocate the null loop space for speed
obs.res.out <- data.frame(GID_0 = character(), obs.kba.mean = numeric(), 
                          obs.10kmbuff.mean = numeric(), obs.diff = numeric())

null.res.out <- data.frame(GID_0 = rep(NA,1000), null.mn.samp = rep(NA,1000), 
                           null.mn.outsamp = rep(NA,1000), null.diff = rep(NA,1000),
                           kba.cells = rep(NA,1000), buff10km.cells = rep(NA,1000))

null.res.all <- data.frame(GID_0 = character(), null.mn.samp = numeric(), 
                           null.mn.outsamp = numeric(), null.diff = numeric(),
                           kba.cells = numeric(), buff10km.cells = numeric())  

## loop iterates per metric, per country, then within each metric per country
## samples that countries KBA area + 10km buffer the same number of cells as the
## KBA area in the country. Then calculates the mean difference between the sampled and
## not sampled area to create a null distribution of differences between KBA sized samples 
## and the not sampled.
## Also calculates the observed difference
m <- "SR"
targ.c <- "RUS"

for(m in c("SR", "Thr")) {
  for (i in 1:length(country.list)) {
    
    targ.c <- country.list[i] #get country
      cat("Reading in ", targ.c, '\n')
    
    metric.kba.c <- rast(paste0("F:/Data/KBA_paper/Country_masks/ResBrNbr/",m,"/",targ.c, ".KBAonly.", m, ".tif" ))
    metric.kba10km.c <- rast(paste0("F:/Data/KBA_paper/Country_masks/ResBrNbr/",m,"/",targ.c, ".KBA10kmbuff.", m, ".tif" ))
  
      cat("Converting ", targ.c, '\n')
   
      metric.kba.c.df <- as.data.table.raster(metric.kba.c, na.rm = TRUE)
   metric.kba10km.c.df <- as.data.table.raster(metric.kba10km.c, na.rm = TRUE)
  
  kba.ncells.c <- nrow(metric.kba.c.df) ## get the number of KBA cells
  kbabuff10km.ncells.c <- nrow(metric.kba10km.c.df) ## get the number of KBA+10kmbuff cells
  buff10km.ncells.c <- kbabuff10km.ncells.c - kba.ncells.c## get the number of 10kmbuff cells
  
  obs.kba.mn.c <- mean(metric.kba.c.df[[1]])
  obs.10kmbuff.mn.c <- (sum(metric.kba10km.c.df[[1]]) - sum(metric.kba.c.df[[1]]))/buff10km.ncells.c
  
  obs.diff.c <- obs.kba.mn.c - obs.10kmbuff.mn.c

  ## results 
  res.c <- data.frame(GID_0 = targ.c, obs.kba.mean = obs.kba.mn.c, 
                      obs.10kmbuff.mean = obs.10kmbuff.mn.c, obs.diff = obs.diff.c)
  
  obs.res.out <- rbind(obs.res.out, res.c)
  
  write.csv(obs.res.out, paste0("Data/Randomization/ResBrNbr/Obs/", m,"_All_10kmbuff_obs.csv"))
  
  ## GET NULL DIFF DISTR ##
  sum_kba10km <- sum(metric.kba10km.c.df[[1]])
  
  #j <- 1
  for (j in 1:1000) {
    cat(i, "out of", length(country.list), "Samp:", j , '\n')
    
    samp.c <- metric.kba10km.c.df[dqsample(kbabuff10km.ncells.c ,kba.ncells.c, replace = FALSE)]
    
    null.mn.samp.c <- mean(samp.c[[1]])
    null.mn.outsamp.c <- (sum_kba10km - sum(samp.c))/buff10km.ncells.c
    
    null.diff.c <- null.mn.samp.c - null.mn.outsamp.c
    
     null.c <- data.frame(GID_0 = targ.c, null.mn.samp = null.mn.samp.c, 
                          null.mn.outsamp = null.mn.outsamp.c, null.diff = null.diff.c,
                          kba.cells = kba.ncells.c, buff10km.cells = buff10km.ncells.c)
    
    #null.res.out <- rbind(null.res.out, null.c)
     null.res.out[j,] <- null.c
     
  }
  null.res.all <- rbind(null.res.all, null.res.out)
  write.csv(null.res.all, paste0("Data/Randomization/ResBrNbr/Nulls/", m, "All_10kmbuff_null.csv"))
  metric.kba10km.c.df <- NULL
  metric.kba10km.c.df <- NULL
  samp.c <- NULL
  null.res.out <- data.frame(GID_0 = character(), null.mn.samp = numeric(), 
                             null.mn.outsamp = numeric(), null.diff = numeric(),
                             kba.cells = numeric(), buff10km.cells = numeric())
  
  }
  obs.res.out <- data.frame(GID_0 = character(), obs.kba.mean = numeric(), 
                            obs.10kmbuff.mean = numeric(), obs.diff = numeric())
  
  null.res.out <- data.frame(GID_0 = rep(NA,1000), null.mn.samp = rep(NA,1000), 
                             null.mn.outsamp = rep(NA,1000), null.diff = rep(NA,1000),
                             kba.cells = rep(NA,1000), buff10km.cells = rep(NA,1000))
  
  null.res.all <- data.frame(GID_0 = character(), null.mn.samp = numeric(), 
                             null.mn.outsamp = numeric(), null.diff = numeric(),
                             kba.cells = numeric(), buff10km.cells = numeric()) 
}


#### Calc significance ####

null.df <- read.csv("Data/Randomization/ResBrNbr/Nulls/SRAll_10kmbuff_null.csv")
obs.df <- read.csv("Data/Randomization/ResBrNbr/Obs/SR_All_10kmbuff_obs.csv")

## Small island and ocean states or territories.
null.df %>% filter(is.na(null.diff)) %>% distinct(GID_0)
obs.df %>% filter(is.na(obs.diff)) %>% distinct(GID_0)

null.df <- null.df %>% filter(!is.na(null.diff))
obs.df <- obs.df %>% filter(!is.na(obs.diff))

c.lst <- unique(obs.df$GID_0)
m.lst <- c("SR", "RW", "Thr")
p.df <- data.frame(metric = character(), GID_0 = character(), p.val = numeric())

for (m in 1:3) {
  m.m <- m.lst[m]
  null.df <- read.csv(paste0("Data/Randomization/ResBrNbr/Nulls/",m.m,"All_10kmbuff_null.csv"))
  obs.df <- read.csv(paste0("Data/Randomization/ResBrNbr/Obs/",m.m,"_All_10kmbuff_obs.csv"))
  null.df <- null.df %>% filter(!is.na(null.diff))
  obs.df <- obs.df %>% filter(!is.na(obs.diff))
  for (i in 1:length(c.lst)) {
    c.i <- c.lst[i]
    null.i <- filter(null.df, GID_0 == c.i)
    obs.i <- filter(obs.df, GID_0 == c.i)
    pValue <- sum(abs(null.i$null.diff) >= abs(obs.i$obs.diff)) / 1000
    out.df <- data.frame(metric= m.m, GID_0 = c.i, p.val = pValue)
    p.df <- rbind(p.df, out.df)
  }
}

write.csv(p.df, "Data/Randomization/ResBrNbr/Sig/Randomization.pvals.csv")




rnaturalearth::ne_countries(type = "countries") %>% as.data.frame()
w2 <- rnaturalearth::ne_countries(returnclass = "sf", scale = 10) %>% as.data.frame()
wc <- w2 %>% ungroup() %>% select(iso_a3, adm0_a3, continent, region_un, region_wb, subregion)
d <- sr.obs.df2 %>% left_join(wc, by = c("GID_0" = "iso_a3"))
