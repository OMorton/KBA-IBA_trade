
library(terra)
library(tidyterra)
library(tidyverse)

terraOptions(tempdir = "I:/TERRA", memfrac = .8)


#### All SR ####
SR_rast <- rast("F:/Data/Spatial/AOH_steps/All_SR/All_SR_ResBrNbr_NOVupd_zero_to_NA.tif")
#plot(SR_rast)

sr_samp <- spatSample(SR_rast, 10^7, "regular")
q95_sr <- quantile(sr_samp, probs=c(0.95), na.rm=TRUE)
q75_sr <- quantile(sr_samp, probs=c(0.75), na.rm=TRUE)

m <- c(0, q75_sr[[1]], NA,
       q75_sr[[1]], q95_sr[[1]], 1,
       q95_sr[[1]], 1000, 2)
pclmat <- matrix(m, ncol=3, byrow=TRUE)
SR_rast_p <- classify(SR_rast, pclmat, include.lowest=TRUE,datatype = "INT1U",
                      filename="F:/Data/Spatial/Processed_AOH/TopProportion/SR_All_ResBrNbr_5_25top_NA_NOVupd.tif")



#### All threatscore ####
thr_rast <- rast("F:/Data/Spatial/AOH_steps/All_ThrScore/Total_ResBrNbr_ThreatScore_NOVupd_zero_to_NA.tif")

#plot(thr_rast)

thr_samp <- spatSample(thr_rast, 10^7, "regular")
q95_thr <- quantile(thr_samp, probs=c(0.95), na.rm=TRUE)
q75_thr <- quantile(thr_samp, probs=c(0.75), na.rm=TRUE)

m <- c(0, q75_thr[[1]], NA,
       q75_thr[[1]], q95_thr[[1]], 1,
       q95_thr[[1]], 100, 2)
pclmat <- matrix(m, ncol=3, byrow=TRUE)
thr_rast_p <- classify(thr_rast, pclmat, include.lowest=TRUE,datatype = "INT1U",
                         filename="F:/Data/Spatial/Processed_AOH/TopProportion/Threatscore_All_ResBrNbr_5_25top_NA_NOVupd.tif")


#### All threatscore ALT ####
thr_rast <- rast("F:/Data/Spatial/AOH_steps/All_ThrScore/Total_ResBrNbr_ThreatScore_ALT_NOVupd_zero_to_NA.tif")

#plot(thr_rast)

thr_samp <- spatSample(thr_rast, 10^7, "regular")
q95_thr <- quantile(thr_samp, probs=c(0.95), na.rm=TRUE)
q75_thr <- quantile(thr_samp, probs=c(0.75), na.rm=TRUE)

m <- c(0, q75_thr[[1]], NA,
       q75_thr[[1]], q95_thr[[1]], 1,
       q95_thr[[1]], 100, 2)
pclmat <- matrix(m, ncol=3, byrow=TRUE)
thr_rast_p <- classify(thr_rast, pclmat, include.lowest=TRUE,datatype = "INT1U",
                       filename="F:/Data/Spatial/Processed_AOH/TopProportion/Threatscore_ALT_All_ResBrNbr_5_25top_NA_NOVupd.tif")


#### All rw ####
rw_rast <- rast("F:/Data/Spatial/AOH_RW_steps/All_RW/All_RW_ResBrNbr_Ave_NOVupd.tif")

#plot(rw_rast)

rw_samp <- spatSample(rw_rast, 10^7, "regular")
q95_rw <- quantile(rw_samp, probs=c(0.95), na.rm=TRUE)
q75_rw <- quantile(rw_samp, probs=c(0.75), na.rm=TRUE)

m <- c(0, q75_rw[[1]], NA,
       q75_rw[[1]], q95_rw[[1]], 1,
       q95_rw[[1]], 100, 2)
pclmat <- matrix(m, ncol=3, byrow=TRUE)
rw_rast_p <- classify(rw_rast, pclmat, include.lowest=TRUE,datatype = "INT1U",
                         filename="F:/Data/Spatial/Processed_AOH/TopProportion/RW_All_ResBrNbr_5_25top_NA_NOVupd.tif")


#### Top Areas in KBAs ####

KBA_spat <- vect("F:/Data/Spatial/AOH-KBA/kbas_minus_nine_invalid.shp")

SR_all <- rast("F:/Data/Spatial/Processed_AOH/TopProportion/SR_All_ResBrNbr_5_25top_NA_NOVupd.tif")
RW_all <- rast("F:/Data/Spatial/Processed_AOH/TopProportion/RW_All_ResBrNbr_5_25top_NA_NOVupd.tif")
Thr_all <- rast("F:/Data/Spatial/Processed_AOH/TopProportion/Threatscore_All_ResBrNbr_5_25top_NA_NOVupd.tif")
Thr_ALT <- rast("F:/Data/Spatial/Processed_AOH/TopProportion/Threatscore_ALT_All_ResBrNbr_5_25top_NA_NOVupd.tif")

AOH_crs <- crs(SR_all)
kba_reproj <- sf::st_transform(sf::st_as_sf(KBA_spat), crs = AOH_crs)

## SR
SR_All_kba <- mask(SR_all, kba_reproj, 
                   filename = "F:/Data/Spatial/Processed_AOH/KBA_mask/SR_All_ResBrNbr_5_25top_NA_kba_NOVupd.tif",
                   datatype = "INT1U")

## RW 
RW_All_kba <- mask(RW_all, kba_reproj, 
                   filename = "F:/Data/Spatial/Processed_AOH/KBA_mask/RW_All_ResBrNbr_5_25top_NA_kba_NOVupd.tif",
                   datatype = "INT1U")


## Thr
Thr_All_kba <- mask(Thr_all, kba_reproj, 
                    filename = "F:/Data/Spatial/Processed_AOH/KBA_mask/Thr_All_ResBrNbr_5_25top_NA_kba_NOVupd.tif",
                   datatype = "INT1U")

## Thr ALT
Thr_All_kba <- mask(Thr_ALT, kba_reproj, 
                    filename = "F:/Data/Spatial/Processed_AOH/KBA_mask/Thr_ALT_All_ResBrNbr_5_25top_NA_kba_NOVupd.tif",
                    datatype = "INT1U")

kba_list <- list.files("F:/Data/Spatial/Processed_AOH/KBA_mask", full.names = TRUE)
all_list <- list.files("F:/Data/Spatial/Processed_AOH/TopProportion/", full.names = TRUE)

full_ls <- c(kba_list[grepl("NA",kba_list)& grepl("_5_",kba_list) & 
                        grepl("All",kba_list) & grepl("ResBrNbr",kba_list) &
                        grepl("NOVupd",kba_list) &
                        !grepl("ALT", kba_list)],
              # all_list[grepl("NA",all_list)& grepl("_5_",all_list)& 
              #            grepl("All",all_list) & grepl("ResBrNbr",all_list) &
              #            grepl("NOVupd",kba_list) &
              #            !grepl("ALT", all_list)]
             "F:/Data/Spatial/Processed_AOH/TopProportion/RW_All_ResBrNbr_5_25top_NA_NOVupd.tif",
             "F:/Data/Spatial/Processed_AOH/TopProportion/Threatscore_All_ResBrNbr_5_25top_NA_NOVupd.tif",
             "F:/Data/Spatial/Processed_AOH/TopProportion/SR_All_ResBrNbr_5_25top_NA_NOVupd.tif" )

freq_df <- data.frame()

for (i in 1:length(full_ls)) {
 print(i)
   pth <- full_ls[[i]]
  
  if (grepl("_kba",pth)) {type <- "KBA"} else {type <- "Global"}
  if (grepl("RW_",pth)) {m <- "RW"}
  if (grepl("SR_",pth)) {m <- "SR"}
  if (grepl("Thr",pth)) {m <- "ThrSc"}
  r <- rast(full_ls[[i]]) %>% freq() %>% as.data.frame() %>%
    mutate(scale = type, metric = m, type = "All")
  
  freq_df <- rbind(freq_df, r)
}

write.csv(freq_df, "Data/KBA_AOH/Hotspot_overlay_ResBrNbr.csv")

