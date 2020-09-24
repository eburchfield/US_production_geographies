source("BS_func.R")
source("BS_load.R")

####################################################################################################################################
#  Make LRR file
####################################################################################################################################

# # Load panels built in WF/Diversity
# panel_list <- Sys.glob(paste0("C:/Users/eburchf/OneDrive - Emory University/WF/Diversity/out/panels/*.RDS"))
# null <- readRDS(panel_list[1])
# 
# # Load LRR data built in WF/RR_Bright_spots_corn
# library(rgdal)
# library(maptools)
# library(raster)
# mlra <- readOGR("C:/Users/eburchf/OneDrive - Emory University/WF/RR_Bright_spots_corn/Bright_spots/data/eco_region/mlra_v42.shp") %>%
#   arrange(LRR_NAME)
# county_sp <- as(county, "Spatial")
# lrr <- unionSpatialPolygons(SpP = mlra, IDs = mlra$LRR_NAME)
# lrr_df <- data.frame(NAME = row.names(lrr),
#                      ID = 1:29, 
#                      row.names = row.names(lrr))
# lrr <- SpatialPolygonsDataFrame(lrr, lrr_df)
# county_sp <- spTransform(county_sp, mlra@proj4string)
# # lrr_crop <- crop(lrr, county_sp)
# 
# # region assigned based on where the county centroid located
# center <- gCentroid(county_sp, byid=T, id = county$GEOID)
# center$GEOID <- as.factor(row.names(center))
# contents <- over(lrr, center, returnList = T)
# 
# 
# for (i in 1:length(contents)) {
#   null$LRR[null$GEOID %in% contents[[i]]$GEOID] <- lrr$ID[i] 
# }
# 
# fix <- unique(null$GEOID[is.na(null$LRR)])
# null$LRR[null$GEOID == "37133"] <- 2
# null$LRR[null$GEOID == "37137"] <- 2
# null$LRR[null$GEOID == "12005"] <- 23
# null$LRR[null$GEOID == "26083"] <- 18
# null$LRR[null$GEOID == "48007"] <- 2
# null$LRR[null$GEOID == "53029"] <- 19
# null$LRR[null$GEOID == "53035"] <- 19
# null$LRR[null$GEOID == "44005"] <- 14
# null$LRR <- as.factor(null$LRR)
# 
# lrr <- null %>% dplyr::select(GEOID, LRR)
# lrr <- distinct(lrr)
# saveRDS(lrr, "./out/lrr.RDS")

# frr_shp <- county %>% group_by(FRR) %>% summarize(N = n())
# saveRDS(frr_shp, "./out/frr_shp.RDS")


####################################################################################################################################
#  Prepare panels
####################################################################################################################################

# # Load panels built in WF/Diversity
# panel_list <- Sys.glob(paste0("C:/Users/eburchf/OneDrive - Emory University/WF/Diversity/out/panels/*.RDS"))
# lrr <- readRDS("./out/lrr.RDS")
# county_df <- county
# st_geometry(county_df) <- NULL
# 
# for (i in 1:length(panel_list)) {
#   print(i)
#   original <- readRDS(panel_list[i])
#   df <- merge(original, lrr, by = "GEOID", all = T)
#   df <- merge(df, county_df %>% select(GEOID), by = "GEOID", all = T)
#   df <- df %>% select(-c(PHASE1, T_USDA_TEX_CLASS))
#   saveRDS(df, paste0("./out/", df$CROP[1], ".RDS"))
# }

####################################################################################################################################
#  Add fertilizer data to panels (8/11/20)
####################################################################################################################################

