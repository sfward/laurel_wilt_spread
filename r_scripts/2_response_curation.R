# Load packages  --------------------------------------------
#
#
#---
source("r_scripts/packages.R")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


# Load in master data LWD distribution data --------------------------------------------
#
#
#---
# from Lynn: If the lw_detection_year is blank, it hasn't been found there yet. 
LWD_cty <- read_excel(path="data/LWD_county_level/Laurel_Wilt_Counties_2_9_2022.xlsx",
                      sheet = "Laurel_Wilt_Counties_2_9_2022",
                      range = "A1:H1264", col_names=T)
table(nchar(LWD_cty$FIPS))
LWD_cty$FIPS <- ifelse(nchar(LWD_cty$FIPS) == 4, paste("0", LWD_cty$FIPS, sep=""), paste(LWD_cty$FIPS))
table(nchar(LWD_cty$FIPS))

summary(LWD_cty)
hist(LWD_cty$lw_detection_year)

# remove 2022 data
LWD_cty <- LWD_cty[which(LWD_cty$lw_detection_year <= 2021),]
hist(LWD_cty$lw_detection_year)



# Load in county shapefile --------------------------------------------
#
#
#---
uscountiesall <- st_read(dsn="data/human_population_density", layer="County")
uscounties <- uscountiesall %>% filter(State %!in% c("Alaska", "Hawaii", "Puerto Rico"))
#plot(st_geometry(uscounties))
table(uscounties$State)
uscounties$FIPS <- uscounties$GEOID
uscounties$tot.pop <- uscounties$B01001_001
uscounties$popden.km <- uscounties$B01001_cal
uscounties <- uscounties %>% dplyr::select(FIPS, State, NAME, tot.pop, popden.km)


# Loop assigning invasion year to county --------------------------------------------
#
#
#---
invaded_counties <- LWD_cty[!is.na(LWD_cty$lw_detection_year),]
nrow(invaded_counties)
invaded_counties
i <- 1
uscounties$Invaded <- NULL
uscounties$YrInv <- NULL
for(i in 1:nrow(uscounties)){
  curr_fip <- uscounties$FIPS[i]
  curr_fip_inv_status <- LWD_cty[which(LWD_cty$FIPS %in% curr_fip),"lw_detection_year"]
  curr_fip_inv_status <- ifelse(nrow(curr_fip_inv_status)==0,NA,curr_fip_inv_status)
  uscounties$Invaded[i] <- ifelse(is.na(curr_fip_inv_status), "no", "yes")
  uscounties$YrInv[i] <- curr_fip_inv_status
}
uscounties$YrInv <- unlist(uscounties$YrInv)
summary(uscounties$YrInv)
table(uscounties$YrInv)
sum(table(uscounties$YrInv))
#plot(uscounties["YrInv"])

uscounties_ALBERS <- st_transform(uscounties, CRS("+init=epsg:5070"))
#uscounties_ALBERS$centroids <- st_centroid(st_geometry(uscounties_ALBERS))

st_write(uscounties_ALBERS, "shapefiles/invasion_data/uscounties_ALBERS.shp", delete_layer =T)


# Loop creating dataset for COXPH model --------------------------------------------
#
#
#---
# create data frame to be populated
lwd_start_inv <- range(na.omit(uscounties$YrInv))[1]
lwd_last_obs <- 2021
lwd_time_series_length <- length(lwd_start_inv:lwd_last_obs)

lwd_df <- as.data.frame(uscounties)
# repeat each row so that it can have one row for each year in the invasion time series
lwd_df_CPH <-   lwd_df[rep(seq_len(nrow(lwd_df)), each = lwd_time_series_length), ]
nrow(lwd_df_CPH) == nrow(lwd_df)*lwd_time_series_length
lwd_counties <- unique(lwd_df_CPH$FIPS)
length(lwd_counties)==nrow(lwd_df)

table(lwd_df_CPH$YrInv)


lwd_pb = txtProgressBar(min = 0, max = length(lwd_counties), initial = 0) 

for(i in lwd_counties){
  #i <- lwd_df_CPH$FIPS[1]
  if(i == lwd_df_CPH$FIPS[1]){lwd_counter <- 1} else {lwd_counter <- lwd_counter +1}
  
  lwd_curr_county_df <- lwd_df_CPH[which(lwd_df_CPH$FIPS %in% i),]
  
  # create rows for CPH model
  lwd_curr_county_df$time0 <- lwd_start_inv:lwd_last_obs
  lwd_curr_county_df$time1 <- (lwd_start_inv+1):(lwd_last_obs+1)
  
  # get year of initial invasion at point AQUI
  lwd_year_of_invasion <- unique(ifelse(is.na(lwd_curr_county_df$YrInv), 2022,lwd_curr_county_df$YrInv))
  
  # assign 0s and 1s (if invaded) based on year of initial invasion
  lwd_curr_county_df$invaded_time1 <- ifelse(lwd_curr_county_df$time1 <= lwd_year_of_invasion, 0, 1)

  # run for loop to assess propagule pressure in each year for each point
  lwd_curr_county_df$prop_press_dist <- NA
  lwd_curr_county_df$prop_press_dist2 <- NA
  #
  for(r in 1:nrow(lwd_curr_county_df)){
    # r <- 13
    lwd_curr_year <- lwd_curr_county_df[r,"time0"]
    
    # get all previous invaded points, and calculate distance to such points
    lwd_curr_FIPS <- uscounties_ALBERS[which(uscounties_ALBERS$FIPS %in% lwd_curr_county_df$FIPS),]
    lwd_previously_invaded_points <- uscounties_ALBERS[which(uscounties_ALBERS$YrInv < lwd_curr_year),]
    
    suppressWarnings(
      lwd_dists_to_all_points <- st_distance(st_centroid(lwd_curr_FIPS), st_centroid(lwd_previously_invaded_points))
    )
    #
    lwd_dists_to_all_points <- lwd_dists_to_all_points[which(as.numeric(lwd_dists_to_all_points)>0)]/1000
    
    #
    lwd_curr_county_df$prop_press_dist[r] <- sum(1/(lwd_dists_to_all_points))
    lwd_curr_county_df$prop_press_dist2[r] <- sum(1/(lwd_dists_to_all_points^2))
  }
  
  #
  invasion_row  <- which(lwd_curr_county_df$invaded_time1 == 1)[1]
  if(is.na(invasion_row)==F){
    lwd_curr_county_df <- lwd_curr_county_df[1:invasion_row,]}
  
  
  if(i == lwd_df_CPH$FIPS[1]){lwd_CPH_fin <- lwd_curr_county_df} else {
    lwd_CPH_fin <- rbind.data.frame(lwd_CPH_fin,lwd_curr_county_df)}
  
  setTxtProgressBar(lwd_pb, lwd_counter)
}
summary(lwd_CPH_fin)
nrow(lwd_CPH_fin)

lwd_CPH_fin$geometry <- NULL
fwrite(lwd_CPH_fin, "data/CPH_data/CPH_LWD.csv", na=".")

