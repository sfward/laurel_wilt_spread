# Load packages  --------------------------------------------
#
#
#---
source("r_scripts/packages.R")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


# Load in shapefile of us counties with invasion information --------------------------------------------
#
#
#---
spread_ALBERS = st_read("shapefiles/invasion_data/uscounties_ALBERS.shp")

location_of_PortWentworth <- cbind.data.frame(long=-81.1629, lat=32.1480)
location_of_PortWentworth$NAME <- "First detection"
crs_longlat <- CRS("+init=epsg:4326") # proj4string of coords

# make the SpatialPointsDataFrame object
pts_PortWW <- SpatialPointsDataFrame(coords = location_of_PortWentworth[,c("long", "lat")],
                                     data = location_of_PortWentworth,
                                     proj4string = crs_longlat)
disc_loc_ALBERs <- spTransform(pts_PortWW, crs(spread_ALBERS))
disc_loc_sf <- st_as_sf(disc_loc_ALBERs)



# invaded counties - create buffer
spread_ALBERS_invaded <- spread_ALBERS[which(spread_ALBERS$Invaded  %in% "yes"), ]
plot(spread_ALBERS_invaded[,"YrInv"])
# 500000/1000
# 500000 == 500*1000

counties_buff_500km <- st_buffer( spread_ALBERS_invaded, 500*1000)%>%  # 500 km
  st_union() %>% # unite to a geometry object
  st_sf() # make the geometry a data frame object

plot(counties_buff_500km)
plot(spread_ALBERS[,"YrInv"], add=T)


library(viridis)
my_LWD_colors <-viridis(17)
library(RColorBrewer)
my_LWD_colors_spec <-brewer.pal(11,"Spectral")
display.brewer.pal(n = 9, name = 'YlGnBu')


table(spread_ALBERS_invaded$YrInv)
spread_ALBERS_cropped <- st_crop(spread_ALBERS, spread_ALBERS_invaded)


counties_within_buffer_pts <- st_intersection(st_centroid(spread_ALBERS), counties_buff_500km)
counties_within_buffer <- spread_ALBERS[which(spread_ALBERS$FIPS %in% counties_within_buffer_pts$FIPS),]
nrow(counties_within_buffer)
st_write(counties_within_buffer, "shapefiles/invasion_data/buffered_invasion.shp", delete_layer = TRUE) 
plot(counties_within_buffer[,"Invaded"])
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#
# number of counties infested by adjoining spread vs isolated spread
#
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
run <- "n"
if(run == "y"){

# create spatial neighborhood for getting neighbors of isolated counties and adjoining counties to isolated counties
sec_order <- poly2nb(spread_ALBERS, queen = T, row.names=spread_ALBERS$FIPS)
sec_order[[1]]

# i <- 2018
for(i in 2004:2021){
    
    invaded_in_i <- spread_ALBERS[which(spread_ALBERS$YrInv %in% i),]
    
    prev_invaded_counties <- spread_ALBERS[which(spread_ALBERS$YrInv %in% (2004:(i-1))),]
    
    if(nrow(invaded_in_i) > 0){
      # sanity check - neighbors list
      for(j in 1:nrow(invaded_in_i)){
        #j <- 1
        # get points in invaded county and in each of the previously invaded counties
        new_invaded_point <- invaded_in_i[j,"FIPS"]
        
        if(i == 2004){
          prev_invaded_points <- disc_loc_ALBERs
        }else{
          prev_invaded_points <-  prev_invaded_counties
        }
        
        prev_invaded_points <- st_as_sf(prev_invaded_points)
        
        # get distance between current county and ground zero
        pts_ground_zero <- suppressWarnings(pointDistance(st_centroid(new_invaded_point),disc_loc_sf,lonlat=F))
        # get minimum distance to ground zero
        invaded_in_i$DtoGZ_km[j] <- min(pts_ground_zero)[1]/1000 
        
        # find minimum distance
        dist_vec_centroids <- suppressWarnings(pointDistance(st_centroid(new_invaded_point),st_centroid(prev_invaded_points),lonlat=F))
        invaded_in_i$D_BDY_km[j] <- min(dist_vec_centroids)[1]/1000
        
        # get closest county
        val <- which(dist_vec_centroids == min(dist_vec_centroids)[1], arr.ind = TRUE)
        closest_county_BNDRY <-  prev_invaded_counties[val,]
        invaded_in_i$NR_CNTY[j] <- as.character(closest_county_BNDRY$FIPS)
        
        # determine whether a neighbor was invaded
        loc_in_vec <- which(spread_ALBERS$FIPS %in% invaded_in_i[j, "FIPS"]) # where is the current county located in data frame
        neighbs <- spread_ALBERS[sec_order[[loc_in_vec]],] # get current county's neighbors
        neighbs_prev <- prev_invaded_counties[which(prev_invaded_counties$FIPS %in% neighbs$FIPS),] # determine if current county had neighbors in the previous year that were infested
        
        if(nrow(neighbs_prev) == 0){invaded_in_i$STAT[j] <- "iso"} else{
          invaded_in_i$STAT[j] <- "adj"
        }}
      
      vec_variables <- c("FIPS", "D_BDY_km", "DtoGZ_km", "STAT", "YrInv", "NR_CNTY")
      
      if(i == 2004){
        invaded_counties <- as.data.frame(invaded_in_i[, paste(vec_variables)])
      } else {
        invaded_counties <- rbind.data.frame(invaded_counties,as.data.frame(invaded_in_i[, paste(vec_variables)]))
      }
    }
    
    cat("Year", i,"out of", 2021,"\n") 
  }
  
  
  # beep(2)
  spread_ALBERS$STAT <- NA
  spread_ALBERS$D_BDY_km <- NA
  spread_ALBERS$DtoGZ_km <- NA
  spread_ALBERS$NR_CNTY <- NA
  
  for(i in 1:nrow(spread_ALBERS)){
    curr_fips <- spread_ALBERS$FIPS[i]
    if( length(which(invaded_counties$FIPS %in% curr_fips)) > 0 ){
      spread_ALBERS$STAT[i] <- invaded_counties[which(invaded_counties$FIPS %in% curr_fips), "STAT"]
      spread_ALBERS$D_BDY_km[i] <- invaded_counties[which(invaded_counties$FIPS %in% curr_fips), "D_BDY_km"]
      spread_ALBERS$DtoGZ_km[i] <- invaded_counties[which(invaded_counties$FIPS %in% curr_fips), "DtoGZ_km"]
      spread_ALBERS$NR_CNTY[i] <- invaded_counties[which(invaded_counties$FIPS %in% curr_fips), "NR_CNTY"]
    }
  }

st_write(obj=spread_ALBERS, dsn="shapefiles/county_status/contigUS_ALBERS_LWD.shp", delete_layer=T) 
} 

# D_BDY - distance to neareast invaded county in previous year
# DtoGZ - distance to ground zero
# STAT - county invaded had adjoining county invaded in previous year or isolated county
# YrInv - year of invasion according to 2002 onwards APHIS data
# NR_CNTY - the nearest invaded county 
contigUS_ALBERS_LWD <- st_read("shapefiles/county_status/contigUS_ALBERS_LWD.shp")
summary(contigUS_ALBERS_LWD)

LWD.adj.iso <- contigUS_ALBERS_LWD[which(contigUS_ALBERS_LWD$STAT %in% c("adj","iso")),]

contigUS_ALBERS_LWD$JumpInv <- ifelse(contigUS_ALBERS_LWD$STAT == "iso", "Non-contiguous", contigUS_ALBERS_LWD$STAT)
contigUS_ALBERS_LWD$JumpInv <- ifelse(contigUS_ALBERS_LWD$JumpInv == "adj", "Contiguous", contigUS_ALBERS_LWD$JumpInv)

contigUS_ALBERS_LWD$color <- NA
contigUS_ALBERS_LWD$color[contigUS_ALBERS_LWD$JumpInv == "Contiguous"] <- "gray"
contigUS_ALBERS_LWD$color[contigUS_ALBERS_LWD$JumpInv == "Non-contiguous"] <- "tomato2"
contigUS_ALBERS_LWD$color[is.na(contigUS_ALBERS_LWD$JumpInv)] <- "white"


isolated_counties <- LWD.adj.iso[which(LWD.adj.iso$STAT == "iso"),]
summary(isolated_counties)
isolated_counties[which(isolated_counties$D_BDY_km %in% min(isolated_counties$D_BDY_km)),]
plot(LWD.adj.iso[which(LWD.adj.iso$FIPS %in% c("13051","13179")),"Invaded"])

nrow(LWD.adj.iso[which(LWD.adj.iso$STAT == "iso"),])
nrow(LWD.adj.iso[which(LWD.adj.iso$STAT == "adj"),])
nrow(LWD.adj.iso[which(LWD.adj.iso$STAT %in% c("iso","adj")),])
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''












# Load in US shapefile --------------------------------------------
#
#
#---
statesUS = st_read("shapefiles/states/cb_2018_us_state_20m.shp")
statesUS_albers <- st_transform(statesUS, st_crs(spread_ALBERS_cropped))
statesUS_albers <- statesUS_albers %>% filter(NAME %!in% c("Puerto Rico", "Alaska", "Hawaii"))
statesUS_albers_crop <- st_crop(statesUS_albers, spread_ALBERS_invaded)

# create shapefile with predictor values for plotting
CPH_df_hosts <- fread("data/CPH_data/CPH_LWD_COMBINED.csv", na=".")
table(CPH_df_hosts$YrInv)
CPH_df_hosts$FIPS <- ifelse(nchar(CPH_df_hosts$FIPS)==4, paste(0,CPH_df_hosts$FIPS,sep=""), paste(CPH_df_hosts$FIPS))
table(nchar(CPH_df_hosts$FIPS))
predictors <- as_tibble(CPH_df_hosts) %>% dplyr::select(FIPS, POPESTIMATE2019, ltot, n_campgrounds, med_income, host_biomass, non_host_biomass, ISOTHERM,  MINMOTEMP, MAP, redbay_biomass, sassafras_biomass) %>% 
  group_by(FIPS) 
predictors <- predictors[match(unique(predictors$FIPS), predictors$FIPS),]
predictors_albers <- merge(spread_ALBERS,predictors,by="FIPS")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''








# Figure: invaded counties --------------------------------------------
#
#
#---
# inset map with bounding box
# https://www.r-bloggers.com/2019/12/inset-maps-with-ggplot2/
library(cowplot)
lwd_bb = st_as_sfc(st_bbox(statesUS_albers_crop))
inset <- ggplot() +
  geom_sf(data = statesUS_albers, fill="white") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data = lwd_bb, fill = NA, color = "black", size = 1.2) +
  coord_sf(datum=st_crs(spread_ALBERS_cropped))
#
invaded_area_continuous <- ggplot() +
  geom_sf(data = spread_ALBERS_cropped, aes(fill = YrInv), color="light gray") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers_crop, fill="transparent", color="black", lwd=1)+
  scale_fill_gradientn(colors = my_LWD_colors_spec, name="", na.value = "white",   guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme(legend.position=c(0.6,0.25),legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8), legend.justification = "left")+
  geom_sf(data=disc_loc_sf, color = "black", size = 3)+
  ggrepel::geom_text_repel(
    data = disc_loc_sf,
    aes(label = NAME, geometry = geometry),
    stat = "sf_coordinates",
    min.segment.length = 0,
    nudge_x = Inf,
    colour = "black",
    segment.colour = "black")+
  coord_sf(datum=st_crs(spread_ALBERS_cropped))

#
resize.win(6.85,5)
gg_inset_map1_cont = ggdraw() +
  draw_plot(invaded_area_continuous) +
  draw_plot(inset, x = 0.16, y = 0.06, width = 0.25, height = 0.25)
gg_inset_map1_cont

tiff("manuscript_spatial_drivers/figures/Fig1-UPDATED.tiff", units="in", width=6.85, height=5, res=600)
gg_inset_map1_cont
dev.off()
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''




# Percent of habitat invaded --------------------------------------------
#
#
#---
host_data <- predictors_albers
host_data$Invaded_number <- ifelse(host_data$Invaded=="yes", 1, 0)
#
counties_with_redbay <- host_data %>% filter(redbay_biomass  > 0)
mean(counties_with_redbay$Invaded_number)
counties_with_redbay %>% group_by(Invaded) %>% summarise(sum_red = sum(redbay_biomass))
10437347837/(6899172641+10437347837)
#
counties_with_sassafras <- host_data %>% filter(sassafras_biomass  > 0)
mean(counties_with_sassafras$Invaded_number)
counties_with_sassafras %>% group_by(Invaded) %>% summarise(sum_sas = sum(sassafras_biomass))
6247875028/(88064386369+6247875028)

#
counties_with_redbay %>% group_by(Invaded) %>% summarise(red_temp = mean(MINMOTEMP))
counties_with_sassafras %>% group_by(Invaded) %>% summarise(sas_temp = mean(MINMOTEMP))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''








# Figure: host distribution --------------------------------------------
#
#
#---
predictors_albers$redbay_biomass_millions <- predictors_albers$redbay_biomass/1e6
predictors_albers$sassafras_biomass_millions <- predictors_albers$sassafras_biomass/1e6
max(predictors_albers$sassafras_biomass_millions)-max(predictors_albers$redbay_biomass_millions)
hist(predictors_albers$redbay_biomass_millions, col="green")
hist(predictors_albers$sassafras_biomass_millions,add=T)


my_breaks <- c(0.1, 3, 20, 148, 1096)

# checking biomass values
predictors_albers$tons <- predictors_albers$sassafras_biomass/907.2
as.data.frame(predictors_albers[predictors_albers$State == "Pennsylvania",])

resize.win(6.85,7)
redbay_p <- ggplot() +
  geom_sf(data = predictors_albers, aes(fill = redbay_biomass_millions), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  #scale_fill_continuous(low="darkgreen", high="yellow", 
  #                     na.value="white",    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  
  scale_fill_viridis(option="viridis", direction=-1, trans="log",  breaks=my_breaks, labels=my_breaks,
                     na.value="white",    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title= bquote("Biomass (kg)")))+
  geom_sf(data=statesUS_albers, fill="transparent", color="black", lwd=1)+
  theme(legend.key.size = unit(0.3, 'cm'), legend.text = element_text(size=7), legend.justification = "left", 
        legend.title = element_blank(), legend.position=c(0.86,0.25))+
  coord_sf(datum=st_crs(spread_ALBERS_cropped))



sassafras_p <- 
  ggplot() +
  geom_sf(data = predictors_albers, aes(fill = sassafras_biomass_millions), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  #scale_fill_continuous(low="darkgreen", high="yellow",  
  #                      na.value="white",    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  scale_fill_viridis(option="viridis", direction=-1, trans="log", breaks=my_breaks, labels=my_breaks,
                     na.value="white",    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", title= bquote("Biomass (kg)")))+
  
  geom_sf(data=statesUS_albers, fill="transparent", color="black", lwd=1)+
  theme(legend.key.size = unit(0.3, 'cm'), legend.text = element_text(size=7), legend.justification = "left", 
        legend.title = element_blank(), legend.position=c(0.86,0.25))+
  coord_sf(datum=st_crs(spread_ALBERS_cropped))



resize.win(6.85,7)
tiff("manuscript_spatial_drivers/figures/Fig2-UPDATED.tiff", units="in", width=6.85, height=7, res=600)
ggarrange(redbay_p, sassafras_p, ncol = 1, nrow = 2,  align = "hv",
          labels="auto", label.x =0.135, label.y =1)
dev.off()
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# Effective range radius --------------------------------------------
#
#
#---
st_write(spread_ALBERS, "shapefiles/data_for_manuscript/spread_albers.shp", delete_layer=T)

spread_ALBERS$county_area <- st_area(spread_ALBERS)

annual_spread <- data.frame(year=2004:2021, area_invaded_km2=NA)
years.to.plot <- sort(unique(spread_ALBERS$YrInv))
i <- years.to.plot[1]
for(i in years.to.plot){
  curr.invaded <- spread_ALBERS[which(spread_ALBERS$YrInv %in% 2004:i),]
  annual_spread[which(annual_spread$year %in% i), "year"] <- i
  annual_spread[which(annual_spread$year==i), "area_invaded_km2"] <- sum(st_area(curr.invaded)/1e6)
}

annual_spread$radius <- sqrt(annual_spread$area_invaded_km2/pi)
annual_spread$radius_try2 <- sqrt(2*annual_spread$area_invaded_km2/pi)
plot(radius~year, data=annual_spread)
fit_spread1 <- lm(radius ~ year, data=annual_spread)
summary(fit_spread1)
32.14^2

# proof of concept (Koch's actual numbers)
r <- c(4398, 18965, 41480)
r_semi <- sqrt(2*r/pi)
t <- c(2004,2005,2006)
fit_Koch_POC <- lm(r_semi ~ t)
summary(fit_Koch_POC)

annual_spread_2006_onwards <- annual_spread[which(annual_spread$year>= 2006),]
summary(annual_spread_2006_onwards)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''






# Figure (one panel only): Effective range radius --------------------------------------------
#
#
#---
lm_eqn <- function(df){
  m <- lm(radius_try2 ~ year, df);
  eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}

set.seed(12)


fit_rangeexpansion <- lm(radius_try2 ~ year, data=annual_spread)
summary(fit_rangeexpansion)
round(summary(fit_rangeexpansion)$coef,2)
round(summary(fit_rangeexpansion)$coef,2)[2,3]^2


# breakpoint regression
library(segmented)
segmented.mod <- segmented(fit_rangeexpansion, control = seg.control(display = FALSE))
summary(segmented.mod)

#
annual_spread_BREAK_BEFORE <- annual_spread[which(annual_spread$year < 2010),]
fit_rangeexpansion_BEFORE <- lm(radius_try2 ~ year, data=annual_spread_BREAK_BEFORE)
summary(fit_rangeexpansion_BEFORE)
round(summary(fit_rangeexpansion_BEFORE)$coef,2)[2,3]^2

#
annual_spread_BREAK_AFTER <- annual_spread[which(annual_spread$year >= 2009),]
fit_rangeexpansion_AFTER <- lm(radius_try2 ~ year, data=annual_spread_BREAK_AFTER)
summary(fit_rangeexpansion_AFTER)
round(summary(fit_rangeexpansion_AFTER)$coef,2)[2,3]^2


#
annual_spread_Koch <- annual_spread[annual_spread$year %in% c(2004:2006),]
fit_spread1_Koch <- lm(radius_try2 ~ year, data=annual_spread_Koch)
summary(fit_spread1_Koch)
round(summary(fit_spread1_Koch)$coef,2)
round(summary(fit_spread1_Koch)$coef,2)[2,3]^2


#
annual_spread_2006_onwards <- annual_spread[which(annual_spread$year>= 2006),]
summary(annual_spread_2006_onwards)
fit_rangeexpansion_post2006 <- lm(radius_try2 ~ year, annual_spread_2006_onwards)
summary(fit_rangeexpansion_post2006)
round(summary(fit_rangeexpansion_post2006)$coef,2)
round(summary(fit_rangeexpansion_post2006)$coef,2)[2,3]^2



y_initial_ERR <- 530
color_dots <- viridis(6)[1]
fig_ERR_COMP <- ggplot(data=annual_spread, aes(x=year, y=radius_try2)) +
  ylab("Radius of invaded area (km)")+
  xlab("Year")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(2005,2022,4), limits=c(2003,2022), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,600,100), limits=c(0,600), expand = c(0,0)) +
  
  geom_point(aes(x = (year), y = radius_try2), col=color_dots, size=3)+
  stat_smooth(data = annual_spread, method = "lm", col = "black", se=F)+
  geom_text(x = 2005, y = y_initial_ERR, label = lm_eqn(annual_spread), parse = TRUE, size=3, check_overlap = TRUE, hjust = 0)+
  geom_segment(aes(x = 2004, y = y_initial_ERR, xend =  2004+0.8, yend = y_initial_ERR), colour = "black",size=0.8)+
  
  stat_smooth(data = annual_spread_Koch, method = "lm", col = "black", se=F, linetype="dashed")+
  geom_text(x = 2005, y = y_initial_ERR-55, label = lm_eqn(annual_spread_Koch), parse = TRUE, size=3,check_overlap = TRUE,  hjust = 0)+
  geom_segment(aes(x = 2004, y = y_initial_ERR-55, xend =  2004+0.8, yend = y_initial_ERR-55), colour = "black",size=0.8, linetype="dashed")+
  
  stat_smooth(data = annual_spread_2006_onwards, method = "lm", col = "dark gray", se=F, linetype="dashed")+
  geom_text(x = 2005, y = y_initial_ERR-55-55, label = lm_eqn(annual_spread_2006_onwards), parse = TRUE, size=3,check_overlap = TRUE, hjust = 0)+
  geom_segment(aes(x = 2004, y = y_initial_ERR-55-55, xend =  2004+0.8, yend = y_initial_ERR-55-55), colour = "dark gray",size=0.8, linetype="dashed")+

  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


seg_pt_col <- rev(viridis(6))[3]
fig_ERR_SEGMENT <- ggplot(data=annual_spread, aes(x=year, y=radius_try2)) +
  ylab("Radius of invaded area (km)")+
  xlab("Year")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(2005,2022,4), limits=c(2003,2022), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,600,100), limits=c(0,600), expand = c(0,0)) +
  
  geom_point(aes(x = (year), y = radius_try2), col=seg_pt_col, size=3)+
  stat_smooth(data = annual_spread_BREAK_BEFORE, method = "lm", col = "black", se=F)+
  geom_text(x = 2005, y = y_initial_ERR, label = lm_eqn(annual_spread_BREAK_BEFORE), parse = TRUE, size=3, check_overlap = TRUE, hjust = 0)+
  geom_segment(aes(x = 2004, y = y_initial_ERR, xend =  2004+0.8, yend = y_initial_ERR), colour = "black",size=0.8)+
  
  stat_smooth(data = annual_spread_BREAK_AFTER, method = "lm", col = "black", se=F, linetype="dashed")+
  geom_text(x = 2005, y = y_initial_ERR-55, label = lm_eqn(annual_spread_BREAK_AFTER), parse = TRUE, size=3,check_overlap = TRUE,  hjust = 0)+
  geom_segment(aes(x = 2004, y = y_initial_ERR-55, xend =  2004+0.8, yend = y_initial_ERR-55), colour = "black",size=0.8, linetype="dashed")+

  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



fig_ERR <- ggplot(data=annual_spread, aes(x=year, y=radius_try2)) +
  ylab("Radius of invaded area (km)")+
  xlab("Year")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(2005,2022,4), limits=c(2003,2022), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,600,100), limits=c(0,600), expand = c(0,0)) +
  
  geom_point(aes(x = (year), y = radius_try2), col=color_dots, size=3)+
  stat_smooth(data = annual_spread, method = "lm", col = "black", se=F)+
  geom_text(x = 2005, y = y_initial_ERR, label = lm_eqn(annual_spread), parse = TRUE, size=3, check_overlap = TRUE, hjust = 0)+
  #geom_segment(aes(x = 2004, y = y_initial_ERR, xend =  2004+0.8, yend = y_initial_ERR), colour = "black",size=0.8)+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Figure (one panel only): Distance regression --------------------------------------------
#
#
#---
spread_ALBERS$dist_from_disc <- as.numeric(st_distance(st_centroid(spread_ALBERS), disc_loc_sf)/1000)

#
plot(st_geometry(spread_ALBERS))
plot(st_geometry(disc_loc_sf), add=T, cex=5, pch=21)

spread_df <-as.data.frame(spread_ALBERS)
spread_df <- spread_df[!is.na(spread_df$YrInv),]
plot(dist_from_disc~YrInv, data=spread_df)
fit_spread2 <- lm(dist_from_disc ~ YrInv, data=spread_df)
summary(fit_spread2)
round(summary(fit_spread2)$coef,2)
confint(fit_spread2)
11.73^2
lm_eqn2 <- function(df){
  m <- lm(dist_from_disc ~ YrInv, spread_df);
  eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}


fig_DR <- ggplot(data=spread_df, aes(x=YrInv, y=dist_from_disc)) +
  ylab("Distance (km)")+
  xlab("Year")+ theme_bw()+
  geom_point(aes(x = YrInv, y = dist_from_disc), col=color_dots, size=2)+
  stat_smooth(data = spread_df, method = "lm", col = "black", se=F)+
  
  scale_x_continuous(breaks = seq(2005,2022,4), limits=c(2003,2022), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1500,300), limits=c(0,1500), expand = c(0,0)) +
  
  #geom_text(x = 2005, y = 1200, label = expression("Spread = 43 ± 4 km/yr (95% CI: 35-51)"), size=6,check_overlap = TRUE, hjust = 0)+
  geom_text(x = 2005, y = 1000, label = lm_eqn2(spread_df), parse = TRUE, size=3,check_overlap = TRUE, hjust = 0)+
  #geom_text(x = 2005, y = 875, label = expression(italic(F)["1,253"]*"= 144.9, "*italic(p)*" < 0.0001"),  size=6,check_overlap = TRUE, hjust = 0)+

theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black"))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''






# Boundary displacement --------------------------------------------
#
#
#---
spread.rad <- data.frame(year=2004:2021)
radii_n <- seq(0,355,5)
length(radii_n )
spread.rad[,paste("radii",radii_n, sep=".")] <- NA

# https://stackoverflow.com/questions/59328707/how-do-i-partition-a-circle-into-equal-multipolygon-slices-with-sf-and-r
st_wedge <- function(x,y,r,start,width,n=20){
  theta = seq(start, start+width, length=n)
  xarc = x + r*sin(theta)
  yarc = y + r*cos(theta)
  xc = c(x, xarc, x)
  yc = c(y, yarc, y)
  st_polygon(list(cbind(xc,yc)))   
}

st_wedges <- function(x, y, r, nsegs){
  width = (2*pi)/nsegs
  starts = (1:nsegs)*width
  polys = lapply(starts, function(s){st_wedge(x,y,r,s,width)})
  mpoly = st_cast(do.call(st_sfc, polys), "MULTIPOLYGON")
  mpoly
}

# 
wedge_length <- 1500 # km
discx <- coordinates(disc_loc_ALBERs)[1]
discy <- coordinates(disc_loc_ALBERs)[2]
wedges = st_wedges(discx,discy,wedge_length*1000,length(radii_n))
st_crs(wedges) <- crs(spread_ALBERS)
plot(wedges)
plot(wedges[1], add=T, col="red")
wedges <- wedges[c(length(radii_n),1:length(radii_n)-1)]
plot(wedges[1], add=T, col="blue")

i <- 2012
for(i in years.to.plot){
  
  curr.invaded.all <- spread_ALBERS[which(spread_ALBERS$YrInv %in% min(na.omit(spread_ALBERS$YrInv)):i),]
  int = suppressWarnings(st_intersects(wedges, st_centroid(curr.invaded.all)))
  for(r in 1:nrow(int)){
    # r <- 2
    curr_wedge_polys <- int[r][[1]]
    curr.invaded.wedge <- curr.invaded.all[curr_wedge_polys,]
    if(nrow(curr.invaded.wedge)>0){
    max_val  <- max(curr.invaded.wedge$dist_from_disc)} else{
      max_val <- 0}
    
    if(r == 1){
      fin_dists <- max_val}else{
      fin_dists <- append(fin_dists,max_val)}
  }
  
  spread.rad[which(spread.rad$year == i), paste("radii",radii_n, sep=".") ] <-  fin_dists
  
}


spread.rad_raw <- spread.rad

i=2
for(i in 2:ncol(spread.rad)){
  curr.col <- spread.rad[,i]
  spread.rad[,i] <- c(NA,((diff(curr.col))))
}

spread.rad <- na.omit(spread.rad)
summary(spread.rad)

sprd.rads.lwd <- spread.rad %>% gather(key = "bearing", value="spread_increment", -year)
sprd.rads.lwd$bearing <- as.numeric(substr(sprd.rads.lwd$bearing,7,nchar(sprd.rads.lwd$bearing)))

# write a for loop to remove years after which the farthest point in an interval was reached
wedge_id <- data.frame(wedge_n = radii_n, year=NA)
int_all = st_intersects(wedges, st_centroid(spread_ALBERS))
# AQUI use boundaries instead?
#int_all = st_intersects(wedges, spread_ALBERS)

for(r in 1:nrow(int_all)){
  # r <- 13 # radii_n[13]
  curr_wedge_polys <- int_all[r][[1]]
  curr.wedge <- spread_ALBERS[curr_wedge_polys,]
  
  if(nrow(curr.wedge)>0){
    max_val_cty <- curr.wedge[which(curr.wedge$dist_from_disc == max(curr.wedge$dist_from_disc)),]
    FIPS_FAR <- max_val_cty$FIPS
  }else{FIPS_FAR <- "00001"}
  
  invaded_cty_wedge <- curr.wedge[which(curr.wedge$Invaded=="yes"),]
  if(nrow(invaded_cty_wedge)>0){
  max_val_invaded_cty <- invaded_cty_wedge[which(invaded_cty_wedge$dist_from_disc == max(invaded_cty_wedge$dist_from_disc)),]
  FIPS_INV_FAR <- max_val_invaded_cty$FIPS
  }else{FIPS_INV_FAR <- "00002"}
  if(FIPS_INV_FAR == FIPS_FAR){
    wedge_id[r,"year"] <- max_val_invaded_cty$YrInv 
  } else{
    wedge_id[r,"year"] <- 2021
  }
  }

  #plot(st_geometry(spread_ALBERS))
  #plot(st_geometry(curr.wedge), col="red", add=T)
bearings <- unique(sprd.rads.lwd$bearing) 
for(b in bearings){
  # b <- bearings[1]
  curr_bearing <- sprd.rads.lwd[which(sprd.rads.lwd$bearing == b),]
  curr_wedge_year_max <- wedge_id[which(wedge_id$wedge_n %in% b), "year"]
  curr_bearing_trunc <- curr_bearing[which(curr_bearing$year <= curr_wedge_year_max),]
if(b == bearings[1]){
  sprd.rads.lwd_proc <- curr_bearing_trunc} else{
    sprd.rads.lwd_proc <- rbind.data.frame(sprd.rads.lwd_proc,curr_bearing_trunc)}
  }
sprd.rads.lwd.g <- sprd.rads.lwd_proc %>% group_by(bearing) %>% summarise(mn_spread_km = mean(spread_increment), max_spread_km = sum(spread_increment))
summary(sprd.rads.lwd.g)
#https://www.r-graph-gallery.com/136-stacked-area-chart
ggplot(sprd.rads.lwd.g, aes(x = bearing, y = mn_spread_km)) +
  geom_area() +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 360, 60))+theme_bw() +
  xlab("Bearing")+
  ylab("Spread distance (km)")

col_bearings <- ggplot(sprd.rads.lwd, aes(x = bearing, y = spread_increment, fill=factor(year))) +
  geom_area() +
  theme_bw()+
  xlab("Bearing")+
  ylab("Spread distance")+ 
  scale_x_continuous(breaks = seq(0, 360, 30), limits = c(0,360), expand=c(0,0))+
  scale_y_continuous(breaks = seq(0, 1500, 300), limits = c(0,1500), expand=c(0,0))+
  theme_bw() +
  scale_fill_manual(values=my_LWD_colors, name="")+
  theme(legend.key.size = unit(0.3, 'cm'), legend.text = element_text(size=7), legend.justification = "left")+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


wedges_12 = st_wedges(discx,discy,wedge_length*1000,12)
st_crs(wedges_12) <- crs(spread_ALBERS)
wedges_4 = st_wedges(discx,discy,(wedge_length+400)*1000,12)
st_crs(wedges_4) <- crs(spread_ALBERS)
#
label_x <- (wedge_length+300)*1000*sin((pi*seq(0,330,30))/180)+discx
label_y <- (wedge_length+300)*1000*cos((pi*seq(0,330,30))/180)+discy

inset_circle <- ggplot() +
  geom_sf(data = wedges_4, fill="transparent", col="transparent") + # polygons filled based on the density value
  geom_sf(data = spread_ALBERS_cropped, aes(fill = YrInv), color="light gray") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers_crop, fill="transparent", color="black", lwd=1)+
  scale_fill_gradientn(colors = my_LWD_colors_spec, name="", na.value = "white",   guide = "none")+
  theme(legend.position="none",legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8), legend.justification = "left")+
  geom_sf(data=disc_loc_sf, color = "black", size = 3)+
  #geom_sf(data = wedges, fill="transparent", col="light gray") + # polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data =wedges_12, fill="transparent", col="black", size=0.5)+# polygons filled based on the density value
  annotate("text", label=seq(0,330,30), x=label_x, y=label_y, size=3)
#

label_x_small <- (wedge_length+300)*900*sin((pi*seq(0,330,30))/180)+discx
label_y_small <- (wedge_length+300)*900*cos((pi*seq(0,330,30))/180)+discy
unit_adj <- c(-0.5,-0.5,-0.5,-0.5)

inset_circle_smaller <- ggplot() +
  geom_sf(data = wedges_4, fill="transparent", col="transparent") + # polygons filled based on the density value
  geom_sf(data = spread_ALBERS_cropped, aes(fill = YrInv), color="light gray") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers_crop, fill="transparent", color="black", lwd=1)+
  scale_fill_gradientn(colors = my_LWD_colors_spec, name="", na.value = "white",   guide = "none")+
  theme(legend.position="none",legend.key.size = unit(0.5, 'cm'),
        legend.text = element_text(size=8), legend.justification = "left")+
  geom_sf(data=disc_loc_sf, color = "black", size = 3)+
  #geom_sf(data = wedges, fill="transparent", col="light gray") + # polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data =wedges_12, fill="transparent", col="black", size=0.5)+# polygons filled based on the density value
  annotate("text", label=seq(0,330,30), x=label_x, y=label_y, size=2.5)+
  theme(plot.margin=unit(unit_adj+2,"cm"))
  
#





# Figure: spread along 5° bearings wedges Morin --------------------------------------------
#
#
#---
resize.win(6.85,4)
gg_inset_bearings = ggdraw() +
  draw_plot(col_bearings) +
  draw_plot(inset_circle, x = 0.1, y = 0.55, width = 0.4, height = 0.4)
gg_inset_bearings

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''






# Figure: spread along 5° bearings wedges --------------------------------------------
#
#
#---
# https://stackoverflow.com/questions/43033628/circular-histogram-in-ggplot2-with-even-spacing-of-bars-and-no-extra-lines
sprd.rads.avg <- sprd.rads.lwd_proc %>% group_by(bearing) %>% summarise(mean_spread = mean(spread_increment))
#summary(sprd.rads.avg)
#sprd.rads.avg <- sprd.rads.avg[which(sprd.rads.avg$bearing %in% seq(0,360, by=5)),]
summary(sprd.rads.avg)
summary(sprd.rads.avg[which(sprd.rads.avg$bearing %!in% 45:157.5),])

sprd.rads.avg[order(sprd.rads.avg$mean_spread, decreasing=T),]
BD_col <- viridis(12)[11]
resize.win(6.85,6.85)

inset_circle_smaller

bearing5_histogram <- ggplot(sprd.rads.avg, aes(x = bearing, y = mean_spread)) +
  coord_polar(theta = "x", start = 0, clip="off") +
  geom_bar(stat = "identity", col= "black", fill = BD_col, width = 5) +
  geom_hline(yintercept = seq(0, 105, by = 15), color = "black", size = 0.2) +
  scale_y_continuous(breaks = seq(0,105,15), limits=c(0,105), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,315+22.5,22.5), expand= c(0,0)) +
  labs(x = "", y = "", title = "") +
  theme_bw()+ theme(panel.grid.minor = element_blank(),
                    #panel.grid.major.x = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major.x = element_line(colour="black", linetype="dashed"),
                    panel.grid.major.y = element_blank(), panel.border = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank())+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=12))+
  theme(plot.margin=unit(unit_adj,"cm"))+
  annotate("text", x=45,y= seq(15,105,15)-5, label= seq(15,105,15), size=2.5, hjust=0, vjust=1, fontface = 2)

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


# only analyze bearings with non-zero spread
#bearings_with_spread <- sprd.rads.lwd.g[which(sprd.rads.lwd.g$max_spread_km>0),]
non_ocean_bearings <- sprd.rads.lwd.g[which(sprd.rads.lwd.g$bearing %in% 45:155),]
spread_BD_year <- sprd.rads.lwd_proc[which(sprd.rads.lwd_proc$bearing %!in% non_ocean_bearings$bearing),]

BD_time_series <- spread_BD_year %>% group_by(year) %>% summarise(mn_spread_km = mean(spread_increment), max_spread_km = sum(spread_increment))
summary(BD_time_series)
round(mean(BD_time_series$mn_spread_km),0)
round(stderr(BD_time_series$mn_spread_km),0)
round(median(BD_time_series$mn_spread_km),0)

fig_BD <- ggplot(spread_BD_year, aes(x=year, y=spread_increment, fill=factor(year))) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=0.3)+
  scale_fill_manual(values=my_LWD_colors, name="")+
  ylab(expression ("Mean spread along intervals (km)"))+
  xlab("Year")+ theme_bw()+
  scale_x_continuous(breaks=seq(2005,2021,2))+
  expand_limits(x = 2004)+
  theme(legend.position="none",panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_text(x = 2005, y = 150, label = "Mean ± SE = 39 ± 7", parse = F, size=3,check_overlap = TRUE,  hjust = 0)+
  geom_text(x = 2005, y = 135, label = "Median = 35", parse = F, size=3,check_overlap = TRUE,  hjust = 0)


# NOTE THAT spread distances jump further here (with BD method) because you are
# just measuring jumps along radii, and not jumps from nearest county
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''




# Figure: spread different ERR segmented --------------------------------------------
#
#
#---
resize.win(6.85,6)
ggarrange(fig_ERR_COMP, fig_ERR_SEGMENT, ncol = 1, nrow = 2,  align = "hv",
          labels="auto", label.x =0.09, label.y =1)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''







# Figure: spread different methods --------------------------------------------
#
#
#---
resize.win(6.85,5)
ggarrange(fig_DR, fig_BD, ncol = 1, nrow =2,  align = "hv",
          labels="auto", label.x =0.12, label.y =1)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''








# Figure: no. counties contiguous vs. non-contiguous spread --------------------------------------------
#
#
#---
adj_col <- rev(viridis(6))[5]
iso_col <- rev(viridis(6))[2]
transp_gray <- alpha("light gray",0.3) 
LWD.adj.iso_fin.crp <- st_crop(contigUS_ALBERS_LWD, spread_ALBERS_cropped)
inv.isoadj_map <- ggplot() +
  geom_sf(data = LWD.adj.iso_fin.crp, aes(fill = JumpInv), color=transp_gray,lwd=0.0001) + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers_crop, fill="transparent", color="black", lwd=1)+
  scale_fill_manual(values= c(adj_col,iso_col),
                    name="", na.translate=FALSE)+
  theme(legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=8), legend.justification = "left",
        legend.position=c(0.3,0.25))



isoadj_cty <- LWD.adj.iso %>%  as.data.frame() %>% group_by(YrInv,STAT) %>% summarize(counties_invaded = n())
add_year0 <- data.frame(YrInv=2004, STAT="iso",  counties_invaded=0)
isoadj_cty <- rbind.data.frame(add_year0,isoadj_cty)

sum(table(LWD.adj.iso$STAT))
table(LWD.adj.iso$STAT)
year_county_summary <- as.data.frame(LWD.adj.iso) %>% group_by(YrInv) %>% summarize(counties_invaded = n())
summary(year_county_summary$counties_invaded)
stderr(year_county_summary$counties_invaded)

inv.isoadj_graph <- 
  ggplot(isoadj_cty,  aes(x=YrInv, y=counties_invaded, fill=STAT)) +
  geom_area( aes(fill=STAT), position="stack")+
  labs(title="",x="Year", y = "Counties invaded")+theme_clean+
  scale_fill_manual(values=c(adj_col, iso_col),labels = c("Contiguous", "Non-contiguous"))+
  scale_x_continuous(expand = c(0, 0))+ # x axis is getting cut-off
  scale_y_continuous(expand = c(0, 0))+ # x axis is getting cut-off
  theme(legend.position="none",legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=8),
        legend.justification = "left", legend.title = element_blank(),  plot.margin = margin(10, 10, 10, 10))+
  annotate('text', x = 2005, y = 30, label =  paste("italic(n)", "~invaded~counties==", nrow(LWD.adj.iso)), parse = T, size=3, check_overlap = TRUE,  hjust = 0)+
  geom_text(x = 2005, y = 26, label = "203 contiguous", parse = F, size=3,check_overlap = TRUE,  hjust = 0)+
  geom_text(x = 2005, y = 22, label = "72 non-contiguous", parse = F, size=3,check_overlap = TRUE,  hjust = 0)


# legend.position=c(0.1,0.8)
bquote(italic("n")~"="~.nrow(LWD.adj.iso$STAT)~"invaded counties")

iso_only <- LWD.adj.iso[which(LWD.adj.iso$STAT == "iso"),]
summary(iso_only$D_BDY_km)
round(mean(iso_only$D_BDY_km),0)
round(stderr(iso_only$D_BDY_km),0)
round(median(iso_only$D_BDY_km),0)
iso_only[order(iso_only$D_BDY_km, decreasing = T),]
quantile(iso_only$D_BDY_km, 0.99)

inv.isoadj_hist <-
  ggplot(iso_only, aes(x=D_BDY_km)) +
  geom_histogram(position="identity", col="black", fill=iso_col, boundary=0, bins=60)+
  #geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
  #           linetype="dashed")+
  labs(title="",x="Distance jumped (km)", y = "Count")+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0), limits= c(0,600))+ # x axis is getting cut-off
  scale_y_continuous(expand = c(0, 0), limits= c(0,15))+ # x axis is getting cut-off
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),  plot.margin = margin(10, 10, 10, 10))+
  geom_text(x = 200, y = 10, label = "Mean ± SE = 164 ± 16 km", parse = F, size=3,check_overlap = TRUE,  hjust = 0)+
  geom_text(x = 200, y = 8, label = "Median = 103 km", parse = F, size=3,check_overlap = TRUE,  hjust = 0)

gc()
memory.limit()
memory.limit(size=1e6)
resize.win(3.30709,7)
ggarrange(inv.isoadj_map, inv.isoadj_graph,inv.isoadj_hist, ncol = 1, nrow = 3, align= "v",
          labels="auto")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''









# Figure: spatial extents considered --------------------------------------------
#
#
#---
host_extent <- contigUS_ALBERS_LWD[which(contigUS_ALBERS_LWD$FIPS %in% c(counties_with_redbay$FIPS,counties_with_sassafras$FIPS)),]
inavded_area_hosts <- ggplot() +
  geom_sf(data =  host_extent, color="light gray",lwd=0.001) + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers, fill="transparent", color="black", lwd=1)
  
counties_within_buffer_extent <- contigUS_ALBERS_LWD[which(contigUS_ALBERS_LWD$FIPS %in% counties_within_buffer$FIPS),]
inavded_area_buffer <- ggplot() +
  geom_sf(data = counties_within_buffer_extent, color="light gray",lwd=0.001) + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers, fill="transparent", color="black", lwd=1)


inavded_area_US <- ggplot() +
  geom_sf(data = contigUS_ALBERS_LWD, color="light gray",lwd=0.001) + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers, fill="transparent", color="black", lwd=1)

resize.win(7,6)
plot(1~1)
inavded_area_hosts
inavded_area_buffer
inavded_area_US
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

























# Sixteen wedges (anisotropy): spread rates in each wedge --------------------------------------------
#
#
#---
spread.rad.16 <- data.frame(year=2004:2021)
radii_n_16 <- seq(1,16,1); length(radii_n_16)
# Create 16 wedges
discx <- coordinates(disc_loc_ALBERs)[1]
discy <- coordinates(disc_loc_ALBERs)[2]
summary(spread_ALBERS[which(spread_ALBERS$YrInv %in% 2004:2021), "dist_from_disc"])
summary(LWD.adj.iso$DtoGZ_km)
summary(isolated_counties$D_BDY_km)
wedges_16 = st_wedges(discx,discy,wedge_length*1000, nsegs=16) # wedge length 

# for adding wedge membership as a predictor
wedges_16_long = st_wedges(discx,discy,wedge_length*3000, nsegs=16) # wedge length 
st_crs(wedges_16_long) <- crs(spread_ALBERS)
st_write(wedges_16_long, "shapefiles/wedges/wedges_16_long.shp", delete_layer=T)
#

st_crs(wedges_16) <- crs(spread_ALBERS)
par(mfrow=c(1,1))
st_crs(wedges_16) <- crs(spread_ALBERS)
plot(wedges_16)
plot(wedges_16[1], add=T, col="red")
wedges_16 <- wedges_16[c(length(radii_n_16),1:length(radii_n_16)-1)]
plot(wedges_16[1], add=T, col="blue")

# set up data frame to be populated
spread.rad.16[,paste("radii",radii_n_16, sep=".")] <- NA

# write a for loop to remove years after which the farthest point in an interval was reached
wedge_id_16 <- data.frame(wedge_n_16 = radii_n_16, year=NA)
int_all_16 = st_intersects(wedges_16, st_centroid(spread_ALBERS))
for(r in 1:nrow(int_all_16)){
  # r <- 13 # radii_n_16[13]
  curr_wedge_polys <- int_all_16[r][[1]]
  curr.wedge <- spread_ALBERS[curr_wedge_polys,]
  
  if(nrow(curr.wedge)>0){
    max_val_cty <- curr.wedge[which(curr.wedge$dist_from_disc == max(curr.wedge$dist_from_disc)),]
    FIPS_FAR <- max_val_cty$FIPS
  }else{FIPS_FAR <- "00001"}
  
  invaded_cty_wedge <- curr.wedge[which(curr.wedge$Invaded=="yes"),]
  if(nrow(invaded_cty_wedge)>0){
    max_val_invaded_cty <- invaded_cty_wedge[which(invaded_cty_wedge$dist_from_disc == max(invaded_cty_wedge$dist_from_disc)),]
    FIPS_INV_FAR <- max_val_invaded_cty$FIPS
  }else{FIPS_INV_FAR <- "00002"}
  if(FIPS_INV_FAR == FIPS_FAR){
    wedge_id_16[r,"year"] <- max_val_invaded_cty$YrInv 
  } else{
    wedge_id_16[r,"year"]<- 2021
  }
}
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



# Invaded area and US map with wedges --------------------------------------------
#
#
#---
which(1:length(wedges) %!in% c(10:31))
wedges_padj <- wedges[which(1:length(wedges) %!in% c(10:31))]
bearings_16_remove <- -c(3:7)
wedges_16_padj <- wedges_16[bearings_16_remove]

invaded_area_discrete <- ggplot() +
  geom_sf(data = spread_ALBERS, aes(fill = (YrInv)), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers, fill="transparent", color="black", lwd=1)+
  scale_fill_viridis(
                    name="", na.value = "white", guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", direction = "horizontal"))+
  theme(legend.position=c(0.49,0.25),legend.key.height =  unit(0.2, 'cm'),
        legend.text = element_text(size=7), legend.justification = "left")+
  geom_sf(data=disc_loc_sf, color = "yellow", fill="black", size = 3, shape=21)+
  ggrepel::geom_text_repel(
    data = disc_loc_sf,
    aes(label = NAME, geometry = geometry),
    stat = "sf_coordinates",
    size=2.8,
    min.segment.length = 0,
    nudge_x = 0.5e6,
    nudge_y = 0.8e5,
    colour = "black",
    segment.colour = "black")+
  geom_sf(data = wedges_16_padj, fill="transparent", color="transparent") + #polygons filled based on the density value
  coord_sf(datum=st_crs(spread_ALBERS_cropped))
invaded_area_discrete


label_x_1 <- (wedge_length+100)*1000*sin((pi*seq(0,337.5,22.5)[-c(4:7)])/180)+discx
label_y_1 <- (wedge_length+100)*1000*cos((pi*seq(0,337.5,22.5)[-c(4:7)])/180)+discy

wedges_p <- ggplot() +
  geom_sf(data = spread_ALBERS, aes(fill = (YrInv)), color="transparent") + #polygons filled based on the density value
  theme_bw()+theme_void()+
  geom_sf(data=statesUS_albers, fill="transparent", color="black", lwd=1)+
  scale_fill_viridis(
    name="", na.value = "white", guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", direction = "horizontal"))+
  theme(legend.position="none")+
  geom_sf(data=statesUS_albers, fill="transparent", color="light gray", lwd=1)+
  geom_sf(data = wedges_padj, fill="transparent", color="gray") + #polygons filled based on the density value
  geom_sf(data = wedges_16_padj, fill="transparent", color="black", size=1) + #polygons filled based on the density value
  geom_sf(data=disc_loc_sf, color = "yellow", fill="black", size = 3, shape=21)+
  coord_sf(datum=st_crs(spread_ALBERS_cropped))+
  annotate("text", label=seq(0,337.5,22.5)[-c(4:7)], x=label_x_1, y=label_y_1, size=3, fontface=2)
#
wedges_p
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# Figure: plot with invaded area and US map with wedges --------------------------------------------
#
#
#---
resize.win(6.85,8)
ggarrange(invaded_area_discrete,wedges_p, ncol=1, nrow=2,
          labels="auto", label.x = 0.13, label.y =1)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



  
# Calculate spread rates using effective range radius and distance regression --------------------------------------------
#
#
#---
for(w in 1:length(wedges_16)){
  # w <- 2
  curr_wedge <- wedges_16[w]# 
  
  ctys.rows.curr.wedge = suppressWarnings(st_intersects(curr_wedge, st_centroid(spread_ALBERS)))
  curr.invaded.wedge <- spread_ALBERS[ctys.rows.curr.wedge[[1]],] # use row locations to extract counties
  curr.invaded.wedge <- curr.invaded.wedge[which(curr.invaded.wedge$Invaded=="yes"),] # use row locations to extract counties
  
  #plot(curr_wedge)
  #plot(curr.invaded.wedge, add=T, col="red")
  curr_wedge_n_16 <- w
  
  #'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  #
  # effective range radius
  #
  #'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  annual_spread_ERR <- data.frame(year=2004:2021, area_invaded_km2=NA)
  years.to.plot_ERR <- sort(unique(spread_ALBERS$YrInv))
  #i <- years.to.plot_ERR[1]
  for(i in years.to.plot){
    curr.invaded_ERR <- curr.invaded.wedge[which(curr.invaded.wedge$YrInv %in% 2004:i),]
    annual_spread_ERR[which(annual_spread_ERR$year %in% i), "year"] <- i
    annual_spread_ERR[which(annual_spread_ERR$year==i), "area_invaded_km2"] <- sum(st_area(curr.invaded_ERR)/1e6)
  }
  
  
  annual_spread_ERR$radius <- sqrt(16*annual_spread_ERR$area_invaded_km2/pi) # AQUI - just checking on 16 :)
  
  if(nrow(curr.invaded.wedge)>1 & length(unique(curr.invaded.wedge$YrInv))>1){
    
    min_yr <- 2004 
    max_yr <- wedge_id_16[wedge_id_16$wedge_n_16==curr_wedge_n_16,"year"]
    
    annual_spread_ERR <- annual_spread_ERR[which(annual_spread_ERR$year %in% min_yr:max_yr),]
    fit_curr_ERR <- lm(radius~year, data=annual_spread_ERR)
    summary(fit_curr_ERR)
    ERR_t_val <- coef(summary(fit_curr_ERR))[2,3]
    ERR_p_val <- coef(summary(fit_curr_ERR))[2,4]
    ERR_spread_rate_mu <- coef(summary(fit_curr_ERR))[2,1]
    ERR_spread_rate_sd <- coef(summary(fit_curr_ERR))[2,2]
    ERR_R_squared <- summary(fit_curr_ERR)$"adj.r.squared"}else{
      ERR_t_val <- 0
      ERR_p_val <- 0
      ERR_spread_rate_mu <- 0
      ERR_spread_rate_sd <- 0
      ERR_R_squared <- 0}
  
  
  
  #'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  #
  # distance regression
  #
  #'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  if(nrow(curr.invaded.wedge)>1 & length(unique(curr.invaded.wedge$YrInv))>1){
  
    min_yr <- 2004 #AQUI: this was min(curr.invaded.wedge$YrInv)...but better to start at initial detection year, right?
    max_yr <- wedge_id_16[wedge_id_16$wedge_n_16==curr_wedge_n_16,"year"]
      
  curr.invaded.wedge_DR <- curr.invaded.wedge[which(curr.invaded.wedge$YrInv %in% min_yr:max_yr),]
    
  fit_curr_DR <- lm(dist_from_disc~YrInv, data=curr.invaded.wedge_DR)
  summary(fit_curr_DR)
  DR_t_val <- coef(summary(fit_curr_DR))[2,3]
  DR_p_val <- coef(summary(fit_curr_DR))[2,4]
  DR_spread_rate_mu <- coef(summary(fit_curr_DR))[2,1]
  DR_spread_rate_sd <- coef(summary(fit_curr_DR))[2,2]
  DR_R_squared <- summary(fit_curr_DR)$"adj.r.squared"}else{
    DR_t_val <- 0
    DR_p_val <- 0
    DR_spread_rate_mu <- 0
    DR_spread_rate_sd <- 0
    DR_R_squared <- 0}

   curr_df_DR <- cbind.data.frame(curr_wedge_n_16,ERR_t_val,ERR_p_val,ERR_spread_rate_mu,ERR_spread_rate_sd,ERR_R_squared,
                                  DR_t_val,DR_p_val,DR_spread_rate_mu,DR_spread_rate_sd,DR_R_squared)
  
  #curr_df_DR <- cbind.data.frame(curr_wedge_n_16,DR_t_val,DR_p_val,DR_spread_rate_mu,DR_spread_rate_sd,DR_R_squared)
  if(w ==1){
    spread_distance_reg_df  <- curr_df_DR}else{
    spread_distance_reg_df  <- rbind.data.frame(spread_distance_reg_df,curr_df_DR)
    }
}

summary(spread_distance_reg_df)



#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#
# calculate spread rates using boundary displacement in each wedge
#
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
i <- 2012
my_wedge_colors <- sample(rainbow(16),16)

par(mfrow=c(4,4))
for(i in years.to.plot){
  curr.invaded.all <- spread_ALBERS[which(spread_ALBERS$YrInv %in% min(na.omit(spread_ALBERS$YrInv)):i),]
  # counties contained in each wedge
  int = suppressWarnings(st_intersects(wedges_16, st_centroid(curr.invaded.all)))

    for(r in 1:nrow(int)){ # nrow(int) should = 16, the number of wedges
    # r <- 1
    curr_wedge_polys <- int[r][[1]] # invaded county (row location) within wedge r
    curr.invaded.wedge <- curr.invaded.all[curr_wedge_polys,] # use row locations to extract counties
    if(nrow(curr.invaded.wedge) > 0){
    max_val  <- max(curr.invaded.wedge$dist_from_disc)} else{ max_val<- 0}
    
    
    #plot(wedges_16)
    # plot(wedges_16[1], col="red", add=T)
    # plot(curr.invaded.wedge, add=T, col="gray")
    if(r == 1){
      par(mfrow=c(4,4))
      plot(wedges_16)
      plot(wedges_16[r], col="gray", add=T)
      plot(st_geometry(curr.invaded.wedge), add=T, col="red")
      fin_dists <- max_val}else{
        fin_dists <- append(fin_dists,max_val)
        plot(wedges_16)
        plot(wedges_16[r], col="gray", add=T)
        plot(st_geometry(curr.invaded.wedge), add=T, col="red")}
  }
    spread.rad.16[which(spread.rad.16$year == i), paste("radii",radii_n_16, sep=".") ] <-  fin_dists
}

spread.rad.16.raw <- spread.rad.16

i=2
for(i in 2:ncol(spread.rad.16)){
  curr.col <- spread.rad.16[,i]
  spread.rad.16[,i] <- c(NA,((diff(curr.col))))
}

spread.rad.16 <- na.omit(spread.rad.16)
summary(spread.rad.16)
sprd.rads16.lwd <- spread.rad.16 %>% gather(key = "bearing", value="spread_increment", -year)
sprd.rads16.lwd$bearing <- as.numeric(substr(sprd.rads16.lwd$bearing,7,nchar(sprd.rads16.lwd$bearing)))



#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# write a for loop to remove years after which the farthest point in an interval was reached
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
table(sprd.rads16.lwd$year)
bearings_16 <- unique(sprd.rads16.lwd$bearing) 
for(b in bearings_16){
  # b <- bearings_16[1]
  curr_bearing <- sprd.rads16.lwd[which(sprd.rads16.lwd$bearing == b),]
  curr_wedge_year_max <- wedge_id_16[which(wedge_id_16$wedge_n_16 %in% b), "year"]
  curr_bearing_trunc <- curr_bearing[which(curr_bearing$year <= curr_wedge_year_max),]
  if(b == bearings_16[1]){
    sprd.rads.lwd_proc_16 <- curr_bearing_trunc} else{
      sprd.rads.lwd_proc_16 <- rbind.data.frame(sprd.rads.lwd_proc_16,curr_bearing_trunc)}
}
table(sprd.rads.lwd_proc_16$year)

sprd.rads16.lwd.g <- sprd.rads.lwd_proc_16 %>% group_by(wedge_n_16=bearing) %>% summarise(BD_mn_spread_km = mean(spread_increment), BD_max_spread_km = sum(spread_increment), BD_se_spread_lm = stderr(spread_increment))
summary(sprd.rads16.lwd.g)


colnames(spread_distance_reg_df)[1]
colnames(spread_distance_reg_df)[1] <- "wedge_n_16"
spread_rates_wedge <- merge(spread_distance_reg_df,sprd.rads16.lwd.g,by="wedge_n_16")
spread_rates_wedge$bearing <- seq(11.25,360-11.25,22.5)


# remove the bearings without spread
spread_rates_wedge <- spread_rates_wedge[bearings_16_remove, ]
spread_rates_wedge


# do spread rates match up between methods - each bearing interval represents a point
plot(DR_spread_rate_mu~BD_mn_spread_km , data=spread_rates_wedge)
fit_spread_comp <- lm(DR_spread_rate_mu~BD_mn_spread_km , data=spread_rates_wedge)
summary(fit_spread_comp)
confint(fit_spread_comp)





ERR_col <- viridis(6)[3]
DR_col <- viridis(6)[5]


resize.win(4,4)
# https://stackoverflow.com/questions/43033628/circular-histogram-in-ggplot2-with-even-spacing-of-bars-and-no-extra-lines
fig_ERR_hist <- 
  ggplot(spread_rates_wedge, aes(x = bearing, y = ERR_spread_rate_mu  )) +
  coord_polar(theta = "x", start = 0, clip="off") +
  geom_bar(stat = "identity", col= "black", fill = ERR_col, width = 22.5) +
  geom_hline(yintercept = seq(0, 90, by = 15), color = "black", size = 0.2) +
  scale_y_continuous(breaks = seq(0,90,15), limits=c(0,90), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,315,45), expand= c(0,0)) +
  labs(x = "", y = "", title = "") +
  theme_bw()+ theme(panel.grid.minor = element_blank(),
                    #panel.grid.major.x = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major.x = element_line(colour="black", linetype="dashed"),
                    panel.grid.major.y = element_blank(), panel.border = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank())+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=12))+
  theme(plot.margin=unit(unit_adj,"cm"))+
  annotate("text", x=45,y= seq(15,90,15)-5, label= seq(15,90,15), size=2.5, hjust=0, vjust=1, fontface =2)
  
fig_DR_hist <- ggplot(spread_rates_wedge, aes(x = bearing, y = DR_spread_rate_mu)) +
  coord_polar(theta = "x", start = 0, clip="off") +
  geom_bar(stat = "identity", col= "black", fill = DR_col, width = 22.5) +
  geom_hline(yintercept = seq(0, 90, by = 15), color = "black", size = 0.2) +
  scale_y_continuous(breaks = seq(0,90,15), limits=c(0,90), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,315,45), expand= c(0,0)) +
  labs(x = "", y = "", title = "") +
  theme_bw()+ theme(panel.grid.minor = element_blank(),
                    #panel.grid.major.x = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major.x = element_line(colour="black", linetype="dashed"),
                    panel.grid.major.y = element_blank(), panel.border = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank())+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=12))+
  theme(plot.margin=unit(unit_adj,"cm"))+
  annotate("text", x=45,y= seq(15,90,15)-5, label= seq(15,90,15), size=2.5, hjust=0, vjust=1, fontface =2)

summary(spread_rates_wedge$BD_mn_spread_km)
fig_BD_hist <- ggplot(spread_rates_wedge, aes(x = bearing, y = BD_mn_spread_km)) +
  coord_polar(theta = "x", start = 0, clip="off") +
  geom_bar(stat = "identity", col= "black", fill = BD_col, width = 22.5) +
  geom_hline(yintercept = seq(0, 90, by = 15), color = "black", size = 0.2) +
  scale_y_continuous(breaks = seq(0,90,15), limits=c(0,90), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0,315,45), expand= c(0,0)) +
  labs(x = "", y = "", title = "") +
  theme_bw()+ theme(panel.grid.minor = element_blank(),
                    #panel.grid.major.x = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major.x = element_line(colour="black", linetype="dashed"),
                    panel.grid.major.y = element_blank(), panel.border = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.text.y = element_blank())+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=12))+
  theme(plot.margin=unit(unit_adj,"cm"))+
  annotate("text", x=45,y= seq(15,90,15)-5, label= seq(15,90,15), size=2.5, hjust=0, vjust=1, fontface =2)

summary(spread_rates_wedge)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



# Figure: spread using different methods (histogram) --------------------------------------------
#
#
#---
resize.win(6.85,6)
ggarrange(bearing5_histogram,  fig_BD_hist,
          fig_ERR_hist,  fig_DR_hist, 
          ncol = 2, nrow = 2, labels = c("a","b","c","d"," "),
          label.x = 0.13, label.y =0.90)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''









# Loop to estimate predictors per wedge --------------------------------------------
# write a for loop that gets the mean of predictors per wedge and the spread rates per wedge using effective range radius
# distance regression, and boundary displacement
#---
CPH_df_hosts <- fread("data/CPH_data/CPH_LWD_COMBINED.csv", na=".")
CPH_df_hosts$FIPS <- ifelse(nchar(CPH_df_hosts$FIPS)==4, paste(0,CPH_df_hosts$FIPS,sep=""), paste(CPH_df_hosts$FIPS))
table(nchar(CPH_df_hosts$FIPS))
summary(CPH_df_hosts)
first_cty <- CPH_df_hosts[match(unique(CPH_df_hosts$FIPS), CPH_df_hosts$FIPS),]
nrow(first_cty)


for(w in 1:length(wedges_16)){
    # w <- 2
    curr_wedge <- wedges_16[w]# 
    ctys.rows.curr.wedge = st_intersects(curr_wedge, st_centroid(spread_ALBERS))
    #plot(wedges_16)
    #plot(st_geometry(spread_ALBERS), add=T)
    curr.invaded.wedge <- spread_ALBERS[ctys.rows.curr.wedge[[1]],] # use row locations to extract counties
    preds_in_wedge <- first_cty[which(first_cty$FIPS %in% curr.invaded.wedge$FIPS),]
    
    preds_df_means <- preds_in_wedge %>% summarize(ln.tot.pop.mu = mean(log(POPESTIMATE2019)),
                                                   med_income.mu = mean(med_income),
                                                   ln.n_campgrounds.mu = mean(log(n_campgrounds+1)),
                                                   ln.ltot.mu = mean(log(ltot+1)),
                                                   MINMOTEMP.mu = mean(MINMOTEMP),
                                                   MAP.mu.m = mean(MAP)/1000,
                                                   ln.host_biomass.mu = mean(log(host_biomass+1)),
                                                   ln.non_host_biomass.mu = mean(log(non_host_biomass+1)),
                                                   ln.redbay_biomass.mu = mean(log(redbay_biomass+1)),
                                                   ln.sassafras_biomass.mu = mean(log(sassafras_biomass+1))) 
    preds_df_means$wedge_n_16 <- w
    
    
    ctys.rows.curr.wedge_jmp = suppressWarnings(st_intersects(curr_wedge, st_centroid(isolated_counties)))
    curr.invaded.wedge_jmp <- isolated_counties[ctys.rows.curr.wedge_jmp[[1]],] # use row locations to extract counties
    
    
    mean_distance_jump <- mean(curr.invaded.wedge_jmp$D_BDY_km)
    max_distance_jump <- max(curr.invaded.wedge_jmp$D_BDY_km)
    
    jumps_in_wedge <- length(ctys.rows.curr.wedge_jmp[[1]])
    
    
    preds_df_means <- cbind.data.frame(preds_df_means,mean_distance_jump,max_distance_jump,jumps_in_wedge)
    if(w ==1){
      preds_df_means_master  <- preds_df_means}else{
        preds_df_means_master  <- rbind.data.frame(preds_df_means_master,preds_df_means)
      }
  }

preds_df_means_master
warnings()

final_aniso <- merge(preds_df_means_master,spread_rates_wedge, by="wedge_n_16")
bearings_16_remove


# Bonferroni correction
p <- 0.05/21
qnorm(p, lower.tail = F)

summary(final_aniso)

plot(ERR_spread_rate_mu~ln.tot.pop.mu, data=final_aniso)
fit_ERR_pop <- lm(ERR_spread_rate_mu~ln.tot.pop.mu, data=final_aniso)
summary(fit_ERR_pop)

#
#
plot(ERR_spread_rate_mu~ln.host_biomass.mu, data=final_aniso)
fit_ERR_host <- lm(ERR_spread_rate_mu~ln.host_biomass.mu, data=final_aniso)
summary(fit_ERR_host)

#
plot(ERR_spread_rate_mu~ln.sassafras_biomass.mu, data=final_aniso)
fit_ERR_sas <- lm(ERR_spread_rate_mu~ln.sassafras_biomass.mu , data=final_aniso)
summary(fit_ERR_sas)

#
plot(ERR_spread_rate_mu~ln.redbay_biomass.mu, data=final_aniso)
fit_ERR_red <- lm(ERR_spread_rate_mu~ln.redbay_biomass.mu, data=final_aniso)
summary(fit_ERR_red)


#
plot(ERR_spread_rate_mu~MINMOTEMP.mu, data=final_aniso)
fit_ERR_temp <- lm(ERR_spread_rate_mu~MINMOTEMP.mu, data=final_aniso)
summary(fit_ERR_temp)

#
plot(ERR_spread_rate_mu~MAP.mu.m, data=final_aniso)
fit_ERR_precip <- lm(ERR_spread_rate_mu~MAP.mu.m, data=final_aniso)
summary(fit_ERR_precip)


#
fit_ERR_tempred <- lm(ERR_spread_rate_mu~MINMOTEMP.mu+ln.redbay_biomass.mu, data=final_aniso)
summary(fit_ERR_tempred)
cor(final_aniso$MINMOTEMP.mu,final_aniso$ln.redbay_biomass.mu)
cor(final_aniso$MINMOTEMP.mu,final_aniso$ln.sassafras_biomass.mu)

plot(ERR_spread_rate_mu~mean_distance_jump, data=final_aniso)
fit_ERR_jump <- lm(ERR_spread_rate_mu~mean_distance_jump, data=final_aniso)
summary(fit_ERR_jump)



plot(DR_spread_rate_mu~ln.tot.pop.mu, data=final_aniso)
fit_DR_pop <- lm(DR_spread_rate_mu~ln.tot.pop.mu, data=final_aniso)
summary(fit_DR_pop)
#
plot(DR_spread_rate_mu~ln.redbay_biomass.mu, data=final_aniso)
fit_DR_red <- lm(DR_spread_rate_mu~ln.redbay_biomass.mu, data=final_aniso)
summary(fit_DR_red)
#
plot(DR_spread_rate_mu~MINMOTEMP.mu, data=final_aniso)
fit_DR_temp <- lm(DR_spread_rate_mu~MINMOTEMP.mu, data=final_aniso)
summary(fit_DR_temp)
#
plot(DR_spread_rate_mu~MAP.mu.m, data=final_aniso)
fit_DR_precip <- lm(DR_spread_rate_mu~MAP.mu.m, data=final_aniso)
summary(fit_DR_precip)
#
plot(DR_spread_rate_mu~mean_distance_jump, data=final_aniso)
fit_DR_jump <- lm(DR_spread_rate_mu~mean_distance_jump, data=final_aniso)
summary(fit_DR_jump)





plot(BD_mn_spread_km ~ln.tot.pop.mu, data=final_aniso)
fit_BD_pop <- lm(BD_mn_spread_km ~ln.tot.pop.mu, data=final_aniso)
summary(fit_BD_pop)
#
plot(BD_mn_spread_km ~ln.redbay_biomass.mu, data=final_aniso)
fit_BD_red <- lm(BD_mn_spread_km~ln.redbay_biomass.mu, data=final_aniso)
summary(fit_BD_red)
#
plot(BD_mn_spread_km ~MINMOTEMP.mu, data=final_aniso)
fit_BD_temp <- lm(BD_mn_spread_km ~ MINMOTEMP.mu, data=final_aniso)
summary(fit_BD_temp)
#
plot(BD_mn_spread_km ~MAP.mu.m, data=final_aniso)
fit_BD_temp <- lm(BD_mn_spread_km ~ MAP.mu.m, data=final_aniso)
summary(fit_BD_temp)
#
#
plot(BD_mn_spread_km ~mean_distance_jump, data=final_aniso)
fit_BD_jump <- lm(BD_mn_spread_km ~mean_distance_jump, data=final_aniso)
summary(fit_BD_jump)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''











# Figure: drivers of anisotropy --------------------------------------------
#
#
#---
lm_eqn_red <- function(df){
  m <- lm(ERR_spread_rate_mu ~ ln.redbay_biomass.mu, df);
  eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}

lm_eqn_temp <- function(df){
  m <- lm(ERR_spread_rate_mu ~ MINMOTEMP.mu, df);
  eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}

lm_eqn_precip <- function(df){
  m <- lm(ERR_spread_rate_mu ~ MAP.mu.m, df);
  eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$adj.r.squared, digits = 2)))
  as.character(as.expression(eq));
}

reg_size_txt <- 2.6


summary(fit_ERR_red)
summary(fit_ERR_temp)
summary(fit_ERR_precip)

# Figure host
fig_ERR_red_aniso <- ggplot(data=final_aniso, aes(x=ln.redbay_biomass.mu, y=ERR_spread_rate_mu)) +
  ylab(expression ("Spread rate (km·"*~year^-1*")"))+
  xlab("Redbay biomass (ln(kg))")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(0,20,5), limits=c(0,20), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,60,10), limits=c(0,60), expand = c(0,0)) +
  
  geom_point(col=ERR_col, size=2)+
  stat_smooth(data = fit_ERR_red, method = "lm", col = "black", se=F)+
  geom_text(x = 4, y = 10, label = lm_eqn_red(final_aniso), parse = TRUE, size=reg_size_txt, check_overlap = TRUE, hjust = 0)+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Figure temp
fig_ERR_temp_aniso <- ggplot(data=final_aniso, aes(x=MINMOTEMP.mu, y=ERR_spread_rate_mu)) +
  ylab(expression (" "))+
  xlab("Mininum temperature (°C)")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(-8,11,4), limits=c(-8,11), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,60,10), limits=c(0,60), expand = c(0,0)) +
  
  geom_point(col=ERR_col, size=2)+
  stat_smooth(data = fit_ERR_temp, method = "lm", col = "black", se=F)+
  geom_text(x = -4, y = 10, label = lm_eqn_temp(final_aniso), parse = TRUE, size=reg_size_txt, check_overlap = TRUE, hjust = 0)+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Figure precip
fig_ERR_precip_aniso <- ggplot(data=final_aniso, aes(x=MAP.mu.m, y=ERR_spread_rate_mu)) +
  ylab(expression (" "))+
  xlab("Total precipitation (m)")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(1,1.5,0.1), limits=c(1,1.5), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,60,10), limits=c(0,60), expand = c(0,0)) +
  
  geom_point(col=ERR_col, size=2)+
  stat_smooth(data = fit_ERR_precip, method = "lm", col = "black", se=F)+
  geom_text(x = 1.12, y = 10, label = lm_eqn_precip(final_aniso), parse = TRUE, size=reg_size_txt, check_overlap = TRUE, hjust = 0)+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



resize.win(6.85,2.3)
ggarrange(fig_ERR_red_aniso, fig_ERR_temp_aniso, fig_ERR_precip_aniso, ncol = 3, nrow = 1,
          labels="auto", label.x = 0.25, align="v")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# SFIWC FIGURE SPREAD --------------------------------------------
#
#
#---
y_initial_ERR <- 530
color_dots <- viridis(6)[1]
ggplot(data=annual_spread, aes(x=year, y=radius_try2)) +
  ylab("Radius of invaded area (km)")+
  xlab("Year")+ theme_bw()+
  
  scale_x_continuous(breaks = seq(2005,2022,4), limits=c(2003,2022), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,600,100), limits=c(0,600), expand = c(0,0)) +
  
  geom_point(aes(x = (year), y = radius_try2), col=seg_pt_col, size=6)+
  stat_smooth(data = annual_spread_BREAK_BEFORE, method = "lm", col = "black", se=F, size=1.2)+
  geom_text(x = 2005, y = y_initial_ERR, label = lm_eqn(annual_spread_BREAK_BEFORE), parse = TRUE, size=4, check_overlap = TRUE, hjust = 0)+
  geom_segment(aes(x = 2004, y = y_initial_ERR, xend =  2004+0.8, yend = y_initial_ERR), colour = "black",size=1.2)+
  
  stat_smooth(data = annual_spread_BREAK_AFTER, method = "lm", col = "black", se=F, linetype="dashed", size=1.2)+
  geom_text(x = 2005, y = y_initial_ERR-55, label = lm_eqn(annual_spread_BREAK_AFTER), parse = TRUE, size=4,check_overlap = TRUE,  hjust = 0)+
  geom_segment(aes(x = 2004, y = y_initial_ERR-55, xend =  2004+0.8, yend = y_initial_ERR-55), colour = "black",size=1.2, linetype="dashed")+
  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"))+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





# SFIWC FIGURE ANISOTROPY --------------------------------------------
#
#
#---

col_bearings <- ggplot(sprd.rads.lwd, aes(x = bearing, y = spread_increment, fill=factor(year))) +
  geom_area() +
  theme_bw()+
  xlab("Bearing")+
  ylab("Spread distance (km)")+ 
  scale_x_continuous(breaks = seq(0, 360, 30), limits = c(0,360), expand=c(0,0))+
  scale_y_continuous(breaks = seq(0, 1500, 300), limits = c(0,1500), expand=c(0,0))+
  theme_bw() +
  scale_fill_manual(values=my_LWD_colors, name="")+
  theme(legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=9), legend.justification = "left")+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"))+
  theme(axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))



gg_inset_bearings = ggdraw() +
  draw_plot(col_bearings) +
  draw_plot(inset_circle, x = 0.1, y = 0.55, width = 0.4, height = 0.4)
gg_inset_bearings

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



# County size through time --------------------------------------------
#
#
#---
spread_ALBERS$county_area_km <- spread_ALBERS$county_area/1e6
county_area_time <- as.data.frame(spread_ALBERS)
county_area_time$county_area_km <- as.numeric(county_area_time$county_area_km)
county_area_time$YrInv <- as.numeric(county_area_time$YrInv)
county_area_time$dist_from_disc <- as.numeric(county_area_time$dist_from_disc)

ggplot(county_area_time, aes(x=YrInv, y=county_area_km)) + 
  geom_point(color="gray")+
  geom_smooth(method=lm, color="black")+
  xlab('Year of invasion')+
  ylab(bquote('County area '(km^2)))+
  theme_classic() 

ggplot(county_area_time, aes(x=dist_from_disc, y=county_area_km)) + 
  geom_point(color="gray")+
  geom_smooth(method=lm, color="black")+
  xlab('Distance from discovery location')+
  ylab(bquote('County area '(km^2)))+
  theme_classic() 


summary(county_area_time)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


beep(1)


