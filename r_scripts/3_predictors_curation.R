# Load packages  --------------------------------------------
#
#
#---
source("r_scripts/packages.R")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


# Load in predictors --------------------------------------------
#
#
#---

#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#
# load shapefile for processing other predictors
#
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
US_ALBERS = st_read("shapefiles/invasion_data/uscounties_ALBERS.shp")
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''





#''''''''''''''''''''''''''''''''''''''''''''
#
# humans
#
#''''''''''''''''''''''''''''''''''''''''''''
humans <- read_excel("data/human_population_density/co-est2019-alldata.xlsx", range=c("A1:S3194"))
humans <- humans[,c("STATE", "COUNTY", "STNAME",  "CTYNAME", "POPESTIMATE2019")]
humans$STATE <- ifelse(nchar(humans$STATE)<2, paste("0",humans$STATE, sep=""), humans$STATE)
humans$COUNTY <- ifelse(nchar(humans$COUNTY)<3, paste("0",humans$COUNTY, sep=""), humans$COUNTY)
humans$COUNTY <- ifelse(nchar(humans$COUNTY)<3, paste("0",humans$COUNTY, sep=""), humans$COUNTY)
humans$FIPS <- paste(humans$STATE, humans$COUNTY, sep="")
humans <- humans[-which(humans$COUNTY == "000"),]
head(humans)
#
table(humans[which(humans$FIPS %!in% US_ALBERS$FIPS),"STNAME"])
nrow(US_ALBERS)
nrow(humans)
US_ALBERS <- merge(x = US_ALBERS, y = humans, by = "FIPS", all.x = TRUE)
nrow(US_ALBERS)
nrow(humans)
#''''''''''''''''''''''''''''''''''''''''''''







#''''''''''''''''''''''''''''''''''''''''''''
#
# income
#
#''''''''''''''''''''''''''''''''''''''''''''
county_incomeAB <- read_excel("data/county_income/Unemployment.xlsx", range=c("A5:B3201"))
county_incomeCM <- read_excel("data/county_income/Unemployment.xlsx", range=c("CM5:CM3201")) # median income 2019
county_income <- cbind.data.frame(county_incomeAB,county_incomeCM)
head(county_income)
#
county_income <- county_income[,-2]
colnames(county_income) <- c("FIPS","med_income")
#
nrow(US_ALBERS)
nrow(county_income)
US_ALBERS <- merge(x = US_ALBERS, y = county_income, by = "FIPS", all.x = TRUE)
nrow(US_ALBERS)
nrow(county_income)

#''''''''''''''''''''''''''''''''''''''''''''




#''''''''''''''''''''''''''''''''''''''''''''
#
# roads
#
#''''''''''''''''''''''''''''''''''''''''''''
roads = st_read("data/road_data/tl_2016_us_primaryroads.shp")
roads_ALBERS <- st_transform(roads,crs(US_ALBERS))
plot(st_geometry(roads))
ints = st_intersection(US_ALBERS,roads_ALBERS)
ltot <- tapply(st_length(ints), ints$FIPS,sum)
US_ALBERS$ltot = rep(0,nrow(US_ALBERS))
US_ALBERS$ltot[match(names(ltot),US_ALBERS$FIPS)] = ltot
plot(US_ALBERS[,"ltot"])
#''''''''''''''''''''''''''''''''''''''''''''






#''''''''''''''''''''''''''''''''''''''''''''
#
# campgrounds
#
#''''''''''''''''''''''''''''''''''''''''''''
us_campgrounds <- read_excel("data/campgrounds/campgrounds_US.xlsx")
head(us_campgrounds)
campgrounds_sf = st_as_sf(us_campgrounds, coords = c("longitude", "latitude"), 
                 crs = 4326)
plot(campgrounds_sf)
campgrounds_sf_ALBERS <- st_transform(campgrounds_sf , crs=crs(US_ALBERS))
US_ALBERS$n_campgrounds <- lengths(st_intersects(US_ALBERS, campgrounds_sf_ALBERS))
#''''''''''''''''''''''''''''''''''''''''''''






#''''''''''''''''''''''''''''''''''''''''''''
#
# host tree data
#
#''''''''''''''''''''''''''''''''''''''''''''
# # https://www.fs.usda.gov/rds/archive/catalog/RDS-2013-0013
# persea_spp <- raster("data/host_tree_data/s720.img")
# plot(persea_spp)
# persea <- raster("data/host_tree_data/s721.img")
# plot(persea)
# sassafras <- raster("data/host_tree_data/s931.img")
# plot(sassafras)
# 
# host_basal_area <- sum(persea,sassafras)
# host_basal_area_r <- rast(host_basal_area)
# host_basal_area_r <- terra::project(host_basal_area_r, crs(US_ALBERS))
# 
# host_basal_area_m <- aggregate(host_basal_area_r, fact=40, fun=mean)
# plot(host_basal_area_m)
# host_basal_area_s <- aggregate(host_basal_area_r, fact=40, fun=sum)
# plot(host_basal_area_s)
# #
# US_ALBERS$host_BA_mean <- extract(raster(host_basal_area_m), US_ALBERS, fun=mean)
# US_ALBERS$host_BA_sum <-  extract(raster(host_basal_area_m), US_ALBERS, fun=sum)
#''''''''''''''''''''''''''''''''''''''''''''


#''''''''''''''''''''''''''''''''''''''''''''
# Data from Songlin and Sandy
#''''''''''''''''''''''''''''''''''''''''''''
tree_data <- read.table("data/tree_biomass_county/SPP_Biomass_by_Cnty.txt", sep="\t", header=T) 
head(tree_data)
nrow(tree_data)

# fix state-level codes
table(nchar(as.character(tree_data$STATECD)))
tree_data$STATECD <- ifelse(nchar(tree_data$STATECD) == 1, paste("0", tree_data$STATECD, sep=""), paste(tree_data$STATECD))
table(nchar(as.character(tree_data$STATECD)))

# fix county-level codes
table(nchar(as.character(tree_data$COUNTYCD)))
tree_data$COUNTYCD <- ifelse(nchar(tree_data$COUNTYCD) == 1, paste("00", tree_data$COUNTYCD, sep=""), paste(tree_data$COUNTYCD))
table(nchar(as.character(tree_data$COUNTYCD)))
tree_data$COUNTYCD <- ifelse(nchar(tree_data$COUNTYCD) == 2, paste("0", tree_data$COUNTYCD, sep=""), paste(tree_data$COUNTYCD))
table(nchar(as.character(tree_data$COUNTYCD)))

# create fips
tree_data$FIPS <- as.factor(paste(tree_data$STATECD,tree_data$COUNTYCD, sep=""))
table(nchar(as.character(tree_data$FIPS)))
tree_data$Tree_Latin <- paste(tree_data$GENUS, tree_data$SPECIES, sep=" ")
head(tree_data)
nrow(tree_data)
summary(tree_data)

host_genus_county <- tree_data %>% filter(Tree_Latin %in% c("Persea borbonia", "Sassafras albidum"))
host_avocado_county <- tree_data %>% filter(Tree_Latin %in% c("Persea americana"))
table(host_genus_county$Tree_Latin)
#
US_ALBERS$host_biomass <- NA
US_ALBERS$non_host_biomass <- NA
i <- 1
for(i in 1:nrow(US_ALBERS)){
  curr_county <- tree_data[which(tree_data$FIPS %in% US_ALBERS$FIPS[i]),]
  #
  curr_county_hosts <- curr_county %>% filter(Tree_Latin %in%  c("Persea borbonia", "Sassafras albidum")) %>% 
    group_by(FIPS) %>% summarise(host_biomass = sum(BIOTOT), host_volume = sum(VOL_ALL)) 
  host_biomass_val <- ifelse(nrow(curr_county_hosts)==0, 0,curr_county_hosts$host_biomass)
  #
  curr_county_non_hosts <- curr_county %>% filter(Tree_Latin %!in% c("Persea borbonia", "Sassafras albidum")) %>% 
    group_by(FIPS) %>% 
    summarise(non_host_biomass = sum(BIOTOT), non_host_volume = sum(VOL_ALL)) 
  non_host_biomass_val <- ifelse(nrow(curr_county_non_hosts)==0, 0, curr_county_non_hosts$non_host_biomass)
  
  curr_county_redbay <- curr_county %>% filter(Tree_Latin %in%  c("Persea borbonia")) %>% 
    group_by(FIPS) %>% summarise(host_biomass = sum(BIOTOT), host_volume = sum(VOL_ALL)) 
  redbay_biomass_val <- ifelse(nrow(curr_county_redbay)==0, 0, curr_county_redbay$host_biomass)
  US_ALBERS$redbay_biomass[i] <- redbay_biomass_val
  
  
  curr_county_sassafras <-curr_county %>% filter(Tree_Latin %in%  c("Sassafras albidum")) %>% 
    group_by(FIPS) %>% summarise(host_biomass = sum(BIOTOT), host_volume = sum(VOL_ALL)) 
  sassafras_biomass_val <- ifelse(nrow(curr_county_sassafras)==0, 0, curr_county_sassafras$host_biomass )
  US_ALBERS$sassafras_biomass[i] <-   sassafras_biomass_val
 
  US_ALBERS$host_biomass[i] <-   host_biomass_val
  US_ALBERS$non_host_biomass[i] <-  non_host_biomass_val
}



US_ALBERS_hostplot <- US_ALBERS
ggplot(data = US_ALBERS_hostplot) +
  geom_sf() +
  geom_sf(data = US_ALBERS, size = 4, shape = 23, fill = "darkred") +
  coord_sf(xlim = c(-88, -78), ylim = c(24.5, 33), expand = FALSE)


US_ALBERS_hostplot$hosts_1 <- log(US_ALBERS_hostplot$host_biomass+1)
# ggplot(data = US_ALBERS_hostplot) +
#   geom_sf() +
#   geom_sf(data = US_ALBERS_hostplot, aes(fill = hosts_1)) +
#   scale_fill_manual(trans = "log", alpha = .4) +
#   coord_sf(xlim = c(-88, -78), ylim = c(24.5, 33), expand = FALSE)

ggplot(US_ALBERS_hostplot) +
  geom_sf(aes(fill = hosts_1))

ggplot() + 
  geom_sf(data = US_ALBERS_hostplot, aes(fill = hosts_1)) + 
  scale_y_continuous(breaks = 0:20)
# US_ALBERS_hostplot <- st_centroid(US_ALBERS_hostplot)
# pts <- do.call(rbind, st_geometry(US_ALBERS_hostplot))
# US_ALBERS_hostplot <- st_centroid(US_ALBERS)
# ggplot() + geom_sf(data = US_ALBERS, color="black",fill="white", size=0.25)+
#   geom_point(data = US_ALBERS,
#              aes(x = X, y = Y, size = host_biomass, color=factor(Invaded)))
#''''''''''''''''''''''''''''''''''''''''''''
#''''''''''''''''''''''''''''''''''''''''''''



#''''''''''''''''''''''''''''''''''''''''''''
#
# climate normals (ppt, temp)
#
#''''''''''''''''''''''''''''''''''''''''''''
climate_r <- worldclim_global("bio", res=5, "data/climate_normals") # res is in minutes of a degree
climate_r <- climate_r[[c(3,6,12)]]
names(climate_r) <- c("ISOTHERM","MINMOTEMP","MAP")
#plot(climate_r)
#
climate_r <- aggregate(climate_r, fact=4, fun=mean, na.rm=T)
climate_r_albers <- terra::project(climate_r,as.character(crs(US_ALBERS)))
summary(climate_r_albers$MINMOTEMP)

#plot(climate_r_albers[[1]])
#plot(US_ALBERS,add=T)
#
climate_r_albers_r <- stack(climate_r_albers)
summary(climate_r_albers_r$MINMOTEMP)

US_ALBERS$ISOTHERM <- extract(na.omit(climate_r_albers_r[["ISOTHERM"]]), US_ALBERS, fun=mean, na.rm=TRUE)
US_ALBERS$MINMOTEMP <- extract(na.omit(climate_r_albers_r[["MINMOTEMP"]]), US_ALBERS, fun=mean, na.rm=TRUE)
US_ALBERS$MAP <- extract(na.omit(climate_r_albers_r[["MAP"]]), US_ALBERS, fun=mean, na.rm=TRUE)
#''''''''''''''''''''''''''''''''''''''''''''


#''''''''''''''''''''''''''''''''''''''''''''
#
# combine all the dataframes
#
#''''''''''''''''''''''''''''''''''''''''''''
lwd_CPH_fin_df <- fread("data/CPH_data/CPH_LWD.csv", na=".")
nrow(lwd_CPH_fin_df)
table(nchar(lwd_CPH_fin_df$FIPS))
lwd_CPH_fin_df$FIPS <- ifelse(nchar(lwd_CPH_fin_df$FIPS)==4, paste(0,lwd_CPH_fin_df$FIPS,sep=""), paste(lwd_CPH_fin_df$FIPS))
table(nchar(lwd_CPH_fin_df$FIPS))

US_ALBERS_df <- as.data.frame(US_ALBERS)
head(US_ALBERS_df)
US_ALBERS_df <- US_ALBERS_df[,c("FIPS","ltot", "POPESTIMATE2019", "n_campgrounds", "med_income", 
                                "host_biomass", "non_host_biomass", "redbay_biomass","sassafras_biomass",
                                "ISOTHERM","MINMOTEMP","MAP")]
head(US_ALBERS_df)
summary(US_ALBERS_df)


fwrite(US_ALBERS_df, "data/CPH_data/US_ALBERS_df.csv", na=".")
lwd_CPH_fin_df <- merge(x = lwd_CPH_fin_df, y = US_ALBERS_df, by = "FIPS", all.x = TRUE)
nrow(lwd_CPH_fin_df)
head(lwd_CPH_fin_df)
summary(lwd_CPH_fin_df)
summary(US_ALBERS_df)
colnames(US_ALBERS_df)

fwrite(lwd_CPH_fin_df, "data/CPH_data/CPH_LWD_COMBINED.csv", na=".")
#''''''''''''''''''''''''''''''''''''''''''''

beep(1)
