library(usethis)
usethis::use_git()
usethis::use_git_remote("origin", url = NULL, overwrite = TRUE)


## install if needed (do this exactly once):
## install.packages("usethis")


use_git_config(user.name = "sfward", user.email = "ward.1792@osu.edu")

usethis::create_from_github(
  "https://github.com/sfward/laurel-wilt-spread.git",
  destdir = "C:/Users/sw2442/Desktop/practice"
)

?gh_token_help()
# Libraries --------------------------------------------
#
#
#---
library(sf)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

getwd()

# Effective range radius --------------------------------------------
#
#
#---
spread_ALBERS$county_area <- st_area(spread_ALBERS)

annual_spread <- data.frame(year=2004:2021, area_invaded_km2=NA)
years.to.plot <- sort(unique(spread_ALBERS$YrInv))
i <- years.to.plot[1]

for(i in years.to.plot){
  curr.invaded <- spread_ALBERS[which(spread_ALBERS$YrInv %in% 2004:i),]
  annual_spread[which(annual_spread$year %in% i), "year"] <- i
  annual_spread[which(annual_spread$year==i), "area_invaded_km2"] <- sum(st_area(curr.invaded)/1e6)
}

annual_spread$radius_2 <- sqrt(2*annual_spread$area_invaded_km2/pi)

plot(radius_2~year, data=annual_spread)

fit_ERR <- lm(radius_2 ~ year, data=annual_spread)
summary(fit_ERR)
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

