library(readxl) # load in database from Lynne Womack
library(sf)
library(tidyverse)
library(raster) # pointDistance()
library(survival)
library(data.table)
library(sp)
library(geosphere) # radii from a fixed point
library(concaveman) # concave hull polygon
library(geodata) # worldclim_global
library(spdep)
library(beepr)
library(terra)
library(ggpubr)
library(beepr)


# model building
library(MuMIn)
library(viridis)
library(magick)



write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}

theme_clean <- theme_bw()+  
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
