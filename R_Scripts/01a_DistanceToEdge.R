# Calculates distance from edge for each nest, 2/4/2019, Shane Somers
# Based on rgis.r file, values are distance from each point (nest) to nearest polygon edge (site boundary)

#rm(list=ls())


# install packages
library(sf)
library(ggplot2)
library(geosphere)
library(dplyr)

# read in nest box shapefile
nestboxes <- st_read("Data/Shapefiles/Nestboxes_2016.shp") # proj=tmerc, units seem to be meters

# read in site perimeters
#sites <- st_read("Data/Shapefiles/Sites/Study_sites.shp") # proj=tmerc, units seem to be meters
#sites <- st_read("Data/Shapefiles/Sites_relaxed/Study_sites_relaxed-polygon.shp") # proj=tmerc, units seem to be meters
sites <- st_read("Data/Shapefiles/Sites_strict/Study_sites_strict.shp")

# init vector and then apply or loop with dictionary lookup
sites.abbrev <- vector("character", nrow(nestboxes))

# correct nestboxes location variable spelling
## for some reason was giving an error if i used double[] in expr part of if statement
levels(nestboxes$Location) <- c(levels(nestboxes$Location), c("Ballinphellic", "Castle Bernard", "Dukes Wood",  "Garrettstown")) ## need to add new factor level first
for (i in 1:nrow(nestboxes)) {
  if (nestboxes[[i,"Location"]]=="Ballinphelic") {nestboxes[i,"Location"] <- "Ballinphellic"}
  if (nestboxes[[i,"Location"]]=="CastleBernar"){nestboxes[i,"Location"] <- "Castle Bernard"}
  if (nestboxes[[i,"Location"]]=="DukesWood"){nestboxes[i,"Location"] <- "Dukes Wood"}
  if (nestboxes[[i,"Location"]]=="Garretstown"){nestboxes[i,"Location"] <- "Garrettstown"}
}

#dput(unique(nestboxes$Location))
## NOTE: Some locations misspelled in nestboxes DF
#sites.names <- c("Ballinphelic", "Carrigeen", "CastleBernar", 
#                "DukesWood", "Dunderrow", "Farran", "Garretstown", "Innishannon", 
#               "Kilbrittain", "Lissarda", "Piercetown", "Shippool")

# corrected names
sites.names <- c("Ballinphellic", "Carrigeen", "Castle Bernard", 
                 "Dukes Wood", "Dunderrow", "Farran", "Garrettstown", "Innishannon", 
                 "Kilbrittain", "Lissarda", "Piercetown", "Shippool")

sites.abb <- c("BP","CG","CB","DW","DD","FR","GT","IN","KB","LS","PT","SP")
names(sites.abb) <- sites.names

#site.dict <- cbind(sites.names,sites.abb)

for (i in 1:nrow(nestboxes)) {
  site.holder <- as.character(nestboxes$Location[i])
  sites.abbrev[i] <- sites.abb[site.holder]
  #print(site.holder)
}

nestboxes <- cbind(nestboxes,sites.abbrev)

# create and attach nest ID column
nest.ID <- paste(nestboxes$sites.abbrev, nestboxes$nestbox, sep = "")
nestboxes <- cbind(nestboxes, nest.ID)


# Transform sf objects to EPSG:4326
pointsWgs84.sf <- st_transform(nestboxes, crs = 4326) # as proj=longlat units are implicitly degrees
lineWgs84.sf <- st_transform(sites, crs = 4326)       # as proj=longlat units are implicitly degrees

# try calc dist with original projection which has units in meters
# dist <- dist2Line(st_coordinates(nestboxes)[,1:2],st_coordinates(sites)[,1:2])

# calc distance of each point to nearest line ## NEED TO VERIFY BY PLOTTING ## VERIFY UNITS
dist <- dist2Line(st_coordinates(pointsWgs84.sf)[,1:2],st_coordinates(lineWgs84.sf)[,1:2])

## FIGURE OUT WHAT UNITS DISTANCE BEING MEASURED IN 
## can verify using Ruler tool on Google earth Pro
## confirmed  units are meters

# create df w/ nest.ID and distance to edge, also need to give nest.ID column correct name
## note: lat lon give position of nearest coordinate on the line for each nest
nestbox.edge <- as.data.frame(cbind(as.character(nestboxes$nest.ID), dist)) %>% 
  mutate_at(c("distance", "lon", "lat"), as.numeric) %>% dplyr::rename(nest = V1, DistanceToEdge = distance)

#head(nestbox.edge)

# remove unneccessary variables and data
rm(list = c("dist", "lineWgs84.sf", "nest.ID", 
            "nestboxes", "pointsWgs84.sf", "site.holder", "sites","i", "sites.abb", 
            "sites.abbrev", "sites.names"))
