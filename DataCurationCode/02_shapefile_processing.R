
# Load libraries and data -------------------------------------------------
getwd()
setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

library(dplyr)
install.packages("lubridate")
library(lubridate)
library(sf)

options(scipen = 999)

dates <- list.files("RawData/")




years <- substr(dates, 0, 4)

uas.dataset <- data.frame()

for(i in 1:length(dates)) {
  filepath <- paste0("RawData/", dates[i], "/polygons.shp")


  uas.dataset <- rbind(uas.dataset, 
                       cbind(st_read(filepath), date = dates[i], year = years[i]))
  
}



# Convert to NAD83 / California zone 3 -------------------------------------------------------------


st_crs(uas.dataset)

uas.dataset = st_transform(uas.dataset,"EPSG:26943") #transform shape file 

st_crs(uas.dataset)

plot(st_geometry(uas.dataset))

# Calculate area ----------------------------------------------------------

uas.dataset$area <- as.numeric(st_area(uas.dataset))


# Calculate seal length and width -----------------------------------------

source("DataCurationCode/Functions/1_Seal_Volume_Function_MHM_fixed.R")



uas.dataset$length <- NA
uas.dataset$width <- NA



for(j in 1:length(st_geometry(uas.dataset))){
  pol1 = uas.dataset[j,]
  lens = curved_length_vol(pol1, plt = TRUE)
  uas.dataset[j,7:8] = lens[1:2]
  print(j)
} 

library(sf)

#transforms data in lat/lon
uas.dataset <- st_transform(uas.dataset, 4326)

#here there is an issue with overlapping/intersecting polygons that s2 can't calculate centroids for so we're going to turn that off for a second
sf_use_s2(FALSE)

# calculate the centroids of each seal polygon
centroids <- st_centroid(uas.dataset)

# gets the coordinates for each centroid of a seal polygon
coords <- st_coordinates(centroids)

# Add lat and lon columns to the dataset
uas.dataset$lon <- coords[, 1]
uas.dataset$lat <- coords[, 2]

#this will turn s2 back on
sf_use_s2(TRUE)



# Save outputs ------------------------------------------------------------

uas.data <- uas.dataset |>
  st_drop_geometry()

write.csv(uas.data, "IntermediateData/uasdata.csv",
          row.names = FALSE)

