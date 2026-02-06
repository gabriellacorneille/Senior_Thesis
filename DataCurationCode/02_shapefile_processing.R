
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


#now get into density
density.df <- data.frame()

for (i in 1:length(dates)) {
  
  survey.subset <- uas.dataset %>%
    filter(date == dates[i],
           is.na(p_complete) | !p_complete %in% c("N", "water"),
           age_sex !="pup")
  
  seal.centroids <- st_centroid(survey.subset)
  
  
  seal.buffer <- st_buffer(seal.centroids, 10)
  
  
  int <- st_intersects(seal.buffer, seal.centroids) 
  
  survey.subset$density <- lengths(int) - 1
  
  density.df <- rbind(density.df,survey.subset )
}

st_write(density.df, "uasdata.density.shp")

plot(st_geometry(seal.buffer[12,]))

plot(st_geometry(survey.subset[12,]), add = TRUE, col = "red")

c <- int[[12]] 

nearest.neighbors <- survey.subset[c,]

plot(st_geometry(nearest.neighbors), add = TRUE)

peak.2016 <-density.df %>%
  filter(date == "20160125")


plot(peak.2016["density"],
     xlim = c(1836725 +680, 1836725 + 820), ylim = c(569768-150, 569768-80))


female.density <-density.df %>%
  mutate(date = as_date(date))


hist(female.density$density)


peak.2016 <-female.density %>%
  filter(year == "2016")

library(ggplot2)
ggplot(data = peak.2016, mapping = aes(x = lon, y = lat, color = density)) +
  geom_point(alpha = 0.7) +
  coord_cartesian(ylim = c(37.11264, 37.114), xlim = c(-122.3295, -122.328)) +
  
  
  write.csv(female.density, "female.density.csv")

#save the density.df thing with this geometry
st_write(density.df, "my_file.shp")


#transforms data into lat/lon
#uas.dataset.check <- st_transform(uas.dataset, 4326)

#here there is an issue with overlapping/intersecting polygons that s2 can't calculate centroids for so we're going to turn that off for a second
#sf_use_s2(FALSE)

# calculate the centroids of each seal polygon
#centroids <- st_centroid(uas.dataset)

# gets the coordinates for each centroid of a seal polygon
#coords <- st_coordinates(centroids)

# Add lat and lon columns to the dataset
#uas.dataset$lon <- coords[, 1]
#uas.dataset$lat <- coords[, 2]

#this will turn s2 back on
#sf_use_s2(TRUE)

st_write(uas.dataset, "IntermediateData/uasdata.shp")

# Save outputs ------------------------------------------------------------

uas.data <- uas.dataset |>
  st_drop_geometry()

write.csv(uas.data, "IntermediateData/uasdata.csv",
          row.names = FALSE)

R.version.string

