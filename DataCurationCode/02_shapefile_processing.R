
# Load libraries and data -------------------------------------------------
getwd()
setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

library(dplyr)
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

#read this in when you want to get the plot without having to run the length vs width code again!
uasdataset <- st_read("IntermediateData/uasdata.full.shp")
#now get into density
density.df <- data.frame()

for (i in 1:length(date)) {
  
  survey.subset <- uasdataset %>%
    filter(date == date[i],
           is.na(p_complete) | !p_complete %in% c("N", "water"),
           age_sex !="pup")
  
  seal.centroids <- st_centroid(survey.subset)
  
  
  seal.buffer <- st_buffer(seal.centroids, 10)
  
  
  int <- st_intersects(seal.buffer, seal.centroids) 
  
  survey.subset$density <- lengths(int) - 1
  
  density.df <- rbind(density.df,survey.subset )
}

#save the density.df thing as shp files
st_write(density.df, "uasdata.density.shp")

plot(st_geometry(seal.buffer[12,]))

plot(st_geometry(survey.subset[12,]), add = TRUE, col = "red")

c <- int[[12]] 

nearest.neighbors <- survey.subset[c,]

plot(st_geometry(nearest.neighbors), add = TRUE)

#------------------------------------------------------------------------------
#trying to figure out the xlim and ylim to use to get full colony

df2016 <- density.df[density.df$year == 2016, ]

xmin <- min(df2016$X, na.rm = TRUE)
xmax <- max(df2016$X, na.rm = TRUE)

ymin <- min(df2016$Y, na.rm = TRUE)
ymax <- max(df2016$Y, na.rm = TRUE)

xmin; xmax; ymin; ymax
#------------------------------------------------------------------------------

peak.2016 <-density.df %>%
  filter(date == "20160125")

#this gives just SP
plot(peak.2016["density"],
     xlim = c(1836725 +680, 1836725 + 820), ylim = c(569768-150, 569768-80))

#this gets you full colony
plot(peak.2016["density"],
     xlim = c(1836832, 1837680), ylim = c(569614.5, 570193.2))

png("densityplot.png", width = 800, height = 600)  # pixels
plot(x, y)
dev.off()

female.density <-density.df %>%
  mutate(date = as_date(date))


hist(female.density$density)


peak.2016 <-female.density %>%
  filter(year == "2016")


#also save the uas.dataset with the geometry (without the density)
write_csv(uas.dataset, "uas.dataset.csv")

#---Here I am trying to extract coords - THIS GAVE ME A SEPARATE FILE WITH COORDS (doesnt have density)
library(sf)

uas.dataset$centroid <- st_centroid(uas.dataset$geometry)

uas.dataset.coords <- uas.dataset %>%
  mutate(
    centroid = st_centroid(geometry),
    lon = st_coordinates(centroid)[, 1],
    lat = st_coordinates(centroid)[, 2]
  )

write_csv(uas.dataset.coords, "uas.dataset.coords.csv")

#ok so now you have mulitple files in your intermediate data folder, "uas.dataset.csv", "uas.dataset.coords.csv", and the shp files for density

uasdataset.density <- st_read("IntermediateData/uasdata.density.shp")
#ok so the shp file with the densities seems to have multipolygons because of the way we did the density calculations? grouping the seals into that circle

#so this extracted coords for the centroids. only thing is im not sure how to validate that these are the right coords since the geometry was MULTIPOLYGON and just POLYGON.
#THIS GIVES ME COORDS WITHIN THE SAME DF AS THE DENSITY CALCULATIONS!
st_geometry_type(uasdataset.density)

centroids <- st_centroid(uasdataset.density)
coords <- st_coordinates(centroids)

uasdataset.density$X <- coords[,1]
uasdataset.density$Y <- coords[,2]

#save the shp file (includes morphs, density, coords, and geometry)
st_write(uasdataset.density, "uasdata.full.shp")

