#this script is not super necessary anymore

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

library(readr)
library(dplyr)
library(lubridate)
library(sf)


#load in csv
dataset <- read_csv("IntermediateData/uasdata.csv")

#clean that one weird male
dataset1 <- dataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

#filter out pups, incompletes, and waters
uasdata <- dataset1 %>%
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup")

#this brought the geometry column back
uasdataset <- st_as_sf(
  uasdata,
  coords = c("lon", "lat"),
  crs = "EPSG:26943"   
)


plot(st_geometry(uasdataset))

options(scipen = 999)
dates <- unique(uasdata$date)

years <- substr(dates, 0, 4)

picterra.output <- data.frame()

for(i in 1:length(dates)) {
  filepathadults <- paste0("../1.RawData/Picterra_outputs/", dates[i], "/adults/adults.shp")
  
  filepathpups <- paste0("../1.RawData/Picterra_outputs/", dates[i], "/pups/pups.shp")
  
  picterra.output <- rbind(picterra.output, 
                           cbind(st_read(filepathadults), date = dates[i], year = years[i], class = "adults"),
                           cbind(st_read(filepathpups), date = dates[i], year = years[i], class = "pups"))
  
}

# Convert to NAD83 / California zone 3 -------------------------------------------------------------

st_crs(picterra.output)

picterra.output = st_transform(picterra.output,"EPSG:26943") #transform shape file 

st_crs(picterra.output)

plot(st_geometry(picterra.output))


uas.data <- picterra.output %>%
  filter(class == "adults" & area_m2 < 3.8 &
           class == "adults" &area_m2 > 1.3 |
           class == "pups" & area_m2 <.6 &
           class == "pups" & area_m2 > .14) %>%
  mutate(class = as.factor(case_when(class == "adults" & area_m2 < 2.1 ~ "female", 
                                     class == "adults" & area_m2 >= 2.1 ~ "male",
                                     class == "pups"~ "pup"))) 

table(uas.data$class)



uasdataset <- st_read("IntermediateData/uasdata.shp")

dates <- unique(uasdataset$date)

#calculate density
density.df <- data.frame()

for (i in 1:length(dates)) {
  
  survey.subset <- uasdataset %>%
    filter(date == dates[i],
           age_sex == "female")
  
  seal.centroids <- st_centroid(survey.subset)
  
  
  seal.buffer <- st_buffer(seal.centroids, 10)
  
  
  int <- st_intersects(seal.buffer, seal.centroids) 
  
  survey.subset$density <- lengths(int) - 1
  
  density.df <- rbind(density.df,survey.subset )
}


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


female.density %>%
  filter(year == "2016")

library(ggplot2)
ggplot(data = peak.2016, mapping = aes(x = lon, y = lat, color = density)) +
  geom_point(alpha = 0.7) +
  coord_cartesian(ylim = c(37.11264, 37.114), xlim = c(-122.3295, -122.328)) +
  
  
  write.csv(female.density, "female.density.csv")

#save the density.df thing with this geometry
write.csv(density.df, "Uasdata with density")



