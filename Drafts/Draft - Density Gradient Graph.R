#here is the code for the density gradient plot

library(dplyr)
library(lubridate)
library(sf)

uasdataset <- st_read("IntermediateData/uasdata.full.shp")

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
#this is 2016 specifically

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
#------------------------------------------------------------------------------
#this is 2017 specifically
peak.2018 <-uasdataset %>%
  filter(date == "20180131")

#this gives just SP
plot(peak.2018["density"],
     xlim = c(1836725 +680, 1836725 + 820), ylim = c(569768-150, 569768-80))
#------------------------------------------------------------------------------
#this is 2019 specifically
peak.2019 <-uasdataset %>%
  filter(date == "20190129")

#this gives just SP
plot(peak.2019["density"],
     xlim = c(1836725 +680, 1836725 + 820), ylim = c(569768-150, 569768-80))




