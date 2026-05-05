#Density Gradient Figure - cite Molly McEntee 2026 (unpublished)
#------------------------------------------------------------------------------
#set up
library(dplyr)
library(lubridate)
library(sf)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

uasdataset <- st_read("IntermediateData/uasdata.full.shp")

#------------------------------------------------------------------------------
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

plot(st_geometry(survey.subset[12,]), add = TRUE, col = "blue")

c <- int[[12]] 

nearest.neighbors <- survey.subset[c,]

plot(st_geometry(nearest.neighbors), add = TRUE)

#------------------------------------------------------------------------------
#Get min and max coords so you can section off the larger area (colony) into smaller areas

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

#This limits it to south point (2016)
plot(peak.2016["density"],
     xlim = c(1836725 +680, 1836725 + 820), ylim = c(569768-150, 569768-80),
     main = "2016 South Point Density Gradient",
     key.pos = 4,
     key.width = lcm(1.6),
)

# create file path to save it
outfile <- file.path(getwd(), "SP_peak_density_2016.png")

#save it as a png
png(outfile, width = 8, height = 6, units = "in", res = 600)

#plot 
plot(peak.2016["density"],
     xlim = c(1836725 + 680, 1836725 + 820),
     ylim = c(569768 - 150, 569768 - 80),
     main = "South Point Density Gradient (2016)",
     key.pos = 4,
     key.width = lcm(1.6),
)

# lose device (VERY IMPORTANT for png saving)
dev.off()

#check that it saved (should say TRUE)
file.exists(outfile)
#---------------------------
#this gets you full colony (This method is better for smaller areas though so this may not be necessary)
plot(peak.2016["density"],
     xlim = c(1836832, 1837680), ylim = c(569614.5, 570193.2))

#------------------------------------------------------------------------------
#this is 2018 specifically
peak.2018 <-uasdataset %>%
  filter(date == "20180131")

#this gives just south point
plot(peak.2018["density"],
     xlim = c(1836725 +680, 1836725 + 820), ylim = c(569768-150, 569768-80))
#------------------------------------------------------------------------------
#this is 2019 specifically
peak.2019 <-uasdataset %>%
  filter(date == "20190129")

#this gives just south point
plot(peak.2019["density"],
     xlim = c(1836725 +680, 1836725 + 820), ylim = c(569768-150, 569768-80))
#------------------------------------------------------------------------------
#repeat for every year from here and reference 2016 section for saving method.

#end.


