#Focal Male ID from density
#-------------------------------------------------------------------------------

#set up
library(readr)
library(sf)
library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

uasdataset <- st_read("IntermediateData/Final Data Files/uasdata.full.shp")

uasdata <- uasdataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

#-------------------------------------------------------------------------------
#separate colony into regions
#find the min and max for lon (x) and lat (y)
lat_lon <- uasdata

xmin <- min(lat_lon$X, na.rm = TRUE)
xmax <- max(lat_lon$X, na.rm = TRUE)

ymin <- min(lat_lon$Y, na.rm = TRUE)
ymax <- max(lat_lon$Y, na.rm = TRUE)

xmin; xmax; ymin; ymax
#xmin=lonmin, xmax=lonmax, ymin=latmin, ymax=latmax

#this creates a new column "region" that has south mid and north
uasdata_with_regions <- uasdata %>%
  mutate(region = case_when(
    Y >= 569593 & Y < 569880 ~ "south",
    Y >= 569875 & Y < 570197.25 ~ "mid",
    Y >= 570197.25 & Y < 571025 ~ "north",
    TRUE ~ NA_character_   # anything outside ranges
  ))

#check
ggplot(uasdata_with_regions, aes(X, Y, color = region)) + geom_point()

ggsave("colony regionalization.png",
       plot = colony_regionalization,#rename the plot above with a title so you can save it
       width = 8,
       height = 6,
       dpi = 600)

#-------------------------------------------------------------------------------
#create a focal status column
uasdata_full <- uasdata_with_regions %>%
  group_by(year, region, age_sex) %>%
  mutate(
    threshold = case_when(
      age_sex == "female" ~ quantile(density, 0.99, na.rm = TRUE),#sets the threshold for focal females
      age_sex == "male"   ~ quantile(density, 0.90, na.rm = TRUE)#sets the threshold for focal males
    ),
    status = case_when(
      age_sex == "female" & density >= threshold ~ "focal female",
      age_sex == "male"   & density >= threshold ~ "focal male",
      age_sex == "female" ~ "female",
      age_sex == "male"   ~ "male"
    )
  ) %>%
  ungroup()

#filtering out the focals so I can see their densities
focals <- uasdata_full %>%
  filter(status %in% c("focal female", "focal male"))

#checking how many focal females and focal males i have
table(focals$status, focals$year, focals$region)

#tells me the min and max of density for my focals
focals %>%
  group_by(status) %>%
  summarise(
    min_density = min(density, na.rm = TRUE),
    max_density = max(density, na.rm = TRUE)
  )

#-------------------------------------------------------------------------------
#this gives us an index for every row so now we have an ID system
uasdata1 <- uasdata_full %>%
  mutate(index = row_number())

#converts to sf
uasdata1_sf <- st_as_sf(uasdata1, coords = c("X", "Y"), crs = 26943)

st_crs(uasdata1_sf)

#separate the statuses
males <- uasdata1_sf %>% filter(status == "male")
females <- uasdata1_sf %>% filter(status == "female")
focal_females   <- uasdata1_sf %>% filter(status == "focal female")
focal_males     <- uasdata1_sf %>% filter(status == "focal male")

#check
table(focal_males$year)

#-------------------------------------------------------------------------------
#calculate closest focal male for every regular female
females2 <- females %>%
  group_by(year, region) %>%
  group_modify(~{
    
    foc <- focal_males %>% filter(year == .y$year, region == .y$region)
    if (nrow(foc) == 0) return(.x)   # skip if none that year
    
    idx <- st_nearest_feature(.x, foc)
    
    .x$closest_focal_male_index <- foc$index[idx]
    .x$distance_focal_male <- st_distance(.x, foc[idx, ], by_element = TRUE)
    
    .x
  }) %>%
  ungroup()

#recombine data so it's all in one place
links2 <- females2 %>%
  left_join(
    focal_males %>%
      st_drop_geometry() %>%
      select(index, X, Y),
    by = c("closest_focal_male_index" = "index"),
    suffix = c("_reg", "_foc")
  ) %>%
  filter(!is.na(closest_focal_male_index))

#so within this "links" df we now have all the original data and the new data (region, lon(x), lat (y), index, closest focal male index, distance from focal male, focal male coords)

#-------------------------------------------------------------------------------
#Below is the graph for 2016, for graphing all other years, change the year and region based on what you want to look at!
#change the region filter (either south, mid, north) or just remove the region filter to get the whole colony

#graphing 2016 - change the region filter (either south, mid, north) or just remove the region filter to get the whole colony
fig_2_pt2 <- ggplot() +
  # Draw lines connecting each regular female to its focal male
  geom_segment(data = links %>% filter(year == 2016, region == "south"),
               aes(x = X_foc, y = Y_foc, xend = X_reg, yend = Y_reg),
               color = "grey", size = 0.5) +
  
  #plot regular females
  geom_point(data = females %>% 
               filter(year == 2016, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y, color = "Female", shape = "Female"),
             size = 1.5) +
  
  # Plot focal females
  geom_point(data = focal_females %>% 
               filter(year == 2016, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y, color = "Female", shape = "Female"),
             size = 1.5) +
  
  #plot regular males
  geom_point(data = males %>% 
               filter(year == 2016, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y, color = "Regular Male", shape = "Regular Male"),
             size = 1.5) +
  
  #Plot focal males
  geom_point(data = focal_males %>% 
               filter(year == 2016, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y, color = "Focal Male", shape = "Focal Male"),
             size = 4) +
  #create a color legend
  scale_color_manual(
    name = "Legend",
    values = c(
      "Female" = "royalblue",
      "Regular Male" = "#190933",
      "Focal Male" = "red"
    )
  ) +
  #create a shape legend
  scale_shape_manual(
    name = "Legend",
    values = c(
      "Focal Male" = 17,
      "Regular Male" = 17,
      "Female" = 16
    )
  ) +
  
  labs(x = "Longitude (meters)", y = "Latitude (meters)") +
       #title = "Focal Males (red) and Associated Regular Females (blue)") +
  theme_bw() ; fig_2_pt2

ggsave("TablesFigures/focal ID -> density/focal male/2016/2016_focal_male_south.png",
       width = 8,
       height = 6,
       dpi = 600)

#-------------------------------------------------------------------------------
#end-------------------------------------------------------------------------------


