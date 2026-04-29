#Using focal males to calculate harem size and harem distance averages

#first you need to start with the same set up as the original Focal Male ID(density) code
library(readr)
library(sf)
library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(gt)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

uasdataset <- st_read("IntermediateData/uasdata.full.shp")

uasdata <- uasdataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

#-------------------------------------------------------------------------------
#population wide summary numbers
#males
uasdata %>%
  sf::st_drop_geometry() %>%
  filter(age_sex == "male") %>%
  count(year)

uasdata %>%
  sf::st_drop_geometry() %>%
  filter(age_sex == "male") %>%
  count(year) %>% 
  summarise(mean_count = mean(n))

#females
uasdata %>% 
  sf::st_drop_geometry() %>%
  filter(age_sex =="female") %>% 
  count(year)

uasdata %>% 
  sf::st_drop_geometry() %>%
  filter(age_sex =="female") %>% 
  count(year) %>% 
  summarise(mean_count = mean(n))

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

females <- females %>%
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
links <- females %>%
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
#here is the code for harem size (# of females)  
harem_size <- links %>%
  count(year, closest_focal_male_index, region)

#no get a average harem size per year
average_harem_size <- harem_size %>%
  group_by(year, region) %>%
  summarise(avg_harem_size = mean(n))

harem_summary <- harem_size %>%
  group_by(year, region) %>%
  summarise(
    avg_harem_size = mean(n),
    min_harem_size = min(n),
    max_harem_size = max(n)
  )

harem_summary_table <- harem_summary %>%
  gt()

gtsave(harem_summary_table, "harem_summary_table.png")
#-------------------------------------------------------------------------------
#here is the code for average harem distance (distance from focal male)

#first you need to make sure the distance is a numeric not categorical value
str(links$distance_focal_male)

harem_spatial <- links %>%
  group_by(year, region, closest_focal_male_index) %>%
  summarise(
    mean_distance = mean(distance_focal_male, na.rm = TRUE),
    min_distance = min(distance_focal_male, na.rm = TRUE),
    max_distance = max(distance_focal_male, na.rm = TRUE)
  )
  
harem_spatial_summary <- links %>%
  group_by(year, region, closest_focal_male_index) %>%
  summarise(
    mean_distance = mean(distance_focal_male, na.rm = TRUE),
    min_distance = min(distance_focal_male, na.rm = TRUE),
    max_distance = max(distance_focal_male, na.rm = TRUE)
  ) %>%
  group_by(year, region) %>%
  summarise(
    avg_mean = mean(mean_distance, na.rm = TRUE),
    avg_min = min(min_distance, na.rm = TRUE),
    avg_max = max(max_distance, na.rm = TRUE),
    .groups = "drop"
  )

harem_spatial_summary_table <- harem_spatial_summary %>%
  gt(groupname_col = "year")

gtsave(harem_spatial_summary_table, "harem_spatial_summary_table.png")

#-------------------------------------------------------------------------------
#number of harems

focal_summary <- links %>%
  group_by(year, region) %>%
  summarise(
    n_focal_males = n_distinct(closest_focal_male_index),
    .groups = "drop"
  )

focal_table <- focal_summary %>%
  arrange(year, region) %>%
  gt(groupname_col = "year") %>%
  tab_header(
    title = "Number of Focal Males by Year and Region"
  ) %>%
  cols_label(
    region = "Region",
    n_focal_males = "Number of Focal Males"
  )

focal_table
gtsave(focal_table, "focal_male_summary.png")


#-------------------------------------------------------------------------------











