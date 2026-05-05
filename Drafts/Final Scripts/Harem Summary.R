#Harem Summary (includes population summary, harem size, harem spatial extent, and number of harems)
#uses focal males as harem centers

#set up
library(readr)
library(sf)
library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)
library(gt)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

uasdataset <- st_read("IntermediateData/Final Data Files/uasdata.full.shp")

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
    if (nrow(foc) == 0) return(.x)  
    
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

#now get a average harem size per year
harem_summary <- harem_size %>%
  group_by(year, region) %>%
  summarise(
    "Average Harem Size" = mean(n),
    "Minimum Harem Size" = min(n),
    "Maximum Harem Size" = max(n)
  ) %>%
  mutate(region = tools::toTitleCase(region))

harem_summary <- harem_summary %>%
  rename(Region = region)

harem_summary_table <- harem_summary %>%
  gt(groupname_col = "year") %>%
  tab_style(
    style = list(
      cell_fill(color = "lightgray"),
      cell_text(weight = "bold")
    ),
    locations = cells_row_groups()
  ) %>%
  cols_label(
    `Average Harem Size` = "Average\nHarem Size",
    `Minimum Harem Size` = "Minimum\nHarem Size",
    `Maximum Harem Size` = "Maximum\nHarem Size"
  ) %>%
  cols_width(
    `Average Harem Size` ~ px(90),
    `Minimum Harem Size` ~ px(90),
    `Maximum Harem Size` ~ px(90)
  )%>%
  fmt_number(
    columns = c(`Average Harem Size`),
    decimals = 0 
  )

#see table
harem_summary_table 

#save
gtsave(harem_summary_table, "TablesFigures/Summary Tables/harem_summary_table.png")
#-------------------------------------------------------------------------------
#here is the code for average harem distance (distance from focal male)

#check value type
str(links$distance_focal_male)
  
harem_spatial_summary <- links %>%
  group_by(year, region, closest_focal_male_index) %>%
  summarise(
    mean_distance = mean(as.numeric(distance_focal_male), na.rm = TRUE),
    min_distance = min(as.numeric(distance_focal_male), na.rm = TRUE),
    max_distance = max(as.numeric(distance_focal_male), na.rm = TRUE)
  ) %>%
  group_by(year, region) %>%
  summarise(
    "Average Radius (m)" = mean(mean_distance, na.rm = TRUE),
    "Minimum Radius (m)" = min(min_distance, na.rm = TRUE),
    "Maximum Radius (m)" = max(max_distance, na.rm = TRUE),
    .groups = "drop"
  )%>%
  mutate(region = tools::toTitleCase(region))

harem_spatial_summary <- harem_spatial_summary %>%
  rename(Region = region)

harem_spatial_summary_table <- harem_spatial_summary %>%
  gt(groupname_col = "year") %>% 
  tab_style(
    style = list(
      cell_fill(color = "lightgray"),
      cell_text(weight = "bold")
    ),
    locations = cells_row_groups()
  ) %>%
  cols_label(
    `Average Radius (m)` = "Average\nRadius (m)",
    `Minimum Radius (m)` = "Minimum\nRadius (m)",
    `Maximum Radius (m)` = "Maximum\nRadius (m)"
  ) %>%
  cols_width(
    `Average Radius (m)` ~ px(90),
    `Minimum Radius (m)` ~ px(90),
    `Maximum Radius (m)` ~ px(90)
  )%>%
  fmt_number(
    columns = c(`Average Radius (m)`, `Minimum Radius (m)`, `Maximum Radius (m)`),
    decimals = 2 
  )

#see table
harem_spatial_summary_table

#save
gtsave(harem_spatial_summary_table, "TablesFigures/Summary Tables/harem_spatial_table.png")

#-------------------------------------------------------------------------------
#number of harems

focal_summary <- links %>%
  group_by(year, region) %>%
  summarise(
    n_focal_males = n_distinct(closest_focal_male_index),
    .groups = "drop"
  ) %>%
  mutate(region = tools::toTitleCase(region))

focal_table <- focal_summary %>%
  arrange(year, region) %>%
  gt(groupname_col = "year") %>%
  tab_style(
    style = list(
      cell_fill(color = "lightgray"),
      cell_text(weight = "bold")
    ),
    locations = cells_row_groups()
  ) %>%
  cols_label(
    region = "Region",
    n_focal_males = "Number of \nFocal Males"
  ) %>%
  cols_width(
    `n_focal_males` ~ px(90))

#see table
focal_table

#save
gtsave(focal_table, "TablesFigures/Summary Tables/focal_male.png")

#-------------------------------------------------------------------------------
#end.









