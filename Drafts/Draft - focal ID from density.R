#here I will try to identify focals through density assessment

library(readr)
library(sf)
library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

uasdataset <- st_read("IntermediateData/uasdata.full.shp")

uasdata <- uasdataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

#---------------------THERE IS A PROBLEM HERE! (the problem is that i'm confused what the columns female99 and male90 mean)
#this gives a table of how many females and males we have per year.
table(uasdata$year, uasdata$age_sex)

counts <- table(uasdata$year, uasdata$age_sex)

colMeans(counts)

uasdata_full <- uasdata %>%
  group_by(year, age_sex) %>%
  mutate(
    female_99 = quantile(density[age_sex == "female"], 0.99, na.rm = TRUE),
    male_90   = quantile(density[age_sex == "male"],   0.90, na.rm = TRUE),
    status = case_when(
      age_sex == "female" & density >= female_99 ~ "focal female",
      age_sex == "male"   & density >= male_90   ~ "focal male",
      age_sex == "female" ~ "female",
      age_sex == "male"   ~ "male"
    )
  ) %>%
  ungroup()

table(uasdata_full$year, uasdata_full$status)

#filtering out the focals so I can see their densities
focals <- uasdata_full %>%
  filter(status %in% c("focal female", "focal male"))

#checking how many focal females and focal males i have
table(focals$status)

#tells me the min and max of density for my focals
focals %>%
  group_by(status) %>%
  summarise(
    min_density = min(density, na.rm = TRUE),
    max_density = max(density, na.rm = TRUE)
  )
#--------------------------------------------------------------------------------------
#here I will try to get the distance from the closest focal and closest focal index
#this gives us an index for every row so now we have an ID system
uasdata1 <- uasdata_full %>%
  mutate(index = row_number())

#converts to sf
uasdata1_sf <- st_as_sf(uasdata1, coords = c("X", "Y"), crs = 26943)

st_crs(uasdata1_sf)

#separate the statuses
females <- uasdata1_sf %>% filter(status == "female")
focal_females   <- uasdata1_sf %>% filter(status == "focal female")
focal_males     <- uasdata1_sf %>% filter(status == "focal male")

table(focal_females$year)
#here we calculate the closest focal female for every regular female
females <- females %>%
  group_by(year) %>%
  group_modify(~{
    
    foc <- focal_females %>% filter(year == .y$year)
    if (nrow(foc) == 0) return(.x)   # skip if none that year
    
    idx <- st_nearest_feature(.x, foc)
    
    .x$closest_focal_fem_index <- foc$index[idx]
    .x$distance_focal_female   <- st_distance(.x, foc[idx, ], by_element = TRUE)
    
    .x
  }) %>%
  ungroup()

#now i will try to graph it to see what's going on
links <- females %>%
  left_join(
    focal_females %>%
      st_drop_geometry() %>%
      select(index, X, Y),
    by = c("closest_focal_fem_index" = "index"),
    suffix = c("_reg", "_foc")
  ) %>%
  filter(!is.na(closest_focal_fem_index))

ggplot() +
  # Draw lines connecting each regular female to its focal female
  geom_segment(data = links %>% filter(year == 2016),
               aes(x = X_foc, y = Y_foc, xend = X_reg, yend = Y_reg),
               color = "gray70", size = 0.5) +
  
  # Plot focal females
  geom_point(data = focal_females %>% st_drop_geometry(),
             aes(x = X, y = Y),
             color = "red", size = 3) +
  
  # Plot regular females
  geom_point(data = females %>% st_drop_geometry(),
             aes(x = X, y = Y),
             color = "blue", size = 1.5) +
  
  labs(x = "X (meters)", y = "Y (meters)",
       title = "Focal Females (red) and Associated Regular Females (blue)") +
  theme_minimal()
