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
#-------------------------------------------------------------------------------
#here I will try to split the colony by latitude

#find the min and max for latitude or Y
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
colony_regionalization <- ggplot(uasdata_with_regions, aes(X, Y, color = region)) + geom_point()

ggsave("colony regionalization.png",
       plot = colony_regionalization,
       width = 8,
       height = 6,
       dpi = 600)
#-------------------------------------------------------------------------------

#this creates the focal status column!
uasdata_full <- uasdata_with_regions %>%
  group_by(year, region, age_sex) %>%
  mutate(
    threshold = case_when(
      age_sex == "female" ~ quantile(density, 0.99, na.rm = TRUE),
      age_sex == "male"   ~ quantile(density, 0.90, na.rm = TRUE)
    ),
    status = case_when(
      age_sex == "female" & density >= threshold ~ "focal female",
      age_sex == "male"   & density >= threshold ~ "focal male",
      age_sex == "female" ~ "female",
      age_sex == "male"   ~ "male"
    )
  ) %>%
  ungroup()
#ok this gives me good numbers for north and south focals, but confusing numbers for mid focals

#check
table(uasdata_full$region, uasdata_full$year, uasdata_full$status)
#-------------------------------------------------------------------------------
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
#here I will try to get the distance from the closest focal and closest focal index
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
table(focal_females$year)

#here we calculate the closest focal female for every regular female
females <- females %>%
  group_by(year, region) %>%
  group_modify(~{
    
    foc <- focal_females %>% filter(year == .y$year, region == .y$region)
    if (nrow(foc) == 0) return(.x)   # skip if none that year
    
    idx <- st_nearest_feature(.x, foc)
    
    .x$closest_focal_fem_index <- foc$index[idx]
    .x$distance_focal_female   <- st_distance(.x, foc[idx, ], by_element = TRUE)
    
    .x
  }) %>%
  ungroup()
#running this helped me realize that i need to split up the colony by latitude. 
#GOOD NEWS THOUGH: THE ASSOCIATION METHOD WORKED SO CONTINUE TO USE THIS CODE!

#-------------------------------------------------------------------------------
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

#graphing 2016 
ggplot() +
  # Draw lines connecting each regular female to its focal female
  geom_segment(data = links %>% filter(year == 2016),
               aes(x = X_foc, y = Y_foc, xend = X_reg, yend = Y_reg),
               color = "gray70", size = 0.5) +
  
  # Plot regular females
  geom_point(data = females %>% 
               filter(year == 2016) %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "blue", size = 1.5) +
 
   # Plot focal females
  geom_point(data = focal_females %>% 
               filter(year == 2016) %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "red", size = 3) +
  
  labs(x = "X (meters)", y = "Y (meters)",
       title = "Focal Females (red) and Associated Regular Females (blue)") +
  theme_minimal()

#graphing 2018 (remove region=south if you want full colony)
ggplot() +
  # Draw lines connecting each regular female to its focal female
  geom_segment(data = links %>% filter(year == 2018, region == "south"),
               aes(x = X_foc, y = Y_foc, xend = X_reg, yend = Y_reg),
               color = "gray70", size = 0.5) +
  
  # Plot regular females
  geom_point(data = females %>% 
               filter(year == 2018, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "blue", size = 1.5) +
  
  # Plot focal females
  geom_point(data = focal_females %>% 
               filter(year == 2018, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "red", size = 3) +
  
  #Plot other males
  geom_point(data = males %>% 
               filter(year == 2018, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "yellow", size = 1.5) +
  
  labs(x = "X (meters)", y = "Y (meters)",
       title = "Focal Females (red) and Associated Regular Females (blue)") +
  theme_minimal()

ggsave("2018_focal_fem_SP.png",
       plot = focal_fem_SP_2018,
       width = 8,
       height = 6,
       dpi = 600)

#-------------------------------------------------------------------------------

#now calculating closest focal males
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

links <- females2 %>%
  left_join(
    focal_males %>%
      st_drop_geometry() %>%
      select(index, X, Y),
    by = c("closest_focal_male_index" = "index"),
    suffix = c("_reg", "_foc")
  ) %>%
  filter(!is.na(closest_focal_male_index))

#graphing 2016 
focal_male_SP_2016 <- ggplot() +
  # Draw lines connecting each regular female to its focal female
  geom_segment(data = links %>% filter(year == 2016, region == "south"),
               aes(x = X_foc, y = Y_foc, xend = X_reg, yend = Y_reg, color = "Association"),
               size = 0.5) +
  
  # Plot regular females
  geom_point(data = females2 %>% 
               filter(year == 2016, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "Female", size = 1.5) +
  
  # Plot focal males
  geom_point(data = focal_males %>% 
               filter(year == 2016, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "Focal Male", size = 3) +
  
  #Plot other males
  geom_point(data = males %>% 
               filter(year == 2016, region == "south") %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "Other Male", size = 1.5) +
  
  scale_color_manual(
    name = "Legend",
    values = c(
      "Association" = "gray70",
      "Female" = "blue",
      "Focal Male" = "red",
      "Other Male" = "yellow"
    )
  ) +
  
  labs(x = "X (meters)", y = "Y (meters)",
       title = "Focal Males (red) and Associated Regular Females (blue)") +
  theme_minimal()

ggsave("2016_focal_male_SP.png",
       plot = focal_male_SP_2016,
       width = 8,
       height = 6,
       dpi = 600)

#graphing 2018 
ggplot() +
  # Draw lines connecting each regular female to its focal female
  geom_segment(data = links %>% filter(year == 2018),
               aes(x = X_foc, y = Y_foc, xend = X_reg, yend = Y_reg),
               color = "gray70", size = 0.5) +
  
  # Plot regular females
  geom_point(data = females %>% 
               filter(year == 2018) %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "blue", size = 1.5) +
  
  # Plot focal males
  geom_point(data = focal_males %>% 
               filter(year == 2018) %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "red", size = 3) +
  
  #Plot other males
  geom_point(data = males %>% 
               filter(year == 2018) %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "yellow", size = 1.5) +
  
  labs(x = "X (meters)", y = "Y (meters)",
       title = "Focal Males (red) and Associated Regular Females (blue)") +
  theme_minimal()

#graphing 2019 
ggplot() +
  # Draw lines connecting each regular female to its focal female
  geom_segment(data = links %>% filter(year == 2019),
               aes(x = X_foc, y = Y_foc, xend = X_reg, yend = Y_reg),
               color = "gray70", size = 0.5) +
  
  # Plot regular females
  geom_point(data = females %>% 
               filter(year == 2019) %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "blue", size = 1.5) +
  
  # Plot focal males
  geom_point(data = focal_males %>% 
               filter(year == 2019) %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "red", size = 3) +
  
  #Plot other males
  geom_point(data = males %>% 
               filter(year == 2019) %>% 
               st_drop_geometry(),
             aes(x = X, y = Y),
             color = "yellow", size = 1.5) +
  
  labs(x = "X (meters)", y = "Y (meters)",
       title = "Focal Males (red) and Associated Regular Females (blue)") +
  theme_minimal()

library(units)

ggplot(links, aes(x = length, y = distance_focal_male, color = year)) +
  geom_smooth(se = FALSE, method = "lm") +  # use lm or loess (for smooth line)
  labs(
    x = "female length",
    y = "distance from focal male",
    color = "Year",
    title = "Relationship between Female Length and Harem Position"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right"
  )

ggsave(
  "TablesFigures/Female Length vs. Harem Position.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
