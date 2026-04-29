#this is where I compare female length to harem position based on the focal males ID'd through density

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

#-------------------------------------------------------------------------------
library(units)

#this is checking the class of each variable (they are different)
class(links$distance_focal_male)
class(links$length)

#now we will try to make the focus variables the same class
links <- links %>%
  mutate(distance_focal_male = as.numeric(distance_focal_male))
# Fit linear model
model1 <- lm(distance_focal_male ~ length, data = links)


model2 <- lm(distance_focal_male ~ length + year, data = links)

#change the "model" based on whether you want to include year or not (model1 or model2)
r2 <- summary(model1)$r.squared
pval <- summary(model1)$coefficients[2,4]
# View stats results
summary(model2)
summary(model1)
summary(model)$coefficients
summary(model1)$r.squared

#here's the plot
ggplot(links, aes(x = length, y = distance_focal_male, color = year)) +
  geom_smooth(se = FALSE, method = "lm") +  # use lm or loess (for smooth line)
  labs(
    x = "Female Length (m)",
    y = "Distance from Focal Male (m)",
    color = "Year",
    title = "Relationship between Female Length and Harem Position"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right") 

#save
ggsave(
  "TablesFigures/female length vs harem position/Female Length vs. Harem Position.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)




