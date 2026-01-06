setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

library(dplyr)
library(ggplot2)

#this line will load in your csv to your environment
dataset <- read_csv("IntermediateData/uasdata.csv")

#there is a typo in row 10543 where it says "malfee" instead of male, this will fix that
dataset1 <- dataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

#fixed it

#now we need to filter out the seals that are categorized as (N & water) not complete in the p_complete column
uasdata <- dataset1 %>%
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup")

#ok so this adds the column that differentiates between females, males, and alpha males
uasdata_main <- uasdata %>%
  group_by(year) %>%#this line communicates that 90th percentiles are calculated for each individual year
  mutate(
    alpha_cutoff = quantile(length[age_sex == "male"], 0.9, na.rm = TRUE),
    dominance_group = case_when(
      age_sex == "female" ~ "female",
      age_sex == "male" & length >= alpha_cutoff ~ "alpha male",
      age_sex == "male" ~ "male",
      age_sex == "pup" ~ "pup",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(-alpha_cutoff)

#this filters out only 2016 seals
uasdata_2016 <- uasdata_main %>% filter(year == 2016)

#ok so now we're getting into the spatial association method

library(geosphere)

females2016 <- uasdata_2016 %>% 
  filter(dominance_group == "female")
alpha_males2016 <- uasdata_2016 %>% 
  filter(dominance_group == "alpha male")

# create matrices of coordinates (lon, lat format)
females2016_coords <- females2016[, c("lon", "lat")]
alpha_males2016_coords <- alpha_males2016[, c("lon", "lat")]

# calculate pairwise distances (in meters by default)
dist_matrix <- distm(females2016_coords, alpha_males2016_coords, fun = distHaversine)

# get the min distance for each female from highlight_male
min_distances2016 <- apply(dist_matrix, 1, min)

# get the index of the closest male as well
closest_indices2016 <- apply(dist_matrix, 1, which.min)

# add results to the females data frame as new columns so we can graph them later
females2016$min_distance_to_alpha_male_m <- min_distances2016
females2016$closest_alpha_male_index_2016 <- closest_indices2016

alpha_males2016$index=1:nrow(alpha_males2016)

# Join to get coordinates of closest male for each female
females_with_males = females2016 %>%
  left_join(alpha_males2016 %>% 
              select(index, male_lat = lat, male_lon = lon),
            by = c("closest_alpha_male_index_2016" = "index"))


my_40_colors <- c(
  "1" = "#E41A1C",  "2" = "#377EB8",  "3" = "#4DAF4A",  "4" = "#984EA3",  "5" = "#FF7F00",
  "6" = "#FFFF33",  "7" = "#A65628",  "8" = "#F781BF",  "9" = "#999999", "10" = "#66C2A5",
  "11" = "#FC8D62", "12" = "#8DA0CB", "13" = "#E78AC3", "14" = "#A6D854", "15" = "#FFD92F",
  "16" = "#E5C494", "17" = "#B3B3B3", "18" = "#1B9E77", "19" = "#D95F02", "20" = "#7570B3",
  "21" = "#66A61E", "22" = "#E7298A", "23" = "#A6761D", "24" = "#666666", "25" = "#1F78B4",
  "26" = "#33A02C", "27" = "#FB9A99", "28" = "#E31A1C", "29" = "#FDBF6F", "30" = "#FF7F00",
  "31" = "#CAB2D6", "32" = "#6A3D9A", "33" = "#B2DF8A", "34" = "#33A02C", "35" = "#FFFF99",
  "36" = "#B15928", "37" = "#FBB4AE", "38" = "#B3CDE3", "39" = "#CCEBC5", "40" = "#DECBE4"
)

highlight_levels <- sort(unique(c(
  females2016$closest_alpha_male_index_2016,
  alpha_males2016$index
)))

#here's the plot
ggplot(data = females2016) +
  geom_segment(data = females_with_males,
               aes(x = lon, y = lat, xend = male_lon, yend = male_lat,
                   group = closest_alpha_male_index_2016),
               color = "gray50", alpha = 0.5) +
  geom_point(aes(x = lon, y = lat,
                 fill = factor(closest_alpha_male_index_2016, levels = highlight_levels)),   # Females with color for their closest male
             shape = 21, color = "black", alpha = 0.7) +
  geom_point(data = alpha_males2016,   # Highlight males with same fill color scale
             aes(x = lon, y = lat, fill = factor(index, levels = highlight_levels)), 
             shape = 22, size = 3, color = "black") +
  geom_text(data = alpha_males2016,
            aes(x = lon, y = lat, label = index), size = 2) +
  scale_fill_manual(values = my_40_colors, name = "Alpha Male Index") +
  scale_y_log10() +
  theme_minimal() +
  labs(x = "Longitude",
       y = "Latitude") +
  guides(fill = "none")
  coord_cartesian(xlim = c(-122.330, -122.3275),  # set longitude limits
                  ylim = c(37.1125, 37.1150))    
  
  
  ggsave(
    "spatialassociation_2016_firstdraft.png",
    width = 8,
    height = 6,
    units = "in",
    dpi = 600
  )

