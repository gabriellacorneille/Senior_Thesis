#this will be looking at the relationship between female length and distance from focal males

library(tidyverse)
library(readr)
library(dplyr)
library(geosphere)
library(ggplot2)
library(patchwork)

#now set your working directory for this specific project
setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

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

#here is my attempt at categorizing alphas as the top ten percent length wise based on each individual year. 
uasdata_top10 <- uasdata %>%
  group_by(year) %>%
  filter(
    age_sex == "male",
    length >= quantile(length[age_sex=="male"], 0.9, na.rm = TRUE)) %>%
  ungroup()
#now this will list out the number of males within the top 10% of length for each individual year (these will be your alphas)
uasdata_top10 %>% 
  count(year)

#now is there a way to reintegrate this info into the main dataset?

#this for top 10 percent and gives me a range of (12-27) males per year
uasdata_main <- uasdata %>%
  group_by(year) %>%#this line communicates that 90th percentiles are calculated for each individual year
  mutate(
    alpha_cutoff = quantile(length[age_sex == "male"], 0.9, na.rm = TRUE),
    dominance_group = case_when(
      age_sex == "female" ~ "female",
      age_sex == "male" & length >= alpha_cutoff ~ "alpha male",
      age_sex == "male" ~ "male",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(-alpha_cutoff)

#this tells me how many alphas there are each year (just to check that its giving me the right amount of alphas per year) 
uasdata_main%>%
  filter(dominance_group == "alpha male") %>%
  count(year)

#this tells you how many total across all years
uasdata_main%>%
  filter(dominance_group == "alpha male") %>%
  count()

#now we're only working with uasdata_main2

#here is the start of creating the linear graph

#split  data up by class
females <- uasdata_main %>%
  filter(dominance_group == "female") 

alpha_males <- uasdata_main %>%
  filter(dominance_group == "alpha male")

# create matrices of coordinates (lon, lat format?)
females_coords <- females[, c("lon", "lat")]
alpha_males_coords <- alpha_males[, c("lon", "lat")]

# calculate pairwise distances (in meters by default)
dist_matrix <- distm(females_coords, alpha_males_coords, fun = distHaversine)

# get the min distance for each female from highlight_male
min_distances <- apply(dist_matrix, 1, min)

# get the index of the closest male as well
closest_indices <- apply(dist_matrix, 1, which.min)

# add results to the females data frame as new columns so we can graph them later
females$min_distance_to_alpha_male_m <- min_distances
females$closest_alpha_male_index <- closest_indices


ggplot(data = females, mapping = aes(x = length,
                                     y = min_distance_to_alpha_male_m,
                                     color = as.factor(year))) +
  geom_point(alpha = 0.0) +
  scale_x_log10() +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1) + 
  labs(color = "Year",
       x = "Female Length (m)",
       y = "Distance from Alpha Male (m)") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 20))

ggsave(
  "linearplot_female_length.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
