install.packages("tidyverse")
library(tidyverse)
install.packages(readr)
library(readr)
library(dplyr)
library(geosphere)
library(ggplot2)
library(patchwork)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

data <- read_csv("IntermediateData/uasdata.csv")

data_clean <- data %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

uasdata <- data_clean %>%
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup") 

#split the data up by class into their own datasets
females <- uasdata %>%
  filter(age_sex == "female") 

males <- uasdata %>%
  filter(age_sex == "male")

#now make matrices of coordinates
females_coords <- females[, c("lon", "lat")]
males_coords <- males[, c("lon", "lat")]

# calculate pairwise distances (in meters by default)
dist_matrix <- distm(females_coords, males_coords, fun = distHaversine)

# get the min distance for each female from highlight_male
min_distances <- apply(dist_matrix, 1, min)

# get the index of the closest male as well
closest_indices <- apply(dist_matrix, 1, which.min)

# add results to the females data frame as new columns so we can graph them later
females$min_distance_to_male_m <- min_distances
females$closest_male_index <- closest_indices

#plot it up - this is the relationship between female body size and distance from male
ggplot(data = females, mapping = aes(x = area,
                                     y = min_distance_to_male_m,
                                     color = as.factor(year))) +
  geom_point(alpha = 0.0) +
  scale_x_log10() +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1) + 
  labs(color = "Year",
       x = "Female Body Size (area)",
       y = "Distance from Male (m)") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 20))

ggsave(
  "linearplot_all_males.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
#-------------------------------------------------------------------------------
#I am now trying something new to make the male szie vs number of females graph

#create a male index
uasdata1 <- uasdata %>%
  mutate(male_index = if_else(age_sex == "male",
                              row_number(),
                              NA_integer_))

#now split the sexes for ease 
females <- uasdata1 %>%
  filter(age_sex == "female") 

males <- uasdata1 %>%
  filter(age_sex == "male")

# matrix of male coordinates
male_coords <- as.matrix(males[, c("lat", "lon")])

# function to find closest male index
find_closest_male <- function(lat, lon, male_coords, male_index) {
  dists <- sqrt((male_coords[,1] - lat)^2 +
                  (male_coords[,2] - lon)^2)
  male_index[which.min(dists)]
}

#now apply this to every female
females$closest_male_index <- mapply(
  find_closest_male,
  females$lat,
  females$lon,
  MoreArgs = list(
    male_coords = male_coords,
    male_index  = males$male_index
  )
)

#now bring the females and males files back together
uasdata_full <- bind_rows(
  males %>% mutate(closest_male_index = NA_integer_),
  females
)

#now calculate the number of females each male is associated with (based on his index)
female_counts <- uasdata_full %>%
  filter(age_sex == "female") %>%
  count(closest_male_index, name = "n_females")

uasdata_with_counts <- uasdata_full %>%
  left_join(
    female_counts,
    by = c("male_index" = "closest_male_index")
  ) %>%
  mutate(
    # Fill NA with 0 for males, keep NA for females
    n_females = if_else(age_sex == "male", n_females, NA_integer_),
    n_females = replace_na(n_females, 0L)
  )


#now try to plot male body size vs. number of females associated
#this gives us a scatterplot (meh)
ggplot(uasdata_with_counts %>% filter(age_sex == "male"),
       aes(x = length, y = n_females)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    x = "Male length (m)",
    y = "Number of associated females") +
  theme_minimal()

male_data <- uasdata_with_counts %>% filter(age_sex == "male")

#this gives us a crazy bar graph (i think we can make it better)
ggplot(male_data,
       aes(x = factor(length), y = n_females)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Male length (m)",
    y = "Number of associated females") +
  theme_minimal()

ggsave(
  "bar graph (all males).png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
#here we go trying to make it better
#make bins for the male lengths (0.5 intervals)
male_data <- male_data %>%
  mutate(length_bin = cut(length,
                          breaks = seq(floor(min(length)), ceiling(max(length)), by = 0.5),
                          include.lowest = TRUE,
                          right = FALSE))
#summarize the total number of females for each bin
bin_summary <- male_data %>%
  group_by(length_bin) %>%
  summarize(total_females = sum(n_females), .groups = "drop")

#now graph
ggplot(bin_summary, aes(x = length_bin, y = total_females)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "Male length (intervals of 0.5)",
    y = "Number of associated females") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#this is kinda crazy i'm not even gonna lie...... 

ggsave(
  "bar graph (binned).png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
