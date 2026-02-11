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

male_data <- uasdata_with_counts %>% filter(age_sex == "male")


#make bins for the male lengths (0.5 intervals)
male_data <- male_data %>%
  mutate(length_bin = cut(length,
                          breaks = seq(floor(min(length)), ceiling(max(length)), by = 0.25),
                          include.lowest = TRUE,
                          right = FALSE))
#summarize the total number of females for each bin
bin_summary <- male_data %>%
  group_by(year, length_bin) %>%
  summarize(total_females = sum(n_females), .groups = "drop")


ggplot(bin_summary,
       aes(x = factor(length_bin), y = total_females)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Male length (m)",
    y = "Number of associated females") +
  theme_minimal() +
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "male length vs harem size - all years binned.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
