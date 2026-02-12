library(dplyr)
library(geosphere)
library(ggplot2)
library(patchwork)
library(tidyr)
library(readr)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

data <- read_csv("IntermediateData/uas.dataset.coords.csv")

data_clean <- data %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

uasdata <- data_clean %>%
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup") 
#I am now trying something new to make the male size vs number of females graph
# create a male index
uasdata1 <- uasdata %>%
  mutate(male_index = if_else(age_sex == "male",
                              row_number(),
                              NA_integer_))

# split sexes for ease
females <- uasdata1 %>%
  filter(age_sex == "female") 

males <- uasdata1 %>%
  filter(age_sex == "male")

# function to find closest male index (only within the same year)
find_closest_male <- function(lat, lon, female_year, males) {
  # filter males to the same year
  males_year <- males %>% filter(year == female_year)
  
  # make matrix of their coordinates
  male_coords <- as.matrix(males_year[, c("lat", "lon")])
  male_index  <- males_year$male_index
  
  # compute distances
  dists <- sqrt((male_coords[,1] - lat)^2 + (male_coords[,2] - lon)^2)
  
  # return the male index of the closest male in the same year
  male_index[which.min(dists)]
}

# apply this to every female
females$closest_male_index <- mapply(
  find_closest_male,
  females$lat,
  females$lon,
  females$year,
  MoreArgs = list(males = males)
)

# combine males and females back together
uasdata_full <- bind_rows(
  males %>% mutate(closest_male_index = NA_integer_),
  females
)

# calculate the number of females associated with each male
female_counts <- uasdata_full %>%
  filter(age_sex == "female") %>%
  count(closest_male_index, name = "n_females")

# join back to full dataset
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

#time to check if there were any males matched with females from different years (aka 2018 male matched with 2016 female)
females_check <- uasdata_with_counts %>%
  filter(age_sex == "female") %>%
  left_join(
    uasdata_with_counts %>%
      select(male_index, year) %>%
      rename(male_year = year),
    by = c("closest_male_index" = "male_index")
  )

mismatched_years <- females_check %>%
  filter(year != male_year)
#no mismatched years... yay!

#now you can do everything else

#-------------------------------------------------------------------------------------------------------
#here is where I make the binned graphs for male length vs harem size:

#make bins for the male lengths (0.25 intervals)
male_data_0.25 <- male_data %>%
  mutate(length_bin = cut(length,
                          breaks = seq(floor(min(length)), ceiling(max(length)), by = 0.25),
                          include.lowest = TRUE,
                          right = FALSE))

#summarize the total number of females for each bin
bin_summary <- male_data_0.25 %>%
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

#now by year 2016 binned by 0.1
male_data_2016 <- uasdata_with_counts %>% filter(age_sex == "male",
                                                 year == "2016")

male_data_2016_binned<- male_data_2016 %>%
  mutate(length_bin = cut(length,
                          breaks = seq(floor(min(length)), ceiling(max(length)), by = 0.10),
                          include.lowest = TRUE,
                          right = FALSE))

bin_summary_2016 <- male_data_2016_binned %>%
  group_by(year, length_bin) %>%
  summarize(total_females = sum(n_females), .groups = "drop")

ggplot(bin_summary_2016,
       aes(x = factor(length_bin), y = total_females)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Male length (m)",
    y = "Number of associated females") +
  theme_minimal() +
  #facet_wrap(~ year, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

male_data_2016 %>% 
  filter(age_sex == "male",
         length >= 3.6,
         length < 3.7) %>% 
  count()
#so there are 11 males with the length 3.6-3.7 (the range with the highest female association)

#then this shows for each of those 11 males exactly how many females are associated with them.
male_data_2016 %>%
  filter(
    age_sex == "male",
    length >= 3.6,
    length < 3.7
  ) %>%
  select(n_females)

male_data_2016 %>%
  filter(age_sex == "male") %>%
  arrange(desc(n_females)) %>%
  select(length, n_females)

#now doing all the years with a bin length of 0.1
male_data_0.1 <- male_data %>%
  mutate(length_bin = cut(length,
                          breaks = seq(floor(min(length)), ceiling(max(length)), by = 0.1),
                          include.lowest = TRUE,
                          right = FALSE))
#summarize the total number of females for each bin
bin_summary_0.1 <- male_data_0.1 %>%
  group_by(year, length_bin) %>%
  summarize(total_females = sum(n_females), .groups = "drop")


ggplot(bin_summary_0.1,
       aes(x = factor(length_bin), y = total_females)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Male length (m)",
    y = "Number of associated females") +
  theme_minimal() +
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "male length vs harem size - all years binned at 0.1.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

#now for all years binned at 0.05
male_data_0.05 <- male_data %>%
  mutate(length_bin = cut(length,
                          breaks = seq(floor(min(length)), ceiling(max(length)), by = 0.05),
                          include.lowest = TRUE,
                          right = FALSE))
#summarize the total number of females for each bin
bin_summary_0.05 <- male_data_0.05 %>%
  group_by(year, length_bin) %>%
  summarize(total_females = sum(n_females), .groups = "drop")


ggplot(bin_summary_0.05,
       aes(x = factor(length_bin), y = total_females)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Male length (m)",
    y = "Number of associated females") +
  theme_minimal() +
  facet_wrap(~ year, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "male length vs harem size - all years binned at 0.1.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

#----------------------------------------------------------------------------------------------
#ok so now I will try to assign the males with the most associated females (for each individual year) as focal males

#something is wrong here - come back to this and check it
uasdata_with_counts_2016 <- uasdata_with_counts %>% 
  filter(year == "2016")

focal_males <- c(444, 346, 606, 725, 351, 518, 719, 288, 798, 950, 832, 164, 686, 863, 860)

uasdata_full_2016<- uasdata_with_counts_2016 %>%
  mutate(
    sex_type = case_when(
      age_sex == "female" ~ "female",
      male_index %in% focal_males ~ "focal male",
      TRUE ~ "male"
    )
  )


library(sf)

females2016 <- uasdata_full_2016 %>% 
  filter(sex_type == "female")
focal_males2016 <- uasdata_full_2016 %>% 
  filter(sex_type == "focal male")

# make sf points
females_sf <- st_as_sf(females2016, coords = c("lon","lat"), crs = "EPSG:26943")
males_sf   <- st_as_sf(focal_males2016, coords = c("lon","lat"), crs = "EPSG:26943")

# compute distance matrix (meters)
dist_matrix <- st_distance(females_sf, males_sf)


















