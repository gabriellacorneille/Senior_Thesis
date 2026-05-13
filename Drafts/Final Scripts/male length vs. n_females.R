#Male length vs. Harem size (number of females)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(readr)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

data <- read_csv("IntermediateData/Final Data Files/uas.dataset.coords.csv")

#------------------------------------------------------------------------------
#initial data clean up
data_clean <- data %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

uasdata <- data_clean %>%
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup") 

#------------------------------------------------------------------------------
#organizing data and creating n_females

# create a male index
uasdata1 <- uasdata %>%
  mutate(male_index = if_else(age_sex == "male",
                              row_number(),
                              NA_integer_))
# split sexes
females <- uasdata1 %>%
  filter(age_sex == "female") 

males <- uasdata1 %>%
  filter(age_sex == "male")

#function to find closest male index (only within the same year)
find_closest_male <- function(lat, lon, female_year, males) {
  #filter males to the same year
  males_year <- males %>% filter(year == female_year)
  
  #make matrix of coordinates
  male_coords <- as.matrix(males_year[, c("lat", "lon")])
  male_index  <- males_year$male_index
  
  #calculate distances
  dists <- sqrt((male_coords[,1] - lat)^2 + (male_coords[,2] - lon)^2)
  
  #return the male index of the closest male in the same year
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

#combine males and females back together
uasdata_full <- bind_rows(
  males %>% mutate(closest_male_index = NA_integer_),
  females
)

#calculate the number of females associated with each male
female_counts <- uasdata_full %>%
  filter(age_sex == "female") %>%
  count(closest_male_index, name = "n_females")

#join back to full dataset
uasdata_with_counts <- uasdata_full %>%
  left_join(
    female_counts,
    by = c("male_index" = "closest_male_index")
  ) %>%
  mutate(
    #fill NA with 0 for males, keep NA for females
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
#mismatched_years should have 0 obs.

#-----------------------
#linear model code

male_data <- uasdata_with_counts %>% filter(age_sex == "male")

male_data$year <- as.factor(male_data$year)

# Fit linear model
model1 <- lm(n_females ~ length, data = male_data)

model2 <- lm(n_females ~ length + year, data = male_data)

#change the "model" based on whether you want to include year or not (model1 or model2)
r2 <- summary(model1)$r.squared
pval <- summary(model1)$coefficients[2,4]

# View stats results
summary(model1)
summary(model2)
#-----------------------

#plot
ggplot(male_data, aes(x = length, y = n_females, color = year)) +
  geom_smooth(se = FALSE, method = "lm") +  # use lm or loess (for smooth line)
  labs(
    x = "Male Length (m)",
    y = "Number of Females",
    color = "Year",
    #title = "Relationship between Male Length and Number of Females by Year"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right")

#save
ggsave(
  "TablesFigures/bigger male bigger harem?/male_len vs n_fem (no title).png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

#------------------------------------------------------------------------------
#end.