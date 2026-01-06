#Here are all the necessary packages to create the Linear Graph

install.packages("tidyverse")
library(tidyverse)
install.packages(readr)
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



#this is the original method of categorizing alphas based on length.
#uasdata_complete$highlight_group <- with(
  uasdata_complete,
  ifelse(age_sex == "male" & length>=4, "alpha_male", age_sex)
)#when i ran this, I realized there just weren't males equal to or above 4m in the year 2016 and only 1 in 2018 so I will now try to come up with a way of getting the alpha males to be the seals within the top 10% of length for each individual year.


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
      age_sex == "pup" ~ "pup",
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


ggplot(data = females, mapping = aes(x = area,
                                     y = min_distance_to_alpha_male_m,
                                     color = as.factor(year))) +
  geom_point(alpha = 0.0) +
  scale_x_log10() +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.1) + 
  labs(color = "Year",
       x = "Female Body Size (area)",
       y = "Distance from Alpha Male (m)") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 20))

ggsave(
  "linearplot_firstdraft.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

#soooo.... this looks significantly different from our original results, so I will run stats models on this to see if the difference is actually significant.

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#here is where the stats analysis starts

#this gives you the number of females associated with each alpha male across all years
female_counts_full <- females %>%   
  count(closest_alpha_male_index, name = "num_females")
#running this gives a list of all 154 alphas and their females, but a lot of the alphas have under 10, even just 1 female associated, so maybe we filter those out of the alphas cuz are they alphas if they so little females?

average_females <- female_counts_full %>%
  summarise(avg_females_per_male = mean(num_females)) %>%
  pull(avg_females_per_male)

print(average_females)
#so over all years there is an average of 79.22973 females per alpha male

#min & max of female body size
min(females$area, na.rm = TRUE)
max(females$area, na.rm = TRUE)

#Linear Model

#splits data by year
females_by_year <- split(females, females$year)

#splits by the linear models by year
models_by_year <- lapply(females_by_year, function(df) {
  lm(min_distance_to_alpha_male_m ~ log10(area), data = df)
})

#calculates r squared by year
r_squared_by_year <- sapply(models_by_year, function(model) summary(model)$r.squared)

#filters the r squared values out in a clean data frame
r2_results <- data.frame(
  year = names(r_squared_by_year),
  r_squared = r_squared_by_year
)
print(r2_results)#this gave us weak relationship strength numbers

#calculates slopea and p-values for the linear model
slopes <- sapply(models_by_year, function(model) coef(summary(model))[2, 1])
p_values <- sapply(models_by_year, function(model) coef(summary(model))[2, 4])

print(format(p_values, scientific = FALSE))#ok so this showed that this relationship is signficant for some years and not for others......


#makes a whole data frame with the stats 
regression_stats <- data.frame(
  year = names(models_by_year),
  slope = slopes,
  r_squared = r_squared_by_year,
  p_value = p_values
)
print(regression_stats)


#Here are all the linear models tested.
all_years = lm(min_distance_to_alpha_male_m ~ area, data = females) # Assumes no difference across years
year_effect = lm(min_distance_to_alpha_male_m ~ area + year, data = females) # Assumes year is main effect
full_model = lm(min_distance_to_alpha_male_m ~ area * year, data = females) # Full model where slopes and intercepts vary across year

summary(all_years)
# Use summary(full_model) to look for the interaction term (body_size:year)
# Use the anova() to test if the slopes are significantly different across years

summary(year_effect)
summary(full_model)
anova(year_effect, full_model) # Put in the model that assumes year is the greatest effect and the full model that tests for the interaction between body_size and year

# Look for the p-value (Pr(>F)), if p<0.05, it means that full_model is significantly  better than year_effect which means your linear models are significantly different across the 10 years of the study (i.e., the slopes and intercepts are not consistent)

#If you want to look at a specific year:
all_years = lm(min_distance_to_highlight_male_m ~ area_, data = females) # Same as before, all data, all together without year 
specific_year = lm(distance ~ body_size, data = uasdata[uasdata$year == specific_year], # This is the model for a single year
                   anova(all_years, specific_year)) # Tells if a specific year is significantly different than the model that describes all the years



