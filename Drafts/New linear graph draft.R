library(tidyverse)
library(readr)
library(dplyr)
library(geosphere)
library(ggplot2)
library(patchwork)
library(sf)

#now set your working directory for this specific project
setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

#this line will load in your csv to your environment
uasdata <- st_read("IntermediateData/uasdata.full.shp")

uasdata <- st_drop_geometry(uasdata)

#------------------------------------------------------------------------------
#This categorizes "alphas" as top ten percenters

uasdata_top10 <- uasdata %>%
  group_by(year) %>%
  filter(
    age_sex == "male",
    length >= quantile(length[age_sex=="male"], 0.9, na.rm = TRUE)) %>%
  ungroup()

#now this will list out the number of males within the top 10% of length for each individual year (these will be your alphas)
uasdata_top10 %>% 
  count(year)

#this integrates the top 10 percent thing into the main df 
uasdata_main <- uasdata %>%
  group_by(year) %>%#this line communicates that 90th percentiles are calculated for each individual year
  mutate(
    focal_cutoff = quantile(length[age_sex == "male"], 0.9, na.rm = TRUE),
    type_group = case_when(
      age_sex == "female" ~ "female",
      age_sex == "male" & length >= focal_cutoff ~ "focal male",
      age_sex == "male" ~ "male",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(-focal_cutoff)

#this tells me how many alphas there are each year (just to check that its giving me the right amount of alphas per year) 
uasdata_main%>%
  filter(type_group == "focal male") %>%
  count(year)

#------------------------------------------------------------------------------
#Now here is my method for associating females with the alphas spatially (taken from 3rd attempt - male length vs. harem (current).R)


# create a male index
uasdata1 <- uasdata_main %>%
  mutate(focal_male_index = if_else(age_sex == "focal male",
                              row_number(),
                              NA_integer_))

# split sexes for ease
females <- uasdata1 %>%
  filter(type_group == "female") 

focal_males <- uasdata1 %>%
  filter(type_group == "focal male")

# function to find closest male index (only within the same year)
find_closest_focal_male <- function(Y, X, female_year, focal_males) {
  # filter males to the same year
  focal_males_year <- focal_males %>% 
    filter(year == female_year)
  
  # make matrix of their coordinates
  focal_male_coords <- as.matrix(focal_males_year[, c("Y", "X")])
  focal_male_index  <- focal_males_year$focal_male_index
  
  # compute distances
  dists <- sqrt((focal_male_coords[,1] - Y)^2 + (focal_male_coords[,2] - X)^2)
  
  # return the male index of the closest male in the same year
  focal_male_index[which.min(dists)]
}

# apply this to every female
females <- females %>% 
  mutate(
    closest_focal_male_index = mapply(
      find_closest_focal_male,
      Y,
      X,
      year,
      MoreArgs = list(focal_males = focal_males)
    )
  )
#--------------------checking issue
nrow(focal_males)
table(focal_males$year)
head(focal_males)


# combine males and females back together
uasdata_full <- bind_rows(
  focal_males %>% mutate(closest_focal_male_index = NA_integer_),
  females
)

# calculate the number of females associated with each male
female_counts <- uasdata_full %>%
  filter(type_group == "female") %>%
  count(closest_focal_male_index, name = "n_females")

#------------------------------------------------------------------------------
#testing an issue right here
table(females$year)
table(focal_males$year)

class(females$year)
class(focal_males$year)

colnames(focal_males)

setdiff(unique(females$year), unique(focal_males$year))

sum(is.na(females$X))
sum(is.na(females$Y))

i <- 1

find_closest_focal_male(
  females$Y[i],
  females$X[i],
  females$year[i],
  focal_males)



# join back to full dataset
uasdata_with_counts <- uasdata_full %>%
  left_join(
    female_counts,
    by = c("focal_male_index" = "closest_focal_male_index")
  ) %>%
  mutate(
    # Fill NA with 0 for males, keep NA for females
    n_females = if_else(type_group == "focal male", n_females, NA_integer_),
    n_females = replace_na(n_females, 0L)
  )

#time to check if there were any males matched with females from different years (aka 2018 male matched with 2016 female)
females_check <- uasdata_with_counts %>%
  filter(type_group == "female") %>%
  left_join(
    uasdata_with_counts %>%
      select(focal_male_index, year) %>%
      rename(focal_male_year = year),
    by = c("closest_focal_male_index" = "focal_male_index")
  )

mismatched_years <- females_check %>%
  filter(year != focal_male_year)















