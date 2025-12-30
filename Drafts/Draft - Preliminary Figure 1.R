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
uasdata <- read_csv("IntermediateData/uasdata.csv")

#we need to filter out the seals that are categorized as (N & water) not complete in the p_complete column
uasdata_complete <- uasdata %>%
  filter(is.na(p_complete) | !p_complete %in% c("N", "water"))

#here I am trying to categorize males above a certain length as alpha males
uasdata_complete$highlight_group <- with(
  uasdata_complete,
  ifelse(age_sex == "male" & length>=4, "alpha_male", age_sex)
)#when i ran this, I realized there just weren't males equal to or above 4m in the year 2016 and only 1 in 2018 so I will now try to come up with a way of getting the alpha males to be the seals within the top 10% of length for each individual year.

#here is my attempt 
uasdata_complete_top10 <- uasdata_complete %>%
  group_by(year) %>%
  filter(
    age_sex == "male",
    length >= quantile(length, 0.9, na.rm = TRUE)) %>%
  ungroup()
#now this will list out the number of males within the top 10% of length for each individual year (these will be your alphas)
uasdata_complete_top10 %>% 
  count(year)
#this gives an range of 108 to 247 alpha males across the years which is probably too many so let's make it the top 5 percent of lengths.
uasdata_complete_top5 <- uasdata_complete %>%
  group_by(year) %>%
  filter(
    age_sex == "male",
    length >= quantile(length, 0.95, na.rm = TRUE)) %>%
  ungroup()

uasdata_complete_top5 %>% 
  count(year)
#using the top 5 percent gives us a range of 104 to 154 alpha males across individual years

#now is there a way to reintegrate this info into the main dataset?




