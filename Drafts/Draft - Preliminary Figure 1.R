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
data <- read_csv("IntermediateData/uasdata.csv")

#there is a typo in row 10543 where it says "malfee" instead of male, this will fix that
uasdata <- data %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

#fixed it

#now we need to filter out the seals that are categorized as (N & water) not complete in the p_complete column
uasdata_complete <- uasdata %>%
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup")



#this is the original method of categorizing alphas based on length.
#uasdata_complete$highlight_group <- with(
  uasdata_complete,
  ifelse(age_sex == "male" & length>=4, "alpha_male", age_sex)
)#when i ran this, I realized there just weren't males equal to or above 4m in the year 2016 and only 1 in 2018 so I will now try to come up with a way of getting the alpha males to be the seals within the top 10% of length for each individual year.


#here is my attempt at categorizing alphas as the top ten percent length wise based on each individual year. 
uasdata_complete_top10 <- uasdata_complete %>%
  group_by(year) %>%
  filter(
    age_sex == "male",
    length >= quantile(length[age_sex=="male"], 0.9, na.rm = TRUE)) %>%
  ungroup()
#now this will list out the number of males within the top 10% of length for each individual year (these will be your alphas)
uasdata_complete_top10 %>% 
  count(year)

#now is there a way to reintegrate this info into the main dataset?

#this for top 10 percent and gives me a range of (12-27) males per year
uasdata_main2 <- uasdata_complete %>%
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
uasdata_main2%>%
  filter(dominance_group == "alpha male") %>%
  count(year)


#now we're only working with uasdata_main2

#here is the start of creating the linear graph

#split  data up by class
females <- uasdata_main2 %>%
  filter(dominance_group == "female") 

alpha_males <- uasdata_main2 %>%
  filter(dominance_group == "alpha male")

#ok now we have run into a problem, this new dataset does not have location data for each seal like the dataset from 128L...
#but luckily once we figure that out I can use the code from 128L project to code the easily. 










