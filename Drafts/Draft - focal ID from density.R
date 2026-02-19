#here I will try to identify focals through density assessment

library(readr)
library(sf)
library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

uasdataset <- st_read("IntermediateData/uasdata.full.shp")

uasdata <- uasdataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

#---------------------THERE IS A PROBLEM HERE! (the problem is that i'm confused what the columns female99 and male90 mean)
#this gives a table of how many females and males we have per year.
table(uasdata$year, uasdata$age_sex)

counts <- table(uasdata$year, uasdata$age_sex)

colMeans(counts)

uasdata_full <- uasdata %>%
  group_by(year, age_sex) %>%
  mutate(
    female_99 = quantile(density[age_sex == "female"], 0.99, na.rm = TRUE),
    male_90   = quantile(density[age_sex == "male"],   0.90, na.rm = TRUE),
    status = case_when(
      age_sex == "female" & density >= female_99 ~ "focal female",
      age_sex == "male"   & density >= male_90   ~ "focal male",
      age_sex == "female" ~ "female",
      age_sex == "male"   ~ "male"
    )
  ) %>%
  ungroup()

table(uasdata_full$year, uasdata_full$status)

