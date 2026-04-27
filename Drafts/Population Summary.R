#calculating average population on the beach

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

#males
uasdata %>%
  sf::st_drop_geometry() %>%
  filter(age_sex == "male") %>%
  count(year)

uasdata %>%
  sf::st_drop_geometry() %>%
  filter(age_sex == "male") %>%
  count(year) %>% 
  summarise(mean_count = mean(n))

#females
uasdata %>% 
  sf::st_drop_geometry() %>%
  filter(age_sex =="female") %>% 
  count(year)

uasdata %>% 
  sf::st_drop_geometry() %>%
  filter(age_sex =="female") %>% 
  count(year) %>% 
  summarise(mean_count = mean(n))
