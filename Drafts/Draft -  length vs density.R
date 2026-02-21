#how does density change with female length?

library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(sf)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

uasdataset <- st_read("IntermediateData/uasdata.full.shp")

uasdata <- uasdataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

females <- uasdata %>% 
  filter(age_sex == "female")

ggplot(females, aes(x = length, y = density, color = year)) +
  geom_smooth(se = FALSE, method = "lm") +  # use lm or loess (for smooth line)
  labs(
    x = "female length",
    y = "surrounding density",
    color = "Year",
    title = "Relationship between Female Length and Surrounding Density"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right"
  )

males <- uasdata %>% 
  filter(age_sex == "male")

ggplot(males, aes(x = length, y = density, color = year)) +
  geom_smooth(se = FALSE, method = "loess") +  # use lm or loess (for smooth line)
  labs(
    x = "male length",
    y = "surrounding density",
    color = "Year",
    title = "Relationship between Male Length and Surrounding Density"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right"
  )


