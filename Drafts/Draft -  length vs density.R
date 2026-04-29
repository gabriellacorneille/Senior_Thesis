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

#-------------------------------------------------------------------------------
#graphing female length vs surrounding density
females <- uasdata %>% 
  filter(age_sex == "female")

#-------------
#run stats

# fit linear model
fem_model1 <- lm(density ~ length, data = females)

fem_model2 <- lm(density ~ length + year, data = females)

fem_r2 <- summary(fem_model1)$r.squared
fem_pval <- summary(fem_model1)$coefficients[2,4]
# View stats results
summary(fem_model1)
summary(fem_model2)
summary(fem_model2)$coefficients
summary(fem_model)$r.squared
#-------------

ggplot(females, aes(x = length, y = density, color = year)) +
  geom_smooth(se = FALSE, method = "lm") +  # use lm or loess (for smooth line)
  labs(
    x = "Female Length (m)",
    y = "Surrounding Density (# of seals)",
    color = "Year",
    title = "Relationship between Female Length and Surrounding Density"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right")
 
ggsave(
  "TablesFigures/female length vs harem position/Female Length vs. Surrounding Density.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

#-------------------------------------------------------------------------------
#graphing male length vs surrounding density
males <- uasdata %>% 
  filter(age_sex == "male")

#-------------
#run stats

# fit linear model
male_model1 <- lm(density ~ length, data = males)

amale_model2 <- lm(density ~ length + year, data = males)

male_r2 <- summary(male_model1)$r.squared
male_pval <- summary(male_model1)$coefficients[2,4]
# View stats results
summary(male_model1)
summary(male_model2)
summary(male_model2)$coefficients
summary(male_model2)$r.squared
#-------------

ggplot(males, aes(x = length, y = density, color = year)) +
  geom_smooth(se = FALSE, method = "lm") +  # use lm or loess (for smooth line)
  labs(
    x = "Male Length (m)",
    y = "Surrounding Density (# of seals)",
    color = "Year",
    title = "Relationship between Male Length and Surrounding Density"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right") 

ggsave(
  "TablesFigures/bigger male bigger harem?/Male Length vs. Surrounding Density.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)
#-------------------------------------------------------------------------------
