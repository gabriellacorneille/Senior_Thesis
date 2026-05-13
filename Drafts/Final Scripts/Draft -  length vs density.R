#Length vs. (surrounding) density (females and males separately)

#set up
library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(sf)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

uasdataset <- st_read("IntermediateData/Final Data Files/uasdata.full.shp")

uasdata <- uasdataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#graphing female length vs surrounding density
females <- uasdata %>% 
  filter(age_sex == "female")

#-------------
# fit linear model
fem_model1 <- lm(density ~ length, data = females)
fem_model2 <- lm(density ~ length + year, data = females) #model includes effect of year

fem_r2 <- summary(fem_model1)$r.squared
fem_pval <- summary(fem_model1)$coefficients[2,4]

#summarize
summary(fem_model1)
summary(fem_model2)
#-------------

#plot
ggplot(females, aes(x = length, y = density, color = year)) +
  geom_smooth(se = FALSE, method = "lm") +  # use lm or loess (for smooth line)
  labs(
    x = "Female Length (m)",
    y = "Surrounding Density (# of seals)",
    color = "Year",
    #title = "Relationship between Female Length and Surrounding Density"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right")

#save
ggsave(
  "TablesFigures/female length vs harem position/fem_len vs density(no title).png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#graphing male length vs surrounding density
males <- uasdata %>% 
  filter(age_sex == "male")

#-------------
# fit linear model
male_model1 <- lm(density ~ length, data = males)
male_model2 <- lm(density ~ length + year, data = males) #model includes effect of year

male_r2 <- summary(male_model1)$r.squared
male_pval <- summary(male_model1)$coefficients[2,4]

# summarize
summary(male_model1)
summary(male_model2)
#-------------

#plot
ggplot(males, aes(x = length, y = density, color = year)) +
  geom_smooth(se = FALSE, method = "lm") +  # use lm or loess (for smooth line)
  labs(
    x = "Male Length (m)",
    y = "Surrounding Density (# of seals)",
    color = "Year",
    #title = "Relationship between Male Length and Surrounding Density"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "right") 

#save
ggsave(
  "TablesFigures/bigger male bigger harem?/male_len vs density(no title).png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)

#-------------------------------------------------------------------------------
#end.