library(readr)
library(dplyr)
library(tidyr)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

resight_data <- read_csv("IntermediateData/Thesis Resights.csv")

library(ggplot2)

ggplot(resight_data, aes(x = "Harem Position", fill = "Age (actual)")) +
  geom_bar(position = "dodge") +
  labs(x = "Age",
       y = "Count",
       fill = "Age (actual)") +
  theme_minimal()
