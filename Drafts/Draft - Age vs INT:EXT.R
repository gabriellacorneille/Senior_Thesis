library(readr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)

setwd("/Users/gabbycorneille/Desktop/Senior Thesis/Senior_Thesis")

resights <- read_csv("IntermediateData/Resight Data.csv")

resights$position <- ifelse(resights$"Harem Position" == "INT", 1,
                        ifelse(resights$"Harem Position" == "EXT", 2, NA))

ggplot(resights, aes(x = `Age (actual)`, y = factor(position))) +
  geom_jitter(height = 0.1, width = 0) +
  labs(x = "Age",
       y = "Location Type",
       title = "INT vs EXT by Age")

resights$position <- ifelse(resights$`Harem Position` == "INT", 0, 1)

#-------------------------------------------------------------------------------
logit_model <- glm(position ~ `Age (actual)`, 
                   data = resights, 
                   family = binomial)
summary(logit_model)

anova(logit_model, test = "Chisq")
#-------------------------------------------------------------------------------      
p_age <- signif(summary(logit_model)$coefficients[2,4], 3)  # p-value
n_obs <- nobs(logit_model)

stat_label <- paste0("p-value = ", p_age,
                     "\nn = ", n_obs)

ggplot(resights, aes(x = `Age (actual)`, y = position)) +
  geom_jitter(height = 0.05, alpha = 0.4) +
  geom_smooth(method = "glm",
              method.args = list(family = "binomial"),
              se = TRUE) +
  scale_y_continuous(
    breaks = c(0, 1),
    labels = c("Center", "Edge")
  ) +
  labs(x = "Age",
       y = "Harem Position",
       title = "Logistic Regression: Age vs Harem Position")

ggsave(
  "TablesFigures/resight data/Resight Age vs Harem Position.png",
  width = 8,
  height = 6,
  units = "in",
  dpi = 600
)






