#here is how to make the density plots

plot(st_geometry(seal.buffer[12,]))

plot(st_geometry(survey.subset[12,]), add = TRUE, col = "red")

c <- int[[12]] 

nearest.neighbors <- survey.subset[c,]

plot(st_geometry(nearest.neighbors), add = TRUE)

#this filters out to 2016 specifically
peak.2016 <-density.df %>%
  filter(date == "20160125")

plot(peak.2016["density"],
     xlim = c(1836725 +680, 1836725 + 820), ylim = c(569768-150, 569768-80))



test.dataset <- uas.dataset %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex)) %>% 
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup")


#fixed it

#now we need to filter out the seals that are categorized as (N & water) not complete in the p_complete column
uasdata <- dataset1 %>%
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup")