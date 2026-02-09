#now i will try to pick out the seals with the highest density and use them as focals

library(sf)
uasdataset.density <- st_read("IntermediateData/uasdata.density.shp")

uasdata.coords <- read_csv("IntermediateData/uas.dataset.coords.csv")

uasdata.part <- uasdata.coords %>%
  mutate(age_sex = if_else(age_sex == "malfee", "male", age_sex))

uasdata.full <- uasdata.part %>%
  filter(
    is.na(p_complete) | !p_complete %in% c("N", "water"),
    age_sex !="pup")

#this is me attempting to combine "uasdata.full" and "uasdataset.density"
nrow(uasdata.full) == nrow(uasdataset.density)

#combined.uasdata <- bind_cols(
  #uasdata.full,
 # uasdataset.density %>% select(-age_sex, -year, -p_complete, -date, -year, -area, -length, -width,)
#)

#ok it worked, but I'm not sure this is a good practice
#this is supposed to test if the rows match, LOL THEY DON'T.....
all(uasdata.full$area == uasdataset.density$area &
      uasdata.full$length == uasdataset.density$length &
      uasdata.full$width == uasdataset.density$width &
      uasdata.full$age_sex == uasdataset.density$age_sex)

