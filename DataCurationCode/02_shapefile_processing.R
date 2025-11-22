
# Load libraries and data -------------------------------------------------


library(dplyr)
install.packages("lubridate")
library(lubridate)
library(sf)

options(scipen = 999)

dates <- list.files("RawData/")




years <- substr(dates, 0, 4)

uas.dataset <- data.frame()

for(i in 1:length(dates)) {
  filepath <- paste0("RawData/", dates[i], "/polygons.shp")


  uas.dataset <- rbind(uas.dataset, 
                       cbind(st_read(filepath), date = dates[i], year = years[i]))
  
}



# Convert to NAD83 / California zone 3 -------------------------------------------------------------


st_crs(uas.dataset)

uas.dataset = st_transform(uas.dataset,"EPSG:26943") #transform shape file 

st_crs(uas.dataset)

plot(st_geometry(uas.dataset))

# Calculate area ----------------------------------------------------------

uas.dataset$area <- as.numeric(st_area(uas.dataset))


# Calculate seal length and width -----------------------------------------

source("DataCurationCode/Functions/1_Seal_Volume_Function_MHM_fixed.R")



uas.dataset$length <- NA
uas.dataset$width <- NA



for(j in 1:length(st_geometry(uas.dataset))){
  pol1 = uas.dataset[j,]
  lens = curved_length_vol(pol1, plt = TRUE)
  uas.dataset[j,7:8] = lens[1:2]
  print(j)
} 



# Save outputs ------------------------------------------------------------

uas.data <- uas.dataset |>
  st_drop_geometry()

write.csv(uas.data, "IntermediateData/uasdata.csv",
          row.names = FALSE)

