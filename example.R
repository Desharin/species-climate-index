# Libraries
library(raster)
source("species_climate_parameter_index.R")

# Data
temp = tempfile()
download.file("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip", temp)
file = unzip(temp, "wc2.1_30s_elev.tif")
elevation_raster = raster(file)
unlink(temp)

# Calculate elevation index for species 4280004 (Polytrichum alpinum)
result = gbif_climate_parameter_index("4280004", elevation_raster, user = "userName", pwd = "somePassword", email = "email@example.com")

result
# species                    parameter     mean       sd    perc5   perc25   perc75   perc95    fraction_NA nObs
# Polytrichastrum alpinum wc2.1_30s_elev 1611.062 465.1332 1139.438 1446.758 1822.856 1987.048           0  298
