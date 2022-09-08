#############################
### Species Climate Index ###
#############################

# This function calculates a climate parameter index for given occurrence lattitudes and longitudes for the given climate raster file.
climate_parameter_index = function(
    decimalLatitude, # Vector of decimal latitudes of observations of a single species
    decimalLongitude, # Vector of decimal longitudes of observations of a single species
    climateRaster, # RasterLayer from the package Raster containing the climate parameter to calculate, this layer needs to be 43200 pixels wide and 21600 pixels high such as the WorldClim 30s rasters
    name = NULL, # Row name used for the results
    nLoop = 100, # Number of times the bootstrap is repeated
    removeDuplicates = T # Whether multiple observations in the same pixel should be removed to reduce observation bias
  ) {
  
  # Load libraries
  require(raster)
  
  # Check if parameters are in the correct format
  stopifnot(is.numeric(decimalLatitude))
  stopifnot(is.numeric(decimalLongitude))
  stopifnot(class(climateRaster)[1] == "RasterLayer")
  stopifnot(is.numeric(nLoop))
  stopifnot(is.logical(removeDuplicates))
  stopifnot(ncol(climateRaster) == 43200)
  stopifnot(nrow(climateRaster) == 21600)
  
  # Convert coordinates to pixels
  pixels = (floor((decimalLongitude + 180) * 2 * 60) + 1) + # Longitude
    floor((-decimalLatitude + 90) * 2 * 60) * # Latitude
    43200 # Number of rows
  # Remove occurrences in the same pixel to reduce observation bias
  if(removeDuplicates) {
    duplications = duplicated(pixels)
    pixels = pixels[!duplications]
  } else {
    duplications = F
  }
  
  # Get values from the raster file for the given pixels
  values = extract(climateRaster, pixels)
  values[!is.na(values) & values == 0] = NA
  
  # Each grid cell is identified by the x and y of the bottom left corner
  # Create 250x250 grid cells
  cell250x = seq(from = -25, to = 42.5, by = 2.5)
  cell250y = seq(from = 30, to = 70, by = 2.5)
  cell250 = expand.grid(x = cell250x, y = cell250y)
  # Assign each observation to a 50x50km grid cell
  x50 = floor(decimalLongitude[!duplications] / 0.5) * 0.5
  y50 = floor(decimalLatitude[!duplications] / 0.5) * 0.5
  c50 = paste0(x50, ", ", y50)
  
  additions = c(0, 0.5, 1, 1.5, 2)
  bootstrap.result = NULL
  
  # The bootstrap loop
  for(i in 1:nLoop) {
    # Get random 50x50 cell for each 250x250 cell
    x.additions = sample(additions, length(cell250$x), replace = T)
    y.additions = sample(additions, length(cell250$y), replace = T)
    # Get all occurrences in the selected cells
    selected_values = values[c50 %in% paste0((cell250$x + x.additions), ", ", (cell250$y + y.additions))]
    
    # Calculate mean, sd and quantiles for the selected occurrences
    result.mean = mean(selected_values, na.rm = T)
    result.sd = sd(selected_values, na.rm = T)
    result.perc5 = quantile(selected_values, probs = 0.05, na.rm = T, names = F)
    result.perc25 = quantile(selected_values, probs = 0.25, na.rm = T, names = F)
    result.perc75 = quantile(selected_values, probs = 0.75, na.rm = T, names = F)
    result.perc95 = quantile(selected_values, probs = 0.95, na.rm = T, names = F)
    result.fraction_NA = sum(is.na(selected_values)) / length(selected_values)
    
    bootstrap.result = rbind.data.frame(bootstrap.result, data.frame(mean = result.mean, sd = result.sd, perc5 = result.perc5, perc25 = result.perc25, perc75 = result.perc75, perc95 = result.perc95, fraction_NA = result.fraction_NA))
  }
  
  # Get mean over the samples
  return(data.frame(mean = mean(bootstrap.result$mean, na.rm = T),
                    sd = mean(bootstrap.result$sd, na.rm = T),
                    perc5 = mean(bootstrap.result$perc5, na.rm = T),
                    perc25 = mean(bootstrap.result$perc25, na.rm = T),
                    perc75 = mean(bootstrap.result$perc75, na.rm = T),
                    perc95 = mean(bootstrap.result$perc95, na.rm = T), 
                    fraction_NA = mean(bootstrap.result$fraction_NA, na.rm = T),
                    row.names = name
                    )
  )
  
}

# This function automatically retrieves occurrence data for given taxa from GBIF 
# and calculates climate parameter indices for the given climate raster files.
gbif_climate_parameter_index = function(
    taxonKeys, # Either a single or a list of GBIF taxonKey character strings for which the occurrence data will be retrieved
    climateRasters, # Either a single or a list of RasterLayers from the package Raster. This layer needs to be 43200 pixels wide and 21600 pixels high
    minYear = 1970, # The minimum year for which occurrences will be retrieved
    nLoop = 1000, # The number of bootstrap loops
    user = NULL, # Your GBIF username used to retrieve the occurrences, also see Authentication below
    pwd = NULL, # Your GBIF password used to retrieve the occurrences, also see Authentication below
    email = NULL # Your GBIF email used to retrieve the occurrences, also see Authentication below
    
    ## Authentication
    # For user, pwd, and email parameters, you can set them in one of three ways:
    # - Set them in your .Rprofile file with the names gbif_user, gbif_pwd, and gbif_email
    # - Set them in your .Renviron/.bash_profile (or similar) file with the names GBIF_USER, GBIF_PWD, and GBIF_EMAIL
    # - Simply pass strings to each of the parameters in the function call
    # We strongly recommend the second option - storing your details as environment variables as it’s the most widely used way to store secrets.
    
  ) {
  
  # Libraries
  require(rgbif)
  require(raster)
  
  # Check parameters
  stopifnot(is.character(taxonKeys))
  if(class(climateRasters) == "RasterLayer") {
    climateRasters = c(climateRasters)
  }
  stopifnot(class(climateRasters) == "list")
  for(climateRaster in climateRasters) {
    stopifnot(class(climateRaster) == "RasterLayer")
  }
  stopifnot(is.numeric(minYear))
  stopifnot(is.numeric(nLoop))
  
  # Prepare GBIF download
  print("Preparing GBIF download...")
  
  res = occ_download(
    pred_in("taxonKey", taxonKeys),
    pred_in("basisOfRecord", c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "OCCURRENCE", "PRESERVED_SPECIMEN", "LIVING_SPECIMEN")),
    pred_gte("coordinateUncertaintyInMeters", 0),
    pred_lte("coordinateUncertaintyInMeters", 10000),
    pred("hasCoordinate", T),
    pred("hasGeospatialIssue", F),
    pred("occurrenceStatus", "PRESENT"),
    pred_gte("year", minYear),
    user = user,
    pwd = pwd,
    email = email
  )
  
  print("Awaiting GBIF download...")
  
  occ_download_wait(res[1])
  download_data <- occ_download_get(res[1]) %>%
    occ_download_import()
  
  # Make simple Dataframe with important columns
  data = data.frame(id = download_data$gbifID, species = download_data$species, lat = download_data$decimalLatitude, lon = download_data$decimalLongitude)
  
  # Calculate climate paramets
  print("Calculating climate parameters...")
  
  result = NULL
  
  for(species in unique(data$species)) {
    for(climateRaster in climateRasters) {
      new_result = climate_parameter_index(data$lat[data$species == species], data$lon[data$species == species], climateRaster, nLoop = nLoop, name = species)
      row.names(new_result) = NULL
      result = rbind.data.frame(result, data.frame(species = species, parameter = names(climateRaster), new_result, nObs = nrow(data[data$species == species,])))
    }
  }
  
  print("Done!")
  
  return(result)
}