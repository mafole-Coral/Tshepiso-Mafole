#Simplified climate projections (beginner level)
#Data source: WorldClim
#June 2019

#*********************************************************************************************************************

################---------------------------------Future climate projection-------------------------------------------

#********************************************************************************************************************* 
# can get aridity index for past and future climates here: https://envirem.github.io/ 

# Define lists for GCMs, SSP scenarios, and periods
gcms <- c("BCC-CSM2-MR", "MPI-ESM1-2-HR", "MRI-ESM2-0")
ssps <- c("ssp245")  # Add more scenarios if needed
periods <- c("2061-2080")  # Add more periods if needed

# Loop over GCMs
for (gcm in gcms) {
  # Loop over SSP scenarios
  for (ssp in ssps) {
    # Loop over periods
    for (period in periods) {
      # Construct the file paths
      file_path <- paste0("GIS/WorldClim/cmip6_2070/", ssp, "/wc2.1_10m_bioc_", gcm, "_", ssp, "_", period, ".tif")
      
      # Read the environmental predictors
      gcm_stack = brick(x = file_path)
      gcm_stack = gcm_stack[[c(3, 4, 8, 9, 12, 18, 19)]]
      
      # Stack with aridity index
      AI_file = paste0("AI_", tolower(gcm))
      AI_stack = raster::stack(paste0("/GIS/WorldClim/aridity_Index/AI_245/envirem/", 
                                      AI_file, ".grd"))
      gcm_stack = raster::stack(gcm_stack, AI_stack[[1]])
      
      # Change raster extent to match with (constant) soil stack
      gcm_stack = crop(gcm_stack, extent(soilStack))
      
      # Combine with constant soil variables
      gcm_env = raster::stack(gcm_stack, soilStack)
      
      # Match layer names with "current climate" layers
      names(gcm_env) = names(envStack)
      
      # Remove correlated variables and crop to Med. regions
      gcm_env_uncor = raster::subset(envStackMR, names(gcm_env)[!names(gcm_env) %in% excl])
      
      # Run the model
      outDir = paste0('R/MaxEnt_Out/future/', ssp, '/', genus, '_', spEpithet, '_', gcm)
      if (!file.exists(outDir)) { dir.create(outDir, showWarnings = FALSE) } 
      
      source("R/Scripts/predictionModel_function.R")
      gcm_prj = model_projection_function(gcm_env_uncor)
      
      # Calculate means for predicted probabilities
      model_stack = stack(gcm_prj[[2]])
      meanPP = calc(model_stack, fun = mean)
      sd = calc(model_stack, fun = sd)
      cv = (sd / meanPP) * 100
    }
  }
}


#*********************************************************************************************************************

################---------------------------------Past climate projection---------------------------------------------

#*********************************************************************************************************************        
# Define the models
models <- list(
  paste0('cclgmbi', c('3', '4', '8', '9', '12', '15')),
  paste0('melgmbi', c('3', '4', '8', '9', '12', '15')),
  paste0('mrlgmbi', c('3', '4', '8', '9', '12', '15'))
)

# Loop through each model
for (i in seq_along(models)) {
  #----------------- 1. get environmental predictors- worldclim
  # stack bioclim variables
  folder <- "GIS/WorldClim/lgmClim/"
  model <- raster::stack(paste0(folder, models[[i]], '.tif'))
  
  # aridity index
  AI <- raster(paste0("GIS/WorldClim/aridity_Index/lgm/", tolower(models[[i]]), "_10arcmin_aridityIndexThornthwaite.tif"))
  model <- raster::stack(model, AI)
  
  # add onto soil stack
  model <- raster::stack(model, soilStack)
  
  # remove correlated variables and crop to Med. Regions
  model_Uncor <- raster::subset(envStackMR, names(model)[!names(model) %in% excl])
  
  #----------------- 2. run model
  outDir <- paste0('R/MaxEnt_Out/lgm/', genus, '_', spEpithet)
  if (!file.exists(outDir)) { dir.create(outDir, showWarnings = FALSE) }
  
  source("R/Scripts/predictionModel_lgm_function.R")
  model_Uncor <- model_projection_function(model_Uncor)
  
  #----------------- 3. get prediction mean
  modelStack <- stack(model_Uncor[[2]])
  meanPP <- calc(modelStack, fun = mean)
  sd <- calc(modelStack, fun = sd)
  cv <- (sd / meanPP) * 100
  plot(meanPP)
}

###----------------------------------------------------------------    