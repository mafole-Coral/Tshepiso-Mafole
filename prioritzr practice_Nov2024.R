#Spatial planning with 'prioritzr'
#November 2024
#Package practice run

sessionInfo()


library(prioritizr)
library(prioritizrdata)
library(Rsymphony)
library(terra)
library(sf)

#A story about Mali, a region in West Africa
#*\Useful resources
#1. Study area: https://www.cbd.int/countries/profile?country=ml
#2. Package paper: https://doi.org/10.1111/cobi.14376

#*******************************************************************************
#load species data
acacia = read.csv('~/Desktop/data/0030794-241107131044228.csv', sep = ',', header = T) 


#Convert occurrence data to an sf object
acacia = acacia[!is.na(acacia$decimalLongitude) & !is.na(acacia$decimalLatitude), ]
occ_sf = st_as_sf(acacia, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)


#*******************************************************************************
#define Mali boundary
mali = vect("~/Desktop/data/mali_shapefile.shp")
planning_units = rast(ext(mali), resolution = 0.1) 
planning_units = rasterize(mali, planning_units)

# #define the extent of Mali and create the planning unit raster
# mali_extent = ext(mali_vect)
# planning_units = rast(mali_extent, resolution = 0.1)
# values(planning_units) = NA 

# #rasterize Mali boundary onto the planning units
# mali_raster = rasterize(mali_vect, planning_units, field = 1, background = NA)
# plot(mali_raster)

#rasterize species occurrences
acacia_raster = rasterize(occ_sf, planning_units, field = 1, background = 0)

#add protected areas
protected_areas = vect("~/Desktop/data/WDPA_WDOECM_Nov2024_Public_MLI_shp_0/WDPA_WDOECM_Nov2024_Public_MLI_shp-polygons.shp")
p_areas_raster = rasterize(protected_areas, planning_units)

#*******************************************************************************
#*https://hub.worldpop.org/geodata/summary?id=42376


#load and process population density (cost layer)
#https://energydata.info/dataset/mali-republic-population-density-2015/resource/5388ca5d-e70f-4998-966c-95f5264b5178
pop_density = terra::rast("~/Desktop/data/mli-popmlipppv2b2015unadj/MLI_ppp_v2b_2015_UNadj.tif")
pop_density_mali = crop(pop_density, mali)
pop_density_mali = mask(pop_density_mali, mali)

#align CRS
pop_density_mali = project(pop_density_mali, crs(planning_units))
acacia_raster = project(acacia_raster, crs(planning_units))

#align resolution and extent
pop_density_mali = resample(pop_density_mali, planning_units, method = "bilinear")
acacia_raster = resample(acacia_raster, planning_units, method = "near")

#crop and mask to ensure exact alignment
pop_density_mali = crop(pop_density_mali, ext(planning_units))
pop_density_mali = mask(pop_density_mali, planning_units)

acacia_raster = crop(acacia_raster, ext(planning_units))
acacia_raster = mask(acacia_raster, planning_units)


#assign population density as the cost layer
values(planning_units) = values(pop_density_mali)
summary(planning_units)
max_cost = max(values(planning_units), na.rm = T)
values(planning_units)[is.na(values(planning_units))] = max_cost * 10

#feature stack
features_stack = c(acacia_raster, p_areas_raster)

#create the problem
p = prioritizr::problem(planning_units, features = features_stack) %>%
  add_min_set_objective() %>%             #minimize costs while meeting targets
  add_relative_targets(c(0.4, 0.6)) %>%   #protect 40% of Acacia senegal habitat
  add_boundary_penalties(penalty = 1) %>% #encourage clustering
  add_rsymphony_solver(gap = 0.1, verbose = T)


solution = solve(p, force=T)
plot(solution, main = "Conservation Prioritisation for Mali")

