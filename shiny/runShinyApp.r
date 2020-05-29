setwd("C:/Users/Cornholio/Desktop/FOSSFlood-master")    # This should look something like C:/Users/.../FOSSFlood-master
basedir <- getwd()
#  _ __ ___  __ _  _____  __
# | '__/ _ \/ _` |/ _ \ \/ /
# | | |  __/ (_| |  __/>  < 
# |_|  \___|\__, |\___/_/\_\
#           |___/           
# Here is where the HTA fired VBS replaces strings for user interactions.  At the moment if you want to run FOSSFlood on another platform you can
# Install R studio
# Point R studio the R-Portable executible in the FOSSFlood installation
# Replace code below with desired inputs
# run the entire script

# -- USER Inputs -------------------------------------------------------------------------
UserZipcode <- "66044, 66046, 66047, 66045, 66049"
cachedCheck <- "##cachedCheck" # Unused
UserAddressess <- "OAAdds"
AddressessFile <- ""
UserRoads <- "TIGER"
RoadsFile <- ""
ForecastSource <- "NWM_SR_C"
ForecastFile <- ""
GridChoice <- "Square"
OutputType <- "BVOSMOTM"
PublishSite <- "##PublishSite"  # Unused
PublishSite <- "##PrintOutput"  # Unused
realProjection <- TRUE  # Unused
useOldMethod <- FALSE  

# -- Dev   Comment/uncomment with ctrl-shift-c
UserZipcode <- "66044, 66046, 66047, 66045, 66049"
cachedCheck <- "##cachedCheck" # Unused
UserAddressess <- "OAAdds"
AddressessFile <- ""
UserRoads <- "TIGER"
RoadsFile <- ""
ForecastSource <- "NWM_SR_C"
ForecastFile <- ""
GridChoice <- "Square"
OutputType <- "BVOSMOTM"
PublishSite <- "##PublishSite"  # Unused
PublishSite <- "##PrintOutput"  # Unused
realProjection <- TRUE  # Unused
useOldMethod <- FALSE  

print(paste("-- Welcome to FOSSFlood - Running FOSSFlood in", basedir)) # This should look something like C:/Users/.../FOSSFlood-master
print("-- Welcome to FOSSFlood - Pre-Preloading constants --")

##  _______       ___      .__   __.   _______  _______ .______          ________    ______   .__   __.  _______ 
## |       \     /   \     |  \ |  |  /  _____||   ____||   _  \        |       /   /  __  \  |  \ |  | |   ____|
## |  .--.  |   /  ^  \    |   \|  | |  |  __  |  |__   |  |_)  |       `---/  /   |  |  |  | |   \|  | |  |__   
## |  |  |  |  /  /_\  \   |  . `  | |  | |_ | |   __|  |      /           /  /    |  |  |  | |  . `  | |   __|  
## |  '--'  | /  _____  \  |  |\   | |  |__| | |  |____ |  |\  \----.     /  /----.|  `--'  | |  |\   | |  |____ 
## |_______/ /__/     \__\ |__| \__|  \______| |_______|| _| `._____|    /________| \______/  |__| \__| |_______|
##  
## Edit below this line at your own risk, well, as risky as coding can be...
## Note: If durring your edits FOSSFlood breaks, the easiest way to recover is simply to download FOSSFlood again

# install.packages('installr')
# install.packages('devtools')
# install.packages('dismo')
# install.packages('ggplot2')
# install.packages('ncdf4')
## install.packages('FedData')
# install.packages('ggmap')
# install.packages('curl')
# install.packages('RCurl')
# install.packages('stringr')
# install.packages('gtools')
# install.packages('rvest')
# install.packages('tigris')
# install.packages('sp')
# install.packages('sf')
# install.packages('noncensus')
# install.packages('maptools')
# install.packages('rgdal')
# install.packages('raster')
# install.packages('dplyr')
# install.packages('leaflet')
# install.packages('shapefiles')
# install.packages('httr')
# install.packages('rgeos')
# install.packages('shiny')
# install.packages('shinydashboard')
# install.packages('data.table')
# install.packages('DT')
# install.packages('RColorBrewer')
# install.packages('lattice')
# install.packages('scales')
# install.packages('rio')
# install.packages('cdlTools')
# install.packages('dataRetrieval')
# install.packages('MazamaSpatialUtils')
# install.packages("stargazer")
# install.packages("sjPlot")
# install.packages("formattable")
# install.packages("DT")
# install.packages("tibble")
# install.packages("flexdashboard")
# install.packages("mapview")
# install.packages("geosphere")
# devtools::install_github("mikejohnson51/AOI")
# devtools::install_github("mikejohnson51/nwm")
# devtools::install_github("mikejohnson51/HydroData")
# devtools::install_github("mikejohnson51/FloodMapping")
suppressMessages(library(dismo))
suppressMessages(library(ggplot2))
suppressMessages(library(ncdf4))
# suppressMessages(library(FedData))
suppressMessages(library(ggmap))
suppressMessages(library(curl))
suppressMessages(library(RCurl))
suppressMessages(library(stringr))
suppressMessages(library(gtools))
suppressMessages(library(rvest))
suppressMessages(library(tigris))
suppressMessages(library(sp))
suppressMessages(library(sf))
suppressMessages(library(noncensus))
suppressMessages(library(maptools))
suppressMessages(library(rgdal))
suppressMessages(library(raster))
suppressMessages(library(dplyr))
suppressMessages(library(leaflet))
suppressMessages(library(shapefiles))
suppressMessages(library(httr))
suppressMessages(library(rgeos))
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(data.table))
suppressMessages(library(AOI))
suppressMessages(library(nwm))
suppressMessages(library(HydroData))
suppressMessages(library(FloodMapping))
suppressMessages(library(DT))
suppressMessages(library(RColorBrewer))
suppressMessages(library(lattice))
suppressMessages(library(scales))
suppressMessages(library(rio))
suppressMessages(library(cdlTools))
suppressMessages(library(dataRetrieval))
suppressMessages(library(stargazer))
suppressMessages(library(sjPlot))
suppressMessages(library(MazamaSpatialUtils))
suppressMessages(library(formattable))
suppressMessages(library(DT))
# suppressMessages(library(tibble))
suppressMessages(library(mapview))
suppressMessages(library(flexdashboard))
suppressMessages(library(geosphere))

rasterOptions(datatype="FLT4S")
options(tigris_use_cache = FALSE)

# Helper functions
make_grid <- function(x, type, cell_width, cell_area, clip = FALSE) {
  if (missing(cell_width)) {
    if (missing(cell_area)) {
      stop("Must provide cell_width or cell_area")
    } else {
      if (type == "square") {
        cell_width <- sqrt(cell_area)
      } else if (type == "hexagonal") {
        cell_width <- sqrt(2 * cell_area / sqrt(3))
      }
    }
  }
  # buffered extent of study area to define cells over
  ext <- as(extent(x) + cell_width, "SpatialPolygons")
  projection(ext) <- projection(x)
  # generate grid
  if (type == "square") {
    g <- raster(ext, resolution = cell_width)
    g <- as(g, "SpatialPolygons")
  } else if (type == "hexagonal") {
    # generate array of hexagon centers
    g <- spsample(ext, type = "hexagonal", cellsize = cell_width, offset = c(0, 0))
    # convert center points to hexagons
    g <- HexPoints2SpatialPolygons(g, dx = cell_width)
  }
  
  # clip to boundary of study area
  if (clip) {
    g <- gIntersection(g, x, byid = TRUE)
  } else {
    g <- g[x, ]
  }
  # clean up feature IDs
  row.names(g) <- as.character(1:length(g))
  
  # Get polygons centroids
  centroids <- as.data.frame(centroid(g))
  colnames(centroids) <- c("lon", "lat") 
  centroids <- centroids %>% 
    mutate_if(is.numeric, round, digits = 5)
  lonGroups <- sort(unique(centroids$lon),decreasing = FALSE)
  lonGDF <- as.data.frame( lonGroups )
  lonGDF$LonID <- seq.int(nrow(lonGDF))
  latGroups <- sort(unique(centroids$lat),decreasing = TRUE)
  latGDF <- as.data.frame( latGroups )
  latGDF$ID <- seq.int(nrow(latGDF))
  
  mergeTable <- read.csv(paste0(basedir,"/data/misc/IndexReclassTable.csv"))  
  
  latIndex <- merge(x = latGDF, y = mergeTable, by.x = "ID", by.y = 'From', all.x = TRUE)
  
  centroids <- merge(x = centroids, y = latIndex, by.x = "lat", by.y = 'latGroups', all.x = TRUE)
  names(centroids)[names(centroids)=="ID"] <- "lanNumID"
  names(centroids)[names(centroids)=="To"] <- "LatID"
  centroids <- merge(x = centroids, y = lonGDF, by.x = "lon", by.y = 'lonGroups', all.x = TRUE)
  
  coordinates(centroids) <- c("lon","lat")
  proj4string(centroids) <- proj4string(x) # assign projection
  
  g <- as(g, "SpatialPolygonsDataFrame")
  PolyTransfer <- over(g, centroids)  # Get district data
  PolyTransfer <- mutate(PolyTransfer, IndexLabel = str_replace_all(paste(LatID, LonID, sep = '-'), "'", ""))
  g <- spCbind(g, PolyTransfer)
  g$SparceID <- seq.int(nrow(g))
  return(g)
}
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
getMostRecentForecast <- function() {
  fulltimestamp <- format(Sys.time(), tz = "GMT", format = "%Y-%m-%d %H:%M:%S") 
  fileDate <- gsub("-", "",format(Sys.time(), tz = "GMT", format = "%Y-%m-%d"))
  fileTime <- format(Sys.time(), tz = "GMT", format = "%H")
  
  # -- Check to make sure I have the right date --------------------------------------------------------------------------- 
  testurl <- paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate, '/short_range')
  if (http_error(testurl)) { 
    fileDate <- gsub("-", "",format(as.Date(fulltimestamp, tz = "GMT")-1, format = "%Y-%m-%d"))
  }
  # -- Guess at the most recent time ---------------------------------------------------------------------------------------
  t <- 25
  repeat {
    t = t-1
    rtime <- sprintf("%02d", t)
    testurl <- paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t', rtime, 'z.short_range.channel_rt.f018.conus.nc')
    if(!http_error(testurl)) {
      break
    }
  }
  # coreFileName = fileDate_rtime_fileDate_rtime=fhour
  return (c(fileDate,rtime))
}

forecastTimestamp <- getMostRecentForecast()

print("-- Welcome to FOSSFlood - Pre-Preloading complete --")

# Run only if this is the first time you have run FOSSFlood, unpacks and creates the zip code file
if (!file.exists(paste0(basedir, "/data/misc/Zipcodes/ZIPCodes.shp"))) {
  print("-- Welcome to FOSSFlood - Unpacking FOSSFlood for use --")
  unzip(paste0(basedir, "/data/misc/HUC6.zip"), exdir = paste0(basedir, "/data/misc"))
  unzip(paste0(basedir, "/data/misc/UTM.zip"), exdir = paste0(basedir, "/data/misc"))
  dir.create(paste0(basedir, "/data/misc/Zipcodes"))
  
  FirstZip <- zctas(cb = TRUE)
  UTMShapefile <- readOGR(paste0(basedir, "/data/misc/UTM/UTM.shp"))
  
  ZipCodeUTM <- intersect(FirstZip, UTMShapefile)
  
  writeOGR(obj=ZipCodeUTM, dsn=paste0(basedir, "/data/misc/Zipcodes"), layer='ZIPCodes', driver="ESRI Shapefile")
  
  #Other option, need to figure this out at some point
  # GET('https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/WBD_National_GDB.zip', write_disk(paste0(basedir, "/data/misc/WBD_National_GDB.zip")), overwrite=TRUE)
  # setSpatialDataDir(paste0(basedir, "/data"))
  # installSpatialData(urlBase = "http://mazamascience.com/RData/Spatial",
  #                 file = "mazama_spatial_files-0.5.tar.gz")
  
}

UserZipCodeList <- as.list(strsplit(UserZipcode, ", ")[[1]])
UserZipCodeFileName <- gsub(", ", "_", UserZipcode)
print(paste("-- Welcome to FOSSFlood - Running FOSSFlood for", UserZipcode))
NeedBaseData = file.exists(paste0(basedir,"/AOI/", UserZipCodeFileName, "/hydro/rating.csv"))

if (!NeedBaseData) {
  #/////////////////////////////////////
  # -- Build file structure
  #/////////////////////////////////////
  # -- Build file structure -----------------------------------------------------------------------------------------------------------------------------------------------------
  print("-- Base data collection - This may take upwards of 2 hours to run - Building File list")
  dir.create(paste0(basedir,"/AOI"))
  dir.create(paste0(basedir,"/AOI/", UserZipCodeFileName))
  
  dir.create(paste0(basedir,"/AOI/",UserZipCodeFileName,"/output/"))
  dir.create(paste0(basedir,"/AOI/",UserZipCodeFileName,"/output/roads/"))
  dir.create(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/"))
  dir.create(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp"))
  dir.create(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses"))
  dir.create(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads"))
  dir.create(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/grids"))
  dir.create(paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/"))
  
  # -- I first need an AOI shapefile to run Mike Johnson packages 
  ZipCodeShapefile <- readOGR(paste0(basedir, "/data/misc/ZIPCodes/ZIPCodes.shp"), verbose = FALSE)
  aoiZIPCodes <- subset(ZipCodeShapefile, (ZipCodeShapefile$ZCTA5CE10 %in% UserZipCodeList))
  ZipcodeCRS <- CRS(paste0("+proj=utm +zone=", aoiZIPCodes$ZONE, "+datum=WGS84"))
  
  # -- Run Mikes code
  setwd(paste0(basedir,"/AOI/",UserZipCodeFileName))
  MyAOI = getAOI(clip = aoiZIPCodes)
  getData(MyAOI, name = "tmp")
  setwd(basedir)
  
  # -- Error check, did his stuff exit correctly?
  NeedOldProcess = file.exists(paste0(basedir,"/AOI/", UserZipCodeFileName, "/LIVINGFLOOD/tmp/hydro/rating.rda"))
  
  #/////////////////////////////////////
  # -- If it fails, backup
  #/////////////////////////////////////
  if(useOldMethod || !NeedOldProcess) {
    ## The code in this loop is old and likely broken.  File names will need to be changed at the very least
    ## This is here more for archival reasons, but should be updated so that it is more robust
    
    
    print("-- Base data collection - This may take 4+ hours to run - Collecting HUC6")
    #Pull in zip code shape file(), subset the list, grab the first zone, and transform the AOI
    print(paste0("-- Base data collection - looking for ", basedir, "/data/misc/ZIPCodes/ZIPCodes.shp"))
    ZipCodeShapefile <- readOGR(paste0(basedir, "/data/misc/ZIPCodes/ZIPCodes.shp"))
    SelectedZipCode <- subset(ZipCodeShapefile, (ZipCodeShapefile$ZCTA5CE10 %in% UserZipCodeList))
    #newcrs <- CRS(paste0("+proj=utm +zone=", SelectedZipCode$ZONE, " +datum=WGS84")) 
    #SelectedZipCode <- spTransform(SelectedZipCode, newcrs)
    
    #pull in huc 6 datasets, transform, subset, and collect desired HUC from ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/Hydrography/WBD/National/GDB/
    HUC6Shapefile <- readOGR(paste0(basedir, "/data/misc/HUC6/HUC6.shp"))
    #HUC6Shapefile <- spTransform(HUC6Shapefile, newcrs)
    HUC6Shapefile <- HUC6Shapefile[SelectedZipCode,]
    HUC6Labels <- as.character(HUC6Shapefile$HUC6)
    
    # -- Build URLs --
    print("-- Base data collection - This may take 4+ hours to run - Downloading Base data")
    hand <- list() 
    catchment <- list()
    rate <- list()
    flodbf <- list()
    floprj <- list()
    floshp <- list()
    floshx <- list()
    for(i in 1:length(HUC6Labels)){
      hand[[i]] <- c(paste0(HUC6Labels[i],"/", HUC6Labels[i], "hand.tif"))
      catchment[[i]] <- c(paste0(HUC6Labels[i],"/",HUC6Labels[i], "catchmask.tif"))
      rate[[i]] <- c(paste0(HUC6Labels[i],"/hydroprop-fulltable-",HUC6Labels[i],".csv"))
      flodbf[[i]] <- c(paste0(HUC6Labels[i],"/",HUC6Labels[i], "-flows.dbf"))
      floprj[[i]] <- c(paste0(HUC6Labels[i],"/",HUC6Labels[i], "-flows.prj"))
      floshp[[i]] <- c(paste0(HUC6Labels[i],"/",HUC6Labels[i], "-flows.shp"))
      floshx[[i]] <- c(paste0(HUC6Labels[i],"/",HUC6Labels[i], "-flows.shx"))
    }
    
    # -- NIFE Data download --
    print("-- Base data collection - This may take 4+ hours to run - Downloading HAND")
    for(i in 1:length(hand)){
      URL <- paste0("https://web.corral.tacc.utexas.edu/nfiedata/HAND/", hand[i])
      GET(URL, write_disk(paste0(getwd(),"/data/HAND/",substr(hand[i],8,nchar(hand[i]))), overwrite=TRUE))
    }
    print("-- Base data collection - This may take 4+ hours to run - Downloading Catchment")
    for(i in 1:length(catchment)){
      URL <- paste0("https://web.corral.tacc.utexas.edu/nfiedata/HAND/", catchment[i])
      GET(URL, write_disk(paste0(getwd(),"/data/Catchments/",substr(catchment[i],8,nchar(catchment[i]))), overwrite=TRUE))
    }
    print("-- Base data collection - This may take 4+ hours to run - Downloading Rating Curves")
    for(i in 1:length(rate)){
      URL <- paste0("https://web.corral.tacc.utexas.edu/nfiedata/HAND/", rate[i])
      GET(URL, write_disk(paste0(getwd(),"/data/rate/",substr(rate[i],8,nchar(rate[i]))), overwrite=TRUE))
    }
    print("-- Base data collection - This may take 4+ hours to run - Downloading Flow Lines")
    for(i in 1:length(flodbf)){
      URL <- paste0("https://web.corral.tacc.utexas.edu/nfiedata/HAND/", flodbf[i])
      GET(URL, write_disk(paste0(getwd(),"/data/shp/flow/",substr(flodbf[i],8,nchar(flodbf[i]))), overwrite=TRUE))
    }
    for(i in 1:length(floprj)){
      URL <- paste0("https://web.corral.tacc.utexas.edu/nfiedata/HAND/", floprj[i])
      GET(URL, write_disk(paste0(getwd(),"/data/shp/flow/",substr(floprj[i],8,nchar(floprj[i]))), overwrite=TRUE))
    }
    for(i in 1:length(floshp)){
      URL <- paste0("https://web.corral.tacc.utexas.edu/nfiedata/HAND/", floshp[i])
      GET(URL, write_disk(paste0(getwd(),"/data/shp/flow/",substr(floshp[i],8,nchar(floshp[i]))), overwrite=TRUE))
    }
    for(i in 1:length(floshx)){
      URL <- paste0("https://web.corral.tacc.utexas.edu/nfiedata/HAND/", floshx[i])
      GET(URL, write_disk(paste0(getwd(),"/data/shp/flow/",substr(floshx[i],8,nchar(floshx[i]))), overwrite=TRUE))
    }
    
    # -- Current Merge Shapefile workflow -----------------------------------------------------------------------------------------------------------------------------------------------------
    print("-- Base data collection - This may take 4+ hours to run - Merging Shapefiles")
    setwd(paste0(basedir,"/AOI/", UserZipCodeFileName,"/data/shp/flow/"))
    file_list <- as.list(dir(getwd(), "*.shp"))
    for (ffile in file_list) {
      # if the merged dataset does exist, append to it
      if (exists("MergedFlow")) {
        temp_dataset <-readOGR(ffile)
        MergedFlow <- rbind(MergedFlow, temp_dataset)
        rm(temp_dataset)
      }# if the merged dataset doesn't exist, create it
      if (!exists("MergedFlow")) { 
        MergedFlow <- readOGR(ffile) 
      }
    }
    MergedFlow <- sp::spTransform(MergedFlow, CRS("+init=epsg:4326"))
    writeOGR(obj=MergedFlow, dsn=getwd(), layer='FF_MergedFlows', driver="ESRI Shapefile")
    FF_MergedFlow <- readOGR(paste0(getwd(), "/FF_MergedFlows.shp"))
    LMergedFlows <- FF_MergedFlow[SelectedZipCode,]
    # -- Might as well pull COMID's while we're at it -------------------------------
    comids <- vector(mode = "numeric",length(LMergedFlows))
    comids <- LMergedFlows@data$COMID
    
    # -- Current Merge Raster workflow -----------------------------------------------------------------------------------------------------------------------------------------------------
    print("-- Base data collection - This may take 4+ hours to run - Merging HAND Rasters")
    setwd(paste0(basedir,"/AOI/", UserZipCodeFileName,"/data/HAND/"))
    HANDFiles <- list.files(paste0(basedir,"/AOI/", UserZipCodeFileName,"/data/HAND/"), pattern = ".tif$", full.names = TRUE)
    if(length(HANDFiles) == 1) {
      HANDLayer <- raster(HANDFiles)
      mHAND <- HANDLayer
    } else {
      HANDLayer <- lapply(HANDFiles, raster)
      mHAND <- do.call(raster::merge, c(HANDLayer, tolerance = 30))
    }
    
    print("-- Base data collection - This may take 4+ hours to run - projecting HAND Rasters")
    mtHANDr <- projectRasterForLeaflet(mHAND, method = "ngb")
    #mtHANDr <- spTransform(mtHANDr, CRS("+init=epsg:3857"))
    
    print("-- Base data collection - This may take 4+ hours to run - Trimming HAND Rasters")
    SelectedZipCodes <- spTransform(SelectedZipCode, CRS("+init=epsg:3857"))
    mcHAND <- raster::crop(mtHANDr, extent(SelectedZipCodes), snap="in", filename="MC_HAND.tif", overwrite=TRUE)     #https://gis.stackexchange.com/questions/61243/clipping-a-raster-in-r
    HANDRaster <- raster(paste0(basedir,"/AOI/", UserZipCodeFileName, "/data/HAND/MC_HAND.tif"))
    
    # -- Current Merge Raster workflow -------------
    rasterOptions(datatype="INT4U")
    setwd(paste0(basedir,"/AOI/", UserZipCodeFileName,"/data/Catchments/"))
    CatchFiles <- list.files(paste0(basedir,"/AOI/", UserZipCodeFileName,"/data/Catchments/"), pattern = ".tif$", full.names = TRUE)
    if(length(CatchFiles) == 1) {
      catchLayer <- raster(CatchFiles)
      mCatch <- catchLayer
    } else {
      catchLayer <- lapply(CatchFiles, raster)
      mCatch <- do.call(raster::merge, c(catchLayer, tolerance = 30))
    }
    #writeRaster(mCatch, filename="Merged_Catchments.tif", format="GTiff")
    
    print("-- Base data collection - This may take 4+ hours to run - projecting Catchment Rasters")
    mtCatch <- projectRasterForLeaflet(mCatch, method = "ngb")
    #writeRaster(mtCatch, filename="P_Merged_Catchments.tif", format="GTiff")
    
    # -- resample catchment to hand raster -----
    print("-- Base data collection - This may take 4+ hours to run - Resampling Catchment Rasters")
    mtCatchNew <- raster::resample(mtCatch, HANDRaster, "ngb")
    #writeRaster(mtCatch, filename="Resamp_Catchments.tif", format="GTiff")
    
    # -- save, crop, save, and open  -----
    print("-- Base data collection - This may take 4+ hours to run - Trimming Catchment Rasters")
    #mtCatchNewRas <- raster(mtCatchNew)
    #inRast <- raster(paste0(basedir,"/AOI/", UserZipCodeFileName, "/data/Catchments/Resamp_Catchments.tif"))
    mcCatch <- raster::crop(mtCatchNew,  extent(SelectedZipCodes), snap="in", filename="MC_Catchments.tif")
    #mcCatchnew = mask(mcCatch, HANDRaster)
    writeRaster(mcCatch, filename="MC_Catchments.tif", format="GTiff", overwrite=TRUE)
    
    CatchmentsRaster <- raster(paste0(basedir,"/AOI/", UserZipCodeFileName, "/data/Catchments/MC_Catchments.tif"))
    #file.remove(paste0(basedir,"/AOI/", UserZipCodeFileName, "/data/Catchments/Resamp_Catchments.tif"))
    # rasterOptions(datatype="FLT4S")
    
    
    # -- Create/consolidate Rating curves -----------------------------------------------------------------------------------------------------------------------------------------------------
    print("-- Base data collection - This may take 4+ hours to run - Creating Ratings")
    rating_curves = list()
    setwd(paste0(basedir,"/AOI/", UserZipCodeFileName))
    for (i in 1:length(dir(paste0(basedir,"/AOI/", UserZipCodeFileName, "/data/rate")))) {
      build <-  read.csv(paste0(basedir,"/AOI/", UserZipCodeFileName, "/data/rate/", dir(paste0(basedir,"/AOI/", UserZipCodeFileName, "/data/rate"))[i]))
      build <- cbind(build[,1], build[,15], build[,2])
      build <- build[comids %in% build[,1],]
      rating_curves[[i]] <- build
    }
    all_rating_curves <- do.call(rbind, rating_curves)
    colnames(all_rating_curves) <- c("COMIDs", "Discharge_cms", "Stage")
    write.csv(all_rating_curves, file <- paste0(basedir,"/AOI/", UserZipCodeFileName, "/data/rate/rating_curves.csv"))
    #ratingCurves <- read.csv(paste0(basedir, "/AOI/", UserZipCodeFileName, "/data/rate/rating_curves.csv"))
  }
  #/////////////////////////////////////
  #/////////////////////////////////////
  #/////////////////////////////////////
  
  
  #/////////////////////////////////////
  # -- Reproject and move to FOSSFlood file structure
  #/////////////////////////////////////
  print("-- Base data collection - This may take upwards of 2 hours to run - reprojecting and restructuring raster files")
  CatchRasProj <- leaflet::projectRasterForLeaflet(raster(paste0(basedir,"/AOI/",UserZipCodeFileName,"/LIVINGFLOOD/tmp/geo/catchmask.tif")), method = "ngb")
  writeRaster(CatchRasProj, filename=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/catchmask.tif"),format="GTiff", overwrite=TRUE)
  
  HandRasProj <- leaflet::projectRasterForLeaflet(raster(paste0(basedir,"/AOI/",UserZipCodeFileName,"/LIVINGFLOOD/tmp/geo/hand.tif")), method = "bilinear")
  writeRaster(HandRasProj, filename=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/hand.tif"),format="GTiff", overwrite=TRUE)
  
  print("-- Base data collection - This may take upwards of 2 hours to run - reprojecting and restructuring vector files")
  nhdFlows <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/LIVINGFLOOD/tmp/geo/nhdflowlines.shp"), verbose = FALSE)
  nhdFlowsProj <- sp::spTransform(nhdFlows, CRS("+init=epsg:4326"))
  writeOGR(obj=nhdFlowsProj, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo"), layer='nhdflowlines', driver="ESRI Shapefile")
  
  AOIshp <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/LIVINGFLOOD/tmp/geo/AOI.shp"), verbose = FALSE)
  AOIshpProj <- sp::spTransform(AOIshp, CRS("+init=epsg:4326"))
  writeOGR(obj=AOIshpProj, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo"), layer='AOI', driver="ESRI Shapefile")
  
  aoiZIPCodesProj <- sp::spTransform(aoiZIPCodes, CRS("+init=epsg:4326"))
  writeOGR(obj=aoiZIPCodesProj, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo"), layer='AOIZIPcodes', driver="ESRI Shapefile")
  
  print("-- Base data collection - This may take upwards of 2 hours to run - Building rating curves")
  file.rename(from=paste0(basedir,"/AOI/",UserZipCodeFileName,"/LIVINGFLOOD/tmp/hydro/rating.rda"), to=paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/rating.rda")) #Maybe this one?
  ratingCurves <- load(paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/rating.rda"))
  write.csv(rating_curves[,c(1,2,3)], file=paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/rating.csv"), row.names=F)
  
  unlink(paste0(basedir,"/AOI/",UserZipCodeFileName,"/LIVINGFLOOD"), recursive=TRUE)
  
  
  #/////////////////////////////////////
  # -- Grab backup roads and points 
  #/////////////////////////////////////
  # TIGRIS Roads from Census
  # Dev notes: here I grab the counties and use roads package to download roads (from Census)
  # This is redundent with the OSM data I download for buildings defaults.  
  # There is also no handling of zip codes which straddle state lines
  print("-- Base data collection - This may take upwards of 2 hours to run - Grabbing TIGER road files")
  CountiesFIPS <- quiet(counties(cb=TRUE))
  CountiesFIPSProj <- sp::spTransform(CountiesFIPS, CRS("+init=epsg:4326"))
  aoiCounties <- CountiesFIPSProj[aoiZIPCodesProj, ]
  # Error checking: Download roads for AOI that interset more than one state or county
  if(length(unique(aoiCounties$COUNTYFP)) > 1) {
    myHoldRoads <- quiet(roads(unique(aoiCounties$STATEFP),unique(aoiCounties$COUNTYFP)[1], year = 2014, refresh = TRUE)[NULL,])
    for(i in unique(aoiCounties$COUNTYFP)) {
      myRoads <- quiet(roads(unique(aoiCounties$STATEFP),i, year = 2014, refresh = TRUE))
      myHoldRoads <- rbind(myHoldRoads, myRoads)
    }
    myRoads <- myHoldRoads
  } else {
    myRoads <- quiet(roads(unique(aoiCounties$STATEFP),unique(aoiCounties$COUNTYFP), year = 2014, refresh = TRUE))
  }
  myRoadsProj <- sp::spTransform(myRoads, CRS("+init=epsg:4326"))
  myRoads_subset <- myRoadsProj[aoiZIPCodesProj, ]
  writeOGR(obj=myRoads_subset, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads"), layer="tigerroads", driver="ESRI Shapefile")
  
  # OSM Roads and points
  print("-- Base data collection - This may take upwards of 2 hours to run - Grabbing OSM data")
  osmDataURL <- paste0("https://download.geofabrik.de/north-america/us/",tolower(gsub(" ", "-", fips(aoiCounties$STATEFP[1], to = "Name"))),"-latest-free.shp.zip")
  download.file(osmDataURL, paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp/",tolower(gsub(" ", "-", fips(aoiCounties$STATEFP[1], to = "Name"))),"-latest-free.shp.zip"))
  unzip(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp/",tolower(gsub(" ", "-", fips(aoiCounties$STATEFP[1], to = "Name"))),"-latest-free.shp.zip"), exdir = paste0(basedir,"/AOI/",UserZipCodeFileName, "/geo/tmp"))
  osmPoints <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp/gis_osm_pois_free_1.shp"), verbose = FALSE)
  osmPointsProj <- sp::spTransform(osmPoints, CRS("+init=epsg:4326"))
  writeOGR(obj=osmPointsProj[aoiZIPCodesProj, ], dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses"), layer='OSMaddresses', driver="ESRI Shapefile")
  osmRoads <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp/gis_osm_roads_free_1.shp"), verbose = FALSE)
  osmRoadsProj <- sp::spTransform(osmRoads, CRS("+init=epsg:4326"))
  writeOGR(obj=osmRoadsProj[aoiZIPCodesProj, ], dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads"), layer='OSMroads', driver="ESRI Shapefile")
  
  # OpenAddresses points
  print("-- Base data collection - This may take upwards of 2 hours to run - Grabbing OpenAddresses data")
  if(aoiCounties$STATEFP %in% c('09','23','25','33','44','50','34','36','42')) {
    oaURL <- "https://data.openaddresses.io/openaddr-collected-us_northeast.zip"
  } else if(aoiCounties$STATEFP %in% c('18','17','26','39','55','19','20','27','29','31','38','46')) {
    oaURL <- "https://data.openaddresses.io/openaddr-collected-us_midwest.zip"
  } else if(aoiCounties$STATEFP %in% c('10','11','12','13','24','37','45','51','54','01','21','28','47','05','22','40','48')) {
    oaURL <- "https://data.openaddresses.io/openaddr-collected-us_south.zip"
  } else if(aoiCounties$STATEFP %in% c('04','08','16','35','30','49','32','56','02','06','15','41','53')) {
    oaURL <- "https://data.openaddresses.io/openaddr-collected-us_west.zip"
  }
  download.file(oaURL, paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp/OpenAddresses.zip"))
  unzip(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp/OpenAddresses.zip"), exdir = paste0(basedir,"/AOI/",UserZipCodeFileName, "/geo/tmp"))
  FullOpenAddresses <- paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp/us/",tolower(fips(aoiCounties$STATEFP[1], to = "Abbreviation")),"/statewide.csv")
  FullOpenAddressesCSV<- read.csv(FullOpenAddresses, header=TRUE, stringsAsFactors = FALSE)
  coordinates(FullOpenAddressesCSV)<- ~LON+LAT
  proj4string(FullOpenAddressesCSV) <- CRS("+init=epsg:4326")
  writeOGR(obj=FullOpenAddressesCSV[aoiZIPCodesProj, ], dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses"), layer='OpenAddresses', driver="ESRI Shapefile")
  
  unlink(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/tmp"), recursive = TRUE)
  
  # -- Make some ancillary products
  print("-- Base data collection - This may take upwards of 2 hours to run - generating fishnets")
  hex_grid <- make_grid(AOIshpProj, type = "hexagonal", cell_area = 0.00005, clip = FALSE)
  rec_grid <- make_grid(AOIshpProj, type = "square", cell_area = 0.00005, clip = FALSE)
  writeOGR(obj=hex_grid, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/grids"), layer='HexFishnet', driver="ESRI Shapefile")
  writeOGR(obj=rec_grid, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/grids"), layer='SquareFishnet', driver="ESRI Shapefile")
  
  # Finally, pull in NWIS gages
  print("-- Base data collection - This may take upwards of 2 hours to run - Pulling and populating NWIS gages")
  # Get active gages  - Double check this
  AOIGages <- MyAOI %>% findNWIS()
  load(file = paste0(basedir,"/data/misc/usgs_filter.rda"))
  GagesWithCOMID <- merge(AOIGages$nwis, usgs_filter %>% select('siteID', 'COMID'), by.x = 'site_no', by.y = 'siteID')
  GagesWithCOMIDsProj <- sp::spTransform(GagesWithCOMID, CRS("+init=epsg:4326"))
  writeOGR(obj=GagesWithCOMIDsProj, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro"), layer='USGSGages', driver="ESRI Shapefile")
  
  print("-- Base data collection finished - Loading in data")
  #/////////////////////////////////////
  #Copy from below!
  #/////////////////////////////////////
  if(UserAddressess == "OAAdds") {
    Addresses <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses/OpenAddresses.shp"), verbose = FALSE)
  } else if (UserAddressess == "OSMAdds") {
    Addresses <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses/OSMaddresses.shp"), verbose = FALSE)
  } else if(UserAddressess == "UPAAdds") {
    tryCatch(Addresses <- readOGR(AddressessFile, verbose = FALSE),
             error = function(e) {
               print("-- !ALERT! Invalid addresses provided, defaulting to OpenAddresses database --")
               Addresses <<- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses/OSMaddresses.shp"), verbose = FALSE)
             }
    )
  }
  if(UserRoads == "OSM") {
    Road <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads/OSMroads.shp"), verbose = FALSE)
  } else if (UserRoads == "TIGER") {
    Road <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads/tigerroads.shp"), verbose = FALSE)
  } else if(UserRoads == "UPRoads") {
    tryCatch(Road <- readOGR(RoadsFile, verbose = FALSE),
             error = function(e) {
               print("-- !ALERT! Invalid addresses provided, defaulting to TIGRIS database --")
               Road <<- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads/tigerroads.shp"), verbose = FALSE)
             }
    )
  }
  if(GridChoice == "Square") {
    MapIndex <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/geo/grids/SquareFishnet.shp"), verbose = FALSE)
  } else if(GridChoice == "Hexagonal") {
    MapIndex <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/geo/grids/HexFishnet.shp"), verbose = FALSE)
  }
  Flowlines <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/geo/nhdflowlines.shp"), verbose = FALSE)
  Gages <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/hydro/USGSGages.shp"), verbose = FALSE)
  AOI <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/geo/AOIZIPcodes.shp"), verbose = FALSE)
  #if(baseData){
  if(TRUE) { 
    ZipCodeShapefile <- readOGR(paste0(basedir, "/data/misc/ZIPCodes/ZIPCodes.shp"), verbose = FALSE)
    HUC6Shapefile <- readOGR(paste0(basedir, "/data/misc/HUC6/HUC6.shp"), verbose = FALSE)
    UTMShapefile <- readOGR(paste0(basedir, "/data/misc/UTM/UTM.shp"))
    catchMask <- raster(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/catchmask.tif"), verbose = FALSE)
    handRas <- raster(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/hand.tif"), verbose = FALSE)
  }
} else {
  #/////////////////////////////////////
  #Loading in base data
  #/////////////////////////////////////
  print("-- Loading in base data --")
  if(UserAddressess == "OAAdds") {
    Addresses <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses/OpenAddresses.shp"), verbose = FALSE)
  } else if (UserAddressess == "OSMAdds") {
    Addresses <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses/OSMaddresses.shp"), verbose = FALSE)
  } else if(UserAddressess == "UPAAdds") {
    tryCatch(Addresses <- readOGR(AddressessFile, verbose = FALSE),
             error = function(e) {
               print("-- !ALERT! Invalid addresses provided, defaulting to OpenAddresses database --")
               Addresses <<- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/addresses/OSMaddresses.shp"), verbose = FALSE)
             }
    )
  }
  if(UserRoads == "OSM") {
    Road <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads/OSMroads.shp"), verbose = FALSE)
  } else if (UserRoads == "TIGER") {
    Road <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads/tigerroads.shp"), verbose = FALSE)
  } else if(UserRoads == "UPRoads") {
    tryCatch(Road <- readOGR(RoadsFile, verbose = FALSE),
             error = function(e) {
               print("-- !ALERT! Invalid addresses provided, defaulting to TIGRIS database --")
               Road <<- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/roads/tigerroads.shp"), verbose = FALSE)
             }
    )
  }
  if(GridChoice == "Square") {
    MapIndex <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/geo/grids/SquareFishnet.shp"), verbose = FALSE)
  } else if(GridChoice == "Hexagonal") {
    MapIndex <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/geo/grids/HexFishnet.shp"), verbose = FALSE)
  }
  Flowlines <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/geo/nhdflowlines.shp"), verbose = FALSE)
  Gages <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/hydro/USGSGages.shp"), verbose = FALSE)
  AOI <- readOGR(paste0(basedir,"/AOI/", UserZipCodeFileName, "/geo/AOIZIPcodes.shp"), verbose = FALSE)
  #if(baseData){
  if(TRUE) { 
    ZipCodeShapefile <- readOGR(paste0(basedir, "/data/misc/ZIPCodes/ZIPCodes.shp"), verbose = FALSE)
    HUC6Shapefile <- readOGR(paste0(basedir, "/data/misc/HUC6/HUC6.shp"), verbose = FALSE)
    UTMShapefile <- readOGR(paste0(basedir, "/data/misc/UTM/UTM.shp"))
    catchMask <- raster(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/catchmask.tif"), verbose = FALSE)
    handRas <- raster(paste0(basedir,"/AOI/",UserZipCodeFileName,"/geo/hand.tif"), verbose = FALSE)
  }
}

#/////////////////////////////////////
#Get forecasts from NWM
#/////////////////////////////////////
#Note: Some of this needs to be cleaned up
setwd(paste0(basedir,"/AOI/", UserZipCodeFileName))
print("-- Base data collected -- Checking for Updated Forecasts")
useNOMADS <<- FALSE
if(ForecastSource == "A&A") {
  # -- Grab A & A ----
  print("-- Grabbing NWM forecasts -- Analysis & Assimilation")
  tryCatch(getFlows(name.dir = getwd(), config = "short_range", n = 1),
           error = function(e) {
             print("-- Grabbing NWM forecasts -- Hydroshare THREDDS unreachable, grabbing forecast from NOMADS")
             useNOMADS <<- TRUE
           }
  )
} else if(ForecastSource == "NWM_SR_C") {
  print("-- Grabbing NWM forecasts -- Short Range Forecast")
  tryCatch(getFlows(name.dir = getwd(), config = "short_range", f = 18, n = 18),
           error = function(e) {
             print("-- Grabbing NWM forecasts -- Hydroshare THREDDS unreachable, grabbing forecast from NOMADS")
             useNOMADS <<- TRUE
           }
  )
} else if(ForecastSource == "MR") {
  # -- Grab Medium range ----
  print("-- Grabbing NWM forecasts -- Medium Range Forecast")
  tryCatch(getFlows(name.dir = getwd(), config = "short_range", date = 20190224, f = 10),
           error = function(e) {
             print("-- Grabbing NWM forecasts -- Hydroshare THREDDS unreachable, grabbing forecast from NOMADS")
             useNOMADS <<- TRUE
           }
  )
} else if(ForecastSource == "UPStage") {
  # called ForecastFile above
  print(paste("-- Mapping User Provided Forecast -- Pulling data from /AOI/", UserZipCodeFileName, "/hydro/flows.csv"))
  setwd(paste0(basedir,"/AOI/", UserZipCodeFileName, "/hydro/"))
  discharge <- as.matrix(read.csv(paste0(basedir, "/AOI/", UserZipCodeFileName, "/hydro/flows.csv")))
  rio::convert(discharge, "flows.rda")
  useNOMADS <<- FALSE
} else if(ForecastSource == "UPDischarge") {
  # called ForecastFile above
  print(paste("-- Mapping User Provided Forecast -- Pulling data from /AOI/", UserZipCodeFileName, "/hydro/flows.csv"))
  setwd(paste0(basedir,"/AOI/", UserZipCodeFileName, "/hydro/"))
  discharge <- as.matrix(read.csv(paste0(basedir, "/AOI/", UserZipCodeFileName, "/hydro/flows.csv")))
  rio::convert(discharge, "flows.rda")
  useNOMADS <<- FALSE
}
setwd(basedir)

#/////////////////////////////////////
#NOMADS backup in case hydroshare is wonk
#/////////////////////////////////////
#Note: Some of this needs to be cleaned up or written...
if(useNOMADS) {
  dir.create(paste0(basedir,"/AOI/NWM"), showWarnings = FALSE)
  
  # -- Pick a discharge source ----
  if(ForecastSource == "A&A") {
    # -- Grab A & A ----
  } else if(ForecastSource == "NWM_SR_C") {
    # # -- Pick out the key timestamp variables -------------------------------------------------------------------------------
    fulltimestamp <- format(Sys.time(), tz = "GMT", format = "%Y-%m-%d %H:%M:%S")
    fileDate <- gsub("-", "",format(Sys.time(), tz = "GMT", format = "%Y-%m-%d"))
    fileTime <- format(Sys.time(), tz = "GMT", format = "%H")
    
    # -- Check to make sure I have the right date ---------------------------------------------------------------------------
    testurl <- paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate, '/short_range')
    if (http_error(testurl)) {
      fileDate <- gsub("-", "",format(as.Date(fulltimestamp, tz = "GMT")-1, format = "%Y-%m-%d"))
    }
    # -- Guess at the most recent time ---------------------------------------------------------------------------------------
    t <- 25
    repeat {
      t = t-1
      rtime <- sprintf("%02d", t)
      testurl <- paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t', rtime, 'z.short_range.channel_rt.f018.conus.nc')
      if(!http_error(testurl)) {
        break
      }
    }
    
    # -- Remove old files --  (------------------------------------------------------------------------------ 
    file.remove(list.files(paste0(basedir,"/AOI/NWM/"), pattern = ".nc", full.names = TRUE))
    file.remove(list.files(paste0(basedir,"/AOI/", UserZipCodeFileName,"/hydro"), pattern = "flows.csv", full.names = TRUE))
    file.remove(list.files(paste0(basedir,"/AOI/", UserZipCodeFileName,"/output"), pattern = ".tif", full.names = TRUE))
    
    # -- Grab NWM ------------------------------------------------------------------------------------------------------------
    print("-- Grabbing NWM forecasts -- Downloading Forecasts: This may take a few minutes")
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  1 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f001.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f001.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  2 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f002.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f002.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  3 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f003.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f003.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  4 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f004.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f004.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  5 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f005.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f005.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  6 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f006.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f006.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  7 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f007.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f007.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  8 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f008.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f008.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour  9 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f009.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f009.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 10 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f010.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f010.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 11 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f011.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f011.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 12 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f012.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f012.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 13 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f013.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f013.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 14 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f014.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f014.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 15 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f015.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f015.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 16 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f016.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f016.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 17 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f017.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f017.conus.nc'), overwrite=TRUE))
    print("-- Grabbing NWM forecasts -- Downloading Forecast for hour 18 of 18")
    GET(paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t',rtime, 'z.short_range.channel_rt.f018.conus.nc'), write_disk(paste0(basedir,"/AOI/NWM/nwm.t",rtime, 'z.short_range.channel_rt.f018.conus.nc'), overwrite=TRUE))
    
    # --Grab the COMIDs ------------------------------------------------------------------------------------------------------------
    print("-- Grabbing NWM forecasts -- Subsetting NWM to selected reaches")
    comids <- vector(mode = "numeric",length(Flowlines))
    comids <- Flowlines@data$comid
    
    # --Subset NWM to reaches ------------------------------------------------------------------------------------------------------------
    folder <- paste0(basedir,"/AOI/NWM/")
    files <- dir(folder)
    cstart.index <- 1
    cstart.date <- gsub("-","", Sys.Date())
    end.index <- length(files)
    end.date <- cstart.date
    cstart.time <- NULL
    files <- files[cstart.index:end.index]
    nc <- nc_open(filename = paste0(folder,"/", files[[1]]))
    vars <- nc$var
    test1test2 <- vars$streamflow$dim[[1]]$vals
    comids_of_value <- comids[comids %in% test1test2]
    cstart <- vector(mode = "numeric", length(comids_of_value))
    for(i in 1:length(comids_of_value)) {
      cstart[i] <- which(test1test2 == comids_of_value[i])
    }
    nwm.flow = matrix(NA, ncol = length(files), nrow = length(cstart))
    nc_close(nc)
    for (i in 1:length(files)) {
      nc = nc_open(filename = paste0(folder,"/", files[i]))
      values = ncvar_get(nc, varid = "streamflow")
      for(j in 1:length(cstart)) {
        nwm.flow[j,i] = values[cstart[j]]
      } 
      nc_close(nc)
    }
    rownames(nwm.flow) <- comids_of_value
    colnames(nwm.flow) <- substr(files, 1, 12)
    nwm.flow = nwm.flow * 35.3147
    write.csv(nwm.flow, paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/flows.csv"))
    discharge <- read.csv(paste0(basedir, "/AOI/", UserZipCodeFileName, "/hydro/flows.csv"))      
    
    colnames(discharge) <- c("COMIDS","timestep 1","timestep 2","timestep 3","timestep 4","timestep 5","timestep 6","timestep 7","timestep 8","timestep 9","timestep 10","timestep 11","timestep 12","timestep 13","timestep 14","timestep 15","timestep 16","timestep 17","timestep 18")
    write.csv(discharge, file=paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/flows.csv"), row.names = FALSE)
    save(discharge, file=paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/flows.rda"))
    
  } else if(ForecastSource == "MR") {
    print("Not added Yet")
  }
  
  unlink(paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/NWM"), recursive=TRUE)
}

#/////////////////////////////////////
# This is a hack around the map function from floodmapping, 
# take forecasts of discharge and rating curves, spits out rasters of flood polygons
#/////////////////////////////////////
print("-- Creating Flood Inundation Polygons --")

name.dir = paste0(basedir,"/AOI/", UserZipCodeFileName)
write = TRUE
add = 0

j = NULL
`%dopar%` = foreach::`%dopar%`

all.files = list.files(name.dir, recursive = T, full.names = TRUE)

flows <- import(paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/flows.rda"))
rating_curves <- import(paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/rating.rda"))

catchmentv = velox::velox(all.files[grepl("catchmask.tif", all.files)])
handv = velox::velox(all.files[grepl("hand.tif", all.files)])

if(length(unique(rating_curves$CatchId)) > length(unique(flows$COMIDS))){
  index = rating_curves$CatchId %in% flows$COMIDS
  rating_curves = rating_curves[index,]
} else if(length(unique(rating_curves$CatchId)) < length(unique(flows$COMIDS))) {
  index = flows$COMIDS %in% rating_curves$CatchId
  flows = flows[index,]
}

comids = unique(flows$COMIDS)
stage = NULL

for(i in seq_along(comids)){
  flow = flows[flows$COMIDS == comids[i],] + add
  curve = rating_curves[rating_curves$CatchId == comids[i],]
  fin = NULL
  
  for(j in 2:dim(flow)[2]){
    tmp = curve$Stage[which.min(abs(curve$`Discharge..m3s.1.` - flow[1,j]))]
    if(length(tmp) <= 0){ tmp = NA }
    fin = append(tmp, fin)
  }
  stage = rbind(stage, fin)
}

stage = cbind(comids, stage)
rownames(stage) = NULL
colnames(stage) = c('COMID', paste0("timestep", 1:(dim(stage)[2]-1)))
stage = as.data.frame(stage, stringsAsFactors = FALSE)
write.csv(stage, file=paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/stage.csv"), row.names = FALSE)
save(stage, file=paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/stage.rda"))

catch.v = as.vector(t(catchmentv$rasterbands[[1]]))
hand.v =  as.vector(t(handv$rasterbands[[1]]))

doParallel::registerDoParallel( parallel::detectCores() - 1 )

a <- foreach::foreach(j = 2:NCOL(stage), .combine = raster::stack) %dopar% {
  val.v = stage[fastmatch::fmatch(catch.v, stage$COMID), j]
  fin.v = val.v - hand.v
  #fin.v[fin.v <= 0] <- NA
  #fin.v[fin.v > 0] <- 1
  f.v = matrix(fin.v, ncol = catchmentv$dim[2], byrow = T)
  f = raster::raster(f.v)
  raster::extent(f) <- catchmentv$extent
  raster::crs(f) = catchmentv$crs
  raster::res(f) <- catchmentv$res
  return(f)
}

names(a) = paste0("timestep", c(1:dim(a)[3]))
setwd(paste0(basedir,"/AOI/", UserZipCodeFileName, '/output'))
unlink(getwd(), recursive = TRUE)
writeRaster(a, filename=names(a), bylayer=TRUE,format="GTiff", overwrite=TRUE)

sumGT <- function(x) sum(x > 0)
fFreq <- calc(a, sumGT)
writeRaster(fFreq, filename="ffreq.tif",format="GTiff", overwrite=TRUE)


#/////////////////////////////////////
# Create visual validation data
#/////////////////////////////////////
print("-- Creating graphs --")
# Bring in NWIS inst forecast
goBack <- as.Date(Sys.time(), format="%Y/%m/%d") - as.difftime(2, unit="days")  # Check timezone
GagesrawData <- readNWISuv(Gages$site_no, c('00065', '00060'), startDate = goBack, endDate = "")
NWMflows <- import(paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/flows.rda"))
NWMstage <- import(paste0(basedir,"/AOI/",UserZipCodeFileName,"/hydro/stage.rda"))

# Make the names human readable
colnames(GagesrawData)[colnames(GagesrawData)=="X_00060_00000"] <- "NWIS discharge - cubic feet per second"
colnames(GagesrawData)[colnames(GagesrawData)=="X_00065_00000"] <- "NWIS gage height - feet"
# 
# baseStirng
# 
# 
# Gages
# 
# # Join NWIS data to points
# 

# p <- tidyr::spread(GagesrawData[c('site_no','dateTime','NWIS Discharge - cubic feet per second')], 'dateTime','NWIS Discharge - cubic feet per second')
# MyGagesSeries <- merge(Gages, p, by.x = 'site_no', by.y = 'site_no')
# MyGagesFullSeries <- merge(MyGagesSeries, p, by.x = 'COMID', by.y = 'COMID')                   


#/////////////////////////////////////
# Intersecting and plotting
#/////////////////////////////////////
print("-- Creating impact outputs --")
FloodInnList <- list.files(path = paste0(basedir,"/AOI/", UserZipCodeFileName, '/output'), pattern =  "timestep*", full.names = TRUE)
FloodInnList <- mixedsort(sort(FloodInnList))
ffreqRas <- raster(paste0(basedir,"/AOI/",UserZipCodeFileName,"/output/ffreq.tif"))

print("-- Intersecting Addresses --")
Addresses <- sp::spTransform(Addresses, proj4string(raster(FloodInnList[1])))
for(i in 1:length(FloodInnList)) {
  Addresses <- raster::extract(raster(FloodInnList[i]), Addresses, sp = TRUE)
}
Addresses <- raster::extract(ffreqRas, Addresses,sp = TRUE)
Addresses <- sp::spTransform(Addresses,  CRS("+init=epsg:4326"))
writeOGR(obj=Addresses, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/output"), layer='AddressImpacts', driver="ESRI Shapefile")
addressImpacts <- readOGR(paste0(basedir,"/AOI/",UserZipCodeFileName,"/output/AddressImpacts.shp"), verbose = FALSE)

print("-- Intersecting Roads --")
Road <- sp::spTransform(Road, proj4string(raster(FloodInnList[1])))
roadPoints <- as(Road, "SpatialPointsDataFrame")
roadImpactList <- vector("list", 18)
for(i in 1:18) {
  roadPoints <- raster::extract(raster(FloodInnList[i]), roadPoints,sp = TRUE)
  var <- paste0("timestep",i)
  UniqueListAdd <- unique(roadPoints[roadPoints[[var]] > 0 & !is.na(roadPoints[[var]]),]$FULLNAME)  
  roadImpactList[[i]] <- UniqueListAdd
  
  roadClossureShape <- gBuffer(roadPoints[roadPoints[[var]] > 0 & !is.na(roadPoints[[var]]),], width = 35)
  roadClossureShapeDissolve <- gUnionCascaded(roadClossureShape)
  roadClossureShapePoly <- as(roadClossureShapeDissolve, "SpatialPolygonsDataFrame")
  roadClossureShapePolyProj <- sp::spTransform(roadClossureShapePoly, CRS("+init=epsg:4326"))
  writeOGR(obj=roadClossureShapePolyProj, dsn=paste0(basedir,"/AOI/",UserZipCodeFileName,"/output/roads"), layer=paste0('RoadBuffer',i), driver="ESRI Shapefile")
}
RoadImpactShapes <- list.files(path = paste0(basedir,"/AOI/", UserZipCodeFileName, '/output/roads/'), pattern =  "RoadBuffer.*.shp$", full.names = TRUE)
RoadImpactShapes <- mixedsort(sort(RoadImpactShapes))
Road <- sp::spTransform(Road,  CRS("+init=epsg:4326"))

## convert it to 'sf'
addressImpactsSF = st_as_sf(addressImpacts)
MapIndexSF = st_as_sf(MapIndex)

## transform into a 'data.frame' by removing the geometry
## intersect polygons with points, keeping the information from both
AddressMapbook = st_intersection(MapIndexSF, addressImpactsSF)  
st_geometry(AddressMapbook) = NULL  

FlowlinesSF = st_as_sf(Flowlines)
FlowlinesMapbook = st_intersection(MapIndexSF, FlowlinesSF)  
st_geometry(FlowlinesMapbook) = NULL  

RoadSF = st_as_sf(Road)
RoadMapbook = st_intersection(MapIndexSF, RoadSF)  
st_geometry(RoadMapbook) = NULL  

#  -- Init ------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("-- FOSSFlood Execution complete - Starting Shiny server --")

options(shiny.port = 7777)
print("-- Launching FOSSFlood Viewer --")

#chrome.portable = file.path(paste0(basedir, '/GoogleChromePortable/GoogleChromePortable.exe'))
chrome.portable = file.path(paste0(basedir, '/GoogleChromePortable/App/Chrome-bin/chrome.exe'))
#chrome.portable = file.path(paste0(basedir, '/FirefoxPortable/App/Firefox64/firefox.exe'))
#message('firefox.portable paths:\n', chrome.portable)

launch.browser = function(appUrl, browser.path=chrome.portable) {
  #message('Browser path: ', browser.path)
  shell(sprintf('"%s" --app=%s', browser.path, appUrl))
}

setwd(paste0(basedir,"/shiny"))




library(ggplot2)
library(tidyverse)
library(lattice)
library(sp)
library(leaflet) 
library(htmlwidgets)
library(htmltools)
library(plotly)
# library(dygraphs)
# library(xts)
# library(mapview)
# library(listviewer)
library(mapedit)
library(leafpop)
library(shiny)
library(shinydashboard)


val = as.numeric(c(0:30))
pal = colorNumeric(c("#deebf7", "#9ecae1", "#3182bd"), val, na.color = "transparent")
icon.glyphicon <- makeAwesomeIcon(icon = "glyphicon-remove-sign", markerColor = "blue",
                                  iconColor = "black", library = "glyphicon",
                                  squareMarker =  TRUE)


colour = c("red", "blue", "red", "red", "black", "black")
group = c("C", "I", "M", "O", "S", "U")

RoadPAL <- colorFactor(
  palette = c('blue', 'red', 'black', 'gray', 'gray', 'black'),
  domain = Road$RTTYP
)

# Basedata exploration
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "FOSSFlood V 1.0: Data Exploration"),
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        "Maps", 
        tabName = "maps", 
        icon = icon("globe"),
        menuSubItem("Overview", tabName = "Overview", icon = icon("map"))
      ),
      menuItem(
        "Charts", 
        tabName = "charts", 
        icon = icon("bar-chart"),
        menuSubItem("Address Pages", tabName = "c_Address", icon = icon("area-chart")),
        menuSubItem("Roads Pages", tabName = "c_Roads", icon = icon("area-chart")),
        menuSubItem("Watershed Pages", tabName = "c_pop", icon = icon("area-chart"))
      )
    )
  ),
  dashboardBody(
    tags$style(type = "text/css", "#currentConditions {height: calc(100vh - 150px) !important;}"
               , "#basedata {height: calc(100vh - 150px) !important;}"),
    tabItems(
      tabItem(
        tabName = "Overview",
        fluidRow(
          infoBoxOutput("AddressTotalCount"),
          infoBoxOutput("AddressPeakDamage"),
          infoBoxOutput("Areas")
        ),
        fluidRow(
          box(
            title = "AOI Map",
            collapsible = TRUE,
            width = "100%",
            height = "100%",
            status = "primary",
            leafletOutput("AOI_Map")
          ),
          box(
            title = "DT",
            collapsible = TRUE,
            width = "100%",
            height = "100%",
            status = "primary",
            DT::dataTableOutput("dataTable")
          )
        )
        
      )
    )
  )
  
)



## ===========================================================================================================
## ===========================================================================================================
## ===========================================================================================================
server <- function(input, output) {
  
  output$AddressTotalCount <- renderInfoBox({
    infoBox(
      "Approval", "80%", icon = icon("thumbs-up", lib = "glyphicon"),
      color = "yellow"
    )
  })
  
  output$AddressPeakDamage <- renderInfoBox({
    infoBox(
      "Approval", "40%", icon = icon("thumbs-up", lib = "glyphicon"),
      color = "red"
    )
  })
  
  output$Areas <- renderInfoBox({
    infoBox(
      "Area in Sq Miles", res(ffreqRas)[1] * res(ffreqRas)[2] * sum(ffreqRas[] >= 1, na.rm=TRUE)/2.59e+6, icon = icon("thumbs-up", lib = "glyphicon"),
      color = "red"
    )
  })
  
  output$AOI_Map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
      addProviderTiles('Esri.WorldImagery', group = "WorldImagery") %>%
      addPolylines(data = AOI, color = "#000000", weight = 2, group = "Zip Codes", label = ~GEOID10) %>%
      addPolylines(data = HUC6Shapefile[AOI,], color = "#7B9EC8", weight = 2, group = "HUC 6", label = ~htmlEscape(NAME)) %>%     
      addPolylines(data = Flowlines, color = "#7B9EC8", weight = 1, group = "Flowlines", label = ~htmlEscape(gnis_name)) %>%
      addPolylines(data = Road, color = ~RoadPAL(RTTYP), weight = 1, group = "Roads") %>%
      addPolylines(data = MapIndex, color = "#000000", weight = 2, group = "Map Index") %>%
      addRasterImage(ffreqRas, colors = pal, project = FALSE, group = "Flood Frequency")  %>%
      addMarkers(data=addressImpacts[!is.na(addressImpacts$ffreq) & addressImpacts$ffreq > 0, ], 
                 label = ~as.character(NUMBER), 
                 icon = icon("glyphicon-remove-sign", lib = "glyphicon"),
                 group = "Impacted Addresses", 
                 clusterOptions = markerClusterOptions()) %>%
      fitBounds(extent(AOI)@xmin, extent(AOI)@ymin, extent(AOI)@xmax,extent(AOI)@ymax) %>%
      addLayersControl(
        baseGroups = c("Toner Lite", "WorldImagery"),
        overlayGroups = c("Zip Codes", "HUC 6", "Roads", "Flowlines", "Impacted Addresses", "Flood Frequency"),
        options = layersControlOptions(collapsed = FALSE)
      )
  }) 

  output$dataTable <- DT::renderDataTable({
    DT::datatable(AddressMapbook)
  })
    
  output$ImpactIndex <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
      
      addProviderTiles('Esri.WorldImagery', group = "WorldImagery") %>%
      addPolylines(data = MapIndex, color = "#7B9EC8", weight = 2, group = "Map Index") %>%
      addRasterImage(ffreqRas, colors = pal, project = FALSE, group = "Flood Frequency") 
      addMarkers(data=addressImpacts[!is.na(addressImpacts$ffreq) & addressImpacts$ffreq > 0, ], 
                 label = ~as.character(NUMBER), 
                 icon = icon("glyphicon-remove-sign", lib = "glyphicon"),
                 group = "Impacted Addresses", 
                 clusterOptions = markerClusterOptions()) %>%
      addLayersControl(
        baseGroups = c("Toner Lite", "WorldImagery"),
        overlayGroups = c("Map Index", "Flood Frequency", "Impacted Addresses", "Roads"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  output$imagecapture <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$OpenStreetMap.BlackAndWhite, group = "OSM B&W") %>%     
      addProviderTiles(providers$OpenStreetMap.Mapnik, group = "OSM Map") %>%
      addProviderTiles(providers$CartoDB.DarkMatter, group = "Carto DM") %>%
      addProviderTiles(providers$CartoDB.Positron, group = "Carto Pos") %>%
      addProviderTiles(providers$CartoDB.PositronNoLabels, group = "Carto Pos No Lab") %>%
      addProviderTiles(providers$Stamen.Watercolor, group = "Watercolor") %>%
      addProviderTiles(providers$Stamen.Toner, group = "Toner") %>%
      addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
      addProviderTiles('Esri.WorldImagery', group = "ESRI Satellite") %>%
      addProviderTiles(providers$Esri.WorldTopoMap, group = "ESRI Topo") %>%
      addProviderTiles(providers$Esri.NatGeoWorldMap, group = "ESRI Nat Geo") %>%
      addProviderTiles(providers$Esri.WorldGrayCanvas, group = "ESRI Gray") %>%
      addProviderTiles(providers$Esri.WorldStreetMap, group = "ESRI Street") %>%
      addProviderTiles(providers$Esri.DeLorme, group = "ESRI DeLormo") %>%
      addProviderTiles(providers$HikeBike.HillShading, group = "Hillshade") %>%
      addRasterImage(ffreqRas, colors = pal, project = FALSE, group = "Flood Frequency") 
      addLayersControl(
        baseGroups = c("OSM B&W","OSM Map","Carto DM","Carto Pos","Carto Pos No Lab","Watercolor","Toner","Toner Lite","ESRI Satellite","ESRI Topo","ESRI Nat Geo","ESRI Gray","ESRI Street","ESRI DeLormo","Hillshade"),
        overlayGroups = c("Flood Frequency", "Impacted Addresses", "Roads"),
        options = layersControlOptions(collapsed = FALSE)
      )
  }) 
  
  output$GeographicIndex <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
      addProviderTiles('Esri.WorldImagery', group = "WorldImagery") %>%
      addPolylines(data = AOI, color = "#7B9EC8", weight = 2, group = "Zip Codes") %>%
      addPolylines(data = HUC6Shapefile[AOI,], color = "#7B9EC8", weight = 2, group = "HUC 6") %>%     
      addPolylines(data = Flowlines, color = "#7B9EC8", weight = 2, group = "Flowlines") %>%
      addPolylines(data = Road, color = "#7B9EC8", weight = 2, group = "Roads", label = ~htmlEscape(FULLNAME)) %>%
      addLayersControl(
        baseGroups = c("Toner Lite", "WorldImagery"),
        overlayGroups = c("Zip Codes", "HUC 6", "Roads", "Flowlines"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  output$HydrologicIndex <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
      addProviderTiles('Esri.WorldImagery', group = "WorldImagery")
  })
  
  output$BaseData <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
      addProviderTiles('Esri.WorldImagery', group = "WorldImagery") %>%
      addRasterImage(catchMask, project = FALSE, group = "Cathment Mask") %>%
      addRasterImage(handRas, project = FALSE, group = "HAND") %>%
      addMarkers(data=addressImpacts[!is.na(addressImpacts$timestep4) & addressImpacts$timestep4 > 0, ], 
                 label = ~as.character(NUMBER), 
                 icon = icon("glyphicon-remove-sign", lib = "glyphicon"), 
                 clusterOptions = markerClusterOptions())
  })
  
  output$Frequencymaps <- renderLeaflet({
    leaflet() %>%
      #addProviderTiles(providers$OpenTopoMap, group = "TonerLite") %>%
      addTiles(group = "OSM (default)") %>%
      addProviderTiles(providers$Stamen.Toner, group = "Toner") %>%
      addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
      addPolylines(data = myRoads, color = "#7B9EC8", weight = 1, group = "Roads", lineJoin='round') %>%
      addPolylines(data = myFlowlines, color = "#7B9EC8", weight = 2, group = "NHDLines") %>%
      addPolygons(data = myHexGrid, color = "#000000", weight = 1, group = "girds") 
  })
  
  output$RoadStyle <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
      addProviderTiles('Esri.WorldImagery', group = "WorldImagery")  %>%
      addPolylines(data = Road, color = ~RoadPAL(RTTYP), weight = 1, group = "Roads") 
  })
  
}

shinyApp(ui, server)
