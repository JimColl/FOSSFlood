#  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.  .----------------.   
# | .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |  
# | |  _________   | || |     ____     | || |    _______   | || |    _______   | || |  _________   | || |   _____      | || |     ____     | || |     ____     | || |  ________    | |  
# | | |_   ___  |  | || |   .'    `.   | || |   /  ___  |  | || |   /  ___  |  | || | |_   ___  |  | || |  |_   _|     | || |   .'    `.   | || |   .'    `.   | || | |_   ___ `.  | |  
# | |   | |_  \_|  | || |  /  .--.  \  | || |  |  (__ \_|  | || |  |  (__ \_|  | || |   | |_  \_|  | || |    | |       | || |  /  .--.  \  | || |  /  .--.  \  | || |   | |   `. \ | |  
# | |   |  _|      | || |  | |    | |  | || |   '.___`-.   | || |   '.___`-.   | || |   |  _|      | || |    | |   _   | || |  | |    | |  | || |  | |    | |  | || |   | |    | | | |  
# | |  _| |_       | || |  \  `--'  /  | || |  |`\____) |  | || |  |`\____) |  | || |  _| |_       | || |   _| |__/ |  | || |  \  `--'  /  | || |  \  `--'  /  | || |  _| |___.' / | |  
# | | |_____|      | || |   `.____.'   | || |  |_______.'  | || |  |_______.'  | || | |_____|      | || |  |________|  | || |   `.____.'   | || |   `.____.'   | || | |________.'  | |  
# | |              | || |              | || |              | || |              | || |              | || |              | || |              | || |              | || |              | |  
# | '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |  
#  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------'   
#
# Welcome to the heart of FOSSFlood, an R script which will do the bulk of the heavy lifting for you.  This header section is where the HTA fired VBS replaces strings for user interactions.  
# While the core of FOSSFlood is platform independent, the HTA was needed to create a user interface and thus FOSSFlood is currently limited to the Windows OS.  However, if you want to
# run FOSSFlood on a different OS, all you need to do is...
# * Install R studio
# * In RStudio, point the R installation to the FOSSFlood R-Portable executable (FOSSFlood-master\R-Portable\App\R-Portable\bin\R.exe)
# * Replace code below with desired inputs,
# * and run the entire body of the code :)

rStarttime <- Sys.time()
# setwd("C:/Users/Cornholio/Desktop/FOSSFlood-master")
# setwd("C:/Users/user/Desktop/FOSSFlood-master")
basedir <- getwd()   # This should look something like C:/Users/.../FOSSFlood-master
print(paste("-- Welcome to FOSSFlood - Running in", basedir)) # This should look something like C:/Users/.../FOSSFlood-master
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- User inputs --------------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
user.aoi.string <- "66044, 66046, 66047, 66045, 66049"
user.aoi.source <- "zctas" 										# "zctas", "huc8", "string"
user.address.source <- "OpenAddresses" 						# "OpenStreetMap_Addresses", "OpenAddresses", "User_Provided_Addresses"
user.address.file <- "" 							# Empty "" or filepath to addresses
user.road.source <- "TIGER_Lines_2018"  							# "TIGER_Lines_2018", "OpenStreetMaps", "User_Provided_Roads"
user.road.file <- "" 									# Empty or filepath to addresses
user.forecast.source <- "NWM_SR_C" 						# "NWM_SR_C", USER_DIS_CMS, USER_STAGE
user.forecast.timesteps <- as.numeric("6") 	# number
user.forecast.file <- "" 							# Empty "" or filepath to flows.fst
user.forecast.start <- "##USERFORECASTSTART"						# unused - DEV (ex: "2017-08-23")
user.forecast.end <- "##USERFORECASTEND"							# unused - DEV (ex: "2017-08-31")
user.forecast.members <- as.list("##USERFORECASTMEMBERS") 			# unused - DEV
user.output.choice <- "impacts" 							# GIS_O  basedata  impacts
user.output.grid <- "Square" 								# "Square", "Hexagon"
user.output.hardclip <- as.logical("True") 					# If TRUE, Hard clip data to aoi shape, defaults to bb
user.output.archive <- as.logical("False") 						# Save flows in output folder of the requested AOI using timestamp as file name.  Can be pointed back to later to regenerate outputs
user.output.serverFiles <- as.logical("False")
# Example
#user.aoi.string <- "77479"
#user.aoi.source <- "zctas" 										# "zctas", "huc8", "string"
#user.address.source <- "OpenStreetMap_Addresses" 						# "OpenStreetMap_Addresses", "OpenAddresses", "User_Provided_Addresses"
#user.address.file <- "" 							# Empty "" or filepath to addresses
#user.road.source <- "TIGER_Lines_2018"  							# "TIGER_Lines_2018", "OpenStreetMaps", "User_Provided_Roads"
#user.road.file <- "" 									# Empty or filepath to addresses
#user.forecast.source <- "NWM_SR_C" 						# "NWM_SR_C", USER_DIS_CMS, USER_STAGE
#user.forecast.timesteps <- as.numeric("6") 	# number
#user.forecast.file <- "" 							# Empty "" or filepath to flows.fst
#user.forecast.start <- "2017-08-23"						# unused - DEV
#user.forecast.end <- "2017-08-31"							# unused - DEV
#user.forecast.members <- as.list("##USERFORECASTMEMBERS") 			# unused - DEV
#user.output.choice <- "impacts" 							# GIS_O  basedata  impacts
#user.output.grid <- "Square" 								# "Square", "Hexagon"
#user.output.hardclip <- as.logical("True") 					# If TRUE, Hard clip data to aoi shape, defaults to bb
#user.output.archive <- as.logical("False") 						# Save flows in output folder of the requested AOI using timestamp as file name.  Can be pointed back to later to regenerate outputs
#user.output.serverFiles <- as.logical("False")

# Input cleanup
user.address.source <- stringr::str_replace_all(user.address.source,"_"," ")
user.road.source <- stringr::str_replace_all(user.road.source,"_"," ")



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Preloads (libs) ---------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
print("-- Pre-Preloading constants --")
# Last build date: 3/14/2020
# install.packages("tidyverse")
# install.packages("installr")
# install.packages("devtools")
# install.packages("geosphere")
# install.packages("sf")
# install.packages("randomcoloR")
# install.packages("leaflet")
# install.packages("leafpop")
## install.packages("plotly")
# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("shinyjs")
# install.packages("shinyalert") 
# install.packages("shinycssloaders") 
# install.packages("shinyWidgets")
# install.packages("dataRetrieval")
# install.packages("webshot")
# install.packages("magick")
# install.packages("animation")
# install.packages("imager")
# install.packages("stars")
# install.packages("mapview")
# install.packages("rgeos")
# install.packages("gdalUtilities")
# install.packages("tigris")
# install.packages("MazamaSpatialUtils")
# install.packages("cdlTools")
# install.packages("openadds")
# install.packages("plainview")
# install.packages("leaflet.opacity")
# install.packages("fst")
# install.packages("dygraphs")
# install.packages("xts")
# install.packages("gridExtra")

suppressMessages(library(installr))
suppressMessages(library(devtools))

# url <- "https://cran.r-project.org/src/contrib/Archive/noncensus/noncensus_0.1.tar.gz"
# pkgFile <- "noncensus_0.1.tar.gz"
# download.file(url = url, destfile = pkgFile)
# install.packages(pkgs=pkgFile, type="source", repos=NULL)
# rm(paste0(basedir,"/noncensus_0.1.tar.gz"))

# install.packages("make")  #needed for velox
# url <- "https://cran.r-project.org/src/contrib/Archive/velox/velox_0.2.0.tar.gz"
# pkgFile <- "velox_0.2.0.tar.gz"
# download.file(url = url, destfile = pkgFile)
# install.packages(pkgs=pkgFile, type="source", repos=NULL)
# rm(paste0(basedir,"/noncensus_0.1.tar.gz"))

# devtools::install_github("mikejohnson51/AOI")  # Option 3
# devtools::install_github("mikejohnson51/HydroData")  # Option 3
# devtools::install_github("mikejohnson51/FloodMapping")  # Option 3
# devtools::install_github("hunzikp/velox")  # installs rtools?
# devtools::install_github("ITSLeeds/geofabrik")
## devtools::install_github("mikejohnson51/nomadsNC",auth_token='\')      --- REMOVE ME BEFORE PUSH ----
# devtools::install_github("mikejohnson51/nwmHistoric") 

# install.packages("patchwork")
# install.packages("hrbrthemes")
# install.packages("timeSeries")
# install.packages("TSstudio")
# install.packages("shinybusy")

suppressMessages(library(tidyverse))
suppressMessages(library(geosphere))
suppressMessages(library(sf))
suppressMessages(library(htmlwidgets))
suppressMessages(library(htmltools))
suppressMessages(library(dataRetrieval))
suppressMessages(library(webshot))
suppressMessages(library(magick))
suppressMessages(library(animation))
suppressMessages(library(imager))
suppressMessages(library(stars))
suppressMessages(library(mapview))
# suppressMessages(library(gdalUtilities))
suppressMessages(library(tigris))
# suppressMessages(library(noncensus))
# suppressMessages(library(MazamaSpatialUtils))
suppressMessages(library(cdlTools))
suppressMessages(library(httr))
suppressMessages(library(openadds))
suppressMessages(library(plainview))
suppressMessages(library(fst))
suppressMessages(library(dygraphs))
suppressMessages(library(xts))
suppressMessages(library(timeSeries))
suppressMessages(library(TSstudio))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

suppressMessages(library(randomcoloR))
suppressMessages(library(leaflet))
suppressMessages(library(leafpop))
suppressMessages(library(leaflet.opacity))
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyalert))
suppressMessages(library(shinycssloaders))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinybusy))

suppressMessages(library(AOI))
suppressMessages(library(HydroData))
suppressMessages(library(FloodMapping))
suppressMessages(library(nwmHistoric))
suppressMessages(library(nomadsNC))
suppressMessages(library(velox))
suppressMessages(library(geofabrik))

# suppressMessages(library(patchwork))
# suppressMessages(library(hrbrthemes))
options(tigris_use_cache = FALSE)
apptitle = "FOSSFlood V 1.21"
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Helper Functions ---------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
make_grid <- function(x, type, cell_width, cell_area, clip = FALSE) {
  # Full credit to this function from: https://rpubs.com/dieghernan/beautifulmaps_I
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
  ext <- as(raster::extent(x) + cell_width, "SpatialPolygons")
  raster::projection(ext) <- raster::projection(x)
  # generate grid
  if (type == "square") {
    g <- raster::raster(ext, resolution = cell_width)
    g <- as(g, "SpatialPolygons")
  } else if (type == "hexagonal") {
    # generate array of hexagon centers
    g <- sp::spsample(ext, type = "hexagonal", cellsize = cell_width, offset = c(0, 0))
    # convert center points to hexagons
    g <- sp::HexPoints2SpatialPolygons(g, dx = cell_width)
  }
  
  # clip to boundary of study area
  # First warning: although coordinates are longitude/latitude, st_intersects assumes that they are planar
  if (clip) {
    g <- sf::st_as_sf(g)[x, ]
  } else { g <- sf::st_as_sf(g) }
  # Second warning: st_centroid does not give correct centroids for longitude/latitude data
  centroids = sf::st_coordinates(sf::st_centroid(g))
  colnames(centroids) <- c("lon", "lat")
  
  lonGroups <- sort(unique(centroids[,1]),decreasing = FALSE)
  lonGDF <- as.data.frame( lonGroups )
  lonGDF$LonID <- seq.int(nrow(lonGDF))
  latGroups <- sort(unique(centroids[,2]),decreasing = TRUE)
  latGDF <- as.data.frame( latGroups )
  latGDF$ID <- seq.int(nrow(latGDF))
  
  # alphArr = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
  # i = Input
  # If(i<26){
  #   Print alphArr[i]
  # }else{
  #   //Consider i=27
  #   count = i/26  (here, count=1)
  #   alphabet = i%26  (here alphabet =1)
  #   print alphArr[count]+””+alphArr[alphabet] // Which will be “AA”
  # }
  mergeTable <- read.csv(paste0(basedir,"/data/misc/IndexReclassTable.csv"))  
  
  latIndex <- merge(x = latGDF, y = mergeTable, by.x = "ID", by.y = 'From', all.x = TRUE)
  
  centroids <- merge(x = centroids, y = latIndex, by.x = "lat", by.y = 'latGroups', all.x = TRUE)
  names(centroids)[names(centroids)=="ID"] <- "lanNumID"
  names(centroids)[names(centroids)=="To"] <- "LatID"
  centroids <- merge(x = centroids, y = lonGDF, by.x = "lon", by.y = 'lonGroups', all.x = TRUE)
  
  sp::coordinates(centroids) <- c("lon","lat")
  gs <- as(as_Spatial(g), "SpatialPolygonsDataFrame")
  raster::crs(centroids) <- raster::crs(gs)
  # gs <- as_Spatial(g)
  # st_crs(centroids) = 4326 # assign projection
  # 
  # proj4string(centroids) <- CRS("+proj=longlat +datum=WGS84")
  # 
  # g <- as(g, "SpatialPolygonsDataFrame")
  # spdf <- as_Spatial(g)
  PolyTransfer <- sp::over(gs, centroids)  # Get district data
  PolyTransfer <- mutate(PolyTransfer, IndexLabel = str_replace_all(paste(LatID, LonID, sep = '-'), "'", ""))
  gf <- maptools::spCbind(gs, PolyTransfer)
  gf$SparceID <- seq.int(nrow(gf))
  return(gf)
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
  
  # -- Check to make sure I have the right date --
  testurl <- paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate, '/short_range')
  if (http_error(testurl)) { 
    fileDate <- gsub("-", "",format(as.Date(fulltimestamp, tz = "GMT")-1, format = "%Y-%m-%d"))
  }
  # -- Guess at the most recent time --
  t <- 25
  repeat {
    t = t-1
    zftime <- sprintf("%02d", t)
    testurl <- paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate,'/short_range/nwm.t', zftime, 'z.short_range.channel_rt.f018.conus.nc')
    if(!http_error(testurl)) {
      break
    }
  }
  # coreFileName = fileDate_rtime_fileDate_rtime=fhour
  return (c(fileDate,zftime))
}
getRuntimeTimestamp <- function() {
  fulltimestamp <- format(Sys.time(), tz = "GMT", format = "%Y-%m-%d %H:%M:%S") 
  currentDate <- gsub("-", "",format(Sys.time(), tz = "GMT", format = "%Y-%m-%d"))
  ztime <- format(Sys.time(), tz = "GMT", format = "%H:%M")
  return (c(currentDate,ztime))
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- first run ----------------------------------------------------------------------------------------------
# Block below runs only if this is the first time you have run FOSSFlood, unpacks and creates the zip code 
# files, installs needed dependencies for printing
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
if (!file.exists(paste0(basedir, "/data/misc/zctas.shp"))) {
  print("-- Unpacking FOSSFlood for first use: Downloading zctas --")
  # Create ZCTA
  zctasShape <- quiet(tigris::zctas(cb = FALSE))
  zctasShapeProj <- sf::st_transform(sf::st_as_sf(zctasShape), sf::st_crs(4326))
  sf::write_sf(sf::st_as_sf(zctasShapeProj), paste0(basedir, "/data/misc/zctas.shp"), delete_layer = TRUE, quiet = TRUE)
}
if (!file.exists(paste0(basedir, "/data/misc/huc8.shp"))) {
  print("-- Welcome to FOSSFlood - Unpacking FOSSFlood for first use: Downloading HUC units (this may take a while)--")
  httr::GET('https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/WBD/National/GDB/WBD_National_GDB.zip', write_disk(paste0(basedir, "/data/misc/WBD_National_GDB.zip")), overwrite=TRUE)
  unzip(paste0(basedir, "/data/misc/WBD_National_GDB.zip"), exdir = paste0(basedir, "/data/misc"))
  HUC8Shape <- sf::read_sf(dsn = paste0(basedir, "/data/misc/WBD_National_GDB.gdb"), layer = "WBDHU8")
  HUC8ShapeProj <- sf::st_transform(HUC8Shape, sf::st_crs(4326))
  sf::write_sf(sf::st_as_sf(HUC8ShapeProj), paste0(basedir, "/data/misc/huc8.shp"), delete_layer = TRUE, quiet = TRUE)
  unlink(paste0(basedir, "/data/misc/WBD_National_GDB.gdb"), recursive = TRUE)
  unlink(paste0(basedir, "/data/misc/WBD_National_GDB.zip"), recursive = TRUE)
  unlink(paste0(basedir, "/data/misc/WBD_National_GDB.gdb"), recursive = TRUE)
  unlink(paste0(basedir, "/data/misc/WBD_National_GDB.xml"))
  unlink(paste0(basedir, "/data/misc/WBD_National_GDB.jpg"))
}
if (is.null(webshot:::find_phantom())) {
  print("-- Unpacking FOSSFlood for first use: Installing dependencies --")
  webshot::install_phantomjs(force = FALSE)
}
if (is.null(rmarkdown::pandoc_available())| !isTRUE(rmarkdown::pandoc_available())) {
  print("-- Pandoc is needed in order to print outputs, please follow the installation prompts as they appear.")
  installr::install.pandoc(use_regex = TRUE, to_restart = FALSE)
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Define Mapping area ------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
print(paste("-- Welcome to FOSSFlood - Running FOSSFlood for", user.aoi.string))
if(user.aoi.source %in% c("zctas", "huc8")) {
  user.aoi.stringlist <- as.list(strsplit(user.aoi.string, ", ")[[1]])
  user.aoi.filepath <- gsub(", ", "_", user.aoi.string)
  if(user.aoi.source=="zctas") {
    user.aoi.filepath <- paste0('zctas_', user.aoi.filepath)
    featurefile <- sf::read_sf(paste0(basedir, "/data/misc/zctas.shp"))
    user.aoi.call <- base::subset(featurefile, (featurefile$ZCTA5CE10 %in% user.aoi.stringlist))
  } else if(user.aoi.source=="hucA8") {
    user.aoi.filepath <- paste0('huc8_', user.aoi.filepath)
    featurefile <- sf::read_sf(paste0(basedir, "/data/misc/huc8.shp"))
    user.aoi.call <- base::subset(featurefile, (featurefile$HUC8 %in% user.aoi.stringlist))
  }
} else if(user.aoi.source == "string") {
  user.aoi.call <- user.aoi.string
  user.aoi.stringlist <- as.list(strsplit(user.aoi.string, ", ")[[1]])
  user.aoi.filepath <- gsub(", ", "_", user.aoi.stringlist)
}

if(length(user.aoi.call$geometry)==0){
  print("-- Requested AOI not found, double check that your features exist and that you entered them in correctly. --")
  stop()
}

# Run AOI and Floodmapping packages
AOI = AOI::aoi_get(user.aoi.call)
tryCatch({quiet(geosphere::areaPolygon(sf::as_Spatial(AOI$geometry))>0)}, 
         error = function(e) {
           print("-- Malformed input: The features were not located in the AOI.  Double check that your features exist and that you entered them in correctly.")
           stop(e)
         })

if(file.exists(paste0(basedir,"/AOI///",user.aoi.filepath,"/rating_",user.aoi.filepath,".fst"))) {
  xx = list(hand.path = paste0(basedir,"/AOI///",user.aoi.filepath,"/hand_",user.aoi.filepath,".tif"),
            catch.path = paste0(basedir,"/AOI///",user.aoi.filepath,"/catchmask_",user.aoi.filepath,".tif"),
            rating.path = paste0(basedir,"/AOI///",user.aoi.filepath,"/rating_",user.aoi.filepath,".fst"))
} else {
  xx = quiet(FloodMapping::getRawData(AOI,paste0(basedir,"/AOI/"),user.aoi.filepath))
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Generate aux basedata ----------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
if (!file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_rec.shp"))) {
  # Setup ----------------------------------------------------------------------------
  memory.size(max = TRUE)
  dir.create(paste0(basedir,"/AOI/", user.aoi.filepath,"/tmp"))
  dir.create(paste0(basedir,"/AOI/", user.aoi.filepath,"/output"))
  dir.create(paste0(basedir,"/AOI/", user.aoi.filepath,"/output/server"))
  
  # Write hardclip file
  if(user.aoi.source %in% c("zctas", "huc8")) {
    sf::write_sf(sf::st_as_sf(user.aoi.call), paste0(basedir,"/AOI/", user.aoi.filepath,"/",user.aoi.filepath,"_hardclip.shp"), delete_layer = TRUE, quiet = TRUE)
  }
  
  # grab bounding box for api calls
  xx$aoi.bb <- sf::read_sf(paste0(basedir,"/AOI/", user.aoi.filepath,"/",user.aoi.filepath,".shp")) %>% 
    sf::st_transform(sf::st_crs(4326))
  south <- sf::st_bbox(xx$aoi.bb)[2]
  north <- sf::st_bbox(xx$aoi.bb)[4]
  west <- sf::st_bbox(xx$aoi.bb)[1]
  east <- sf::st_bbox(xx$aoi.bb)[3]
  # ---------------------------------------------------------------------------------
  
  # NHD Flowlines ----------------------------------------------------------------------------
  print("-- Downloading NHD Flowlines --")
  URL <- paste0("https://labs.waterdata.usgs.gov/geoserver/wmadata/ows?service=WFS&version=2.0.0&request=GetFeature&typeNames=wmadata:nhdflowline_network&srsName=EPSG:4326&bbox=",
                south,",",west,",",north,",",east,
                "&outputFormat=SHAPE-ZIP")
  tryCatch({
   quiet(httr::GET(URL, write_disk(paste0(basedir, "/AOI/",user.aoi.filepath,"/tmp/Flowlines.zip"), overwrite=TRUE)))
  },
  error=function(cond) {
    unlink(paste0(basedir,"/AOI/", user.aoi.filepath), recursive=TRUE)
    print("-- Waterdata.usgs.gov threw an error, make sure you're connected to the internet and run the AOI again --")
    print("-- Here's the original error message:")
    message(cond)
  })
  
  unzip(paste0(basedir, "/AOI/",user.aoi.filepath,"/tmp/Flowlines.zip"), exdir = paste0(basedir, "/AOI/",user.aoi.filepath))
  
  # needed to plot correctly (pull and drop feature geom, project and rewrite, clean up)
  tmplines = sf::read_sf(paste0(basedir,"/AOI/",user.aoi.filepath,"/nhdflowline_network.shp")) 
  goodlines <- suppressWarnings(sf::st_zm(tmplines, drop = T, what = "ZM") %>% 
                                  sf::st_set_crs(4326))
  var.out.bool <- names(goodlines) %in% c("comid", "gnis_name", "ftype","streamorde", "streamleve", "streamcalc")
  suppressWarnings(sf::write_sf(goodlines[,var.out.bool], paste0(basedir,"/AOI/",user.aoi.filepath,"/nhdflowline_network.shp"), delete_layer = TRUE, quiet = TRUE))
  quiet(file.remove(paste0(basedir, "/AOI/",user.aoi.filepath,"/wfsrequest.txt"))) 
  # ---------------------------------------------------------------------------------
  
  # gages ----------------------------------------------------------------------------
  print("-- Downloading NWIS stations --")
  skipNWISflag <- FALSE
  URL <- paste0("https://waterdata.usgs.gov/nwis/inventory?nw_longitude_va=",
                west,"&nw_latitude_va=",north,"&se_longitude_va=",east,"&se_latitude_va=",south,
                "&coordinate_format=decimal_degrees&group_key=NONE&format=sitefile_output&sitefile_output_format=rdb_file&column_name=agency_cd&column_name=site_no&column_name=station_nm&column_name=dec_lat_va&column_name=dec_long_va&list_of_search_criteria=lat_long_bounding_box")
  tryCatch(quiet(httr::GET(URL, write_disk(paste0(basedir, "/AOI/",user.aoi.filepath,"/inventory"), overwrite=TRUE))),
           error = function(e) {
             print("-- !ALERT! No NWIS stations found in AOI. --")
             skipNWISflag <<- TRUE
           }
  )
  if(!skipNWISflag) {
    NWISgagestxt <- read.table(paste0(basedir, "/AOI/",user.aoi.filepath,"/inventory"), sep="\t", header=TRUE)[-1,]  # Remove colume data def
    NWISgagestxt$Latitude <- as.numeric(as.character(NWISgagestxt$dec_lat_va))
    NWISgagestxt$Longitude <- as.numeric(as.character(NWISgagestxt$dec_long_va))
    sp::coordinates(NWISgagestxt) <- ~Longitude + Latitude
    NWISgages <- sf::st_as_sf(subset(NWISgagestxt, select=-c(dec_long_va,dec_lat_va))) %>% sf::st_set_crs(4326)
    NWISgages_proj = sf::st_transform(NWISgages, sf::st_crs(3857))
    # Shoehorn in comid
    comidsource = raster::raster(paste0(basedir,"/AOI/", user.aoi.filepath,"/catchmask_",user.aoi.filepath,".tif"))
    comidv = velox::velox(comidsource)
    NWISgageCOMID = comidv$extract_points(NWISgages_proj)
    colnames(NWISgageCOMID) <- "comid"
    NWISgages <- do.call(cbind, list(NWISgages, NWISgageCOMID))
    suppressWarnings(sf::write_sf(NWISgages, paste0(basedir,"/AOI/",user.aoi.filepath,"/NWISgages.shp"), delete_layer = TRUE, quiet = TRUE))
    quiet(file.remove(paste0(basedir, "/AOI/",user.aoi.filepath,"/inventory"))) 
  }
  # ---------------------------------------------------------------------------------
  
  # TIGER ----------------------------------------------------------------------------
  print("-- Downloading TIGER roads --")
  CountiesFIPS <- quiet(tigris::counties(cb=TRUE))
  CountiesFIPS_Proj <- sf::st_transform(sf::st_as_sf(CountiesFIPS), sf::st_crs(4326))
  aoiCounties <- suppressMessages(CountiesFIPS_Proj[xx$aoi.bb, ])
  # Error checking: Download roads for AOI that interset more than one state or county
  if( nrow(aoiCounties) > 1) {
    tigerroads_tmp <- quiet(tigris::roads(aoiCounties[1,]$STATEFP,aoiCounties[1,]$COUNTYFP, year = 2018, refresh = TRUE))
    for (i in 1:nrow(aoiCounties)-1) {
      j = i+1
      tigerroads_tmp1 <- quiet(tigris::roads(aoiCounties[j,]$STATEFP,aoiCounties[j,]$COUNTYFP, year = 2018, refresh = TRUE))
      tigerroads_tmp <- rbind(tigerroads_tmp, tigerroads_tmp1)
    }
    tigerroads <- tigerroads_tmp
  } else {
    tigerroads <- quiet(tigris::roads(unique(aoiCounties$STATEFP),unique(aoiCounties$COUNTYFP), year = 2018, refresh = TRUE))
  }
  tigerroads_Proj <- sf::st_transform(sf::st_as_sf(tigerroads), sf::st_crs(4326))
  tigerroads_Proj_Sub <- suppressMessages(tigerroads_Proj[xx$aoi.bb, ])
  var.out.bool <- names(tigerroads_Proj_Sub) %in% c("FULLNAME", "RTTYP")
  tigerout <- tigerroads_Proj_Sub[,var.out.bool] %>%
    dplyr::rename('name' = FULLNAME) %>%
    dplyr::rename('type' = RTTYP)
  sf::write_sf(tigerout, paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_tiger.shp"), delete_layer = TRUE, quiet = TRUE)
  quiet(gc())
  # ---------------------------------------------------------------------------------
  for(t in 1:length(unique(aoiCounties$STATEFP))){
    print(paste0("-- Downloading OSM Dataset for ",cdlTools::fips( unique(aoiCounties$STATEFP)[t], to = "Abbreviation")," --"))
    tryCatch({
       quiet(httr::GET(paste0("http://download.geofabrik.de/north-america/us/",tolower(gsub(" ", "-", cdlTools::fips( unique(aoiCounties$STATEFP)[t], to = "Name"))),"-latest.osm.pbf"), 
                    write_disk(paste0(basedir,"/AOI/",user.aoi.filepath,"/tmp/",tolower(gsub(" ", "-", cdlTools::fips(unique(aoiCounties$STATEFP)[t], to = "Name"))),"-latest.osm.pbf"), 
                               overwrite=TRUE)))
    },
      error=function(cond) {
        unlink(paste0(basedir,"/AOI/", user.aoi.filepath), recursive=TRUE)
        print("-- OSM servers threw an error, make sure you're connected to the internet and run the AOI again --")
        print("-- Here's the original error message:")
        message(cond)
    })
    
    print(paste0("-- Building OSM Roads for ",cdlTools::fips(unique(aoiCounties$STATEFP)[t], to = "Abbreviation")," --"))
    OSMlines <- geofabrik::read_pbf(paste0(basedir,"/AOI/",user.aoi.filepath,"/tmp/",
                                           tolower(gsub(" ", "-", cdlTools::fips(unique(aoiCounties$STATEFP)[t], to = "Name"))),
                                           "-latest.osm.pbf"),
                                    layer = "lines")
    OSMlines_Proj <- sf::st_transform(sf::st_as_sf(OSMlines), sf::st_crs(4326))
    OSMlines_Proj_Sub <- suppressMessages(OSMlines_Proj[xx$aoi.bb, ])
    var.out.bool <- names(OSMlines_Proj_Sub) %in% c("name", "highway")
    OSMlinesout <- OSMlines_Proj_Sub[,var.out.bool] %>%
      dplyr::rename('name' = name) %>%
      dplyr::rename('type' = highway) %>%
      mutate(type=recode(type, 
                         'track'=as.character(NA),
                         'residential'="M",
                         'service'="O",
                         'path'=as.character(NA),
                         'unclassified'=as.character(NA),
                         'footway'=as.character(NA),
                         'primary'="M",
                         'secondary'="C",
                         'tertiary'="O")) %>%
      drop_na(type) %>%
      sf::st_as_sf()
    # ---------------------------------------------------------------------------------
    
    # OSM addresses ---------------------------------------------------------------------------------
    print(paste0("-- Building OSM Addresses for ",cdlTools::fips(unique(aoiCounties$STATEFP)[t], to = "Abbreviation")," --"))
    OSMpoints <- geofabrik::read_pbf(paste0(basedir,"/AOI/",user.aoi.filepath,"/tmp/",
                                            tolower(gsub(" ", "-", cdlTools::fips(unique(aoiCounties$STATEFP)[t], to = "Name"))),
                                            "-latest.osm.pbf"),
                                     layer = "points")
    OSMpoints_Proj <- sf::st_transform(sf::st_as_sf(OSMpoints), sf::st_crs(4326))
    OSMpoints_Proj_Sub <- suppressMessages(OSMpoints_Proj[xx$aoi.bb, ])
    # layer_options = "ENCODING=UTF-8", delete_layer = TRUE, quiet = TRUE)
    var.out.bool <- names(OSMpoints_Proj_Sub) %in% c("name", "address")
    OSMpointsout <- OSMpoints_Proj_Sub[,var.out.bool]
    OSMpointsout$address <- ifelse(is.na(OSMpointsout$address), OSMpointsout$name, OSMpointsout$address)
    OSMpointsout$address <- ifelse(is.na(OSMpointsout$address), "NA No CAMA Data Avail", OSMpointsout$address)
    OSMpointsout$dataset="OpenStreetMaps"
    OSMpointsout <- OSMpointsout[ , !(names(OSMpointsout) %in% c("name"))]
    
    # Open Addresses ---------------------------------------------------------------------------------
    quiet(gc())
    print(paste0("-- Downloading OpenAddresses for ",cdlTools::fips(unique(aoiCounties$STATEFP)[t], to = "Abbreviation")," --"))
    tryCatch({
      validOA <- openadds::oa_list() %>% dplyr::filter(str_detect(processed, paste0("/us/",tolower(cdlTools::fips(unique(aoiCounties$STATEFP)[t], to = "Abbreviation")))))
      urls = validOA$processed 
      out = list()
      for( l in 1:length(validOA$processed)) {
        out[[l]] <- tryCatch({
          openadds::oa_get(validOA$processed[[l]])},
          error   = function(e){NULL})
        message(l)
      }
      quiet(gc())
      c = parse(text = paste0("mergedStateData <- openadds::oa_combine(",str_c("out[[", c(1:(length(out))), "]]", sep = "",  collapse = ", "),")"))
      mergedStateData <- eval(c)
      quiet(gc())
      mergedStateData <- na.omit(mergedStateData)
      print(paste0("-- Georeferencing OpenAddressess - may be slow depending on the size of the database --"))
      sp::coordinates(mergedStateData) <- ~lon+lat
      sfoadata <- sf::st_as_sf(mergedStateData) %>% sf::st_set_crs(4326)
      print(paste0("-- Trimming OpenAddressess - may be slow depending on the size of the database --"))
      aoiPointBounds <- suppressWarnings(suppressMessages(sfoadata[xx$aoi.bb, ]))
      print(paste0("-- Cleaning OpenAddressess - may be slow depending on the size of the database --"))
      aoiPoints <- suppressWarnings(suppressMessages(aoiPointBounds[!duplicated(aoiPointBounds),]))
    },
    error = function(e){
      quiet(gc())
      print(paste0("-- !ALERT! OpenAddressess failed, generating empty dataset --"))
      lon <<- c(0)
      lat <<- c(0)
      address <<- c("test")
      dataset <<- c("test")
      mergedStateData <<- data.frame(lon, lat, address, dataset)
      sp::coordinates(mergedStateData) <<- ~lon+lat
      sfoadata <<- sf::st_as_sf(mergedStateData) %>% sf::st_set_crs(4326)
      aoiPoints <<- suppressWarnings(suppressMessages(sfoadata[sfoadata$dataset=="empty", ]))
    })
    quiet(gc())
    
    if(t==1) {
      OSMlinesout_hold <- OSMlinesout
      OSMpointsout_hold <- OSMpointsout
      aoiPoints_hold <- aoiPoints
    }
    if(t>=2) {
      OSMlinesout_hold <- rbind(OSMlinesout_hold, OSMlinesout)
      OSMpointsout_hold <- rbind(OSMpointsout_hold, OSMpointsout)
      aoiPoints_hold <- rbind(aoiPoints_hold, aoiPoints)
    }
  }
  
  sf::write_sf(OSMlinesout_hold, paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_osm.shp"), delete_layer = TRUE, quiet = TRUE)
  sf::write_sf(OSMpointsout_hold, paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_osm.shp"), delete_layer = TRUE, quiet = TRUE)
  sf::write_sf(aoiPoints_hold, paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_oa.shp"), delete_layer = TRUE, quiet = TRUE)
  quiet(gc())
  
  # -- Make ancillary products ---------------------------------------------------------------------------------
  print("-- Generating Index Fishnets --")
  suppressWarnings(suppressMessages(hex_grid <- make_grid(xx$aoi.bb, type = "hexagonal", cell_area = 0.00006, clip = FALSE)))
  suppressWarnings(suppressMessages(rec_grid <- make_grid(xx$aoi.bb, type = "square", cell_area = 0.00006, clip = FALSE)))
  sf::write_sf(sf::st_as_sf(hex_grid[,names(hex_grid) %in% c("IndexLabel")]), paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_hex.shp"), delete_layer = TRUE, quiet = TRUE)
  sf::write_sf(sf::st_as_sf(rec_grid[,names(rec_grid) %in% c("IndexLabel")]), paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_rec.shp"), delete_layer = TRUE, quiet = TRUE)
  # ---------------------------------------------------------------------------------
  
  unlink(paste0(basedir,"/AOI/",user.aoi.filepath,"/tmp"), recursive=TRUE)
  
  # -- quickstats ---------------------------------------------------------------------------------
  print("--Data prep finished: report--")
  print(paste0("AOI generated: ",user.aoi.string))
  print(paste0("File path: ",basedir,"/AOI/",user.aoi.filepath))
  print(paste0("Area processed: ",sf::st_area(xx$aoi.bb)," [m^2] => ",toString(sf::st_area(xx$aoi.bb)*3.86102e-7)," [sq mi]"))
  if(all(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/catchmask_",user.aoi.filepath,".tif")),
         file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/hand_",user.aoi.filepath,".tif")),
         file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/rating_",user.aoi.filepath,".fst")))) {
    print("HAND/SRC mapping data downloaded sucessfully") } else { print("AOI package failed") }
  if(all(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_osm.shp")),
         file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_osm.shp")))) {
    print(paste0("No errors generating osm data: ",nrow(OSMpoints_Proj_Sub)," unique OSM addresses and ",nrow(OSMlines_Proj_Sub),
                 " road segments equivalent to ",toString(sum(sf::st_length(OSMlines_Proj_Sub))*0.000621371)," miles of roads in AOI")) } else { print("Issues generating osm data") }
  if(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_oa.shp"))) {
    if(length(aoiPoints$geometry)==0) {
      print(paste0("Openaddress database empty (bad AOI or too large to process in one go)"))
    } else {
      print(paste0("No issues generating openaddress database: ",nrow(aoiPoints)," unique points")) }
  } else { print("Issues generating openaddresses") }
  if(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_tiger.shp"))) {
    print(paste0("No errors generating TIGER data: ",nrow(tigerroads_Proj_Sub)," road segments equivalent to ",toString(sum(sf::st_length(tigerroads_Proj_Sub))*0.000621371)," miles of roads in AOI")) } else { print("Issues generating TIGER roads") }
  if(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/NWISgages.shp"))) {
    print(paste0("No issues downloading NWIS data: ",nrow(NWISgages)," potential points in the AOI")) } else { print("Issues generating NWIS data") }
  if(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/nhdflowline_network.shp"))) {
    print(paste0("No errors generating NHD flowline data: ",nrow(goodlines)," segments equivalent to ",toString(sum(sf::st_length(goodlines))*0.000621371)," miles of waterways in AOI")) } else { print("Issues downloading flow lines") }
  if(all(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_hex.shp")),
         file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_rec.shp")))) {
    print("No errors generating grids") } else { print("Issues generating grids") }
  if(all(
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_hex.shp")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_rec.shp")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/nhdflowline_network.shp")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/NWISgages.shp")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/catchmask_",user.aoi.filepath,".tif")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/hand_",user.aoi.filepath,".tif")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/rating_",user.aoi.filepath,".fst")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_osm.shp")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_osm.shp")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_tiger.shp")),
    file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_oa.shp"))
  )) {
    paste("-- No errors generating base data - Map on! --")
  }
  
   # -- Write out report ---------------------------------------------------------------------------------
  capture.output(
    cat(
      "--Data prep finished: report--\n",
      "AOI generated:",user.aoi.string,"\n",
      "Built on:",date(),"\n",
      "Area processed:",sf::st_area(xx$aoi.bb),"[m^2] => ",toString(sf::st_area(xx$aoi.bb)*3.86102e-7),"[sq mi]\n",
      "Number of HAND Cells:", raster::ncell(raster::raster(xx$hand.path)),"\n\n",
      "Address points:\n",
      "OpenAddresses:",nrow(aoiPoints),"address points\n",
      "Open Street Maps points:",nrow(OSMpoints_Proj_Sub),"address points\n\n",
      "Roads:\n",
      "TIGER Roads:",nrow(tigerroads_Proj_Sub),"road segments equivalent to",toString(sum(sf::st_length(tigerroads_Proj_Sub))*0.000621371),"miles of roads in AOI\n",
      "Open Street Maps Roads:",nrow(OSMlines_Proj_Sub),"road segments equivalent to",toString(sum(sf::st_length(OSMlines_Proj_Sub))*0.000621371),"miles of roads\n\n",  
      "NHD flowline data:\n",
      nrow(goodlines),"segments equivalent to",toString(sum(sf::st_length(goodlines))*0.000621371),"miles of waterways in AOI"
      ),
    file = paste0(basedir,"/AOI/",user.aoi.filepath,"/AOI_Report.txt")
  )
  
  quiet(gc())
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Load in base data and paths ----------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
print("-- Welcome to FOSSFlood - Loading in data")

#/////////////////////////////////////
# Paths
#/////////////////////////////////////
xx$aoibb.path <- paste0(basedir,"/AOI/", user.aoi.filepath,"/",user.aoi.filepath,".shp")
if(user.aoi.source %in% c("zctas", "huc8")) {
  xx$aoi.path <- paste0(basedir,"/AOI/", user.aoi.filepath,"/",user.aoi.filepath,"_hardclip.shp")
} else {
  xx$aoi.path <- paste0(basedir,"/AOI/", user.aoi.filepath,"/",user.aoi.filepath,".shp")
}

#/////////////////////////////////////
# CFIM
#/////////////////////////////////////
xx$aoi.shp <- sf::read_sf(xx$aoi.path)
xx$aoi.shp_t <- sf::st_transform(xx$aoi.shp, sf::st_crs(3857))
xx$catch.path = paste0(basedir,"/AOI/", user.aoi.filepath,"/catchmask_",user.aoi.filepath,".tif")
xx$catch.grid = raster::raster(xx$catch.path)
# raster::projection(xx$catch.grid) <- raster::crs("+init=epsg:3857")
xx$hand.path = paste0(basedir,"/AOI/", user.aoi.filepath,"/hand_",user.aoi.filepath,".tif")
xx$hand.grid = raster::raster(xx$hand.path)
# raster::projection(xx$hand.grid) <- raster::crs("+init=epsg:3857")
xx$rating.path = paste0(basedir,"/AOI/", user.aoi.filepath,"/rating_",user.aoi.filepath,".fst")
xx$rating.file <- fst::read_fst(xx$rating.path)
#----------
SkipNWISFlows <<- FALSE
tryCatch({xx$gage.point = sf::read_sf(paste0(basedir,"/AOI/",user.aoi.filepath,"/NWISgages.shp"))}, 
         error = function(e) {
           SkipNWISFlows <<- TRUE
         })
#----------
xx$flow.line = sf::read_sf(paste0(basedir,"/AOI/",user.aoi.filepath,"/nhdflowline_network.shp"))
xx$flow.point = suppressWarnings(sf::st_cast(xx$flow.line,"POINT"))
xx$tigerroads.path <- paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_tiger.shp")
xx$osmroads.path <- paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_osm.shp")
xx$osmadd.path <- paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_osm.shp")
xx$oaadd.path <- paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_oa.shp")

#/////////////////////////////////////
# Addresses
#/////////////////////////////////////
if(user.address.source == "OpenAddresses") { 
  xx$address.path <- xx$oaadd.path 
} else if (user.address.source == "OpenStreetMap Addresses") { 
  xx$address.path <- xx$osmadd.path 
} else if (user.address.source == "User Provided Addresses") { 
  xx$address.path <- user.address.file 
}
tryCatch(xx$address.point <- sf::read_sf(xx$address.path),
         error = function(e) {
           print("-- !ALERT! Invalid addresses provided, defaulting to OpenAddresses database --")
           xx$address.path <<- xx$oaadd.path
           xx$address.point <<- sf::read_sf(xx$address.path)
           user.address.source <<- "OpenAddresses"
         }
)
if(length(xx$address.point$geometry)==0){
  if(user.address.source == "OpenAddresses") {
    user.address.source <- "OpenStreetMap Addresses"
    xx$address.path <- xx$osmadd.path
    xx$address.point <- quiet(sf::read_sf(xx$address.path))
  } else {
    user.address.source <- "OpenAddresses"
    xx$address.path <- xx$oaadd.path
    xx$address.point <- quiet(sf::read_sf(xx$address.path))
  }
  print(paste0("-- !ALERT! Address database is empty, swapping to ",user.address.source," --"))
}

#/////////////////////////////////////
# Roads
#/////////////////////////////////////
if(user.road.source == "OpenStreetMaps") { 
  xx$road.path <- xx$osmroads.path 
} else if (user.road.source == "TIGER Lines 2018") { 
  xx$road.path <- xx$tigerroads.path 
} else if(user.road.source == "UserRoads") { 
  xx$road.path <- user.road.file 
}  
tryCatch(xx$road.line <- sf::read_sf(xx$road.path),
         error = function(e) {
           print("-- !ALERT! Invalid addresses provided, defaulting to TIGER database --")
           xx$road.path <<- xx$tigerroads.path
           xx$road.line <<- sf::read_sf(xx$road.path)
           user.road.source <<- "TIGER Lines 2018"
         }
)
xx$road.point = suppressWarnings(sf::st_cast(xx$road.line,"POINT"))

# Code in full names
xx$road.line$type_full <- xx$road.line$type %>%
  recode("C"="C - County", "I"="I - Interstate", "M"="M - Common Name", "O"="O - Other","S"="S - State recognized","U"="U - U.S.")

#/////////////////////////////////////
# Grids
#/////////////////////////////////////
if(user.output.grid == "Square") {
  xx$mapindex.poly <- sf::read_sf(paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_rec.shp"))
} else if(user.output.grid == "Hexagon") {
  xx$mapindex.poly <- sf::read_sf(paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_hex.shp"))
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Load CFIM flow format ----------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
print("-- Grabbing flows --")
# Flows in cms
if(substr(user.forecast.source, 1, 4)=="USER") {
  user.forecast.gen <- character()
  user.forecast.gen[1] <- "User"
  user.forecast.gen[2] <- "Provided"
  fileext <- substring(user.forecast.file,nchar(user.forecast.file)-2)
  if(fileext=="fst") {
    xx$nwm.flow <- user.forecast.file
  } #else if(fileext=="csv") {
    #xx$nwm.flow <- user.forecast.file
    #xx$nwm.discharge <- read.csv(xx$nwm.flow)
  #} 0.028316846592
  
} else if(user.forecast.source=="NWM_SR_C") {
  user.forecast.gen <- getMostRecentForecast()
  xx$nwm.flow = quiet(nomadsNC::create_nomads_fst(type = "short_range", num = user.forecast.timesteps, dstfile = paste0(basedir,"/AOI/",user.aoi.filepath,"/flows.fst")))
  if(user.output.archive) { 
    write.fst(xx$nwm.flow, paste0(basedir,"/AOI/",user.aoi.filepath,"/output/",user.forecast.source,getMostRecentForecast()[1],"_",getMostRecentForecast()[2],"_flows.fst")) 
  }
} else if(user.forecast.source=="NWM_ARC") {
  nwmArchiveFlow = nwmHistoric::readNWMdata(comid = unique(raster::getValues(xx$catch.grid)), startDate = user.forecast.start, endDate = user.forecast.end) %>%
    select(-model) %>%
    tidyr::pivot_wider(id_cols = comid, names_from = time, values_from = flow) %>%
    rename(COMID = comid)
   write.fst(nwmArchiveFlow, paste0(basedir,"/AOI/",user.aoi.filepath,"/output/archive_",user.forecast.start,"_",user.forecast.end,"flows.fst")) 
   xx$nwm.flow <- paste0(basedir,"/AOI/",user.aoi.filepath,"/output/archive_",user.forecast.start,"_",user.forecast.end,"flows.fst")
}

print("-- Generating discharge --")
xx$nwm.discharge <- fst::read.fst(xx$nwm.flow)
user.forecast.timesteps <- ncol(xx$nwm.discharge) - 1
if(user.forecast.source=="NWM_SR_C") {
  xx$nwm.discharge.dateTime <- names(xx$nwm.discharge)[-1] %>% 
    stringr::str_sub(3, 21)
} else {
  xx$nwm.discharge.dateTime <- names(xx$nwm.discharge)[-1]
}
xx$nwm.discharge.dateTimeZoneLocal <- as.POSIXct(xx$nwm.discharge.dateTime,tz=Sys.timezone())
xx$nwm.discharge.dateTimeZoneZ <- xx$nwm.discharge.dateTimeZoneLocal 
attr(xx$nwm.discharge.dateTimeZoneZ, "tzone") <- "UTC"
xx$nwm.discharge.dateTimeZoneLocalHR <- format(xx$nwm.discharge.dateTimeZoneLocal, '%m/%d/%Y %I:%M %p')
# xx$nwm.discharge.dateTimeZoneZHR <- sapply(xx$nwm.discharge.dateTimeZoneZ, as.character) 
xx$nwm.discharge.dateTimeZoneZHR <- format(xx$nwm.discharge.dateTimeZoneZ, "%Y-%m-%d %H:%M:%S %Z")

if(!SkipNWISFlows) {
  firstdate <- lubridate::date(xx$nwm.discharge.dateTimeZoneZ[1])
  lastdate <-  lubridate::date(xx$nwm.discharge.dateTimeZoneZ[user.forecast.timesteps])
  tryCatch({
    xx$nwis.discharge <- readNWISuv(xx$gage.point$site_no, "00060", as.Date(firstdate)-2, as.Date(lastdate)+1, tz="UTC")
  }, 
  error = function(e) {
    SkipNWISFlows <<- TRUE
  })
  tryCatch({
    xx$nwis.stage <- readNWISuv(xx$gage.point$site_no, "00065", as.Date(firstdate)-2, as.Date(lastdate)+1, tz="UTC")
  }, 
  error = function(e) {
    SkipNWISFlows <<- TRUE
  })
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- CFIM Flood raster generation ---------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
print("-- Mapping inundation depths --")
xx$flood.grid = FloodMapping::map_flood(hand.path = xx$hand.path,
                                        catchment.path = xx$catch.path,
                                        rating.path = xx$rating.path,
                                        flows.path = xx$nwm.flow,
                                        add = 0)
xx$flood.grid <- round(xx$flood.grid, digits = 3)
xx$nwm.timestep = names(xx$flood.grid)

quiet(gc()) # to avoid Error in x$.self$finalize() : attempt to apply non-function
if(user.output.hardclip) {
  print("-- Cropping to AOI --")
  xx$aoi.shp_d <- sf::st_union(xx$aoi.shp) %>% sf::st_as_sf()
  xx$aoi.shp_td <- sf::st_union(xx$aoi.shp_t) %>% sf::st_as_sf()
  
  xx$catch.grid <- xx$catch.grid %>% raster::mask(xx$aoi.shp_td)
  xx$hand.grid <- xx$hand.grid %>% raster::mask(xx$aoi.shp_td)
  xx$flood.grid <- xx$flood.grid %>% raster::mask(xx$aoi.shp_td)
  
  xx$mapindex.poly <- suppressWarnings(suppressMessages(xx$mapindex.poly[xx$aoi.shp_d, ]))
  xx$road.line <- suppressWarnings(suppressMessages(sf::st_intersection(xx$road.line, xx$aoi.shp_d))) %>% sf::st_collection_extract("LINESTRING")
  xx$flow.line <- suppressWarnings(suppressMessages(sf::st_intersection(xx$flow.line, xx$aoi.shp_d)))# %>% sf::st_collection_extract("LINESTRING")
  xx$address.point <- suppressWarnings(suppressMessages(xx$address.point[xx$aoi.shp_d, ]))
}
# Calc and push ffreq layer
xx$flood.grid <- raster::addLayer(xx$flood.grid, sum(xx$flood.grid>=0, na.rm = TRUE))
names(xx$flood.grid)[raster::nlayers(xx$flood.grid)] <- 'ffreq'
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Generate GIS data --------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
if(user.output.choice=="GIS_O") {
  print("--Writing flood rasters --")
  exportnames <- as.character(xx$nwm.discharge.dateTimeZoneZ) %>% stringr::str_sub(6, 13) %>% paste0("Z.tif")
  exportnames[user.forecast.timesteps+1] <- "ffreq.tif"
  setwd(paste0(basedir,"/AOI/",user.aoi.filepath,"/output"))
  raster::writeRaster(xx$flood.grid, filename=exportnames, bylayer=TRUE,format="GTiff")
  print(paste0("-- Flood inundation rasters written to: ",getwd()," --"))
  setwd(basedir)
  if(user.output.choice=="GIS_O") {
    print("-- FOSSFlood execution complete --")
    stop()
  }
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Generate impacts ---------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
print("-- Intersecting addresses --")
v = velox::velox(xx$flood.grid)
xx$address.velox = v$extract_points(sf::st_transform(xx$address.point, raster::crs(xx$flood.grid)))
colnames(xx$address.velox) <- names(xx$flood.grid)
xx$address.point <- do.call(cbind, list(xx$address.point, xx$address.velox))

xx$flow.point <- suppressWarnings(suppressMessages(sf::st_intersection(xx$mapindex.poly, xx$flow.point))) %>% 
  rename('Index Label' = IndexLabel)
# xx$address.point[names(xx$flood.grid)[1]] <- NA_real_     # testing

processskipflag <- TRUE
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- No impacts or just want to see data? -------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
if(all(xx$address.point$ffreq==0,na.rm = TRUE) | user.output.choice=="basedata") {
  if(all(xx$address.point$ffreq==0,na.rm = TRUE)) { 
    print("-- !Alert!: No addresses forecasted to be impacted, defaulting to base data viewer --") 
    user.output.choice="basedata"
  }
  
  # Build basedata titleblock
  TitleBlock <- magick::image_read(paste0(basedir,"/data/misc/EmptyTitleBlock.png"))
  if(user.aoi.source=="string") {
    TitleBlock <- magick::image_annotate(TitleBlock, user.aoi.string, font = 'Georgia', location = "+398+302", size = 20)
  } else {
    TitleBlock <- magick::image_annotate(TitleBlock, paste(user.aoi.source,user.aoi.string), font = 'Georgia', location = "+398+302", size = 20)
  }
  if(substr(user.forecast.source, 1, 8)=="USER_DIS" || substr(user.forecast.source, 1, 10)=="USER_STAGE") {
    if(substr(user.forecast.source, 1, 8)=="USER_DIS") { 
      TitleBlock <- magick::image_annotate(TitleBlock, "User Provided Discharge", font = 'Georgia', location = "+305+328", size = 20)
    } else { 
      TitleBlock <- magick::image_annotate(TitleBlock, "User Provided Stage", font = 'Georgia', location = "+305+328", size = 20) 
    }
    magick::image_annotate(TitleBlock, Sys.time(), font = 'Georgia', location = "+303+355", size = 20)
  } else if(user.forecast.source=="NWM_SR_C") {
    if(user.forecast.source=="NWM_SR_C") {
      TitleBlock <- magick::image_annotate(TitleBlock, "National Water Model Short Range Forecast", font = 'Georgia', location = "+305+328", size = 20)
    }
    forecastZ <- ISOdatetime(year=substr(user.forecast.gen[1], 0, 4),
                             month=substr(user.forecast.gen[1], 5, 6),
                             day=substr(user.forecast.gen[1], 7, 8),
                             hour=user.forecast.gen[2],
                             min="00",
                             sec="00",
                             tz="GMT") 
    forecastT <- lubridate::with_tz(forecastZ, tz=Sys.timezone())
    TitleBlock <- magick::image_annotate(TitleBlock, paste(format(forecastZ,format='%Y-%m-%d %H:%M') ,"Z =",format(forecastT,format='%Y-%m-%d %H:%M') ,Sys.timezone()), font = 'Georgia', location = "+303+355", size = 20)
  }
  ifelse(user.address.source == "OpenStreetMap Addresses",addString <<- "OSM Addresses",
    ifelse(user.address.source == "User Provided Addresses",addString <<- "User Addresses",
    addString <<- user.address.source
    ))
  TitleBlock <- magick::image_annotate(TitleBlock, format(file.mtime(xx$aoi.path),format='%m/%d/%Y'), font = 'Georgia', location = "+418+383", size = 20) %>% 
    magick::image_annotate(addString, font = 'Georgia', location = "+301+410", size = 20) %>%
    magick::image_annotate(user.road.source, font = 'Georgia', location = "+627+410", size = 20) %>% 
    imager::magick2cimg(alpha = "rm") %>%
    imager::autocrop("white") %>%
    imager::cimg2magick(rotate = T) %>%
    magick::image_flop()
  
  xx$address.mapbook = suppressWarnings(suppressMessages(sf::st_intersection(xx$mapindex.poly, xx$address.point)))
  xx$road.point.mapbook = suppressWarnings(suppressMessages(sf::st_intersection(xx$mapindex.poly, xx$road.point)))
  xx$mapindex.point <- suppressWarnings(suppressMessages(xx$mapindex.poly %>% st_centroid()))
  
  # Human readable field names
  xx$mapindex.poly <- xx$mapindex.poly %>% 
    rename('Index Label' = IndexLabel)
  xx$mapindex.point <- xx$mapindex.point %>% 
    rename('Index Label' = IndexLabel)
  xx$address.mapbook <- xx$address.mapbook %>% 
    rename('Index Label' = IndexLabel)
  xx$road.point.mapbook <- xx$road.point.mapbook %>% 
    rename('Index Label' = IndexLabel)
  # Set process skip flag
  processskipflag <- FALSE
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Generate rest of flood impacts -------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
if(processskipflag) {
  print("-- Intersecting roads --")
  xx$road.depth = v$extract_points(sf::st_transform(xx$road.point, raster::crs(xx$flood.grid)))
  colnames(xx$road.depth) <- names(xx$flood.grid)
  xx$road.point =  do.call(cbind, list(xx$road.point, xx$road.depth))
  xx$address.mapbook = suppressWarnings(suppressMessages(sf::st_intersection(xx$mapindex.poly, xx$address.point)))
  # sf::st_geometry(xx$address.mapbook) = NULL  
  xx$road.point.mapbook = suppressWarnings(suppressMessages(sf::st_intersection(xx$mapindex.poly, xx$road.point)))
  # sf::st_geometry(xx$road.points.mapbook) = NULL  
  
  #/////////////////////////////////////
  # summary tables
  #/////////////////////////////////////
  print("-- Summarizing results --")
  for (i in 1:length(xx$nwm.timestep)) {
    xx$floodsum.x[[i]] <- i
    xx$floodsum.nwmts[[i]] <- xx$nwm.timestep[[i]]
    xx$floodsum.area[[i]] <- raster::cellStats(xx$flood.grid[[i]]>0, 'sum', na.rm=TRUE)*raster::res(xx$flood.grid[[i]])[1]*raster::res(xx$flood.grid[[i]])[2]*1e-6
    xx$floodsum.addimpact[[i]] <- sum(xx$address.point[[xx$nwm.timestep[[i]]]]>0, na.rm = TRUE)
    xx$floodsum.roadimpact[[i]] <- sum(xx$road.point[[xx$nwm.timestep[[i]]]]>0, na.rm = TRUE)
    # Slows the whole process down
    # xx$road.poly[[i]] <- suppressMessages(suppressWarnings(
    #   subset(xx$road.point,xx$nwm.timestep[[i]]>0) %>%
    #     sf::st_buffer(0.000098) %>% 
    #     summarise(area=sum(i))
    # ))
    xx$road.poly[[i]] <- suppressMessages(suppressWarnings(xx$road.point[which(xx$road.point[[xx$nwm.timestep[[i]]]]>0),] %>% 
                                                             sf::st_buffer(0.000098) %>% 
                                                             sf::st_union() %>% 
                                                             sf::st_as_sf()))
  }
  xx$chartdata <- data.frame(xx$floodsum.x,xx$floodsum.nwmts,xx$floodsum.area,xx$floodsum.addimpact,xx$nwm.discharge.dateTimeZoneZ)
  #add freq to swiper
  xx$floodsum.x[[length(xx$nwm.timestep)+1]] <- length(xx$nwm.timestep)+1
  xx$floodsum.nwmts[[length(xx$nwm.timestep)+1]] <- "Frequency"
  xx$floodsum.area[[length(xx$nwm.timestep)+1]] <- raster::cellStats(xx$flood.grid[[length(xx$nwm.timestep)+1]]>0, 'sum', na.rm=TRUE)*raster::res(xx$flood.grid[[length(xx$nwm.timestep)+1]])[1]*raster::res(xx$flood.grid[[length(xx$nwm.timestep)+1]])[2]*1e-6
  xx$road.poly[[length(xx$nwm.timestep)+1]] <- suppressMessages(suppressWarnings(xx$road.point[which(xx$road.point$ffreq>0),] %>% 
                                                                                   sf::st_buffer(0.000098) %>% 
                                                                                   sf::st_union() %>% 
                                                                                   sf::st_as_sf()))
  xx$floodsum.addimpact[[length(xx$nwm.timestep)+1]] <- sum(xx$address.point$ffreq>0, na.rm = TRUE)
  
  print("-- Generating cartography objects --")
  xx$mapindex.impacts <- xx$address.mapbook %>% 
    group_by(IndexLabel) %>% 
    summarise(wethousehours = sum(ffreq),
              impacthousecount = sum(ffreq>0),
              maximumdepth = suppressWarnings(max(eval(parse(text=xx$nwm.timestep)), na.rm=T))) %>%
    mutate(maximumdepth = ifelse(is.infinite(maximumdepth), NA, maximumdepth))
  
  xx$mapindex.poly <- xx$mapindex.poly %>% left_join(as.data.frame(xx$mapindex.impacts)[c("IndexLabel","wethousehours","impacthousecount","maximumdepth")],by="IndexLabel")
  xx$mapindex.poly$maximumdepth <- round(xx$mapindex.poly$maximumdepth, digits = 2)
  xx$mapindex.point <- suppressWarnings(suppressMessages(xx$mapindex.poly %>% st_centroid()))
  # Subset here for impacts as well, since subsetting in app is undesireable
  # xx$impactindex.poly <- subset(xx$mapindex.poly, impacthousecount>0)
  # xx$impactindex.point <- suppressWarnings(suppressMessages(xx$impactindex.poly %>% st_centroid()))
  
  # Human readable field names
  xx$mapindex.poly <- xx$mapindex.poly %>% 
    rename('Wet House Hours' = wethousehours) %>%
    rename('Uniquely Impacted Addresses' = impacthousecount) %>%
    rename('Maximum Impact Depth' = maximumdepth) %>%
    rename('Index Label' = IndexLabel)
  xx$mapindex.point <- xx$mapindex.point %>% 
    rename('Wet House Hours' = wethousehours) %>%
    rename('Uniquely Impacted Addresses' = impacthousecount) %>%
    rename('Maximum Impact Depth' = maximumdepth) %>%
    rename('Index Label' = IndexLabel)
  xx$address.mapbook <- xx$address.mapbook %>% 
    rename('Index Label' = IndexLabel)
  xx$road.point.mapbook <- xx$road.point.mapbook %>% 
    rename('Index Label' = IndexLabel)
  
  # Chart generation
  # xx$chart.addts<- timeSeries::timeSeries(xx$floodsum.addimpact[1:user.forecast.timesteps], xx$nwm.discharge.dateTimeZoneZ, units = "Count of Impacted Addresses")
  # xx$chart.areats <-timeSeries::timeSeries(xx$floodsum.area[1:user.forecast.timesteps], xx$nwm.discharge.dateTimeZoneZ, units = "Square miles of Inundated Land")
  # xx$chart.summary = cbind(xx$chart.addts, xx$chart.areats) 
  
  xx$chart.addts <- xts::xts(x=xx$floodsum.addimpact[1:user.forecast.timesteps],
                             order.by= as.POSIXct(xx$nwm.discharge.dateTimeZoneZ, tz = "UTC"))
  xx$chart.areats <- xts::xts(x=xx$floodsum.area[1:user.forecast.timesteps],
                              order.by= as.POSIXct(xx$nwm.discharge.dateTimeZoneZ, tz = "UTC"))
  xx$chart.summary = cbind(xx$chart.addts, xx$chart.areats) 
  names(xx$chart.summary) <- c("Count of Impacted Addresses","Square miles of Inundated Land")
  
  # Create titleblock
  TitleBlock <- magick::image_read(paste0(basedir,"/data/misc/EmptyTitleBlock.png"))
  if(user.aoi.source=="string") {
    TitleBlock <- magick::image_annotate(TitleBlock, user.aoi.string, font = 'Georgia', location = "+398+302", size = 20)
  } else {
    TitleBlock <- magick::image_annotate(TitleBlock, paste(user.aoi.source,user.aoi.string), font = 'Georgia', location = "+398+302", size = 20)
  }
  if(substr(user.forecast.source, 1, 8)=="USER_DIS" || substr(user.forecast.source, 1, 10)=="USER_STAGE") {
    if(substr(user.forecast.source, 1, 8)=="USER_DIS") { 
      TitleBlock <- magick::image_annotate(TitleBlock, "User Provided Discharge", font = 'Georgia', location = "+305+328", size = 20)
    } else { 
      TitleBlock <- magick::image_annotate(TitleBlock, "User Provided Stage", font = 'Georgia', location = "+305+328", size = 20) 
    }
    magick::image_annotate(TitleBlock, Sys.time(), font = 'Georgia', location = "+303+355", size = 20)
  } else if(user.forecast.source=="NWM_SR_C") {
    if(user.forecast.source=="NWM_SR_C") {
      TitleBlock <- magick::image_annotate(TitleBlock, "National Water Model Short Range Forecast", font = 'Georgia', location = "+305+328", size = 20)
    }
    forecastZ <- ISOdatetime(year=substr(user.forecast.gen[1], 0, 4),
                             month=substr(user.forecast.gen[1], 5, 6),
                             day=substr(user.forecast.gen[1], 7, 8),
                             hour=user.forecast.gen[2],
                             min="00",
                             sec="00",
                             tz="GMT") 
    forecastT <- lubridate::with_tz(forecastZ, tz=Sys.timezone())
    TitleBlock <- magick::image_annotate(TitleBlock, paste(format(forecastZ,format='%Y-%m-%d %H:%M') ,"Z =",format(forecastT,format='%Y-%m-%d %H:%M') ,Sys.timezone()), font = 'Georgia', location = "+303+355", size = 20)
  }
  ifelse(user.address.source == "OpenStreetMap Addresses",addString <<- "OSM Addresses",
    ifelse(user.address.source == "User Provided Addresses",addString <<- "User Addresses",
    addString <<- user.address.source
    ))
  TitleBlock <- magick::image_annotate(TitleBlock, format(file.mtime(xx$aoi.path),format='%Y-%m-%d'), font = 'Georgia', location = "+418+383", size = 20) %>% 
    magick::image_annotate(addString, font = 'Georgia', location = "+301+410", size = 20) %>%
    magick::image_annotate(user.road.source, font = 'Georgia', location = "+627+410", size = 20) %>% 
    imager::magick2cimg(alpha = "rm") %>%
    imager::autocrop("white") %>%
    imager::cimg2magick(rotate = T) %>%
    magick::image_flop() 
  
  magick::image_write(TitleBlock, path = paste0(basedir,"/AOI/",user.aoi.filepath,"/output/titleblock.png"), format = "png")
  
  # impact specific vis
  xx$mapindex.poly.WHHval <- as.numeric(c(0:max(xx$mapindex.poly$'Wet House Hours', na.rm = TRUE)))
  xx$mapindex.poly.WHHpal <- colorBin("Reds", xx$mapindex.poly$'Wet House Hours', 4, pretty = TRUE)
  xx$mapindex.poly.IHCval <- as.numeric(c(0:max(xx$mapindex.poly$'Uniquely Impacted Addresses', na.rm = TRUE)))
  xx$mapindex.poly.IHCpal <- colorBin("Reds", xx$mapindex.poly$'Uniquely Impacted Addresses', 4, pretty = TRUE)
  xx$mapindex.poly.MDval <- as.numeric(c(0:max(xx$mapindex.poly$'Maximum Impact Depth', na.rm = TRUE)))
  xx$mapindex.poly.MDpal <- colorBin("Reds", xx$mapindex.poly$'Maximum Impact Depth', 4, pretty = TRUE)
  
} # End of impact processing
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Loading remaining cart variables -----------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# Index view lock
xx$aoi.bb.view <- suppressMessages(suppressWarnings(sf::st_buffer(
  sf::st_combine(xx$aoi.shp),
  0.08,
  nQuadSegs = 1,
  endCapStyle = "SQUARE",
  joinStyle = "ROUND",
  mitreLimit = 1
)))

# vis params
xx$hand.grid.val <- as.numeric(c(0:raster::maxValue(xx$hand.grid)))
xx$hand.grid.pal <- colorBin("Blues", xx$hand.grid.val, 10,na.color = "transparent")
xx$catch.grid.val <- raster::unique(xx$catch.grid)
xx$catch.grid.pal <- colorFactor(palette = randomcoloR::distinctColorPalette(length(raster::unique(xx$catch.grid))),
                                 domain = xx$catch.grid.val)
xx$flood.grid.val <- as.numeric(c(0:10))
xx$flood.grid.pal <- colorNumeric(c("#deebf7", "#9ecae1", "#3182bd"), xx$flood.grid.val, na.color = "transparent")
xx$flood.grid.ffreqval <- as.numeric(c(1:length(xx$nwm.timestep)))
xx$flood.grid.ffreqpal <- colorNumeric(c("#deebf7", "#9ecae1", "#3182bd"), xx$flood.grid.ffreqval, na.color = "transparent")

xx$vis_roadpal <- colorFactor(palette = c('yellow', 'red', 'black', 'gray', 'black', 'red'), 
                              domain = c("C","I","M","O","S","U"))
xx$vis_roadpalFull <- colorFactor(palette = c('yellow', 'red', 'black', 'gray', 'black', 'red'), 
                                  domain = c("C - County","I - Interstate","M - Common Name","O - Other","S - State recognized","U - U.S."))
xx$flow.line$popup_text <- 
  paste0('<strong>Stream name:', xx$flow.line$gnis_name, '</strong>','<br/>', 
         'Stream order: ', xx$flow.line$streamorde) %>% 
  lapply(htmltools::HTML)
xx$road.line$popup_text <-
  paste0('<strong>Road name:', xx$road.line$name, '</strong>') %>% 
  lapply(htmltools::HTML)
xx$address.point$popup_text <-
  paste0('<strong>Road name:', xx$address.point$address, '</strong>') %>% 
  lapply(htmltools::HTML)

xx$gage.points.standard <- leaflet::awesomeIcons(icon = 'home',lib = 'fa',markerColor = 'transparent',iconColor = "#000000")
xx$address.point.standard <- leaflet::makeAwesomeIcon(icon = 'home',lib = 'fa',iconColor = "black",markerColor = 'transparent')
xx$address.point.yellow <- leaflet::awesomeIcons(icon = 'home',lib = 'fa',markerColor = 'transparent',iconColor = "Yellow")
xx$address.point.orange <- leaflet::awesomeIcons(icon = 'home',lib = 'fa',markerColor = 'transparent',iconColor = "Orange")
xx$address.point.red <- leaflet::awesomeIcons(icon = 'home',lib = 'fa',markerColor = 'transparent',iconColor = "Red")

providers$none=NULL
basemapchoices <- list(
  "Minimal map with lables" = c('Hydda'=providers$Hydda.Full,'TonerLite'=providers$Stamen.TonerLite,'Toner'=providers$Stamen.Toner),
  "Minimal maps" = c('Empty'=providers$none,'Hydda'=providers$Hydda.Base,'TonerBackground'=providers$Stamen.TonerBackground),
  "Just lables" = c('Black Labels'=providers$CartoDB.DarkMatterOnlyLabels,'CartoDB Labels'=providers$CartoDB.VoyagerOnlyLabels,'Terrain Labels'=providers$Stamen.TerrainLabels),
  "Streetmaps" = c('OpenStreetMap'=providers$OpenStreetMap.Mapnik,'ESRI'=providers$Esri.WorldStreetMap,'NatGeo'=providers$Esri.NatGeoWorldMap,'CartoDB Voyager'=providers$CartoDB.Voyager),
  "Natural" = c('ESRI Imagery'=providers$Esri.WorldImagery,'ESRI World Topo'=providers$Esri.WorldTopoMap,'ESRI World Terrain'=providers$Esri.WorldTerrain,'World Shaded Relief'=providers$Esri.WorldShadedRelief),
  "Fun" = c('Watercolor'=providers$Stamen.Watercolor,'Night Mode'=providers$CartoDB.DarkMatter))

DISCLAIMERString =
  "SOURCES: TIGER 2018, OpenStreetMaps, OpenAddresses, and the Free and Open Source Community <br> DISCLAIMER OF WARRANTY: The University disclaims all warranties, including all implied warranties of accuracy, merchantability, and fitness for a particular purpose. User expressly waives any claim that it may have 
against the University, its employees or agents and shall not be liable for any direct, indirect, consequential, special, punitive or other damages arising from or related to the use of the FOSSFlood application and 
the decisions made based on its outputs."

dyCSScoolPrint <- function(dygraph){
  dygraph$x$css <- '
  .dygraph-axis-label {font-size: 25px;}
  .dygraph-title {
  font-size: 25px;
  }
  .dygraph-legend {
  width: auto !important;
  min-width: 150px;
  color: white;
  background-color: #BABABA !important;
  padding-left:5px;
  border-color:#BABABA;
  border-style:solid;
  border-width:thin;
  transition:0s 4s;
  z-index: 80 !important;
  box-shadow: 2px 2px 5px rgba(0, 0, 0, .3);
  border-radius: 3px;
  }
  .dygraph-legend:hover{
  transform: translate(-110%);
  transition: 0s;
  }
  .dygraph-legend > span {
  color: black;
  padding-left:5px;
  padding-right:2px;
  margin-left:-5px;
  background-color: white !important;
  display: block;
  }
  .dygraph-legend > span:first-child {
  margin-top:2px;
  }
  .dygraph-legend > span > span{
  display: inline;
  }
  .highlight {
  border-left: 2px solid #BABABA;
  padding-left:3px !important;
  }
  '
  dygraph
}
dyCSScoolWeb <- function(dygraph){
  dygraph$x$css <- '
  .dygraph-title {
  font-size: 35px;
  padding-bottom:7px !important;
  }
  .dygraph-legend {
  width: auto !important;
  min-width: 150px;
  color: white;
  background-color: #BABABA !important;
  padding-left:5px;
  border-color:#BABABA;
  border-style:solid;
  border-width:thin;
  transition:0s 4s;
  z-index: 80 !important;
  box-shadow: 2px 2px 5px rgba(0, 0, 0, .3);
  border-radius: 3px;
  }
  .dygraph-legend:hover{
  transform: translate(-110%);
  transition: 0s;
  }
  .dygraph-legend > span {
  color: black;
  padding-left:5px;
  padding-right:2px;
  margin-left:-5px;
  background-color: white !important;
  display: block;
  }
  .dygraph-legend > span:first-child {
  margin-top:2px;
  }
  .dygraph-legend > span > span{
  display: inline;
  }
  .highlight {
  border-left: 2px solid #BABABA;
  padding-left:3px !important;
  }
  '
  dygraph
}
dyCSScoolMin <- function(dygraph){
  dygraph$x$css <- '
  .dygraph-title {
  font-size: 20px;
  }
  '
  dygraph
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////




#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Create shiny server files ------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
if(user.output.serverFiles) {
  print("-- Generating shiny server files --")
  unlink(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/*"),recursive=TRUE,force=TRUE)
  magick::image_write(TitleBlock, path = paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/titleblock.png"), format = "png")
  raster::writeRaster(xx$flood.grid,paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/floodstack.grd"),overwrite=TRUE, format="raster")
  saveRDS(xx$aoi.shp, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/aoi.RData"))
  saveRDS(xx$road.poly, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/road_poly.RData"))
  saveRDS(xx$road.line, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/road_line.RData"))
  saveRDS(xx$address.point, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/address.RData"))
  saveRDS(xx$mapindex.poly, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/index_poly.RData")) 
  saveRDS(xx$mapindex.point, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/index_point.RData"))
  saveRDS(xx$address.mapbook, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/address_maps.RData"))
  saveRDS(xx$road.point.mapbook, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/road_maps.RData"))
  saveRDS(xx$flow.line, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/flow_line.RData"))
  saveRDS(xx$chart.summary, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chartdata.RData"))
  saveRDS(xx$floodsum.addimpact, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chartimpact.RData"))
  saveRDS(xx$floodsum.area, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chartarea.RData"))
  saveRDS(xx$nwm.discharge.dateTimeZoneZ, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chart_DTZ.RData"))
  saveRDS(xx$nwm.discharge.dateTimeZoneZHR, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chart_DTZHR.RData"))
  saveRDS(xx$nwm.discharge.dateTimeZoneLocalHR, file=paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chart_DTLHR.RData"))
  file.copy(list.files(paste0(basedir,"/data/misc/liveserver/")), paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/"))
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Base data viewer ---------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
basedataUI <- dashboardPage(
  skin = "blue",
  # dashboardHeader(title="FOSSFlood V 1.2",tags$li(class="dropdown",actionButton("print", "Print"))),
  dashboardHeader(title=apptitle),
  dashboardSidebar(
    sidebarMenu(id="sidebar",
                h3("Map controls"),
                pickerInput(inputId='userbasemap',label='Select a base map:',choices=basemapchoices,selected=basemapchoices[3]$`Just lables`[1]),
                shinyWidgets::switchInput(inputId="gridColor",label="Grid color",onLabel="Black",offLabel="White",value=TRUE),
                shinyWidgets::switchInput(inputId="aoiColor",label="AOI shading",onLabel="Black",offLabel="White",value=TRUE),
                menuItem("basedata", tabName = "basedata", icon = icon("info-circle")),
                menuItem("mapbook", tabName = "mapbook", icon = icon("arrows-h")),
                add_busy_spinner(spin="orbit", color = "red",position='bottom-left', margins = c(40, 25))
    )
  ),
  dashboardBody(
    useShinyalert(),
    useShinyjs(),
    tabItems(
      tabItem(tabName="basedata",value="basedata",
              tags$head(
                tags$style(type = "text/css", "#basemap {height: calc(60vh - 80px) !important;}"),
                tags$style(type = "text/css", "#mapbookmap {height: calc(60vh - 20px) !important;}"),
                tags$style(type="text/css", "#titleblock img {max-width: 100%; width: 100%; height: auto;}")
              ),
              fluidRow(
                column(width = 8,
                       box(width = NULL, solidHeader = FALSE, status = "warning",leafletOutput("basemap")),
                       tabBox(title = tagList(shiny::icon("gear"), "Summary Tables"), width = NULL,
                              tabPanel(title = tagList(shiny::icon("home"), "Addresses"),
                                       DT::dataTableOutput("addressDataset")),
                              tabPanel(title = tagList(shiny::icon("road"), "Roads"),
                                       DT::dataTableOutput("roadsDataset")),
                              tabPanel(title = tagList(shiny::icon("water"), "Streams"),
                                       DT::dataTableOutput("streamsDataset"))
                       )
                ),
                column(width = 4,
                       box(width = NULL,height = 287,imageOutput("titleblock")),
                       box(title = "Map controls", width = NULL, status = "primary","Map controls"),
                       box(title = "Warnings", width = NULL, status = "primary",HTML(DISCLAIMERString))
                )
              )
      ),
      tabItem(tabName = "mapbook",value="mapbook",
              tags$head(
                tags$style(type="text/css", "#swipermaintitleblock img {max-width: 100%; width: 100%; height: auto;}"),
                tags$style(type="text/css", "#titleblock1 img {max-width: 100%; width: 100%; height: auto;}")
              ),
              fluidRow(
                column(width = 8,
                       box(title = "map book map", width = NULL, status = "primary",leafletOutput("mapbookmap")),
                       box(title = "features", width = NULL, status = "primary",tabBox(title = tagList(shiny::icon("gear"), "Summary Tables"), width = NULL,
                                                                                       tabPanel(title = tagList(shiny::icon("home"), "Addresses"),
                                                                                                DT::dataTableOutput("filteredaddressesmapbook")),
                                                                                       tabPanel(title = tagList(shiny::icon("road"), "Roads"),
                                                                                                DT::dataTableOutput("filteredroadsmapbook")),
                                                                                       tabPanel(title = tagList(shiny::icon("water"), "Streams"),
                                                                                                DT::dataTableOutput("filteredstreamsmapbook")))
                       )
                ),
                column(width = 4,
                       imageOutput("titleblock1"),
                       box(title = "Map index", width = NULL, status = "primary",
                           leafletOutput("indexmap"),
                           fluidRow(width = NULL,align="center",
                                    actionBttn(inputId="backbutton", label = "Back-Previous", style = "pill", color = "danger"),
                                    actionBttn(inputId="forwardbutton", label = "Forward-Next", style = "pill", color = "danger"),
                                    textOutput("activeIndexLabel"))
                       ),
                       box(title = "Warnings", width = NULL, status = "primary",HTML(DISCLAIMERString))
                )
              )
      )
    )
  )
)
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
basedataSERVER <- function(input, output, session) {
  # TitleBlocks ========================================
  output$titleblock <- renderImage({
    tmpfile <- TitleBlock %>%
      magick::image_write(tempfile(fileext='jpg'), format = 'jpg')
    list(
      src = tmpfile,
      filetype = "image/jpg",
      deleteFile = TRUE
    )
  })
  output$titleblock1 <- renderImage({
    tmpfile <- TitleBlock %>%
      magick::image_write(tempfile(fileext='jpg'), format = 'jpg')
    list(
      src = tmpfile,
      filetype = "image/jpg",
      deleteFile = TRUE
    )
  })
  output$activeIndexLabel <- renderText({
    paste0("Viewing:",values$currentIndexLabel)
  })
  
  # InforBoxes ========================================
  values <- reactiveValues()
  values$count <- 1
  values$currentIndexLabel <- "None"
  values$currentIndexAdds <- NULL
  values$currentIndexRoads <- NULL
  values$currentIndexStreams <- NULL
  
  # Maps ========================================
  basemap_reactive <- reactive({
    gridcolor <- "Black"
    aoiColor <- "Black"
    sumvisfill <- T
    
    leaflet() %>%
      addProviderTiles("CartoDB.DarkMatterOnlyLabels", group = "Base Map") %>%
      addMapPane("bg", zIndex = 410) %>%
      # addMapPane("flood", zIndex = 412) %>%
      addMapPane("hand", zIndex = 411) %>%
      addMapPane("catch", zIndex = 412) %>%
      addMapPane("flow", zIndex = 414) %>%
      addMapPane("gage", zIndex = 416) %>%
      addMapPane("roadimpact", zIndex = 418) %>%
      addMapPane("road", zIndex = 420) %>%
      addMapPane("add", zIndex = 422) %>%
      addMapPane("index", zIndex = 424) %>%
      addMapPane("indexlabels", zIndex = 426) %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
      addPolylines(data=xx$flow.line,weight = ~streamorde/2,color = "Blue",options=pathOptions(pane="flow"),group="Flowlines") %>%
      addRasterImage(x=xx$hand.grid,colors=xx$hand.grid.pal,layerId="hand",group="hand",project=FALSE,maxBytes=20*1024*1024) %>%
      addRasterImage(x=xx$catch.grid,colors=xx$catch.grid.pal,maxBytes=20*1024*1024,layerId="catch",group="catch") %>%
      addRasterImage(x=xx$flood.grid$ffreq,colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
      addPolylines(data=xx$road.poly[[length(xx$nwm.timestep)+1]], color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpalFull(type_full),weight = 0.8,options=pathOptions(pane="road"),group="Roads") %>%
      addAwesomeMarkers(data=xx$address.point,
                        icon=xx$address.point.standard, 
                        options=pathOptions(pane="add"),
                        clusterOptions=markerClusterOptions(),
                        group="Addresses") %>%
      addPolygons(data=xx$mapindex.poly, fill=F,color=gridcolor,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addLabelOnlyMarkers(data=xx$mapindex.point, label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = T,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","hand","catch","Flooding","Road Impacts","Addresses","Borders","Index","Index Labels"),
        options = layersControlOptions(collapsed = FALSE)
      ) %>% 
      addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
      addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
      addLegend("bottomleft", title = "Flood Frequency",pal = xx$flood.grid.ffreqpal, values = xx$flood.grid.ffreqval, group = "Flooding") %>%
      hideGroup("Addresses") %>%
      hideGroup("Road Impacts") %>%
      hideGroup("catch") %>%
      hideGroup("hand") %>%
      hideGroup("Index Labels")
  })
  indexmap_reactive <- reactive({
    gridcolor <- "Black"
    aoiColor <- "Black"
    sumvisfill <- T
    
    leaflet() %>%
      addProviderTiles("Stamen.TonerBackground", group = "Base Map") %>%
      addMapPane("bg", zIndex = 410) %>%
      # addMapPane("flood", zIndex = 412) %>%
      addMapPane("flow", zIndex = 414) %>%
      addMapPane("gage", zIndex = 416) %>%
      addMapPane("roadimpact", zIndex = 418) %>%
      addMapPane("road", zIndex = 420) %>%
      addMapPane("add", zIndex = 422) %>%
      addMapPane("index", zIndex = 424) %>%
      addMapPane("indexlabels", zIndex = 426) %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders")# %>%
    # addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
    # addPolygons(data=xx$mapindex.poly, fill=F,color=gridcolor,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
    # addLabelOnlyMarkers(data=xx$mapindex.point, label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = T,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
    # addLayersControl(
    #   overlayGroups = c("Flowlines","Roads","Borders","Index","Index Labels"),
    #   options = layersControlOptions(collapsed = FALSE)
    # ) %>%
    # hideGroup("Index Labels")
  })
  mapbookmap_reactive <- reactive({
    gridcolor <- "Black"
    aoiColor <- "Black"
    sumvisfill <- T
    
    leaflet() %>%
      addProviderTiles("CartoDB.DarkMatterOnlyLabels", group = "Base Map") %>%
      addMapPane("bg", zIndex = 410) %>%
      # addMapPane("flood", zIndex = 412) %>%
      addMapPane("hand", zIndex = 411) %>%
      addMapPane("catch", zIndex = 412) %>%
      addMapPane("flow", zIndex = 414) %>%
      addMapPane("gage", zIndex = 416) %>%
      addMapPane("roadimpact", zIndex = 418) %>%
      addMapPane("road", zIndex = 420) %>%
      addMapPane("add", zIndex = 422) %>%
      addMapPane("index", zIndex = 424) %>%
      addMapPane("indexlabels", zIndex = 426) %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
      # addRasterImage(x=xx$flood.grid[[1]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
      addPolylines(data=xx$flow.line,weight = ~streamorde/2,color = "Blue",options=pathOptions(pane="flow"),group="Flowlines") %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpalFull(type_full),weight = 0.8,options=pathOptions(pane="road"),group="Roads") %>%
      addAwesomeMarkers(data=xx$address.point,
                        icon=xx$address.point.standard,
                        options=pathOptions(pane="add"),
                        clusterOptions=markerClusterOptions(),
                        group="Addresses") %>%
      addPolygons(data=xx$mapindex.poly, fill=F,color=gridcolor,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addLabelOnlyMarkers(data=xx$mapindex.point, label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = T,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Addresses","Borders","Index","Index Labels"),
        options = layersControlOptions(collapsed = FALSE)
      ) %>%
      addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
      addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow")
  })
  
  output$basemap <- renderLeaflet({
    basemap_reactive()
  })
  output$indexmap <- renderLeaflet({
    indexmap_reactive()
  })
  output$mapbookmap <- renderLeaflet({
    mapbookmap_reactive()
  })
  
  observeEvent(input$basemap_zoom, {
    if(as.numeric(input$basemap_zoom) > 11) {
      leafletProxy("basemap") %>%
        showGroup("Index Labels")
    } else {
      leafletProxy("basemap") %>%
        hideGroup("Index Labels")
    }
  })
  observeEvent(input$indexmap_zoom, {
    if(as.numeric(input$indexmap_zoom) > 12) {
      leafletProxy("indexmap") %>%
        showGroup("Index Labels")
    } else {
      leafletProxy("indexmap") %>%
        hideGroup("Index Labels")
    }
  })
  observeEvent(input$mapbookmap_zoom, {
    if(as.numeric(input$mapbookmap_zoom) > 12) {
      leafletProxy("mapbookmap") %>%
        showGroup("Index Labels")
    } else {
      leafletProxy("mapbookmap") %>%
        hideGroup("Index Labels")
    }
  })
  
  observeEvent((input$gridColor | input$aoiColor), {
    shiny::req(input$timeslider)
    shiny::req(input$uislidertime)
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    
    leafletProxy("basemap") %>%
      clearGroup("Borders") %>%
      clearGroup("Index") %>%
      addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
      addPolygons(data=xx$mapindex.poly, fill=F,color=gridcolor,weight = 1,options=pathOptions(pane="index"),group="Index")
    
    leafletProxy("mapbookmap") %>%
      clearGroup("Borders") %>%
      clearGroup("Index") %>%
      addPolygons(data=xx$mapindex.poly, fill=F,color=gridcolor,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders")
    
  })
  observeEvent(input$userbasemap, {
    
    leafletProxy("basemap") %>%
      clearGroup("Base Map") %>%
      addProviderTiles(input$userbasemap, group = "Base Map")
    
    leafletProxy("mapbookmap") %>%
      clearGroup("Base Map") %>%
      addProviderTiles(input$userbasemap, group = "Base Map")
    
  })
  observeEvent(input$backbutton, {
    shiny::req(input$gridColor)
    shiny::req(input$aoiColor)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(values$count != 1) {
      values$count <- values$count - 1
    } else {
      values$count <- length(xx$mapindex.poly)
    }
    
    values$currentIndexLabel <- xx$mapindex.point[values$count,]$`Index Label`
    values$currentIndexAdds <- subset(xx$address.mapbook, xx$address.mapbook$`Index Label`==values$currentIndexLabel) %>% select(address,`Index Label`,dataset)
    values$currentIndexRoads <- subset(xx$road.point.mapbook, xx$road.point.mapbook$`Index Label`==values$currentIndexLabel) %>% select(`Index Label`,name)
    values$currentIndexStreams <- subset(xx$flow.point, xx$flow.point$`Index Label`==values$currentIndexLabel) %>% select(`Index Label`,comid,gnis_name,streamorde)
    
    mybbox <- AOI::bbox_coords(sf::st_buffer(
      x=xx$mapindex.poly[values$count,],
      dist=0.0002,
      nQuadSegs = 1,
      endCapStyle = "SQUARE",
      joinStyle = "ROUND",
      mitreLimit = 1,
      singleSide = FALSE
    ))
    
    leafletProxy("indexmap") %>%
      clearGroup("Index") %>%
      addPolygons(data=xx$mapindex.poly, fill=F,color=gridcolor,fillOpacity=1,weight = 2,options=pathOptions(pane="index"),group="Index") %>%
      addPolygons(data=xx$mapindex.poly[values$count,], fill=F,color="Red",weight = 3,options=pathOptions(pane="index"),group="Index")
    
    leafletProxy("mapbookmap") %>%
      flyToBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
    
  })
  observeEvent(input$forwardbutton, {
    shiny::req(input$gridColor)
    shiny::req(input$aoiColor)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(values$count != nrow(xx$mapindex.poly)) {
      values$count <- values$count + 1
    } else {
      values$count <- 1
    }
    
    values$currentIndexLabel <- xx$mapindex.point[values$count,]$`Index Label`
    values$currentIndexAdds <- subset(xx$address.mapbook, xx$address.mapbook$`Index Label`==values$currentIndexLabel) %>% select(address,`Index Label`,dataset)
    values$currentIndexRoads <- subset(xx$road.point.mapbook, xx$road.point.mapbook$`Index Label`==values$currentIndexLabel) %>% select(`Index Label`,name)
    values$currentIndexStreams <- subset(xx$flow.point, xx$flow.point$`Index Label`==values$currentIndexLabel) %>% select(`Index Label`,comid,gnis_name,streamorde)
    
    mybbox <- AOI::bbox_coords(sf::st_buffer(
      x=xx$mapindex.poly[values$count,],
      dist=0.0002,
      nQuadSegs = 1,
      endCapStyle = "SQUARE",
      joinStyle = "ROUND",
      mitreLimit = 1,
      singleSide = FALSE
    ))
    
    leafletProxy("indexmap") %>%
      clearGroup("Index") %>%
      addPolygons(data=xx$mapindex.poly, fill=F,color=gridcolor,fillOpacity=1,weight = 2,options=pathOptions(pane="index"),group="Index") %>%
      addPolygons(data=xx$mapindex.poly[values$count,], fill=F,color="Red",weight = 3,options=pathOptions(pane="index"),group="Index")
    
    leafletProxy("mapbookmap") %>%
      flyToBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
    
  })
  
  observeEvent({input$indexmap_click}, {
    coords <- input$indexmap_click
    # clicked <- sf::st_as_sf(sp::SpatialPoints(matrix(c(coords$lng, coords$lat),nrow = 1)),crs = 4326)
    clicked <- sf::st_sfc(sf::st_point(c(coords$lng, coords$lat)), crs = 4326)
    # clicked <- sp::SpatialPoints(matrix(c(coords$lng, coords$lat),nrow = 1))
    # sp::proj4string(clicked) <- sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    compare <- sf::st_intersects(
      # rgeos::gBuffer(clicked, width = 0.00005), 
      clicked,
      xx$mapindex.poly)[[1]]
    # isolate(values$count <- as.numeric(match(compare$`Index Label`,xx$mapindex.poly$`Index Label`) - 1 ))
    isolate(values$count <- as.numeric(compare - 1 ))
    # shinyalert(title = paste0("vale:",coords$lng, coords$lat), type = "success")
    click("forwardbutton")
  })
  
  # Tables ========================================
  output$addressDataset <- DT::renderDataTable({
    DT::datatable(as.data.frame(xx$address.mapbook) %>% select(address,`Index Label`,dataset),
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$roadsDataset <- DT::renderDataTable({
    DT::datatable(as.data.frame(xx$road.point.mapbook) %>% select(`Index Label`,name),
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$streamsDataset <- DT::renderDataTable({
    DT::datatable(as.data.frame(xx$flow.line) %>% select(comid,gnis_name,streamorde),
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$filteredaddressesmapbook <- DT::renderDataTable({
    DT::datatable(as.data.frame(values$currentIndexAdds),
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$filteredroadsmapbook <- DT::renderDataTable({
    DT::datatable(as.data.frame(values$currentIndexRoads),
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$filteredstreamsmapbook <- DT::renderDataTable({
    DT::datatable(as.data.frame(values$currentIndexStreams),
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  
  # Events ========================================
  observeEvent(input$print, {
    # Take a screenshot of the map
    mapshot(user_created_map(), file=paste0(getwd(), '/exported_map.png'))
    shinyalert(title = "Print complete, see /shiny!", type = "success")
  })
  
  session$onSessionEnded(stopApp)
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Impact viewer ------------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
impactUI <- dashboardPage(
  skin = "blue",
  dashboardHeader(title=apptitle,tags$li(class="dropdown",actionButton("save","Save raster data"),actionButton("print", "Print"))),
  # dashboardHeader(title=apptitle,tags$li(class="dropdown",actionButton("relaunch","Relaunch UI"),actionButton("save","Save raster data"),actionButton("print", "Print"))),
  dashboardSidebar(
    sidebarMenu(id="sidebar",
                h3("Map controls"),
                pickerInput(inputId='userbasemap',label='Select a base map:',choices=basemapchoices,selected=basemapchoices[3]$`Just lables`[1]),
                shinyWidgets::switchInput(inputId="gridColor",label="Grid color",onLabel="Black",offLabel="White",value=TRUE),
                shinyWidgets::switchInput(inputId="aoiColor",label="AOI shading",onLabel="Black",offLabel="White",value=TRUE),
                menuItem("summarytab", tabName = "Summary", icon = icon("info-circle")),
                menuItem("swipertab", icon = icon("arrows-h"), tabName = "Swiper"),
                # add_busy_gif(src = "https://jeroen.github.io/images/banana.gif",height = 70, width = 70)
                # use_busy_spinner(spin="hollow-dots",choices = epic_spinners())
                add_busy_spinner(spin="orbit", color = "red",position='bottom-left', margins = c(40, 25))
    )
  ),
  dashboardBody(
    useShinyalert(),
    useShinyjs(),
    tabItems(
      tabItem(tabName = "Summary",value="Summary",
              tags$head(
                tags$style(type = "text/css", "#summary_map {height: calc(60vh - 80px) !important;}"),
                tags$style(type="text/css", "#summarytitleblock img {max-width: 100%; width: 100%; height: auto;}")
              ),
              fluidRow(
                column(width = 4,
                       box(width = NULL,height = 287,imageOutput("summarytitleblock")),
                       tabBox(title = tagList(shiny::icon("gear"), "Summary Tables"), width = NULL,
                              tabPanel(title = tagList(shiny::icon("th"), "Map index"),
                                       DT::dataTableOutput("mapindextable")),
                              tabPanel(title = tagList(shiny::icon("home"), "Addresses impacted"),
                                       DT::dataTableOutput("addressimpacts")),
                              tabPanel(title = tagList(shiny::icon("road"), "Roads impacted"),
                                       DT::dataTableOutput("roadimpacts"))
                       ),
                       radioGroupButtons(inputId="summaryviewlayer",
                                         label="View impacts as:",
                                         direction = "vertical",
                                         choices=c("Wet House Hours","Count of Impacted Addresses","Maximum Impact Depth","Individual Points"),individual = TRUE,
                                         checkIcon = list(yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),no = tags$i(class = "fa fa-circle-o", style = "color: steelblue"))),
                       box(title = "Credits & Warnings", width = NULL, status = "primary",HTML(DISCLAIMERString))),
                column(width = 8,
                       box(title = tagList(shiny::icon("chart-line"), "Time series"), width = NULL,dygraphOutput("summarychart",height = 287)),
                       box(width = NULL, solidHeader = FALSE,status = "warning",leafletOutput("summary_map")))
              )
      ),
      tabItem(tabName = "Swiper",value="Swiper",
              tags$head(
                tags$style(type = "text/css", "#summary_map {height: calc(60vh - 80px) !important;}"),
                tags$style(type="text/css", "#swipermaintitleblock img {max-width: 100%; width: 100%; height: auto;}")
              ),
              fluidRow(
                column(width = 8,
                       box(title = "time series", width = NULL,status = "primary", 
                           dygraphOutput("dygraphswiperchart"),
                           fluidRow(
                             box(width = 10, sliderTextInput(inputId = "timeslider", grid = TRUE,force_edges = TRUE,label = "Timestep to view:",choices = xx$nwm.discharge.dateTimeZoneZHR,selected = xx$nwm.discharge.dateTimeZoneZHR[1])),
                             box(width = 2, radioGroupButtons(inputId = "uislidertime",label = "Slider time labels:",choices = c("Zulu", "Local")))
                           )
                       ),
                       box(title = "SwiperMap", width = NULL, status = "primary",leafletOutput("swiper_map")),
                       box(title = "Impacted features", width = NULL, status = "primary",tabBox(title = tagList(shiny::icon("gear"), "Summary Tables"), width = NULL,
                                                                                                tabPanel(title = tagList(shiny::icon("home"), "Addresses impacted"),
                                                                                                         DT::dataTableOutput("filteraddressimpacts")),
                                                                                                tabPanel(title = tagList(shiny::icon("road"), "Roads impacted"),
                                                                                                         DT::dataTableOutput("filterroadimpacts"))
                       ))
                ),
                column(width = 4,
                       imageOutput("swipermaintitleblock"),
                       box(title = "Map index", width = NULL, status = "primary",
                           leafletOutput("index_map"),
                           fluidRow(width = NULL,align="center",
                                    actionBttn(inputId = "backbutton", label = "Back-Previous", style = "pill", color = "danger"),
                                    actionBttn(inputId = "forwardbutton", label = "Forward-Next", style = "pill", color = "danger"),
                                    textOutput("activeIndexLabel"))
                       ),
                       box(title = "Warnings", width = NULL, status = "primary",HTML(DISCLAIMERString))
                )
              )
      )
    )
  )
)
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
impactSERVER <- function(input, output, session) {
  # TitleBlocks ========================================
  output$summarytitleblock <- renderImage({
    tmpfile <- TitleBlock %>%
      magick::image_write(tempfile(fileext='jpg'), format = 'jpg')
    list(
      src = tmpfile,
      filetype = "image/jpg",
      deleteFile = TRUE
    )
  })
  output$swipermaintitleblock <- renderImage({
    tmpfile <- TitleBlock %>%
      magick::image_write(tempfile(fileext='jpg'), format = 'jpg')
    list(
      src = tmpfile,
      filetype = "image/jpg",
      deleteFile = TRUE
    )
  })
  output$activeIndexLabel <- renderText({
    paste0("Viewing:",values$currentIndexLabel)
  })
  
  # InforBoxes ========================================
  values <- reactiveValues()
  values$count <- 1
  values$currentIndexLabel <- "None"
  values$currentIndexAdds <- NULL
  values$currentIndexRoads <- NULL
  
  # Maps ========================================
  summary_map_reactive <- reactive({
    gridcolor <- "Black"
    aoiColor <- "Black"
    sumvisfill <- T
    
    leaflet() %>%
      addProviderTiles("CartoDB.DarkMatterOnlyLabels", group = "Base Map") %>%
      addMapPane("bg", zIndex = 410) %>%
      # addMapPane("flood", zIndex = 412) %>%
      addMapPane("flow", zIndex = 414) %>%
      addMapPane("gage", zIndex = 416) %>%
      addMapPane("roadimpact", zIndex = 418) %>%
      addMapPane("road", zIndex = 420) %>%
      addMapPane("add", zIndex = 422) %>%
      addMapPane("index", zIndex = 424) %>%
      addMapPane("indexlabels", zIndex = 426) %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
      addPolylines(data=xx$flow.line,weight = ~streamorde/2,color = "Blue",options=pathOptions(pane="flow"),group="Flowlines") %>%
      addRasterImage(x=xx$flood.grid$ffreq,colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
      # addMarkers(gages)
      addPolylines(data=xx$road.poly[[length(xx$nwm.timestep)+1]], color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpalFull(type_full),weight = 0.8,options=pathOptions(pane="road"),group="Roads") %>%
      addAwesomeMarkers(data=subset(xx$address.point, ffreq>0),
                        icon=xx$address.point.standard, 
                        options=pathOptions(pane="add"),
                        clusterOptions=markerClusterOptions(),
                        group="Addresses") %>%
      addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addLabelOnlyMarkers(data=subset(xx$mapindex.point, `Wet House Hours`>0), label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = T,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Road Impacts","Addresses","Borders","Index","Index Labels"),
        options = layersControlOptions(collapsed = FALSE)
      ) %>% 
      addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
      addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
      addLegend("bottomleft", title = "Flood Frequency",pal = xx$flood.grid.ffreqpal, values = xx$flood.grid.ffreqval, group = "Flooding") %>%
      addLegend("bottomleft", title = "Wet House Hours", pal = xx$mapindex.poly.WHHpal, values = subset(xx$mapindex.poly, `Wet House Hours`>0)$`Wet House Hours`, group = "Index",layerId = "Index_legend") %>%
      hideGroup("Addresses") %>%
      hideGroup("Road Impacts")
  })
  swiper_map_reactive <- reactive({
    gridcolor <- "Black"
    aoiColor <- "Black"
    sumvisfill <- T
    
    leaflet() %>%
      addProviderTiles("CartoDB.DarkMatterOnlyLabels", group = "Base Map") %>%
      addMapPane("bg", zIndex = 410) %>%
      # addMapPane("flood", zIndex = 412) %>%
      addMapPane("flow", zIndex = 414) %>%
      addMapPane("gage", zIndex = 416) %>%
      addMapPane("roadimpact", zIndex = 418) %>%
      addMapPane("road", zIndex = 420) %>%
      addMapPane("add", zIndex = 422) %>%
      addMapPane("index", zIndex = 424) %>%
      addMapPane("indexlabels", zIndex = 426) %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
      addRasterImage(x=xx$flood.grid[[1]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
      addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
      # addMarkers(gages)
      # addPolylines(data=xx$road.poly[[1]], color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),options=pathOptions(pane="road"),group="Roads") %>%
      # addAwesomeMarkers(data=subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),
      #                   icon=xx$address.point.standard,
      #                   options=pathOptions(pane="add"),
      #                   clusterOptions=markerClusterOptions(),
      #                   group="Addresses") %>%
      # addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),], fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      # addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`)), label = subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`))$`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") #%>%
      # addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`)), labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") #%>%
      addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
      addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
      # addLegend("bottomleft", title = "Flood Depth",pal = xx$flood.grid.pal, values = xx$flood.grid.val, group = "Flooding") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Road Impacts","Addresses","Borders","Index","Index Labels"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  index_map_reactive <- reactive({
    # May need to test for empty first timestep?    "Current index","Count of Impacted Addresses","Maximum Impact Depth","Individual Points"
    gridcolor <- "Black"
    aoiColor <- "Black"
    sumvisfill <- T
    
    leaflet() %>%
      addProviderTiles("Stamen.TonerBackground", group = "Base Map") %>%
      addMapPane("bg", zIndex = 410) %>%
      # addMapPane("flood", zIndex = 412) %>%
      addMapPane("flow", zIndex = 414) %>%
      addMapPane("gage", zIndex = 416) %>%
      addMapPane("roadimpact", zIndex = 418) %>%
      addMapPane("road", zIndex = 420) %>%
      addMapPane("add", zIndex = 422) %>%
      addMapPane("index", zIndex = 424) %>%
      addMapPane("indexlabels", zIndex = 426) %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
      addRasterImage(x=xx$flood.grid$ffreq,colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
      addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
      # addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      # addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),][1,],fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      # addLabelOnlyMarkers(data=subset(xx$mapindex.point,`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`),label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Borders","Index"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  # Map objects ========================================
  output$summary_map <- renderLeaflet({
    summary_map_reactive()
  })
  output$swiper_map <- renderLeaflet({
    swiper_map_reactive()
  })
  output$index_map <- renderLeaflet({
    index_map_reactive()
  })
  
  # Maps for mapshot ========================================
  mapshot_summary_map <- function() {
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    mybbox <- suppressWarnings(suppressMessages(AOI::bbox_coords(sf::st_buffer(x=xx$aoi.shp,dist=0.00005,nQuadSegs=1,endCapStyle="SQUARE",joinStyle="ROUND",mitreLimit=1,singleSide=FALSE))))
    
    m = leaflet(height=1000, width=1000) %>%
      addProviderTiles(input$userbasemap, group = "Base Map") %>%
      addMapPane("bg", zIndex = 410) %>%
      # addMapPane("flood", zIndex = 412) %>%
      addMapPane("flow", zIndex = 414) %>%
      addMapPane("gage", zIndex = 416) %>%
      addMapPane("roadimpact", zIndex = 418) %>%
      addMapPane("road", zIndex = 420) %>%
      addMapPane("add", zIndex = 422) %>%
      addMapPane("index", zIndex = 424) %>%
      addMapPane("indexlabels", zIndex = 426) %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
      addPolylines(data=xx$flow.line,weight = ~streamorde/2,color = "Blue",options=pathOptions(pane="flow"),group="Flowlines") %>%
      addRasterImage(x=xx$flood.grid$ffreq,colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpalFull(type_full),weight = 0.8,options=pathOptions(pane="road"),group="Roads") %>%
      addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addLabelOnlyMarkers(data=subset(xx$mapindex.point, `Wet House Hours`>0), label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = T,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
      addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
      addLegend("bottomleft", title = "Flood Frequency",pal = xx$flood.grid.ffreqpal, values = xx$flood.grid.ffreqval, group = "Flooding") %>%
      addLegend("bottomleft", title = "Wet House Hours", pal = xx$mapindex.poly.WHHpal, values = subset(xx$mapindex.poly, `Wet House Hours`>0)$`Wet House Hours`, group = "Index",layerId = "Index_legend") %>%
      fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
    
    return(m)
  }
  mapshot_index_map <- function(timestep,index,bboxobj){
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(is.null(index)) {
      m = leaflet(height=1000, width=1000) %>%
        addProviderTiles("Stamen.TonerBackground", group = "Base Map") %>%
        addMapPane("bg", zIndex = 410) %>%
        # addMapPane("flood", zIndex = 412) %>%
        addMapPane("flow", zIndex = 414) %>%
        addMapPane("gage", zIndex = 416) %>%
        addMapPane("roadimpact", zIndex = 418) %>%
        addMapPane("road", zIndex = 420) %>%
        addMapPane("add", zIndex = 422) %>%
        addMapPane("index", zIndex = 424) %>%
        addMapPane("indexlabels", zIndex = 426) %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
        addRasterImage(x=xx$flood.grid[[timestep]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        fitBounds(bboxobj$xmax,bboxobj$ymin,bboxobj$xmin,bboxobj$ymax)
      
    } else {
      
      m = leaflet(height=1000, width=1000) %>%
        addProviderTiles("Stamen.TonerBackground", group = "Base Map") %>%
        addMapPane("bg", zIndex = 410) %>%
        # addMapPane("flood", zIndex = 412) %>%
        addMapPane("flow", zIndex = 414) %>%
        addMapPane("gage", zIndex = 416) %>%
        addMapPane("roadimpact", zIndex = 418) %>%
        addMapPane("road", zIndex = 420) %>%
        addMapPane("add", zIndex = 422) %>%
        addMapPane("index", zIndex = 424) %>%
        addMapPane("indexlabels", zIndex = 426) %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
        addRasterImage(x=xx$flood.grid[[timestep]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),][xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`==index,],
                    fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`)),
                            label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
        fitBounds(bboxobj$xmax,bboxobj$ymin,bboxobj$xmin,bboxobj$ymax)
    }
    
    return(m)
  }
  mapshot_swiper_map <- function(timestep,index){
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(is.null(index)) {
      mybbox <- AOI::bbox_coords(sf::st_buffer(
        x=xx$aoi.shp,
        dist=0.0002,
        nQuadSegs = 1,
        endCapStyle = "SQUARE",
        joinStyle = "ROUND",
        mitreLimit = 1,
        singleSide = FALSE
      ))
      m = leaflet(height=1000, width=1000) %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addMapPane("bg", zIndex = 410) %>%
        # addMapPane("flood", zIndex = 412) %>%
        addMapPane("flow", zIndex = 414) %>%
        addMapPane("gage", zIndex = 416) %>%
        addMapPane("roadimpact", zIndex = 418) %>%
        addMapPane("road", zIndex = 420) %>%
        addMapPane("add", zIndex = 422) %>%
        addMapPane("index", zIndex = 424) %>%
        addMapPane("indexlabels", zIndex = 426) %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
        addRasterImage(x=xx$flood.grid[[timestep]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),options=pathOptions(pane="road"),group="Roads") %>%
        addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
        addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
        addLegend("bottomleft", title = "Flood Depth",pal = xx$flood.grid.pal, values = xx$flood.grid.val, group = "Flooding") %>%
        fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
      
    } else {
      mybbox <- AOI::bbox_coords(sf::st_buffer(
        x=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),][xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`==index,],
        dist=0.0002,
        nQuadSegs = 1,
        endCapStyle = "SQUARE",
        joinStyle = "ROUND",
        mitreLimit = 1,
        singleSide = FALSE
      ))
      m = leaflet(height=1000, width=1000) %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addMapPane("bg", zIndex = 410) %>%
        # addMapPane("flood", zIndex = 412) %>%
        addMapPane("flow", zIndex = 414) %>%
        addMapPane("gage", zIndex = 416) %>%
        addMapPane("roadimpact", zIndex = 418) %>%
        addMapPane("road", zIndex = 420) %>%
        addMapPane("add", zIndex = 422) %>%
        addMapPane("index", zIndex = 424) %>%
        addMapPane("indexlabels", zIndex = 426) %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=1,options=pathOptions(pane="bg"),group="Borders") %>%
        addRasterImage(x=xx$flood.grid[[timestep]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),options=pathOptions(pane="road"),group="Roads") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        # addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),][xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`==index,],
        # fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addAwesomeMarkers(data=subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),
                          icon=xx$address.point.standard, 
                          options=pathOptions(pane="add"),
                          # clusterOptions=markerClusterOptions(),
                          group="Addresses") %>%
        addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`)),
                            label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
        addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
        addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
        addLegend("bottomleft", title = "Flood Depth",pal = xx$flood.grid.pal, values = xx$flood.grid.val, group = "Flooding") %>%
        fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
    }
    
    return(m)
  }
  
  # Tables ========================================
  output$summarychart <- renderDygraph({
    dygraphs::dygraph(data = xx$chart.summary,
                      main = paste0("Forecast impacts for ", user.aoi.string),
                      xlab = "Date/Time in UTC") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = FALSE) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      # dyLegend(width = 650) %>%
      dyCSScoolWeb()
  })
  swiperchartdata <- reactive({
    req(input$uislidertime)
    ifelse(input$uislidertime=="Zulu", timeseriesOrder <- as.POSIXct(xx$nwm.discharge.dateTimeZoneZ, tz = "GMT"),timeseriesOrder <- as.POSIXct(xx$nwm.discharge.dateTimeZoneLocal, tz = Sys.timezone()))
    xx$chartXTS.addts <- xts::xts(x=xx$floodsum.addimpact[1:user.forecast.timesteps], order.by= timeseriesOrder)
    xx$chartXTS.areats <- xts::xts(x=xx$floodsum.area[1:user.forecast.timesteps], order.by= timeseriesOrder)
    xx$chartXTS.summary = cbind(xx$chartXTS.addts, xx$chartXTS.areats) 
    names(xx$chartXTS.summary) <- c("Count of Impacted Addresses","Square miles of Inundated Land")
    return(xx$chartXTS.summary)
  })
  output$dygraphswiperchart <- renderDygraph({
    req(input$uislidertime)
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    dygraphs::dygraph(data = swiperchartdata(),
                      main = paste0("Forecast impacts for ", user.aoi.string),
                      xlab = "Date/Time") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = FALSE) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      dyShading(from = lubridate::force_tz(xx$nwm.discharge.dateTimeZoneZ[match(input$timeslider, timestepseries)]-180,tzone = "UTC"), 
                to = lubridate::force_tz(xx$nwm.discharge.dateTimeZoneZ[match(input$timeslider, timestepseries)]+180,tzone = "UTC"), axis = "x", color = "#FFEAEA") %>%
      dyShading(from = Sys.time()-180, to = Sys.time()+180, axis = "x", color = "#CCEBD6") %>%
      # dyLegend(width = 600)
      dyCSScoolWeb()
  })
  
  # webshot for charts
  summararychart_image <- function(workingdir,filename) {
    setwd(workingdir)
    chart <- dygraphs::dygraph(data = xx$chart.summary,
                               main = paste0("Forecast impacts for ", user.aoi.string),
                               xlab = "Date/Time in UTC") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = FALSE, labelWidth = 120) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      dyLegend(width = 650) %>%
      dyCSScoolMin()
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1100
    height <- 400
    webshot::webshot("temp.html", file = filename,cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height)
    unlink(paste0(workingdir,"/temp_files"), recursive=TRUE)
    unlink(paste0(workingdir,"/temp.html"))
    setwd(basedir)
  }
  swiperchart_image <- function(printTimestep) {
    ifelse(input$uislidertime=="Zulu", timeseriesOrder <- as.POSIXct(xx$nwm.discharge.dateTimeZoneZ, tz = "GMT"),timeseriesOrder <- as.POSIXct(xx$nwm.discharge.dateTimeZoneLocal, tz = Sys.timezone()))
    xx$chartXTS.addts <- xts::xts(x=xx$floodsum.addimpact[1:user.forecast.timesteps], order.by= timeseriesOrder)
    xx$chartXTS.areats <- xts::xts(x=xx$floodsum.area[1:user.forecast.timesteps], order.by= timeseriesOrder)
    xx$chartXTS.summary = cbind(xx$chartXTS.addts, xx$chartXTS.areats) 
    names(xx$chartXTS.summary) <- c("Count of Impacted Addresses","Square miles of Inundated Land")
    
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    chart <- dygraphs::dygraph(data = xx$chartXTS.summary,
                               main = paste0("Forecast impacts for ", user.aoi.string),
                               xlab = "Date/Time") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = FALSE) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      dyShading(from = lubridate::force_tz(xx$nwm.discharge.dateTimeZoneZ[printTimestep]-180,tzone = "UTC"), 
                to = lubridate::force_tz(xx$nwm.discharge.dateTimeZoneZ[printTimestep]+180,tzone = "UTC"), axis = "x", color = "#FFEAEA") %>%
      dyShading(from = Sys.time()-180, to = Sys.time()+180, axis = "x", color = "#CCEBD6") %>%
      dyLegend(width = 650) %>%
      dyCSScoolMin()
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1100
    height <- 400
    webshot::webshot("temp.html", file = filename,cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height )
  }
  swiperchart_print_image <- function(printTimestep,workingdir,filename) {
    setwd(workingdir)
    ifelse(input$uislidertime=="Zulu", timeseriesOrder <- as.POSIXct(xx$nwm.discharge.dateTimeZoneZ, tz = "GMT"),timeseriesOrder <- as.POSIXct(xx$nwm.discharge.dateTimeZoneLocal, tz = Sys.timezone()))
    xx$chartXTS.addts <- xts::xts(x=xx$floodsum.addimpact[1:user.forecast.timesteps], order.by= timeseriesOrder)
    xx$chartXTS.areats <- xts::xts(x=xx$floodsum.area[1:user.forecast.timesteps], order.by= timeseriesOrder)
    xx$chartXTS.summary = cbind(xx$chartXTS.addts, xx$chartXTS.areats) 
    names(xx$chartXTS.summary) <- c("Count of Impacted Addresses","Square miles of Inundated Land")
    
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    chart <- dygraphs::dygraph(data = xx$chartXTS.summary,
                               main = paste0("Forecast impacts for ", user.aoi.string),
                               xlab = "Date/Time") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = FALSE) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      dyShading(from = lubridate::force_tz(xx$nwm.discharge.dateTimeZoneZ[printTimestep]-180,tzone = "UTC"), 
                to = lubridate::force_tz(xx$nwm.discharge.dateTimeZoneZ[printTimestep]+180,tzone = "UTC"), axis = "x", color = "#FFEAEA") %>%
      dyLegend(width = 650) %>%
      dyCSScoolMin()
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1100
    height <- 400
    webshot::webshot("temp.html", file = filename,cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height )
    unlink(paste0(workingdir,"/temp_files"), recursive=TRUE)
    unlink(paste0(workingdir,"/temp.html"))
    setwd(basedir)
  }
  
  # Tables ========================================
  output$mapindextable <- DT::renderDataTable({
    DT::datatable(as.data.frame(xx$mapindex.point[which(xx$mapindex.point$`Wet House Hours`>0),]) %>% select(`Index Label`,`Wet House Hours`,`Uniquely Impacted Addresses`,`Maximum Impact Depth`), 
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$addressimpacts <- DT::renderDataTable({
    DT::datatable(as.data.frame(xx$address.mapbook[which(xx$address.mapbook$ffreq>0),]) %>% select(address,`Index Label`,dataset), 
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$roadimpacts <- DT::renderDataTable({
    DT::datatable(as.data.frame(xx$road.point.mapbook[which(xx$road.point.mapbook$ffreq>0),]) %>% select(`Index Label`,name), 
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$filteraddressimpacts <- DT::renderDataTable({
    DT::datatable(as.data.frame(values$currentIndexAdds), 
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  output$filterroadimpacts <- DT::renderDataTable({
    DT::datatable(as.data.frame(values$currentIndexRoads), 
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  
  # Events ========================================
  observeEvent(input$summaryviewlayer, {
    shiny::req(input$gridColor)
    shiny::req(input$aoiColor)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(input$summaryviewlayer=="Individual Points") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        showGroup("Addresses")
    } else if(input$summaryviewlayer=="Wet House Hours") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Wet House Hours", pal = xx$mapindex.poly.WHHpal, values = subset(xx$mapindex.poly, `Wet House Hours`>0)$`Wet House Hours`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Count of Impacted Addresses", pal = xx$mapindex.poly.IHCpal, values = subset(xx$mapindex.poly, `Uniquely Impacted Addresses`>0)$`Uniquely Impacted Addresses`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Maximum Impact Depth") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Maximum Impact Depth", pal = xx$mapindex.poly.MDpal, values = subset(xx$mapindex.poly, `Maximum Impact Depth`>0)$`Maximum Impact Depth`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    }
  })
  observeEvent(input$Indexviewlayer, {
    shiny::req(input$gridColor)
    shiny::req(input$aoiColor)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(input$Indexviewlayer=="Individual Points") {
      leafletProxy("index_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        showGroup("Addresses")
    } else if(input$Indexviewlayer=="Current index") {
      leafletProxy("index_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~"Red",fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Indexviewlayer")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("index_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$Indexviewlayer=="Maximum Impact Depth") {
      leafletProxy("index_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    }
  })
  
  observeEvent((input$gridColor | input$aoiColor), {
    shiny::req(input$summaryviewlayer)
    shiny::req(input$timeslider)
    shiny::req(input$uislidertime)
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    
    if(input$summaryviewlayer=="Individual Points") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        showGroup("Addresses")
    } else if(input$summaryviewlayer=="Wet House Hours") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Maximum Impact Depth") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        hideGroup("Addresses")
    }
    leafletProxy("swiper_map") %>%
      clearGroup("Borders") %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
      hideGroup("Addresses")
  })
  observeEvent(input$userbasemap, {
    leafletProxy("summary_map") %>%
      clearGroup("Base Map") %>%
      addProviderTiles(input$userbasemap, group = "Base Map")
    
    leafletProxy("swiper_map") %>%
      clearGroup("Base Map") %>%
      addProviderTiles(input$userbasemap, group = "Base Map")
    
  })
  observeEvent(input$backbutton, {
    shiny::req(input$timeslider)
    shiny::req(input$uislidertime)
    shiny::req(input$gridColor)
    shiny::req(input$aoiColor)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    
    # if( nrow(subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0))>0 & (values$currentIndexLabel!="None") ) {
    if( nrow(subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0))>0) {
      if(values$count != 1) {
        values$count <- values$count - 1
      } else {
        values$count <- length(xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),])
      }
      values$currentIndexLabel <- subset(xx$mapindex.point,`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),]$`Index Label`)[values$count,]$`Index Label`
      values$currentIndexAdds <- subset(xx$address.mapbook, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0 & xx$address.mapbook$`Index Label`==values$currentIndexLabel) %>% select(address,`Index Label`,dataset)
      values$currentIndexRoads <- subset(xx$road.point.mapbook, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0 & xx$road.point.mapbook$`Index Label`==values$currentIndexLabel) %>% select(`Index Label`,name)
      mybbox <- AOI::bbox_coords(sf::st_buffer(
        x=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),][values$count,],
        dist=0.0002,
        nQuadSegs = 1,
        endCapStyle = "SQUARE",
        joinStyle = "ROUND",
        mitreLimit = 1,
        singleSide = FALSE
      ))
      leafletProxy("index_map") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 2,options=pathOptions(pane="index"),group="Index") %>% 
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),][values$count,], fill=F,color="Red",weight = 3,options=pathOptions(pane="index"),group="Index")
      leafletProxy("swiper_map") %>%
        flyToBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
    } else {
      values$currentIndexLabel = "No impacts in timestep"
    }
  })
  observeEvent(input$forwardbutton, {
    shiny::req(input$timeslider)
    shiny::req(input$uislidertime)
    shiny::req(input$gridColor)
    shiny::req(input$aoiColor)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    
    if( nrow(subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0))>0 ) {
      if(values$count != nrow(xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),])) {
        values$count <- values$count + 1
      } else {
        values$count <- 1
      }
      values$currentIndexLabel <- subset(xx$mapindex.point,`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),]$`Index Label`)[values$count,]$`Index Label`
      values$currentIndexAdds <- subset(xx$address.mapbook, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0 & xx$address.mapbook$`Index Label`==values$currentIndexLabel) %>% select(address,`Index Label`,dataset)
      values$currentIndexRoads <- subset(xx$road.point.mapbook, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0 & xx$road.point.mapbook$`Index Label`==values$currentIndexLabel) %>% select(`Index Label`,name)
      mybbox <- AOI::bbox_coords(sf::st_buffer(
        x=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),][values$count,],
        dist=0.0002,
        nQuadSegs = 1,
        endCapStyle = "SQUARE",
        joinStyle = "ROUND",
        mitreLimit = 1,
        singleSide = FALSE
      ))
      leafletProxy("index_map") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 2,options=pathOptions(pane="index"),group="Index") %>% 
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),][values$count,], fill=F,color="Red",weight = 3,options=pathOptions(pane="index"),group="Index")
      leafletProxy("swiper_map") %>%
        flyToBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
    } else {
      values$currentIndexLabel = "No impacts in timestep"
    }
  })
  observeEvent(input$uislidertime, {
    choices <- switch(
      input$uislidertime,
      "Zulu" = xx$nwm.discharge.dateTimeZoneZHR,
      "Local" = xx$nwm.discharge.dateTimeZoneLocalHR
    )
    updateSliderTextInput(
      session = session,
      inputId = "timeslider",
      choices = choices
    )
  })
  observeEvent(input$timeslider, {
    shiny::req(input$uislidertime)
    values$count <- 1
    values$currentIndexLabel <- "None"
    values$currentIndexAdds <- NULL
    values$currentIndexRoads <- NULL
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    
    if(nrow(subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0))==0) {
      leafletProxy("swiper_map") %>%
        clearGroup("Flooding") %>%
        clearGroup("Road Impacts") %>%
        clearGroup("Addresses") %>% 
        clearGroup("Index") %>%
        clearGroup("Index Labels") %>%
        addRasterImage(x=xx$flood.grid[[match(input$timeslider, timestepseries)]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addLegend("bottomleft", title = "Flood Depth",pal = xx$flood.grid.pal, values = xx$flood.grid.val, group = "Flooding") 
      # addPolylines(data=xx$road.poly[[match(input$timeslider, timestepseries)]],color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts")
      
      leafletProxy("index_map") %>%
        clearGroup("Flooding") %>%
        clearGroup("Road Impacts") %>%
        clearGroup("Addresses") %>% 
        clearGroup("Index") %>%
        clearGroup("Index Labels")
    } else {
      leafletProxy("swiper_map") %>%
        clearGroup("Flooding") %>%
        clearGroup("Road Impacts") %>%
        clearGroup("Addresses") %>% 
        clearGroup("Index") %>%
        clearGroup("Index Labels") %>%
        addRasterImage(x=xx$flood.grid[[match(input$timeslider, timestepseries)]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$road.poly[[match(input$timeslider, timestepseries)]],color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
        addAwesomeMarkers(data=subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),
                          icon=xx$address.point.standard, 
                          options=pathOptions(pane="add"),
                          # clusterOptions=markerClusterOptions(),
                          group="Addresses") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),], 
                    fill=F,
                    color=gridcolor,
                    fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),
                    # fillOpacity=1,
                    weight = 1,
                    options=pathOptions(pane="index"),
                    group="Index") %>%
        addLabelOnlyMarkers(data=subset(xx$mapindex.point,`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),]$`Index Label`), 
                            label = ~`Index Label`, 
                            labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),
                            options=pathOptions(pane="indexlabels"),group="Index Labels")
      
      leafletProxy("index_map") %>%
        clearGroup("Flooding") %>%
        clearGroup("Road Impacts") %>%
        clearGroup("Addresses") %>% 
        clearGroup("Index") %>%
        clearGroup("Index Labels") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),][1,],fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),]$`Index Label`)),label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index")
      
      values$count <- 0
      input$forwardbutton
    }
  })
  observeEvent(input$sidebar, {
    if(input$sidebar == "Swiper"){
      shiny::req(input$uislidertime)
      ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
      updateSliderTextInput(session = session,inputId="timeslider", selected = timestepseries[2])
      updateSliderTextInput(session = session,inputId="timeslider", selected = timestepseries[1])
    }
  })
  
  observeEvent(input$mapindextable_rows_selected, {
    leafletProxy("summary_map") %>%
      flyTo(xx$mapindex.point[which(xx$mapindex.point$`Uniquely Impacted Addresses`>0),][input$mapindextable_rows_selected,]$geometry[[1]][1],
            xx$mapindex.point[which(xx$mapindex.point$`Uniquely Impacted Addresses`>0),][input$mapindextable_rows_selected,]$geometry[[1]][2],17)
  })
  observeEvent(input$addressimpacts_rows_selected, {
    leafletProxy("summary_map") %>%
      flyTo(xx$address.mapbook[which(xx$address.mapbook$ffreq>0),][input$addressimpacts_rows_selected,]$geometry[[1]][1],
            xx$address.mapbook[which(xx$address.mapbook$ffreq>0),][input$addressimpacts_rows_selected,]$geometry[[1]][2],17)
  })
  observeEvent(input$roadimpacts_rows_selected, {
    leafletProxy("summary_map") %>%
      flyTo(xx$road.point.mapbook[which(xx$road.point.mapbook$ffreq>0),][input$roadimpacts_rows_selected,]$geometry[[1]][1],
            xx$road.point.mapbook[which(xx$road.point.mapbook$ffreq>0),][input$roadimpacts_rows_selected,]$geometry[[1]][2],17)
  })
  # observeEvent(input$relaunch, {
  #   session$reload()
  # })
  observeEvent(input$save, {
    exportnames <- as.character(xx$nwm.discharge.dateTimeZoneZ) %>% stringr::str_sub(6, 13) %>% paste0("Z.tif")
    exportnames[user.forecast.timesteps+1] <- "ffreq.tif"
    setwd(paste0(basedir,"/AOI/",user.aoi.filepath,"/output"))
    raster::writeRaster(xx$flood.grid, filename=exportnames,overwrite=TRUE,bylayer=TRUE,format="GTiff")
    setwd(basedir)
    shinyalert(title = paste0("-- Flood innundation rasters written to: ",basedir,"/AOI/",user.aoi.filepath,"/output --"), type = "success")
  })
  observeEvent(input$print, {
    outputdir <- paste0(basedir,"/AOI/",user.aoi.filepath,"/output")
    
    # Title Page 
    TitlePageTitleBlock <- magick::image_read(paste0(basedir,"/data/misc/EmptyFrontPageBlock.png"))
    if(user.aoi.source=="string") {
      TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, user.aoi.string, font = 'Georgia', location = "+158+328", size = 20)
    } else {
      TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, paste(user.aoi.source,user.aoi.string), font = 'Georgia', location = "+158+328", size = 20)
    }
    if(user.output.choice=="basedata") {
      TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, "Base data map books", font = 'Georgia', location = "+215+357", size = 20)
    } else {
      TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, "Flood impact mapping", font = 'Georgia', location = "+215+357", size = 20)
    }
    TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, format(file.mtime(xx$aoi.path),format='%m/%d/%Y'), font = 'Georgia', location = "+418+383", size = 20) %>% 
      magick::image_annotate(user.address.source, font = 'Georgia', location = "+301+410", size = 20) %>%
      magick::image_annotate(user.road.source, font = 'Georgia', location = "+627+410", size = 20) %>%
      magick::image_trim()
    mybbox <- suppressWarnings(suppressMessages(AOI::bbox_coords(sf::st_buffer(x=xx$aoi.shp,dist=0.0001,nQuadSegs=1,endCapStyle="SQUARE",joinStyle="ROUND",mitreLimit=1,singleSide=FALSE))))
    # mybigbbox <- suppressWarnings(suppressMessages(AOI::bbox_coords(sf::st_buffer(x=xx$aoi.shp,dist=0.004,nQuadSegs=1,endCapStyle="SQUARE",joinStyle="ROUND",mitreLimit=1,singleSide=FALSE))))
    boundrydatapoly <- suppressWarnings(suppressMessages(featurefile[sf::st_buffer(xx$aoi.shp,dist=0.002,nQuadSegs=1,endCapStyle="SQUARE",joinStyle="ROUND",mitreLimit=1,singleSide=FALSE),]))
    boundrydatapoint <- suppressWarnings(suppressMessages(boundrydatapoly %>% st_centroid()))
    quiet(ifelse(user.aoi.source=="zctas",insetlablefield<-boundrydatapoint$ZCTA5CE10,boundrydatapoint<-featurefile$HUC8))
    mapshot(leaflet(height=500, width=500) %>%
              addProviderTiles("Esri.WorldStreetMap") %>%
              addPolygons(data=boundrydatapoly,color="black",fillColor="black",fillOpacity=0.08,weight=1,labelOptions=labelOptions(noHide=T,direction="bottom")) %>%
              addPolygons(data=xx$aoi.shp,color="red",fill=FALSE,weight=5) %>%
              # addLegend("bottomleft", title = "Processed Area", pal = colorFactor(palette = "red",domain = "AOI"), values = "AOI") %>%
              addLabelOnlyMarkers(data=boundrydatapoint, label =insetlablefield, labelOptions = labelOptions(noHide=T,sticky=T,textOnly=TRUE,direction="bottom")) %>%
              fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
            ,file=paste0(outputdir, '/FrontPageInset.png'))   
    mapshot(leaflet(height=1000, width=1000) %>%
              addProviderTiles("Esri.WorldImagery") %>%
              addPolygons(data=xx$aoi.shp,color="red",fill=FALSE,weight=3) %>%
              addPolylines(data=xx$flow.line,weight = ~streamorde/2,group="Flowlines") %>%
              # fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax) %>%
              addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines") %>%
              addLegend("bottomright", title = "Processed Area",pal = colorFactor(palette = "Red",domain = "AOI"), values = "AOI") %>%
              # addScaleBar(position = "bottomleft", options = scaleBarOptions(maxWidth = 500,imperial = TRUE)) %>%
              fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
            ,file=paste0(outputdir, '/FrontPageMap.png'))   
    blankpage <- magick::image_read(paste0(basedir,'/data/misc/Empty8x11.png'))
    AOIMapInset <- magick::image_read(paste0(outputdir, '/FrontPageInset.png')) %>%
      magick::image_trim() %>%
      magick::image_scale("750x")
    frontpagemap <- magick::image_read(paste0(outputdir, '/FrontPageMap.png'))
    AOIMapImage <- magick::image_composite(blankpage, magick::image_scale(frontpagemap,"2250x"), offset = "+150+900")
    AOIMapImage <- magick::image_composite(AOIMapImage, AOIMapInset, offset = "+1650+150")
    AOIMapImage <- magick::image_composite(AOIMapImage, magick::image_scale(TitlePageTitleBlock,"1500x"), offset = "+150+150")
    unlink(paste0(outputdir,'/FrontPageInset.png'))
    unlink(paste0(outputdir,'/FrontPageMap.png'))
    magick::image_write(AOIMapImage, path = paste0(outputdir,'/00_TitlePage.png'), format = "png")
    
    # Summary page
    mapshot(mapshot_summary_map(), file=paste0(outputdir, '/summary_map.png'))   
    summarymap <- magick::image_read(paste0(outputdir,'/summary_map.png'))
    summararychart_image(outputdir,'summary_graph_table.png')  
    # summarygraph <- magick::image_read(paste0(outputdir,'/summary_graph_table.png'))
    summarygraph <- magick::image_read(paste0(outputdir,'/summary_graph_table.png')) %>%
      magick::image_trim()
    ggsave(
      filename=paste0(outputdir, '/summary_map_table.png'),
      plot = ggplot() +
        annotation_custom(tableGrob(as.matrix(as.data.frame(xx$mapindex.point[which(xx$mapindex.point$`Wet House Hours`>0),]) %>% select(`Index Label`,`Wet House Hours`,`Uniquely Impacted Addresses`,`Maximum Impact Depth`)), vp = viewport(width=10))) +
        theme_void(),width=10,units="in",dpi = 300)
    summarytable <- magick::image_read(paste0(outputdir,'/summary_map_table.png')) %>%
      magick::image_trim()
    magick::image_write(AOIMapImage, path = paste0(outputdir,'/00_TitlePage.png'), format = "png")
    titleblockimage <- magick::image_read(paste0(outputdir,'/titleblock.png'))
    page1 <- magick::image_composite(blankpage, magick::image_scale(titleblockimage,"975x"), offset = "+150+150")
    page1 <- magick::image_composite(page1,
                                     magick::image_scale(summarymap, "2250x") %>% image_border("black", "4x4"), 
                                     offset = "+150+900")
    page1 <- magick::image_composite(page1,
                                     magick::image_scale(summarygraph, "967x") %>% image_border("black", "4x4"), 
                                     offset = "+150+533")
    page1 <- magick::image_composite(page1,
                                     magick::image_scale(summarytable, "1000x") %>% image_border("black", "4x4"), 
                                     offset = "+1400+150")
    
    unlink(paste0(outputdir,'/summary_graph_table.png'))
    unlink(paste0(outputdir,'/summary_map.png'))
    unlink(paste0(outputdir,'/summary_map_table.png'))
    magick::image_write(page1, path = paste0(outputdir,'/01_summarypage.png'), format = "png")
    
    for(i in 1:length(xx$nwm.discharge.dateTimeZoneLocalHR)) {
      workingpage <- magick::image_composite(blankpage, magick::image_scale(titleblockimage,"975x"), offset = "+150+150")
      swiperchart_print_image(i,outputdir,paste0("swiper_",i,"_graph.png"))
      swiperchart <- magick::image_read(paste0(outputdir,"/swiper_",i,"_graph.png")) %>%
        magick::image_trim()
      workingpage <- magick::image_composite(workingpage, magick::image_scale(swiperchart, "967x") %>% image_border("black", "4x4"), offset = "+150+533")
      
      if(nrow(subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0))>0) {
        
        for(j in 1:length(xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),])) {
          mapshot(mapshot_index_map(i,xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),][j,]$'Index Label',mybbox), file=paste0(outputdir, '/swiper_index_',i,'_',j,'_.png'))
          pindexmap <- magick::image_read(paste0(outputdir,'/swiper_index_',i,'_',j,'_.png')) %>%
            magick::image_scale("750x")
          
          mapshot(mapshot_swiper_map(i,xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),][j,]$'Index Label'), file=paste0(outputdir, '/swiper_',i,'_',j,'_.png'))
          pswipermap <- magick::image_read(paste0(outputdir,'/swiper_',i,'_',j,'_.png')) %>%
            magick::image_trim() %>%
            magick::image_scale("x731")
          
          tg = gridExtra::tableGrob(as.matrix(
            as.data.frame(
              subset(xx$address.mapbook, eval(parse(text = names(xx$flood.grid)[i]))>0 & xx$address.mapbook$`Index Label`==xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),][j,]$'Index Label')  ) %>%
              select(address)))
          h = grid::convertHeight(sum(tg$heights), "in", TRUE)
          w = grid::convertWidth(sum(tg$widths), "in", TRUE)
          ggplot2::ggsave(paste0(outputdir, '/swiper_impact_',i,'_',j,'_.png'), tg, width=w, height=h)
          swipertable <- magick::image_read(paste0(outputdir,'/swiper_impact_',i,'_',j,'_.png')) %>%
            magick::image_trim()
          
          workingpage_new <- magick::image_composite(workingpage, magick::image_scale(pswipermap, "2250x") %>% image_border("black", "4x4"), offset = "+150+900")
          workingpage_new <- magick::image_composite(workingpage_new, pindexmap %>% image_border("black", "4x4"),offset = "+1125+150")
          workingpage_new <- magick::image_composite(workingpage_new,magick::image_scale(swipertable, "270x") %>% image_border("black", "4x4"), offset = "+2114+248")
          workingpage_new <- magick::image_annotate(workingpage_new, paste0("Impacts for: ",xx$nwm.discharge.dateTimeZoneZHR[i]), font = 'Georgia', location = "+1891+150", size = 30)
          workingpage_new <- magick::image_annotate(workingpage_new, paste0("Index: ",xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),][j,]$'Index Label'), font = 'Georgia', location = "+2180+182", size = 30)
          workingpage_new <- magick::image_annotate(workingpage_new, paste0("Impacted Addresses:"), font = 'Georgia', location = "+2094+214", size = 30)
          magick::image_write(workingpage_new, path = paste0(outputdir,'/02_impactpage_',i,'_',j,'_.png'), format = "png")
          
          unlink(paste0(outputdir,'/swiper_impact_',i,'_',j,'_.png'))
          unlink(paste0(outputdir,'/swiper_index_',i,'_',j,'_.png'))
          unlink(paste0(outputdir,'/swiper_',i,'_',j,'_.png'))
        }
        
      } else {
        # no impact in timestep
        mapshot(mapshot_index_map(i,NULL,mybbox), file=paste0(outputdir, '/swiper_',i,'_0.png'))
        pswipermap <- magick::image_read(paste0(outputdir,'/swiper_',i,'_0.png'))
        
        mapshot(mapshot_swiper_map(i,NULL), file=paste0(outputdir, '/swiper_index_',i,'_0.png'))
        pindexmap <- magick::image_read(paste0(outputdir,'/swiper_index_',i,'_0.png')) %>%
          magick::image_trim() %>%
          magick::image_scale("x731")
        
        workingpage_new <- magick::image_composite(workingpage,
                                                   magick::image_scale(pswipermap, "2250x") %>% image_border("black", "4x4"), 
                                                   offset = "+150+900")
        workingpage_new <- magick::image_composite(workingpage_new, pindexmap %>% image_border("black", "4x4"), offset = "+1125+150")
        workingpage_new <- magick::image_annotate(workingpage_new, paste0("Impacts for: ",xx$nwm.discharge.dateTimeZoneZHR[i]), font = 'Georgia', location = "+1891+150", size = 30)
        workingpage_new <- magick::image_annotate(workingpage_new, paste0("Index: NONE"), font = 'Georgia', location = "+2180+182", size = 30)
        magick::image_write(workingpage_new, path = paste0(outputdir,'/02_impactpage_',i,'_.png'), format = "png")
        
        unlink(paste0(outputdir,'/swiper_',i,'_0.png'))
        unlink(paste0(outputdir,'/swiper_index_',i,'_0.png'))
      }
      unlink(paste0(outputdir,"/swiper_",i,"_graph.png"))
    }
    
    all_images = list.files(outputdir, full.names = TRUE, pattern = '.png')
    all_images <- all_images[-length(all_images)]
    all_images_1 <- purrr::reduce(purrr::map(all_images,image_read),c)
    image_write(all_images_1 , format = "pdf", paste0(outputdir,"/FOSSFlood_Impact_Output_",user.forecast.gen[1],"_",user.forecast.gen[2],"Z.pdf"))
    file.remove(all_images)
    
    # Alert to finish :)
    shinyalert(title = paste0("See ",user.aoi.filepath,"/output"), type = "success")
  })
  
  cleanFiles <- function() {
    unlink(paste0(basedir,"/temp_files"), recursive=TRUE)
    unlink(paste0(basedir,"/temp.html"))
  }
  
  # session$allowReconnect(TRUE)
  session$onSessionEnded(function() {
    shiny::stopApp()
  })
}
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////



#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
# ----- Server launch ------------------------------------------------------------------------------------------
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////
rRuntime <- difftime(Sys.time(), rStarttime, units = "mins")
print(paste0("-- FOSSFlood execution complete in ",round(rRuntime, digits = 2), " minutes - Starting Shiny server --"))
# setwd(paste0(basedir,'/shiny'))  # needed for two and three file apps

browser.chromeium = file.path(paste0(basedir, '/ChromiumPortable/App/Chromium/64/chrome.exe'))
launch.browser = function(appUrl, browser.path=browser.chromeium) {
  system(sprintf('"%s" --disable-gpu --app="data:text/html,<html>
                 <head>
                 <title>System Configuration</title>
                 </head>
                 <body>
                 <script>window.resizeTo(830,675);window.location=\'%s\';</script>
                 </body></html>" &', browser.path, appUrl), wait=FALSE)
}

# shinyWidgets::shinyWidgetsGallery()
if(user.output.choice=="basedata") {
  finalUI <- basedataUI
  finalSERVER <- basedataSERVER
} else {
  finalUI <- impactUI
  finalSERVER <- impactSERVER
}

shinyApp(
  ui = finalUI,
  server = finalSERVER,
  options = list(launch.browser=launch.browser)
)
#grep("^ChromiumPortable",readLines(textConnection(system('tasklist',intern=TRUE))),value=TRUE)
