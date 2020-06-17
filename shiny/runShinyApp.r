rStarttime <- Sys.time()
# setwd("C:/Users/Cornholio/Desktop/FOSSFlood-master")
basedir <- getwd()   # This should look something like C:/Users/.../FOSSFlood-master

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

# Welcome to the heart of FOSSFlood, R script which will do the bulk of the heavy lifting for you.  This header section is where the HTA fired VBS replaces strings for user interactions.  
# While the core of FOSSFlood is platform independent, the HTA was needed to create a user interface and thus FOSSFlood is currently limited to the Windows OS.  However, if you want to
# run FOSSFlood on a different OS, all you need to do is...
# * Install R studio
# * In RStudio, point the R installation to the FOSSFlood R-Portable executable (FOSSFlood-master\R-Portable\App\R-Portable\bin\R.exe)
# * Replace code below with desired inputs,
# * and run the entire body of the code :)

# -- USER Inputs -------------------------------------------------------------------------
user.aoi.string <- "66044, 66046, 66047, 66045, 66049"
user.aoi.source <- "zctas" 										# "zctas", "huc8", "string"
user.address.source <- "OpenAddresses" 						# "OpenStreetMap_Addresses", "OpenAddresses", "User_Provided_Addresses"
user.address.file <- "" 							# Empty "" or filepath to addresses
user.road.source <- "TIGER_Lines_2018"  							# "TIGER_Lines_2018", "OpenStreetMaps", "User_Provided_Roads"
user.road.file <- "" 									# Empty or filepath to addresses
user.forecast.source <- "NWM_SR_C" 						# "NWM_SR_C", USER_DIS, USER_STAGE
user.forecast.timesteps <- as.numeric("6") 	# number
user.forecast.file <- "" 							# Empty "" or filepath to flows.fst
user.forecast.members <- as.list("##USERFORECASTMEMBERS") 			# unused - DEV
user.output.choice <- "impacts" 							# GIS_O  basedata  impacts
user.output.grid <- "Square" 								# "Square", "Hexagon"
user.output.hardclip <- as.logical("False") 					# If TRUE, Hard clip data to aoi shape, defaults to bb
user.output.archive <- as.logical("False") 						# Save flows in output folder of the requested AOI using timestamp as file name.  Can be pointed back to later to regenerate outputs

# Input cleanup
user.address.source <- stringr::str_replace_all(user.address.source,"_"," ")
user.road.source <- stringr::str_replace_all(user.road.source,"_"," ")

# -- Dev Comment/uncomment with ctrl-shift-c
# user.aoi.string <- "66044, 66046, 66047, 66045, 66049"
# # user.aoi.string <- "03216"
# user.aoi.source <- "zctas"
# user.address.source <- "OpenAddresses"
# user.address.file <- ""
# user.road.source <- "TIGER Lines 2018"
# user.road.file <- ""
# user.forecast.source <- "NWM_SR_C"
# user.forecast.timesteps <- as.numeric("8")
# user.forecast.file <- ""
# user.forecast.members <- as.list(1)
# user.output.choice <- "impacts"
# user.output.grid <- "Square"
# user.output.hardclip <- TRUE
# user.output.archive <- FALSE

print(paste("-- Welcome to FOSSFlood - Running FOSSFlood in", basedir)) # This should look something like C:/Users/.../FOSSFlood-master
print("-- Pre-Preloading constants --")

##  _______       ___      .__   __.   _______  _______ .______          ________    ______   .__   __.  _______ 
## |       \     /   \     |  \ |  |  /  _____||   ____||   _  \        |       /   /  __  \  |  \ |  | |   ____|
## |  .--.  |   /  ^  \    |   \|  | |  |  __  |  |__   |  |_)  |       `---/  /   |  |  |  | |   \|  | |  |__   
## |  |  |  |  /  /_\  \   |  . `  | |  | |_ | |   __|  |      /           /  /    |  |  |  | |  . `  | |   __|  
## |  '--'  | /  _____  \  |  |\   | |  |__| | |  |____ |  |\  \----.     /  /----.|  `--'  | |  |\   | |  |____ 
## |_______/ /__/     \__\ |__| \__|  \______| |_______|| _| `._____|    /________| \______/  |__| \__| |_______|
##  
## Edit below this line at your own risk, well, as risky as coding can be...
## Note: If during your edits FOSSFlood breaks, the easiest way to recover is simply to download FOSSFlood again.

# Last build date: 3/14/2020
# install.packages("tidyverse")
# install.packages("installr")
# install.packages("devtools")
# install.packages("geosphere")
# install.packages("sf")
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
suppressMessages(library(nomadsNC))
suppressMessages(library(velox))
suppressMessages(library(geofabrik))

# suppressMessages(library(patchwork))
# suppressMessages(library(hrbrthemes))
options(tigris_use_cache = FALSE)


#/////////////////////////////////////
# Helper Functions                                     
#/////////////////////////////////////
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
  
  # -- Check to make sure I have the right date --------------------------------------------------------------------------- 
  testurl <- paste0('http://nomads.ncep.noaa.gov/pub/data/nccf/com/nwm/prod/nwm.', fileDate, '/short_range')
  if (http_error(testurl)) { 
    fileDate <- gsub("-", "",format(as.Date(fulltimestamp, tz = "GMT")-1, format = "%Y-%m-%d"))
  }
  # -- Guess at the most recent time ---------------------------------------------------------------------------------------
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
#/////////////////////////////////////
#/////////////////////////////////////


#/////////////////////////////////////
# Block below runs only if this is the first time you have run FOSSFlood, unpacks and creates the zip code files, installs needed dependencies for printing
#/////////////////////////////////////
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
if (is.null(rmarkdown::pandoc_available())) {
  print("-- Pandoc is needed in order to print outputs, please follow the installation prompts as they appear.")
  installr::install.pandoc(use_regex = TRUE, to_restart = FALSE)
}
#/////////////////////////////////////
#/////////////////////////////////////


#/////////////////////////////////////
# Parse user inputs
#/////////////////////////////////////
print(paste("-- Welcome to FOSSFlood - Running FOSSFlood for", user.aoi.string))
if(user.aoi.source %in% c("zctas", "huc8")) {
  user.aoi.stringlist <- as.list(strsplit(user.aoi.string, ", ")[[1]])
  user.aoi.filepath <- gsub(", ", "_", user.aoi.string)
  if(user.aoi.source=="zctas") {
    user.aoi.filepath <- paste0('zctas_', user.aoi.filepath)
    featurefile <- sf::read_sf(paste0(basedir, "/data/misc/zctas.shp"))
    user.aoi.call <- base::subset(featurefile, (featurefile$ZCTA5CE10 %in% user.aoi.stringlist))
  } else if(user.aoi.source=="huc8") {
    user.aoi.filepath <- paste0('huc8_', user.aoi.filepath)
    featurefile <- sf::read_sf(paste0(basedir, "/data/misc/huc8.shp"))
    user.aoi.call <- base::subset(featurefile, (featurefile$HUC8 %in% user.aoi.stringlist))
  }
} else if(user.aoi.source == "string") {
  user.aoi.call <- user.aoi.string
  user.aoi.stringlist <- as.list(strsplit(user.aoi.string, ", ")[[1]])
  user.aoi.filepath <- gsub(", ", "_", user.aoi.stringlist)
}

# Run AOI and Floodmapping packages
AOI = AOI::aoi_get(user.aoi.call)
tryCatch({quiet(geosphere::areaPolygon(sf::as_Spatial(AOI$geometry))>0)}, 
         error = function(e) {
           print("-- Malformed input: The features were not located in the AOI.  Double check that your features exist and that you entered them in correctly.")
           stop(e)
         })
xx = quiet(FloodMapping::getRawData(AOI,paste0(basedir,"/AOI/"),user.aoi.filepath))

# generate base data as needed
if (!file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/grid_rec.shp"))) {
  # Setup ----------------------------------------------------------------------------
  dir.create(paste0(basedir,"/AOI/", user.aoi.filepath,"/tmp"))
  dir.create(paste0(basedir,"/AOI/", user.aoi.filepath,"/output"))
  
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
  URL <- paste0("https://cida.usgs.gov/nwc/geoserver/nhdplus/ows?service=WFS&version=2.0.0&request=GetFeature&typeNames=nhdplus:nhdflowline_network&srsName=EPSG:4326&bbox=",
                south,",",west,",",north,",",east,
                "&outputFormat=SHAPE-ZIP")
  quiet(httr::GET(URL, write_disk(paste0(basedir, "/AOI/",user.aoi.filepath,"/tmp/Flowlines.zip"), overwrite=TRUE)))
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
  # Error correcting: not set for areas that cross state boundaries
  print("-- Downloading TIGER roads --")
  CountiesFIPS <- quiet(tigris::counties(cb=TRUE))
  CountiesFIPS_Proj <- sf::st_transform(sf::st_as_sf(CountiesFIPS), sf::st_crs(4326))
  aoiCounties <- suppressMessages(CountiesFIPS_Proj[xx$aoi.bb, ])
  # Error checking: Download roads for AOI that interset more than one state or county
  if(length(unique(aoiCounties$COUNTYFP)) > 1) {
    tigerroads_tmp <- quiet(tigris::roads(unique(aoiCounties$STATEFP),unique(aoiCounties$COUNTYFP)[1], year = 2018, refresh = TRUE)[NULL,])
    for(i in unique(aoiCounties$COUNTYFP)) {
      tigerroads <- quiet(tigris::roads(unique(aoiCounties$STATEFP),i, year = 2018, refresh = TRUE))
      tigerroads_tmp <- rbind(tigerroads_tmp, tigerroads)
    }
    tigerroads <- tigerroads_tmp
  } else
    tigerroads <- quiet(tigris::roads(unique(aoiCounties$STATEFP),unique(aoiCounties$COUNTYFP), year = 2018, refresh = TRUE))
  
  tigerroads_Proj <- sf::st_transform(sf::st_as_sf(tigerroads), sf::st_crs(4326))
  tigerroads_Proj_Sub <- suppressMessages(tigerroads_Proj[xx$aoi.bb, ])
  var.out.bool <- names(tigerroads_Proj_Sub) %in% c("FULLNAME", "RTTYP")
  tigerout <- tigerroads_Proj_Sub[,var.out.bool] %>%
    dplyr::rename('name' = FULLNAME) %>%
    dplyr::rename('type' = RTTYP)
  sf::write_sf(tigerout, paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_tiger.shp"), delete_layer = TRUE, quiet = TRUE)
  # ---------------------------------------------------------------------------------
  
  # OSM Downloads ---------------------------------------------------------------------------------
  # Error correcting: not set for areas that cross state boundaries
  print("-- Downloading OSM Dataset --")
  quiet(httr::GET(paste0("http://download.geofabrik.de/north-america/us/",tolower(gsub(" ", "-", cdlTools::fips(aoiCounties$STATEFP[1], to = "Name"))),"-latest.osm.pbf"), 
                  write_disk(paste0(basedir,"/AOI/",user.aoi.filepath,"/tmp/",tolower(gsub(" ", "-", cdlTools::fips(aoiCounties$STATEFP[1], to = "Name"))),"-latest.osm.pbf"), 
                             overwrite=TRUE)))
  
  # OSM Roads ---------------------------------------------------------------------------------
  print("-- Building OSM Roads --")
  OSMlines <- geofabrik::read_pbf(paste0(basedir,"/AOI/",user.aoi.filepath,"/tmp/",
                                         tolower(gsub(" ", "-", cdlTools::fips(aoiCounties$STATEFP[1], to = "Name"))),
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
  sf::write_sf(OSMlinesout, paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_osm.shp"), delete_layer = TRUE, quiet = TRUE)
  # ---------------------------------------------------------------------------------
  
  # OSM addresses ---------------------------------------------------------------------------------
  print("-- Building OSM Addresses --")
  OSMpoints <- geofabrik::read_pbf(paste0(basedir,"/AOI/",user.aoi.filepath,"/tmp/",
                                          tolower(gsub(" ", "-", cdlTools::fips(aoiCounties$STATEFP[1], to = "Name"))),
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
  sf::write_sf(OSMpointsout, paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_osm.shp"), delete_layer = TRUE, quiet = TRUE)
  # ---------------------------------------------------------------------------------
  
  # Open Addresses ---------------------------------------------------------------------------------
  # Error correcting: not set for areas that cross state boundaries
  print("-- Downloading OpenAddresses --")
  validOA <- openadds::oa_list() %>% dplyr::filter(str_detect(processed, paste0("/us/",tolower(cdlTools::fips(aoiCounties$STATEFP[1], to = "Abbreviation")))))
  urls = validOA$processed 
  out = list()
  for( i in 1:length(validOA$processed)) {
    out[[i]] <- tryCatch({
      openadds::oa_get(validOA$processed[[i]])},
      error   = function(e){NULL})
    message(i)
  }
  # out <- out[-which(sapply(out, is.null))]
  c = parse(text = paste0("mergedStateData <- openadds::oa_combine(",str_c("out[[", c(1:(length(out))), "]]", sep = "",  collapse = ", "),")"))
  mergedStateData <- eval(c)
  sp::coordinates(mergedStateData) <- ~lon+lat
  sfoadata <- sf::st_as_sf(mergedStateData) %>% sf::st_set_crs(4326)
  aoiPoints <- suppressWarnings(suppressMessages(sfoadata[!duplicated(sfoadata),][xx$aoi.bb, ]))
  sf::write_sf(aoiPoints, paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_oa.shp"), delete_layer = TRUE, quiet = TRUE)
  # ---------------------------------------------------------------------------------
  
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
                 " road segments equivelent to ",toString(sum(sf::st_length(OSMlines_Proj_Sub))*0.000621371)," miles of roads in AOI")) } else { print("Issues generating osm data") }
  if(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/addresses_oa.shp"))) {
    print(paste0("No issues generating openaddress database: ",nrow(aoiPoints)," unique points")) } else { print("Issues generating openaddresses") }
  if(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/roads_tiger.shp"))) {
    print(paste0("No errors generating TIGER data: ",nrow(tigerroads_Proj_Sub)," road segments equivelent to ",toString(sum(sf::st_length(tigerroads_Proj_Sub))*0.000621371)," miles of roads in AOI")) } else { print("Issues generating TIGER roads") }
  if(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/NWISgages.shp"))) {
    print(paste0("No issues downloading NWIS data: ",nrow(NWISgages)," potential points in the AOI")) } else { print("Issues generating NWIS data") }
  if(file.exists(paste0(basedir,"/AOI/",user.aoi.filepath,"/nhdflowline_network.shp"))) {
    print(paste0("No errors generating NHD flowline data: ",nrow(goodlines)," segments equivelent to ",toString(sum(sf::st_length(goodlines))*0.000621371)," miles of waterways in AOI")) } else { print("Issues downloading flow lines") }
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
}

print("-- Welcome to FOSSFlood - Loading in data")
#/////////////////////////////////////
# AOI - load in base data and paths
#/////////////////////////////////////
xx$aoibb.path <- paste0(basedir,"/AOI/", user.aoi.filepath,"/",user.aoi.filepath,".shp")
if(user.aoi.source %in% c("zctas", "huc8")) {
  xx$aoi.path <- paste0(basedir,"/AOI/", user.aoi.filepath,"/",user.aoi.filepath,"_hardclip.shp")
} else {
  xx$aoi.path <- paste0(basedir,"/AOI/", user.aoi.filepath,"/",user.aoi.filepath,".shp")
}

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
#/////////////////////////////////////
#/////////////////////////////////////


#/////////////////////////////////////
# Grab and create Flows                                  
#/////////////////////////////////////
print("-- Downloading flows --")
# Flows in cms
if(user.forecast.source=="USER_DIS") {
  xx$nwm.flow = fst::read.fst(user.forecast.file)
  user.forecast.timesteps <- 6
} else if(user.forecast.source=="USER_STAGE") {
  xx$nwm.flow = fst::read.fst(xx$nwm.flow)
  user.forecast.timesteps <- 6
} else if(user.forecast.source=="NWM_SR_C") {
  user.forecast.gen <- getMostRecentForecast()
  xx$nwm.flow = quiet(nomadsNC::create_nomads_fst(type = "short_range", num = user.forecast.timesteps, dstfile = paste0(basedir,"/AOI/",user.aoi.filepath,"/flows.fst")))
}
if(user.output.archive) { write.fst(xx$nwm.flow, paste0(basedir,"/AOI/",user.aoi.filepath,"output/",user.forecast.source,getMostRecentForecast()[1],"_",getMostRecentForecast()[2],"_flows.fst")) }
xx$nwm.discharge <- fst::read.fst(xx$nwm.flow)
xx$nwm.discharge.dateTime <- names(xx$nwm.discharge)[-1] %>% 
  stringr::str_sub(3, 21)
xx$nwm.discharge.dateTimeZoneLocal <- as.POSIXct(xx$nwm.discharge.dateTime,tz=Sys.timezone())
xx$nwm.discharge.dateTimeZoneZ <- xx$nwm.discharge.dateTimeZoneLocal 
attr(xx$nwm.discharge.dateTimeZoneZ, "tzone") <- "UTC"
xx$nwm.discharge.dateTimeZoneLocalHR <- format(xx$nwm.discharge.dateTimeZoneLocal, '%m/%d/%Y %I:%M %p')
xx$nwm.discharge.dateTimeZoneZHR <- sapply(xx$nwm.discharge.dateTimeZoneZ, as.character) 

if(!SkipNWISFlows) {
  firstdate <- lubridate::date(xx$nwm.discharge.dateTimeZoneZ[1])
  lastdate <-  lubridate::date(xx$nwm.discharge.dateTimeZoneZ[user.forecast.timesteps])
  xx$nwis.discharge <- readNWISuv(xx$gage.point$site_no, "00060", as.Date(firstdate)-2, as.Date(lastdate)+1, tz="UTC")
  xx$nwis.stage <- readNWISuv(xx$gage.point$site_no, "00065", as.Date(firstdate)-2, as.Date(lastdate)+1, tz="UTC")
}

#/////////////////////////////////////
# Map inpacts                                 
#/////////////////////////////////////
print("-- Mapping Inundation Depths --")
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
print("-- Generating impacts --")
xx$flood.grid <- raster::addLayer(xx$flood.grid, sum(xx$flood.grid>=0, na.rm = TRUE))
names(xx$flood.grid)[raster::nlayers(xx$flood.grid)] <- 'ffreq'

#/////////////////////////////////////
#/////////////////////////////////////
if(user.output.choice=="GIS_O") {
  exportnames <- as.character(xx$nwm.discharge.dateTimeZoneZ) %>% stringr::str_sub(6, 13) %>% paste0("Z.tif")
  exportnames[user.forecast.timesteps+1] <- "ffreq.tif"
  setwd(paste0(basedir,"/AOI/",user.aoi.filepath,"/output"))
  raster::writeRaster(xx$flood.grid, filename=exportnames, bylayer=TRUE,format="GTiff")
  print(paste0("-- Flood innundation rasters written to: ",getwd()," --"))
  setwd(basedir)
  if(user.output.choice=="GIS_O") {stop()}
}
#/////////////////////////////////////
#/////////////////////////////////////

# Pull Addresses
v = velox::velox(xx$flood.grid)
xx$address.velox = v$extract_points(sf::st_transform(xx$address.point, raster::crs(xx$flood.grid)))
colnames(xx$address.velox) <- names(xx$flood.grid)
xx$address.point <- do.call(cbind, list(xx$address.point, xx$address.velox))

xx$address.point[names(xx$flood.grid)[1]] <- NA_real_              # testing
xx$address.point[names(xx$flood.grid)[4]] <- NA_real_              # testing

processskipflag <- FALSE
if(all(xx$address.point$ffreq==0,na.rm = TRUE) | user.output.choice=="basedata") {
  if(all(xx$address.point$ffreq==0,na.rm = TRUE)) { 
    print("-- Alert: No addresses forecasted to be impacted, defaulting to base data viewer --") 
    user.output.choice="basedata"
  }
  
  # Build basedata titleblock
  TitleBlock <- magick::image_read(paste0(basedir,"/data/misc/EmptyTitleBlock.png"))
  if(user.aoi.source=="string") {
    TitleBlock <- magick::image_annotate(TitleBlock, user.aoi.string, font = 'Georgia', location = "+398+302", size = 20)
  } else {
    TitleBlock <- magick::image_annotate(TitleBlock, paste(user.aoi.source,user.aoi.string), font = 'Georgia', location = "+398+302", size = 20)
  }
  if(user.forecast.source=="USER_DIS" || user.forecast.source=="USER_STAGE") {
    if(user.forecast.source=="USER_DIS") { 
      TitleBlock <- magick::image_annotate(TitleBlock, "User Provided Discharge", font = 'Georgia', location = "+305+328", size = 20)
    } else { TitleBlock <- magick::image_annotate(TitleBlock, "User Provided Stage", font = 'Georgia', location = "+305+328", size = 20) }
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
  TitleBlock <- magick::image_annotate(TitleBlock, format(file.mtime(xx$aoi.path),format='%m/%d/%Y'), font = 'Georgia', location = "+418+383", size = 20) %>% 
    magick::image_annotate(user.address.source, font = 'Georgia', location = "+301+410", size = 20) %>%
    magick::image_annotate(user.road.source, font = 'Georgia', location = "+627+410", size = 20) %>% 
    imager::magick2cimg(alpha = "rm") %>%
    imager::autocrop("white") %>%
    imager::cimg2magick(rotate = T) %>%
    magick::image_flop()
  
  # Set process skip flag
  processskipflag <- TRUE
}

if(!processskipflag) {
  # Pull roads
  xx$road.depth = v$extract_points(sf::st_transform(suppressWarnings(sf::st_cast(xx$road.line,"POINT")), raster::crs(xx$flood.grid)))
  colnames(xx$road.depth) <- names(xx$flood.grid)
  xx$road.point =  do.call(cbind, list(suppressWarnings(sf::st_cast(xx$road.line,"POINT")), xx$road.depth))
  xx$address.mapbook = suppressWarnings(suppressMessages(sf::st_intersection(xx$mapindex.poly, xx$address.point)))
  # sf::st_geometry(xx$address.mapbook) = NULL  
  xx$road.point.mapbook = suppressWarnings(suppressMessages(sf::st_intersection(xx$mapindex.poly, xx$road.point)))
  # sf::st_geometry(xx$road.points.mapbook) = NULL  
  
  #/////////////////////////////////////
  # summary tables
  #/////////////////////////////////////
  print("-- Summerizing results --")
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
  
  # Map index prep
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
  
  print("-- Generating cartography objects --")
  
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
  # plot(sbuxmsft.ts, col="blue", lwd=2, ylab="Adjusted close", main="Monthly closing price of SBUX") 
  # ff <- xx$nwis.discharge %>%
  #   group_by(site_no) %>%
  #   group_split()
  # 
  # 
  # ff[[2]]
  # TS3 <- timeSeries::timeSeries(ff[[2]]$X_00060_00000, as.POSIXct(ff[[2]]$dateTime, tz = "UTC"), units = "RAND")
  # 
  # 
  # TS2 <- timeSeries::timeSeries(xx$nwis.discharge)
  # sbuxmsft.ts = cbind(TS, TS1, TS3) 
  # sbuxmsft.ts
  # 
  # 
  # as.POSIXct(xx$nwm.discharge.dateTime, tz = "UTC")
  # 
  # tst <- as.ts(data=xx$floodsum1.addimpact,as.POSIXct(xx$nwm.discharge.dateTime, tz = "UTC")) 
  # tst[1:5] sbux.ts
  # plot(sbux.ts, col="blue", lwd=2, ylab="Adjusted close",main="Monthly closing price of SBUX") 
  # 
  # 
  # 
  # xx$chartdata <- data.frame(xx$floodsum1.x,xx$floodsum1.nwmts,xx$floodsum1.area,xx$floodsum1.addimpact,xx$nwm.discharge.dateTimeZone)
  # xx$nwis.discharge <- readNWISuv(xx$gage.point$site_no, "00060", as.Date(firstdate)-2, as.Date(lastdate)+1, tz="UTC")
  # xx$nwis.stage <- readNWISuv(xx$gage.point$site_no, "00065", as.Date(firstdate)-2, as.Date(lastdate)+1, tz="UTC")
  # 
  # test <- data.frame(xx$nwm.discharge.dateTimeZone,xx$floodsum1.addimpact)
  # plot(test, plot.type="s")
  # 
  # xts(xx$chartdata, order.by=as.Date(xx$chartdata$xx.nwm.discharge.dateTimeZone))
  # 
  # par(mfrow=c(1, 1))
  # plot(xx$chartdata, plot.type="s")
  
  # Create titleblock
  TitleBlock <- magick::image_read(paste0(basedir,"/data/misc/EmptyTitleBlock.png"))
  if(user.aoi.source=="string") {
    TitleBlock <- magick::image_annotate(TitleBlock, user.aoi.string, font = 'Georgia', location = "+398+302", size = 20)
  } else {
    TitleBlock <- magick::image_annotate(TitleBlock, paste(user.aoi.source,user.aoi.string), font = 'Georgia', location = "+398+302", size = 20)
  }
  if(user.forecast.source=="USER_DIS" || user.forecast.source=="USER_STAGE") {
    if(user.forecast.source=="USER_DIS") { 
      TitleBlock <- magick::image_annotate(TitleBlock, "User Provided Discharge", font = 'Georgia', location = "+305+328", size = 20)
    } else { TitleBlock <- magick::image_annotate(TitleBlock, "User Provided Stage", font = 'Georgia', location = "+305+328", size = 20) }
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
  TitleBlock <- magick::image_annotate(TitleBlock, format(file.mtime(xx$aoi.path),format='%Y-%m-%d'), font = 'Georgia', location = "+418+383", size = 20) %>% 
    magick::image_annotate(user.address.source, font = 'Georgia', location = "+301+410", size = 20) %>%
    magick::image_annotate(user.road.source, font = 'Georgia', location = "+627+410", size = 20) %>% 
    imager::magick2cimg(alpha = "rm") %>%
    imager::autocrop("white") %>%
    imager::cimg2magick(rotate = T) %>%
    magick::image_flop() 
  # vis params
  xx$flood.grid.val <- as.numeric(c(1:10))
  xx$flood.grid.pal <- colorNumeric(c("#deebf7", "#9ecae1", "#3182bd"), xx$flood.grid.val, na.color = "transparent")
  xx$flood.grid.ffreqval <- as.numeric(c(1:length(xx$nwm.timestep)))
  xx$flood.grid.ffreqpal <- colorNumeric(c("#deebf7", "#9ecae1", "#3182bd"), xx$flood.grid.ffreqval, na.color = "transparent")
  xx$mapindex.poly.WHHval <- as.numeric(c(0:max(xx$mapindex.poly$'Wet House Hours', na.rm = TRUE)))
  xx$mapindex.poly.WHHpal <- colorBin("Reds", xx$mapindex.poly$'Wet House Hours', 4, pretty = TRUE)
  xx$mapindex.poly.IHCval <- as.numeric(c(0:max(xx$mapindex.poly$'Uniquely Impacted Addresses', na.rm = TRUE)))
  xx$mapindex.poly.IHCpal <- colorBin("Reds", xx$mapindex.poly$'Uniquely Impacted Addresses', 4, pretty = TRUE)
  xx$mapindex.poly.MDval <- as.numeric(c(0:max(xx$mapindex.poly$'Maximum Impact Depth', na.rm = TRUE)))
  xx$mapindex.poly.MDpal <- colorBin("Reds", xx$mapindex.poly$'Maximum Impact Depth', 4, pretty = TRUE)
} # End of impact processing

# Index view lock
xx$aoi.bb.view <- suppressMessages(suppressWarnings(sf::st_buffer(
  sf::st_combine(xx$aoi.shp),
  0.08,
  nQuadSegs = 1,
  endCapStyle = "SQUARE",
  joinStyle = "ROUND",
  mitreLimit = 1
)))

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

# This app deployment format is less than ideal, but the other alternative I was able to implement was a disgusting file-replace using vbs like so... 
#
# strText = Replace(strText, "##useForecastFile", Replace(userFlowFile.Value, "\", "/"))
# strText = Replace(strText, "##UserOutputT", UserOutputType.Value)
# Set objFile = objFSO.OpenTextFile(shinyPathFilename, ForWriting)
# objFile.WriteLine strText
# objFile.Close
# dim shinyUIPathFilename
# shinyUIPathFilename = sCurPath & "\shiny\ui.r"
# dim shinyServerPathFilename
# shinyServerPathFilename = sCurPath & "\shiny\server.r"
# dim oldUIpath
# dim oldServerPath
# set newUI = CreateObject("Scripting.FileSystemObject")
# set newServer = CreateObject("Scripting.FileSystemObject")
# If UserOutputType.Value="BVCARTPos" Then
#   oldUIpath = sCurPath & "\shiny\ui_bview.r"
#   oldServerPath =  sCurPath & "\shiny\server_BVCARTPos.r"
#   newUI.CopyFile oldUIpath, shinyUIPathFilename
#   newServer.CopyFile oldServerPath, shinyServerPathFilename
# ElseIf UserOutputType.Value="BVOSMBW" Then
#   ...
# End If
# dim RGTXT 
# Set RGTXT = CreateObject("Scripting.FileSystemObject") 
# RGTXT.DeleteFile(shinyRPathFilename) 
# strComputer = "."
# Set objWMIService = GetObject("winmgmts:\\" & strComputer & "\root\cimv2")
# dim shinyPath2
# shinyPath2 = sCurPath & "\shiny"
# Set colFiles = objWMIService.ExecQuery _
#   ("ASSOCIATORS OF {Win32_Directory.Name='"& shinyPath2 &"'} Where " _
#     & "ResultClass = CIM_DataFile")
# For Each objFile In colFiles
#   strExtension = objFile.Extension
#   strExtension = Replace(strExtension, "txt", "r")
#   strNewName = objFile.Drive & objFile.Path & objFile.FileName & "." & strExtension
#   errResult = objFile.Rename(strNewName)
# Next
# dim shinyRAPPPath
# shinyRAPPPath = sCurPath & "\shiny\runShinyApp.r"
# 

basedataUI <- dashboardPage(
  skin = "blue",
  dashboardHeader(title="FOSSFlood V 1.1",tags$li(class="dropdown",actionButton("home","Home"),actionButton("print", "Print"))),
  dashboardSidebar(
    sidebarMenu(id = "sidebar",
                h3("Map controls"),
                pickerInput(inputId='userbasemap',label='Select a base map:',choices=basemapchoices,selected=basemapchoices[2]$`Minimal maps`[2]),
                shinyWidgets::switchInput(inputId="indexColor",label="Grid color",onLabel="Black",offLabel="White",value=TRUE)
    )
  ),
  dashboardBody(
    tags$head(tags$script('
                          // Define function to set height of "map" and "map_container"
                          setHeight = function() {
                          var window_height = $(window).height();
                          var header_height = $(".main-header").height();
                          
                          var boxHeight = window_height - header_height - 80;
                          
                          $("#map").height(boxHeight);
                          $("#summary_map").height(boxHeight - 20);
                          };
                          
                          // Set input$box_height when the connection is established
                          $(document).on("shiny:connected", function(event) {
                          setHeight();
                          });
                          
                          // Refresh the box height on every window resize event    
                          $(window).on("resize", function(){
                          setHeight();
                          });
                          '),
              tags$style(
                type="text/css",
                "#titleblock img {max-width: 100%; width: 100%; height: auto;}"
              )),
    fluidRow(
      column(width = 8,
             box(width = NULL, solidHeader = FALSE, status = "warning",leafletOutput("basemap")),
             tabBox(title = tagList(shiny::icon("gear"), "Summary Tables"), width = NULL,
                    tabPanel(title = tagList(shiny::icon("th"), "Map index"),
                             DT::dataTableOutput("mapindextable")),
                    tabPanel(title = tagList(shiny::icon("home"), "Addresses"),
                             DT::dataTableOutput("addressesmapbook")),
                    tabPanel(title = tagList(shiny::icon("road"), "Roads"),
                             DT::dataTableOutput("roadsmapbook")),
                    tabPanel(title = tagList(shiny::icon("water"), "Streams"),
                             DT::dataTableOutput("streamsmapbook"))
             )
             
      ),
      column(width = 4,
             box(width = NULL, status = "primary",imageOutput("titleblock")),
             box(title = "Map controls", width = NULL, status = "primary","Map controls"),
             box(title = "Warnings/Credits", width = NULL, status = "primary",DISCLAIMERString)
      )
      
    )
    )
)
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
  
  # InforBoxes ========================================
  output$selected_var <- renderText({
    paste0("t_",gsub(":",".",gsub("-",".",gsub(" ",".",as.character(as.POSIXlt(input$SwiperTimestep, origin="1970-01-01"))))))
  })
  
  # Maps ========================================
  summary_map_reactive <- reactive({
    gridcolor <- "Black"
    sumvisfill <- T
    
    leaflet() %>%
      addProviderTiles("Stamen.TonerBackground", group = "Base Map") %>%
      addPolygons(data=xx$aoi.shp, color=gridcolor, weight=1,group="Borders") %>%
      addPolylines(data=xx$flow.line,weight=1,group="Flowlines") %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),group="Roads") %>%
      addRasterImage(x=xx$flood.grid$ffreq,colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
      addPolylines(data=xx$road.poly[[length(xx$nwm.timestep)+1]], color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,group="Road Impacts") %>%
      addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1,group="Index") %>%
      addLabelOnlyMarkers(data=subset(xx$mapindex.point, `Wet House Hours`>0), label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T), group="Index Lables") %>%
      addAwesomeMarkers(data=subset(xx$address.point, ffreq>0),
                        icon=xx$address.point.standard, 
                        clusterOptions=markerClusterOptions(),
                        group="Addresses") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Road Impacts","Addresses","Borders","Index","Index Lables"),
        options = layersControlOptions(collapsed = FALSE)
      ) %>% 
      hideGroup("Roads") %>% 
      hideGroup("Addresses")  
  })
  output$summary_map <- renderLeaflet({
    summary_map_reactive()
  })
  observe({
    gridcolor <- ifelse(input$indexColor==TRUE, "Black", "White")
    
    if(input$summaryviewlayer=="Individual Points") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor, weight=1,group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,group="Index") %>% 
        showGroup("Addresses")
    } else if(input$summaryviewlayer=="Wet House Hours") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor, weight=1,group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1.3,group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor, weight=1,group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Maximum Impact Depth") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor, weight=1,group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,group="Index") %>% 
        hideGroup("Addresses")
    }
    
  })
  
  # Tables ========================================
  output$summarychart <- renderPlot({
    timeSeries::plot(xx$chart.summary, plot.type="s", col="blue", lwd=2, ylab="Number", main="Summary") 
  })
  
  # Tables ========================================
  output$mapindextable<-DT::renderDataTable({
    DT::datatable(as.data.frame(xx$mapindex.point[which(xx$mapindex.point$`Wet House Hours`>0),]) %>% select(`Index Label`,`Wet House Hours`,`Uniquely Impacted Addresses`,`Maximum Impact Depth`), 
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  observe({
    if (length(input$mapindextable_rows_selected)) {
      leafletProxy("summary_map") %>%
        flyTo(xx$mapindex.point[which(xx$mapindex.point$`Uniquely Impacted Addresses`>0),][input$mapindextable_rows_selected,]$geometry[[1]][1],
              xx$mapindex.point[which(xx$mapindex.point$`Uniquely Impacted Addresses`>0),][input$mapindextable_rows_selected,]$geometry[[1]][2],17)
    }
  })
  output$addressimpacts<-DT::renderDataTable({
    DT::datatable(as.data.frame(xx$address.mapbook[which(xx$address.mapbook$ffreq>0),]) %>% select(address,`Index Label`,dataset), 
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  observe({
    if (length(input$addressimpacts_rows_selected)) {
      leafletProxy("summary_map") %>%
        flyTo(xx$address.mapbook[which(xx$address.mapbook$ffreq>0),][input$addressimpacts_rows_selected,]$geometry[[1]][1],
              xx$address.mapbook[which(xx$address.mapbook$ffreq>0),][input$addressimpacts_rows_selected,]$geometry[[1]][2],17)
    }
  })
  # There are better ways to do this...
  output$roadimpacts<-DT::renderDataTable({
    DT::datatable(as.data.frame(xx$road.point.mapbook[which(xx$road.point.mapbook$ffreq>0),]) %>% select(`Index Label`,name), 
                  editable = FALSE, class = 'compact hover stripe order-column row-border stripe', options = list(scrollX = TRUE), selection = "single")
  })
  observe({
    if (length(input$roadimpacts_rows_selected)) {
      leafletProxy("summary_map") %>%
        flyTo(xx$road.point.mapbook[which(xx$road.point.mapbook$ffreq>0),][input$roadimpacts_rows_selected,]$geometry[[1]][1],
              xx$road.point.mapbook[which(xx$road.point.mapbook$ffreq>0),][input$roadimpacts_rows_selected,]$geometry[[1]][2],17)
    }
  })
  
  # Events ========================================
  observeEvent(input$home, {
    updateTabItems(session, "sidebar", "home")
  })
  # Print the map to the working directory
  observeEvent(input$print, {
    # Take a screenshot of the map
    mapshot(user_created_map(), file=paste0(getwd(), '/exported_map.png'))
    shinyalert(title = "Print complete, see /shiny!", type = "success")
  })
  
  session$onSessionEnded(stopApp)
}

##///////////////////////////////////////////////////////////////////////////////////////////
##///////////////////////////////////////////////////////////////////////////////////////////
##///////////////////////////////////////////////////////////////////////////////////////////
impactUI <- dashboardPage(
  skin = "blue",
  dashboardHeader(title="FOSSFlood V 1.1",tags$li(class="dropdown",actionButton("relaunch","Relaunch UI"),actionButton("save","Save raster data"),actionButton("print", "Print"))),
  dashboardSidebar(
    sidebarMenu(id = "sidebar",
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
      tabItem(tabName = "Summary",
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
      tabItem(tabName = "Swiper",
              tags$head(
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
      addLabelOnlyMarkers(data=subset(xx$mapindex.point, `Wet House Hours`>0), label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
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
      addPolylines(data=xx$road.poly[[1]], color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),options=pathOptions(pane="road"),group="Roads") %>%
      addAwesomeMarkers(data=subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),
                        icon=xx$address.point.standard, 
                        options=pathOptions(pane="add"),
                        clusterOptions=markerClusterOptions(),
                        group="Addresses") %>%
      addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),], fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addLabelOnlyMarkers(data=subset(xx$mapindex.point,`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`),label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Road Impacts","Addresses","Borders","Index","Index Labels"),
        options = layersControlOptions(collapsed = FALSE)
      ) %>% 
      addLegend("bottomright", pal = xx$vis_roadpal, values = xx$road.line$type, group = "Roads") %>%
      addLegend("bottomright", pal = xx$vis_roadpal, values = xx$road.line$type, group = "flow") %>%
      hideGroup("Roads") %>%
      hideGroup("Addresses") 
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
      addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),][1,],fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addLabelOnlyMarkers(data=subset(xx$mapindex.point,`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`),label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Road Impacts","Addresses","Borders","Index","Index Labels"),
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
  mapshot_summary_map <- reactive({
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(input$summaryviewlayer=="Individual Points") {
      m = leaflet() %>%
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
        addRasterImage(x=xx$flood.grid$ffreq,colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        # addMarkers(gages)
        addPolylines(data=xx$road.poly[[length(xx$nwm.timestep)+1]], color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
        addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),options=pathOptions(pane="road"),group="Roads") %>%
        addAwesomeMarkers(data=subset(xx$address.point, ffreq>0),
                          icon=xx$address.point.standard, 
                          options=pathOptions(pane="add"),
                          clusterOptions=markerClusterOptions(),
                          group="Addresses") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        addLabelOnlyMarkers(data=subset(xx$mapindex.point, `Wet House Hours`>0), label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
        addLegend("bottomright", pal = xx$vis_roadpal, values = xx$road.line$type, group = "Roads") %>%
        addLegend("bottomright", pal = xx$vis_roadpal, values = xx$road.line$type, group = "flow") 
    } else if(input$summaryviewlayer=="Wet House Hours") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Maximum Impact Depth") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    }
    m
  })
  mapshot_index_map <- reactive({
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    m = summary_map_reactive()
    m
  })
  mapshot_swiper_map <- reactive({
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    m = summary_map_reactive()
    m
  })
  
  
  # Tables ========================================
  output$summarychart <- renderDygraph({
    dygraphs::dygraph(data = xx$chart.summary,
                      main = paste0("Forecast impacts for ", user.aoi.string),
                      xlab = "Date/Time in UTC") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = TRUE) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      dyLegend(width = 600)
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
      dyLegend(width = 600)
  })
  
  # webshot for charts
  summararychart_image <- function() {
    chart <- dygraphs::dygraph(data = xx$chart.summary,
                               main = paste0("Forecast impacts for ", user.aoi.string),
                               xlab = "Date/Time in UTC") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = FALSE) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      dyLegend(width = 600)
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1080
    height <- 610
    webshot::webshot("temp.html", file = "Rplot.png",cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height )
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
      dyLegend(width = 600)
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1080
    height <- 610
    webshot::webshot("temp.html", file = "Rplot.png",cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height )
  }
  swiperchart_print_image <- function(printTimestep) {
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
      dyLegend(width = 600)
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1080
    height <- 610
    webshot::webshot("temp.html", file = "Rplot.png",cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height )
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
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(input$summaryviewlayer=="Individual Points") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addProviderTiles(input$userbasemap,group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        showGroup("Addresses")
    } else if(input$summaryviewlayer=="Wet House Hours") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Wet House Hours", pal = xx$mapindex.poly.WHHpal, values = subset(xx$mapindex.poly, `Wet House Hours`>0)$`Wet House Hours`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Count of Impacted Addresses", pal = xx$mapindex.poly.IHCpal, values = subset(xx$mapindex.poly, `Uniquely Impacted Addresses`>0)$`Uniquely Impacted Addresses`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Maximum Impact Depth") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Maximum Impact Depth", pal = xx$mapindex.poly.MDpal, values = subset(xx$mapindex.poly, `Maximum Impact Depth`>0)$`Maximum Impact Depth`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    }
  })
  observeEvent(input$Indexviewlayer, {
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(input$Indexviewlayer=="Individual Points") {
      leafletProxy("index_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap,group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        showGroup("Addresses")
    } else if(input$Indexviewlayer=="Current index") {
      leafletProxy("index_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~"Red",fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Indexviewlayer")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("index_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$Indexviewlayer=="Maximum Impact Depth") {
      leafletProxy("index_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
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
  
  observeEvent(input$gridColor | input$aoiColor, {
    shiny::req(input$summaryviewlayer)
    shiny::req(input$index_symbology)
    shiny::req(input$userbasemap)
    shiny::req(input$timeslider)
    shiny::req(input$uislidertime)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    
    if(input$summaryviewlayer=="Individual Points") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap,group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        showGroup("Addresses")
    } else if(input$summaryviewlayer=="Wet House Hours") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Maximum Impact Depth") {
      leafletProxy("summary_map") %>%
        clearGroup("Base Map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addProviderTiles(input$userbasemap, group = "Base Map") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    }
    
    leafletProxy("swiper_map") %>%
      clearGroup("Base Map") %>%
      clearGroup("Borders") %>%
      clearGroup("Index") %>%
      addProviderTiles(input$userbasemap, group = "Base Map") %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
      addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
      hideGroup("Addresses")    
    
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
    
    if( nrow(subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0))>0) {
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
    shiny::req(input$timeslider)
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
        addRasterImage(x=xx$flood.grid[[match(input$timeslider, timestepseries)]],colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$road.poly[[match(input$timeslider, timestepseries)]],color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts")
      
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
        addRasterImage(x=xx$flood.grid[[match(input$timeslider, timestepseries)]],colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$road.poly[[match(input$timeslider, timestepseries)]],color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
        addAwesomeMarkers(data=subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),
                          icon=xx$address.point.standard, 
                          options=pathOptions(pane="add"),
                          clusterOptions=markerClusterOptions(),
                          group="Addresses") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),], 
                    fill=T,
                    color=gridcolor,
                    fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),
                    fillOpacity=1,
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
        clearGroup("Index Labels")
      
      values$count <- 0
      input$forwardbutton
    }
  })  
  
  observeEvent(input$relaunch, {
    session$reload()
  })
  observeEvent(input$save, {
    exportnames <- as.character(xx$nwm.discharge.dateTimeZoneZ) %>% stringr::str_sub(6, 13) %>% paste0("Z.tif")
    exportnames[user.forecast.timesteps+1] <- "ffreq.tif"
    setwd(paste0(basedir,"/AOI/",user.aoi.filepath,"/output"))
    raster::writeRaster(xx$flood.grid, filename=exportnames,overwrite=TRUE,bylayer=TRUE,format="GTiff")
    setwd(basedir)
    shinyalert(title = paste0("-- Flood innundation rasters written to: ",basedir,"/AOI/",user.aoi.filepath,"/output --"), type = "success")
  })
  observeEvent(input$print, {
    swiperchart_image(3)  # creates temp_files folder and temp.html
    cleanFiles()
    mapshot_summary_map()
    # for (i in 1:user.forecast.timesteps)
    # 
    # # Take a screenshot of the map
    # mapshot(user_created_map(), file=paste0(getwd(), '/exported_map.png'))
    # 
    # # Generate summery table images
    # ggsave(
    #   "SummeryTable.png",
    #   plot = ggplot() + 
    #     annotation_custom(tableGrob(as.matrix(as.data.frame(table(subset(as.data.frame(dataInstance$Data),select=input$variabledisplay),dnn=input$variabledisplay))), vp = viewport(width=10))) + 
    #     theme_void(),
    #   width=10,
    #   dpi = 300
    # )
    # 
    # # Generate table images
    # maxrow = 30
    # npages = ceiling(nrow(dataInstance$Data)/maxrow)
    # for (i in 1:npages) {
    #   idx = seq(1+((i-1)*maxrow), i*maxrow)
    #   if(i*maxrow >= length(dataInstance$Data)){
    #     idx = seq(1+((i-1)*maxrow), length(dataInstance$Data))
    #   }
    #   imageName = paste0("rawTable",i,".png")
    #   ggsave(
    #     imageName,
    #     plot = ggplot() + 
    #       annotation_custom(tableGrob(as.matrix(as.data.frame(dataInstance$Data[idx, ])), vp = viewport(width=10))) + 
    #       theme_void(),
    #     width=10,
    #     dpi = 300
    #   )
    # }
    # 
    # # Modify FirstPage with dates
    # MapPage <- image_read(paste0(basedir,"/data/ClusteR_MapPageBase.png"))
    # MapPageChage <- image_annotate(MapPage, paste("Generated on:", Sys.Date(), " ",format(Sys.time(), "%H:%M:%S")), font = 'Helvetica', location = "+30+1100", size = 77)
    # 
    # # and map image
    # mappaste <- image_read(paste0(getwd(),"/exported_map.png"))
    # FirstPage <- image_composite(MapPageChage, mappaste, offset = "+30+1200")
    # 
    # # and summery table (probably do this on export...)
    # sumtablarge <- load.image(paste0(getwd(),"/SummeryTable.png"))
    # imager::save.image(autocrop(sumtablarge,"white"),paste0(getwd(),"/SummeryTableSamll.png"))
    # sumtabSmall <- image_read(paste0(getwd(),"/SummeryTableSamll.png"))
    # HeadPage <- image_composite(FirstPage, sumtabSmall, offset = "+1000+1200")
    # image_write(HeadPage, path = "HeadPage.png", format = "png")
    # 
    # # remove extra images
    # file.remove('exported_map.png')
    # file.remove('SummeryTable.png')
    # file.remove('SummeryTableSamll.png')
    # file.remove('cropsum.png')
    # 
    # # Flatten images into pdf
    # all_images = list.files(getwd(), full.names = TRUE, pattern = '.png')
    # all_images_1 <- purrr::reduce(
    #   purrr::map(all_images,image_read),
    #   c
    # )
    # image_write(all_images_1 , format = "pdf", "check.pdf")
    # file.remove(all_images)
    
    # Alert to finish :)
    shinyalert(title = "Print complete, see /shiny!", type = "success")
  })
  
  cleanFiles <- function() {
    unlink(paste0(basedir,"/temp_files"), recursive=TRUE)
    unlink(paste0(basedir,"/temp.html"))
  }
  
  session$allowReconnect(TRUE)
}







#/////////////////////////////////////
# Launch shiny
#/////////////////////////////////////
rRuntime <- Sys.time() - rStarttime
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
print("also hit")
