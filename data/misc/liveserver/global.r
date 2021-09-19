setwd("C:/Users/User/Downloads/FOSSFlood-master")
basedir <- getwd()   # This should look something like C:/Users/.../FOSSFlood-master
user.aoi.filepath <- "zctas_03216"
user.aoi.string <- "03216"

suppressMessages(library(installr))
suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
suppressMessages(library(geosphere))
suppressMessages(library(sf))
suppressMessages(library(htmlwidgets))
suppressMessages(library(htmltools))
suppressMessages(library(webshot))
suppressMessages(library(magick))
suppressMessages(library(animation))
suppressMessages(library(imager))
suppressMessages(library(stars))
suppressMessages(library(mapview))
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
options(tigris_use_cache = FALSE)
apptitle = "FOSSFlood V 1.21"

xx = list()
TitleBlock <- magick::image_read(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/titleblock.png"))
xx$flood.grid <- raster::brick(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/floodstack.grd"))
xx$nwm.timestep <- names(xx$flood.grid)[1:length(names(xx$flood.grid)) - 1]
xx$aoi.shp = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/aoi.RData"))
xx$road.poly = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/road_poly.RData"))
xx$road.line = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/road_line.RData"))
xx$address.point = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/address.RData"))
xx$mapindex.poly = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/index_poly.RData")) 
xx$mapindex.point = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/index_point.RData"))
xx$address.mapbook = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/address_maps.RData"))
xx$road.point.mapbook = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/road_maps.RData"))
xx$flow.line = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/flow_line.RData"))
xx$chart.summary = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chartdata.RData"))
xx$floodsum.addimpact = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chartimpact.RData"))
xx$floodsum.area = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chartarea.RData"))
xx$nwm.discharge.dateTimeZoneZ = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chart_DTZ.RData"))
xx$nwm.discharge.dateTimeZoneZHR = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chart_DTZHR.RData"))
xx$nwm.discharge.dateTimeZoneLocalHR = readRDS(paste0(basedir,"/AOI/",user.aoi.filepath,"/output/server/chart_DTLHR.RData"))
user.forecast.timesteps <- length(xx$nwm.discharge.dateTimeZoneZHR)

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
xx$flood.grid.val <- as.numeric(c(0:10))
xx$flood.grid.pal <- colorNumeric(c("#deebf7", "#9ecae1", "#3182bd"), xx$flood.grid.val, na.color = "transparent")
xx$flood.grid.ffreqval <- as.numeric(c(1:length(xx$nwm.timestep)))
xx$flood.grid.ffreqpal <- colorNumeric(c("#deebf7", "#9ecae1", "#3182bd"), xx$flood.grid.ffreqval, na.color = "transparent")
xx$mapindex.poly.WHHval <- as.numeric(c(0:max(xx$mapindex.poly$'Wet House Hours', na.rm = TRUE)))
xx$mapindex.poly.WHHpal <- colorBin("Reds", xx$mapindex.poly$'Wet House Hours', 4, pretty = TRUE)
xx$mapindex.poly.IHCval <- as.numeric(c(0:max(xx$mapindex.poly$'Uniquely Impacted Addresses', na.rm = TRUE)))
xx$mapindex.poly.IHCpal <- colorBin("Reds", xx$mapindex.poly$'Uniquely Impacted Addresses', 4, pretty = TRUE)
xx$mapindex.poly.MDval <- as.numeric(c(0:max(xx$mapindex.poly$'Maximum Impact Depth', na.rm = TRUE)))
xx$mapindex.poly.MDpal <- colorBin("Reds", xx$mapindex.poly$'Maximum Impact Depth', 4, pretty = TRUE)
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