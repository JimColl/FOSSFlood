#' Read geojson or other formats from a local file or a URL
#'
#' @export
#'
#' @param x (character) Path to a local file or a URL.
#' @param what (character) What to return. One of "list" or "sp" (for 
#' Spatial class). Default: "list". If "sp" chosen, forced to 
#' `method="local"` 
#' @template read
#' 
#' @seealso [topojson_read()], [geojson_write()]
#' 
#' @return various, depending on what's chosen in `what` parameter
#' 
#' @details Uses [file_to_geojson()] internally to give back geojson, 
#' and other helper functions when returning spatial classes.
#' 
#' This function supports various geospatial file formats from a URL, as well 
#' as local kml, shp, and geojson file formats.
#' 
#' @section File size:
#' When using `method="web"`, be aware of file sizes.
#' https://ogre.adc4gis.com that we use for this option does not document 
#' what file size is too large, but you should get an error message like 
#' "maximum file length exceeded" when that happens. `method="local"`
#' shouldn't be sensitive to file sizes.
#'
#' @examples \dontrun{
#' # From a file
#' file <- system.file("examples", "california.geojson", package = "geojsonio")
#' (out <- geojson_read_old(file))
#'
#' # From a URL
#' url <- "https://raw.githubusercontent.com/glynnbird/usstatesgeojson/master/california.geojson"
#' geojson_read_old(url, method = "local")
#' 
#' # Use as.location first if you want
#' geojson_read_old(as.location(file))
#' 
#' # use jsonlite to parse to data.frame structures where possible
#' geojson_read_old(url, method = "local", parse = TRUE)
#' 
#' # output a SpatialClass object
#' ## read kml
#' file <- system.file("examples", "norway_maple.kml", package = "geojsonio")
#' geojson_read_old(as.location(file), what = "sp")
#' ## read geojson
#' file <- system.file("examples", "california.geojson", package = "geojsonio")
#' geojson_read_old(as.location(file), what = "sp")
#' ## read geojson from a url
#' url <- "https://raw.githubusercontent.com/glynnbird/usstatesgeojson/master/california.geojson"
#' geojson_read_old(url, what = "sp")
#' ## read from a shape file
#' file <- system.file("examples", "bison.zip", package = "geojsonio")
#' dir <- tempdir()
#' unzip(file, exdir = dir)
#' shpfile <- list.files(dir, pattern = ".shp", full.names = TRUE)
#' geojson_read_old(shpfile, what = "sp")
#' 
#' x <- "https://raw.githubusercontent.com/johan/world.geo.json/master/countries.geo.json"
#' geojson_read_old(x, method = "local", what = "sp")
#' geojson_read_old(x, method = "local", what = "list")
#' 
#' utils::download.file(x, destfile = basename(x))
#' geojson_read_old(basename(x), method = "local", what = "sp")
#' 
#' # doesn't work right now
#' ## file <- system.file("examples", "feature_collection.geojson", 
#' ##   package = "geojsonio")
#' ## geojson_read_old(file, what = "sp")
#' }
geojson_read_old <- function(x, method = "web", parse = FALSE, what = "list", stringsAsFactors = FALSE, ...) {
  UseMethod("geojson_read_old")
}

#' @export
geojson_read_old.default <- function(x, method = "web", parse = FALSE, what = "list", stringsAsFactors = FALSE, ...) { 
  stop("no 'geojson_read_old' method for ", class(x), call. = FALSE)
}

#' @export
geojson_read_old.character <- function(x, method = "web", parse = FALSE, what = "list", stringsAsFactors = FALSE, ...) { 
  read_json_old(as.location(x), method, parse, what, stringsAsFactors, ...)
}

#' @export
geojson_read_old.location_ <- function(x, method = "web", parse = FALSE, what = "list", stringsAsFactors = FALSE, ...) {
  read_json_old(x, method, parse, what, stringsAsFactors, ...)
}

read_json_old <- function(x, method, parse, what, stringsAsFactors = FALSE, ...) {
  what <- match.arg(what, c("list", "sp"))
  switch(what, 
         list = file_to_geojson(x, method, output = ":memory:", parse, ...), 
         sp = file_to_sp_old(x, stringsAsFactors = stringsAsFactors, ...)
  )
}

file_to_sp_old <- function(input, stringsAsFactors = FALSE, ...) {
  fileext <- ftype(input)
  fileext <- match.arg(fileext, c("shp", "kml", "geojson", "json"))
  input <- handle_remote(input)
  switch(
    fileext, 
    kml = rgdal::readOGR(input, rgdal::ogrListLayers(input)[1], 
                         drop_unsupported_fields = TRUE, verbose = FALSE, ...),
    # kml = tosp(input, stringsAsFactors, ...),
    shp = rgdal::readOGR(input, rgdal::ogrListLayers(input), verbose = FALSE, ...),
    # shp = tosp(input, stringsAsFactors, ...),
    geojson = rgdal::readOGR(input, rgdal::ogrListLayers(input), verbose = FALSE, ...),
    # geojson = tosp(input, stringsAsFactors, ...),
    json = rgdal::readOGR(input, rgdal::ogrListLayers(input), verbose = FALSE, ...)
    # json = tosp(input, stringsAsFactors, ...)
  )
}
