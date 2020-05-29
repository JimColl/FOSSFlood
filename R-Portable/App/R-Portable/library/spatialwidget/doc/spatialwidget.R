## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "# "
)

library(spatialwidget)
library(sfheaders)
library(geojsonsf)


## -----------------------------------------------------------------------------
head( widget_capitals )

## -----------------------------------------------------------------------------
js <- spatialwidget::widget_point(
  data = widget_capitals
  , fill_colour = "country"
  , legend = TRUE
  )

substr( js$data, 1, 200 )

## -----------------------------------------------------------------------------
substr( js$legend, 1, 100 )

## -----------------------------------------------------------------------------
l <- widget_point(
  widget_capitals[1:2, ]
  , fill_colour = "country"
  , legend = T
  )

substr( l$data, 1, 200 )

## -----------------------------------------------------------------------------
l <- widget_line(
  widget_roads[1:2, ]
  , stroke_colour = "ROAD_NAME"
  , legend = T
  )

substr( l$data, 1, 200 )

## -----------------------------------------------------------------------------
l <- widget_polygon(
  widget_melbourne[1:2, ]
  , fill_colour = "AREASQKM16"
  , legend = F
  )

substr( l$data, 1, 200 )

## -----------------------------------------------------------------------------
feat1 <- '{"type":"Feature","properties":{"id":1},"geometry":{"type":"Point","coordinates":[0,0]}}'
feat2 <- '{"type":"Feature","properties":{"id":2},"geometry":{"type":"Point","coordinates":[1,1]}}'
geojson <- paste0('[{"type":"FeatureCollection","features":[',feat1,',',feat2,']}]')
sf <- geojsonsf::geojson_sf( geojson )
sf

## -----------------------------------------------------------------------------
geo <- geojsonsf::sf_geojson( sf )
geo

## -----------------------------------------------------------------------------
geojsonsf::sf_geojson( sf, atomise = TRUE )

## -----------------------------------------------------------------------------
geojson <- spatialwidget:::rcpp_geojson_sf(sf = widget_arcs, geometries = c("origin","destination"))
substr( geojson, 1, 500)

## -----------------------------------------------------------------------------
geojson <- spatialwidget:::rcpp_geojson( sf = widget_capitals, geometry = "geometry")
substr( geojson, 1, 300)

## -----------------------------------------------------------------------------
df <- sfheaders::sf_to_df( widget_capitals )

geojson <- spatialwidget:::rcpp_geojson_df(df = df, list(geometry = c("x","y")) )
substr( geojson, 1, 500 )

## -----------------------------------------------------------------------------
df$z <- sample(1:500, size = nrow(df), replace = TRUE )
geojson <- spatialwidget:::rcpp_geojson_dfz( df, geometries = list(geometry = c("x","y","z") ) )
substr( geojson, 1, 500 )


