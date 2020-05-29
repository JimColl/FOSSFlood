## ----echo=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE, 
  message = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("geojsonlint")

## ----eval=FALSE---------------------------------------------------------------
#  remotes::install_github("ropensci/geojsonlint")

## -----------------------------------------------------------------------------
library("geojsonlint")

## -----------------------------------------------------------------------------
geojson_hint(x = '{"type": "Point", "coordinates": [-100, 80]}')

## -----------------------------------------------------------------------------
geojson_validate(x = '{"type": "Point", "coordinates": [-100, 80]}')

## -----------------------------------------------------------------------------
geojson_hint('{"type": "FooBar"}')

## -----------------------------------------------------------------------------
geojson_validate('{ "type": "FeatureCollection" }')

## -----------------------------------------------------------------------------
geojson_hint('{"type": "FooBar"}', inform = TRUE)

## -----------------------------------------------------------------------------
geojson_validate('{ "type": "FeatureCollection" }', inform = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  geojson_hint('{"type": "FooBar"}', error = TRUE)
#  #> Error: Line 1
#  #>    - The type FooBar is unknown

## ----eval=FALSE---------------------------------------------------------------
#  geojson_validate('{ "type": "FeatureCollection" }', error = TRUE)
#  #> Error: 1 error validating json:
#  #> 	- data: no (or more than one) schemas match

