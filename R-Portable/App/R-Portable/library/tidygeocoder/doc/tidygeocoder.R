## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(collapse = T, comment = "#>")
options(tibble.print_min = 4L, tibble.print_max = 4L)
set.seed(42)

## ----setup, warning = FALSE, message = FALSE----------------------------------
library(tibble)
library(dplyr)
library(tidygeocoder)

address_single <- tibble(singlelineaddress = c(
  "11 Wall St, NY, NY",
  "600 Peachtree Street NE, Atlanta, Georgia"
))
address_components <- tribble(
  ~street, ~cty, ~st,
  "11 Wall St", "NY", "NY",
  "600 Peachtree Street NE", "Atlanta", "GA"
)

## -----------------------------------------------------------------------------
census_s1 <- address_single %>%
  geocode(address = singlelineaddress, method = "census", verbose = TRUE)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(census_s1)

## -----------------------------------------------------------------------------
osm_s1 <- geo(
  address = address_single$singlelineaddress, method = "osm",
  lat = latitude, long = longitude
)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(osm_s1)

## -----------------------------------------------------------------------------
census_c1 <- address_components %>%
  geocode(street = street, city = cty, state = st, method = "census")

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(census_c1)

## -----------------------------------------------------------------------------
addr_comp1 <- address_components %>%
  bind_rows(
    tibble(
      cty = c("Toronto", "Tokyo"), 
      country = c("Canada", "Japan")
    )
  )

cascade1 <- addr_comp1 %>% geocode(
  street = street, state = st, city = cty,
  country = country, method = "cascade"
)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(cascade1)

## -----------------------------------------------------------------------------
census_full1 <- address_single %>% geocode(
  address = singlelineaddress,
  method = "census", full_results = TRUE, return_type = "geographies"
)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(census_full1)

## -----------------------------------------------------------------------------
salz <- geo("Salzburg, Austria", method = "osm", full_results = TRUE) %>%
  select(-licence)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(salz)

## -----------------------------------------------------------------------------
lat_longs1 <- tibble(
  latitude = c(38.895865, 43.6534817),
  longitude = c(-77.0307713, -79.3839347)
)

rev1 <- lat_longs1 %>%
  reverse_geocode(lat = latitude, long = longitude, address = addr, method = "osm")

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(rev1)

## -----------------------------------------------------------------------------
rev2 <- reverse_geo(
  lat = lat_longs1$latitude,
  long = lat_longs1$longitude,
  method = "osm",
  full_results = TRUE
)

glimpse(rev2)

## -----------------------------------------------------------------------------
# create a dataset with duplicate and NA addresses
duplicate_addrs <- address_single %>%
  bind_rows(address_single) %>%
  bind_rows(tibble(singlelineaddress = rep(NA, 3)))

duplicates_geocoded <- duplicate_addrs %>%
  geocode(singlelineaddress, verbose = TRUE)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(duplicates_geocoded)

## -----------------------------------------------------------------------------
uniqueonly1 <- duplicate_addrs %>%
  geocode(singlelineaddress, unique_only = TRUE)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(uniqueonly1)

## -----------------------------------------------------------------------------
geo_limit <- geo(
  c("Lima, Peru", "Cairo, Egypt"),
  method = "osm",
  limit = 3, full_results = TRUE
)

glimpse(geo_limit)

## -----------------------------------------------------------------------------
cairo_geo <- geo("Cairo, Egypt",
  method = "osm", full_results = TRUE,
  custom_query = list(polygon_geojson = 1), verbose = TRUE
)

glimpse(cairo_geo)

## -----------------------------------------------------------------------------
noquery1 <- geo(c("Vancouver, Canada", "Las Vegas, NV"),
  no_query = TRUE,
  method = "arcgis"
)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(noquery1)

