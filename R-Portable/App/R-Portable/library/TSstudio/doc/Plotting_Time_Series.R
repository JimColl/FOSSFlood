## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
fig.width=8, 
fig.height=5,
collapse = TRUE,
comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(TSstudio)

## ----message=FALSE, warning=FALSE---------------------------------------------
data(USgas)

ts_info(USgas)

ts_plot(USgas)

## -----------------------------------------------------------------------------
ts_plot(USgas,
        title = "US Monthly Natural Gas Consumption",
        Xtitle = "Time",
        Ytitle = "Billion Cubic Feet")

## -----------------------------------------------------------------------------
ts_plot(USgas,
        title = "US Monthly Natural Gas Consumption",
        Xtitle = "Time",
        Ytitle = "Billion Cubic Feet",
        slider = TRUE)

## -----------------------------------------------------------------------------
ts_plot(USgas,
        title = "US Monthly Natural Gas Consumption",
        Xtitle = "Time",
        Ytitle = "Billion Cubic Feet",
        color = "black",
        width = 3)

## -----------------------------------------------------------------------------
ts_plot(USgas,
        title = "US Monthly Natural Gas Consumption",
        Xtitle = "Time",
        Ytitle = "Billion Cubic Feet",
        dash = "dash")

## -----------------------------------------------------------------------------
ts_plot(USgas,
        title = "US Monthly Natural Gas Consumption",
        Xtitle = "Time",
        Ytitle = "Billion Cubic Feet",
        line.mode =  "lines+markers")

## -----------------------------------------------------------------------------
class(ts_plot(USgas))

## ----message=FALSE, warning=FALSE---------------------------------------------
library(plotly)

ts_plot(USgas,
        title = "US Monthly Natural Gas Consumption",
        Xtitle = "Time",
        Ytitle = "Billion Cubic Feet",
        color =  "pink",
        Xgrid = TRUE,
        Ygrid = TRUE) %>%
  layout(paper_bgcolor = "black",
         plot_bgcolor = "black",
         font = list(color = "white"),
         yaxis = list(linecolor = "#6b6b6b",
                      zerolinecolor = "#6b6b6b",
                      gridcolor= "#444444"),
         xaxis = list(linecolor = "#6b6b6b",
                      zerolinecolor = "#6b6b6b",
                      gridcolor= "#444444"))

## -----------------------------------------------------------------------------
data("Coffee_Prices")

ts_info(Coffee_Prices)

## -----------------------------------------------------------------------------
ts_plot(Coffee_Prices)

## -----------------------------------------------------------------------------
ts_plot(Coffee_Prices,
        type = "multiple")

## -----------------------------------------------------------------------------
data("Michigan_CS")

ts_info(Michigan_CS)

## -----------------------------------------------------------------------------
ts_plot(Michigan_CS)

## -----------------------------------------------------------------------------
USgas_df <- ts_to_prophet(USgas)

str(USgas_df)


ts_plot(USgas_df)

