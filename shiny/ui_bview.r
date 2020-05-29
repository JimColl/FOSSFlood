library(shiny)
library(leaflet)

navbarPage("FOSSFlood", id = "nav",
	tabPanel("Interactive Map:",
		div(class = "outer",
			tags$head(
				includeCSS("styles.css"),
				includeScript("gomap.js")
			),

			leafletOutput("OverviewMap", width="100%", height="100%"),

			absolutePanel(id = "controls",
				class = "panel panel-default", 
				fixed = TRUE,
				draggable = TRUE, 
				top = 60, 
				left = "auto", 
				right = 20, 
				bottom = "auto",
				width = 330, 
				height = "auto",
				h1("FOSSFlood V 0.9"),
				checkboxInput("showStreamLines", "Show stream lines", value = FALSE, width = NULL),
				checkboxInput("showFloods", "Show Flooding", value = FALSE, width = NULL),
				checkboxInput("showGrids", "Show Grids", value = TRUE, width = NULL),
				sliderInput("fHourInt", "Forecast hour", min = 1, max = 18, value = 1, animate = animationOptions(interval = 500, loop = TRUE))
			),
			tags$div(id="cite", 'FOSSFlood V0.9 - Jim Coll and Mike Johnson')
		)
	),
	
	conditionalPanel("false", icon("crosshair"))
)