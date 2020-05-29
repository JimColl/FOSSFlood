library(shiny)
library(leaflet)

# Choices for drop-downs
vars <- c(
	"Is SuperZIP?" = "superzip",
	"Centile score" = "centile",
	"College education" = "college",
	"Median income" = "income",
	"Population" = "adultpop"
)

navbarPage("FOSSFlood", id = "nav",
	tabPanel("Interactive Map:",
		div(class = "outer",
			tags$head(
				includeCSS("styles.css"),
				includeScript("gomap.js")
			),

			leafletOutput("map", width="100%", height="100%"),

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
				h4(UserZipcode),
				sliderInput("fHourInt", "Forecast hour:",	min = 1,	max = 18, value = 1, animate = animationOptions(interval = 500, loop = TRUE)),
				sliderInput("gridInt", "GridZoom:", min = 1, max = 18, value = 1),
				actionButton("button", "Forward grid")
			),
			tags$div(id="cite", 'FOSSFlood V0.9 - Jim Coll and Mike Johnson')
		)
	),
	
	tabPanel("Impact Summery:",
		fluidRow(
			column(3, selectInput("Roads", "Roads", c("All states"="", structure(state.abb, names=state.name), "Washington, DC"="DC"), multiple=TRUE))
		),
		fluidRow(
			column(1, numericInput("minScore", "Min score", min=0, max=100, value=0)),
			column(1,	numericInput("maxScore", "Max score", min=0, max=100, value=100))
		),
		hr(),
		DT::dataTableOutput("ziptable")
	),
	
	conditionalPanel("false", icon("crosshair"))
)