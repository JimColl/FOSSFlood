library(leaflet)
library(RColorBrewer)
library(scales)
library(lattice)
library(dplyr)
library(shiny)

val = as.numeric(c(0:30))
pal = colorNumeric(c("#deebf7", "#9ecae1", "#3182bd"), val, na.color = "transparent")

function(input, output, session) {
	
	output$OverviewMap <- renderLeaflet({
		leaflet() %>%
			#addProviderTiles(providers$OpenTopoMap, group = "TonerLite") %>%
			addTiles(group = "OSM (default)") %>%
			addProviderTiles(providers$Stamen.Toner, group = "Toner") %>%
			addProviderTiles(providers$Stamen.TonerLite, group = "Toner Lite") %>%
			addPolylines(data = myFlowlines, color = "#7B9EC8", weight = 2, group = "NHDLines") %>%
			addPolygons(data = myHexGrid, color = "#000000", weight = 1, group = "girds")
			addLayersControl(
				baseGroups = c("OSM (default)", "Toner", "Toner Lite"),
				overlayGroups = c("NHDLines", "girds"),
				options = layersControlOptions(collapsed = FALSE)
			)
	})
	
	observe({
		leafletProxy("OverviewMap") %>% 
			clearImages() %>%
			addRasterImage(x = raster(FloodInnList[input$fHourInt]), colors = pal, project = FALSE, group = "FloodRas")
	})
	
	observe({
		if (input$showStreamLines) {
			leafletProxy("OverviewMap") %>% 
				showGroup(group = "NHDLines")
		} else {
			leafletProxy("OverviewMap") %>% 
				hideGroup(group = "NHDLines")
		}
	})
	
	observe({
		if (input$showFloods) {
			leafletProxy("OverviewMap") %>% 
				showGroup(group = "FloodRas")
		} else {
			leafletProxy("OverviewMap") %>% 
				hideGroup(group = "FloodRas")
		}
	})
	
	observe({
		if (input$showGrids) {
			leafletProxy("OverviewMap") %>% 
				showGroup(group = "girds")
		} else {
			leafletProxy("OverviewMap") %>% 
				hideGroup(group = "girds")
		}
	})
	
	
	session$onSessionEnded(stopApp)
}