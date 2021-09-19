function(input, output, session) {
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
      addLabelOnlyMarkers(data=subset(xx$mapindex.point, `Wet House Hours`>0), label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = T,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
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
      # addPolylines(data=xx$road.poly[[1]], color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),options=pathOptions(pane="road"),group="Roads") %>%
      # addAwesomeMarkers(data=subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),
      #                   icon=xx$address.point.standard,
      #                   options=pathOptions(pane="add"),
      #                   clusterOptions=markerClusterOptions(),
      #                   group="Addresses") %>%
      # addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),], fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      # addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`)), label = subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`))$`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") #%>%
      # addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`)), labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") #%>%
      addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
      addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
      # addLegend("bottomleft", title = "Flood Depth",pal = xx$flood.grid.pal, values = xx$flood.grid.val, group = "Flooding") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Road Impacts","Addresses","Borders","Index","Index Labels"),
        options = layersControlOptions(collapsed = FALSE)
      )
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
      # addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      # addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),][1,],fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      # addLabelOnlyMarkers(data=subset(xx$mapindex.point,`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[1]))>0),]$`Index Label`),label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLayersControl(
        overlayGroups = c("Flowlines","Roads","Flooding","Borders","Index"),
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
  mapshot_summary_map <- function() {
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    mybbox <- suppressWarnings(suppressMessages(AOI::bbox_coords(sf::st_buffer(x=xx$aoi.shp,dist=0.00005,nQuadSegs=1,endCapStyle="SQUARE",joinStyle="ROUND",mitreLimit=1,singleSide=FALSE))))
    
    m = leaflet(height=1000, width=1000) %>%
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
      addPolylines(data=xx$flow.line,weight = ~streamorde/2,color = "Blue",options=pathOptions(pane="flow"),group="Flowlines") %>%
      addRasterImage(x=xx$flood.grid$ffreq,colors=xx$flood.grid.ffreqpal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
      addPolylines(data=xx$road.line,color= ~xx$vis_roadpalFull(type_full),weight = 0.8,options=pathOptions(pane="road"),group="Roads") %>%
      addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
      addLabelOnlyMarkers(data=subset(xx$mapindex.point, `Wet House Hours`>0), label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = T,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
      addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
      addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
      addLegend("bottomleft", title = "Flood Frequency",pal = xx$flood.grid.ffreqpal, values = xx$flood.grid.ffreqval, group = "Flooding") %>%
      addLegend("bottomleft", title = "Wet House Hours", pal = xx$mapindex.poly.WHHpal, values = subset(xx$mapindex.poly, `Wet House Hours`>0)$`Wet House Hours`, group = "Index",layerId = "Index_legend") %>%
      fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
    
    return(m)
  }
  mapshot_index_map <- function(timestep,index,bboxobj){
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(is.null(index)) {
      m = leaflet(height=1000, width=1000) %>%
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
        addRasterImage(x=xx$flood.grid[[timestep]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        fitBounds(bboxobj$xmax,bboxobj$ymin,bboxobj$xmin,bboxobj$ymax)
      
    } else {
      
      m = leaflet(height=1000, width=1000) %>%
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
        addRasterImage(x=xx$flood.grid[[timestep]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),][xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`==index,],
                    fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`)),
                            label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
        fitBounds(bboxobj$xmax,bboxobj$ymin,bboxobj$xmin,bboxobj$ymax)
    }
    
    return(m)
  }
  mapshot_swiper_map <- function(timestep,index){
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(is.null(index)) {
      mybbox <- AOI::bbox_coords(sf::st_buffer(
        x=xx$aoi.shp,
        dist=0.0002,
        nQuadSegs = 1,
        endCapStyle = "SQUARE",
        joinStyle = "ROUND",
        mitreLimit = 1,
        singleSide = FALSE
      ))
      m = leaflet(height=1000, width=1000) %>%
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
        addRasterImage(x=xx$flood.grid[[timestep]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),options=pathOptions(pane="road"),group="Roads") %>%
        addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
        addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
        addLegend("bottomleft", title = "Flood Depth",pal = xx$flood.grid.pal, values = xx$flood.grid.val, group = "Flooding") %>%
        fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
      
    } else {
      mybbox <- AOI::bbox_coords(sf::st_buffer(
        x=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),][xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`==index,],
        dist=0.0002,
        nQuadSegs = 1,
        endCapStyle = "SQUARE",
        joinStyle = "ROUND",
        mitreLimit = 1,
        singleSide = FALSE
      ))
      m = leaflet(height=1000, width=1000) %>%
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
        addRasterImage(x=xx$flood.grid[[timestep]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$flow.line,weight=1,options=pathOptions(pane="flow"),group="Flowlines") %>%
        addPolylines(data=xx$road.line,color= ~xx$vis_roadpal(type),options=pathOptions(pane="road"),group="Roads") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        # addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),][xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`==index,],
        # fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addAwesomeMarkers(data=subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),
                          icon=xx$address.point.standard, 
                          options=pathOptions(pane="add"),
                          # clusterOptions=markerClusterOptions(),
                          group="Addresses") %>%
        addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[timestep]))>0),]$`Index Label`)),
                            label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index Labels") %>%
        addLegend("bottomright",title = "Road Class", pal = xx$vis_roadpalFull, values = xx$road.line$type_full, group = "Roads") %>%
        addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines", group = "Flow") %>%
        addLegend("bottomleft", title = "Flood Depth",pal = xx$flood.grid.pal, values = xx$flood.grid.val, group = "Flooding") %>%
        fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
    }
    
    return(m)
  }
  
  # Tables ========================================
  output$summarychart <- renderDygraph({
    dygraphs::dygraph(data = xx$chart.summary,
                      main = paste0("Forecast impacts for ", user.aoi.string),
                      xlab = "Date/Time in UTC") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = FALSE) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      # dyLegend(width = 650) %>%
      dyCSScoolWeb()
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
      # dyLegend(width = 600)
      dyCSScoolWeb()
  })
  
  # webshot for charts
  summararychart_image <- function(workingdir,filename) {
    setwd(workingdir)
    chart <- dygraphs::dygraph(data = xx$chart.summary,
                               main = paste0("Forecast impacts for ", user.aoi.string),
                               xlab = "Date/Time in UTC") %>%
      dyOptions(useDataTimezone = TRUE) %>%
      dyAxis("y", label = "Count of Impacted Addresses", independentTicks = TRUE) %>%
      dyAxis("y2", label = "Inundated Area (Sq. miles)", independentTicks = FALSE, labelWidth = 120) %>%
      dySeries("Square miles of Inundated Land", axis = 'y2') %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3)) %>%
      dyLegend(width = 650) %>%
      dyCSScoolMin()
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1100
    height <- 400
    webshot::webshot("temp.html", file = filename,cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height)
    unlink(paste0(workingdir,"/temp_files"), recursive=TRUE)
    unlink(paste0(workingdir,"/temp.html"))
    setwd(basedir)
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
      dyLegend(width = 650) %>%
      dyCSScoolMin()
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1100
    height <- 400
    webshot::webshot("temp.html", file = filename,cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height )
  }
  swiperchart_print_image <- function(printTimestep,workingdir,filename) {
    setwd(workingdir)
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
      dyLegend(width = 650) %>%
      dyCSScoolMin()
    htmlwidgets::saveWidget(chart, "temp.html", selfcontained = FALSE)
    width<- 1100
    height <- 400
    webshot::webshot("temp.html", file = filename,cliprect = c(10,30,width+50,height+50) ,vwidth = width, vheight = height )
    unlink(paste0(workingdir,"/temp_files"), recursive=TRUE)
    unlink(paste0(workingdir,"/temp.html"))
    setwd(basedir)
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
    shiny::req(input$gridColor)
    shiny::req(input$aoiColor)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(input$summaryviewlayer=="Individual Points") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        showGroup("Addresses")
    } else if(input$summaryviewlayer=="Wet House Hours") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Wet House Hours", pal = xx$mapindex.poly.WHHpal, values = subset(xx$mapindex.poly, `Wet House Hours`>0)$`Wet House Hours`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Count of Impacted Addresses", pal = xx$mapindex.poly.IHCpal, values = subset(xx$mapindex.poly, `Uniquely Impacted Addresses`>0)$`Uniquely Impacted Addresses`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Maximum Impact Depth") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        removeControl("Index_legend") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addLegend("bottomleft", title = "Maximum Impact Depth", pal = xx$mapindex.poly.MDpal, values = subset(xx$mapindex.poly, `Maximum Impact Depth`>0)$`Maximum Impact Depth`, group = "Index",layerId = "Index_legend") %>%
        hideGroup("Addresses")
    }
  })
  observeEvent(input$Indexviewlayer, {
    shiny::req(input$gridColor)
    shiny::req(input$aoiColor)
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    
    if(input$Indexviewlayer=="Individual Points") {
      leafletProxy("index_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        showGroup("Addresses")
    } else if(input$Indexviewlayer=="Current index") {
      leafletProxy("index_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~"Red",fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Indexviewlayer")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("index_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    } else if(input$Indexviewlayer=="Maximum Impact Depth") {
      leafletProxy("index_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>% 
        hideGroup("Addresses")
    }
  })
  
  observeEvent((input$gridColor | input$aoiColor), {
    shiny::req(input$summaryviewlayer)
    shiny::req(input$timeslider)
    shiny::req(input$uislidertime)
    ifelse(input$aoiColor==TRUE,aoiColor <- "Black",aoiColor <- "White")
    ifelse(input$gridColor==TRUE,gridcolor <- "Black",gridcolor <- "White")
    ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
    
    if(input$summaryviewlayer=="Individual Points") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp, color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=F,color=gridcolor,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        showGroup("Addresses")
    } else if(input$summaryviewlayer=="Wet House Hours") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Count of Impacted Addresses") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.IHCpal(`Uniquely Impacted Addresses`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        hideGroup("Addresses")
    } else if(input$summaryviewlayer=="Maximum Impact Depth") {
      leafletProxy("summary_map") %>%
        clearGroup("Borders") %>%
        clearGroup("Index") %>%
        addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
        addPolygons(data=subset(xx$mapindex.poly, `Wet House Hours`>0), fill=T,color=gridcolor,fillColor=~xx$mapindex.poly.MDpal(`Maximum Impact Depth`),fillOpacity=1,weight = 1.3,options=pathOptions(pane="index"),group="Index") %>%
        hideGroup("Addresses")
    }
    leafletProxy("swiper_map") %>%
      clearGroup("Borders") %>%
      addPolygons(data=xx$aoi.shp,color=gridcolor,fillColor=aoiColor,weight=2,options=pathOptions(pane="bg"),group="Borders") %>%
      hideGroup("Addresses")
  })
  observeEvent(input$userbasemap, {
    
    leafletProxy("summary_map") %>%
      clearGroup("Base Map") %>%
      addProviderTiles(input$userbasemap, group = "Base Map")
    
    leafletProxy("swiper_map") %>%
      clearGroup("Base Map") %>%
      addProviderTiles(input$userbasemap, group = "Base Map")
    
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
    
    if( nrow(subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0))>0 ) {
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
        addRasterImage(x=xx$flood.grid[[match(input$timeslider, timestepseries)]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addLegend("bottomleft", title = "Flood Depth",pal = xx$flood.grid.pal, values = xx$flood.grid.val, group = "Flooding") 
      # addPolylines(data=xx$road.poly[[match(input$timeslider, timestepseries)]],color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts")
      
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
        addRasterImage(x=xx$flood.grid[[match(input$timeslider, timestepseries)]],colors=xx$flood.grid.pal,layerId="Flooding",group="Flooding",project=FALSE,maxBytes=20*1024*1024) %>%
        addPolylines(data=xx$road.poly[[match(input$timeslider, timestepseries)]],color="red",fill = T,fillColor = "red", fillOpacity = 0.75, weight = 1.2,options=pathOptions(pane="roadimpact"),group="Road Impacts") %>%
        addAwesomeMarkers(data=subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),
                          icon=xx$address.point.standard, 
                          options=pathOptions(pane="add"),
                          # clusterOptions=markerClusterOptions(),
                          group="Addresses") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),], 
                    fill=F,
                    color=gridcolor,
                    fillColor=~xx$mapindex.poly.WHHpal(`Wet House Hours`),
                    # fillOpacity=1,
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
        clearGroup("Index Labels") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),], fill=F,color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addPolygons(data=xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),][1,],fill=T,fillColor="red",color=gridcolor,fillOpacity=1,weight = 1,options=pathOptions(pane="index"),group="Index") %>%
        addLabelOnlyMarkers(data=subset(xx$mapindex.point,(`Index Label` %in% xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[match(input$timeslider, timestepseries)]))>0),]$`Index Label`)),label = ~`Index Label`, labelOptions = labelOptions(noHide = T, sticky = F,textOnly = T),options=pathOptions(pane="indexlabels"),group="Index")
      
      values$count <- 0
      input$forwardbutton
    }
  })
  observeEvent(input$sidebar, {
    if(input$sidebar == "Swiper"){
      shiny::req(input$uislidertime)
      ifelse(input$uislidertime=="Zulu",timestepseries <- xx$nwm.discharge.dateTimeZoneZHR,timestepseries <- xx$nwm.discharge.dateTimeZoneLocalHR)
      updateSliderTextInput(session = session,inputId="timeslider", selected = timestepseries[2])
      updateSliderTextInput(session = session,inputId="timeslider", selected = timestepseries[1])
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
  # observeEvent(input$relaunch, {
  #   session$reload()
  # })
  observeEvent(input$save, {
    exportnames <- as.character(xx$nwm.discharge.dateTimeZoneZ) %>% stringr::str_sub(6, 13) %>% paste0("Z.tif")
    exportnames[user.forecast.timesteps+1] <- "ffreq.tif"
    setwd(paste0(basedir,"/AOI/",user.aoi.filepath,"/output"))
    raster::writeRaster(xx$flood.grid, filename=exportnames,overwrite=TRUE,bylayer=TRUE,format="GTiff")
    setwd(basedir)
    shinyalert(title = paste0("-- Flood innundation rasters written to: ",basedir,"/AOI/",user.aoi.filepath,"/output --"), type = "success")
  })
  observeEvent(input$print, {
    outputdir <- paste0(basedir,"/AOI/",user.aoi.filepath,"/output")
    
    # Title Page 
    TitlePageTitleBlock <- magick::image_read(paste0(basedir,"/data/misc/EmptyFrontPageBlock.png"))
    if(user.aoi.source=="string") {
      TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, user.aoi.string, font = 'Georgia', location = "+158+328", size = 20)
    } else {
      TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, paste(user.aoi.source,user.aoi.string), font = 'Georgia', location = "+158+328", size = 20)
    }
    if(user.output.choice=="basedata") {
      TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, "Base data map books", font = 'Georgia', location = "+215+357", size = 20)
    } else {
      TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, "Flood impact mapping", font = 'Georgia', location = "+215+357", size = 20)
    }
    TitlePageTitleBlock <- magick::image_annotate(TitlePageTitleBlock, format(file.mtime(xx$aoi.path),format='%m/%d/%Y'), font = 'Georgia', location = "+418+383", size = 20) %>% 
      magick::image_annotate(user.address.source, font = 'Georgia', location = "+301+410", size = 20) %>%
      magick::image_annotate(user.road.source, font = 'Georgia', location = "+627+410", size = 20) %>%
      magick::image_trim()
    mybbox <- suppressWarnings(suppressMessages(AOI::bbox_coords(sf::st_buffer(x=xx$aoi.shp,dist=0.0001,nQuadSegs=1,endCapStyle="SQUARE",joinStyle="ROUND",mitreLimit=1,singleSide=FALSE))))
    # mybigbbox <- suppressWarnings(suppressMessages(AOI::bbox_coords(sf::st_buffer(x=xx$aoi.shp,dist=0.004,nQuadSegs=1,endCapStyle="SQUARE",joinStyle="ROUND",mitreLimit=1,singleSide=FALSE))))
    boundrydatapoly <- suppressWarnings(suppressMessages(featurefile[sf::st_buffer(xx$aoi.shp,dist=0.002,nQuadSegs=1,endCapStyle="SQUARE",joinStyle="ROUND",mitreLimit=1,singleSide=FALSE),]))
    boundrydatapoint <- suppressWarnings(suppressMessages(boundrydatapoly %>% st_centroid()))
    quiet(ifelse(user.aoi.source=="zctas",insetlablefield<-boundrydatapoint$ZCTA5CE10,boundrydatapoint<-featurefile$HUC8))
    mapshot(leaflet(height=500, width=500) %>%
              addProviderTiles("Esri.WorldStreetMap") %>%
              addPolygons(data=boundrydatapoly,color="black",fillColor="black",fillOpacity=0.08,weight=1,labelOptions=labelOptions(noHide=T,direction="bottom")) %>%
              addPolygons(data=xx$aoi.shp,color="red",fill=FALSE,weight=5) %>%
              # addLegend("bottomleft", title = "Processed Area", pal = colorFactor(palette = "red",domain = "AOI"), values = "AOI") %>%
              addLabelOnlyMarkers(data=boundrydatapoint, label =insetlablefield, labelOptions = labelOptions(noHide=T,sticky=T,textOnly=TRUE,direction="bottom")) %>%
              fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
            ,file=paste0(outputdir, '/FrontPageInset.png'))   
    mapshot(leaflet(height=1000, width=1000) %>%
              addProviderTiles("Esri.WorldImagery") %>%
              addPolygons(data=xx$aoi.shp,color="red",fill=FALSE,weight=3) %>%
              addPolylines(data=xx$flow.line,weight = ~streamorde/2,group="Flowlines") %>%
              # fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax) %>%
              addLegend("bottomright", pal = colorFactor(palette = "Blue",domain = "Stream Lines"), values = "Stream Lines") %>%
              addLegend("bottomright", title = "Processed Area",pal = colorFactor(palette = "Red",domain = "AOI"), values = "AOI") %>%
              # addScaleBar(position = "bottomleft", options = scaleBarOptions(maxWidth = 500,imperial = TRUE)) %>%
              fitBounds(mybbox$xmax,mybbox$ymin,mybbox$xmin,mybbox$ymax)
            ,file=paste0(outputdir, '/FrontPageMap.png'))   
    blankpage <- magick::image_read(paste0(basedir,'/data/misc/Empty8x11.png'))
    AOIMapInset <- magick::image_read(paste0(outputdir, '/FrontPageInset.png')) %>%
      magick::image_trim() %>%
      magick::image_scale("750x")
    frontpagemap <- magick::image_read(paste0(outputdir, '/FrontPageMap.png'))
    AOIMapImage <- magick::image_composite(blankpage, magick::image_scale(frontpagemap,"2250x"), offset = "+150+900")
    AOIMapImage <- magick::image_composite(AOIMapImage, AOIMapInset, offset = "+1650+150")
    AOIMapImage <- magick::image_composite(AOIMapImage, magick::image_scale(TitlePageTitleBlock,"1500x"), offset = "+150+150")
    unlink(paste0(outputdir,'/FrontPageInset.png'))
    unlink(paste0(outputdir,'/FrontPageMap.png'))
    magick::image_write(AOIMapImage, path = paste0(outputdir,'/00_TitlePage.png'), format = "png")
    
    # Summary page
    mapshot(mapshot_summary_map(), file=paste0(outputdir, '/summary_map.png'))   
    summarymap <- magick::image_read(paste0(outputdir,'/summary_map.png'))
    summararychart_image(outputdir,'summary_graph_table.png')  
    # summarygraph <- magick::image_read(paste0(outputdir,'/summary_graph_table.png'))
    summarygraph <- magick::image_read(paste0(outputdir,'/summary_graph_table.png')) %>%
      magick::image_trim()
    ggsave(
      filename=paste0(outputdir, '/summary_map_table.png'),
      plot = ggplot() +
        annotation_custom(tableGrob(as.matrix(as.data.frame(xx$mapindex.point[which(xx$mapindex.point$`Wet House Hours`>0),]) %>% select(`Index Label`,`Wet House Hours`,`Uniquely Impacted Addresses`,`Maximum Impact Depth`)), vp = viewport(width=10))) +
        theme_void(),width=10,units="in",dpi = 300)
    summarytable <- magick::image_read(paste0(outputdir,'/summary_map_table.png')) %>%
      magick::image_trim()
    magick::image_write(AOIMapImage, path = paste0(outputdir,'/00_TitlePage.png'), format = "png")
    titleblockimage <- magick::image_read(paste0(outputdir,'/titleblock.png'))
    page1 <- magick::image_composite(blankpage, magick::image_scale(titleblockimage,"975x"), offset = "+150+150")
    page1 <- magick::image_composite(page1,
                                     magick::image_scale(summarymap, "2250x") %>% image_border("black", "4x4"), 
                                     offset = "+150+900")
    page1 <- magick::image_composite(page1,
                                     magick::image_scale(summarygraph, "967x") %>% image_border("black", "4x4"), 
                                     offset = "+150+533")
    page1 <- magick::image_composite(page1,
                                     magick::image_scale(summarytable, "1000x") %>% image_border("black", "4x4"), 
                                     offset = "+1400+150")
    
    unlink(paste0(outputdir,'/summary_graph_table.png'))
    unlink(paste0(outputdir,'/summary_map.png'))
    unlink(paste0(outputdir,'/summary_map_table.png'))
    magick::image_write(page1, path = paste0(outputdir,'/01_summarypage.png'), format = "png")
    
    for(i in 1:length(xx$nwm.discharge.dateTimeZoneLocalHR)) {
      workingpage <- magick::image_composite(blankpage, magick::image_scale(titleblockimage,"975x"), offset = "+150+150")
      swiperchart_print_image(i,outputdir,paste0("swiper_",i,"_graph.png"))
      swiperchart <- magick::image_read(paste0(outputdir,"/swiper_",i,"_graph.png")) %>%
        magick::image_trim()
      workingpage <- magick::image_composite(workingpage, magick::image_scale(swiperchart, "967x") %>% image_border("black", "4x4"), offset = "+150+533")
      
      if(nrow(subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0))>0) {
        
        for(j in 1:length(xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),])) {
          mapshot(mapshot_index_map(i,xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),][j,]$'Index Label',mybbox), file=paste0(outputdir, '/swiper_index_',i,'_',j,'_.png'))
          pindexmap <- magick::image_read(paste0(outputdir,'/swiper_index_',i,'_',j,'_.png')) %>%
            magick::image_scale("750x")
          
          mapshot(mapshot_swiper_map(i,xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),][j,]$'Index Label'), file=paste0(outputdir, '/swiper_',i,'_',j,'_.png'))
          pswipermap <- magick::image_read(paste0(outputdir,'/swiper_',i,'_',j,'_.png')) %>%
            magick::image_trim() %>%
            magick::image_scale("x731")
          
          tg = gridExtra::tableGrob(as.matrix(
            as.data.frame(
              subset(xx$address.mapbook, eval(parse(text = names(xx$flood.grid)[i]))>0 & xx$address.mapbook$`Index Label`==xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),][j,]$'Index Label')  ) %>%
              select(address)))
          h = grid::convertHeight(sum(tg$heights), "in", TRUE)
          w = grid::convertWidth(sum(tg$widths), "in", TRUE)
          ggplot2::ggsave(paste0(outputdir, '/swiper_impact_',i,'_',j,'_.png'), tg, width=w, height=h)
          swipertable <- magick::image_read(paste0(outputdir,'/swiper_impact_',i,'_',j,'_.png')) %>%
            magick::image_trim()
          
          workingpage_new <- magick::image_composite(workingpage, magick::image_scale(pswipermap, "2250x") %>% image_border("black", "4x4"), offset = "+150+900")
          workingpage_new <- magick::image_composite(workingpage_new, pindexmap %>% image_border("black", "4x4"),offset = "+1125+150")
          workingpage_new <- magick::image_composite(workingpage_new,magick::image_scale(swipertable, "270x") %>% image_border("black", "4x4"), offset = "+2114+248")
          workingpage_new <- magick::image_annotate(workingpage_new, paste0("Impacts for: ",xx$nwm.discharge.dateTimeZoneZHR[i]), font = 'Georgia', location = "+1891+150", size = 30)
          workingpage_new <- magick::image_annotate(workingpage_new, paste0("Index: ",xx$mapindex.poly[subset(xx$address.point, eval(parse(text = names(xx$flood.grid)[i]))>0),][j,]$'Index Label'), font = 'Georgia', location = "+2180+182", size = 30)
          workingpage_new <- magick::image_annotate(workingpage_new, paste0("Impacted Addresses:"), font = 'Georgia', location = "+2094+214", size = 30)
          magick::image_write(workingpage_new, path = paste0(outputdir,'/02_impactpage_',i,'_',j,'_.png'), format = "png")
          
          unlink(paste0(outputdir,'/swiper_impact_',i,'_',j,'_.png'))
          unlink(paste0(outputdir,'/swiper_index_',i,'_',j,'_.png'))
          unlink(paste0(outputdir,'/swiper_',i,'_',j,'_.png'))
        }
        
      } else {
        # no impact in timestep
        mapshot(mapshot_index_map(i,NULL,mybbox), file=paste0(outputdir, '/swiper_',i,'_0.png'))
        pswipermap <- magick::image_read(paste0(outputdir,'/swiper_',i,'_0.png'))
        
        mapshot(mapshot_swiper_map(i,NULL), file=paste0(outputdir, '/swiper_index_',i,'_0.png'))
        pindexmap <- magick::image_read(paste0(outputdir,'/swiper_index_',i,'_0.png')) %>%
          magick::image_trim() %>%
          magick::image_scale("x731")
        
        workingpage_new <- magick::image_composite(workingpage,
                                                   magick::image_scale(pswipermap, "2250x") %>% image_border("black", "4x4"), 
                                                   offset = "+150+900")
        workingpage_new <- magick::image_composite(workingpage_new, pindexmap %>% image_border("black", "4x4"), offset = "+1125+150")
        workingpage_new <- magick::image_annotate(workingpage_new, paste0("Impacts for: ",xx$nwm.discharge.dateTimeZoneZHR[i]), font = 'Georgia', location = "+1891+150", size = 30)
        workingpage_new <- magick::image_annotate(workingpage_new, paste0("Index: NONE"), font = 'Georgia', location = "+2180+182", size = 30)
        magick::image_write(workingpage_new, path = paste0(outputdir,'/02_impactpage_',i,'_.png'), format = "png")
        
        unlink(paste0(outputdir,'/swiper_',i,'_0.png'))
        unlink(paste0(outputdir,'/swiper_index_',i,'_0.png'))
      }
      unlink(paste0(outputdir,"/swiper_",i,"_graph.png"))
    }
    
    all_images = list.files(outputdir, full.names = TRUE, pattern = '.png')
    all_images <- all_images[-length(all_images)]
    all_images_1 <- purrr::reduce(purrr::map(all_images,image_read),c)
    image_write(all_images_1 , format = "pdf", paste0(outputdir,"/FOSSFlood_Impact_Output_",user.forecast.gen[1],"_",user.forecast.gen[2],"Z.pdf"))
    file.remove(all_images)
    
    # Alert to finish :)
    shinyalert(title = paste0("See ",user.aoi.filepath,"/output"), type = "success")
  })
  
  cleanFiles <- function() {
    unlink(paste0(basedir,"/temp_files"), recursive=TRUE)
    unlink(paste0(basedir,"/temp.html"))
  }
  
  # session$allowReconnect(TRUE)
  session$onSessionEnded(function() {
    shiny::stopApp()
  })
}