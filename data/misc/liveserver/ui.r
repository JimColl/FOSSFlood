dashboardPage(
  skin = "blue",
  dashboardHeader(title=apptitle,tags$li(class="dropdown",actionButton("save","Save raster data"),actionButton("print", "Print"))),
  # dashboardHeader(title=apptitle,tags$li(class="dropdown",actionButton("relaunch","Relaunch UI"),actionButton("save","Save raster data"),actionButton("print", "Print"))),
  dashboardSidebar(
    sidebarMenu(id="sidebar",
                h3("Map controls"),
                pickerInput(inputId='userbasemap',label='Select a base map:',choices=basemapchoices,selected=basemapchoices[3]$`Just lables`[1]),
                shinyWidgets::switchInput(inputId="gridColor",label="Grid color",onLabel="Black",offLabel="White",value=TRUE),
                shinyWidgets::switchInput(inputId="aoiColor",label="AOI shading",onLabel="Black",offLabel="White",value=TRUE),
                menuItem("summarytab", tabName = "Summary", icon = icon("info-circle")),
                menuItem("swipertab", icon = icon("arrows-h"), tabName = "Swiper"),
                # add_busy_gif(src = "https://jeroen.github.io/images/banana.gif",height = 70, width = 70)
                # use_busy_spinner(spin="hollow-dots",choices = epic_spinners())
                add_busy_spinner(spin="orbit", color = "red",position='bottom-left', margins = c(40, 25))
    )
  ),
  dashboardBody(
    useShinyalert(),
    useShinyjs(),
    tabItems(
      tabItem(tabName = "Summary",value="Summary",
              tags$head(
                tags$style(type = "text/css", "#summary_map {height: calc(60vh - 80px) !important;}"),
                tags$style(type="text/css", "#summarytitleblock img {max-width: 100%; width: 100%; height: auto;}")
              ),
              fluidRow(
                column(width = 4,
                       box(width = NULL,height = 287,imageOutput("summarytitleblock")),
                       tabBox(title = tagList(shiny::icon("gear"), "Summary Tables"), width = NULL,
                              tabPanel(title = tagList(shiny::icon("th"), "Map index"),
                                       DT::dataTableOutput("mapindextable")),
                              tabPanel(title = tagList(shiny::icon("home"), "Addresses impacted"),
                                       DT::dataTableOutput("addressimpacts")),
                              tabPanel(title = tagList(shiny::icon("road"), "Roads impacted"),
                                       DT::dataTableOutput("roadimpacts"))
                       ),
                       radioGroupButtons(inputId="summaryviewlayer",
                                         label="View impacts as:",
                                         direction = "vertical",
                                         choices=c("Wet House Hours","Count of Impacted Addresses","Maximum Impact Depth","Individual Points"),individual = TRUE,
                                         checkIcon = list(yes = tags$i(class = "fa fa-circle", style = "color: steelblue"),no = tags$i(class = "fa fa-circle-o", style = "color: steelblue"))),
                       box(title = "Credits & Warnings", width = NULL, status = "primary",HTML(DISCLAIMERString))),
                column(width = 8,
                       box(title = tagList(shiny::icon("chart-line"), "Time series"), width = NULL,dygraphOutput("summarychart",height = 287)),
                       box(width = NULL, solidHeader = FALSE,status = "warning",leafletOutput("summary_map")))
              )
      ),
      tabItem(tabName = "Swiper",value="Swiper",
              tags$head(
                tags$style(type = "text/css", "#summary_map {height: calc(60vh - 80px) !important;}"),
                tags$style(type="text/css", "#swipermaintitleblock img {max-width: 100%; width: 100%; height: auto;}")
              ),
              fluidRow(
                column(width = 8,
                       box(title = "time series", width = NULL,status = "primary", 
                           dygraphOutput("dygraphswiperchart"),
                           fluidRow(
                             box(width = 10, sliderTextInput(inputId = "timeslider", grid = TRUE,force_edges = TRUE,label = "Timestep to view:",choices = xx$nwm.discharge.dateTimeZoneZHR,selected = xx$nwm.discharge.dateTimeZoneZHR[1])),
                             box(width = 2, radioGroupButtons(inputId = "uislidertime",label = "Slider time labels:",choices = c("Zulu", "Local")))
                           )
                       ),
                       box(title = "SwiperMap", width = NULL, status = "primary",leafletOutput("swiper_map")),
                       box(title = "Impacted features", width = NULL, status = "primary",tabBox(title = tagList(shiny::icon("gear"), "Summary Tables"), width = NULL,
                                                                                                tabPanel(title = tagList(shiny::icon("home"), "Addresses impacted"),
                                                                                                         DT::dataTableOutput("filteraddressimpacts")),
                                                                                                tabPanel(title = tagList(shiny::icon("road"), "Roads impacted"),
                                                                                                         DT::dataTableOutput("filterroadimpacts"))
                       ))
                ),
                column(width = 4,
                       imageOutput("swipermaintitleblock"),
                       box(title = "Map index", width = NULL, status = "primary",
                           leafletOutput("index_map"),
                           fluidRow(width = NULL,align="center",
                                    actionBttn(inputId = "backbutton", label = "Back-Previous", style = "pill", color = "danger"),
                                    actionBttn(inputId = "forwardbutton", label = "Forward-Next", style = "pill", color = "danger"),
                                    textOutput("activeIndexLabel"))
                       ),
                       box(title = "Warnings", width = NULL, status = "primary",HTML(DISCLAIMERString))
                )
              )
      )
    )
  )
)