shinyUI(fluidPage(
  theme = shinytheme("sandstone"),
  #shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
  titlePanel("SpANOVA"),
  sidebarLayout(
    sidebarPanel(
      selectInput("spem",label = "Spatial Error Modelling",
                  choices = c("Geostatistical Approach" = "geoest",
                      "Spatial Simultaneous Autoregressive Approach" = "arm"),
                       selected = "geoest"),
      uiOutput("expd"),
      tags$hr(),
      tabsetPanel(
        tabPanel("Input Variables",
                 #radioButtons("typeA","Analysis of",c("Single-factor","Factorial"),
                  #            selected = "Single-factor", inline = TRUE),
                 bsAlert("alert"), #alerta under construction
                 uiOutput("variables"),value = "variables"),
        tabPanel("Modelling Controls", uiOutput("geoUI"),
                 value = "modelCont"), id = "sidepanel")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data",
                 fluidRow(
          column(5,
                 selectInput("dataType","Data Input",
                             c("Import data from drive" = 'ifd',
                               "Pre-loaded dataset" = 'dkb')),
                 uiOutput("dataTBS")),
          column(7,dataTableOutput("table"))
        ),value = "data"),
        tabPanel("Plot Output",
         uiOutput("PlotGeo"),
            value = "panel1"),
        tabPanel("Analysis",
                 tags$br(),
                 #busyIndicator(text = "Calculation in progress..",wait = 500),
                 strong(textOutput("txt1")),
                 uiOutput("anovaTable"),
                 tags$br(),
                 strong(textOutput("txt2")),
                 fluidRow(
                   column(6, plotOutput("histogramaRes")),
                   column(6, plotOutput("qqplotNorm"))
                 ),
                 fluidRow(plotOutput("scatter")),
                 uiOutput("tests"),
                 uiOutput("mct"),
                 #busyIndicator(text = "Calculation in progress..",wait = 500),
                 uiOutput("testeComp"),
                 uiOutput("dld"),
                 value = "panel2"), id = "inTabset"
      )

    )
  )
))


