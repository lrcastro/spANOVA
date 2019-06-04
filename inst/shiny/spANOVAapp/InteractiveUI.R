#Output que retorna os controles da parte de modelagem
#geoestatistica. Obs: falta por os graficos
output$geoUI <- renderUI({
  if(input$spem == "geoest"){
  return(
    tagList(
      #radioButtons("pe","Parameter Estimation",c("Ordinary Least Squares","EyeFit"),
      #            selected = "Ordinary Least Squares", inline = TRUE),
  selectInput("Thmodel","Correlation Function",
              choices = c("Spherical" = 'spherical',"Exponential" = 'exponential',
                          "Gaussian" = 'gaussian',"Wave" = 'wave',
                          "Circular" = 'circular',"Matern" = 'matern',"Cubic" = 'cubic',
                          "Cauchy" = 'cauchy',"Pure.nugget" = 'pure.nugget'),
              selected = "spherical"),
    radioButtons("trend", "Spatial Trend", c("Cte"='cte',"1st"),
                 selected = "cte", inline = T),
    sliderInput("dm", "Cutoff: % of Maximum Distance", 0.3, 1, step = 0.01,
                value = 0.50),
    #strong("Initial Parameters Values"),
    #checkboxInput("defaultSearch", "Default Search", TRUE),
    #uiOutput("iniParams"),
    strong("Parameter Estimates"),
    busyIndicator(text = "Calculation in progress..",wait = 500),
    textOutput("psillEst"),
    textOutput("RangeEst"),
    textOutput("NuggetEst"),
    textOutput("SillEst"),
    tags$hr(),
    div(style="display: inline-block;vertical-align:top; width: 230px;",bsButton("runS", "Compute Semivariogram",icon=icon("caret-right"), style="default")),
    div(style="display: inline-block; vertical-align:top; width: 10px;",bsButton("runA", "Run Analysis",icon=icon("ban"), style="danger",disabled = T))
    #div(style="display: inline-block;vertical-align:top; width: 150px;",downloadButton("report", "Download report"))
    #radioButtons('format', 'Save as', c('PDF', 'DOCX'),
    #           inline = TRUE,selected = 'PDF')
  )
  )}else{
    return(
      tagList(
        #strong("Sequence Radius"),
        #checkboxInput("sequence", label = "Default Sequence", value = TRUE),
        strong("Parameter Estimate"),
        textOutput("radius"),
        textOutput("rho"),
        tags$br(),
        strong("AIC"),
        textOutput("aic"),
        tags$hr(),
        div(style="display: inline-block;vertical-align:top; width: 150px;",bsButton("runA", "Run analysis",icon=icon("caret-right"), style="primary"))
        #div(style="display: inline-block;vertical-align:top; width: 150px;",downloadButton("report1", "Generate report")),
        # radioButtons('format1', 'Save as', c('PDF', 'DOCX'),
        #              inline = TRUE,selected = 'PDF')
      )
    )
  }
  })


#Output que retorna o tipo de delineamento a ser escolhido,
#para GEO as opcoes sao DIC e DBC, para AR as opcoes sao DIC, DBC

output$expd<-renderUI({
    selectInput("expd",label = "Choose an Experimental Design",
                choices = c("Completely Randomized Design" = "crd",
                    "Randomized Block Design" = "rbd"), selected = "crd")
})


#Output que retorna a selecao das variaveis do modelo
#condicionado ao tipo de delineamento escolhido e
#ao tipo de analise
output$variables <- renderUI({
    req(input$expd)
    if(input$expd == "crd"){
      return(
        tagList(
          selectInput("crdR","Response Variable",""),
          selectInput("crdF","Factor",""),
          strong("Match the columns corresponding to coordinates in your dataset"),
          div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput(inputId="xcoordCrd", label="X coord", choices = "")),
          div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput(inputId="ycoordCrd", label="Y coord", choices = "")),
          checkboxInput(inputId="scaleCoordCRD", label="Scale coordinates", value = FALSE),
          div(style="display: inline-block;vertical-align:top; width: 150px;",bsButton("buttonNext", "Next step",icon=icon("angle-double-right"),
                                                                                       style="default",type = "action"))#,
         # bsAlert("alert1") #alerta done
        )
      )
    }
    #variable selection rbd
    if(input$expd == "rbd"){
      return(
        tagList(
          selectInput("rbdR","Response Variable",""),
          selectInput("rbdF","Factor",""),
          selectInput("rbdB","Block",""),
          strong("Match the columns corresponding to coordinates in your dataset"),
          div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput(inputId="xcoordRbd", label="X coord", choices = "")),
          div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput(inputId="ycoordRbd", label="Y coord", choices = "")),
          checkboxInput(inputId="scaleCoordRCBD", label="Scale coordinates", value = FALSE),
          div(style="display: inline-block;vertical-align:top; width: 150px;",bsButton("buttonNext", "Next step",icon=icon("angle-double-right"),
                                                                                       style="default",type = "action"))
        )
      )
    }

})

#To show the table
output$dataTBS<-renderUI({
  if(input$dataType=='ifd'){
    return(
      tagList(
      helpText("Only '.csv' and '.txt' files are supported"),
      fileInput('file', 'Choose a file',
                accept = c('text/csv',
                           'text/comma-separated-values',
                           'text/tab-separated-values',
                           'text/plain',
                           '.csv',
                           '.tsv')),
     # helpText("Note: Your data file must contain X and Y coordinates in the 1st and 2nd column
      #         (in this order)"),
      checkboxInput('header', 'Header', TRUE),
      helpText("Check this box if your file has a header in the first line"),
      radioButtons('sep', 'Field separator character',
                   c("semicolon"=';',
                     "comma"=',',
                     "tab"='\t',
                     "white space"= ""),
                   ';'),
      radioButtons('dec', 'Character used for decimal points',
                   c("dot"='.',
                     "comma"=','),'.')
      ))} else {
      return(tagList(
        selectInput("loadedData", choices = c("crd_simulated", "candeia"),
                    selected = "crd_simulated", label = "Choose a dataset:"),
        div(style = "display: inline-block;vertical-align:top; width: 230px;",
            bsButton("runD", "Ok",icon = icon("caret-right"),
                     style="default"))
      ))
    }
})


#Conditions to show the plot menu
output$PlotGeo<-renderUI({
  #Geostatistics Ordinary Least Squares
  if(input$spem=="geoest"){
    return(
      tagList(
        fluidRow(
          column(6,
                 busyIndicator(text = "Calculation in progress..",wait = 500),
                 plotOutput("semivariog")),
          column(6,
                 busyIndicator(text = "Calculation in progress..",wait = 500),
                 plotOutput("chartQt")
          )
        ),
        fluidRow(
          column(6,
                 busyIndicator(text = "Calculation in progress..",wait = 500),
                 plotOutput("chart2")),
          column(6,
                 busyIndicator(text = "Calculation in progress..",wait = 500),
                 plotOutput("chart3")
          )
        )
      )
    )} else {
      return(
        tagList(
          fluidRow(
            column(12,
                   busyIndicator(text = "Calculation in progress..",wait = 500),
                   plotOutput("chartQt")
            )
          ),
          fluidRow(
            column(6,
                   busyIndicator(text = "Calculation in progress..",wait = 500),
                   plotOutput("chart2")),
            column(6,
                   busyIndicator(text = "Calculation in progress..",wait = 500),
                   plotOutput("chart3")
            )
          )
        )
      )
    }
})

