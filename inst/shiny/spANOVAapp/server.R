server <- function(input, output, session) {
  # sourcing the interactive part of ui
  source("InteractiveUI.R",local = T)

  # load the data set
  myData <- reactive({
    if(input$dataType=='ifd'){
      req(input$file)
      read.table(input$file$datapath,header = input$header,
                 sep = input$sep, dec = input$dec)
    }
  })

  # Pre-loaded dataset
  myData2 <- eventReactive(input$runD, {
    if(input$dataType!='ifd'){
      get(input$loadedData)
    }
  })


  # Output the data table
  output$table <-renderDataTable({
    if(input$dataType=='ifd'){
      datatable(myData(),filter = 'none')
    } else {
      datatable(myData2(),filter = 'none')
    }
  })


  # Response variable choice
  # Feed blank fields after chosing CRD
  observe({
    if(input$dataType=='ifd'){
      updateSelectInput(session, "crdR",choices = names(myData()))
    } else {
      updateSelectInput(session, "crdR",choices = names(myData2()))
    }
  })

  # Feed blank fields after chosing RCBD
  observe({
    if(input$dataType=='ifd'){
      updateSelectInput(session, "rbdR",choices = names(myData()))
    } else {
      updateSelectInput(session, "rbdR",choices = names(myData2()))
    }
  })

  # Explanatory variable choice
  # When the choice is made the name of the variable is saved, and thus
  # it will not appear in the next field


  # Factor levels to choose in crd
  observe({
    if(input$dataType=='ifd'){
      updateSelectInput(session, "crdF",
                        choices = names(myData())[!(names(myData()) %in% input$crdR)])
    } else {
      updateSelectInput(session, "crdF",
                        choices = names(myData2())[!(names(myData2()) %in% input$crdR)])
    }
  })

  # Factor levels to choose in RCBD
  observe({
    if(input$dataType=='ifd'){
    updateSelectInput(session, "rbdF",
                      choices = names(myData())[!(names(myData()) %in% input$rbdR)])
    } else {
      updateSelectInput(session, "rbdF",
                        choices = names(myData2())[!(names(myData2()) %in% input$rbdR)])
    }
  })

  # Blocking factor choice in RCBD
  observe({
    if(input$dataType=='ifd'){
      updateSelectInput(session, "rbdB",
                        choices = names(myData())[!(names(myData()) %in% input$rbdR) & !(names(myData()) %in% input$rbdF)])
    } else {
      updateSelectInput(session, "rbdB",
                        choices = names(myData2())[!(names(myData2()) %in% input$rbdR) & !(names(myData2()) %in% input$rbdF)])
    }
  })

  # Escolha das colunas restantes para ser as coordenadas em coordX e coordY
  observe({
    if(input$dataType=='ifd'){
      updateSelectInput(session, "xcoordCrd",
                        choices = names(myData())[!(names(myData()) %in% input$crdF) & !(names(myData()) %in% input$crdR)])
    } else {
      updateSelectInput(session, "xcoordCrd",
                        choices = names(myData2())[!(names(myData2()) %in% input$crdF) & !(names(myData2()) %in% input$crdR)])
    }
  })

  observe({
    if(input$dataType=='ifd'){
      updateSelectInput(session, "ycoordCrd",choices
                        = names(myData())[!(names(myData()) %in% input$crdF) & !(names(myData()) %in% input$crdR) & !(names(myData()) %in% input$xcoordCrd)])
    } else {
      updateSelectInput(session, "ycoordCrd",choices
                        = names(myData2())[!(names(myData2()) %in% input$crdF) & !(names(myData2()) %in% input$crdR) & !(names(myData2()) %in% input$xcoordCrd)])
    }
  })

  # Escolha das colunas restantes para ser as coordenadas em coordX e coordY no RBD
  observe({
    if(input$dataType=='ifd'){
      updateSelectInput(session, "xcoordRbd",choices
                        = names(myData())[!(names(myData()) %in% input$rbdF) &
                                            !(names(myData()) %in% input$rbdR) &  !(names(myData()) %in% input$rbdB)])
    } else {
      updateSelectInput(session, "xcoordRbd",choices
                        = names(myData2())[!(names(myData2()) %in% input$rbdF) &
                                             !(names(myData2()) %in% input$rbdR) &  !(names(myData2()) %in% input$rbdB)])
    }
  })

  observe({
    if(input$dataType=='ifd'){
      updateSelectInput(session, "ycoordRbd",choices
                        = names(myData())[!(names(myData()) %in% input$rbdF) & !(names(myData()) %in% input$rbdR) &
                                            !(names(myData()) %in% input$xcoordRbd) & !(names(myData()) %in% input$rbdB)])
    } else {
      updateSelectInput(session, "ycoordRbd",choices
                        = names(myData2())[!(names(myData2()) %in% input$rbdF) & !(names(myData2()) %in% input$rbdR) &
                                             !(names(myData2()) %in% input$xcoordRbd) & !(names(myData2()) %in% input$rbdB)])
    }
  })


  # Quando clicar em next mudar para o controles de modelagem
  observeEvent(input$buttonNext,({
    updateTabsetPanel(session, "sidepanel",
                      selected = "modelCont")
    updateTabsetPanel(session, "inTabset",
                      selected = "panel1")
  }))


  # Criacao do objeto geodata
  geodados <- reactive({
    if(input$spem == "geoest"){
      if(input$crdR != ""){
        if(input$dataType=='ifd'){
          x <- as.geodata(myData(), coords.col = c(which(names(myData()) == input$xcoordCrd),
                                                   which(names(myData()) == input$ycoordCrd)),
                          data.col = which(names(myData()) == input$crdR),
                          covar.col = which(names(myData()) == input$crdF))
        } else {
          x <- as.geodata(myData2(), coords.col = c(which(names(myData2()) == input$xcoordCrd),
                                                    which(names(myData2()) == input$ycoordCrd)),
                          data.col = which(names(myData2()) == input$crdR),
                          covar.col = which(names(myData2()) == input$crdF))
        }
      } else {
        if(input$dataType=='ifd'){
          x <- as.geodata(myData(), coords.col = c(which(names(myData()) == input$xcoordRbd),
                                                   which(names(myData()) == input$ycoordRbd)),
                          data.col = which(names(myData()) == input$rbdR),
                          covar.col = c(which(names(myData()) == input$rbdF),
                                        which(names(myData()) == input$rbdB)))
        } else {
          x <- as.geodata(myData2(), coords.col = c(which(names(myData2()) == input$xcoordRbd),
                                                    which(names(myData2()) == input$ycoordRbd)),
                          data.col = which(names(myData2()) == input$rbdR),
                          covar.col = c(which(names(myData2()) == input$rbdF),
                                        which(names(myData2()) == input$rbdB)))
        }
      }
      return(x)
    }
  })

  # Determinacao da distancia maxima
  dist.max <- reactive({
    if(input$spem == "geoest"){
      h_max <- summary(geodados())[[3]][[2]]
      dist_mx <- input$dm*h_max
      return(dist_mx)
    }
  })


  # Esse trecho cria o objeto spvariog
  semivar <- reactive({
    if(input$spem == "geoest"){
      validate(
        need(input$dm != "", "Please configure modelling controls")
      )
      if(input$crdR != ""){
        b <- spVariog(geodata = geodados(), trend = input$trend, max.dist = dist.max(),
                      design = "crd", scale = input$scaleCoordCRD)
      } else {
        b <- spVariog(geodata = geodados(), trend = input$trend, max.dist = dist.max(),
                      design = "rcbd", scale = input$scaleCoordRCBD)
      }
      return(b)
    }
  })


  # Esse trecho faz o ajuste do modelo via OLS
  fit <- reactive({
    if(input$spem == "geoest"){
      if(input$iniPar != "default"){
        req(input$Nug)
        spVariofit(semivar(), max.dist = dist.max(), cov.model = input$Thmodel,
                   weights = input$estMethod, nugget = as.numeric(input$Nug),
                   ini.cov.pars = c(as.numeric(input$Sil), as.numeric(input$Ran)))
      } else {
        spVariofit(semivar(), max.dist = dist.max(), cov.model = input$Thmodel,
                   weights = input$estMethod)
      }
    }
  })

  # Grafico semivariograma experimental e ajuste
  output$semivariog <- renderPlot({
    if(input$spem == "geoest"){
      plot(semivar(),main = "Initial Semivariogram")
      lines(fit(),col=1)
    }
  })


  source("grafico1.R",local = T)

  # Segundo grafico (apos o semivariog)
  output$chartQt<-renderPlot({
    if(input$crdR != ""){
      if(input$dataType=='ifd'){
        ch1 <- chart1(myData(),resp = input$crdR,
                      trat = input$crdF, coordX = input$xcoordCrd, coordY = input$ycoordCrd)
      } else {
        ch1 <- chart1(myData2(),resp = input$crdR,
                      trat = input$crdF, coordX = input$xcoordCrd, coordY = input$ycoordCrd)
      }
    } else {
      if(input$dataType=='ifd'){
        ch1 <- chart1(myData(),resp = input$rbdR,
                      trat = input$rbdF, coordX = input$xcoordRbd, coordY = input$ycoordRbd)
      } else {
        ch1 <- chart1(myData2(),resp = input$rbdR,
                      trat = input$rbdF, coordX = input$xcoordRbd, coordY = input$ycoordRbd)
      }
    }
    return(ch1)
  })

  # Terceiro grafico
  output$chart2 <- renderPlot({
    if(input$crdR != ""){
      if(input$dataType=='ifd'){
        ch2 <- chart2(myData(),resp = input$crdR,
                      trat = input$crdF, coordX = input$xcoordCrd, coordY = input$ycoordCrd)
      } else {
        ch2 <- chart2(myData2(),resp = input$crdR,
                      trat = input$crdF, coordX = input$xcoordCrd, coordY = input$ycoordCrd)
      }
    } else {
      if(input$dataType=='ifd'){
        ch2 <- chart2(myData(),resp = input$rbdR,
                      trat = input$rbdF, coordX = input$xcoordRbd, coordY = input$ycoordRbd)
      } else {
        ch2 <- chart2(myData2(),resp = input$rbdR,
                      trat = input$rbdF, coordX = input$xcoordRbd, coordY = input$ycoordRbd)
      }
    }
    return(ch2)
  })

  # Quarto grafico
  output$chart3 <- renderPlot({
    if(input$crdR != ""){
      if(input$dataType=='ifd'){
        ch3 <- chart3(myData(), resp = input$crdR,
                      trat = input$crdF, coordX = input$xcoordCrd, coordY = input$ycoordCrd)
      } else {
        ch3 <- chart3(myData2(), resp = input$crdR,
                      trat = input$crdF, coordX = input$xcoordCrd, coordY = input$ycoordCrd)
      }
    }else{
      if(input$dataType=='ifd'){
        ch3 <- chart3(myData(), resp = input$rbdR,
                      trat = input$rbdF, coordX = input$xcoordRbd, coordY = input$ycoordRbd)
      } else {
        ch3 <- chart3(myData2(), resp = input$rbdR,
                      trat = input$rbdF, coordX = input$xcoordRbd, coordY = input$ycoordRbd)
      }
    }
    return(ch3)
  })

  # Text with the parameter estimates
  # psill
  output$psillEst<-renderText({
    req(input$runA)
    if(input$spem == "geoest"){
      return(paste("Partial Sill:", round(fit()$mod$cov.pars[1],3)))
    }
  })

  # sill
  output$SillEst <- renderText({
    req(input$runA)
    if(input$spem == "geoest"){
      return(paste("Sill:",round(fit()$mod$cov.pars[1]+fit()$mod$nugget,3)))
    }
  })

  # Range
  output$RangeEst <- renderText({
    req(input$runA)
    if(input$spem == "geoest"){
      return(paste("Range:",round(fit()$mod$cov.pars[2],3)))
    }
  })

  # nugget
  output$NuggetEst <- renderText({
    req(input$runA)
    if(input$spem == "geoest"){
      return(paste("Nugget:",round(fit()$mod$nugget,3)))
    }
  })


  # Observe event para o botao run analysis
  # observeEvent(input$runS,({
  #   updateButton(session, "runA", style = "success",disabled = FALSE ,icon = icon("check"))
  # }))

  # # Quando clicar em compute semivariogram mudar para os graficos
  # observeEvent(input$runA,({
  #   updateTabsetPanel(session, "inTabset",
  #                     selected = "panel1")
  #   # write.table(input$runA,"run.txt")
  # }))

  # Quando clicar em run analysis mudar para o controles de modelagem
  observeEvent(input$runA,({
    updateTabsetPanel(session, "inTabset",
                      selected = "panel2")
    # write.table(input$runA,"run.txt")
  }))

  # Textos
  output$txt1 <- renderText({
    req(input$runA)
    if(input$runA > 0){
      return("Analysis of Variance Table")
    }
  })

  output$txt2 <- renderText({
    req(input$runA)
    if(input$runA > 0){
      return("Checking the Residuals")
    }
  })

  # Modelo Geo e AR
  an <- eventReactive(input$runA,{
    # Modelo geoestatistico crd ou rcbd
    if(input$spem == "geoest"){
      ant <- aovGeo(fit(), cutoff = input$dm)
    } else {
      if(input$dataType=='ifd'){
        coordenadasCRD <- cbind(myData()[,which(names(myData()) == input$xcoordCrd)],
                                myData()[,which(names(myData()) == input$ycoordCrd)])
        coordenadsRCBD <- cbind(myData()[,which(names(myData()) == input$xcoordRbd)],
                                myData()[,which(names(myData()) == input$ycoordRbd)])
        responseCRD <- as.numeric(myData()[,which(names(myData()) == input$crdR)])
        responseRCBD <- as.numeric(myData()[,which(names(myData()) == input$rbdR)])
        factorCRD <- as.numeric(myData()[,which(names(myData()) == input$crdF)])
        factorRCBD <- as.numeric(myData()[,which(names(myData()) == input$rbdF)])
        blocoRCBD <- as.numeric(myData()[,which(names(myData()) == input$rbdB)])
      } else {
        coordenadasCRD <- cbind(myData2()[,which(names(myData2()) == input$xcoordCrd)],
                                myData2()[,which(names(myData2()) == input$ycoordCrd)])
        coordenadsRCBD <- cbind(myData2()[,which(names(myData2()) == input$xcoordRbd)],
                                myData2()[,which(names(myData2()) == input$ycoordRbd)])
        responseCRD <- as.numeric(myData2()[,which(names(myData2()) == input$crdR)])
        responseRCBD <- as.numeric(myData2()[,which(names(myData2()) == input$rbdR)])
        factorCRD <- as.numeric(myData2()[,which(names(myData2()) == input$crdF)])
        factorRCBD <- as.numeric(myData2()[,which(names(myData2()) == input$rbdF)])
        blocoRCBD <- as.numeric(myData2()[,which(names(myData2()) == input$rbdB)])
      }
      if(input$expd == "crd"){
        ant <- aovSar.crd(responseCRD,
                          factorCRD,
                          coordenadasCRD)
      } else {
        ant <- aovSar.rcbd(responseRCBD,
                           factorRCBD,
                           blocoRCBD,
                           coordenadsRCBD)
      }

      return(ant)
    }
  })

  # Summary para o modelo AR
  smry <- reactive({
    if(input$spem == "arm"){
      z <- summary(an())
    }
    return(z)
  })

  # Parameters for AR model: radius
  output$radius <- renderText({
    if(input$spem == "arm"){
      return(paste("Radius:", round(smry()[1],3)))
    }
  })

  # Parameters for AR model: rho
  output$rho <- renderText({
    if(input$spem == "arm"){
      return(paste("Rho:", round(smry()[2],3)))
    }
  })

  # AIC
  output$aic <- renderText({
    if(input$spem == "arm"){
      return(as.character(round(smry()[3],3)))
    }
  })

  antable <- eventReactive(input$runA,{
    # Tabela de analise de variancia
    tabl <- anova(an())
    return(tabl)
  })

  output$anovaTable<-renderUI({
    M <- print(xtable(antable(), caption = "Analysis of Variance"),
               floating = FALSE, tabular.environment = "array", comment = FALSE,
               print.results = FALSE, include.rownames = TRUE,
               sanitize.text.function = identity)

    # Colocar o $ no inicio e no fim do codigo (modo matematico)
    html <- paste0("$$", M, "$$")

    # Exportar a tabela para MathJax em html
    list(withMathJax(HTML(html)))

  })


  # Teste de Shapiro
  Tabela <- eventReactive(input$runA, {
    # Testes
    sht <- shapiro.test(an()$residuals)$p.value

    # Teste Moran.I
    # obtendo as coordenadas
    if(input$dataType=='ifd'){
      coordenadasCRD <- cbind(myData()[,which(names(myData()) == input$xcoordCrd)],
                              myData()[,which(names(myData()) == input$ycoordCrd)])
      coordenadsRCBD <- cbind(myData()[,which(names(myData()) == input$xcoordRbd)],
                              myData()[,which(names(myData()) == input$ycoordRbd)])
      if(input$expd == "crd"){
        Xcoordinates <-(coordenadasCRD[,1])-min(coordenadasCRD[,1])
        Ycoordinates <-(coordenadasCRD[,2])-min(coordenadasCRD[,2])
      } else {
        Xcoordinates <-(coordenadsRCBD[,1])-min(coordenadsRCBD[,1])
        Ycoordinates <-(coordenadsRCBD[,2])-min(coordenadsRCBD[,2])
      }
    } else {
      coordenadasCRD <- cbind(myData2()[,which(names(myData2()) == input$xcoordCrd)],
                              myData2()[,which(names(myData2()) == input$ycoordCrd)])
      coordenadsRCBD <- cbind(myData2()[,which(names(myData2()) == input$xcoordRbd)],
                              myData2()[,which(names(myData2()) == input$ycoordRbd)])
      if(input$expd == "crd"){
        Xcoordinates <-(coordenadasCRD[,1])-min(coordenadasCRD[,1])
        Ycoordinates <-(coordenadasCRD[,2])-min(coordenadasCRD[,2])
      } else {
        Xcoordinates <-(coordenadsRCBD[,1])-min(coordenadsRCBD[,1])
        Ycoordinates <-(coordenadsRCBD[,2])-min(coordenadsRCBD[,2])
      }
    }

    # Calculando o teste
    coordinates <-cbind(Xcoordinates,Ycoordinates)
    distance <- as.matrix(dist(coordinates))
    distance.inv <- distance^(-1)
    diag(distance.inv) <- 0
    moran.pvalue <- Moran.I(an()$residuals, distance.inv, alt = "g")$p.value


    # Tabela
    tbla <- data.frame(Assumption = c("Normality", "Independence"),
                       Test = c("Shapiro-Wilk", "Moran-I"),
                       P.value=c(sht, moran.pvalue))
    return(tbla)
  })

  # Histograma com densidade dos residuos
  output$histogramaRes <- renderPlot({
    hist(an()$residuals, xlab = "Residuals", main = "", freq = FALSE)
    lines(density(an()$residuals))
  })

  # QQnorm
  output$qqplotNorm <- renderPlot({
    qqnorm(an()$residuals)
    qqline(an()$residuals)
  })

  # Scatter-plot residuals
  output$scatter <- renderPlot({
    plot(an()$residuals, ylab = "Residuals")
  })


  # Exportando em HTML
  output$tests <- renderUI({

    M <- print(xtable(Tabela(), caption = "Checking Model Assumptions"),
               floating = FALSE, tabular.environment = "array", comment=FALSE,
               print.results = FALSE, include.rownames = FALSE,
               sanitize.text.function = identity)

    # Colocar o $$ no inicio e no fim do codigo (modo matematico)
    html <- paste0("$$", M, "$$")

    # Exportar a tabela para MathJax em html
    list(withMathJax(HTML(html)))

  })

  # Escolha do teste de comparacao multipla
  output$mct<-renderUI({
      tagList(
        div(style="display: inline-block;vertical-align:center; width: 230px",
            selectInput("testeCM","Multiple-comparison procedure:",
                        choices = c("Tukey","Mutivariate T", "Scott-Knott"),
                        selected = "Scott-Knott")),

        div(style="display: inline-block;vertical-align:14px;width:230px",
            bsButton("runTes", "OK",icon=icon("caret-right"),
                     style="primary",size = "small"))
      )
  })


  # Multicomparisons tests

  mcompTests <- eventReactive(input$runTes,{
    #write.table(trd,"trend.txt")
    if(input$testeCM=="Scott-Knott"){
      w <- spScottKnott(an())
    }
    if(input$testeCM=="Mutivariate T"){
      w <- spMVT(an())
    }
    if(input$testeCM=="Tukey"){
      w <- spTukey(an())
    }
    return(w)
  })

  output$testeComp<-renderUI({
    M <- print(xtable(mcompTests(), caption = "Multicomparison Test"),
               floating = FALSE, tabular.environment = "array", comment = FALSE,
               print.results = FALSE, include.rownames = TRUE,
               sanitize.text.function = identity)

    # Colocar o $$ no inicio e no fim do codigo (modo matematico)
    html <- paste0("$$", M, "$$")

    # Exportar a tabela para MathJax em html
    list(withMathJax(HTML(html)))
  })


  # To pass as parameter to rmrkdown
  semivPar <- reactive({
    if(input$spem == "geoest"){
      d <- data.frame(Partial_Sill = round(fit()$mod$cov.pars[1],3),
                      Range = round(fit()$mod$cov.pars[2],3),
                      Nugget = round(fit()$mod$nugget,3),
                      Sill = round(fit()$mod$cov.pars[1]+fit()$mod$nugget,3))
      return(d)
    }
  })


  # Report
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = function() {
      paste('SpAnova-report', sep = '.', 'docx')
    },
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      if(input$spem == "geoest"){
        src <- normalizePath('reportGeo.Rmd')

        # temporarily switch to the temp dir, in case you do not have write
        # permission to the current working directory

        tempReport <- file.path(tempdir(), "reportGeo.Rmd")
        file.copy("reportGeo.Rmd", tempReport, overwrite = TRUE)

        #owd <- setwd(tempdir())
        #owd<-tempdir()
        #on.exit(setwd(owd))
        #arq<-'reportGeo.Rmd'
        #fl<-paste(owd,arq,sep = "/")
        #file.copy(src, fl, overwrite = TRUE)
        # Set up parameters to pass to Rmd document
        params <- list(tab = antable(),
                       modelVariog = semivar(), table = Tabela(),
                       semiPar = semivPar(), CorF = input$Thmodel,
                       McomP = mcompTests(), McompName = input$testeCM,
                       modelVariofit = fit(), modelGeo = an())
      } else {
        src <- normalizePath('reportAr.Rmd')

        # temporarily switch to the temp dir, in case you do not have write
        # permission to the current working directory

        tempReport <- file.path(tempdir(), "reportAr.Rmd")
        file.copy("reportAr.Rmd", tempReport, overwrite = TRUE)

        #owd <- setwd(tempdir())
        #owd<-tempdir()
        #on.exit(setwd(owd))
        #arq<-'reportGeo.Rmd'
        #fl<-paste(owd,arq,sep = "/")
        #file.copy(src, fl, overwrite = TRUE)
        # Set up parameters to pass to Rmd document
        params <- list(tab = antable(),
                       paramEst = smry(),
                       table = Tabela(),
                       modelAR = an(),
                       McomP = mcompTests(), McompName = input$testeCM)
      }

      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      render(tempReport,output_format = word_document(),
             output_file = file, params = params,
      envir = new.env(parent = globalenv())
      )
    }
  )

  # Botao de Download do report geo

  output$dld <- renderUI({
    req(input$runTes)
    tagList(
      div(style="display: inline-block;vertical-align:top; width: 150px;",
          downloadButton("report", "Download report"))
    )
  })

  # Selecao do metodo de estimacao
  output$iniParams <- renderUI({
    if(input$iniPar == "values"){
      tagList(
        sliderInput("Nug", "Nugget", min = 0, max = round(2*max(semivar()$vario.res$v),1),
                    0.1*max(semivar()$vario.res$v),0.1),
        sliderInput("Sil", "Sill",  min = 0, max = round(2*max(semivar()$vario.res$v),1),
                    0.8*max(semivar()$vario.res$v),0.1),
        sliderInput("Ran", "Range", min = 0, max = round(2*max(semivar()$vario.res$u),1),
                    max(semivar()$vario.res$u)/3,0.1)
      )
    }
  })

}
