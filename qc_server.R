library(tidyverse)
library(shiny)

qc_server <- function(input, output, session) { 
  mylist <- reactiveVal()
  observeEvent(
    c(
      input$includeBaltimore,
      input$includeBrazil,
      input$includeChicago,
      input$includeDenver,
      input$includeNigeria,
      input$includeWashington_DC
    ),
    {
      if (input$includeDenver != TRUE &
          input$includeChicago != TRUE &
          input$includeNigeria != TRUE) {
        updateCheckboxInput(session = session,
                            "includeChild", "Children", FALSE)
      }
      if (input$includeDenver == TRUE |
          input$includeChicago == TRUE |
          input$includeNigeria == TRUE) {
        updateCheckboxInput(session = session,
                            "includeChild", "Children", TRUE)
      }
    }
  )
  observeEvent(c(input$includeChild), {
    if (input$includeDenver != TRUE &
        input$includeChicago != TRUE &
        input$includeNigeria != TRUE) {
      updateCheckboxInput(session = session,
                          "includeChild", "Children", FALSE)
    }
  })
  
  # Force Adults Selection
  observeEvent(
    c(
      input$includeBaltimore,
      input$includeBrazil,
      input$includeChicago,
      input$includeDenver,
      input$includeNigeria,
      input$includeWashington_DC
    ),
    {
      if (input$includeBaltimore != TRUE &
          input$includeBrazil != TRUE &
          input$includeNigeria != TRUE &
          input$includeWashington_DC != TRUE &
          input$includeChicago != TRUE) {
        updateCheckboxInput(session = session,
                            "includeAdult", "Adults", FALSE)
      }
      if (input$includeBaltimore == TRUE |
          input$includeBrazil == TRUE |
          input$includeNigeria == TRUE |
          input$includeWashington_DC == TRUE |
          input$includeChicago == TRUE) {
        updateCheckboxInput(session = session,
                            "includeAdult", "Adults", TRUE)
      }
      if ((as.character(input$CellType) != "Nasal")) {
        updateCheckboxInput(session = session,
                            "includeBaltimore", "Baltimore", FALSE)
      }
    }
  )
  observeEvent(c(input$includeAdult), {
    if (input$includeBaltimore != TRUE &
        input$includeBrazil != TRUE &
        input$includeNigeria != TRUE &
        input$includeWashington_DC != TRUE &
        input$includeChicago != TRUE) {
      updateCheckboxInput(session = session,
                          "includeAdult", "Adults", FALSE)
    }
  })
  
  # Update inputs for cell type selection
  observeEvent(input$CellType, {
    if ((as.character(input$CellType) != "Nasal")) {
      updateSelectInput(session = session,
                        inputId = 'PASS_ONLY',
                        label = "Remove samples that failed PBMC_cellcounts?")
      updateSelectInput(session = session,
                        inputId = 'Include_Boxplots',
                        label = "Should the correlation plots feture an extra row & column illustrating PASS/FAIL for PBMC_cellcounts?")
      updateCheckboxInput(session = session,
                          "includeBaltimore", "Baltimore", FALSE)
    }
    if ((as.character(input$CellType) != "PBMC")) {
      updateSelectInput(session = session,
                        inputId = 'PASS_ONLY',
                        label = "Remove samples that failed Slide_status?")
      updateSelectInput(session = session,
                        inputId = 'Include_Boxplots',
                        label = "Should the correlation plots feture an extra row & column illustrating PASS/FAIL for Slide_status?")
      updateCheckboxInput(session = session,
                          "includeBaltimore", "Baltimore", TRUE)
    }
  })
  
  # Update inputs for DNA or RNA
  observeEvent(input$DNAorRNA, {
    if ((as.character(input$DNAorRNA) != "DNA")) {
      updateSliderInput(
        session = session,
        inputId = 'LUTHRESH_Nanodrop_260_280',
        value = c(1.70, 2.20)
      )
      updateSliderInput(session = session,
                        inputId = 'LTHRESH_TOTAL_Nanodrop_ug',
                        value = c(0.60))
      updateSliderInput(session = session,
                        inputId = 'LTHRESH_TOTAL_Qubit_ug',
                        value = c(0.60))
    }
    if ((as.character(input$DNAorRNA) != "RNA")) {
      updateSliderInput(
        session = session,
        inputId = 'LUTHRESH_Nanodrop_260_280',
        value = c(1.40, 2.15)
      )
      updateSliderInput(session = session,
                        inputId = 'LTHRESH_TOTAL_Nanodrop_ug',
                        value = c(0.75))
      updateSliderInput(session = session,
                        inputId = 'LTHRESH_TOTAL_Qubit_ug',
                        value = c(0.75))
    }
  })
  
  observe({
    mylist(
      list(
        # Fname = as.character(input$Fname),
        filterMissing = input$filterMissing,
        includeBaltimore = input$includeBaltimore,
        includeBrazil = input$includeBrazil,
        includeChicago = input$includeChicago,
        includeDenver = input$includeDenver,
        includeNigeria = input$includeNigeria,
        includeAdult = input$includeAdult,
        includeChild = input$includeChild,
        includeAsthma = input$includeAsthma,
        includeNonAsthma = input$includeNonAsthma,
        includeMales = input$includeMales,
        includeFemales = input$includeFemales,
        includeWashington_DC = input$includeWashington_DC,
        PBMC_PASS_ONLY = as.character(input$PASS_ONLY),
        Slide_status_PASS_ONLY = as.character(input$PASS_ONLY),
        Include_Boxplots = as.character(input$Include_Boxplots),
        CellType = as.character(input$CellType),
        DNAorRNA = as.character(input$DNAorRNA),
        # LTHRESH_Nanodrop_260_280 = as.numeric(input$LTHRESH_Nanodrop_260_280),
        # UTHRESH_Nanodrop_260_280 = as.numeric(input$UTHRESH_Nanodrop_260_280),
        LTHRESH_Nanodrop_260_280 = as.numeric(input$LUTHRESH_Nanodrop_260_280[1]),
        UTHRESH_Nanodrop_260_280 = as.numeric(input$LUTHRESH_Nanodrop_260_280[2]),
        LTHRESH_Agilent_DIN = as.numeric(input$LTHRESH_Agilent_DIN),
        # UTHRESH_Agilent_DIN = as.numeric(input$UTHRESH_Agilent_DIN),
        LTHRESH_Agilent_RINe = as.numeric(input$LTHRESH_Agilent_RINe),
        # UTHRESH_Agilent_RINe = as.numeric(input$UTHRESH_Agilent_RINe),
        LTHRESH_TOTAL_Nanodrop_ug = as.numeric(input$LTHRESH_TOTAL_Nanodrop_ug),
        # UTHRESH_TOTAL_Nanodrop_ug = as.numeric(input$UTHRESH_TOTAL_Nanodrop_ug),
        LTHRESH_TOTAL_Qubit_ug = as.numeric(input$LTHRESH_TOTAL_Qubit_ug)
        # UTHRESH_TOTAL_Qubit_ug = as.numeric(input$UTHRESH_TOTAL_Qubit_ug)
      )
    )
  })
  
  # data <- reactive({
  #   req(input$file1)
  #
  #   inFile <- input$file1
  #
  #
  #   df <-
  #     read.delim(
  #       inFile$datapath,
  #       header = input$header,
  #       sep = input$sep,
  #       quote = input$quote
  #     )
  #
  #
  #   return(df)
  # })
  #
  # output$contents <- renderTable({
  #   data()
  # })
  
  observeEvent(input$runScript, {
    env_list <<- mylist()
    save(env_list, file = "env_list.RData")
    source("./source_plots_shiny.R", local = list2env(mylist()))
    load("savePlots.RData")
    
    output$cGram <- renderPlot({
      cGram
    }, width = "auto", height = 520)
    
    output$mytable <- renderUI(print(summaryTable))
    
    output$myGridPlots1 <- renderPlotly({
      ggplotly(myGridPlots[[1]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    output$myGridPlots2 <- renderPlotly({
      ggplotly(myGridPlots[[2]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    output$myGridPlots3 <- renderPlotly({
      ggplotly(myGridPlots[[3]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    output$myGridPlots4 <- renderPlotly({
      ggplotly(myGridPlots[[4]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    # output$myGridPlots5 <- renderPlotly({
    #   ggplotly(myGridPlots[[5]])  %>% layout(legend = list(orientation = "V", x = 1, y =1))
    # })
    
    output$myGridPlots6 <- renderPlotly({
      ggplotly(myGridPlots[[6]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    output$myGridPlots7 <- renderPlotly({
      ggplotly(myGridPlots[[7]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    output$myGridPlots8 <- renderPlotly({
      ggplotly(myGridPlots[[8]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    # output$myGridPlots9 <- renderPlotly({
    #   ggplotly(myGridPlots[[9]])  %>% layout(legend = list(orientation = "V", x = 1, y =1))
    # })
    #
    # output$myGridPlots10 <- renderPlotly({
    #   ggplotly(myGridPlots[[10]])  %>% layout(legend = list(orientation = "V", x = 1, y =1))
    # })
    
    output$myGridPlots11 <- renderPlotly({
      ggplotly(myGridPlots[[11]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    output$myGridPlots12 <- renderPlotly({
      ggplotly(myGridPlots[[12]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
    
    # output$myGridPlots13 <- renderPlotly({
    #   ggplotly(myGridPlots[[13]])  %>% layout(legend = list(orientation = "V", x = 1, y =1))
    # })
    #
    # output$myGridPlots14 <- renderPlotly({
    #   ggplotly(myGridPlots[[14]])  %>% layout(legend = list(orientation = "V", x = 1, y =1))
    # })
    #
    # output$myGridPlots15 <- renderPlotly({
    #   ggplotly(myGridPlots[[15]])  %>% layout(legend = list(orientation = "V", x = 1, y =1))
    # })
    
    output$myGridPlots16 <- renderPlotly({
      ggplotly(myGridPlots[[16]])  %>% layout(legend = list(
        orientation = "V",
        x = 1,
        y = 1
      ))
    })
  })
  
  # output$MyPlot <- renderPlot({
  #
  #   x <- data()
  #   ggpairs(x, columns = c(6,12,18,19),lower = list(continuous = wrap("smooth",method="pearson",linetype="blank", alpha = 0.5, size=2)), ggplot2::aes(colour=Slide_status))
  #
  # })
  
  
  # observeEvent(input$plateSamples, {
  #   tableForCorey <- get_manifest()
  #   output$tableForCorey <- DT::renderDataTable({ tableForCorey })
  # })
}