library(shiny)
library(datasets)
# library(ggpairs)
require(GGally)
require(ggplot2)
library(tidyverse)
library(gridExtra)
library(plotly)
library(DT)
library(reshape2)
library(table1)
library(pals)
library(ggthemes)

# install.packages("shiny",
#                  "datasets",
#                  "GGally",
#                  "ggplot2",
#                  "tidyverse",
#                  "gridExtra",
#                  "plotly",
#                  "DT")

source("manifest_shiny.R")

ui <- shinyUI(fluidPage(
  titlePanel("CAAPA2 QC Metrics"),
  tabsetPanel(
    id = "tabs",
    tabPanel(
      titlePanel("Data Setup & QC"),
      sidebarPanel(
        tags$h4("Testing Sites"),
        checkboxInput("includeBaltimore", "Baltimore", TRUE),
        checkboxInput("includeBrazil", "Brazil", TRUE),
        checkboxInput("includeChicago", "Chicago", TRUE),
        checkboxInput("includeDenver", "Denver", TRUE),
        checkboxInput("includeNigeria", "Nigeria", TRUE),
        checkboxInput("includeWashington_DC", "Washington DC", TRUE),
        tags$h4("Age Groups"),
        checkboxInput("includeAdult", "Adults", TRUE),
        checkboxInput("includeChild", "Children", TRUE),
        tags$h4("Asthma Status"),
        checkboxInput("includeAsthma", "Asthmatics", TRUE),
        checkboxInput("includeNonAsthma", "Non-Asthmatics", TRUE),
        tags$h4("Gender"),
        checkboxInput("includeMales", "Males", TRUE),
        checkboxInput("includeFemales", "Females", TRUE),
        
        selectInput(
          "CellType",
          label = "Select cell type",
          c("PBMC",
            "Nasal"),
          selected = "Nasal"
        ),
        
        selectInput(
          "DNAorRNA",
          label = "Select DNA or RNA",
          c("DNA",
            "RNA"),
          selected = "DNA"
        ),
        
        selectInput(
          "PASS_ONLY",
          label = "Remove samples that failed PBMC_cellcounts?",
          c("Yes", "No"),
          selected = "Yes"
        ),
        
        selectInput(
          "Include_Boxplots",
          label = "Should the correlation plots feture an extra row & column illustrating PASS/FAIL for PBMC_cellcounts?",
          #or Slide_status (for Nasal)?
          c("Yes", "No"),
          selected = "No"
        ),
        
        tags$head(
          tags$style(
            type = "text/css",
            '
            .js-irs-1 .irs-line-mid{
            background: #428bca ;
            border: 1px solid #428bca ;
            }
            .js-irs-1 .irs-line-right{
            background: #428bca ;
            }
            .js-irs-1 .irs-bar {
            background: linear-gradient(to bottom, #DDD -50%, #FFF 150%);
            border-top: 1px solid #CCC ;
            border-bottom: 1px solid #CCC ;
            }
            .js-irs-1 .irs-bar-edge {
            background: inherit ;
            border: inherit ;
            }
            
            '
            )
            ),
        
        tags$head(
          tags$style(
            type = "text/css",
            '
            .js-irs-2 .irs-line-mid{
            background: #428bca ;
            border: 1px solid #428bca ;
            }
            .js-irs-2 .irs-line-right{
            background: #428bca ;
            }
            .js-irs-2 .irs-bar {
            background: linear-gradient(to bottom, #DDD -50%, #FFF 150%);
            border-top: 1px solid #CCC ;
            border-bottom: 1px solid #CCC ;
            }
            .js-irs-2 .irs-bar-edge {
            background: inherit ;
            border: inherit ;
            }
            
            '
            )
            ),
        
        tags$head(
          tags$style(
            type = "text/css",
            '
            .js-irs-3 .irs-line-mid{
            background: #428bca ;
            border: 1px solid #428bca ;
            }
            .js-irs-3 .irs-line-right{
            background: #428bca ;
            }
            .js-irs-3 .irs-bar {
            background: linear-gradient(to bottom, #DDD -50%, #FFF 150%);
            border-top: 1px solid #CCC ;
            border-bottom: 1px solid #CCC ;
            }
            .js-irs-3 .irs-bar-edge {
            background: inherit ;
            border: inherit ;
            }
            
            '
            )
            ),
        
        tags$head(
          tags$style(
            type = "text/css",
            '
            .js-irs-4 .irs-line-mid{
            background: #428bca ;
            border: 1px solid #428bca ;
            }
            .js-irs-4 .irs-line-right{
            background: #428bca ;
            }
            .js-irs-4 .irs-bar {
            background: linear-gradient(to bottom, #DDD -50%, #FFF 150%);
            border-top: 1px solid #CCC ;
            border-bottom: 1px solid #CCC ;
            }
            .js-irs-4 .irs-bar-edge {
            background: inherit ;
            border: inherit ;
            }
            
            '
            )
            ),
        
        sliderInput(
          "LUTHRESH_Nanodrop_260_280",
          label = HTML(
            "Lower & Upper Thresholds for Nanodrop_260_280,<br/>DNA ≈ 1.4 to 2.15<br/>RNA ≈ 1.7 to 2.2"
          ),
          min = 1,
          max = 3,
          value = c(1.40, 2.15),
          step = 0.01
        ),
        
        
        # textInput("LTHRESH_Agilent_DIN", label="Lower Threshold for Agilent_DIN", value = 6),
        sliderInput(
          "LTHRESH_Agilent_DIN",
          label = HTML("Lower Threshold for Agilent_DIN"),
          min = 0,
          max = 10,
          value = 6.00,
          step = 0.01
        ),
        
        # textInput("UTHRESH_Agilent_DIN", label="Upper Threshold for Agilent_DIN", value = Inf),
        # place holder for slider
        
        # textInput("LTHRESH_Agilent_RINe", label="Lower Threshold for Agilent_RINe", value = 6),
        sliderInput(
          "LTHRESH_Agilent_RINe",
          label = HTML("Lower Threshold for Agilent_RINe"),
          min = 0,
          max = 10,
          value = 6.00,
          step = 0.01
        ),
        
        # textInput("UTHRESH_Agilent_RINe", label="Upper Threshold for Agilent_RINe", value = Inf),
        # place holder for slider
        
        # textInput("LTHRESH_TOTAL_Nanodrop_ug", label="Lower Threshold for TOTAL_Nanodrop_ug", value = 0.75),
        sliderInput(
          "LTHRESH_TOTAL_Nanodrop_ug",
          label = HTML(
            "Lower Threshold for TOTAL_Nanodrop_ug<br/>DNA ≈ 0.75<br/>RNA ≈ 0.60"
          ),
          min = 0,
          max = 1,
          value = 0.75,
          step = 0.01
        ),
        
        # textInput("UTHRESH_TOTAL_Nanodrop_ug", label="Upper Threshold for TOTAL_Nanodrop_ug", value = Inf),
        # place holder for slider
        
        # textInput("LTHRESH_TOTAL_Qubit_ug", label="Lower Threshold for TOTAL_Qubit_ug", value = 0.75),
        sliderInput(
          "LTHRESH_TOTAL_Qubit_ug",
          label = HTML(
            "Lower Threshold for TOTAL_Qubit_ug<br/>DNA ≈ 0.75<br/>RNA ≈ 0.60"
          ),
          min = 0,
          max = 1,
          value = 0.75,
          step = 0.01
        ),
        
        # textInput("UTHRESH_TOTAL_Qubit_ug", label="Upper Threshold for TOTAL_Qubit_ug", value = Inf),
        # place holder for slider
        
        tags$h4("Missing Data"),
        checkboxInput("filterMissing", "Filter Missing", FALSE),
        
        actionButton("runScript", "Run QC & Generate Plots"),
        tags$h6(
          "Email christopher.arehart@cuanschutz.edu to report bugs or provide suggestions"
        )
            ),
      
      titlePanel("Summary Figures"),
      mainPanel(
        plotOutput("cGram", width = "auto", height = "520px"),
        uiOutput("mytable")
        # DT::dataTableOutput("mytable")
        
        # plotlyOutput("myGridPlots1"),
        # plotlyOutput("myGridPlots2"),
        # plotlyOutput("myGridPlots3"),
        # plotlyOutput("myGridPlots4"),
        # # plotlyOutput("myGridPlots5"),
        # plotlyOutput("myGridPlots6"),
        # plotlyOutput("myGridPlots7"),
        # plotlyOutput("myGridPlots8"),
        # # plotlyOutput("myGridPlots9"),
        # # plotlyOutput("myGridPlots10"),
        # plotlyOutput("myGridPlots11"),
        # plotlyOutput("myGridPlots12"),
        # # plotlyOutput("myGridPlots13"),
        # # plotlyOutput("myGridPlots14"),
        # # plotlyOutput("myGridPlots15"),
        # plotlyOutput("myGridPlots16")
      ),
      
      fluidRow(
        # column(5,
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots1")
        ),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots2")
        ),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots3")
        ),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots4")
        ),
        # plotlyOutput("myGridPlots5"),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots6")
        ),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots7")
        ),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots8")
        ),
        # plotlyOutput("myGridPlots9"),
        # plotlyOutput("myGridPlots10"),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots11")
        ),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots12")
        ),
        # plotlyOutput("myGridPlots13"),
        # plotlyOutput("myGridPlots14"),
        # plotlyOutput("myGridPlots15"),
        column(
          width = 6,
          offset = 0,
          plotlyOutput("myGridPlots16")
        )
        
        # plotlyOutput("myGridPlots2"),
        # plotlyOutput("myGridPlots3"),
        # plotlyOutput("myGridPlots4"),
        # # plotlyOutput("myGridPlots5"),
        # plotlyOutput("myGridPlots6"),
        # plotlyOutput("myGridPlots7"),
        # plotlyOutput("myGridPlots8"),
        # # plotlyOutput("myGridPlots9"),
        # # plotlyOutput("myGridPlots10"),
        # plotlyOutput("myGridPlots11"),
        # plotlyOutput("myGridPlots12"),
        # # plotlyOutput("myGridPlots13"),
        # # plotlyOutput("myGridPlots14"),
        # # plotlyOutput("myGridPlots15"),
        # plotlyOutput("myGridPlots16")
        # )
    )),
    

    # tabPanel(
    #   titlePanel("Plate Samples"),
    #   mainPanel(
    #     tags$h3("This is a tab for corey's plating script"),
    #     actionButton("plateSamples", "Reload Passed Samples"),
    #     downloadButton("downloadManifest", "Download Manifest"),
    #     # Corey, the passed samples table is saved in this Rdata file:
    #     # load("savePassed.RData"),
    #     DT::dataTableOutput("tableForCorey")
    #   )
    # )
    manifest_ui(input, output, session)
  )
))


server <- shinyServer(function(input, output, session) {
  mylist <- reactiveVal()
  
  # checkboxInput("includeBaltimore", "Baltimore", TRUE),
  # checkboxInput("includeBrazil", "Brazil", TRUE),
  # checkboxInput("includeChicago", "Chicago", TRUE),
  # checkboxInput("includeDenver", "Denver", TRUE),
  # checkboxInput("includeNigeria", "Nigeria", TRUE),
  # checkboxInput("includeWashington_DC", "Washington DC", TRUE),
  # checkboxInput("includeAdult", "Adults", TRUE),
  # checkboxInput("includeChild", "Children", TRUE),
  # Force Children Selection
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
  
  manifest_server(input, output, session)
})

shinyApp(ui, server)