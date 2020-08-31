library(shiny)
library(datasets)
# library(ggpairs)
require(GGally)
require(ggplot2)
library(tidyverse)
library(gridExtra)
library(plotly)
library(DT)

ui <- shinyUI(fluidPage(
  
                      titlePanel("QC metrics"),
                      tabsetPanel(
                        tabPanel(titlePanel("Data Setup"),
                         # selectInput("Fname", label="File Name",
                         #             c("./InputData/pbmc_rna_Rep1_data_for_plots.txt",
                         #               "./InputData/Nasal_RNA_Final_for_plots_062520.txt",
                         #               "./InputData/Nasal_DNA_Final_for_plots_062520.txt",
                         #               "./InputData/PBMC_RNA_Final_for_plots_062520.txt",
                         #               "./InputData/PBMC_DNA_Final_for_plots_062520.txt"), selected = "./InputData/pbmc_rna_Rep1_data_for_plots.txt"),
                                 
                        selectInput("Fname", label="File Name",
                                    c("./InputData/pbmc_rna_Rep1_data_for_plots.txt",
                                      "./InputData/pbmc_DNA_Rep1_data_for_plots.txt",
                                      "./InputData/Nasal_rna_Rep1_data_for_plots.txt",
                                      "./InputData/Nasal_DNA_Rep1_data_for_plots.txt"), selected = "./InputData/pbmc_rna_Rep1_data_for_plots.txt"),
                        # textInput("Fname", label="File Name", value = "./InputData/pbmc_rna_Rep1_data_for_plots.txt"),
                        selectInput("PBMC_PASS_ONLY", label = "Include only the samples that passed PBMC_status",
                                    c("TRUE", "FALSE"), selected = "TRUE"),
                        selectInput("Rep1_Slide_PASS_ONLY", label = "Include only the samples that passed Rep1_Slide",
                                    c("TRUE", "FALSE"), selected = "TRUE"),
                        
                        selectInput("Include_Boxplots", label = "If TRUE, the correlation plots will feture a 5th row and 5th column showing initial pass and fail based on PBMC_statys or Rep1_Slide",
                                    c("TRUE", "FALSE"), selected = "FALSE"),
                        
                        # sliderInput("UTHRESH_Rep1_X260.280", label = HTML("Upper Threshold for Rep1_X260.280,<br/>DNA ≈ 2.1 and RNA ≈ 2.2"), 
                        #             min = 1.8, max = 3, value = 2.10, step = 0.01),
                        # # textInput("UTHRESH_Rep1_X260.280", label=HTML("Upper Threshold for Rep1_X260.280,<br/>DNA ≈ 2.1 and RNA ≈ 2.2"), value = 2.2),
                        # sliderInput("LTHRESH_Rep1_X260.280", label = HTML("Lower Threshold for Rep1_X260.280,<br/>DNA ≈ 1.4 and RNA ≈ 1.7"), 
                        #             min = 1, max = 2, value = 1.40, step = 0.01),
                        # # textInput("LTHRESH_Rep1_X260.280", label=NULL, value = 1.7),
                        
                        tags$head( tags$style( type = "text/css", '
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
                                               
                                               ')), 
                        
                        tags$head( tags$style( type = "text/css", '
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
                                               
                                               ')), 
                        
                        tags$head( tags$style( type = "text/css", '
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
                                               
                                               ')), 
                        
                        tags$head( tags$style( type = "text/css", '
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
                                               
                                               ')), 
                        
                        sliderInput("LUTHRESH_Rep1_X260.280", label = HTML("Lower & Upper Thresholds for Rep1_X260.280,<br/>DNA ≈ 1.4 to 2.1<br/>RNA ≈ 1.7 to 2.2"), 
                                    min = 1, max = 3, value = c(1.40,2.10), step = 0.01),
                        
                        
                        # textInput("LTHRESH_Rep1_Agilent.DIN", label="Lower Threshold for Rep1_Agilent.DIN", value = 6),
                        sliderInput("LTHRESH_Rep1_Agilent.DIN", label = HTML("Lower Threshold for Rep1_Agilent.DIN"), 
                                    min = 0, max = 10, value = 6.00, step = 0.01),
                        
                        # textInput("UTHRESH_Rep1_Agilent.DIN", label="Upper Threshold for Rep1_Agilent.DIN", value = Inf),
                        # place holder for slider
                        
                        # textInput("LTHRESH_Rep1_Agilent.RINe", label="Lower Threshold for Rep1_Agilent.RINe", value = 6),
                        sliderInput("LTHRESH_Rep1_Agilent.RINe", label = HTML("Lower Threshold for Rep1_Agilent.RINe"), 
                                    min = 0, max = 10, value = 6.00, step = 0.01),
                        
                        # textInput("UTHRESH_Rep1_Agilent.RINe", label="Upper Threshold for Rep1_Agilent.RINe", value = Inf),
                        # place holder for slider
                        
                        # textInput("LTHRESH_Rep1_TOTAL_ug", label="Lower Threshold for Rep1_TOTAL_ug", value = 0.75),
                        sliderInput("LTHRESH_Rep1_TOTAL_ug", label = HTML("Lower Threshold for Rep1_TOTAL_ug<br/>DNA ≈ 0.75 and RNA ≈ 0.60"), 
                                    min = 0, max = 1, value = 0.75, step = 0.01),
                        
                        # textInput("UTHRESH_Rep1_TOTAL_ug", label="Upper Threshold for Rep1_TOTAL_ug", value = Inf),
                        # place holder for slider
                        
                        # textInput("LTHRESH_Rep1_TOTAL_Qubit_ug", label="Lower Threshold for Rep1_TOTAL_Qubit_ug", value = 0.75),
                        sliderInput("LTHRESH_Rep1_TOTAL_Qubit_ug", label = HTML("Lower Threshold for Rep1_TOTAL_Qubit_ug<br/>DNA ≈ 0.75 and RNA ≈ 0.60"), 
                                    min = 0, max =1, value = 0.75, step = 0.01),
                        
                        # textInput("UTHRESH_Rep1_TOTAL_Qubit_ug", label="Upper Threshold for Rep1_TOTAL_Qubit_ug", value = Inf),
                        # place holder for slider
                        
                        actionButton("runScript", "Run & Generate Plots")
                      )),
                      titlePanel("Summary Figures"),
                      mainPanel(
                        plotOutput("cGram"),
                        DT::dataTableOutput("mytable"),
                        plotlyOutput("myGridPlots1"),
                        plotlyOutput("myGridPlots2"),
                        plotlyOutput("myGridPlots3"),
                        plotlyOutput("myGridPlots4"),
                        # plotlyOutput("myGridPlots5"),
                        plotlyOutput("myGridPlots6"),
                        plotlyOutput("myGridPlots7"),
                        plotlyOutput("myGridPlots8"),
                        # plotlyOutput("myGridPlots9"),
                        # plotlyOutput("myGridPlots10"),
                        plotlyOutput("myGridPlots11"),
                        plotlyOutput("myGridPlots12"),
                        # plotlyOutput("myGridPlots13"),
                        # plotlyOutput("myGridPlots14"),
                        # plotlyOutput("myGridPlots15"),
                        plotlyOutput("myGridPlots16")
                        
                      )

                        )
              )


server <- shinyServer(function(input, output, session) {
  
  mylist <- reactiveVal()
  
  # autoUpdateSlider(textName = "LTHRESH_Rep1_X260.280", sliderName = "LTHRESH_Rep1_X260.280_slider")
  
  # observeEvent(input$LTHRESH_Rep1_X260.280, {
  #   print(input$LTHRESH_Rep1_X260.280)
  #   if ((as.numeric(input$LTHRESH_Rep1_X260.280) != input$LTHRESH_Rep1_X260.280_slider) &
  #       input$LTHRESH_Rep1_X260.280 != "" &  input$LTHRESH_Rep1_X260.280_slider != "")
  #   {
  #     updateSliderInput(
  #       session = session,
  #       inputId = 'LTHRESH_Rep1_X260.280_slider',
  #       value = input$LTHRESH_Rep1_X260.280
  #     )
  #   } else {
  #     if (input$LTHRESH_Rep1_X260.280 == "") {
  #       updateSliderInput(session = session,
  #                         inputId = 'LTHRESH_Rep1_X260.280_slider',
  #                         value = 0)
  #     }
  #   }
  # })
  # 
  # observeEvent(input$LTHRESH_Rep1_X260.280_slider, {
  #   if ((as.numeric(input$LTHRESH_Rep1_X260.280) != input$LTHRESH_Rep1_X260.280_slider) &
  #       input$LTHRESH_Rep1_X260.280_slider != "" & input$LTHRESH_Rep1_X260.280 != "")
  #   {
  #     updateTextInput(
  #       session = session,
  #       inputId = 'LTHRESH_Rep1_X260.280',
  #       value = input$LTHRESH_Rep1_X260.280_slider
  #     )
  # 
  #   }
  # })
  
  observe({
    mylist(list(
      Fname = as.character(input$Fname),
      PBMC_PASS_ONLY = as.character(input$PBMC_PASS_ONLY),
      Rep1_Slide_PASS_ONLY = as.character(input$Rep1_Slide_PASS_ONLY),
      Include_Boxplots = as.character(input$Include_Boxplots),
      # LTHRESH_Rep1_X260.280 = as.numeric(input$LTHRESH_Rep1_X260.280),
      # UTHRESH_Rep1_X260.280 = as.numeric(input$UTHRESH_Rep1_X260.280),
      LTHRESH_Rep1_X260.280 = as.numeric(input$LUTHRESH_Rep1_X260.280[1]),
      UTHRESH_Rep1_X260.280 = as.numeric(input$LUTHRESH_Rep1_X260.280[2]),
      LTHRESH_Rep1_Agilent.DIN = as.numeric(input$LTHRESH_Rep1_Agilent.DIN),
      # UTHRESH_Rep1_Agilent.DIN = as.numeric(input$UTHRESH_Rep1_Agilent.DIN),
      LTHRESH_Rep1_Agilent.RINe = as.numeric(input$LTHRESH_Rep1_Agilent.RINe),
      # UTHRESH_Rep1_Agilent.RINe = as.numeric(input$UTHRESH_Rep1_Agilent.RINe),
      LTHRESH_Rep1_TOTAL_ug = as.numeric(input$LTHRESH_Rep1_TOTAL_ug),
      # UTHRESH_Rep1_TOTAL_ug = as.numeric(input$UTHRESH_Rep1_TOTAL_ug),
      LTHRESH_Rep1_TOTAL_Qubit_ug = as.numeric(input$LTHRESH_Rep1_TOTAL_Qubit_ug)
      # UTHRESH_Rep1_TOTAL_Qubit_ug = as.numeric(input$UTHRESH_Rep1_TOTAL_Qubit_ug)
    ))
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
    }, width="auto", height = "auto")
    
    output$mytable <- DT::renderDataTable({
      summaryTable
    })
    
    output$myGridPlots1 <- renderPlotly({
      ggplotly(myGridPlots[[1]])
    })
    
    output$myGridPlots2 <- renderPlotly({
      ggplotly(myGridPlots[[2]])
    })

    output$myGridPlots3 <- renderPlotly({
      ggplotly(myGridPlots[[3]])
    })

    output$myGridPlots4 <- renderPlotly({
      ggplotly(myGridPlots[[4]])
    })

    # output$myGridPlots5 <- renderPlotly({
    #   ggplotly(myGridPlots[[5]])
    # })

    output$myGridPlots6 <- renderPlotly({
      ggplotly(myGridPlots[[6]])
    })

    output$myGridPlots7 <- renderPlotly({
      ggplotly(myGridPlots[[7]])
    })

    output$myGridPlots8 <- renderPlotly({
      ggplotly(myGridPlots[[8]])
    })

    # output$myGridPlots9 <- renderPlotly({
    #   ggplotly(myGridPlots[[9]])
    # })
    # 
    # output$myGridPlots10 <- renderPlotly({
    #   ggplotly(myGridPlots[[10]])
    # })

    output$myGridPlots11 <- renderPlotly({
      ggplotly(myGridPlots[[11]])
    })

    output$myGridPlots12 <- renderPlotly({
      ggplotly(myGridPlots[[12]])
    })

    # output$myGridPlots13 <- renderPlotly({
    #   ggplotly(myGridPlots[[13]])
    # })
    # 
    # output$myGridPlots14 <- renderPlotly({
    #   ggplotly(myGridPlots[[14]])
    # })
    # 
    # output$myGridPlots15 <- renderPlotly({
    #   ggplotly(myGridPlots[[15]])
    # })

    output$myGridPlots16 <- renderPlotly({
      ggplotly(myGridPlots[[16]])
    })
  })
  
  # output$MyPlot <- renderPlot({
  #
  #   x <- data()
  #   ggpairs(x, columns = c(6,12,18,19),lower = list(continuous = wrap("smooth",method="pearson",linetype="blank", alpha = 0.5, size=2)), ggplot2::aes(colour=Rep1_Slide))
  #
  # })
})

shinyApp(ui, server)