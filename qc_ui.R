qc_panel <- div(
  style = "min-width: 240px; max-width: 240px;",
  wellPanel(
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
    
    
    sliderInput(
      "LTHRESH_Agilent_DIN",
      label = HTML("Lower Threshold for Agilent_DIN"),
      min = 0,
      max = 10,
      value = 6.00,
      step = 0.01
    ),
    
    sliderInput(
      "LTHRESH_Agilent_RINe",
      label = HTML("Lower Threshold for Agilent_RINe"),
      min = 0,
      max = 10,
      value = 6.00,
      step = 0.01
    ),
    
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
    
    tags$h4("Missing Data"),
    checkboxInput("filterMissing", "Filter Missing", FALSE),
    
    actionButton("runScript", "Run QC & Generate Plots"),
    downloadButton("qcReport", "Download QC Report"),
    tags$h6(
      "Email Corey.Cox@CUAnschutz.edu to report bugs or provide suggestions"
)))
  
  ### TODO: Move figures and tables to live figures within app
qc_main <- div(
  style="display:inline-block; min-width: 400px; padding-left:25px; padding-top:10px",
    titlePanel("Quality Control Figures"),
    plotOutput("cGram", width = "auto", height = "520px"),
    uiOutput("mytable"),
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
      column(
        width = 6,
        offset = 0,
        plotlyOutput("myGridPlots16")
      )
    )
  )

qc_ui <- tabPanel(titlePanel("Quality Control"), flowLayout(qc_panel, qc_main))
