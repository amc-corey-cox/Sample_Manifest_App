upload_ctrl <- fileInput("files", "Upload Phenotype File", multiple = FALSE,
  # fileInput("files", "Upload Phenotype File", multiple = TRUE,
  accept = c(".csv", ".xls", ".xlsx", "text/csv",
             "text/comma-separated-values,text/plain",
             "application/vnd.ms-excel",
             "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"))

file_ctrl <- conditionalPanel( condition = "output.fileUploaded",
                    checkboxInput("col_names", "Column names in header", TRUE),
                    numericInput("skip", "Skip lines:", 0, min = 0))

excel_ctrl <- conditionalPanel( condition = "output.fileExcel",
                    selectInput("sheet", "Sheet", choices = NULL))

txt_ctrl <- conditionalPanel(
  condition = "output.fileUploaded && ! output.fileExcel",
  radioButtons("delim", "Separator", selected = "\t",
    choices = c(Comma = ",", Semicolon = ";", Tab = "\t")),
  radioButtons("quote", "Quote", selected = '"',
    choices = c(None = "", "Double Quote" = '"', "Single Quote" = "'")))

id_ctrl <- conditionalPanel(
  condition = "output.fileUploaded",
  # radioButtons("id_radio", "Sample ID Column", choices = c(1))
  selectInput("d_id_col", "ID Column Name", choices = NULL))

data_panel <- div(
  style = "min-width: 240px; max-width: 240px;",
  wellPanel(
    conditionalPanel(
      condition = "input.d_tabs == 'Phenotype'",
      upload_ctrl,
      file_ctrl,
      excel_ctrl,
      txt_ctrl),
    conditionalPanel(condition = "input.d_tabs == 'Phenotype' || input.d_tabs == 'Cleanup'", id_ctrl),
    conditionalPanel(
      condition = "input.d_tabs == 'Cleanup'",
      uiOutput("ui_clean_cols")
)))

data_tabs <- div(
  style="display:inline-block; min-width: 400px; padding-left:25px; padding-top:10px",
  tabsetPanel(id = "d_tabs", type = "pills",
    # tabPanel( "Phenotype", uiOutput('uploadUI') ),
    tabPanel("Phenotype", tableOutput("pheno")),
    tabPanel("Cleanup", tableOutput("cleaned_pheno"))
    # tabPanel("Filter", tableOutput("filtered_pheno"))
))

data_ui <- tabPanel(titlePanel("Import Data"),
  splitLayout(
    cellWidths = c("55%", "20%", "25%"),
    titlePanel("Load and Prepare Data"),
    radioButtons("d_disp", "Display", selected = "head", inline = TRUE, choices = c(Head = "head", All = "all")),
    downloadButton("downloadCleanData", "Get Cleaned Data"),
    tags$style(type = 'text/css', "#d_disp { margin-top: 15px; }"),
    tags$style(type = 'text/css', "#downloadCleanData { margin-top: 20px; }")),
  flowLayout(data_panel, data_tabs))