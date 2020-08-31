# setwd("/Users/chris/Documents/CARGO_QC")

# library(tidyverse)
# library(GGally)
# library(ggplot2)
# library(gridExtra)
# Fname = "./InputData/pbmc_rna_Rep1_data_for_plots.txt"
# PBMC_PASS_ONLY = "FALSE"
# Rep1_Slide_PASS_ONLY = "FALSE"
# LTHRESH_Rep1_X260.280 = 1.7
# UTHRESH_Rep1_X260.280 = 2.2
# LTHRESH_Rep1_Agilent.DIN = 6
# UTHRESH_Rep1_Agilent.DIN = Inf
# LTHRESH_Rep1_Agilent.RINe = 6
# UTHRESH_Rep1_Agilent.RINe = Inf
# LTHRESH_Rep1_TOTAL_ug = 0.75
# UTHRESH_Rep1_TOTAL_ug = Inf
# LTHRESH_Rep1_TOTAL_Qubit_ug = 0.75
# UTHRESH_Rep1_TOTAL_Qubit_ug = Inf

UTHRESH_Rep1_Agilent.RINe = Inf
UTHRESH_Rep1_Agilent.DIN = Inf
UTHRESH_Rep1_TOTAL_ug = Inf
UTHRESH_Rep1_TOTAL_Qubit_ug = Inf


# Load id original data
# F1 <- read.table("./InputData/pbmc_rna_Rep1_data_for_plots.txt", sep = "\t", header=T)
# F2 <- read.table("./InputData/pbmc_DNA_Rep1_data_for_plots.txt", sep = "\t", header=T)
# F3 <- read.table("./InputData/Nasal_rna_Rep1_data_for_plots.txt", sep = "\t", header=T)
# F4 <- read.table("./InputData/Nasal_DNA_Rep1_data_for_plots.txt", sep = "\t", header=T)

F1 <- read.table(Fname, sep = "\t", header=T)


# PBMC_PASS_ONLY = T
# Rep1_Slide_PASS_ONLY = T

# Subset data given choice of subset
dataFocus <- as_tibble(F1)

# Select meaningful columns
dataFocusSelect <-  dataFocus %>% select(starts_with("Rep1_TOTAL_ug"), 
                                         starts_with("Rep1_TOTAL_Qubit_ug"), 
                                         starts_with("Rep1_X260.280"), 
                                         starts_with("Rep1_Agilent.DIN"), 
                                         starts_with("Rep1_Agilent.RINe"),
                                         starts_with("PBMC_status"),
                                         starts_with("Rep1_Slide")
                                         )
glimpse(dataFocusSelect)
summary(dataFocusSelect)

c1 <- "red4"
c2 <- "coral"
my_density <- function(data,mapping,...){ggplot(data=data,mapping=mapping)+geom_density(...,lwd=1, colour = c1, fill = c2)}

cGramFunc <- function(inputData, colorColumn){
  if(is.null(colorColumn)){
    cGramPlot <- ggpairs(inputData, title="Correlation of QC Metrics", 
                         mapping=ggplot2::aes_string(alpha = 0.5),
                         lower = list(continuous = wrap("points", alpha = 0.5, colour = c1, size=0.2), 
                                      combo = wrap("dot", alpha = 0.5, colour = c1, size=0.2) ),
                         diag = list(continuous = my_density)
    )
  }else{
    cGramPlot <- ggpairs(inputData, title="Correlation of QC Metrics", 
                         mapping=ggplot2::aes_string(colour = colorColumn, alpha = 0.5),
                         lower = list(continuous = wrap("points", alpha = 0.5,    size=0.2), 
                                      combo = wrap("dot", alpha = 0.5,            size=0.2) )
    )
  }
  
  
  
  index_Rep1_TOTAL_ug <- which(names(inputData) == "Rep1_TOTAL_ug")
  index_Rep1_TOTAL_Qubit_ug <- which(names(inputData) == "Rep1_TOTAL_Qubit_ug")
  index_Rep1_X260.280 <- which(names(inputData) == "Rep1_X260.280")
  index_Rep1_Agilent.DIN <- which(names(inputData) == "Rep1_Agilent.DIN")
  index_Rep1_Agilent.RINe <- which(names(inputData) == "Rep1_Agilent.RINe")
  
  toKeep <- which(c("Rep1_TOTAL_ug", "Rep1_TOTAL_Qubit_ug", "Rep1_X260.280", "Rep1_Agilent.DIN",
                    "Rep1_Agilent.RINe")
                  %in% cGramPlot$xAxisLabels)
  
  indexVec <- c(index_Rep1_TOTAL_ug, index_Rep1_TOTAL_Qubit_ug, index_Rep1_X260.280, index_Rep1_Agilent.DIN,
                index_Rep1_Agilent.RINe)
  
  LTHRESHVec <- c(LTHRESH_Rep1_TOTAL_ug, LTHRESH_Rep1_TOTAL_Qubit_ug, LTHRESH_Rep1_X260.280, LTHRESH_Rep1_Agilent.DIN,
                  LTHRESH_Rep1_Agilent.RINe)[toKeep]
  
  UTHRESH_Vec <- c(UTHRESH_Rep1_TOTAL_ug, UTHRESH_Rep1_TOTAL_Qubit_ug, UTHRESH_Rep1_X260.280, UTHRESH_Rep1_Agilent.DIN,
                  UTHRESH_Rep1_Agilent.RINe)[toKeep]
  
  thresholdColors1 <- "blue"
  thresholdSize <- 0.55
  
  for(i in 1:length(indexVec)){
    for(j in 1:length(indexVec)){
      if(i > j){
        cGramPlot[indexVec[i],indexVec[j]] = cGramPlot[indexVec[i],indexVec[j]] + 
          geom_vline(xintercept = LTHRESHVec[j],color=thresholdColors1,size=thresholdSize)+
          geom_vline(xintercept = UTHRESH_Vec[j],color=thresholdColors1,size=thresholdSize)+
          geom_hline(yintercept = LTHRESHVec[i],color=thresholdColors1,size=thresholdSize)+
          geom_hline(yintercept = UTHRESH_Vec[i],color=thresholdColors1,size=thresholdSize)
      }else if(i == j){
        cGramPlot[indexVec[i],indexVec[j]] = cGramPlot[indexVec[i],indexVec[j]] + 
          geom_vline(xintercept = LTHRESHVec[i],color=thresholdColors1,size=thresholdSize)+
          geom_vline(xintercept = UTHRESH_Vec[i],color=thresholdColors1,size=thresholdSize)
      }
      # if(is.null(colorColumn)){
      #   cGramPlot[indexVec[i],indexVec[j]] = cGramPlot[indexVec[i],indexVec[j]] + 
      #     scale_fill_manual(values=c("green")) +
      #     scale_color_manual(values=c("coral"))  
      # }
    }
  }
  
  return(cGramPlot)
}

if("PBMC_status" %in% names(dataFocusSelect)){
  if(PBMC_PASS_ONLY == "FALSE"){
    dataFocusSelectSub <- subset(dataFocusSelect, PBMC_status %in% c("PASS","FAIL"))
    if(Include_Boxplots == "FALSE"){
      dataFocusSelectSub$PBMC_status <- NULL
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = NULL)
    }else{
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "PBMC_status")
    }
  }else if(PBMC_PASS_ONLY == "TRUE"){
    dataFocusSelectSub <- subset(dataFocusSelect, PBMC_status %in% c("PASS"))
    if(Include_Boxplots == "FALSE"){
      dataFocusSelectSub$PBMC_status <- NULL
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = NULL)
    }else{
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "PBMC_status")
    }
  }
}else if("Rep1_Slide" %in% names(dataFocusSelect)){
  if(Rep1_Slide_PASS_ONLY == "FALSE"){
    dataFocusSelectSub <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS","FAIL"))
    if(Include_Boxplots == "FALSE"){
      dataFocusSelectSub$Rep1_Slide <- NULL
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = NULL)
    }else{
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "Rep1_Slide")
    }
  }else if(Rep1_Slide_PASS_ONLY == "TRUE"){
    dataFocusSelectSub <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS"))
    if(Include_Boxplots == "FALSE"){
      dataFocusSelectSub$Rep1_Slide <- NULL
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = NULL)
    }else{
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "Rep1_Slide")
    }
  }
}

cGram <- cGram + theme(legend.position = "none", 
                         panel.grid.major = element_blank(), 
                         axis.ticks = element_blank(), 
                         panel.border = element_rect(linetype = "solid", colour = "white", fill = NA))



# print(cGram)

thresholdPlotter <- function(inputData, var1, var2, thresh1_lower, thresh1_upper, thresh2_lower, thresh2_upper){
  inputData$var1Pass <- "FAIL"
  inputData$var1Pass[inputData[[var1]] >= thresh1_lower & inputData[[var1]] <= thresh1_upper] <- "PASS"
  inputData$var2Pass <- "FAIL"
  inputData$var2Pass[inputData[[var2]] >= thresh2_lower & inputData[[var2]] <= thresh2_upper] <- "PASS"
  inputData$var12Pass <- "FAIL"
  inputData$var12Pass[inputData$var1Pass == "PASS" & inputData$var2Pass == "PASS"] <- "PASS"
  
  inputData$var12Pass <- as.character(inputData$var12Pass)
  inputData$FAIL_ANY <- as.character(inputData$FAIL_ANY)
  
  inputData$aesCol <- "PASSES ALL THESHOLDS"
  inputData$aesCol[inputData$FAIL_ANY == "FAIL ANY" & inputData$var12Pass == "PASS"] <- "FAILS OTHER THRESHOLD"
  inputData$aesCol[inputData$var12Pass == "FAIL"] <- "FAILS THESE THRESHOLDS"
  
  inputData$aesCol[inputData$aesCol == "PASSES ALL THESHOLDS"] <- paste0("PASSES ALL THESHOLDS, N = ",length(inputData$aesCol[inputData$aesCol == "PASSES ALL THESHOLDS"]))
  inputData$aesCol[inputData$aesCol == "FAILS OTHER THRESHOLD"] <- paste0("FAILS OTHER THRESHOLD, N = ",length(inputData$aesCol[inputData$aesCol == "FAILS OTHER THRESHOLD"]))
  inputData$aesCol[inputData$aesCol == "FAILS THESE THRESHOLDS"] <- paste0("FAILS THESE THRESHOLDS, N = ",length(inputData$aesCol[inputData$aesCol == "FAILS THESE THRESHOLDS"]))
  inputData$aesCol <- as.factor(inputData$aesCol)
  
  thresholdColors2 <- "blue"
  
  myPlot <- ggplot()+
    geom_point(data = inputData, aes_string(x = var1, y = var2, color = "aesCol"), size = 1)+
    # geom_point(data = inputData, aes_string(x = var1, y = var2, color = "var12Pass", shape = "FAIL_ANY"), size = 1)+
    geom_vline(xintercept = thresh1_lower, color = thresholdColors2)+
    geom_vline(xintercept = thresh1_upper, color = thresholdColors2)+
    geom_hline(yintercept = thresh2_lower, color = thresholdColors2)+
    geom_hline(yintercept = thresh2_upper, color = thresholdColors2)+
    # scale_color_manual(NULL, breaks = c("FAIL", "PASS"), labels = c("FAIL", "PASS"), values = c("firebrick3", "forestgreen"))+
    scale_color_manual(paste0("Total N = ",length(inputData$aesCol)), 
                       # labels = c(paste0("FAILS THESE THRESHOLDS, N = ",length(inputData$aesCol[inputData$aesCol == "FAILS THESE THRESHOLDS"])), 
                       #            paste0("FAILS OTHER THRESHOLD, N = ",length(inputData$aesCol[inputData$aesCol == "FAILS OTHER THRESHOLD"])),  
                       #            paste0("PASSES ALL THESHOLDS, N = ",length(inputData$aesCol[inputData$aesCol == "PASSES ALL THESHOLDS"]))
                       #            ), 
                       # breaks = c("FAILS THESE THRESHOLDS", "FAILS OTHER THRESHOLD", "PASSES ALL THESHOLDS"), 
                       values = c("orange", "firebrick3", "forestgreen"))+
    # scale_shape_manual(NULL, breaks = c("FAIL ANY", "PASS ALL"), labels = c("FAIL ANY", "PASS ALL"), values = c(4, 16))
  # scale_fill_manual(labels = c("PASS ALL", "FAIL ANY"), values = c("black", "pink"))+
  theme_bw(base_size = 8)
  return(myPlot)
}

gridMaker <- function(inputData, VarVec){
  
  myplots <- vector('list', length(VarVec)^2)
  cc <- 1
  for(i in 1:length(VarVec)){
  # for(i in 1:2){
    for(j in 1:length(VarVec)){
    # for(j in 1:1){
      myplots[[cc]] <- local({
        i <- i
        j <- j
        idat <- dataFocusSelectSubContinue %>% select(VarVec[i], VarVec[j], FAIL_ANY)
        p1 <- thresholdPlotter(inputData = idat, 
                               var1 = VarVec[i], var2 = VarVec[j], 
                               thresh1_lower = eval(parse(text = paste0("LTHRESH_",VarVec[i]))), thresh1_upper = eval(parse(text = paste0("UTHRESH_",VarVec[i]))), 
                               thresh2_lower = eval(parse(text = paste0("LTHRESH_",VarVec[j]))), thresh2_upper = eval(parse(text = paste0("UTHRESH_",VarVec[j]))))
        p1
      })
      cc <- cc+1
      }
  }
  
  return(myplots)

}


# SET THE THRESHOLDS ----
pdf("testing.pdf", width = 20, height = 15)
if("Rep1_Agilent.DIN" %in% names(dataFocusSelect)){
  # DNA THRESHOLDS
  # LTHRESH_Rep1_X260.280 <- 1.4
  # UTHRESH_Rep1_X260.280 <- 2.1
  # LTHRESH_Rep1_Agilent.DIN <- 6
  # UTHRESH_Rep1_Agilent.DIN <- Inf
  # LTHRESH_Rep1_TOTAL_ug <- 0.75
  # UTHRESH_Rep1_TOTAL_ug <- Inf
  # LTHRESH_Rep1_TOTAL_Qubit_ug <- 0.75
  # UTHRESH_Rep1_TOTAL_Qubit_ug <- Inf
  
  # Subset to those that pass
  if("PBMC_status" %in% names(dataFocusSelect)){
    if(PBMC_PASS_ONLY == "FALSE"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_status %in% c("PASS","FAIL"))
    }else if(PBMC_PASS_ONLY == "TRUE"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_status %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- "PASS ALL"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_X260.280"]] <= LTHRESH_Rep1_X260.280 | dataFocusSelectSubContinue[["Rep1_X260.280"]] >= UTHRESH_Rep1_X260.280 ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_Agilent.DIN"]] <= LTHRESH_Rep1_Agilent.DIN | dataFocusSelectSubContinue[["Rep1_Agilent.DIN"]] >= UTHRESH_Rep1_Agilent.DIN ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] <= LTHRESH_Rep1_TOTAL_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] >= UTHRESH_Rep1_TOTAL_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] <= LTHRESH_Rep1_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] >= UTHRESH_Rep1_TOTAL_Qubit_ug ] <- "FAIL ANY"
    # dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
    # run grid plotter
    myGridPlots <- gridMaker(inputData = dataFocusSelectSubContinue, VarVec = c("Rep1_X260.280", "Rep1_Agilent.DIN", "Rep1_TOTAL_ug", "Rep1_TOTAL_Qubit_ug") )
    plotDF <- data.frame(matrix(unlist(1:4^2), nrow=4, byrow=T),stringsAsFactors=FALSE)
    plotDFmat <- as.matrix(plotDF)
    gridPlot <- grid.arrange(grobs = myGridPlots, layout_matrix = plotDFmat)
    print(gridPlot)
    
    summaryTable <- data.frame(matrix(ncol = 0, nrow = 3))
    summaryTable$`_____` <- c("N", "PASS ALL", "FAIL ANY")
    summaryTable$Total <- c(nrow(dataFocusSelectSubContinue),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Pass_PBMC_status <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "PASS")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Fail_PBMC_status <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "FAIL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    # summaryTable
    
  }else if("Rep1_Slide" %in% names(dataFocusSelect)){
    if(Rep1_Slide_PASS_ONLY == "FALSE"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS","FAIL"))
    }else if(Rep1_Slide_PASS_ONLY == "TRUE"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- "PASS ALL"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_X260.280"]] <= LTHRESH_Rep1_X260.280 | dataFocusSelectSubContinue[["Rep1_X260.280"]] >= UTHRESH_Rep1_X260.280 ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_Agilent.DIN"]] <= LTHRESH_Rep1_Agilent.DIN | dataFocusSelectSubContinue[["Rep1_Agilent.DIN"]] >= UTHRESH_Rep1_Agilent.DIN ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] <= LTHRESH_Rep1_TOTAL_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] >= UTHRESH_Rep1_TOTAL_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] <= LTHRESH_Rep1_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] >= UTHRESH_Rep1_TOTAL_Qubit_ug ] <- "FAIL ANY"
    # dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
    # run grid plotter
    myGridPlots <- gridMaker(inputData = dataFocusSelectSubContinue, VarVec = c("Rep1_X260.280", "Rep1_Agilent.DIN", "Rep1_TOTAL_ug", "Rep1_TOTAL_Qubit_ug") )
    plotDF <- data.frame(matrix(unlist(1:4^2), nrow=4, byrow=T),stringsAsFactors=FALSE)
    plotDFmat <- as.matrix(plotDF)
    gridPlot <- grid.arrange(grobs = myGridPlots, layout_matrix = plotDFmat)
    print(gridPlot)
    
    summaryTable <- data.frame(matrix(ncol = 0, nrow = 3))
    summaryTable$`_____` <- c("N", "PASS ALL", "FAIL ANY")
    summaryTable$Total <- c(nrow(dataFocusSelectSubContinue),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Pass_Rep1_Slide <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "PASS")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Fail_Rep1_Slide <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "FAIL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    # summaryTable
  }
}else if("Rep1_Agilent.RINe" %in% names(dataFocusSelect)){
  # RNA THRESHOLDS
  # LTHRESH_Rep1_X260.280 <- 1.7
  # UTHRESH_Rep1_X260.280 <- 2.2
  # LTHRESH_Rep1_Agilent.RINe <- 6
  # UTHRESH_Rep1_Agilent.RINe <- Inf
  # LTHRESH_Rep1_TOTAL_ug <- 0.75
  # UTHRESH_Rep1_TOTAL_ug <- Inf
  # LTHRESH_Rep1_TOTAL_Qubit_ug <- 0.75
  # UTHRESH_Rep1_TOTAL_Qubit_ug <- Inf
  
  # Subset to those that pass
  if("PBMC_status" %in% names(dataFocusSelect)){
    if(PBMC_PASS_ONLY == "FALSE"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_status %in% c("PASS","FAIL"))
    }else if(PBMC_PASS_ONLY == "TRUE"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_status %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- "PASS ALL"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_X260.280"]] <= LTHRESH_Rep1_X260.280 | dataFocusSelectSubContinue[["Rep1_X260.280"]] >= UTHRESH_Rep1_X260.280 ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_Agilent.RINe"]] <= LTHRESH_Rep1_Agilent.RINe | dataFocusSelectSubContinue[["Rep1_Agilent.RINe"]] >= UTHRESH_Rep1_Agilent.RINe ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] <= LTHRESH_Rep1_TOTAL_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] >= UTHRESH_Rep1_TOTAL_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] <= LTHRESH_Rep1_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] >= UTHRESH_Rep1_TOTAL_Qubit_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
    # run grid plotter
    myGridPlots <- gridMaker(inputData = dataFocusSelectSubContinue, VarVec = c("Rep1_X260.280", "Rep1_Agilent.RINe", "Rep1_TOTAL_ug", "Rep1_TOTAL_Qubit_ug") )
    plotDF <- data.frame(matrix(unlist(1:4^2), nrow=4, byrow=T),stringsAsFactors=FALSE)
    plotDFmat <- as.matrix(plotDF)
    gridPlot <- grid.arrange(grobs = myGridPlots, layout_matrix = plotDFmat)
    print(gridPlot)
    
    summaryTable <- data.frame(matrix(ncol = 0, nrow = 3))
    summaryTable$`_____` <- c("N", "PASS ALL", "FAIL ANY")
    summaryTable$Total <- c(nrow(dataFocusSelectSubContinue),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Pass_PBMC_status <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "PASS")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Fail_PBMC_status <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "FAIL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_status == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    # summaryTable
    
  }else if("Rep1_Slide" %in% names(dataFocusSelect)){
    if(Rep1_Slide_PASS_ONLY == "FALSE"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS","FAIL"))
    }else if(Rep1_Slide_PASS_ONLY == "TRUE"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- "PASS ALL"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_X260.280"]] <= LTHRESH_Rep1_X260.280 | dataFocusSelectSubContinue[["Rep1_X260.280"]] >= UTHRESH_Rep1_X260.280 ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_Agilent.RINe"]] <= LTHRESH_Rep1_Agilent.RINe | dataFocusSelectSubContinue[["Rep1_Agilent.RINe"]] >= UTHRESH_Rep1_Agilent.RINe ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] <= LTHRESH_Rep1_TOTAL_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] >= UTHRESH_Rep1_TOTAL_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] <= LTHRESH_Rep1_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] >= UTHRESH_Rep1_TOTAL_Qubit_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
    # run grid plotter
    myGridPlots <- gridMaker(inputData = dataFocusSelectSubContinue, VarVec = c("Rep1_X260.280", "Rep1_Agilent.RINe", "Rep1_TOTAL_ug", "Rep1_TOTAL_Qubit_ug") )
    plotDF <- data.frame(matrix(unlist(1:4^2), nrow=4, byrow=T),stringsAsFactors=FALSE)
    plotDFmat <- as.matrix(plotDF)
    gridPlot <- grid.arrange(grobs = myGridPlots, layout_matrix = plotDFmat)
    print(gridPlot)
    
    summaryTable <- data.frame(matrix(ncol = 0, nrow = 3))
    summaryTable$`_____` <- c("N", "PASS ALL", "FAIL ANY")
    summaryTable$Total <- c(nrow(dataFocusSelectSubContinue),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Pass_Rep1_Slide <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "PASS")),
                                      nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                      nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Fail_Rep1_Slide <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "FAIL")),
                                      nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                      nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Rep1_Slide == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    # summaryTable
  }
}

dev.off()

rownames(summaryTable) <- summaryTable$`_____`
summaryTable$`_____` <- NULL
print(summaryTable)

save(summaryTable, cGram, myGridPlots, gridPlot, file="savePlots.RData")

print("SCRIPT HAS FINISHED")


print(PBMC_PASS_ONLY)
print(Rep1_Slide_PASS_ONLY)




