setwd("/Users/chris/Documents/CARGO_QC")

library(tidyverse)
library(GGally)
library(ggplot2)
library(gridExtra)

# Load id original data
F1 <- read.table("./InputData/pbmc_rna_Rep1_data_for_plots.txt", sep = "\t", header=T)
F2 <- read.table("./InputData/pbmc_DNA_Rep1_data_for_plots.txt", sep = "\t", header=T)
F3 <- read.table("./InputData/Nasal_rna_Rep1_data_for_plots.txt", sep = "\t", header=T)
F4 <- read.table("./InputData/Nasal_DNA_Rep1_data_for_plots.txt", sep = "\t", header=T)

PBMC_PASS_ONLY = T
Rep1_Slide_PASS_ONLY = T

# Subset data given choice of subset
dataFocus <- as_tibble(F2)

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

cGramFunc <- function(inputData, colorColumn){
  cGramPlot <- ggpairs(inputData, title="Correlation of QC Metrics", 
                   mapping=ggplot2::aes_string(colour = colorColumn, alpha = 0.5),
                   lower = list(continuous = wrap("points", alpha = 0.5,    size=0.2), 
                                combo = wrap("dot", alpha = 0.5,            size=0.2) )
  )
  return(cGramPlot)
}

if("PBMC_status" %in% names(dataFocusSelect)){
  if(PBMC_PASS_ONLY == F){
    dataFocusSelectSub <- subset(dataFocusSelect, PBMC_status %in% c("PASS","FAIL"))
    cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "PBMC_status")
  }else if(PBMC_PASS_ONLY == T){
    dataFocusSelectSub <- subset(dataFocusSelect, PBMC_status %in% c("PASS"))
    cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "PBMC_status")
  }
}else if("Rep1_Slide" %in% names(dataFocusSelect)){
  if(Rep1_Slide_PASS_ONLY == F){
    dataFocusSelectSub <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS","FAIL"))
    cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "Rep1_Slide")
  }else if(Rep1_Slide_PASS_ONLY == T){
    dataFocusSelectSub <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS"))
    cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "Rep1_Slide")
  }
}

cGram <- cGram + theme(legend.position = "none", 
                         panel.grid.major = element_blank(), 
                         axis.ticks = element_blank(), 
                         panel.border = element_rect(linetype = "solid", colour = "white", fill = NA))

print(cGram)

thresholdPlotter <- function(inputData, var1, var2, thresh1_lower, thresh1_upper, thresh2_lower, thresh2_upper){
  inputData$var1Pass <- "FAIL"
  inputData$var1Pass[inputData[[var1]] >= thresh1_lower & inputData[[var1]] <= thresh1_upper] <- "PASS"
  inputData$var2Pass <- "FAIL"
  inputData$var2Pass[inputData[[var2]] >= thresh2_lower & inputData[[var2]] <= thresh2_upper] <- "PASS"
  inputData$var12Pass <- "FAIL"
  inputData$var12Pass[inputData$var1Pass == "PASS" & inputData$var2Pass == "PASS"] <- "PASS"
  
  inputData$var12Pass <- as.factor(inputData$var12Pass)
  inputData$FAIL_ANY <- as.factor(inputData$FAIL_ANY)
  
  myPlot <- ggplot()+
    geom_point(data = inputData, aes_string(x = var1, y = var2, color = "var12Pass", shape = "FAIL_ANY"), size = .7)+
    geom_vline(xintercept = thresh1_lower, color = "darkorchid4")+
    geom_vline(xintercept = thresh1_upper, color = "darkorchid4")+
    geom_hline(yintercept = thresh2_lower, color = "lightblue4")+
    geom_hline(yintercept = thresh2_upper, color = "lightblue4")+
    scale_color_manual(NULL, breaks = c("FAIL", "PASS"), labels = c("FAIL", "PASS"), values = c("firebrick3", "forestgreen"))+
    scale_shape_manual(NULL, breaks = c("FAIL ANY", "PASS ALL"), labels = c("FAIL ANY", "PASS ALL"), values = c(4, 16))
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
  LTHRESH_Rep1_X260.280 <- 1.4
  UTHRESH_Rep1_X260.280 <- 2.1
  LTHRESH_Rep1_Agilent.DIN <- 6
  UTHRESH_Rep1_Agilent.DIN <- Inf
  LTHRESH_Rep1_TOTAL_ug <- 0.75
  UTHRESH_Rep1_TOTAL_ug <- Inf
  LTHRESH_Rep1_TOTAL_Qubit_ug <- 0.75
  UTHRESH_Rep1_TOTAL_Qubit_ug <- Inf
  
  # Subset to those that pass
  if("PBMC_status" %in% names(dataFocusSelect)){
    if(PBMC_PASS_ONLY == F){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_status %in% c("PASS","FAIL"))
    }else if(PBMC_PASS_ONLY == T){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_status %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- "PASS ALL"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_X260.280"]] <= LTHRESH_Rep1_X260.280 | dataFocusSelectSubContinue[["Rep1_X260.280"]] >= UTHRESH_Rep1_X260.280 ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_Agilent.DIN"]] <= LTHRESH_Rep1_Agilent.DIN | dataFocusSelectSubContinue[["Rep1_Agilent.DIN"]] >= UTHRESH_Rep1_Agilent.DIN ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] <= LTHRESH_Rep1_TOTAL_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] >= UTHRESH_Rep1_TOTAL_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] <= LTHRESH_Rep1_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] >= UTHRESH_Rep1_TOTAL_Qubit_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
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
    if(Rep1_Slide_PASS_ONLY == F){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS","FAIL"))
    }else if(Rep1_Slide_PASS_ONLY == T){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- "PASS ALL"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_X260.280"]] <= LTHRESH_Rep1_X260.280 | dataFocusSelectSubContinue[["Rep1_X260.280"]] >= UTHRESH_Rep1_X260.280 ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_Agilent.DIN"]] <= LTHRESH_Rep1_Agilent.DIN | dataFocusSelectSubContinue[["Rep1_Agilent.DIN"]] >= UTHRESH_Rep1_Agilent.DIN ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] <= LTHRESH_Rep1_TOTAL_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_ug"]] >= UTHRESH_Rep1_TOTAL_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] <= LTHRESH_Rep1_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["Rep1_TOTAL_Qubit_ug"]] >= UTHRESH_Rep1_TOTAL_Qubit_ug ] <- "FAIL ANY"
    dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
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
  LTHRESH_Rep1_X260.280 <- 1.7
  UTHRESH_Rep1_X260.280 <- 2.2
  LTHRESH_Rep1_Agilent.RINe <- 6
  UTHRESH_Rep1_Agilent.RINe <- Inf
  LTHRESH_Rep1_TOTAL_ug <- 0.75
  UTHRESH_Rep1_TOTAL_ug <- Inf
  LTHRESH_Rep1_TOTAL_Qubit_ug <- 0.75
  UTHRESH_Rep1_TOTAL_Qubit_ug <- Inf
  
  # Subset to those that pass
  if("PBMC_status" %in% names(dataFocusSelect)){
    if(PBMC_PASS_ONLY == F){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_status %in% c("PASS","FAIL"))
    }else if(PBMC_PASS_ONLY == T){
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
    if(Rep1_Slide_PASS_ONLY == F){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Rep1_Slide %in% c("PASS","FAIL"))
    }else if(Rep1_Slide_PASS_ONLY == T){
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

summaryTable





