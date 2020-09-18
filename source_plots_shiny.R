# setwd("/Users/chris/Documents/CARGO_QC")
#
library(tidyverse)
library(GGally)
library(ggplot2)
library(gridExtra)
library(reshape2)

# CellType = "Nasal"
# DNAorRNA = "DNA"
# PBMC_PASS_ONLY = "Yes"
# Slide_status_PASS_ONLY = "Yes"
# Include_Boxplots = "No"
# LTHRESH_Nanodrop_260_280 = 1.7
# UTHRESH_Nanodrop_260_280 = 2.2
# LTHRESH_Agilent_DIN = 6
# UTHRESH_Agilent_DIN = Inf
# LTHRESH_Agilent_RINe = 6
# UTHRESH_Agilent_RINe = Inf
# LTHRESH_TOTAL_Nanodrop_ug = 0.75
# UTHRESH_TOTAL_Nanodrop_ug = Inf
# LTHRESH_TOTAL_Qubit_ug = 0.75
# UTHRESH_TOTAL_Qubit_ug = Inf
# includeBaltimore = T
# includeBrazil = T
# includeChicago = T
# includeDenver = T
# includeNigeria = T
# includeWashington_DC = T
# includeAdult = T
# includeChild = T
# includeMales = T
# includeFemales = T
# includeAsthma = T
# includeNonAsthma = T

UTHRESH_Agilent_RINe = Inf
UTHRESH_Agilent_DIN = Inf
UTHRESH_TOTAL_Nanodrop_ug = Inf
UTHRESH_TOTAL_Qubit_ug = Inf

# Load id original data
# F1 <- read.table("./InputData/pbmc_rna_Rep1_data_for_plots.txt", sep = "\t", header=T)
# F2 <- read.table("./InputData/pbmc_DNA_Rep1_data_for_plots.txt", sep = "\t", header=T)
# F3 <- read.table("./InputData/Nasal_rna_Rep1_data_for_plots.txt", sep = "\t", header=T)
# F4 <- read.table("./InputData/Nasal_DNA_Rep1_data_for_plots.txt", sep = "\t", header=T)

nameDate <- "072320_updated"
inputDir <- "./InputData/"
inputDir <- if_else(filterMissing, "./FilteredData/", "./InputData/")
# inputDir
# nameDate <- "062520"

# nameDate <- "091720"
# inputDir <- "./InputData2/"
# inputDir <- if_else(filterMissing, "./FilteredData/", "./InputData2/")

if(CellType == "PBMC" & DNAorRNA == "DNA"){
  Fname = paste0(inputDir,"PBMC_DNA_Final_for_plots_",nameDate,".csv")
}else if(CellType == "PBMC" & DNAorRNA == "RNA"){
  Fname = paste0(inputDir,"PBMC_RNA_Final_for_plots_",nameDate,".csv")
}else if(CellType == "Nasal" & DNAorRNA == "DNA"){
  Fname = paste0(inputDir,"Nasal_DNA_Final_for_plots_",nameDate,".csv")
}else if(CellType == "Nasal" & DNAorRNA == "RNA"){
  Fname = paste0(inputDir,"Nasal_RNA_Final_for_plots_",nameDate,".csv")
}

F1 <- read.table(Fname, sep = ",", header=T, na.strings = c("",NA))
# F1_all <- read.table(Fname, sep = ",", header=T, na.strings = c("",NA))
# F1 <- F1_all %>% select(subject, site, type, Age_category, Asthma, Gender, TOTAL_Nanodrop_ug,
#                         TOTAL_Qubit_ug, Nanodrop_260_280, Agilent_RINe, Slide_status)
if("Agilent_DIN" %in% names(F1)){
  F1$Agilent_DIN <- as.character(F1$Agilent_DIN)
  F1$Agilent_DIN <- as.numeric(F1$Agilent_DIN)
}
if("Agilent_RINe" %in% names(F1)){
  F1$Agilent_RINe <- as.character(F1$Agilent_RINe)
  F1$Agilent_RINe <- as.numeric(F1$Agilent_RINe)
}
F1$Asthma <- as.character(F1$Asthma)
F1$Gender <- as.character(F1$Gender)


# Subset data given choice of subset
sitesToInclude <- c("Baltimore", "Brazil", "Chicago", "Denver", "Nigeria", "Washington DC")[c(includeBaltimore,includeBrazil,includeChicago,includeDenver,includeNigeria,includeWashington_DC)]
F1 <- subset(F1, F1$site %in% sitesToInclude)

AgeToInclude <- c("Adult", "Child")[c(includeAdult, includeChild)]
F1 <- subset(F1, F1$Age_category %in% AgeToInclude)

# For Asthma: 1=control, 2=case
AsthmaToInclude <- c("1", "2")[c(includeNonAsthma, includeAsthma)]
F1 <- subset(F1, F1$Asthma %in% AsthmaToInclude)

# For Age: 1=male, 2=female
GenderToInclude <- c("1", "2")[c(includeMales, includeFemales)]
F1 <- subset(F1, F1$Gender %in% GenderToInclude)

# Select meaningful columns
dataFocus <- as_tibble(F1)
dataFocusSelect <-  dataFocus %>% select(starts_with("TOTAL_Nanodrop_ug"), 
                                         starts_with("TOTAL_Qubit_ug"), 
                                         starts_with("Nanodrop_260_280"), 
                                         starts_with("Agilent_DIN"), 
                                         starts_with("Agilent_RINe"),
                                         starts_with("PBMC_cellcounts"),
                                         starts_with("Slide_status")
                                         )

glimpse(dataFocusSelect)
summary(dataFocusSelect)

c1 <- "slategrey"
# c2 <- "wheat2"
c2 <- "slategray1"
my_density <- function(data,mapping,...){ggplot(data=data,mapping=mapping)+geom_density(...,lwd=1, colour = c1, fill = c2)}

cGramFunc <- function(inputData, colorColumn){
  if(is.null(colorColumn)){
    cGramPlot <- ggpairs(inputData, title="Correlation of QC Metrics", 
                         mapping=ggplot2::aes_string(alpha = 0.65),
                         lower = list(continuous = wrap("points", alpha = 0.65, colour = c1, size=0.8), 
                                      combo = wrap("dot", alpha = 0.65, colour = c1, size=0.8) ),
                         diag = list(continuous = my_density)
    )
  }else{
    cGramPlot <- ggpairs(inputData, title="Correlation of QC Metrics", 
                         mapping=ggplot2::aes_string(colour = colorColumn, alpha = 0.65),
                         lower = list(continuous = wrap("points", alpha = 0.65,    size=0.8), 
                                      combo = wrap("dot", alpha = 0.65,            size=0.8) )
    )
  }
  
  
  
  index_TOTAL_Nanodrop_ug <- which(names(inputData) == "TOTAL_Nanodrop_ug")
  index_TOTAL_Qubit_ug <- which(names(inputData) == "TOTAL_Qubit_ug")
  index_Nanodrop_260_280 <- which(names(inputData) == "Nanodrop_260_280")
  index_Agilent_DIN <- which(names(inputData) == "Agilent_DIN")
  index_Agilent_RINe <- which(names(inputData) == "Agilent_RINe")
  
  toKeep <- which(c("TOTAL_Nanodrop_ug", "TOTAL_Qubit_ug", "Nanodrop_260_280", "Agilent_DIN",
                    "Agilent_RINe")
                  %in% cGramPlot$xAxisLabels)
  
  indexVec <- c(index_TOTAL_Nanodrop_ug, index_TOTAL_Qubit_ug, index_Nanodrop_260_280, index_Agilent_DIN,
                index_Agilent_RINe)
  
  LTHRESHVec <- c(LTHRESH_TOTAL_Nanodrop_ug, LTHRESH_TOTAL_Qubit_ug, LTHRESH_Nanodrop_260_280, LTHRESH_Agilent_DIN,
                  LTHRESH_Agilent_RINe)[toKeep]
  
  UTHRESH_Vec <- c(UTHRESH_TOTAL_Nanodrop_ug, UTHRESH_TOTAL_Qubit_ug, UTHRESH_Nanodrop_260_280, UTHRESH_Agilent_DIN,
                  UTHRESH_Agilent_RINe)[toKeep]
  
  thresholdColors1 <- "lightgoldenrod3"
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

if("PBMC_cellcounts" %in% names(dataFocusSelect)){
  if(PBMC_PASS_ONLY == "No"){
    dataFocusSelectSub <- subset(dataFocusSelect, PBMC_cellcounts %in% c("PASS","FAIL"))
    if(Include_Boxplots == "No"){
      dataFocusSelectSub$PBMC_cellcounts <- NULL
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = NULL)
    }else{
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "PBMC_cellcounts")
    }
  }else if(PBMC_PASS_ONLY == "Yes"){
    dataFocusSelectSub <- subset(dataFocusSelect, PBMC_cellcounts %in% c("PASS"))
    if(Include_Boxplots == "No"){
      dataFocusSelectSub$PBMC_cellcounts <- NULL
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = NULL)
    }else{
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "PBMC_cellcounts")
    }
  }
}else if("Slide_status" %in% names(dataFocusSelect)){
  if(Slide_status_PASS_ONLY == "No"){
    dataFocusSelectSub <- subset(dataFocusSelect, Slide_status %in% c("PASS","FAIL"))
    if(Include_Boxplots == "No"){
      dataFocusSelectSub$Slide_status <- NULL
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = NULL)
    }else{
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "Slide_status")
    }
  }else if(Slide_status_PASS_ONLY == "Yes"){
    dataFocusSelectSub <- subset(dataFocusSelect, Slide_status %in% c("PASS"))
    if(Include_Boxplots == "No"){
      dataFocusSelectSub$Slide_status <- NULL
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = NULL)
    }else{
      cGram <- cGramFunc(inputData = dataFocusSelectSub, colorColumn = "Slide_status")
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
  inputData$FAIL_ANY <- as.numeric(inputData$FAIL_ANY)
  
  inputData$aesCol <- "UNCLASSIFIED"
  inputData$aesCol[inputData$FAIL_ANY == 0] <- "Passes all thresholds"
  inputData$aesCol[inputData$FAIL_ANY == 1 & inputData$var12Pass == "PASS"] <- "Passes here, fails 1 other threshold"
  inputData$aesCol[inputData$FAIL_ANY >= 2 & inputData$var12Pass == "PASS"] <- "Passes here, fails >1 other threshold"
  inputData$aesCol[inputData$FAIL_ANY == 1 & inputData$var12Pass == "FAIL"] <- "Fails only these thresholds"
  inputData$aesCol[inputData$FAIL_ANY >= 2 & inputData$var12Pass == "FAIL"] <- "Fails these and other thresholds"
  inputData$aesCol[as.numeric(inputData[["TOTAL_Nanodrop_ug"]]) == -5 | as.numeric(inputData[["TOTAL_Qubit_ug"]]) == -5] <- "Not yet processed (missing)"
  inputData$aesCol[as.numeric(inputData[["Nanodrop_260_280"]]) == 0] <- "Not yet processed (missing)"
  if("Agilent_DIN" %in% names(inputData)){
    inputData$aesCol[as.numeric(inputData[["Agilent_DIN"]]) == -1] <- "Not yet processed (missing)"
  }
  if("Agilent_RINe" %in% names(inputData)){
    inputData$aesCol[as.numeric(inputData[["Agilent_RINe"]]) == -1] <- "Not yet processed (missing)"
  }
  
  
  # LTHRESH_Nanodrop_260_280 = 1.7
  # UTHRESH_Nanodrop_260_280 = 2.2
  
  for(q in 1:nrow(inputData)){
    if("Agilent_DIN" %in% names(inputData)){
      if(as.numeric(inputData$Nanodrop_260_280[q]) < 1.4 | as.numeric(inputData$Nanodrop_260_280[q]) > 2.15){
        if(as.numeric(inputData$TOTAL_Nanodrop_ug[q]) >= 0.75){
          if(as.numeric(inputData$TOTAL_Qubit_ug[q]) >= 0.75){
            if(as.numeric(inputData$Agilent_DIN[q]) >= 7){
              inputData$aesCol[q] <- "Manual Pass"
            }
          }
        }
      }
    }
    if("Agilent_RINe" %in% names(inputData)){
      if(as.numeric(inputData$Nanodrop_260_280[q]) < 1.7 | as.numeric(inputData$Nanodrop_260_280[q]) > 2.2){
        if(as.numeric(inputData$TOTAL_Nanodrop_ug[q]) >= 0.6){
          if(as.numeric(inputData$TOTAL_Qubit_ug[q]) >= 0.6){
            if(as.numeric(inputData$Agilent_RINe[q]) >= 7){
              inputData$aesCol[q] <- "Manual Pass"
            }
          }
        }
      }
    }
  }
  
  inputData_old <- inputData
  
  inputData$aesCol <- as.character(inputData$aesCol)
  inputData$aesCol[inputData$aesCol == "Passes all thresholds"] <- paste0("Passes all thresholds; N = ",length(inputData$aesCol[inputData$aesCol == "Passes all thresholds"]))
  inputData$aesCol[inputData$aesCol == "Passes here, fails 1 other threshold"] <- paste0("Passes here, fails 1 other threshold; N = ",length(inputData$aesCol[inputData$aesCol == "Passes here, fails 1 other threshold"]))
  inputData$aesCol[inputData$aesCol == "Passes here, fails >1 other threshold"] <- paste0("Passes here, fails >1 other threshold; N = ",length(inputData$aesCol[inputData$aesCol == "Passes here, fails >1 other threshold"]))
  inputData$aesCol[inputData$aesCol == "Fails only these thresholds"] <- paste0("Fails only these thresholds; N = ",length(inputData$aesCol[inputData$aesCol == "Fails only these thresholds"]))
  inputData$aesCol[inputData$aesCol == "Fails these and other thresholds"] <- paste0("Fails these and other thresholds; N = ",length(inputData$aesCol[inputData$aesCol == "Fails these and other thresholds"]))
  inputData$aesCol[inputData$aesCol == "Not yet processed (missing)"] <- paste0("Not yet processed (missing); N = ",length(inputData$aesCol[inputData$aesCol == "Not yet processed (missing)"]))
  inputData$aesCol[inputData$aesCol == "Manual Pass"] <- paste0("Manual Pass; N = ",length(inputData$aesCol[inputData$aesCol == "Manual Pass"]))
  
  orderedFact <- c(paste0("Passes all thresholds; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Passes all thresholds"])),
                   paste0("Manual Pass; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Manual Pass"])),
                   paste0("Passes here, fails 1 other threshold; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Passes here, fails 1 other threshold"])),
                   paste0("Passes here, fails >1 other threshold; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Passes here, fails >1 other threshold"])),
                   paste0("Fails only these thresholds; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Fails only these thresholds"])),
                   paste0("Fails these and other thresholds; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Fails these and other thresholds"])),
                   paste0("Not yet processed (missing); N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Not yet processed (missing)"]))
  )
  # print(orderedFact)
  thresh1_lowerQ <- thresh1_lower
  thresh1_upperQ <- thresh1_upper
  if(thresh1_upperQ == Inf){thresh1_upperQ <- max(inputData[[var1]])}
  thresh2_lowerQ <- thresh2_lower
  thresh2_upperQ <- thresh2_upper
  if(thresh2_upperQ == Inf){thresh2_upperQ <- max(inputData[[var2]])}
  
  inputData$aesCol <- as.character(inputData$aesCol)
  inputData$aesCol <- factor(inputData$aesCol, levels = orderedFact, ordered = T)
  
  topLeft <- length(which(inputData[[var1]] < thresh1_lower & inputData[[var2]] > thresh2_upper))
  topLeftX <- thresh1_lowerQ - (max(inputData[[var1]]) - min(inputData[[var1]]))*.05
  topLeftY <- thresh2_upperQ + (max(inputData[[var2]]) + min(inputData[[var2]]))*.05
  
  topMiddle <- length(which(inputData[[var1]] >= thresh1_lower & inputData[[var1]] <= thresh1_upper & inputData[[var2]] > thresh2_upper))
  topMiddleX <- mean(c(thresh1_lowerQ,thresh1_upperQ))
  topMiddleY <- thresh2_upperQ + (max(inputData[[var2]]) + min(inputData[[var2]]))*.05
  
  topRight <- length(which(inputData[[var1]] > thresh1_upper & inputData[[var2]] > thresh2_upper))
  topRightX <- thresh1_upperQ + (max(inputData[[var1]]) - min(inputData[[var1]]))*.05
  topRightY <- thresh2_upperQ + (max(inputData[[var2]]) + min(inputData[[var2]]))*.05
  
  middleRight <- length(which(inputData[[var1]] > thresh1_upper & inputData[[var2]] >= thresh2_lower & inputData[[var2]] <= thresh2_upper))
  middleRightX <- thresh1_upperQ + (max(inputData[[var1]]) - min(inputData[[var1]]))*.05
  middleRightY <- mean(c(thresh2_lowerQ,thresh2_upperQ))
  
  bottomRight <- length(which(inputData[[var1]] > thresh1_upper & inputData[[var2]] < thresh2_lower))
  bottomRightX <- thresh1_upperQ + (max(inputData[[var1]]) - min(inputData[[var1]]))*.05
  bottomRightY <- thresh2_lowerQ - (max(inputData[[var2]]) + min(inputData[[var2]]))*.05
  
  bottomMiddle <- length(which(inputData[[var1]] >= thresh1_lower & inputData[[var1]] <= thresh1_upper & inputData[[var2]] < thresh2_lower))
  bottomMiddleX <- mean(c(thresh1_lowerQ,thresh1_upperQ))
  bottomMiddleY <- thresh2_lowerQ - (max(inputData[[var2]]) + min(inputData[[var2]]))*.05
  
  bottomLeft <- length(which(inputData[[var1]] < thresh1_lower & inputData[[var2]] < thresh2_lower))
  bottomLeftX <- thresh1_lowerQ - (max(inputData[[var1]]) - min(inputData[[var1]]))*.05
  bottomLeftY <- thresh2_lowerQ - (max(inputData[[var2]]) + min(inputData[[var2]]))*.05
  
  middleLeft <- length(which(inputData[[var1]] < thresh1_lower & inputData[[var2]] <= thresh2_upper & inputData[[var2]] >= thresh2_lower))
  middleLeftX <- thresh1_lowerQ - (max(inputData[[var1]]) - min(inputData[[var1]]))*.05
  middleLeftY <- mean(c(thresh2_lowerQ,thresh2_upperQ))
  
  middleMiddle <- length(which(inputData[[var1]] >= thresh1_lower & inputData[[var1]] <= thresh1_upper & inputData[[var2]] >= thresh2_lower & inputData[[var2]] <= thresh2_upper))
  middleMiddleX <- mean(c(thresh1_lowerQ,thresh1_upperQ))
  middleMiddleY <- mean(c(thresh2_lowerQ,thresh2_upperQ))
  
  # a data frame with all the annotation info
  annotation <- data.frame(
    x = c(topLeftX, topMiddleX, topRightX, middleRightX, bottomRightX, bottomMiddleX, bottomLeftX, middleLeftX, middleMiddleX),
    y = c(topLeftY, topMiddleY, topRightY, middleRightY, bottomRightY, bottomMiddleY, bottomLeftY, middleLeftY, middleMiddleY),
    label = c(topLeft, topMiddle, topRight, middleRight, bottomRight, bottomMiddle, bottomLeft, middleLeft, middleMiddle)
  )
  annotation[annotation==0] <- NA
  
  thresholdColors2 <- "lightgoldenrod3"
  
  myPlot <- ggplot()+
    geom_point(data = inputData, aes_string(x = var1, y = var2, color = "aesCol"), size = 1, alpha = 0.8)+
    # geom_point(data = inputData, aes_string(x = var1, y = var2, color = "var12Pass", shape = "FAIL_ANY"), size = 1)+
    geom_vline(xintercept = thresh1_lower, color = thresholdColors2)+
    geom_vline(xintercept = thresh1_upper, color = thresholdColors2)+
    geom_hline(yintercept = thresh2_lower, color = thresholdColors2)+
    geom_hline(yintercept = thresh2_upper, color = thresholdColors2)+
    # scale_color_manual(NULL, breaks = c("FAIL", "PASS"), labels = c("FAIL", "PASS"), values = c("firebrick3", "forestgreen"))+
    scale_color_manual(paste0("Total N = ",length(inputData$aesCol)), 
                       # labels = c(paste0("FAILS THESE THRESHOLDS; N = ",length(inputData$aesCol[inputData$aesCol == "FAILS THESE THRESHOLDS"])), 
                       #            paste0("FAILS OTHER THRESHOLD; N = ",length(inputData$aesCol[inputData$aesCol == "FAILS OTHER THRESHOLD"])),  
                       #            paste0("PASSES ALL THESHOLDS; N = ",length(inputData$aesCol[inputData$aesCol == "PASSES ALL THESHOLDS"]))
                       #            ), 
                       breaks = c(paste0("Passes all thresholds; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Passes all thresholds"])),
                                  paste0("Passes here, fails 1 other threshold; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Passes here, fails 1 other threshold"])),
                                  paste0("Passes here, fails >1 other threshold; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Passes here, fails >1 other threshold"])),
                                  paste0("Fails only these thresholds; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Fails only these thresholds"])),
                                  paste0("Fails these and other thresholds; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Fails these and other thresholds"])),
                                  paste0("Not yet processed (missing); N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Not yet processed (missing)"])),
                                  paste0("Manual Pass; N = ",length(inputData_old$aesCol[inputData_old$aesCol == "Manual Pass"]))
                       ),
                       values = c("green", "cyan2", "purple", "orange", "darkred", "honeydew3", "turquoise4"))+
                        # values = c("orange", "darkred", "green", "purple", "deepskyblue"))+
                       # values = c("#44AA99","#332288", "#D55E00", ))+ #blue red green
    # scale_shape_manual(NULL, breaks = c("FAIL ANY", "PASS ALL"), labels = c("FAIL ANY", "PASS ALL"), values = c(4, 16))
  # scale_fill_manual(labels = c("PASS ALL", "FAIL ANY"), values = c("black", "pink"))+
  theme_bw(base_size = 8)+
    geom_text(data=annotation, aes( x=x, y=y, label=label),
               color="gray22", 
               size=3, angle=0, fontface="italic" )
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
        # idat <- dataFocusSelectSubContinue %>% select(VarVec[i], VarVec[j], FAIL_ANY)
        idat <- dataFocusSelectSubContinue
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
# pdf("testing.pdf", width = 20, height = 15)
if("Agilent_DIN" %in% names(dataFocusSelect)){
  # DNA THRESHOLDS
  # LTHRESH_Nanodrop_260_280 <- 1.4
  # UTHRESH_Nanodrop_260_280 <- 2.1
  # LTHRESH_Agilent_DIN <- 6
  # UTHRESH_Agilent_DIN <- Inf
  # LTHRESH_TOTAL_Nanodrop_ug <- 0.75
  # UTHRESH_TOTAL_Nanodrop_ug <- Inf
  # LTHRESH_TOTAL_Qubit_ug <- 0.75
  # UTHRESH_TOTAL_Qubit_ug <- Inf
  
  # Subset to those that pass
  if("PBMC_cellcounts" %in% names(dataFocusSelect)){
    if(PBMC_PASS_ONLY == "No"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_cellcounts %in% c("PASS","FAIL"))
    }else if(PBMC_PASS_ONLY == "Yes"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_cellcounts %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- 0
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Nanodrop_260_280"]] < LTHRESH_Nanodrop_260_280 | dataFocusSelectSubContinue[["Nanodrop_260_280"]] > UTHRESH_Nanodrop_260_280 ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Nanodrop_260_280"]] < LTHRESH_Nanodrop_260_280 | dataFocusSelectSubContinue[["Nanodrop_260_280"]] > UTHRESH_Nanodrop_260_280 ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Agilent_DIN"]] < LTHRESH_Agilent_DIN | dataFocusSelectSubContinue[["Agilent_DIN"]] > UTHRESH_Agilent_DIN ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Agilent_DIN"]] < LTHRESH_Agilent_DIN | dataFocusSelectSubContinue[["Agilent_DIN"]] > UTHRESH_Agilent_DIN ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] < LTHRESH_TOTAL_Nanodrop_ug | dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] > UTHRESH_TOTAL_Nanodrop_ug ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] < LTHRESH_TOTAL_Nanodrop_ug | dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] > UTHRESH_TOTAL_Nanodrop_ug ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] < LTHRESH_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] > UTHRESH_TOTAL_Qubit_ug ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] < LTHRESH_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] > UTHRESH_TOTAL_Qubit_ug ] + 1
    
    # dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
    # run grid plotter
    myGridPlots <- gridMaker(inputData = dataFocusSelectSubContinue, VarVec = c("Nanodrop_260_280", "Agilent_DIN", "TOTAL_Nanodrop_ug", "TOTAL_Qubit_ug") )
    # plotDF <- data.frame(matrix(unlist(1:4^2), nrow=4, byrow=T),stringsAsFactors=FALSE)
    # plotDFmat <- as.matrix(plotDF)
    # gridPlot <- grid.arrange(grobs = myGridPlots, layout_matrix = plotDFmat)
    # print(gridPlot)
    
    summaryTable <- data.frame(matrix(ncol = 0, nrow = 3))
    summaryTable$`_____` <- c("N", "PASS ALL", "FAIL ANY")
    summaryTable$Total <- c(nrow(dataFocusSelectSubContinue),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Pass_PBMC_cellcounts <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "PASS")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Fail_PBMC_cellcounts <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "FAIL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    # summaryTable
    
  }else if("Slide_status" %in% names(dataFocusSelect)){
    if(Slide_status_PASS_ONLY == "No"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Slide_status %in% c("PASS","FAIL"))
    }else if(Slide_status_PASS_ONLY == "Yes"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Slide_status %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- 0
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Nanodrop_260_280"]] < LTHRESH_Nanodrop_260_280 | dataFocusSelectSubContinue[["Nanodrop_260_280"]] > UTHRESH_Nanodrop_260_280 ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Nanodrop_260_280"]] < LTHRESH_Nanodrop_260_280 | dataFocusSelectSubContinue[["Nanodrop_260_280"]] > UTHRESH_Nanodrop_260_280 ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Agilent_DIN"]] < LTHRESH_Agilent_DIN | dataFocusSelectSubContinue[["Agilent_DIN"]] > UTHRESH_Agilent_DIN ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Agilent_DIN"]] < LTHRESH_Agilent_DIN | dataFocusSelectSubContinue[["Agilent_DIN"]] > UTHRESH_Agilent_DIN ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] < LTHRESH_TOTAL_Nanodrop_ug | dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] > UTHRESH_TOTAL_Nanodrop_ug ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] < LTHRESH_TOTAL_Nanodrop_ug | dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] > UTHRESH_TOTAL_Nanodrop_ug ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] < LTHRESH_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] > UTHRESH_TOTAL_Qubit_ug ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] < LTHRESH_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] > UTHRESH_TOTAL_Qubit_ug ] + 1
    
    # dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
    # run grid plotter
    myGridPlots <- gridMaker(inputData = dataFocusSelectSubContinue, VarVec = c("Nanodrop_260_280", "Agilent_DIN", "TOTAL_Nanodrop_ug", "TOTAL_Qubit_ug") )
    # plotDF <- data.frame(matrix(unlist(1:4^2), nrow=4, byrow=T),stringsAsFactors=FALSE)
    # plotDFmat <- as.matrix(plotDF)
    # gridPlot <- grid.arrange(grobs = myGridPlots, layout_matrix = plotDFmat)
    # print(gridPlot)
    
    summaryTable <- data.frame(matrix(ncol = 0, nrow = 3))
    summaryTable$`_____` <- c("N", "PASS ALL", "FAIL ANY")
    summaryTable$Total <- c(nrow(dataFocusSelectSubContinue),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Pass_Slide_status <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "PASS")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Fail_Slide_status <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "FAIL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    # summaryTable
  }
}else if("Agilent_RINe" %in% names(dataFocusSelect)){
  # RNA THRESHOLDS
  # LTHRESH_Nanodrop_260_280 <- 1.7
  # UTHRESH_Nanodrop_260_280 <- 2.2
  # LTHRESH_Agilent_RINe <- 6
  # UTHRESH_Agilent_RINe <- Inf
  # LTHRESH_TOTAL_Nanodrop_ug <- 0.75
  # UTHRESH_TOTAL_Nanodrop_ug <- Inf
  # LTHRESH_TOTAL_Qubit_ug <- 0.75
  # UTHRESH_TOTAL_Qubit_ug <- Inf
  
  # Subset to those that pass
  if("PBMC_cellcounts" %in% names(dataFocusSelect)){
    if(PBMC_PASS_ONLY == "No"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_cellcounts %in% c("PASS","FAIL"))
    }else if(PBMC_PASS_ONLY == "Yes"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, PBMC_cellcounts %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- 0
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Nanodrop_260_280"]] < LTHRESH_Nanodrop_260_280 | dataFocusSelectSubContinue[["Nanodrop_260_280"]] > UTHRESH_Nanodrop_260_280 ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Nanodrop_260_280"]] < LTHRESH_Nanodrop_260_280 | dataFocusSelectSubContinue[["Nanodrop_260_280"]] > UTHRESH_Nanodrop_260_280 ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Agilent_RINe"]] < LTHRESH_Agilent_RINe | dataFocusSelectSubContinue[["Agilent_RINe"]] > UTHRESH_Agilent_RINe ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Agilent_RINe"]] < LTHRESH_Agilent_RINe | dataFocusSelectSubContinue[["Agilent_RINe"]] > UTHRESH_Agilent_RINe ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] < LTHRESH_TOTAL_Nanodrop_ug | dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] > UTHRESH_TOTAL_Nanodrop_ug ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] < LTHRESH_TOTAL_Nanodrop_ug | dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] > UTHRESH_TOTAL_Nanodrop_ug ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] < LTHRESH_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] > UTHRESH_TOTAL_Qubit_ug ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] < LTHRESH_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] > UTHRESH_TOTAL_Qubit_ug ] + 1
    
    # dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
    # run grid plotter
    myGridPlots <- gridMaker(inputData = dataFocusSelectSubContinue, VarVec = c("Nanodrop_260_280", "Agilent_RINe", "TOTAL_Nanodrop_ug", "TOTAL_Qubit_ug") )
    # plotDF <- data.frame(matrix(unlist(1:4^2), nrow=4, byrow=T),stringsAsFactors=FALSE)
    # plotDFmat <- as.matrix(plotDF)
    # gridPlot <- grid.arrange(grobs = myGridPlots, layout_matrix = plotDFmat)
    # print(gridPlot)
    
    summaryTable <- data.frame(matrix(ncol = 0, nrow = 3))
    summaryTable$`_____` <- c("N", "PASS ALL", "FAIL ANY")
    summaryTable$Total <- c(nrow(dataFocusSelectSubContinue),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Pass_PBMC_cellcounts <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "PASS")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Fail_PBMC_cellcounts <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "FAIL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                       nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$PBMC_cellcounts == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    # summaryTable
    
  }else if("Slide_status" %in% names(dataFocusSelect)){
    if(Slide_status_PASS_ONLY == "No"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Slide_status %in% c("PASS","FAIL"))
    }else if(Slide_status_PASS_ONLY == "Yes"){
      dataFocusSelectSubContinue <- subset(dataFocusSelect, Slide_status %in% c("PASS"))
    }
    
    dataFocusSelectSubContinue$FAIL_ANY <- 0
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Nanodrop_260_280"]] < LTHRESH_Nanodrop_260_280 | dataFocusSelectSubContinue[["Nanodrop_260_280"]] > UTHRESH_Nanodrop_260_280 ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Nanodrop_260_280"]] < LTHRESH_Nanodrop_260_280 | dataFocusSelectSubContinue[["Nanodrop_260_280"]] > UTHRESH_Nanodrop_260_280 ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Agilent_RINe"]] < LTHRESH_Agilent_RINe | dataFocusSelectSubContinue[["Agilent_RINe"]] > UTHRESH_Agilent_RINe ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["Agilent_RINe"]] < LTHRESH_Agilent_RINe | dataFocusSelectSubContinue[["Agilent_RINe"]] > UTHRESH_Agilent_RINe ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] < LTHRESH_TOTAL_Nanodrop_ug | dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] > UTHRESH_TOTAL_Nanodrop_ug ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] < LTHRESH_TOTAL_Nanodrop_ug | dataFocusSelectSubContinue[["TOTAL_Nanodrop_ug"]] > UTHRESH_TOTAL_Nanodrop_ug ] + 1
    dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] < LTHRESH_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] > UTHRESH_TOTAL_Qubit_ug ] <- 
      dataFocusSelectSubContinue$FAIL_ANY[dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] < LTHRESH_TOTAL_Qubit_ug | dataFocusSelectSubContinue[["TOTAL_Qubit_ug"]] > UTHRESH_TOTAL_Qubit_ug ] + 1
    
    # dataFocusSelectSubContinue$FAIL_ANY <- as.factor(dataFocusSelectSubContinue$FAIL_ANY)
    
    # run grid plotter
    myGridPlots <- gridMaker(inputData = dataFocusSelectSubContinue, VarVec = c("Nanodrop_260_280", "Agilent_RINe", "TOTAL_Nanodrop_ug", "TOTAL_Qubit_ug") )
    # plotDF <- data.frame(matrix(unlist(1:4^2), nrow=4, byrow=T),stringsAsFactors=FALSE)
    # plotDFmat <- as.matrix(plotDF)
    # gridPlot <- grid.arrange(grobs = myGridPlots, layout_matrix = plotDFmat)
    # print(gridPlot)
    
    summaryTable <- data.frame(matrix(ncol = 0, nrow = 3))
    summaryTable$`_____` <- c("N", "PASS ALL", "FAIL ANY")
    summaryTable$Total <- c(nrow(dataFocusSelectSubContinue),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                            nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Pass_Slide_status <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "PASS")),
                                      nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                      nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "PASS" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    summaryTable$Fail_Slide_status <- c(nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "FAIL")),
                                      nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "PASS ALL")),
                                      nrow(subset(dataFocusSelectSubContinue, dataFocusSelectSubContinue$Slide_status == "FAIL" & dataFocusSelectSubContinue$FAIL_ANY == "FAIL ANY"))
    )
    # summaryTable
  }
}
print(summaryTable)
# dev.off()

pre_summaryTable <- dataFocus
pre_summaryTable_old <- dataFocus
pre_summaryTable$Nanodrop_260_280[pre_summaryTable[["Nanodrop_260_280"]] >= LTHRESH_Nanodrop_260_280 & pre_summaryTable[["Nanodrop_260_280"]] <= UTHRESH_Nanodrop_260_280 ] <- "PASS"
pre_summaryTable$Nanodrop_260_280[pre_summaryTable$Nanodrop_260_280 != "PASS"] <- "FAIL"
if("Agilent_RINe" %in% names(pre_summaryTable)){
  pre_summaryTable$Agilent_RINe[pre_summaryTable[["Agilent_RINe"]] >= LTHRESH_Agilent_RINe & pre_summaryTable[["Agilent_RINe"]] <= UTHRESH_Agilent_RINe ] <- "PASS"
  pre_summaryTable$Agilent_RINe[pre_summaryTable$Agilent_RINe != "PASS"] <- "FAIL"
}
if("Agilent_DIN" %in% names(pre_summaryTable)){
  pre_summaryTable$Agilent_DIN[pre_summaryTable[["Agilent_DIN"]] >= LTHRESH_Agilent_DIN & pre_summaryTable[["Agilent_DIN"]] <= UTHRESH_Agilent_DIN ] <- "PASS"
  pre_summaryTable$Agilent_DIN[pre_summaryTable$Agilent_DIN != "PASS"] <- "FAIL"
}
pre_summaryTable$TOTAL_Nanodrop_ug[pre_summaryTable[["TOTAL_Nanodrop_ug"]] >= LTHRESH_TOTAL_Nanodrop_ug & pre_summaryTable[["TOTAL_Nanodrop_ug"]] <= UTHRESH_TOTAL_Nanodrop_ug ] <- "PASS"
pre_summaryTable$TOTAL_Nanodrop_ug[pre_summaryTable$TOTAL_Nanodrop_ug != "PASS"] <- "FAIL"
pre_summaryTable$TOTAL_Qubit_ug[pre_summaryTable[["TOTAL_Qubit_ug"]] >= LTHRESH_TOTAL_Qubit_ug & pre_summaryTable[["TOTAL_Qubit_ug"]] <= UTHRESH_TOTAL_Qubit_ug ] <- "PASS"
pre_summaryTable$TOTAL_Qubit_ug[pre_summaryTable$TOTAL_Qubit_ug != "PASS"] <- "FAIL"
# For Age: 1=male, 2=female
pre_summaryTable$Gender[pre_summaryTable$Gender == 1] <- "Male"
pre_summaryTable$Gender[pre_summaryTable$Gender == 2] <- "Female"
# For Asthma: 1=control, 2=case
pre_summaryTable$Asthma[pre_summaryTable$Asthma == 1] <- "Control"
pre_summaryTable$Asthma[pre_summaryTable$Asthma == 2] <- "Case"
# pass all filters?
pre_summaryTable$pass_all <- "FAIL ANY"
if("Agilent_RINe" %in% names(pre_summaryTable)){
  pre_summaryTable$pass_all[pre_summaryTable$Nanodrop_260_280 == "PASS" & pre_summaryTable$Agilent_RINe == "PASS" & pre_summaryTable$TOTAL_Nanodrop_ug == "PASS" & pre_summaryTable$TOTAL_Qubit_ug == "PASS"] <- "PASS ALL"
}
if("Agilent_DIN" %in% names(pre_summaryTable)){
  pre_summaryTable$pass_all[pre_summaryTable$Nanodrop_260_280 == "PASS" & pre_summaryTable$Agilent_DIN == "PASS" & pre_summaryTable$TOTAL_Nanodrop_ug == "PASS" & pre_summaryTable$TOTAL_Qubit_ug == "PASS"] <- "PASS ALL"
}
for(q in 1:nrow(pre_summaryTable_old)){
  if("Agilent_DIN" %in% names(pre_summaryTable_old)){
    if(as.numeric(pre_summaryTable_old$Nanodrop_260_280[q]) < 1.4 | as.numeric(pre_summaryTable_old$Nanodrop_260_280[q]) > 2.15){
      if(as.numeric(pre_summaryTable_old$TOTAL_Nanodrop_ug[q]) >= 0.75){
        if(as.numeric(pre_summaryTable_old$TOTAL_Qubit_ug[q]) >= 0.75){
          if(as.numeric(pre_summaryTable_old$Agilent_DIN[q]) >= 7){
            pre_summaryTable$pass_all[q] <- "MANUAL PASS"
          }
        }
      }
    }
  }
  if("Agilent_RINe" %in% names(pre_summaryTable_old)){
    if(as.numeric(pre_summaryTable_old$Nanodrop_260_280[q]) < 1.7 | as.numeric(pre_summaryTable_old$Nanodrop_260_280[q]) > 2.2){
      if(as.numeric(pre_summaryTable_old$TOTAL_Nanodrop_ug[q]) >= 0.6){
        if(as.numeric(pre_summaryTable_old$TOTAL_Qubit_ug[q]) >= 0.6){
          if(as.numeric(pre_summaryTable_old$Agilent_RINe[q]) >= 7){
            pre_summaryTable$pass_all[q] <- "MANUAL PASS"
          }
        }
      }
    }
  }
}

# pre_summaryTable$subject <- NULL
pre_summaryTable$site <- as.factor(pre_summaryTable$site)
pre_summaryTable$type <- NULL
pre_summaryTable$site <- as.factor(pre_summaryTable$site)
pre_summaryTable$Age_category <- as.factor(pre_summaryTable$Age_category)
pre_summaryTable$Asthma <- as.factor(pre_summaryTable$Asthma)
pre_summaryTable$Gender <- as.factor(pre_summaryTable$Gender)
pre_summaryTable$TOTAL_Nanodrop_ug <- as.factor(pre_summaryTable$TOTAL_Nanodrop_ug)
pre_summaryTable$TOTAL_Qubit_ug <- as.factor(pre_summaryTable$TOTAL_Qubit_ug)
pre_summaryTable$Nanodrop_260_280 <- as.factor(pre_summaryTable$Nanodrop_260_280)
if("Agilent_RINe" %in% names(pre_summaryTable)){
  pre_summaryTable$Agilent_RINe <- as.factor(pre_summaryTable$Agilent_RINe)
}
if("Agilent_DIN" %in% names(pre_summaryTable)){
  pre_summaryTable$Agilent_DIN <- as.factor(pre_summaryTable$Agilent_DIN)
}

# Subset to only the primary passing samples...
if("PBMC_cellcounts" %in% names(pre_summaryTable)){
  if(PBMC_PASS_ONLY == "Yes"){
    pre_summaryTable <- subset(pre_summaryTable, PBMC_cellcounts %in% c("PASS"))
  }
  pre_summaryTable$PBMC_cellcounts <- NULL
}
if("Slide_status" %in% names(pre_summaryTable)){
  if(PBMC_PASS_ONLY == "Yes"){
    pre_summaryTable <- subset(pre_summaryTable, Slide_status %in% c("PASS"))
  }
  pre_summaryTable$Slide_status <- NULL
}
pre_summaryTable$pass_all <- as.factor(pre_summaryTable$pass_all)

# summary(pre_summaryTable)
if("Agilent_RINe" %in% names(pre_summaryTable)){
  pre_summaryTable_melt <- melt(pre_summaryTable, id=c("Nanodrop_260_280","Agilent_RINe","TOTAL_Nanodrop_ug","TOTAL_Qubit_ug","pass_all"))
}
if("Agilent_DIN" %in% names(pre_summaryTable)){
  pre_summaryTable_melt <- melt(pre_summaryTable, id=c("Nanodrop_260_280","Agilent_DIN","TOTAL_Nanodrop_ug","TOTAL_Qubit_ug","pass_all"))
}
pre_summaryTable_melt_total <- pre_summaryTable
pre_summaryTable_melt_total$site <- NULL
pre_summaryTable_melt_total$Age_category <- NULL
pre_summaryTable_melt_total$Asthma <- NULL
pre_summaryTable_melt_total$Gender <- NULL

# pre_summaryTable_melt$value <- factor(pre_summaryTable_melt$value, levels=c("Baltimore","Brazil","Chicago","Denver","Nigeria","Washington DC","Adult","Child","Case","Control","Male","Female"), ordered=T)
AsthmaToIncludeTable1 <- c("Control", "Case")[c(includeNonAsthma, includeAsthma)]
GenderToIncludeTable1 <- c("Male", "Female")[c(includeMales, includeFemales)]
pre_summaryTable_melt$value <- factor(pre_summaryTable_melt$value, levels=c(sitesToInclude,AgeToInclude,AsthmaToIncludeTable1,GenderToIncludeTable1), ordered=T)

if("Agilent_DIN" %in% names(pre_summaryTable)){
  labels <- list(
    variables=list(Nanodrop_260_280="Nanodrop_260_280",
                   Agilent_DIN="Agilent_DIN",
                   TOTAL_Nanodrop_ug="TOTAL_Nanodrop_ug",
                   TOTAL_Qubit_ug="TOTAL_Qubit_ug",
                   pass_all="All Filters"
    ),
    groups=list("","Recruitment Site", "Age Group", "Asthma Status", "Gender"))
  strata <- c(list(Overall=pre_summaryTable_melt_total), split(pre_summaryTable_melt, pre_summaryTable_melt$value))
  rndr <- function(x, name, ...) {
    if (!is.numeric(x)) return(render.categorical.default(x))
    what <- switch(name,
                   TOTAL_Nanodrop_ug = "FREQ, PCT",
                   Nanodrop_260_280 = "FREQ, PCT",
                   Agilent_DIN = "FREQ, PCT",
                   PBMC_cellcounts = "FREQ, PCT",
                   pass_all = "FREQ, PCT"
    )
    parse.abbrev.render.code(c("", what))(x)
  }
  # summaryTable <- table1(strata, labels, groupspan=c(1,6,2,2,2), render=rndr)
  summaryTable <- table1(strata, labels, groupspan=c(1,length(sitesToInclude),length(AgeToInclude),length(AsthmaToInclude),length(GenderToInclude)), render=rndr)
}else if("Agilent_RINe" %in% names(pre_summaryTable)){
  labels <- list(
    variables=list(Nanodrop_260_280="Nanodrop_260_280",
                   Agilent_RINe="Agilent_RINe",
                   TOTAL_Nanodrop_ug="TOTAL_Nanodrop_ug",
                   TOTAL_Qubit_ug="TOTAL_Qubit_ug",
                   pass_all="Pass All"
    ),
    groups=list("","Recruitment Site", "Age Group", "Asthma Status", "Gender"))
  strata <- c(list(Overall=pre_summaryTable_melt_total), split(pre_summaryTable_melt, pre_summaryTable_melt$value))
  rndr <- function(x, name, ...) {
    if (!is.numeric(x)) return(render.categorical.default(x))
    what <- switch(name,
                   TOTAL_Nanodrop_ug = "FREQ, PCT",
                   Nanodrop_260_280 = "FREQ, PCT",
                   Agilent_RINe = "FREQ, PCT",
                   PBMC_cellcounts = "FREQ, PCT",
                   pass_all = "FREQ, PCT"
    )
    parse.abbrev.render.code(c("", what))(x)
  }
  # summaryTable <- table1(strata, labels, groupspan=c(1,6,2,2,2), render=rndr)
  summaryTable <- table1(strata, labels, groupspan=c(1,length(sitesToInclude),length(AgeToInclude),length(AsthmaToInclude),length(GenderToInclude)), render=rndr)
}

# rownames(summaryTable) <- summaryTable$`_____`
# summaryTable$`_____` <- NULL
# print(summaryTable)

# save(summaryTable, cGram, myGridPlots, gridPlot, file="savePlots.RData")

preForCorey <- dataFocus
idsToKeep <- subset(pre_summaryTable, pass_all == "MANUAL PASS" | pass_all == "PASS ALL")
forCorey <- subset(preForCorey, subject %in% idsToKeep$subject)
# For Age: 1=male, 2=female
forCorey$Gender[forCorey$Gender == 1] <- "Male"
forCorey$Gender[forCorey$Gender == 2] <- "Female"
# For Asthma: 1=control, 2=case
forCorey$Asthma[forCorey$Asthma == 1] <- "Control"
forCorey$Asthma[forCorey$Asthma == 2] <- "Case"

# F1_others <- F1_all %>% select(subject, Container_ID, Volume, Conc_Nanodrop, Conc_Qubit)

forCorey <- left_join(forCorey, F1_others, by = "subject")

save(summaryTable, cGram, myGridPlots, file="savePlots.RData")
save(forCorey, file="savePassed.RData")

print("SCRIPT HAS FINISHED")




