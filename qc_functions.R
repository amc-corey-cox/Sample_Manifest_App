gg_density <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(..., lwd=1, colour = c1, fill = c2)
}

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
    }
  }
  
  return(cGramPlot)
}