
# Setting .on.public.web
.on.public.web <- FALSE;

# OTHER SETTINGS
# NA

# Used to defined the parallel namespace for peak picking
.peak_function_list <- list("PerformPeakPicking",
                           "PeakPicking_centWave_slave",
                           "PerformPeakGrouping",
                           "Densitygrouping_slave",
                           "PerformRTcorrection",
                           "PerformPeakAlignment",
                           "RT.Adjust_peakGroup",
                           "adjustRtime_obiwarp",
                           "adjustRtime_peakGroup",
                           "mSet.obiwarp",
                           "PerformPeakFiling",
                           "mSet2xcmsSet",
                           "updateRawSpectraParam",
                           "continuousPtsAboveThreshold",
                           "getLocalNoiseEstimate",
                           "continuousPtsAboveThresholdIdx",
                           "MSW.cwt",
                           "MSW.extendNBase",
                           "MSW.extendLength",
                           "MSW.getLocalMaximumCWT",
                           "MSW.localMaximum",
                           "MSW.getRidge",
                           "descendMin",
                           "descendMinTol",
                           "joinOverlappingPeaks",
                           "trimm",
                           "findEqualGreaterM",
                           "na.flatfill",
                           "SSgauss",
                           "rectUnique",
                           ".narrow_rt_boundaries",
                           ".rawMat",
                           ".getPeakGroupsRtMatrix",
                           ".peakIndex",
                           ".applyRtAdjToChromPeaks",
                           ".applyRtAdjustment",
                           ".getChromPeakData",
                           ".feature_values",
                           ".groupval",
                           ".createProfileMatrix",
                           "profMat",
                           "binYonX",
                           "imputeLinInterpol",
                           "colMax",
                           "filtfft",
                           "descendZero",
                           "which.colMax",
                           "breaks_on_nBins",
                           ".aggregateMethod2int",
                           "findmzROI",
                           "getEIC",
                           "getMZ",
                           "weighted.mean"
                           
                           
)

# Used to defined the parallel namespace for parameters optimization
.optimize_function_list <- c(list("PeakPicking_prep",
                                 "Statistic_doe",
                                 "SlaveCluster_doe",
                                 "calculateSet_doe",
                                 "SlaveCluster_doe",
                                 "calculateSet_doe",
                                 "calculatePPKs",
                                 "calculateGPRT",
                                 "calcPPS2",
                                 "calcCV",
                                 "resultIncreased_doe",
                                 "calcRCS_GSValues",
                                 "calcGaussianS",
                                 "peaks_IPO",
                                 "getClusterType",
                                 "calcMaximumCarbon",
                                 "getMaximumLevels",
                                 "getMaxSettings",
                                 "expand.grid.subset",
                                 "typeCastParams",
                                 "getCcdParameter",
                                 "combineParams",
                                 "getRGTVValues",
                                 "findIsotopes.IPO",
                                 "creatPeakTable",
                                 "createModel",
                                 "decode",
                                 "decodeAll",
                                 "encode",
                                 "attachList",
                                 "checkParams"),
                            .peak_function_list)

####Scheduler

SetRawFileNames <- function(filenms){
  cat(filenms)
  rawfilenms.vec <<- filenms
  return(1);
}

####Raw Spectra Upload
ReadRawMeta<-function(fileName){
  if(grepl(".txt", fileName, fixed=T)){
    tbl=read.table(fileName,header=TRUE, stringsAsFactors = F);
  }else if(grepl(".csv", fileName, fixed=T)){
    tbl = read.csv(fileName,header=TRUE, stringsAsFactors = F);
  }else{
    cat("wrongfiletype\n")
  }
  
  rawFileNms<-as.vector(tbl[,1])
  rawClassNms<-as.vector(tbl[,2])
  rawFileNms <- sapply(strsplit(rawFileNms, "\\."), function(x) paste0(head(x,-1), collapse="."));
  clsTable = table(rawClassNms)
  #check replicate number
  clsTypes = names(table(rawClassNms))
  for(name in clsTypes){
    if(toupper(name) !="QC"){
      replicateNum = clsTable[[name]]
      cat(replicateNum,"\n")
    }
  }
  
  rawFileNms<<-rawFileNms
  rawClassNms<<-rawClassNms
  return(1);
}

GetRawFileNms <- function(){
  return(rawFileNms)
}

GetRawClassNms <- function(){
  return(rawClassNms)
}


#' Title
#'
#' @param mes 
#' @param ecol 
#' @param progress 
#'
#' @return
#'
#' @examples
MessageOutput <- function(mes, ecol, progress) {
  if (!is.null(mes)) {
    if (.on.public.web) {
      # write down message
      write.table(
        mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = ecol
      )
    } else {
      # print message
      if(ecol == "\n"){
        message(mes)
      } else {
        message(mes,appendLF = FALSE)
      }
      
    }
  }
  
  # write down progress
  if (.on.public.web & !is.null(progress)) {
    progress <- as.numeric(progress)
    
    write.table(
      progress,
      file = paste0(fullUserPath, "log_progress.txt"),
      row.names = F,
      col.names = F
    )
  }
  
}

fast.write.csv <- function(dat, file, row.names=TRUE){
  tryCatch(
    {
      if(is.data.frame(dat)){
        # there is a rare bug in data.table (R 3.6) which kill the R process in some cases 
        data.table::fwrite(dat, file, row.names=row.names);
      }else{
        write.csv(dat, file, row.names=row.names);  
      }
    }, error=function(e){
      print(e);
      write.csv(dat, file, row.names=row.names);   
    }, warning=function(w){
      print(w);
      write.csv(dat, file, row.names=row.names); 
    });
}



