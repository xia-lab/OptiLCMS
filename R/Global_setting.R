
# Setting .on.public.web
.on.public.web <- FALSE;

# OTHER SETTINGS

# Used to defined the parallel namespace for peak picking
peak_function_list <- list("PerformPeakPicking",
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
                           ".aggregateMethod2int"
                           
                           
)

# Used to defined the parallel namespace for parameters optimization
optimize_function_list <- c(list("PeakPicking_prep",
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
                            peak_function_list)

####Scheduler

SetRawFileNames <- function(filenms){
  print(filenms)
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
    print("wrongfiletype")
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
      print(replicateNum)
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


