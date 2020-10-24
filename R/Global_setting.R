# Setting .on.public.web (if on.public.web, the version number will be even, local package is odd)
.on.public.web <- FALSE;

# Setting the global Variable to avoid notes in R CMD Check
utils::globalVariables(c(".SwapEnv"))

# OTHER SETTINGS
#' @references Gatto L, Gibb S, Rainer J (2020). “MSnbase, efficient and elegant R-based processing and visualisation of raw mass spectrometry data.” bioRxiv.
.MSnExpReqFvarLabels <- c("fileIdx", "spIdx", "acquisitionNum",
                          "retentionTime", "msLevel", "precursorScanNum")


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
                           "GaussModel",
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
                                 "checkParams",
                                 "SSgaussStats"
                                 ),
                            .peak_function_list)

#' MessageOutput
#' @noRd
MessageOutput <- function(mes = NULL, ecol = "\n", progress =NULL, SuppressWeb = FALSE) {
  if (!is.null(mes)) {
    if (.on.public.web & !SuppressWeb) {
      # write down message
      write.table(
        mes,
        file = "metaboanalyst_spec_proc.txt",
        append = TRUE,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE,
        eol = ecol
      )
    } else {
      # print message
      if(ecol == "\n"){
        message(mes)
      } else {
        message(mes, appendLF = FALSE)
      }
      
    }
  }
  
  # write down progress
  if (.on.public.web & !is.null(progress)) {
    progress <- as.numeric(progress)
    
    write.table(
      progress,
      file = "log_progress.txt",
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
}

#' fast.write.csv
#' @author Jeff xia
#' @noRd
#' @importFrom data.table fwrite
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
      warning(e);
      write.csv(dat, file, row.names=row.names);   
    }, warning=function(w){
      warning(w);
      write.csv(dat, file, row.names=row.names); 
    });
}

#' SetGlobalParallel
#' @description SetGlobalParallel used to set the global core numbers
#' @param ncore Numeric, used to set the global core numbers, default is 1
#' @export
#' @import BiocParallel
#' @examples
#' SetGlobalParallel(1);
#' register(bpstop());
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)

SetGlobalParallel <- function(ncore = 1){
  
  if (.Platform$OS.type == "unix") {
    register(bpstart(MulticoreParam(ceiling(ncore))))
  } else if (.Platform$OS.type == "windows") {
    register(bpstart(SnowParam(ceiling(ncore))))
  }
}

