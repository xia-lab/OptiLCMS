#' @title testOptiLCMS
#' @description testOptiLCMS
#' @author zhiqiang Pang
#' @export
testOptiLCMS <- function(){
  print("Running Successfully !")
}

# MessageOutput(
#   mes = paste0("Run isotope peak annotation.."),
#   ecol = "\n",
#   progress = NULL
# )
# 
# .running.as.plan <- .optimize_switch <<- .on.public.web <<- F;
# plan_count <<- 1;
# Params <- SetPeakParam();
# library(BiocParallel);library(MetaboAnalystR);library(progress);library(MSnbase)
# mSet<-InitDataObjects("spec", "raw", FALSE)
# rawfilenmsIncluded <- c('0-5_16_1_vq_002.mzXML','0-5_16_1_vq_006.mzXML','0-5_16_1_vq_004.mzXML','0-13_100q3_002.mzXML','0-1_100v8_006.mzXML','0-1_100v8_004.mzXML','0-13_100q3_006.mzXML','0-13_100q3_004.mzXML','0-1_100v8_002.mzXML','0-11_1_256_vq_006.mzXML','0-11_1_256_vq_002.mzXML','0-11_1_256_vq_004.mzXML','0-6_4_1_vq_004.mzXML','0-6_4_1_vq_002.mzXML','0-6_4_1_vq_006.mzXML','0-10_1_64_vq_002.mzXML','0-10_1_64_vq_006.mzXML','0-10_1_64_vq_004.mzXML','0-7_1_1_vq_004.mzXML','0-7_1_1_vq_002.mzXML','0-7_1_1_vq_006.mzXML','0-2_1024_1_vq_004.mzXML','0-2_1024_1_vq_006.mzXML','0-2_1024_1_vq_002.mzXML','0-9_1_16_vq_006.mzXML','0-9_1_16_vq_004.mzXML','0-9_1_16_vq_002.mzXML','0-4_64_1_vq_002.mzXML','0-4_64_1_vq_004.mzXML','0-4_64_1_vq_006.mzXML','0-3_256_1_vq_004.mzXML','0-3_256_1_vq_002.mzXML','0-3_256_1_vq_006.mzXML','0-12_1_1024_vq_006.mzXML','0-12_1_1024_vq_002.mzXML','0-12_1_1024_vq_004.mzXML','0-8_1_4_vq_004.mzXML','0-8_1_4_vq_002.mzXML','0-8_1_4_vq_006.mzXML')
# mSet <- updateSpectraFiles(mSet, "/home/glassfish/projects/MetaboDemoRawData");
# 
# mSet <- ImportRawMSData(mSet, "/home/glassfish/projects/MetaboDemoRawData/upload/",ncores = 10, plotSettings = SetPlotParam(Plot = F))
# mSet@params <- updateRawSpectraParam (Params)
# PeakPicking_centWave_slave <- MetaboAnalystR:::PeakPicking_centWave_slave
# 
# mSet <- PerformPeakPicking(mSet = mSet)
# findEqualGreaterM <- MetaboAnalystR:::findEqualGreaterM;
# Densitygrouping_slave <- MetaboAnalystR:::Densitygrouping_slave;
# rectUnique <- MetaboAnalystR:::rectUnique;
# na.flatfill <- MetaboAnalystR:::na.flatfill
# 
# mSet <- PerformPeakGrouping(mSet = mSet);
# mSet<-PerformRTcorrection(mSet);
# mSet <- PerformPeakGrouping(mSet = mSet);
# mSet <- PerformPeakFiling(mSet = mSet);
# 
# mSet <- PerformPeakProfiling(mSet,ncore = 10,plotSettings = SetPlotParam(Plot = F))
# 
# save(mSet, file = "mSet.rda")
# 
# annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.015)
# annotPeaks2 <- PerformPeakAnnotation(mSet, annParams,running.controller=NULL)
# 
# load("mSet.rda");library(OptiLCMS)
# mSet <- PerformROIExtraction("/home/glassfish/projects/MetaboDemoRawData/upload/QC/",rt.idx = 0.9);
# param_initial <- SetPeakParam(platform = "UPLC-Q/E");
# param_optimized <- PerformParamsOptimization(mSet, param = SetPeakParam(), ncore = 8)

