#' @title testOptiLCMS
#' @description testOptiLCMS
#' @author zhiqiang Pang
#' @export
testOptiLCMS <- function(){
  print("Running Successfully !")
}

MessageOutput(
  mes = paste0("Run isotope peak annotation.."),
  ecol = "\n",
  progress = NULL
)

.optimize_switch <<- .on.public.web <<- F;
Params <- SetPeakParam();
library(BiocParallel)
mSet<-InitDataObjects("spec", "raw", FALSE)
rawfilenmsIncluded <- c('0-5_16_1_vq_002.mzXML','0-5_16_1_vq_006.mzXML','0-5_16_1_vq_004.mzXML','0-13_100q3_002.mzXML','0-1_100v8_006.mzXML','0-1_100v8_004.mzXML','0-13_100q3_006.mzXML','0-13_100q3_004.mzXML','0-1_100v8_002.mzXML','0-11_1_256_vq_006.mzXML','0-11_1_256_vq_002.mzXML','0-11_1_256_vq_004.mzXML','0-6_4_1_vq_004.mzXML','0-6_4_1_vq_002.mzXML','0-6_4_1_vq_006.mzXML','0-10_1_64_vq_002.mzXML','0-10_1_64_vq_006.mzXML','0-10_1_64_vq_004.mzXML','0-7_1_1_vq_004.mzXML','0-7_1_1_vq_002.mzXML','0-7_1_1_vq_006.mzXML','0-2_1024_1_vq_004.mzXML','0-2_1024_1_vq_006.mzXML','0-2_1024_1_vq_002.mzXML','0-9_1_16_vq_006.mzXML','0-9_1_16_vq_004.mzXML','0-9_1_16_vq_002.mzXML','0-4_64_1_vq_002.mzXML','0-4_64_1_vq_004.mzXML','0-4_64_1_vq_006.mzXML','0-3_256_1_vq_004.mzXML','0-3_256_1_vq_002.mzXML','0-3_256_1_vq_006.mzXML','0-12_1_1024_vq_006.mzXML','0-12_1_1024_vq_002.mzXML','0-12_1_1024_vq_004.mzXML','0-8_1_4_vq_004.mzXML','0-8_1_4_vq_002.mzXML','0-8_1_4_vq_006.mzXML')
mSet <- updateSpectraFiles(mSet, "/home/glassfish/projects/MetaboDemoRawData");
mSet <- ImportRawMSData(mSet, "/home/glassfish/projects/MetaboDemoRawData/upload/",ncores = 4, plotSettings = SetPlotParam(Plot = F))
mSet@params <- updateRawSpectraParam (Params)
mSet <- PerformPeakPicking(mSet = mSet)
findEqualGreaterM <- MetaboAnalystR:::findEqualGreaterM;
Densitygrouping_slave <- MetaboAnalystR:::Densitygrouping_slave;
rectUnique <- MetaboAnalystR:::rectUnique;
na.flatfill <- MetaboAnalystR:::na.flatfill
mSet <- PerformPeakGrouping(mSet = mSet)
mSet<-PerformRTcorrection(mSet)

