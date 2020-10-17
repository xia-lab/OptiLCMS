
# 
# library(OptiLCMS);
# mSet<-InitDataObjects("spec", "raw", FALSE)
# param_initial <- SetPeakParam(mzdiff = 0.01,bw = 5,
#                               min_peakwidth = 5,max_peakwidth = 20,
#                               ppm = 15,noise = 0,snthresh = 6,
#                               prefilter = 3,value_of_prefilter = 100,
#                               minFraction = 0.85,RT_method = "obiwarp"
# );
# mSet <- PerformROIExtraction("/home/glassfish/projects/MetaboDemoRawData/upload/QC/",rt.idx = 0.9);
# param_optimized <- PerformParamsOptimization(mSet, ncore = 8,param = param_initial)
# mSet <- ImportRawMSData(mSet, "/home/glassfish/projects/MetaboDemoRawData/upload/",ncores = 6, plotSettings = SetPlotParam(Plot = F))
# mSet <- PerformPeakProfiling(mSet,ncore = 10,Params = param_optimized, plotSettings = SetPlotParam(Plot = F))
# annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.015)
# mSet <- PerformPeakAnnotation(mSet, annParams)
# maPeaks <- FormatPeakList(mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1)


# library(OptiLCMS);
# mSet<-InitDataObjects("spec", "raw", FALSE)
# param_initial <- SetPeakParam(Peak_method = "Massifquant",RT_method = "obiwarp");
# 
# mSet <- PerformROIExtraction("/home/glassfish/projects/MetaboDemoRawData/upload/QC/",rt.idx = 0.9);
# param_optimized <- PerformParamsOptimization(mSet, ncore = 8,param = param_initial)
# 
# mSet <- ImportRawMSData(mSet, "/home/qiang/Downloads/cdf_example/test_covid/",ncores = 6, plotSettings = SetPlotParam(Plot = F))
# mSet <- PerformPeakProfiling(mSet,ncore = 6,Params = param_initial, plotSettings = SetPlotParam(Plot = F))
# annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.015)
# mSet <- PerformPeakAnnotation(mSet, annParams)
# maPeaks <- FormatPeakList(mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1)


# tempZip <- tempfile(fileext = ".zip");
# download.file("https://www.dropbox.com/s/kabienyoadmzdpm/SpectraData.zip", destfile = tempZip, method = "wget");
# dataDir <- paste0(tempdir(),"/SpectraData");
# out <- unzip(tempZip, exdir = dataDir);
# 
# library(OptiLCMS);
# setwd("/home/glassfish/projects/MetaboDemoRawData/")
# plan <- InitializaPlan("raw_opt")
# plan <- running.plan(plan,
#                      data_folder_QC <- paste0(dataDir,"/QC"),
#                      mSet <- PerformROIExtraction(datapath = data_folder_QC, rt.idx = 0.5, plot = F, rmConts = F, running.controller = rc),
#                      param_initial <- SetPeakParam(),
#                      best_parameters <- PerformParamsOptimization(mSet = mSet, param_initial, ncore = 1, running.controller = rc),
#                      param <- best_parameters,
#                      plotSettings1 <- SetPlotParam(Plot=T),
#                      plotSettings2 <- SetPlotParam(Plot=T),
#                      mSet <- ImportRawMSData(mSet = mSet, foldername = dataDir, plotSettings = plotSettings1, running.controller = rc),
#                      mSet <- PerformPeakProfiling(mSet = mSet, Params = param, plotSettings = plotSettings2, running.controller = rc),
#                      annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.025),
#                      mSet <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams, ncore =1, running.controller = rc),
#                      mSet <- FormatPeakList(mSet = mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1));
# result <- ExecutePlan(plan);
# 
# plan <- running.plan(plan,
#                      data_folder_QC <- paste0(dataDir,"/QC"),
#                      mSet <- PerformROIExtraction(datapath = data_folder_QC, rt.idx = 0.5, plot = F, rmConts = F, running.controller = rc),
#                      param_initial <- SetPeakParam(),
#                      best_parameters <- PerformParamsOptimization(mSet = mSet, param_initial, ncore = 1, running.controller = rc),
#                      param <- best_parameters,
#                      plotSettings1 <- SetPlotParam(Plot=T),
#                      plotSettings2 <- SetPlotParam(Plot=T),
#                      mSet <- ImportRawMSData(mSet = mSet, foldername = dataDir, plotSettings = plotSettings1, running.controller = rc),
#                      mSet <- PerformPeakProfiling(mSet = mSet, Params = param, plotSettings = plotSettings2, running.controller = rc),
#                      annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.035),
#                      mSet <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams, ncore =1, running.controller = rc),
#                      mSet <- FormatPeakList(mSet = mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1));
# result <- ExecutePlan(plan);
# #

#' ## Download the raw spectra example data
#' tempZip <- tempfile(fileext = ".zip");
#' download.file("https://www.dropbox.com/s/kabienyoadmzdpm/SpectraData.zip",
#'               destfile = tempZip, method = "wget");
#' dataDir <- paste0(tempdir(),"/SpectraData");
#' out <- unzip(tempZip, exdir = dataDir);
#' 
#' ## Load OptiLCMS
#' library(OptiLCMS);
#' 
#' ## Initialize your plan
#' mSet <- ImSet<-InitDataObjects("spec", "raw", FALSE);
#' 
#' ## Import raw data
#' mSet <- ImportRawMSData(mSet = mSet, foldername = dataDir,
#'                         plotSettings = SetPlotParam(Plot=T));
#' 
#' ## Perform raw spectra profiling
#' mSet <- PerformPeakProfiling(mSet = mSet, Params = SetPeakParam(),
#'                              plotSettings = SetPlotParam(Plot=T), 
#'                              ncore = 1);
#' 
#' ## Set annotation parameters and run
#' annParams <- SetAnnotationParam(polarity = 'negative',
#'                                 mz_abs_add = 0.025);
#' mSet <- PerformPeakAnnotation(mSet = mSet,
#'                               annotaParam = annParams, 
#'                               ncore =1);
#' ## Format the PeakList
#' mSet <- FormatPeakList(mSet = mSet, 
#'                        annParams,
#'                        filtIso =F, 
#'                        filtAdducts = FALSE,
#'                        missPercent = 1)
#' 
#' ## Export the annotation result
#' Export.Annotation(mSet);
#' 
#' ## Export the Peak Table
#' Export.PeakTable(mSet);
#' 
#' ## Export the Peak summary
#' Export.PeakSummary(mSet)




# importFrom("grDevices", "boxplot.stats", "dev.off", "jpeg")
# importFrom("graphics", "abline", "boxplot", "contour", "grid",
#            "legend", "par", "points")
# importFrom("stats", "approx", "approxfun", "as.formula", "convolve",
#            "cor", "cor.test", "cutree", "deriv3", "dist", "dnorm",
#            "fft", "fitted", "hclust", "kmeans", "lm", "loess", "lsfit",
#            "na.omit", "nextn", "nls", "prcomp", "predict", "sd",
#            "smooth.spline", "weighted.mean")
# 






