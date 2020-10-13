# Reguler part ----
# library(OptiLCMS);
# mSet<-InitDataObjects("spec", "raw", FALSE)
# mSet <- PerformROIExtraction("/home/glassfish/projects/MetaboDemoRawData/upload/QC/",rt.idx = 0.9);
# param_initial <- SetPeakParam(platform = "UPLC-Q/E");
# param_optimized <- PerformParamsOptimization(mSet, ncore = 2)

# mSet <- ImportRawMSData(mSet, "/home/glassfish/projects/MetaboDemoRawData/upload/",ncores = 10, plotSettings = SetPlotParam(Plot = F))
# mSet@params <- updateRawSpectraParam (param_optimized$best_parameters);
# mSet <- PerformPeakProfiling(mSet,ncore = 10,plotSettings = SetPlotParam(Plot = F))
# annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.015)
# mSet <- PerformPeakAnnotation(mSet, annParams)
# maPeaks <- FormatPeakList(mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1)

# remove.packages("OptiLCMS", lib="~/R/x86_64-pc-linux-gnu-library/4.0")

# Resumming Running ----
# MessageOutput <- OptiLCMS:::MessageOutput
# .on.public.web <- F
# setwd("/home/qiang/Data_IBD/test/upload/trimmed/")
# plot.MS_3D <- OptiLCMS:::plot.MS_3D
# library(lattice)
# PerformDataInspect("data/QC/Trimmed_0020a_XAV_iHMP2_FFA_PREFA02.mzML",res = 35)
# a1 <- PerformROIExtraction("/home/qiang/Data_IBD/test/upload/",mode = "rt_specific",write = T,rt = seq(480, 660, 0.2),
#                            rtdiff = 0.101,rmConts = F,plot = F)
# # 
# library(OptiLCMS);
# plan <- InitializaPlan("raw_opt","/home/glassfish/projects/MetaboDemoRawData/")
# plan <- running.plan(plan,
#                      data_folder_QC <- 'upload/QC/',
#                      mSet <- PerformROIExtraction(datapath = data_folder_QC, rt.idx = 0.95, plot = F, rmConts = F, running.controller = rc),
#                      param_initial <- SetPeakParam(),
#                      best_parameters <- PerformParamsOptimization(mSet = mSet, param_initial, ncore = 2, running.controller = rc),
#                      data_folder_Sample <- 'upload/',
#                      param <- best_parameters,
#                      plotSettings1 <- SetPlotParam(Plot=T),
#                      plotSettings2 <- SetPlotParam(Plot=T),
#                      mSet <- ImportRawMSData(mSet = mSet, foldername = data_folder_Sample, plotSettings = plotSettings1, running.controller = rc),
#                      mSet <- PerformPeakProfiling(mSet = mSet, Params = param, plotSettings = plotSettings2, running.controller = rc),
#                      annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.035),
#                      mSet <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams, ncore =1, running.controller = rc),
#                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1));
# ExecutePlan(plan);
# #
# plan <- InitializaPlan("raw_ms","/home/glassfish/projects/MetaboDemoRawData/")
# plan <- running.plan(plan,
#                      data_folder_Sample <- 'upload/',
#                      param <- SetPeakParam(),
#                      plotSettings1 <- SetPlotParam(Plot=T),
#                      plotSettings2 <- SetPlotParam(Plot=T),
#                      mSet <- ImportRawMSData(foldername = data_folder_Sample, plotSettings = plotSettings1, running.controller = rc),
#                      mSet <- PerformPeakProfiling(mSet = mSet, Params = param, plotSettings = plotSettings2, running.controller = rc),
#                      annParams <- SetAnnotationParam(polarity = 'positive', mz_abs_add = 0.035),
#                      mSet <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams, ncore =1, running.controller = rc),
#                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1));
# 
# ExecutePlan(plan)


#Or download from this link: https://drive.google.com/file/d/1CjEPed1WZrwd5T3Ovuic1KVF-Uz13NjO/view?usp=sharing
# load("mSet.rda") # load mSet for further analysis, e.g. Statistics/ Mummichog etc.



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
# param_initial <- SetPeakParam(Peak_method = "matchedFilter",RT_method = "obiwarp");
# 
# mSet <- PerformROIExtraction("/home/glassfish/projects/MetaboDemoRawData/upload/QC/",rt.idx = 0.9);
# param_optimized <- PerformParamsOptimization(mSet, ncore = 8,param = param_initial)
# 
# mSet <- ImportRawMSData(mSet, "/home/glassfish/projects/MetaboDemoRawData/upload/",ncores = 6, plotSettings = SetPlotParam(Plot = F))
# mSet <- PerformPeakProfiling(mSet,ncore = 10,Params = param_initial, plotSettings = SetPlotParam(Plot = F))
# annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.015)
# mSet <- PerformPeakAnnotation(mSet, annParams)
# maPeaks <- FormatPeakList(mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1)

# importFrom("grDevices", "boxplot.stats", "dev.off", "jpeg")
# importFrom("graphics", "abline", "boxplot", "contour", "grid",
#            "legend", "par", "points")
# importFrom("stats", "approx", "approxfun", "as.formula", "convolve",
#            "cor", "cor.test", "cutree", "deriv3", "dist", "dnorm",
#            "fft", "fitted", "hclust", "kmeans", "lm", "loess", "lsfit",
#            "na.omit", "nextn", "nls", "prcomp", "predict", "sd",
#            "smooth.spline", "weighted.mean")
# 






