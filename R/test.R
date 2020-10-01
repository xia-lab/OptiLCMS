# Reguler part ----
# library(OptiLCMS);
# mSet<-InitDataObjects("spec", "raw", FALSE)
# mSet <- PerformROIExtraction("/home/glassfish/projects/MetaboDemoRawData/upload/QC/",rt.idx = 0.9);
# param_initial <- SetPeakParam(platform = "UPLC-Q/E");
# param_optimized <- PerformParamsOptimization(mSet, ncore = 8)

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
# 
# library(OptiLCMS);
# plan <- InitializaPlan("raw_opt","/home/qiang/Data_IBD/test/upload/trimmed/")
# plan <- running.plan(plan,
#                      data_folder_QC <- 'data/QC/',
#                      mSet <- PerformROIExtraction(datapath = data_folder_QC, rt.idx = 0.95, plot = F, rmConts = F, running.controller = rc),
#                      param_initial <- SetPeakParam(),
#                      best_parameters <- PerformParamsOptimization(mSet = mSet, param_initial, ncore = 8, running.controller = rc),
#                      data_folder_Sample <- 'data/',
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

# 
# library(OptiLCMS);
# plan <- InitializaPlan("raw_opt","/home/qiang/Downloads/Data_IBD/")
# plan <- running.plan(plan,
#                      data_folder_QC <- 'QC/',
#                      mSet <- PerformROIExtraction(datapath = data_folder_QC, rt.idx = 0.95, plot = F, rmConts = T, running.controller = rc),
#                      param_initial <- SetPeakParam(),
#                      best_parameters <- PerformParamsOptimization(mSet = mSet, param_initial, ncore = 4, running.controller = rc),
#                      data_folder_Sample <- '',
#                      param <- best_parameters,
#                      plotSettings1 <- SetPlotParam(Plot=T),
#                      plotSettings2 <- SetPlotParam(Plot=T),
#                      mSet <- ImportRawMSData(mSet = mSet, foldername = data_folder_Sample, plotSettings = plotSettings1, running.controller = rc),
#                      mSet <- PerformPeakProfiling(mSet = mSet, Params = param, plotSettings = plotSettings2, running.controller = rc),
#                      annParams <- SetAnnotationParam(polarity = 'negative', mz_abs_add = 0.025),
#                      mSet <- PerformPeakAnnotation(mSet = mSet, annotaParam = annParams, ncore =1, running.controller = rc),
#                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F, filtAdducts = FALSE,missPercent = 1));
# ExecutePlan(plan)
# 


# load("mSet.rda") # load mSet for further analysis, e.g. Statistics/ Mummichog etc.