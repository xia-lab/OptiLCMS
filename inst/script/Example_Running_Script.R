# Example Running Tutorail ----------------------------------------------

## load googledrive package to download example data
# library("googledrive");

# Set data folder
# data_folder_Sample <- "~/Data_IBD";
# temp <- tempfile(fileext = ".zip");

# Please authorize the package to download the data from web
# dl <- drive_download(as_id("1CjEPed1WZrwd5T3Ovuic1KVF-Uz13NjO"), path = temp, overwrite = TRUE);
# out <- unzip(temp, exdir = data_folder_Sample);
# out;
# #### Running as resumable procedure: seamless pipeline
## Initialize running plan
# plan <- InitializaPlan("raw_opt","~/Data_IBD/")
## define/set running plan
# plan <- running.plan(plan,
#                      data_folder_QC <- "/home/glassfish/projects/MetaboDemoRawData/upload/QC/",
#                      mSet <- PerformROIExtraction(datapath = data_folder_QC,
#                                                   rt.idx = 0.9, plot = F,
#                                                   rmConts = F,
#                                                   running.controller = rc),
#                      param_initial <- SetPeakParam(),
#                      best_parameters <- PerformParamsOptimization(mSet = mSet,
#                                                   param_initial, ncore = 8,
#                                                   running.controller = rc),
#                      data_folder_Sample <- '/home/glassfish/projects/MetaboDemoRawData/upload/',
#                      param <- best_parameters,
#                      plotSettings1 <- SetPlotParam(Plot=T),
#                      plotSettings2 <- SetPlotParam(Plot=T),
#                      mSet <- ImportRawMSData(mSet = mSet,
#                                              path = data_folder_Sample,
#                                              plotSettings = plotSettings1,
#                                              running.controller = rc),
#                      mSet <- PerformPeakProfiling(mSet = mSet,
#                                              Params = param,
#                                              plotSettings = plotSettings2,
#                                              running.controller = rc),
#                      annParams <- SetAnnotationParam(polarity = 'negative',
#                                              mz_abs_add = 0.025),
#                      mSet <- PerformPeakAnnotation(mSet = mSet,
#                                              annotaParam = annParams,
#                                              ncore =1,
#                                              running.controller = rc),
#                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F,
#                                              filtAdducts = FALSE, missPercent = 1));
# # Execute the defined plan
# ExecutePlan(plan)

#' # revise running plan, for example, revise mz_abs_add as 0.030
# plan <- running.plan(plan,
#                      data_folder_QC <- "~/Data_IBD/QC",
#                      mSet <- PerformROIExtraction(datapath = data_folder_QC,
#                                                   rt.idx = 0.95, plot = F,
#                                                   rmConts = F,
#                                                   running.controller = rc),
#                      param_initial <- SetPeakParam(),
#                      best_parameters <- PerformParamsOptimization(mSet = mSet,
#                                                   param_initial, ncore = 2,
#                                                   running.controller = rc),
#                      data_folder_Sample <- '',
#                      param <- best_parameters,
#                      plotSettings1 <- SetPlotParam(Plot=T),
#                      plotSettings2 <- SetPlotParam(Plot=T),
#                      mSet <- ImportRawMSData(mSet = mSet,
#                                              foldername = data_folder_Sample,
#                                              plotSettings = plotSettings1,
#                                              running.controller = rc),
#                      mSet <- PerformPeakProfiling(mSet = mSet,
#                                              Params = param,
#                                              plotSettings = plotSettings2,
#                                              running.controller = rc),
#                      annParams <- SetAnnotationParam(polarity = 'negative',
#                                              mz_abs_add = 0.030),
#                      mSet <- PerformPeakAnnotation(mSet = mSet,
#                                              annotaParam = annParams,
#                                              ncore =1,
#                                              running.controller = rc),
#                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F,
#                                              filtAdducts = FALSE,
#                                              missPercent = 1));
## Re-execute the defined plan, unnecessary steps will be skipped
# ExecutePlan(plan)