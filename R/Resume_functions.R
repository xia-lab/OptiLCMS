## Perform Running Pipeline - script

#' Initializing running plan
#' @param type Initialized plan type for a resumable running mode. Can be "raw_opt" for automated optimization option, or "raw_ms" for customized pipeline.
#' @param specDir this parameter is used to introduce the spectral directory (where the raw spectra data exists).
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
#' @export
#' @examples 
#' ## load googledrive package to download example data
#' # library("googledrive");
#'
#' # Set data folder
#' # data_folder_Sample <- "~/Data_IBD";
#' # temp <- tempfile(fileext = ".zip");
#'
#' # Please authorize the package to download the data from web
#' # dl <- drive_download(as_id("1CjEPed1WZrwd5T3Ovuic1KVF-Uz13NjO"), path = temp, overwrite = TRUE);
#' # out <- unzip(temp, exdir = data_folder_Sample);
#' # out;
# 
#' #### Running as resumable procedure: seamless pipeline
#' ## Initialize running plan
#' # plan <- InitializaPlan("raw_opt","~/Data_IBD/")
#' ## define/set running plan
#' # plan <- running.plan(plan,
#' #                      data_folder_QC <- "~/Data_IBD/QC",
#' #                      mSet <- PerformROIExtraction(datapath = data_folder_QC, 
#' #                                                   rt.idx = 0.95, plot = F, 
#' #                                                   rmConts = F, 
#' #                                                   running.controller = rc),
#' #                      param_initial <- SetPeakParam(),
#' #                      best_parameters <- PerformParamsOptimization(mSet = mSet, 
#' #                                                   param_initial, ncore = 2, 
#' #                                                   running.controller = rc),
#' #                      data_folder_Sample <- '',
#' #                      param <- best_parameters,
#' #                      plotSettings1 <- SetPlotParam(Plot=T),
#' #                      plotSettings2 <- SetPlotParam(Plot=T),
#' #                      mSet <- ImportRawMSData(mSet = mSet, 
#' #                                              foldername = data_folder_Sample, 
#' #                                              plotSettings = plotSettings1, 
#' #                                              running.controller = rc),
#' #                      mSet <- PerformPeakProfiling(mSet = mSet, 
#' #                                              Params = param, 
#' #                                              plotSettings = plotSettings2, 
#' #                                              running.controller = rc),
#' #                      annParams <- SetAnnotationParam(polarity = 'negative', 
#' #                                              mz_abs_add = 0.025),
#' #                      mSet <- PerformPeakAnnotation(mSet = mSet, 
#' #                                              annotaParam = annParams, 
#' #                                              ncore =1, 
#' #                                              running.controller = rc),
#' #                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F, 
#' #                                              filtAdducts = FALSE, missPercent = 1));
#' ## Execute the defined plan
#' # ExecutePlan(plan)

InitializaPlan <- function(type="spec", specDir=NULL){
  
  if(!is.null(specDir)){
    #setwd(WorkingDir);
    fullUserPath <- specDir;
  } else {
    fullUserPath <- getwd();
  };
  
  if(!exists(".SwapEnv")){
    .SwapEnv <<- new.env(parent = .GlobalEnv);
    .SwapEnv$.optimize_switch <- FALSE;
    .SwapEnv$count_current_sample <- 0;
    .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
    .SwapEnv$envir <- new.env();
  }
  
  MessageOutput(paste0("Running Status -- Plan Initialized Successfully at: ", Sys.time(), "\nPlease define your running plan ..."), "\n", 0);
  
  temp.path <- paste0(getwd(), "/temp/");
  
  if(dir.exists(temp.path)){
    unlink(temp.path, recursive = TRUE, force = TRUE);
  }
  
  #=============== raw_opt
  if (type == "raw_opt") {
    
    plan.path <- paste0(getwd(), "/temp/plan");
    
    if (!dir.exists(plan.path)) {
      dir.create(paste0(getwd(), "/temp/plan"),
                 recursive = TRUE)
    }
    
    #.running.as.plan <- .SwapEnv$.running.as.plan <- TRUE;
    plan <- new("ResumingePlan");
    .plan_count <- plan@PlanNumber <- 0;
    plan@WorkingDir <- fullUserPath;
    
    saveRDS(plan, file = paste0(plan.path, "/plan.rds"));
    #saveRDS(.running.as.plan, file = paste0(plan.path, "/running.as.plan.rds"));
    saveRDS(.plan_count, file = paste0(plan.path, "/plan_count.rds"));
    #----------------------
    .optimize_switch <- .SwapEnv$.optimize_switch <- TRUE;
    switch.path <- paste0(getwd(), "/temp/plan")
    saveRDS(.optimize_switch,
            file = paste0(switch.path, "/optimize_switch_", .plan_count, ".rds"))
    
    #----------------------
    envir.path <- paste0(getwd(), "/temp/envir")
    if (!dir.exists(envir.path)) {
      dir.create(paste0(getwd(), "/temp/envir"),
                 recursive = T)
    }
    
    envir <- .SwapEnv$envir <- new.env()
    saveRDS(envir, file = paste0(envir.path, "/envir.rds"))
    #----------------------
    # if (.on.public.web) {
    #   rawFileNames <- paste0(getwd(), "/temp/plan")
    #   saveRDS(rawfilenms,
    #           file = paste0(rawFileNames, "/rawfilenms_", .plan_count, ".rds"))
    # }
    #----------------------
    
    if (exists(".plan_count") & .plan_count > 0) {
      return(plan)
    }
  }
  
  #================== raw_ms
  if (type == "raw_ms") {
    plan.path <- paste0(getwd(), "/temp/plan")
    
    if (!dir.exists(plan.path)) {
      dir.create(paste0(getwd(), "/temp/plan"),
                 recursive = T)
    }
    
    #---------------
    #.running.as.plan <- .SwapEnv$.running.as.plan <- TRUE;
    plan <- new("ResumingePlan");
    .plan_count <- plan@PlanNumber <- 0;
    plan@WorkingDir <- fullUserPath;
    
    saveRDS(plan, file = paste0(plan.path, "/plan.rds"));
    #saveRDS(.running.as.plan, file = paste0(plan.path, "/running.as.plan.rds"));
    saveRDS(.plan_count, file = paste0(plan.path, "/plan_count.rds"));
    #----------------------
    .optimize_switch <- .SwapEnv$.optimize_switch <- FALSE;
    switch.path <- paste0(getwd(), "/temp/plan");
    saveRDS(.optimize_switch,
            file = paste0(switch.path, "/optimize_switch_", .plan_count, ".rds"));
    
    #----------------------
    envir.path <- paste0(getwd(), "/temp/envir")
    if (!dir.exists(envir.path)) {
      dir.create(paste0(getwd(), "/temp/envir"),
                 recursive = T)
    }
    
    envir <- .SwapEnv$envir <- new.env()
    saveRDS(envir, file = paste0(envir.path, "/envir.rds"))
    #----------------------
    # rawFileNames <- paste0(getwd(), "/temp/plan")
    # 
    # if (.on.public.web) {
    #   saveRDS(rawfilenms,
    #           file = paste0(rawFileNames, "/rawfilenms_", .plan_count, ".rds"))
    # } else {
    #   # do nothing for local
    # }
    # 
    #----------------------
    if (exists(".plan_count") & .plan_count > 0) {
      return(plan)
    }
  }
  
  ## Recording cache file information
  record.info <- matrix(nrow = 30, ncol = 2)
  record.path <- paste0(getwd(), "/temp/records")
  
  if (!dir.exists(record.path)) {
    dir.create(paste0(getwd(), "/temp/records"),
               recursive = T)
  }
  
  saveRDS(record.info, file = paste0(record.path, "/records.rds"))
  
  return(plan)
  #return(.set.mSet(plan))
}

#' running.plan
#' @param plan ResummingPlan object. The object is generated by 'InitializaPlan' function.
#' @param ... Multiple Processing commands can be input here.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
#' @export
#' @examples 
#' ## load googledrive package to download example data
#' # library("googledrive");
#'
#' # Set data folder
#' # data_folder_Sample <- "~/Data_IBD";
#' # temp <- tempfile(fileext = ".zip");
#'
#' # Please authorize the package to download the data from web
#' # dl <- drive_download(as_id("1CjEPed1WZrwd5T3Ovuic1KVF-Uz13NjO"), path = temp, overwrite = TRUE);
#' # out <- unzip(temp, exdir = data_folder_Sample);
#' # out;
# 
#' #### Running as resumable procedure: seamless pipeline
#' ## Initialize running plan
#' # plan <- InitializaPlan("raw_opt","~/Data_IBD/")
#' ## define/set running plan
#' # plan <- running.plan(plan,
#' #                      data_folder_QC <- "~/Data_IBD/QC",
#' #                      mSet <- PerformROIExtraction(datapath = data_folder_QC, 
#' #                                                   rt.idx = 0.95, plot = F, 
#' #                                                   rmConts = F, 
#' #                                                   running.controller = rc),
#' #                      param_initial <- SetPeakParam(),
#' #                      best_parameters <- PerformParamsOptimization(mSet = mSet, 
#' #                                                   param_initial, ncore = 2, 
#' #                                                   running.controller = rc),
#' #                      data_folder_Sample <- '',
#' #                      param <- best_parameters,
#' #                      plotSettings1 <- SetPlotParam(Plot=T),
#' #                      plotSettings2 <- SetPlotParam(Plot=T),
#' #                      mSet <- ImportRawMSData(mSet = mSet, 
#' #                                              foldername = data_folder_Sample, 
#' #                                              plotSettings = plotSettings1, 
#' #                                              running.controller = rc),
#' #                      mSet <- PerformPeakProfiling(mSet = mSet, 
#' #                                              Params = param, 
#' #                                              plotSettings = plotSettings2, 
#' #                                              running.controller = rc),
#' #                      annParams <- SetAnnotationParam(polarity = 'negative', 
#' #                                              mz_abs_add = 0.025),
#' #                      mSet <- PerformPeakAnnotation(mSet = mSet, 
#' #                                              annotaParam = annParams, 
#' #                                              ncore =1, 
#' #                                              running.controller = rc),
#' #                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F, 
#' #                                              filtAdducts = FALSE, missPercent = 1));
#' ## Execute the defined plan
#' # ExecutePlan(plan)
#' 
#' #' # revise running plan, for example, revise mz_abs_add as 0.030
#' # plan <- running.plan(plan,
#' #                      data_folder_QC <- "~/Data_IBD/QC",
#' #                      mSet <- PerformROIExtraction(datapath = data_folder_QC, 
#' #                                                   rt.idx = 0.95, plot = F, 
#' #                                                   rmConts = F, 
#' #                                                   running.controller = rc),
#' #                      param_initial <- SetPeakParam(),
#' #                      best_parameters <- PerformParamsOptimization(mSet = mSet, 
#' #                                                   param_initial, ncore = 2, 
#' #                                                   running.controller = rc),
#' #                      data_folder_Sample <- '',
#' #                      param <- best_parameters,
#' #                      plotSettings1 <- SetPlotParam(Plot=T),
#' #                      plotSettings2 <- SetPlotParam(Plot=T),
#' #                      mSet <- ImportRawMSData(mSet = mSet, 
#' #                                              foldername = data_folder_Sample, 
#' #                                              plotSettings = plotSettings1, 
#' #                                              running.controller = rc),
#' #                      mSet <- PerformPeakProfiling(mSet = mSet, 
#' #                                              Params = param, 
#' #                                              plotSettings = plotSettings2, 
#' #                                              running.controller = rc),
#' #                      annParams <- SetAnnotationParam(polarity = 'negative', 
#' #                                              mz_abs_add = 0.030),
#' #                      mSet <- PerformPeakAnnotation(mSet = mSet, 
#' #                                              annotaParam = annParams, 
#' #                                              ncore =1, 
#' #                                              running.controller = rc),
#' #                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F, 
#' #                                              filtAdducts = FALSE, 
#' #                                              missPercent = 1));
#' ## Re-execute the defined plan, unnecessary steps will be skipped
#' # ExecutePlan(plan)

running.plan <- function(plan=NULL,...){
  
  #plan <- .get.current.plan(plan);
  commands <- match.call(expand.dots = FALSE)$...
  
  ## Declare controller
  plan@running.controller <- controller.resetter();
  ##
  
  if (!length(commands)) {
    stop("No command provided to run !");
  }
  
  plan@PlanNumber <- .plan_count <- plan@PlanNumber + 1;
  
  CommandsVerified <- CommandsVerify(commands);
  MessageOutput("Commands Origanization Finished!", ecol = "\n", NULL);
  
  plan@CommandSet[[paste0("command_set_",.plan_count)]] <- CommandsVerified;
  
  plan.path <- paste0(getwd(), "/temp/plan");
  saveRDS(plan, file = paste0(plan.path, "/plan.rds"));
  saveRDS(.plan_count, file = paste0(plan.path, "/plan_count.rds"));
  
  recordMarker_resetter(.plan_count);
  
  return(plan)
}

#' @noRd
.get.current.plan <- function(plan){
  
  if (is.null(plan)){
    return(list())
  } else {
    return (plan)
  }
  
}

#' ExecutePlan
#' @param plan ResummingPlan object. The object is generated by running.plan() function.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
#' @export
#' @examples 
#' ## load googledrive package to download example data
#' # library("googledrive");
#'
#' # Set data folder
#' # data_folder_Sample <- "~/Data_IBD";
#' # temp <- tempfile(fileext = ".zip");
#'
#' # Please authorize the package to download the data from web
#' # dl <- drive_download(as_id("1CjEPed1WZrwd5T3Ovuic1KVF-Uz13NjO"), path = temp, overwrite = TRUE);
#' # out <- unzip(temp, exdir = data_folder_Sample);
#' # out;
# 
#' #### Running as resumable procedure: seamless pipeline
#' ## Initialize running plan
#' # plan <- InitializaPlan("raw_opt","~/Data_IBD/")
#' ## define/set running plan
#' # plan <- running.plan(plan,
#' #                      data_folder_QC <- "~/Data_IBD/QC",
#' #                      mSet <- PerformROIExtraction(datapath = data_folder_QC, 
#' #                                                   rt.idx = 0.95, plot = F, 
#' #                                                   rmConts = F, 
#' #                                                   running.controller = rc),
#' #                      param_initial <- SetPeakParam(),
#' #                      best_parameters <- PerformParamsOptimization(mSet = mSet, 
#' #                                                   param_initial, ncore = 2, 
#' #                                                   running.controller = rc),
#' #                      data_folder_Sample <- '',
#' #                      param <- best_parameters,
#' #                      plotSettings1 <- SetPlotParam(Plot=T),
#' #                      plotSettings2 <- SetPlotParam(Plot=T),
#' #                      mSet <- ImportRawMSData(mSet = mSet, 
#' #                                              foldername = data_folder_Sample, 
#' #                                              plotSettings = plotSettings1, 
#' #                                              running.controller = rc),
#' #                      mSet <- PerformPeakProfiling(mSet = mSet, 
#' #                                              Params = param, 
#' #                                              plotSettings = plotSettings2, 
#' #                                              running.controller = rc),
#' #                      annParams <- SetAnnotationParam(polarity = 'negative', 
#' #                                              mz_abs_add = 0.025),
#' #                      mSet <- PerformPeakAnnotation(mSet = mSet, 
#' #                                              annotaParam = annParams, 
#' #                                              ncore =1, 
#' #                                              running.controller = rc),
#' #                      maPeaks <- FormatPeakList(mSet = mSet, annParams, filtIso =F, 
#' #                                              filtAdducts = FALSE, 
#' #                                              missPercent = 1));
#' ## Execute the defined plan
#' # ExecutePlan(plan)

ExecutePlan <- function(plan=NULL){
  
  MessageOutput(paste0("Current Working Directory: ",getwd(),"\n"));
  if (is.null(plan)){
    stop("No running plan input. Please design you plan first with 'running.plan' function !")
  };
  
  .plan_count <- plan@PlanNumber;
  
  # Reset running.controller to make sure everything is normal at beginning
  plan@running.controller <- controller.resetter();
  
  if (length(plan@CommandSet) == 1){
    
    .SwapEnv$envir$rc <- plan@running.controller;
    perform.plan(plan@CommandSet[["command_set_1"]]);
    
    envir <- .SwapEnv$envir;
    envir.path <- paste0(getwd(),"/temp/envir");
    saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
    
  } else if (length(plan@CommandSet) > 1) {
    
    ## This is the most important part for the whole pipeline
    ## Module 1-6: Dectect whether the parameters changed or not so as to operate the whole pipeline
    
    last_command_set <- plan@CommandSet[[.plan_count-1]];
    new_command_set <- plan@CommandSet[[.plan_count]];
    
    if(class(last_command_set) != class(new_command_set)){
      
      plan@running.controller <- controller.resetter();
      
    } else {
      plan <-
        controller.modifier(new_command_set,
                            last_command_set,
                            plan)
    }

    ## Module 7 - Detect whether some steps have been run in last plan excuting process (Secondary Decision switch)
    plan <- recording_identifier(plan);
    
    ## Module 8 - Dectect whether current plan type (raw_ms or raw_pre) is different from the last one (Final Decision switch)
    # plan <- planType_identifier(plan);
    
    .SwapEnv$envir$rc <- plan@running.controller
    # define.plan.controller <- str2lang('rc <- plan$running.controller');
    # eval (define.plan.controller,envir = envir);
    #perform.plan(new_command_set);
    
    if(class(new_command_set) == "CustCommandSet"){
      MessageOutput(mes = "Step 1/6: Ready to start the customized pipeline !",
                    ecol = "\n",
                    progress = 10.0)
    }
    
    mSetInfo <- tryCatch(perform.plan(new_command_set), error = function(e){e});
    
    if (class(mSetInfo)[1]=="simpleError"){
      print_mes <- paste0("<font color=\"red\">","\nERROR:",mSetInfo$message,"</font>");
      write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
      stop(paste0("EXCEPTION POINT CODE: ", mSetInfo$message));
    }
    
    envir <- .SwapEnv$envir;
    envir.path <- paste0(getwd(),"/temp/envir");
    saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
    
  } else {
    stop("Wrong plan has been prepared !");
  }
  
}

#' @noRd
recording_identifier <- function(plan) {
  
  .plan_count <- plan@PlanNumber;
  record.path <- paste0(getwd(),"/temp/records");
  
  if(file.exists(paste0(record.path,"/records_marker_",.plan_count-1,".rds"))){
    record.marker_last <- readRDS(paste0(record.path,"/records_marker_",.plan_count-1,".rds"));
  } else {
    recordMarker_resetter(.plan_count-1);
    record.marker_last <- readRDS(paste0(record.path,"/records_marker_",.plan_count-1,".rds"));
  }
  
  
  # If things were not run last time, forcedly start running the slot this time
  if(!as.logical(record.marker_last[1,2])){plan@running.controller@ROI_extract[1] <- TRUE};
  if(!as.logical(record.marker_last[2,2])){plan@running.controller@ROI_extract[2] <- TRUE};
  if(!as.logical(record.marker_last[3,2])){plan@running.controller@ROI_extract[3] <- TRUE};
  if(!as.logical(record.marker_last[4,2])){plan@running.controller@ROI_extract[4] <- TRUE};
  
  if(!as.logical(record.marker_last[5,2])){plan@running.controller@data_import[1] <- TRUE}
  if(!as.logical(record.marker_last[6,2])){plan@running.controller@data_import[2] <- TRUE}
  if(!as.logical(record.marker_last[7,2])){plan@running.controller@data_import[3] <- TRUE}
  if(!as.logical(record.marker_last[8,2])){plan@running.controller@data_import[4] <- TRUE}
  
  if(!as.logical(record.marker_last[9,2])){plan@running.controller@peak_profiling[1] <- TRUE}
  if(!as.logical(record.marker_last[10,2])){plan@running.controller@peak_profiling[2] <- TRUE}
  if(!as.logical(record.marker_last[11,2])){plan@running.controller@peak_profiling[3] <- TRUE}
  if(!as.logical(record.marker_last[12,2])){plan@running.controller@peak_profiling[4] <- TRUE}
  
  if(!as.logical(record.marker_last[13,2])){plan@running.controller@peak_annotation[1] <- TRUE}
  if(!as.logical(record.marker_last[14,2])){plan@running.controller@peak_annotation[2] <- TRUE}
  if(!as.logical(record.marker_last[15,2])){plan@running.controller@peak_annotation[3] <- TRUE}
  if(!as.logical(record.marker_last[16,2])){plan@running.controller@peak_annotation[4] <- TRUE}
  
  if(!as.logical(record.marker_last[17,2])){plan@running.controller@others_1[1] <- TRUE}
  if(!as.logical(record.marker_last[18,2])){plan@running.controller@others_1[2] <- TRUE}
  if(!as.logical(record.marker_last[19,2])){plan@running.controller@others_1[3] <- TRUE}
  if(!as.logical(record.marker_last[20,2])){plan@running.controller@others_1[4] <- TRUE}
  
  if(!as.logical(record.marker_last[21,2])){plan@running.controller@operators[1] <- TRUE}
  if(!as.logical(record.marker_last[22,2])){plan@running.controller@operators[2] <- FALSE} # Switch off this operator at the end of plan define
  if(!as.logical(record.marker_last[23,2])){plan@running.controller@operators[3] <- FALSE} # Switch off this operator at the end of plan define
  if(!as.logical(record.marker_last[24,2])){plan@running.controller@operators[4] <- FALSE}
  if(!as.logical(record.marker_last[25,2])){plan@running.controller@operators[5] <- FALSE}
  if(!as.logical(record.marker_last[26,2])){plan@running.controller@operators[6] <- FALSE}
  if(!as.logical(record.marker_last[27,2])){plan@running.controller@operators[7] <- FALSE}
  if(!as.logical(record.marker_last[28,2])){plan@running.controller@operators[8] <- FALSE}
  
  plan@RunningHistory@setNO <- .plan_count-1;
  plan@RunningHistory@FunishedPosition <- which(as.logical(record.marker_last[,2]))
  
  return(plan)
  
}

#' @noRd
recordMarker_resetter <- function(.plan_count){
  
  ## Recording initial markers about being ran or not
  record.marker <- matrix(nrow = 28,ncol = 2);
  record.marker[,1] <- c("ROI_extract_1","ROI_extract_2","ROI_extract_3","ROI_extract_4",
                         "data_import_1","data_import_2","data_import_3","data_import_4",
                         "profiling_1","profiling_2","profiling_3","profiling_4",
                         "peak_annotation_1","peak_annotation_2","peak_annotation_3","peak_annotation_4",
                         "others_1","others_2","others_3","others_4",
                         "operators_1","operators_2","operators_3","operators_4",
                         "operators_5","operators_6","operators_7","operators_8");
  record.marker[,2] <- rep(FALSE, 28);
  record.path <- paste0(getwd(),"/temp/records");
  
  saveRDS(record.marker,file = paste0(record.path,"/records_marker_",.plan_count,".rds"));
  
}

#' @noRd
marker_record <- function(functionNM){
  
  record.path <- paste0(getwd(),"/temp/records");
  plan.path <- paste0(getwd(),"/temp/plan/");
  
  .plan_count <- readRDS(paste0(plan.path,"plan_count.rds"))
  
  if (!file.exists(paste0(record.path,"/records_marker_",.plan_count,".rds"))){
    record.marker <- matrix(nrow = 28,ncol = 2);
    record.marker[,1] <- c("ROI_extract_1","ROI_extract_2","ROI_extract_3","ROI_extract_4",
                           "data_import_1","data_import_2","data_import_3","data_import_4",
                           "profiling_1","profiling_2","profiling_3","profiling_4",
                           "peak_annotation_1","peak_annotation_2","peak_annotation_3","peak_annotation_4",
                           "others_1","others_2","others_3","others_4",
                           "operators_1","operators_2","operators_3","operators_4",
                           "operators_5","operators_6","operators_7","operators_8");
    record.marker[,2] <- rep(FALSE, 28);
    saveRDS(record.marker,file = paste0(record.path,"/records_marker_",.plan_count,".rds"))
  } else {
    record.marker <- readRDS(paste0(record.path,"/records_marker_",.plan_count,".rds"));
  };
  
  # If this step has been run at ".plan_count" time, is will be marked as T
  if(functionNM=="ROI_extract_c1"){record.marker[1,2] <- TRUE};
  if(functionNM=="ROI_extract_c2"){record.marker[2,2] <- TRUE};
  if(functionNM=="ROI_extract_c3"){record.marker[3,2] <- TRUE};
  if(functionNM=="ROI_extract_c4"){record.marker[4,2] <- TRUE};
  
  if(functionNM=="data_import_c1"){record.marker[5,2] <- TRUE};
  if(functionNM=="data_import_c2"){record.marker[6,2] <- TRUE};
  if(functionNM=="data_import_c3"){record.marker[7,2] <- TRUE}; # RETAINED 
  if(functionNM=="data_import_c4"){record.marker[8,2] <- TRUE}; # RETAINED 
  
  if(functionNM=="peak_profiling_c1"){record.marker[9,2] <- TRUE};
  if(functionNM=="peak_profiling_c2"){record.marker[10,2] <- TRUE};
  if(functionNM=="peak_profiling_c3"){record.marker[11,2] <- TRUE};
  if(functionNM=="peak_profiling_c4"){record.marker[12,2] <- TRUE};
  
  if(functionNM=="peak_annotation_c1"){record.marker[13,2] <- TRUE};
  if(functionNM=="peak_annotation_c2"){record.marker[14,2] <- TRUE}; # RETAINED 
  if(functionNM=="peak_annotation_c3"){record.marker[15,2] <- TRUE}; # RETAINED 
  if(functionNM=="peak_annotation_c4"){record.marker[16,2] <- TRUE}; # RETAINED 
  
  if(functionNM=="others_c1"){record.marker[17,2] <- TRUE};
  if(functionNM=="others_c2"){record.marker[18,2] <- TRUE}; # RETAINED 
  if(functionNM=="others_c3"){record.marker[19,2] <- TRUE}; # RETAINED 
  if(functionNM=="others_c4"){record.marker[19,2] <- TRUE}; # RETAINED 
  
  if(functionNM=="operators_c1"){record.marker[13,2] <- TRUE};
  if(functionNM=="operators_c2"){record.marker[14,2] <- TRUE};
  if(functionNM=="operators_c3"){record.marker[15,2] <- TRUE};
  if(functionNM=="operators_c4"){record.marker[16,2] <- TRUE};
  if(functionNM=="operators_c5"){record.marker[17,2] <- FALSE}; # RETAINED 
  if(functionNM=="operators_c6"){record.marker[18,2] <- FALSE}; # RETAINED 
  if(functionNM=="operators_c7"){record.marker[19,2] <- FALSE}; # RETAINED 
  if(functionNM=="operators_c8"){record.marker[20,2] <- FALSE}; # RETAINED 
  
  saveRDS(record.marker,file = paste0(record.path,"/records_marker_",.plan_count,".rds"));
}

#' @noRd
controller.resetter <- function() {
  running.controller <- new("controller");
  
  points <- rep(T, 5);
  names(points) <- c("c1", "c2", "c3", "c4", "c5");
  
  running.controller@ROI_extract <-
    running.controller@data_import <-
    running.controller@peak_profiling <-
    running.controller@peak_annotation <-
    running.controller@others_1 <- 
    points;
  
  running.controller@operators <- c(F, F, F, F, F, F, F, F);
  
  names(running.controller@operators) <-
    c(
      "operators_1",
      "operators_2",
      "operators_3",
      "operators_4",
      "operators_5",
      "operators_6",
      "operators_7",
      "operators_8"
    );
  
  return(running.controller)
}

controller.modifier <- function(new_command_set, last_command_set, plan){
  
  ###-------------Operators definition ------------//
  
  # operators_1 : PerformDataTrimming ~ PerformParamsOptimization;
  # operators_2 : PerformParamsOptimization ~ peak_profiling;  #' 
  # operators_3 : PerformPeakProfiling ~ PerformPeakAnnotation;  #' 
  # operators_4 : ;
  
  
  # others_1 c1 : Control Optimization or not;
  # others_1 c2 : ;
  # others_1 c3 : ;
  # others_1 c4 : .
  
  ##----------------------------------------------//
  
  if(class(new_command_set) == "OptiCommandSet"){
    
    # 1. ROIExtraction: -----
    # 1.1 Note on controller: c1, read; c2, trim; c3, write; c4, plot; C5, rmConts.
    new_ROIExtraction <- new_command_set@ROIExtraction[[3]];
    last_ROIExtraction <- last_command_set@ROIExtraction[[3]];
    ChangedArugsArray <- NULL;
    
    for (i in 2:length(new_ROIExtraction)){
      for (j in 2:length(last_ROIExtraction)){
        if (names(new_ROIExtraction[i]) == names(last_ROIExtraction[j]) & 
            new_ROIExtraction[[i]] != last_ROIExtraction[[j]]) {
          ChangedArugsArray <- c(ChangedArugsArray, names(new_ROIExtraction[i]))
        } 
      }
    }
    
    newArugsNMs <- names(new_ROIExtraction)[-1];
    lastArugsNMs <- names(last_ROIExtraction)[-1];
    NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
    RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
    
    ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray)
    
    if(is.null(ChangedArugsArray) | length(ChangedArugsArray) == 0){
      # Not changed !
      plan@running.controller@ROI_extract[c(1:5)] <- rep(FALSE, 5);
    }
    
    if("datapath" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[c(1:5)] <- rep(TRUE,5);
      plan@running.controller@operators[1] <- TRUE;
    } else {
      plan@running.controller@ROI_extract[1] <- FALSE;
    }
    
    if("rmConts" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[5] <- TRUE;
      plan@running.controller@operators[1] <- TRUE;
    }else {
      plan@running.controller@ROI_extract[5] <- FALSE;
    }
    
    if(any(c("mode", "mz", "mzdiff", "rt", "rtdiff", "rt.idx") %in% ChangedArugsArray)){
      # if any of these params changed, need to run the following steps
      plan@running.controller@ROI_extract[c(2:4)] <- rep(TRUE,3);
      plan@running.controller@operators[1] <- TRUE;
    } else if (any(plan@running.controller@ROI_extract[c(1,5)])){
      # if datapath changed, still need to run the following steps
      plan@running.controller@ROI_extract[c(2:4)] <- rep(TRUE,3);
      plan@running.controller@operators[1] <- TRUE;
    } else {
      # if all these params and datapath did not changed, skipped trimming step
      plan@running.controller@ROI_extract[2] <- FALSE;
    }
    
    if("write" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[3] <- TRUE;
    }
    
    if("plot" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[4] <- TRUE;
    }
    
    # 2. ParamsOptimization: -----
    # 2.1 Note on controller: others c1, run or not.
    new_ParamsOptimization <- new_command_set@ParamsOptimization[[3]];
    last_ParamsOptimization <- last_command_set@ParamsOptimization[[3]];
    ChangedArugsArray <- NULL;
    
    for (i in 2:length(new_ParamsOptimization)){
      for (j in 2:length(last_ParamsOptimization)){
        if (names(new_ParamsOptimization[i]) == names(last_ParamsOptimization[j]) & 
            new_ParamsOptimization[[i]] != last_ParamsOptimization[[j]]) {
          ChangedArugsArray <- c(ChangedArugsArray, names(new_ParamsOptimization[i]))
        } 
      }
    }
    
    newArugsNMs <- names(new_ParamsOptimization)[-1];
    lastArugsNMs <- names(last_ParamsOptimization)[-1];
    NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
    RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
    
    ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray)
    
    if(is.null(ChangedArugsArray) | length(ChangedArugsArray) == 0){
      # Not changed !
      plan@running.controller@others_1[1] <- FALSE;
    }
    
    if("param" %in% ChangedArugsArray){
      
      if(length(RmArugsArray == "param") != 0){
        if(RmArugsArray == "param"){
          # refer to: the omitted Argus is param
          if(last_ParamsOptimization$param == "SetPeakParam()"){
            # handle: the omitted param is SetPeakParam()
            plan@running.controller@others_1[1] <- FALSE;
          } else {
            # handle: the omitted param is not SetPeakParam()
            plan@running.controller@others_1[1] <- TRUE;
          }
        }
      }
      
      if (length(NewArugsArray == "param") != 0) {
        if (NewArugsArray == "param") {
          # refer to: the added Argus is param
          if (new_ParamsOptimization$param == "SetPeakParam()") {
            # handle: the added param is SetPeakParam()
            plan@running.controller@others_1[1] <- FALSE;
          } else {
            # handle: the added param is not SetPeakParam()
            plan@running.controller@others_1[1] <- TRUE;
          }
        }
      }
      
      if (length(NewArugsArray == "param") == 0 & 
          length(RmArugsArray == "param") == 0) {
        plan@running.controller@others_1[1] <- TRUE;
      }
      
    }
    
    if(plan@running.controller@operators[1]){
      # Data ROI changed ! Optimization has to be re-run !
      plan@running.controller@others_1[1] <- TRUE;
    }
    
    if(plan@running.controller@others_1[1]){
      # the peak profiling also need to be re-run
      plan@running.controller@operators[2] <- TRUE;
    }
  }
  
  # 3. ImportRawMSData: -----
  # 3.1 Note on controller: others c2, data import. c3, plot.
  
  new_ImportRawMSData <- new_command_set@ImportRawMSData[[3]];
  last_ImportRawMSData <- last_command_set@ImportRawMSData[[3]];
  ChangedArugsArray <- NULL;
  
  for (i in 2:length(new_ImportRawMSData)){
    for (j in 2:length(last_ImportRawMSData)){
      if (names(new_ImportRawMSData[i]) == names(last_ImportRawMSData[j]) & 
          new_ImportRawMSData[[i]] != last_ImportRawMSData[[j]]) {
        ChangedArugsArray <- c(ChangedArugsArray, names(new_ImportRawMSData[i]))
      } 
    }
  }
  
  newArugsNMs <- names(new_ImportRawMSData)[-1];
  lastArugsNMs <- names(last_ImportRawMSData)[-1];
  NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
  RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
  
  ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray);
  
  if(is.null(ChangedArugsArray) | length(ChangedArugsArray) == 0){
    # Not changed !
    plan@running.controller@data_import[c(1:2)] <- rep(TRUE, 2);
  }
  
  if("foldername" %in% ChangedArugsArray){
    plan@running.controller@data_import[1] <- TRUE;
  }
  
  if("plotSettings" %in% ChangedArugsArray){
    plan@running.controller@data_import[2] <- TRUE;
  }
  
  # 4. PeakProfiling: -----
  # 4.1 Note on controller peak_profiling: c1, peak picking; c2, peak alignment; c3, peak filing; c4, plotting.
  
  new_PeakProfiling <- new_command_set@PeakProfiling[[3]];
  last_PeakProfiling <- last_command_set@PeakProfiling[[3]];
  ChangedArugsArray <- NULL;
  
  for (i in 2:length(new_PeakProfiling)){
    for (j in 2:length(last_PeakProfiling)){
      if (names(new_PeakProfiling[i]) == names(last_PeakProfiling[j])) {
        if(new_PeakProfiling[[i]] != last_PeakProfiling[[j]]){
          ChangedArugsArray <- c(ChangedArugsArray, names(new_PeakProfiling[i]))
        }
      } 
    }
  }
  
  newArugsNMs <- names(new_PeakProfiling)[-1];
  lastArugsNMs <- names(last_PeakProfiling)[-1];
  NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
  RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
  
  ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray);
  
  if(is.null(ChangedArugsArray) | length(ChangedArugsArray) == 0){
    # Not changed !
    plan@running.controller@peak_profiling <- rep(FALSE, 4);
    names(plan@running.controller@peak_profiling) <- c("c1","c2","c3","c4");
  }
  
  if(plan@running.controller@operators[2]){
    # OPtimized parameters changed ! re-run everything in peak profiling.
    plan@running.controller@peak_profiling <- rep(TRUE, 4);
    names(plan@running.controller@peak_profiling) <- c("c1","c2","c3","c4");
    
    plan@running.controller@operators[3] <- TRUE;
  }
  
  if("plotSettings" %in% ChangedArugsArray){
    plan@running.controller@peak_profiling[4] <- TRUE;
  }
  
  if(.on.public.web){
    peakParams <- NULL;
    # load params.rda and envir.rds and do a comparison
    envir.path <- paste0(getwd(),"/temp/envir");
    envir_tmp <- readRDS(paste0(envir.path,"/envir.rds"));
    last_param <- envir_tmp[["param"]];
    load("params.rda");
    new_param <- peakParams;
    # TODO: need to configure with the web; handle the annotate params at the same time
    
  } else {
    
    if("Params" %in% ChangedArugsArray){
      
      newParamArgus <- new_PeakProfiling$Params;
      lastParamArgus <- last_PeakProfiling$Params;
      
      ChangedParamArgus <- ParamsChanged(lastParamArgus, newParamArgus);
      
      if(is.null(ChangedParamArgus)){
        plan@running.controller@peak_profiling[c(1:3)] <- rep(FALSE,3);
      } else if (any(ChangedParamArgus %in% c("min_peakwidth","max_peakwidth","mzdiff","ppm","noise","prefilter","value_of_prefilter",
                                              "Peak_method","snthresh","fwhm","sigma","steps"))){
        plan@running.controller@peak_profiling[c(1:3)] <- rep(TRUE,3);
        plan@running.controller@operators[3] <- TRUE;
      } else if (any(ChangedParamArgus %in% c("bw","RT_method","minFraction","minSamples","maxFeatures","family","smooth",
                                              "span","integrate","mzCenterFun","fitgauss"))){
        plan@running.controller@peak_profiling[c(1:3)] <- c(FALSE, TRUE, TRUE);
        plan@running.controller@operators[3] <- TRUE;
      } else if (ChangedParamArgus == "all"){
        plan@running.controller@peak_profiling[c(1:3)] <- rep(TRUE,3);
        plan@running.controller@operators[3] <- TRUE;
      }
    }
  }
  
  # 5. PeakAnnotation: -----
  # 5.1 Note on controller PeakAnnotation: c1, peak annotation; others not used for now.
  
  new_PeakAnnotation <- new_command_set@PeakAnnotation[[3]];
  last_PeakAnnotation <- last_command_set@PeakAnnotation[[3]];
  ChangedArugsArray <- NULL;
  
  for (i in 2:length(new_PeakAnnotation)){
    for (j in 2:length(last_PeakAnnotation)){
      if (names(new_PeakAnnotation[i]) == names(last_PeakAnnotation[j])) {
        if(new_PeakAnnotation[[i]] != last_PeakAnnotation[[j]]){
          ChangedArugsArray <- c(ChangedArugsArray, names(new_PeakAnnotation[i]))
        }
      } 
    }
  }
  
  newArugsNMs <- names(new_PeakAnnotation)[-1];
  lastArugsNMs <- names(last_PeakAnnotation)[-1];
  NewArugsArray <- setdiff(newArugsNMs,lastArugsNMs);
  RmArugsArray <- setdiff(lastArugsNMs, newArugsNMs);
  
  ChangedArugsArray <- c(ChangedArugsArray, NewArugsArray, RmArugsArray);
  
  if(.on.public.web){
    # load params.rda and envir.rds and do a comparison
    # TODO: need to configure with the web and profiling
    
  } else {
    
    if ("annotaParam" %in% ChangedArugsArray) {
      newAnnoParamArgus <- new_PeakAnnotation$annotaParam;
      lastAnnoParamArgus <- last_PeakAnnotation$annotaParam;
      
      ChangedParamArgus <-
        ParamsChanged(lastAnnoParamArgus, newAnnoParamArgus);
      
      if (is.null(ChangedParamArgus)) {
        plan@running.controller@peak_annotation[1] <- FALSE;
      } else {
        plan@running.controller@peak_annotation[1] <- TRUE;
      }
    } else {
      plan@running.controller@peak_annotation[1] <- FALSE;
    }
  }
  
  if(plan@running.controller@operators[3]){
    plan@running.controller@peak_annotation[1] <- TRUE;
  }
  
  # 6. FormatPeakList: -----
  # 6.1 Note on controller FormatPeakList.
  # No need to resume for this function
  
  return(plan)
}

#' @noRd
perform.plan <- function(plan.set){
  
  if(class(plan.set) == "OptiCommandSet"){
    perform.command(plan.set@ROIExtraction);
    perform.command(plan.set@ParamsOptimization);
    perform.command(plan.set@ImportRawMSData);
    perform.command(plan.set@PeakProfiling);
    perform.command(plan.set@PeakAnnotation);
    perform.command(plan.set@FormatPeakList);
  }
  
  if(class(plan.set) == "CustCommandSet"){
    perform.command(plan.set@ImportRawMSData);
    perform.command(plan.set@PeakProfiling);
    perform.command(plan.set@PeakAnnotation);
    perform.command(plan.set@FormatPeakList);
  }
}

#' @noRd
perform.command <- function(command){
  
  if(!exists(".SwapEnv")){
    .SwapEnv <<- new.env(parent = .GlobalEnv);
    .SwapEnv$.optimize_switch <- FALSE;
    .SwapEnv$count_current_sample <- 0;
    .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
    .SwapEnv$envir <- new.env();
  }
  
  envir <- .SwapEnv$envir;
  
  if (as.character(command[[1]])=="<-"){
    
    eval(command,envir = envir);
    
    envir.path <- paste0(getwd(),"/temp/envir");
    saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
  } else {
    
    eval(command,envir = envir);
    
    envir.path <- paste0(getwd(),"/temp/envir");
    saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
  }
  
}

#' @noRd
cache.save <- function(obj, funpartnm){
  
  tmp_path <- paste0(getwd(),"/temp/cache");
  if (!dir.exists(tmp_path)){dir.create(tmp_path,recursive = T)};
  temp <- tempfile(tmpdir=tmp_path,fileext = ".rds");
  saveRDS(obj,file = temp);
  
  tmp_path_r <- paste0(getwd(),"/temp/records");
  if (!dir.exists(tmp_path_r)){dir.create(tmp_path,recursive = T)};
  if (!file.exists(paste0(getwd(),"/temp/records/file_record.rds"))){
    
  };
  temp_file_name <- basename(temp)
  info.save(funpartnm, tmp_path_r, temp_file_name)
}

#' @noRd
info.save <- function(funpartnm, tmp_path_r, temp_file_name){
  
  info.matrix <- readRDS(paste0(tmp_path_r,"/records.rds"));
  if (identical(which(info.matrix[,1] == funpartnm),integer(0))){
    record.row <- which(is.na(info.matrix[,1]))[1];
  } else {
    record.row <- which(info.matrix[,1]==funpartnm)[1];
  }
  
  info.matrix[record.row,1] <- funpartnm;
  info.matrix[record.row,2] <- temp_file_name;
  
  saveRDS(info.matrix,file = paste0(tmp_path_r,"/records.rds"))
}

#' @noRd
cache.read <- function(function.name, point){
  
  info.matrix <- readRDS(paste0(getwd(),"/temp/records/records.rds"));
  temp_point <- paste0(function.name,"_",point);
  temp_file <- info.matrix[match(temp_point,info.matrix[,1]),2];
  
  obj <- readRDS(paste0(getwd(),"/temp/cache/",temp_file));
  
  return(obj)
  
}

#' @noRd
profiling_param_identifier <- function(new_command,last_command){
  
  if(!exists(".SwapEnv")){
    .SwapEnv <<- new.env(parent = .GlobalEnv);
    .SwapEnv$.optimize_switch <- FALSE;
    .SwapEnv$count_current_sample <- 0;
    .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
    .SwapEnv$envir <- new.env();
  }
  
  envir <- .SwapEnv$envir
  
  new <- eval(new_command,envir = envir);
  last <- eval(last_command,envir = envir);
  diff.names <- character();
  
  for (i in 1:length(new)){
    for (j in 1:length(last)){
      
      if ((names(new[i]) == names(last[j]))){
        if (unlist(unname(new[i])) != unlist(unname(last[j]))) {
          diff.names <- c(diff.names,names(new[i]))
        }
      }
    }
  };
  return(diff.names)
  
}

#' @noRd
CommandsVerify <- function(commands){
  
  CustPipeline <- OptiPipeline <- FALSE;

  # Identify if the commandset is automated optimization
  OptiPipeline <- any(unlist(lapply(
    commands,
    FUN = function(x) {
      if (class(x[[3]]) == "call") {
        if (x[[3]][[1]] == "PerformParamsOptimization") {
          return(TRUE)
        }
      }
    }
  )))
  
  CustPipeline <- all(unlist(lapply(
    commands,
    FUN = function(x) {
      if (class(x[[3]]) == "call") {
        if (x[[3]][[1]] == "PerformParamsOptimization") {
          return(FALSE)
        } else {
          return(TRUE)
        }
      }
    }
  )))
  
  # Organize the commands as a standard opti/cust commandset for next step
  # including folowing parts:
  
  # 1. PerformROIExtraction / PerformDataTrimming;
  # 2. PerformParamsOptimization;
  # 3. ImportRawMSData;
  # 4. PerformPeakProfiling;
  # 5. PerformPeakAnnotation;
  # 6. FormatPeakList;
  
  if(OptiPipeline){
    
    VARCommandArray <- NULL;
    FUNCommandArray <- NULL;
    StandCommand <- new("OptiCommandSet");
    
    for (i in seq_along(commands)){

      if(class(commands[[i]][[3]]) == "call"){
        
        if(commands[[i]][[3]][[1]] == "PerformDataTrimming" | commands[[i]][[3]][[1]] == "PerformROIExtraction"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("datapath", "mode", "write", "mz", "mzdiff", "rt", "rtdiff", "rt.idx", "plot", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@ROIExtraction <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformParamsOptimization"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "param", "method", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@ParamsOptimization <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "ImportRawMSData"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "foldername", "mode", "ncores", "plotSettings", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@ImportRawMSData <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformPeakProfiling"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "Params", "plotSettings", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@PeakProfiling <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformPeakAnnotation"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "annotaParam", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@PeakAnnotation <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "FormatPeakList"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "annParams", "filtIso", "filtAdducts", "missPercent");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@FormatPeakList <- tmp_command;
        }
        
      }
      
    }
    
    if(length(commands) > 6){
      MessageOutput("More functions than standard OptiPipeline were detected !\n")
      warning(paste0("NOTE: Only ",paste(scales::ordinal(FUNCommandArray),collapse = ", "), 
                           " functions in this plan and their direct defination on the argument will be included !\n"))
    }
  }
  
  if(CustPipeline){
    
    VARCommandArray <- NULL;
    FUNCommandArray <- NULL;
    StandCommand <- new("CustCommandSet");
    
    for (i in seq_along(commands)){
      
      if(class(commands[[i]][[3]]) == "call"){
        
        if(commands[[i]][[3]][[1]] == "ImportRawMSData"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "foldername", "mode", "ncores", "plotSettings", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@ImportRawMSData <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformPeakProfiling"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "Params", "plotSettings", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@PeakProfiling <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformPeakAnnotation"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "annotaParam", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@PeakAnnotation <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "FormatPeakList"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "annParams", "filtIso", "filtAdducts", "missPercent");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- which(ArgNM == ArguList);
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@FormatPeakList <- tmp_command;
        }
        
      }
      
    }
    
    if(length(commands) > 4){
      MessageOutput("More functions than standard OptiPipeline were detected !\n")
      MessageOutput(paste0("NOTE: Only ",paste(scales::ordinal(FUNCommandArray),collapse = ", "), 
                   " functions in this plan and their direct defination on the argument will be included !\n"))
    }
    
  }
  
  return(StandCommand)
}

#' @noRd
CommandOrganize <- function(command, commands, FUNPos){
  
  varNMs <- as.character(command[[3]]);
  if(varNMs[1] == "PerformROIExtraction"){
    Vars <- varNMs[-1];
  } 
  
  if(varNMs[2] != "mSet"){
    Vars <- varNMs[-1];
  } else {
    Vars <- varNMs[c(-1,-2)];
  }
  
  for (i in Vars){
    
    VARDefinedPos <- NULL;
    VarPos <- which(varNMs == i);
    
    for (j in seq_len(FUNPos)){
      if (commands[[j]][[2]] == i){
        VARDefinedPos <- c(VARDefinedPos, j);
      }
    }
    
    if (is.null(VARDefinedPos)) {
      next()
    } else if(length(VARDefinedPos) == 1) {
      command[[3]][[VarPos]]<- commands[[VARDefinedPos]][[3]];
    } else if(length(VARDefinedPos) > 1) {
      VARDefinedPos <- VARDefinedPos[which.min(abs(VARDefinedPos - FUNPos))];
      command[[3]][[VarPos]]<- commands[[VARDefinedPos]][[3]];
    }
  }
  
  return(command)
}

#' #' @noRd
#' CreateRawRscript <- function(guestName, planString, planString2, rawfilenms.vec){
#'   
#'   guestName <<- guestName;
#'   planString <<- planString;
#'   
#'   if(dir.exists("/home/glassfish/payara5/glassfish/domains/")){
#'     
#'     scr.path <- "/home/glassfish/payara5/glassfish/domains/domain1/applications/MetaboAnalyst/resources/rscripts/metaboanalystr/"
#'     users.path <-paste0("/data/glassfish/projects/metaboanalyst//", guestName)
#'     
#'   } else if(dir.exists("/media/zzggyy/disk/")){
#'     
#'     users.path <-getwd();
#'     scr.path <-"/media/zzggyy/disk/MetaboAnalyst/target/MetaboAnalyst-5-18/resources/rscripts/metaboanalystr/"
#'     
#'   }else {
#'     
#'     users.path <-getwd();
#'     scr.path <-paste0(strsplit(users.path, "users")[[1]][1], "rscripts/metaboanalystr/")
#'     
#'   }
#'   
#'   ## Prepare the script for running
#'   str <- paste0('scripts.path <- ','"', scr.path,'"');
#'   str <- paste0(str, '\n', 'general_files <- c("general_data_utils.R","general_misc_utils.R","general_lib_utils.R","generic_c_utils.R");')
#'   str <- paste0(str, '\n', 'spectra_files <- c("spectra_generic_utils.R","data_trimming.R","parameters_db.R","parameters_optimization.R","preproc_utils.R","spectra_resume.R","spectra_utils.R");')
#'   str <- paste0(str, '\n', 'rawfilenms <<-', rawfilenms.vec)
#'   str <- paste0(str, '\n', 'err.vec <<-', '""')
#'   str <- paste0(str, '\n', '.on.public.web <<- TRUE')
#'   str <- paste0(str, '\n', 'file.sources = c(general_files, spectra_files)')
#'   str <- paste0(str, '\n', 'for (f in file.sources) {source(paste0(scripts.path, f))}')
#'   
#'   ## Construct the opt pipeline
#'   if(planString2 == "opt"){
#'     str <- paste0(str, '\n', 'plan <- InitializaPlan("raw_opt","' , users.path,  '/")')
#'     str <- paste0(str, '\n', 'data_folder_QC <- "',users.path , '/upload/QC/"')
#'     str <- paste0(str, '\n',  planString)
#'     str <- paste0(str, '\n',  "ExecutePlan(plan)")
#'   }
#'   
#'   ## Construct the default pipeline
#'   if(planString2 == "default"){
#'     str <- paste0(str, '\n', 'plan <- InitializaPlan("raw_ms","' , users.path,  '/")')
#'     str <- paste0(str, '\n',  planString)
#'     str <- paste0(str, '\n',  "ExecutePlan(plan)")
#'   }
#'   
#'   # sink command for running
#'   sink("ExecuteRawSpec.R");
#'   cat(str);
#'   sink();
#'   return(1)
#'   
#' }
#' 
#' #' @noRd
#' CreateRawRscript2 <- function(guestName, meth){
#'   
#'   guestName <<- guestName;  
#'   
#'   if(dir.exists("/home/glassfish/payara5/glassfish/domains/")){
#'     
#'     scr.path <- "/home/glassfish/payara5/glassfish/domains/domain1/applications/MetaboAnalyst/resources/rscripts/metaboanalystr/"
#'     users.path <-paste0("/data/glassfish/projects/metaboanalyst/", guestName)
#'     
#'   } else {
#'     
#'     users.path <-getwd();
#'     scr.path <-paste0(strsplit(users.path, "users")[[1]][1], "rscripts/metaboanalystr/")
#'     
#'   }
#'   
#'   str <- paste0('scripts.path <- ','"', scr.path,'"');
#'   str <- paste0(str, '\n', 'general_files <- c("general_data_utils.R","general_misc_utils.R","general_lib_utils.R","generic_c_utils.R");')
#'   str <- paste0(str, '\n', 'spectra_files <- c("spectra_generic_utils.R","data_trimming.R","parameters_db.R","parameters_optimization.R","preproc_utils.R","spectra_resume.R","spectra_utils.R");')
#'   str <- paste0(str, '\n', '.on.public.web <<- TRUE')
#'   str <- paste0(str, '\n', 'file.sources = c(general_files, spectra_files)')
#'   str <- paste0(str, '\n', 'for (f in file.sources) {source(paste0(scripts.path, f))}')
#'   
#'   
#'   ## Prepare the script for running
#'   if(meth == "auto"){
#'     str <- paste0(str, '\n', 'FastRunningShow_auto ("',users.path,'")');
#'   } else {
#'     str <- paste0(str, '\n', 'FastRunningShow_customized ("',users.path,'")');
#'   }
#'   
#'   # sink command for running
#'   sink("ExecuteRawSpec.R");
#'   cat(str);
#'   sink();
#'   return(1)
#'   
#' }

#' @noRd
ParamsChanged <- function(lastParamArgus, newParamArgus){
  
  if (is.call(newParamArgus) & is.call(lastParamArgus)){
    new_param <- eval(newParamArgus);
    last_param <- eval(lastParamArgus);
    
    identifiers <- NULL;
    
    for (i in 1:length(new_param)){
      for (j in 1:length(last_param)){
        if (names(new_param[i])==names(last_param[j]) & new_param[[i]] != last_param[[j]]){
          identifiers <- c(identifiers, names(new_param[i]))
        }
      }
    }
    
    return(identifiers)
    
  } else {
    return("all")
  }
  
}


## End of the resuming script ---
