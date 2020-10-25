#' @title Initializing running plan
#' @description Initialize a running plan
#' @param type Character, Initialized plan type for a resumable running mode. Can be "raw_opt" for 
#' automated optimization option, or "raw_ms" for customized pipeline.
#' @export
#' @seealso \code{\link{ExecutePlan}} for the this resumable running pipeline.
#' @examples 
#' library(OptiLCMS);
#' plan <- InitializaPlan("raw_opt")
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)

InitializaPlan <- function(type = "raw_ms"){
  
  if(.on.public.web){
    tmp_cache_path <- getwd();
  }else {
    tmp_cache_path <- tempdir();
  }
  
  plan <- NULL;
  
  if(!exists(".SwapEnv")){
    .SwapEnv <<- new.env(parent = .GlobalEnv);
    .SwapEnv$.optimize_switch <- FALSE;
    .SwapEnv$count_current_sample <- 0;
    .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
    .SwapEnv$envir <- new.env();
    .SwapEnv$PlanWorkingDir <- paste0(tmp_cache_path, "/specTemp/");
  }
  
  temp.path <- paste0(tmp_cache_path, "/specTemp/");
  
  if(dir.exists(temp.path) & !.on.public.web){
    unlink(temp.path, recursive = TRUE, force = TRUE);
  }
  
  .SwapEnv$PlanWorkingDir <- temp.path;

  #=============== raw_opt
  if (type == "raw_opt") {
    
    plan.path <- paste0(temp.path, "/plan");
    
    if (!dir.exists(plan.path)) { # THIS MEANS: if no following stuff, will create them (no matter local or web); If local, always recreate them.
      dir.create(paste0(temp.path, "/plan"),
                 recursive = TRUE, showWarnings = FALSE);
      
      #.running.as.plan <- .SwapEnv$.running.as.plan <- TRUE;
      plan <- new("ResumingePlan");
      .plan_count <- plan@PlanNumber <- 0;
      plan@WorkingDir <- temp.path;
      
      saveRDS(plan, file = paste0(plan.path, "/plan.rds"));
      saveRDS(.plan_count, file = paste0(plan.path, "/plan_count.rds"));
      #----------------------/
      .optimize_switch <- .SwapEnv$.optimize_switch <- TRUE;
      switch.path <- paste0(temp.path, "/plan");
      saveRDS(.optimize_switch,
              file = paste0(switch.path, "/optimize_switch_", .plan_count, ".rds"))
      
    } else {
      if(file.exists(paste0(plan.path, "/plan.rds"))){
        plan <- readRDS(paste0(plan.path, "/plan.rds"));
        plan@WorkingDir <- temp.path;
      } else {
        stop(paste0("Cache file has been damaged ! Please remove your cache folder from ", tmp_cache_path, "!"))
      }
    }

    #----------------------
    envir.path <- paste0(temp.path, "/envir")
    if (!dir.exists(envir.path)) {
      dir.create(paste0(temp.path, "/envir"),
                 recursive = TRUE, showWarnings = FALSE)
      
      envir <- .SwapEnv$envir <- new.env()
      saveRDS(envir, file = paste0(envir.path, "/envir.rds"))
    }
  }
  
  #================== raw_ms
  if (type == "raw_ms") {
    plan.path <- paste0(temp.path, "/plan")
    
    if (!dir.exists(plan.path)) {
      dir.create(paste0(temp.path, "/plan"),
                 recursive = TRUE, showWarnings = FALSE)
      
      MessageOutput(paste0("Running Status -- Plan Initialized Successfully at: ", 
                           Sys.time(), 
                           "\nCurrent OptiLCMS version is ",
                           packageVersion("OptiLCMS"),
                           "\nPlease define your running plan ..."), "\n", 0);
      #---------------
      #.running.as.plan <- .SwapEnv$.running.as.plan <- TRUE;
      plan <- new("ResumingePlan");
      .plan_count <- plan@PlanNumber <- 0;
      plan@WorkingDir <- temp.path;
      
      saveRDS(plan, file = paste0(plan.path, "/plan.rds"));
      saveRDS(.plan_count, file = paste0(plan.path, "/plan_count.rds"));
      #----------------------
      .optimize_switch <- .SwapEnv$.optimize_switch <- FALSE;
      switch.path <- paste0(temp.path, "/plan");
      saveRDS(.optimize_switch,
              file = paste0(switch.path, "/optimize_switch_", .plan_count, ".rds"));
    } else {
      if(file.exists(paste0(plan.path, "/plan.rds"))){
        plan <- readRDS(paste0(plan.path, "/plan.rds"));
        plan@WorkingDir <- temp.path;
      } else {
        stop(paste0("Cache file has been damaged ! Please remove your cache folder from ", tmp_cache_path, "!"))
      }
    }
    
    #----------------------
    envir.path <- paste0(temp.path, "/envir")
    if (!dir.exists(envir.path)) {
      dir.create(paste0(temp.path, "/envir"),
                 recursive = TRUE, showWarnings = FALSE)
      
      envir <- .SwapEnv$envir <- new.env()
      saveRDS(envir, file = paste0(envir.path, "/envir.rds"))
    }

  }

  #============ Recording cache file information
  record.path <- paste0(temp.path, "/records");
  
  if (!dir.exists(record.path)) {
    record.info <- matrix(nrow = 30, ncol = 2);
    
    dir.create(paste0(temp.path, "/records"),
               recursive = TRUE, showWarnings = FALSE);
    
    saveRDS(record.info, file = paste0(record.path, "/records.rds"))
  }

  return(plan)
}

#' @title running.plan
#' @description define a plan for resumalbe running
#' @param plan ResummingPlan object. The object is generated by 'InitializaPlan' function.
#' @param ... Multiple Processing commands can be input here.
#' @export
#' @seealso \code{\link{ExecutePlan}} for the this resumable running pipeline.
#' @examples 
#' ##' Download the raw spectra data
#' DataFiles <- dir(system.file("mzData", package = "mtbls2"), full.names = TRUE,
#'             recursive = TRUE)[c(10:12, 14:16)]
#' ##' Create a phenodata data.frame
#' pd <- data.frame(sample_name = sub(basename(DataFiles), pattern = ".mzData", 
#'                                    replacement = "", fixed = TRUE),
#'                  sample_group = c(rep("col0", 3), rep("cyp79", 3)),
#'                  stringsAsFactors = FALSE)
#' 
#' ##' Load OptiLCMS
#' library(OptiLCMS);
#' 
#' ##' Initialize your plan
#' plan <- InitializaPlan("raw_opt")
#' 
#' ##' Define your plan
#' plan <- running.plan(plan,
#'                      mSet <- PerformROIExtraction(datapath = DataFiles[c(1:2)], rt.idx = 0.025,
#'                                                   plot = FALSE, rmConts = FALSE,
#'                                                   running.controller = rc),
#'                      param_initial <- SetPeakParam(),
#'                      best_parameters <- PerformParamsOptimization(mSet = mSet, param_initial,
#'                                                                   ncore = 1,
#'                                                                   running.controller = rc),
#'                      param <- best_parameters,
#'                      plotSettings1 <- SetPlotParam(Plot=TRUE),
#'                      plotSettings2 <- SetPlotParam(Plot=TRUE),
#'                      mSet <- ImportRawMSData(mSet = mSet, path = DataFiles, 
#'                                              metadata = pd,
#'                                              plotSettings = plotSettings1,
#'                                              running.controller = rc),
#'                      mSet <- PerformPeakProfiling(mSet = mSet, Params = param,
#'                                                   plotSettings = plotSettings2, ncore = 1,
#'                                                   running.controller = rc),
#'                      annParams <- SetAnnotationParam(polarity = 'negative',
#'                                                      mz_abs_add = 0.025),
#'                      mSet <- PerformPeakAnnotation(mSet = mSet,
#'                                                    annotaParam = annParams, ncore =1,
#'                                                    running.controller = rc),
#'                      mSet <- FormatPeakList(mSet = mSet, annParams,
#'                                             filtIso =FALSE, filtAdducts = FALSE,
#'                                             missPercent = 1));
#' ##' Run it!
#' # result <- ExecutePlan(plan);
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)

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
  
  plan.path <- paste0(plan@WorkingDir, "/plan");
  saveRDS(plan, file = paste0(plan.path, "/plan.rds"));
  saveRDS(.plan_count, file = paste0(plan.path, "/plan_count.rds"));
  
  recordMarker_resetter(.plan_count);
  
  return(plan)
}

#' .get.current.plan
#' @noRd
.get.current.plan <- function(plan){
  
  if (is.null(plan)){
    return(list())
  } else {
    return (plan)
  }
  
}

#' @title ExecutePlan
#' @param plan ResummingPlan object. The object is generated by running.plan() function.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
#' @export
#' @examples 
#' ##' Download the raw spectra data
#' DataFiles <- dir(system.file("mzData", package = "mtbls2"), full.names = TRUE,
#'             recursive = TRUE)[c(10:12, 14:16)]
#' ##' Create a phenodata data.frame
#' pd <- data.frame(sample_name = sub(basename(DataFiles), pattern = ".mzData", 
#'                                    replacement = "", fixed = TRUE),
#'                  sample_group = c(rep("col0", 3), rep("cyp79", 3)),
#'                  stringsAsFactors = FALSE)
#' 
#' ##' Load OptiLCMS
#' library(OptiLCMS);
#' 
#' ##' Initialize your plan
#' plan <- InitializaPlan("raw_opt")
#' 
#' ##' Define your plan
#' plan <- running.plan(plan,
#'                      mSet <- PerformROIExtraction(datapath = DataFiles[c(1:2)], rt.idx = 0.025,
#'                                                   plot = FALSE, rmConts = FALSE,
#'                                                   running.controller = rc),
#'                      param_initial <- SetPeakParam(),
#'                      best_parameters <- PerformParamsOptimization(mSet = mSet, param_initial,
#'                                                                   ncore = 1,
#'                                                                   running.controller = rc),
#'                      param <- best_parameters,
#'                      plotSettings1 <- SetPlotParam(Plot=TRUE),
#'                      plotSettings2 <- SetPlotParam(Plot=TRUE),
#'                      mSet <- ImportRawMSData(mSet = mSet, path = DataFiles, 
#'                                              metadata = pd,
#'                                              plotSettings = plotSettings1,
#'                                              running.controller = rc),
#'                      mSet <- PerformPeakProfiling(mSet = mSet, Params = param,
#'                                                   plotSettings = plotSettings2, ncore = 1,
#'                                                   running.controller = rc),
#'                      annParams <- SetAnnotationParam(polarity = 'negative',
#'                                                      mz_abs_add = 0.025),
#'                      mSet <- PerformPeakAnnotation(mSet = mSet,
#'                                                    annotaParam = annParams, ncore =1,
#'                                                    running.controller = rc),
#'                      mSet <- FormatPeakList(mSet = mSet, annParams,
#'                                             filtIso =FALSE, filtAdducts = FALSE,
#'                                             missPercent = 1));
#' ##' Run it!
#' # result <- ExecutePlan(plan);
#' 
#' ##' Re-define your plan with a change on mz_abs_add from 0.025 to 0.035
#' plan <- running.plan(plan,
#'                      mSet <- PerformROIExtraction(datapath = DataFiles[c(1:2)], rt.idx = 0.025,
#'                                                   plot = FALSE, rmConts = FALSE,
#'                                                   running.controller = rc),
#'                      param_initial <- SetPeakParam(),
#'                      best_parameters <- PerformParamsOptimization(mSet = mSet, param_initial,
#'                                                                   ncore = 1,
#'                                                                   running.controller = rc),
#'                      param <- best_parameters,
#'                      plotSettings1 <- SetPlotParam(Plot=TRUE),
#'                      plotSettings2 <- SetPlotParam(Plot=TRUE),
#'                      mSet <- ImportRawMSData(mSet = mSet, 
#'                                              path = DataFiles, 
#'                                              metadata = pd,
#'                                              plotSettings = plotSettings1,
#'                                              running.controller = rc),
#'                      mSet <- PerformPeakProfiling(mSet = mSet, Params = param,
#'                                                   plotSettings = plotSettings2, ncore = 1,
#'                                                   running.controller = rc),
#'                      annParams <- SetAnnotationParam(polarity = 'negative',
#'                                                      mz_abs_add = 0.035),
#'                      mSet <- PerformPeakAnnotation(mSet = mSet,
#'                                                    annotaParam = annParams, ncore =1,
#'                                                    running.controller = rc),
#'                      mSet <- FormatPeakList(mSet = mSet, annParams,
#'                                             filtIso =FALSE, filtAdducts = FALSE,
#'                                             missPercent = 1));
#' 
#' ##' Re-run it! Most steps will be resumed from cache and save your time!
#' # result <- ExecutePlan(plan);

ExecutePlan <- function(plan=NULL){
  
  MessageOutput("This plan is being excuted !\n", SuppressWeb = TRUE);
  MessageOutput(paste0("Working Directory: ",getwd()), SuppressWeb = TRUE);
  MessageOutput(paste0("Cache Directory: ",plan@WorkingDir,"\n"), SuppressWeb = TRUE);
  
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
    envir.path <- paste0(plan@WorkingDir,"/envir");
    envir.save(envir, path = envir.path);
    
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
    save(plan, file = "plan_1.rda")
    ## Module 7 - Detect whether some steps have been run in last plan excuting process (Secondary Decision switch)
    plan <- recording_identifier(plan);
    save(plan, file = "plan_2.rda")
    ## Module 8 - Dectect whether current FilesInclusion is different from the last one (Final Decision switch)
    plan <- FilesInclusion_identifier(plan);
    save(plan, file = "plan_3.rda")
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
      if(.on.public.web){
        print_mes <- paste0("<font color=\"red\">","\nERROR:",mSetInfo$message,"</font>");
        write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
      } else {
        stop(paste0("EXCEPTION POINT CODE: ", mSetInfo$message));
      }
    }
    
    envir <- .SwapEnv$envir;
    envir.path <- paste0(plan@WorkingDir,"/envir");
    envir.save(envir, path = envir.path);
    
  } else {
    stop("Wrong plan has been prepared !");
  }
  
  return(envir)
}

#' @noRd
recording_identifier <- function(plan) {
  
  .plan_count <- plan@PlanNumber;
  record.path <- paste0(plan@WorkingDir,"/records");
  
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
FilesInclusion_identifier <- function(plan) {
  
  if(.on.public.web){
    
    envOld <- mSet <- NULL;
    
    load("mSet.rda");
    currentfiles <- mSet@rawfiles;
    envOld <- readRDS("specTemp/envir/envir.rds");
    lastfiles <- envOld$mSet@rawfiles;
    
  } else {
    
    if(file.exists(paste0(tempdir(), "/envir/envir.rds"))){
      
      envOld <- readRDS(paste0(tempdir(), "/envir/envir.rds"));
      lastfiles <- envOld$mSet@rawfiles;
      
      currentfiles <- list.files(plan@CommandSet[[plan@PlanNumber]]@ImportRawMSData[[3]][["path"]], 
                                 recursive = TRUE, full.names = TRUE);
      
    } else {
      stop("Cache missing! Please redo the InitializePlan() from the begining!")
    }
  }
  
  currentfiles <- basename(currentfiles);
  lastfiles <- basename(lastfiles);
  
  res <- c(setdiff(currentfiles, lastfiles), setdiff(lastfiles, currentfiles));

  if(length(res) != 0){
    plan@running.controller <- controller.resetter();
  }
  
  return(plan)
}

#' @noRd
recordMarker_resetter <- function(.plan_count){
  
  tmp_cache_path <- tempdir();
  
  if(!exists(".SwapEnv")){
    .SwapEnv <<- new.env(parent = .GlobalEnv);
    .SwapEnv$.optimize_switch <- FALSE;
    .SwapEnv$count_current_sample <- 0;
    .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
    .SwapEnv$envir <- new.env();
    .SwapEnv$PlanWorkingDir <- paste0(tmp_cache_path, "/specTemp/");
  }
  
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
  record.path <- paste0(.SwapEnv$PlanWorkingDir,"/records");
  saveRDS(record.marker,file = paste0(record.path,"/records_marker_",.plan_count,".rds"));
}

#' @noRd
marker_record <- function(functionNM){
  
  tmp_cache_path <- tempdir();
  
  if(!exists(".SwapEnv")){
    .SwapEnv <<- new.env(parent = .GlobalEnv);
    .SwapEnv$.optimize_switch <- FALSE;
    .SwapEnv$count_current_sample <- 0;
    .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
    .SwapEnv$envir <- new.env();
    .SwapEnv$PlanWorkingDir <- paste0(tmp_cache_path, "/specTemp/");
  }
  
  record.path <- paste0(.SwapEnv$PlanWorkingDir,"/records");
  plan.path <- paste0(.SwapEnv$PlanWorkingDir,"/plan/");
  .plan_count <- readRDS(paste0(plan.path,"plan_count.rds"));
  
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

#' @noRd
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
    plan@running.controller@data_import[c(1:2)] <- rep(FALSE, 2);
  }
  
  if("path" %in% ChangedArugsArray){
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
    envir.path <- paste0(plan@WorkingDir,"/envir");
    envir_tmp <- readRDS(paste0(envir.path,"/envir.rds"));
    last_param <- envir_tmp$mSet@params;
    load("params.rda");
    new_param <- updateRawSpectraParam(peakParams);
    
    # Compare the difference
    if(class(new_command_set) == "OptiCommandSet"){
      # For auto web pipeline
      if(plan@running.controller@operators[2]){
        plan@running.controller@peak_profiling[c(1:3)] <- rep(TRUE,3);
        plan@running.controller@operators[3] <- TRUE;
      }
      
    } else {
      # For non-optimized web pipeline
      ChangedParamArgus <- names(new_param)[
      which(unlist(
        sapply(c(1:25), function(i){
          return(last_param[[i]] != new_param[[i]])
        })
      ))]
     
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
      }
    }
    
    
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
    
    envir.path <- paste0(.SwapEnv$PlanWorkingDir,"/envir");
    envir.save(envir, path = envir.path);
  } else {
    
    eval(command,envir = envir);
    
    envir.path <- paste0(.SwapEnv$PlanWorkingDir,"/envir");
    envir.save(envir, path = envir.path);
  }
  
}

#' @noRd
envir.save <- function(envir, path){
  # this function is used to save the envir cache slowly 
  # and avoid error from reading due to the interupt last time
  
  saveRDS(envir, file = paste0(path,"/envir_tmp.rds"));
  
  tmmenvir <- try(readRDS(paste0(path,"/envir_tmp.rds")));
  
  # To ensure the envir cache can be read
  if(class(tmmenvir) != "try-error"){
    file.rename(paste0(path,"/envir_tmp.rds"), paste0(path,"/envir.rds"))
  } else {
    # do nothing
  }
}

#' @noRd
cache.save <- function(obj, funpartnm){
  
  tmp_cache_path <- tempdir();
  
  if(!exists(".SwapEnv")){
    .SwapEnv <<- new.env(parent = .GlobalEnv);
    .SwapEnv$.optimize_switch <- FALSE;
    .SwapEnv$count_current_sample <- 0;
    .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
    .SwapEnv$envir <- new.env();
    .SwapEnv$PlanWorkingDir <- paste0(tmp_cache_path, "/specTemp/");
  }
  
  tmp_path <- paste0(.SwapEnv$PlanWorkingDir,"/cache");
  if (!dir.exists(tmp_path)){dir.create(tmp_path,recursive = TRUE)};
  temp <- tempfile(tmpdir=tmp_path,fileext = ".rds");
  saveRDS(obj,file = temp);
  
  tmp_path_r <- paste0(.SwapEnv$PlanWorkingDir,"/records");
  if (!dir.exists(tmp_path_r)){dir.create(tmp_path,recursive = TRUE)};
  if (!file.exists(paste0(.SwapEnv$PlanWorkingDir,"/records/file_record.rds"))){
    
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
  
  tmp_cache_path <- tempdir();
  
  if(!exists(".SwapEnv")){
    .SwapEnv <<- new.env(parent = .GlobalEnv);
    .SwapEnv$.optimize_switch <- FALSE;
    .SwapEnv$count_current_sample <- 0;
    .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
    .SwapEnv$envir <- new.env();
    .SwapEnv$PlanWorkingDir <- paste0(tmp_cache_path, "/specTemp/");
  }
  
  info.matrix <- readRDS(paste0(.SwapEnv$PlanWorkingDir,"/records/records.rds"));
  temp_point <- paste0(function.name,"_",point);
  temp_file <- info.matrix[match(temp_point,info.matrix[,1]),2];
  
  obj <- readRDS(paste0(.SwapEnv$PlanWorkingDir,"/cache/",temp_file));
  
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
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@ROIExtraction <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformParamsOptimization"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "param", "method", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@ParamsOptimization <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "ImportRawMSData"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "foldername", "mode", "ncores", "plotSettings", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@ImportRawMSData <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformPeakProfiling"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "Params", "plotSettings", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@PeakProfiling <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformPeakAnnotation"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "annotaParam", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@PeakAnnotation <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "FormatPeakList"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "annParams", "filtIso", "filtAdducts", "missPercent");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@FormatPeakList <- tmp_command;
        }
        
      }
      
    }
    
    if(length(commands) > 6){
      MessageOutput("More functions than standard OptiPipeline were detected !\n", SuppressWeb = TRUE)
      MessageOutput(paste0("NOTE: Only ",paste(scales::ordinal(FUNCommandArray),collapse = ", "), 
                     " functions in this plan and their direct defination on the argument will be included !\n"), SuppressWeb = TRUE)
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
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@ImportRawMSData <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformPeakProfiling"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "Params", "plotSettings", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@PeakProfiling <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "PerformPeakAnnotation"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "annotaParam", "ncore", "running.controller");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@PeakAnnotation <- tmp_command;
        }
        
        if(commands[[i]][[3]][[1]] == "FormatPeakList"){
          FUNCommandArray <- c(FUNCommandArray, i);
          tmp_command <- CommandOrganize(commands[[i]], commands, i);
          
          ArguList <- c("mSet", "annParams", "filtIso", "filtAdducts", "missPercent");
          ArgNM <- ArguList[which(names(tmp_command[[3]])[-1] == "")];
          ArgPos <- unlist(sapply(ArgNM, function(x){which(x == ArguList)}))
          names(tmp_command[[3]])[ArgPos + 1] <- ArgNM;
          
          StandCommand@FormatPeakList <- tmp_command;
        }
        
      }
      
    }
    
    if(length(commands) > 4){
      MessageOutput("More functions than standard OptiPipeline were detected !\n", SuppressWeb = TRUE)
      MessageOutput(paste0("NOTE: Only ",paste(scales::ordinal(FUNCommandArray),collapse = ", "), 
                           " functions in this plan and their direct defination on the argument will be included !\n"), 
                    SuppressWeb = TRUE)
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

