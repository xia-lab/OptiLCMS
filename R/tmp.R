
#' @export
controller.modifier <- function(new_command_set, last_command_set, plan){
  
  ###-------------Operators definition ------------//
  
  #' operators_1 : PerformDataTrimming ~ PerformParamsOptimization;
  #' operators_2 : PerformParamsOptimization ~ peak_profiling;  #' 
  #' operators_3 : PerformPeakProfiling ~ PerformPeakAnnotation;  #' 
  #' operators_4 : ;
  
  
  #' others_1 c1 : Control Optimization or not;
  #' others_1 c2 : ;
  #' others_1 c3 : ;
  #' others_1 c4 : .
  
  ##----------------------------------------------//
  
  if(class(new_command_set) == "OptiCommandSet"){
    
    # 1. ROIExtraction: -----
    # 1.1 Note on controller: c1, read; c2, trim; c3, write; c4, plot.
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
      plan@running.controller@ROI_extract[c(1:4)] <- rep(FALSE, 4);
    }
    
    if("datapath" %in% ChangedArugsArray){
      plan@running.controller@ROI_extract[c(1:4)] <- rep(TRUE,4);
      plan@running.controller@operators[1] <- TRUE;
    }
    
    if(any(c("mode", "mz", "mzdiff", "rt", "rtdiff", "rt.idx") %in% ChangedArugsArray)){
      plan@running.controller@ROI_extract[c(2:4)] <- rep(TRUE,3);
      plan@running.controller@operators[1] <- TRUE;
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
                                              "span","integrate","mzCenterFun","verbose.columns","fitgauss"))){
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
    
    if("Params" %in% ChangedArugsArray){
      
      newAnnoParamArgus <- new_PeakAnnotation$annotaParam;
      lastAnnoParamArgus <- last_PeakAnnotation$annotaParam;
      
      ChangedParamArgus <- ParamsChanged(lastAnnoParamArgus, newAnnoParamArgus);
      
      if(is.null(ChangedParamArgus)){
        plan@running.controller@peak_annotation[1] <- FALSE;
      } else {
        plan@running.controller@peak_annotation[1] <- TRUE;
      }
    }
    
  }
  
  if(plan@running.controller@operators[3]){
    plan@running.controller@peak_annotation[1] <- TRUE;
  }
  
  # 6. FormatPeakList: -----
  # 6.1 Note on controller FormatPeakList.
  # No need to resume for this function
  
  
  
  
  
  
  #' 
  #' ############### --------------- old code
  #' 
  #' functions <- sapply(1:length(plan[[length(plan)]]),FUN=function(i){plan[[length(plan)]][[i]][[3]][[1]]});
  #' data_trim_order <- which(functions=="PerformDataTrimming" | functions=="PerformROIExtraction");
  #' dara_ms_order <- which(functions=="ImportRawMSData");
  #' 
  #' ###-------------Operators definition ------------//
  #' 
  #' #' operators_1 : SetPeakParam ~ PerformParamsOptimization;
  #' #' operators_2 : data_folder_trainning ~ PerformDataTrimming;
  #' #' operators_3 : data_folder_sample ~ peak_profiling;
  #' #' operators_4 : PerformPeakProfiling ~ PerformPeakAnnotation;
  #' 
  #' 
  #' #' others_1 c1 : SetPeakParam;
  #' #' others_1 c2 : ImportRawData - reading part;
  #' #' others_1 c3 : ImportRawData - plotting part;
  #' #' others_1 c4 : Start to do annotation.
  #' 
  #' ##----------------------------------------------//
  #' 
  #' # Module 1 - Detect whether the tranning data folder has been changed
  #' if (!identical(data_trim_order, integer(0))){
  #'   
  #'   if ((class(new_command[[3]]) == "character") & (plan[[length(plan)]][[data_trim_order]][[3]][[2]] == new_command[[2]])){
  #'     if (new_command[[3]] != last_command[[3]]){
  #'       # to deal with the case the data folder did change
  #'       plan$running.controller$data_trim <- c(T,T,T,T);
  #'       names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'       plan$running.controller$operators[["operators_2"]] <- T;
  #'     } else if(OptiFileChanged(last_command[[3]])){
  #'       # to deal with the case some files are deleted or added for the same folder
  #'       plan$running.controller$data_trim <- c(T,T,T,T);
  #'       names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'       plan$running.controller$operators[["operators_2"]] <- T;
  #'     }
  #'   }
  #' }
  #' 
  #' # Module 2 - Detect whether the params of datatrimming has been changed
  #' if (new_command[[3]][[1]] == "PerformDataTrimming" | new_command[[3]][[1]] == "PerformROIExtraction"){
  #'   
  #'   data_folder_funciton_order <- which(sapply(1:length(plan[[length(plan)]]),
  #'                                              FUN = function(x){
  #'                                                plan[[length(plan)]][[x]][[2]]})==new_command[[3]][[2]] & 
  #'                                         as.character(sapply(1:length(plan[[length(plan)]]),
  #'                                                             FUN = function(x){plan[[length(plan)]][[x]][[1]]})) == "`<-`")
  #'   
  #'   if ((plan[[length(plan)]][[data_folder_funciton_order]] != plan[[length(plan)-1]][[data_folder_funciton_order]]) | 
  #'       plan$running.controller$operators[["operators_2"]]) {
  #'     # To make sure the data did change or not. If the data changed, re-run all.
  #'     plan$running.controller$data_trim <- c(T,T,T,T);
  #'     names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'     
  #'   } else if (new_command == last_command) {
  #'     # If the setting did not change, skip the trim step
  #'     if(.on.public.web){
  #'       # retained for further process 
  #'       # If on web, to detect whether the rmConts (remove contaminants param changed or not !)
  #'       print("run here : web version - param change detection !-")
  #'       
  #'       envir.path <- paste0(getwd(),"/temp/envir");
  #'       envir_tmp <- readRDS(paste0(envir.path,"/envir.rds"));
  #'       last_param <- envir_tmp[["param"]];
  #'       load("params.rda");
  #'       new_param <- peakParams;
  #'       
  #'       if(is.null(last_param)){ # This is designed for the case that param have not been finsihed when killed
  #'         
  #'         print("run here : param phase killed last time ! trim phase-")
  #'         plan$running.controller$data_trim <- c(T,T,T,T);
  #'         names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'         
  #'       } else if (last_param[["rmConts"]] != new_param[["rmConts"]]) {
  #'         
  #'         print("run here : param rmConts changed ! trim phase-")
  #'         plan$running.controller$data_trim <- c(F,T,T,T);
  #'         names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'       }
  #'       
  #'     } else {
  #'       # Otherwise, no need to detect
  #'       plan$running.controller$data_trim <- c(F,F,F,F);
  #'       names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'     }
  #'     
  #'   } else {
  #'     # If the setting did change, skip some steps in different cases
  #'     for (i in 2:length(new_command[[3]])){
  #'       for (j in 2:length(last_command[[3]])){
  #'         if ((names(new_command[[3]])[i] == names(last_command[[3]])[j]) && 
  #'             (new_command[[3]][[i]] != last_command[[3]][[j]])){
  #'           
  #'           if (names(new_command[[3]])[i] == "datapath" | (names(new_command[[3]])[i] == "" & i ==1)){
  #'             plan$running.controller$data_trim <- c(T,T,T,T);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "mode" | (names(new_command[[3]])[i] == "" & i ==2)){
  #'             plan$running.controller$data_trim <- c(F,T,T,T);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "mz" | (names(new_command[[3]])[i] == "" & i ==4)){
  #'             plan$running.controller$data_trim <- c(F,T,T,T);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "mzdiff" | (names(new_command[[3]])[i] == "" & i ==5)){
  #'             plan$running.controller$data_trim <- c(F,T,T,T);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "rt" | (names(new_command[[3]])[i] == "" & i ==6)){
  #'             plan$running.controller$data_trim <- c(F,T,T,T);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "rtdiff" | (names(new_command[[3]])[i] == "" & i ==7)){
  #'             plan$running.controller$data_trim <- c(F,T,T,T);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "rt.idx" | (names(new_command[[3]])[i] == "" & i ==8)){
  #'             plan$running.controller$data_trim <- c(F,T,T,T);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "write" | (names(new_command[[3]])[i] == "" & i ==3)){
  #'             plan$running.controller$data_trim <- c(F,F,T,F);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "plot" | (names(new_command[[3]])[i] == "" & i ==9)){
  #'             plan$running.controller$data_trim <- c(F,F,F,T);
  #'             names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
  #'           };
  #'         };
  #'       };
  #'     };
  #'   };
  #'   
  #'   
  #' };
  #' 
  #' if (new_command[[3]][[1]] == "SetPeakParam"){
  #'   
  #'   if (new_command != last_command) {
  #'     
  #'     plan$running.controller$others_1[1] <- T;
  #'     plan$running.controller$operators[["operators_1"]] <- T;
  #'     
  #'   } else if(new_command == last_command & .on.public.web) {
  #'     
  #'     #load("params.rda");
  #'     plan$running.controller$operators[["operators_1"]] <- F; # switch of the operator for optimization
  #'     plan$running.controller$others_1[1] <- F; # switch of the operator for optimization
  #'     
  #'   };
  #'   
  #' };
  #' 
  #' if (new_command[[3]][[1]] == "PerformParamsOptimization"){
  #'   
  #'   if (any(plan$running.controller$data_trim[c(1,2)])){
  #'     plan$running.controller$others_1[1] <- T;
  #'     
  #'   } else if (plan$running.controller$operators[["operators_1"]]) {
  #'     plan$running.controller$others_1[1] <- T;
  #'     
  #'   } else {
  #'     plan$running.controller$others_1[1] <- F;
  #'   };
  #'   
  #' };
  #' 
  #' # Module 3 - Dectect whether the raw ms data folder has been changed
  #' if (!identical(dara_ms_order, integer(0))){
  #'   if ((class(new_command[[3]]) == "character") & (plan[[length(plan)]][[dara_ms_order]][[3]][[2]] == new_command[[2]])){
  #'     if (new_command[[3]] != last_command[[3]] | ProcessFileChanged(last_command[[3]])){ 
  #'       # if the data folder changed at local or the files changed inside the folder
  #'       plan$running.controller$peak_profiling <- c(T,T,T,T);
  #'       names(plan$running.controller$peak_profiling) <- c("c1","c2","c3","c4"); 
  #'       plan$running.controller$operators[["operators_3"]] <- T;
  #'       
  #'     } else if(new_command[[3]] == last_command[[3]] & .on.public.web){ # to identify if the data folder changed at local
  #'       # TODO: need to correct this issue.
  #'       rawFileNames <- paste0(getwd(),"/temp/plan");
  #'       rawfilenms_last <- readRDS(paste0(rawFileNames,"/rawfilenms_",.plan_count-1,".rds"));
  #'       rawfilenms_new <- readRDS(paste0(rawFileNames,"/rawfilenms_",.plan_count,".rds"));
  #'       
  #'       if(identical(setdiff(rawfilenms_last,rawfilenms_new), character(0))){ 
  #'         # if files included didn't change
  #'         plan$running.controller$operators[["operators_3"]] <- F;
  #'         
  #'       } else {
  #'         # if files included did change
  #'         
  #'         #plan$running.controller$data_trim[["c1"]] <- T; 
  #'         #plan$running.controller$others_1[["c1"]] <- T; 
  #'         # TODO: to avoid re-do the trimming when QC did not change;
  #'         
  #'         plan$running.controller$operators[["operators_3"]] <- T;
  #'         plan$running.controller$peak_profiling <- c(T,T,T,T);
  #'         names(plan$running.controller$peak_profiling) <- c("c1","c2","c3","c4"); 
  #'         
  #'       }
  #'       
  #'     }
  #'   }
  #' };
  #' 
  #' # Module 4 - Detect whether the params of data Import has been changed
  #' if (new_command[[3]][[1]] == "ImportRawMSData"){
  #'   
  #'   data_folder_sample_order <- which(sapply(1:length(plan[[length(plan)]]),
  #'                                            FUN = function(x){plan[[length(plan)]][[x]][[2]]})==new_command[[3]][[2]] & 
  #'                                       as.character(sapply(1:length(plan[[length(plan)]]),
  #'                                                           FUN = function(x){plan[[length(plan)]][[x]][[1]]})) == "`<-`")
  #'   
  #'   if (!identical(which(names(new_command[[3]])=="plotSettings"),integer(0))){
  #'     plot_function1_order <-
  #'       which(new_command[[3]][[which(names(new_command[[3]]) == "plotSettings")]] == sapply(
  #'         plan[[length(plan)]],
  #'         FUN = function(x) {
  #'           x[[2]]
  #'         }
  #'       ))
  #'     
  #'   } else {
  #'     plot_function1_order <-
  #'       which(new_command[[3]][[4]] == sapply(
  #'         plan[[length(plan)]],
  #'         FUN = function(x) {
  #'           x[[2]]
  #'         }
  #'       ))
  #'   }
  #'   
  #'   if ((plan[[length(plan)]][[data_folder_sample_order]] != plan[[length(plan)-1]][[data_folder_sample_order]]) | 
  #'       plan$running.controller$operators[["operators_3"]]) {
  #'     # To make sure the data did change or not. If the data changed, re-run all.
  #'     # c2 in others_1 is used to control the ImportRawMSData - Reading Part (c3 is used for plotting part)
  #'     
  #'     plan[["running.controller"]][["others_1"]][["c2"]] <- T;
  #'     plan[["running.controller"]][["others_1"]][["c3"]] <- T;
  #'     
  #'   } else if (new_command == last_command) {
  #'     # If the setting did not change, skip the import step - Reading Part (c3 is used for plotting part)
  #'     
  #'     plan[["running.controller"]][["others_1"]][["c2"]] <- F;
  #'     
  #'   } else {
  #'     # If the setting did change, skip some steps in different cases
  #'     for (i in 2:length(new_command[[3]])){
  #'       for (j in 2:length(last_command[[3]])){
  #'         if ((names(new_command[[3]])[i] == names(last_command[[3]])[j]) && 
  #'             (new_command[[3]][[i]] != last_command[[3]][[j]])){
  #'           
  #'           if (names(new_command[[3]])[i] == "foldername" | (names(new_command[[3]])[i] == "" & i ==1)){
  #'             
  #'             plan[["running.controller"]][["others_1"]][["c2"]] <- T;
  #'             
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "mode" | (names(new_command[[3]])[i] == "" & i ==2)){
  #'             
  #'             plan[["running.controller"]][["others_1"]][["c2"]] <- T;
  #'             
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "ncores" | (names(new_command[[3]])[i] == "" & i ==3)){
  #'             
  #'             plan[["running.controller"]][["others_1"]][["c2"]] <- F;
  #'             
  #'           };
  #'           
  #'           if (names(new_command[[3]])[i] == "plotSettings" | (names(new_command[[3]])[i] == "" & i ==4)){
  #'             
  #'             plan[["running.controller"]][["others_1"]][["c2"]] <- F;
  #'             plan[["running.controller"]][["others_1"]][["c3"]] <- T;
  #'             
  #'           };
  #'         }
  #'       }
  #'     };
  #'     
  #'   }
  #'   
  #'   # others_1: c3 is the plotting part
  #'   if (plan[[length(plan)]][[plot_function1_order]] != plan[[length(plan)-1]][[plot_function1_order]]) {
  #'     
  #'     plan[["running.controller"]][["others_1"]][["c3"]] <- T;
  #'     
  #'   } else if(plan$running.controller$operators[["operators_3"]]){
  #'     
  #'     plan[["running.controller"]][["others_1"]][["c3"]] <- T;
  #'     
  #'   } else {
  #'     
  #'     plan[["running.controller"]][["others_1"]][["c3"]] <- F;
  #'     
  #'   }
  #'   
  #' }
  #' 
  #' # Module 5 - Detect whether the params of peak profiling +  annotation has been changed
  #' if (new_command[[3]][[1]] == "PerformPeakProfiling"){
  #'   print("run here PerformPeakProfiling + annotation parame changes detection !-");
  #'   
  #'   # TO get the position of functions defined the 'param'
  #'   profiling_param_order <- which(sapply(1:length(plan[[length(plan)]]),FUN = function(x){plan[[length(plan)]][[x]][[2]]})==new_command[[3]][[3]] & 
  #'                                    as.character(sapply(1:length(plan[[length(plan)]]),FUN = function(x){plan[[length(plan)]][[x]][[1]]})) == "`<-`")
  #'   
  #'   
  #'   if (!identical(which(names(new_command[[3]])=="plotSettings"),integer(0))){
  #'     plot_function2_order <- which(new_command[[3]][[which(names(new_command[[3]])=="plotSettings")]] == sapply(plan[[length(plan)]], FUN=function(x){x[[2]]}));
  #'   } else {
  #'     plot_function2_order <- which(new_command[[3]][[4]] == sapply(plan[[length(plan)]], FUN=function(x){x[[2]]}));
  #'   };
  #'   
  #'   if (plan$running.controller$others_1[[2]]){ 
  #'     # If data Import step (reading) was excecuted, the profiling step also has to be run/re-run;
  #'     plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
  #'     names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'     
  #'   } else if (.on.public.web & new_command == last_command) {
  #'     # If the param setting did change, run some profiling steps (web version)
  #'     
  #'     print("run here : web version - param change detection !-")
  #'     
  #'     envir.path <- paste0(getwd(),"/temp/envir");
  #'     envir_tmp <- readRDS(paste0(envir.path,"/envir.rds"));
  #'     last_param <- envir_tmp[["param"]];
  #'     load("params.rda");
  #'     new_param <- peakParams;
  #'     
  #'     if(is.null(last_param)){ # This is designed for the case that param have not been finsihed when killed
  #'       
  #'       print("run here : param phase killed last time !-")
  #'       
  #'       plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
  #'       names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'       
  #'       plan[["running.controller"]][["operators"]][[4]] <- T;
  #'       names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
  #'       
  #'     } else {
  #'       
  #'       identifiers <- NULL;
  #'       
  #'       for (i in 1:length(new_param)){
  #'         for (j in 1:length(last_param)){
  #'           if (names(new_param[i])==names(last_param[j]) & new_param[[i]] != last_param[[j]]){
  #'             identifiers <- c(identifiers, names(new_param[i]))
  #'           }
  #'         }
  #'       }
  #'       
  #'       switch.path <- paste0(getwd(),"/temp/plan");      
  #'       new_optimize_switch <- readRDS(paste0(switch.path,"/optimize_switch_",.plan_count,".rds"));
  #'       last_optimize_switch <- readRDS(paste0(switch.path,"/optimize_switch_",.plan_count-1,".rds"));
  #'       
  #'       if(new_optimize_switch == T & 
  #'          last_optimize_switch == T & 
  #'          !plan$running.controller$operators[["operators_2"]] & 
  #'          !plan$running.controller$operators[["operators_3"]]){
  #'         
  #'         identifiers <- NULL;
  #'         
  #'       } else if(plan$running.controller$operators[["operators_3"]] | plan$running.controller$operators[["operators_2"]]){
  #'         
  #'         # if data files included changed, re-run everything!
  #'         print("run here : data files included changed !")
  #'         
  #'         plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
  #'         names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'         
  #'         plan[["running.controller"]][["operators"]][[4]] <- T;
  #'         names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
  #'         
  #'       }
  #'       
  #'       if(is.null(identifiers)){
  #'         #0. No parameters changed
  #'         plan[["running.controller"]][["peak_profiling"]] <- c(F,F,F,F);
  #'         names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'         
  #'       } else if (any(identifiers %in% c("min_peakwidth","max_peakwidth","mzdiff","ppm","noise","prefilter","value_of_prefilter",
  #'                                         "Peak_method","snthresh","fwhm","sigma","steps"))){
  #'         
  #'         print("run here : picking parameters change found !-")
  #'         # 1. change picking parameters
  #'         plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
  #'         names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'         
  #'         plan[["running.controller"]][["operators"]][[4]] <- T;
  #'         names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
  #'         
  #'       } else if (any(identifiers %in% c("bw","RT_method","minFraction","minSamples","maxFeatures","family","smooth",
  #'                                         "span","integrate","mzCenterFun","verbose.columns","fitgauss"))){
  #'         
  #'         print("run here : alignment parameters change found !-")
  #'         # 2. change alignment parameters
  #'         plan[["running.controller"]][["peak_profiling"]] <- c(F,T,T,T);
  #'         names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'         
  #'         plan[["running.controller"]][["operators"]][[4]] <- T;
  #'         names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
  #'         
  #'       } else if(any(identifiers %in% c("polarity","perc_fwhm","mz_abs_iso","max_charge","max_iso","corr_eic_th","mz_abs_add"))){
  #'         
  #'         print("run here : Annotation parameters change found !-");
  #'         
  #'         # 3. change Annotation parameters
  #'         plan[["running.controller"]][["peak_profiling"]] <- c(F,F,F,F);
  #'         names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'         
  #'         plan[["running.controller"]][["operators"]][[4]] <- T;
  #'         names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
  #'         
  #'       };
  #'       
  #'     } 
  #'     
  #'   } else if (plan[[length(plan)]][[profiling_param_order]] != plan[[length(plan)-1]][[profiling_param_order]]) {
  #'     # If the param setting did change, run some profiling steps (Package version)
  #'     
  #'     identifiers <- profiling_param_identifier(plan[[length(plan)]][[profiling_param_order]],
  #'                                               plan[[length(plan)-1]][[profiling_param_order]]);
  #'     print("run here -8------------------------------")
  #'     # 1. change picking parameters
  #'     
  #'     if (any(identifiers %in% c("min_peakwidth","max_peakwidth","mzdiff","ppm","noise","prefilter","value_of_prefilter",
  #'                                "Peak_method","snthresh","fwhm","sigma","steps"))){
  #'       plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
  #'       names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'     };
  #'     
  #'     # 2. change alignment parameters
  #'     
  #'     if (any(identifiers %in% c("bw","RT_method","minFraction","minSamples","maxFeatures","family","smooth",
  #'                                "span","integrate","mzCenterFun","verbose.columns","fitgauss"))){
  #'       plan[["running.controller"]][["peak_profiling"]] <- c(F,T,T,T);
  #'       names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'     };
  #'     
  #'     
  #'   } else if (new_command == last_command) {
  #'     # If the setting did not change, skip the profiling step
  #'     
  #'     plan[["running.controller"]][["peak_profiling"]] <- c(F,F,F,F);
  #'     names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
  #'     
  #'   } else {
  #'     # Some unexcepted cases appears. Re-run everything to make sure correctness
  #'     print("run here : unexcepted cases appears - Re-run everything !-")
  #'     plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
  #'   };
  #'   
  #'   # peak_profiling: c4 is the plotting part
  #'   if ((plan[[length(plan)]][[plot_function2_order]] != plan[[length(plan)-1]][[plot_function2_order]]) & 
  #'       (!any(plan[["running.controller"]][["peak_profiling"]][c(1:3)]))) {
  #'     
  #'     print("run here -10------------------------------")
  #'     plan[["running.controller"]][["peak_profiling"]][[4]] <- T;
  #'     
  #'   } else if(any(plan[["running.controller"]][["peak_profiling"]][c(1:3)])) {
  #'     
  #'     print("run here -11------------------------------")
  #'     plan[["running.controller"]][["peak_profiling"]][[4]] <- T;
  #'     
  #'   } else {
  #'     plan[["running.controller"]][["peak_profiling"]][[4]] <- F;
  #'   };
  #'   
  #' }
  #' 
  #' # Module 6 - Detect whether annotation need to be re-run
  #' if (new_command[[3]][[1]] == "SetAnnotationParam"){
  #'   
  #'   if (new_command != last_command) {
  #'     plan$running.controller$operators[[4]] <- T; # operators_4 for annotations
  #'   };
  #'   
  #'   if (any(plan$running.controller$peak_profiling[1:3])){
  #'     plan$running.controller$operators[[4]] <- T;
  #'   };
  #'   
  #' }
  #' 
  #' # Return the prepared plan for excuting
  return(plan)
  
  # 
}
