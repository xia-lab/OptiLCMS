## Perform Running Pipeline - script

#' Initializing running plan
#' @param type 
#'
#' @param path 
#'
#' @export
InitializaPlan <- function(type="spec", path){
  if(path != ""){
    setwd(path);
    fullUserPath <<- path;
  } else {
    fullUserPath <<- "";
  };
  
  if (.on.public.web){
    
    if (!file.exists(paste0(fullUserPath, "log_progress.txt"))){
      write.table(0.0, file = paste0(fullUserPath, "log_progress.txt"),row.names = F,col.names = F);
    }
    if (!file.exists(paste0(fullUserPath, "opt_param.txt"))){
      write.table(0.0, file = paste0(fullUserPath, "opt_param.txt"),row.names = T,col.names = F);
    }
    
    if (!file.exists("metaboanalyst_spec_proc.txt")){
      print_mes <- paste0("Running Status -- Job Submitted Successfully at: ", Sys.time(), "\n Waiting in queue to start...")
      write.table(print_mes, file = paste0(fullUserPath, "metaboanalyst_spec_proc.txt"),row.names = F,col.names = F);
    }
    
  }
  
  #================= raw_pre
  if (type=="spec"){
    
    plan.path <- paste0(getwd(),"/temp/plan");
    if (!dir.exists(plan.path)){dir.create(paste0(getwd(),"/temp/plan"),recursive = T)};
    
    if (!file.exists(paste0(plan.path,"/plan.rds"))){
      
      running.as.plan <<- T;
      plan <- list();
      plan_count <<- 0;
      
      
      saveRDS(plan, file = paste0(plan.path,"/plan.rds"));
      saveRDS(running.as.plan, file = paste0(plan.path,"/running.as.plan.rds"));
      saveRDS(plan_count, file = paste0(plan.path,"/plan_count.rds"));
      
      
    } else {
      
      plan <<- readRDS(paste0(plan.path,"/plan.rds"));
      running.as.plan <<- readRDS(paste0(plan.path,"/running.as.plan.rds"));
      plan_count <<- readRDS(paste0(plan.path,"/plan_count.rds"));
      plan_count <<- plan_count + 1;
      
    };
    
    #----------------------
    
    envir.path <- paste0(getwd(),"/temp/envir");
    if (!dir.exists(envir.path)){dir.create(paste0(getwd(),"/temp/envir"),recursive = T)};
    
    if (!file.exists(paste0(envir.path,"/envir.rds"))){
      envir <<- new.env();
      saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
    } else {
      envir <<- readRDS(paste0(envir.path,"/envir.rds"));
    }
    
    if(exists("plan_count") & plan_count > 0){
      return(plan)
    }
    
  };
  
  
  #=============== raw_opt
  if (type=="raw_opt"){
    
    plan.path <- paste0(getwd(),"/temp/plan");
    if (!dir.exists(plan.path)){dir.create(paste0(getwd(),"/temp/plan"),recursive = T)};
    
    if (!file.exists(paste0(plan.path,"/plan.rds"))){
      
      running.as.plan <<- T;
      plan <- list();
      plan_count <<- 0;
      
      
      saveRDS(plan, file = paste0(plan.path,"/plan.rds"));
      saveRDS(running.as.plan, file = paste0(plan.path,"/running.as.plan.rds"));
      saveRDS(plan_count, file = paste0(plan.path,"/plan_count.rds"));
      
      
    } else {
      
      plan <<- readRDS(paste0(plan.path,"/plan.rds"));
      running.as.plan <<- readRDS(paste0(plan.path,"/running.as.plan.rds"));
      plan_count <<- readRDS(paste0(plan.path,"/plan_count.rds"));
      plan_count <<- plan_count + 1;
      
    };
    
    #----------------------
    
    optimize_switch <<- T;
    switch.path <- paste0(getwd(),"/temp/plan");
    saveRDS(optimize_switch, file = paste0(switch.path,"/optimize_switch_",plan_count,".rds"));  
    
    #----------------------
    
    envir.path <- paste0(getwd(),"/temp/envir");
    if (!dir.exists(envir.path)){dir.create(paste0(getwd(),"/temp/envir"),recursive = T)};
    
    if (!file.exists(paste0(envir.path,"/envir.rds"))){
      envir <<- new.env();
      saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
    } else {
      envir <<- readRDS(paste0(envir.path,"/envir.rds"));
    }
    
    #----------------------
    
    rawFileNames <- paste0(getwd(),"/temp/plan");
    saveRDS(rawfilenms, file=paste0(rawFileNames,"/rawfilenms_",plan_count,".rds"))
    
    #----------------------
    
    if(exists("plan_count") & plan_count > 0){
      return(plan)
    }
    
  };
  
  #================== raw_ms
  if (type=="raw_ms"){
    
    plan.path <- paste0(getwd(),"/temp/plan");
    
    if (!dir.exists(plan.path)){dir.create(paste0(getwd(),"/temp/plan"),recursive = T)};
    
    #---------------
    if (!file.exists(paste0(plan.path,"/plan.rds"))){
      
      running.as.plan <<- T;
      plan <- list();
      plan_count <<- 0;
      
      saveRDS(plan, file = paste0(plan.path,"/plan.rds"));
      saveRDS(running.as.plan, file = paste0(plan.path,"/running.as.plan.rds"));
      saveRDS(plan_count, file = paste0(plan.path,"/plan_count.rds"));
      
      
    } else {
      
      plan <<- readRDS(paste0(plan.path,"/plan.rds"));
      running.as.plan <<- readRDS(paste0(plan.path,"/running.as.plan.rds"));
      plan_count <<- readRDS(paste0(plan.path,"/plan_count.rds"));
      plan_count <<- plan_count + 1;
      
    };
    
    #----------------------
    
    optimize_switch <<- F;
    switch.path <- paste0(getwd(),"/temp/plan");
    saveRDS(optimize_switch, file = paste0(switch.path,"/optimize_switch_",plan_count,".rds"));
    
    #----------------------
    
    envir.path <- paste0(getwd(),"/temp/envir");
    if (!dir.exists(envir.path)){dir.create(paste0(getwd(),"/temp/envir"),recursive = T)};
    if (!file.exists(paste0(envir.path,"/envir.rds"))){
      envir <<- new.env();
      saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
    } else {
      envir <<- readRDS(paste0(envir.path,"/envir.rds"));
    }
    
    #----------------------
    
    rawFileNames <- paste0(getwd(),"/temp/plan");
    saveRDS(rawfilenms, file=paste0(rawFileNames,"/rawfilenms_",plan_count,".rds"))
    
    #----------------------
    
    if(exists("plan_count") & plan_count > 0){
      return(plan)
    }
    
  };
  
  ## Recording cache file information
  record.info <- matrix(nrow = 180,ncol = 2);record.path <- paste0(getwd(),"/temp/records");
  if (!dir.exists(record.path)){dir.create(paste0(getwd(),"/temp/records"),recursive = T)};
  
  if (!file.exists(paste0(record.path,"/records.rds"))){
    saveRDS(record.info,file = paste0(record.path,"/records.rds"))
  } else {
    records <- readRDS(paste0(record.path,"/records.rds"));
    if (sum(is.na(records[,1])) < 1){
      AddErrMsg("Too many caches for data proceesing operation !")
    };
  };
  
  
  return(plan)
  #return(.set.mSet(plan))
}

#' running.plan
#' @param plan 
#'
#' @param ... 
#'
#' @export
running.plan <- function(plan=NULL,...){
  
  #plan <- .get.current.plan(plan);
  commands <- match.call(expand.dots = FALSE)$...
  
  ## Declare controller
  plan$running.controller <- controller.resetter()
  ##
  
  if (!length(commands)) {
    stop("No command provided to run !")
  }
  
  plan[[paste0("command_set_",plan_count)]] <- commands
  
  if(.on.public.web){
    plan.path <- paste0(getwd(),"/temp/plan");
    saveRDS(plan, file = paste0(plan.path,"/plan.rds"));
    saveRDS(plan_count, file = paste0(plan.path,"/plan_count.rds"));
  }
  
  recordMarker_resetter(plan_count);
  
  return(plan)
  
}
.get.current.plan <- function(plan){
  
  if (is.null(plan)){
    return(list())
  } else {
    return (plan)
  }
  
}

#' ExecutePlan
#' @param plan 
#'
#' @export
ExecutePlan <- function(plan=NULL){
  print(getwd());
  if (is.null(plan)){
    stop("No running plan input. Please design you plan first with 'running.plan' function !")
  };
  
  
  # Reset running.controller to make sure everything is normal at beginning
  plan$running.controller <- controller.resetter()
  
  if (length(plan)==2){
    
    envir$rc <<- plan$running.controller
    perform.plan(plan[["command_set_1"]]);
    
    envir.path <- paste0(getwd(),"/temp/envir");
    saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
    
  } else if (length(plan) > 2) {
    
    ## This is the most important part for the whole pipeline
    ## Module 1-6: Dectect whether the parameters changed or not so as to operate the whole pipeline
    
    last_command_set <- plan[[names(plan)[length(plan)-1]]];
    new_command_set <- plan[[names(plan)[length(plan)]]];
    
    for (i in 1:length(new_command_set)){
      for (j in 1:length(last_command_set)){
        
        if (((new_command_set[[i]][[1]] == "<-") & (last_command_set[[j]][[1]] == "<-")) & 
            (new_command_set[[i]][[3]][[1]] == last_command_set[[j]][[3]][[1]])){
          plan <- controller.modifier(new_command_set[[i]],last_command_set[[j]],plan)
        };
        
        if (((new_command_set[[i]][[1]] != "<-") & (last_command_set[[j]][[1]] != "<-")) &
            (new_command_set[[i]][[1]] == last_command_set[[j]][[1]])){
          plan <- controller.modifier(new_command_set[[i]],last_command_set[[j]],plan)
        };
        
      }
    }
    
    ## Module 7 - Detect whether some steps have been run in last plan excuting process (Secondary Decision switch)
    plan <- recording_identifier(plan)
    
    
    ## Module 8 - Dectect whether current plan type (raw_ms or raw_pre) is different from the last one (Final Decision switch)
    plan <- planType_identifier(plan)
    
    
    
    envir$rc <<- plan$running.controller
    # define.plan.controller <- str2lang('rc <- plan$running.controller');
    # eval (define.plan.controller,envir = envir);
    #perform.plan(new_command_set);
    
    mSetInfo <- tryCatch(perform.plan(new_command_set), error = function(e){e})
    
    if (class(mSetInfo)[1]=="simpleError"){
      #write.table(0.0, file = paste0(fullUserPath, "log_progress.txt"),row.names = F,col.names = F);
      print_mes <- paste0("<font color=\"red\">","\nERROR:",mSetInfo$message,"</font>");
      write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
      stop("EXCEPTION POINT CODE: EX1");
    }
    
    envir.path <- paste0(getwd(),"/temp/envir");
    saveRDS(envir, file = paste0(envir.path,"/envir.rds"));
    
  } else {
    stop("Wrong plan has been prepared !")
  }
  
}

controller.modifier <- function(new_command,last_command,plan){
  
  functions <- sapply(1:length(plan[[length(plan)]]),FUN=function(i){plan[[length(plan)]][[i]][[3]][[1]]});
  data_trim_order <- which(functions=="PerformDataTrimming");
  dara_ms_order <- which(functions=="ImportRawMSData");
  
  ###-------------Operators definition ------------//
  
  #' operators_1 : SetPeakParam ~ PerformParamsOptimization;
  #' operators_2 : data_folder_trainning ~ PerformDataTrimming;
  #' operators_3 : data_folder_sample ~ peak_profiling;
  #' operators_4 : PerformPeakProfiling ~ PerformPeakAnnotation;
  
  
  #' others_1 c1 : SetPeakParam;
  #' others_1 c2 : ImportRawData - reading part;
  #' others_1 c3 : ImportRawData - plotting part;
  #' others_1 c4 : Start to do annotation.
  
  ##----------------------------------------------//
  
  # Module 1 - Detect whether the tranning data folder has been changed
  if (!identical(data_trim_order, integer(0))){
    if ((class(new_command[[3]]) == "character") & (plan[[length(plan)]][[data_trim_order]][[3]][[2]] == new_command[[2]])){
      if (new_command[[3]] != last_command[[3]]){
        plan$running.controller$data_trim <- c(T,T,T,T);
        names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
        plan$running.controller$operators[["operators_2"]] <- T;
      }
    }
  }
  
  # Module 2 - Detect whether the params of datatrimming has been changed
  if (new_command[[3]][[1]] == "PerformDataTrimming"){
    
    data_folder_funciton_order <- which(sapply(1:length(plan[[length(plan)]]),FUN = function(x){plan[[length(plan)]][[x]][[2]]})==new_command[[3]][[2]] & 
                                          as.character(sapply(1:length(plan[[length(plan)]]),FUN = function(x){plan[[length(plan)]][[x]][[1]]})) == "`<-`")
    
    if ((plan[[length(plan)]][[data_folder_funciton_order]] != plan[[length(plan)-1]][[data_folder_funciton_order]]) | 
        plan$running.controller$operators[["operators_2"]]){
      # To make sure the data did change or not. If the data changed, re-run all.
      plan$running.controller$data_trim <- c(T,T,T,T);
      names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
      
    } else if (new_command == last_command) {
      # If the setting did not change, skip the trim step
      if(.on.public.web){
        # If on web, to detect whether the rmConts (remove contaminants param changed or not !)
        print("run here : web version - param change detection !-")
        
        envir.path <- paste0(getwd(),"/temp/envir");
        envir_tmp <- readRDS(paste0(envir.path,"/envir.rds"));
        last_param <- envir_tmp[["param"]];
        load("params.rda");
        new_param <- peakParams;
        
        if(is.null(last_param)){ # This is designed for the case that param have not been finsihed when killed
          
          print("run here : param phase killed last time ! trim phase-")
          plan$running.controller$data_trim <- c(T,T,T,T);
          names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
          
        } else if (last_param[["rmConts"]] != new_param[["rmConts"]]) {
          
          print("run here : param rmConts changed ! trim phase-")
          plan$running.controller$data_trim <- c(F,T,T,T);
          names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
        }
        
        
        
      } else {
        # Otherwise, no need to detect
        plan$running.controller$data_trim <- c(F,F,F,F);
        names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
      }
      
      
    } else {
      # If the setting did change, skip some steps in different cases
      for (i in 2:length(new_command[[3]])){
        for (j in 2:length(last_command[[3]])){
          if ((names(new_command[[3]])[i] == names(last_command[[3]])[j]) && 
              (new_command[[3]][[i]] != last_command[[3]][[j]])){
            
            if (names(new_command[[3]])[i] == "datapath" | (names(new_command[[3]])[i] == "" & i ==1)){
              plan$running.controller$data_trim <- c(T,T,T,T);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
            
            if (names(new_command[[3]])[i] == "mode" | (names(new_command[[3]])[i] == "" & i ==2)){
              plan$running.controller$data_trim <- c(F,T,T,T);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
            
            if (names(new_command[[3]])[i] == "mz" | (names(new_command[[3]])[i] == "" & i ==4)){
              plan$running.controller$data_trim <- c(F,T,T,T);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
            
            if (names(new_command[[3]])[i] == "mzdiff" | (names(new_command[[3]])[i] == "" & i ==5)){
              plan$running.controller$data_trim <- c(F,T,T,T);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
            
            if (names(new_command[[3]])[i] == "rt" | (names(new_command[[3]])[i] == "" & i ==6)){
              plan$running.controller$data_trim <- c(F,T,T,T);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
            
            if (names(new_command[[3]])[i] == "rtdiff" | (names(new_command[[3]])[i] == "" & i ==7)){
              plan$running.controller$data_trim <- c(F,T,T,T);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
            
            if (names(new_command[[3]])[i] == "rt.idx" | (names(new_command[[3]])[i] == "" & i ==8)){
              plan$running.controller$data_trim <- c(F,T,T,T);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
            
            if (names(new_command[[3]])[i] == "write" | (names(new_command[[3]])[i] == "" & i ==3)){
              plan$running.controller$data_trim <- c(F,F,T,F);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
            
            if (names(new_command[[3]])[i] == "plot" | (names(new_command[[3]])[i] == "" & i ==9)){
              plan$running.controller$data_trim <- c(F,F,F,T);
              names(plan$running.controller$data_trim) <- c("c1","c2","c3","c4");
            };
          };
        };
      };
    };
    
    
  };
  
  if (new_command[[3]][[1]] == "SetPeakParam"){
    
    if (new_command != last_command) {
      
      plan$running.controller$others_1[1] <- T;
      plan$running.controller$operators[["operators_1"]] <- T;
      
    } else if(new_command == last_command & .on.public.web) {
      
      #load("params.rda");
      plan$running.controller$operators[["operators_1"]] <- F; # switch of the operator for optimization
      plan$running.controller$others_1[1] <- F; # switch of the operator for optimization
      
    };
    
  };
  
  if (new_command[[3]][[1]] == "PerformParamsOptimization"){
    
    if (any(plan$running.controller$data_trim[c(1,2)])){
      plan$running.controller$others_1[1] <- T;
      
    } else if (plan$running.controller$operators[["operators_1"]]) {
      plan$running.controller$others_1[1] <- T;
      
    } else {
      plan$running.controller$others_1[1] <- F;
    };
    
  };
  
  # Module 3 - Dectect whether the raw ms data folder has been changed
  if (!identical(dara_ms_order, integer(0))){
    if ((class(new_command[[3]]) == "character") & (plan[[length(plan)]][[dara_ms_order]][[3]][[2]] == new_command[[2]])){
      
      if (new_command[[3]] != last_command[[3]]){ # if the data folder changed at local
        
        plan$running.controller$peak_profiling <- c(T,T,T,T);
        names(plan$running.controller$peak_profiling) <- c("c1","c2","c3","c4"); 
        plan$running.controller$operators[["operators_3"]] <- T;
        
      } else if(new_command[[3]] == last_command[[3]] & .on.public.web){ # to identify if the data folder changed at local
        
        rawFileNames <- paste0(getwd(),"/temp/plan");
        rawfilenms_last <- readRDS(paste0(rawFileNames,"/rawfilenms_",plan_count-1,".rds"));
        rawfilenms_new <- readRDS(paste0(rawFileNames,"/rawfilenms_",plan_count,".rds"));
        
        if(identical(setdiff(rawfilenms_last,rawfilenms_new), character(0))){ # if files included didn't change
          
          plan$running.controller$operators[["operators_3"]] <- F;
          
        } else {# if files included did change
          
          print("run here -0-0------------------------------");
          
          #plan$running.controller$data_trim[["c1"]] <- T; 
          #plan$running.controller$others_1[["c1"]] <- T; # TO DO: to avoid re-do the trimming when QC did not change;
          
          plan$running.controller$operators[["operators_3"]] <- T;
          plan$running.controller$peak_profiling <- c(T,T,T,T);
          names(plan$running.controller$peak_profiling) <- c("c1","c2","c3","c4"); 
          
        }
        
      }
    }
  };
  
  # Module 4 - Detect whether the params of data Import has been changed
  if (new_command[[3]][[1]] == "ImportRawMSData"){
    
    data_folder_sample_order <- which(sapply(1:length(plan[[length(plan)]]),FUN = function(x){plan[[length(plan)]][[x]][[2]]})==new_command[[3]][[2]] & 
                                        as.character(sapply(1:length(plan[[length(plan)]]),FUN = function(x){plan[[length(plan)]][[x]][[1]]})) == "`<-`")
    
    if (!identical(which(names(new_command[[3]])=="plotSettings"),integer(0))){
      print("run here -0-1------------------------------")
      plot_function1_order <- which(new_command[[3]][[which(names(new_command[[3]])=="plotSettings")]] == sapply(plan[[length(plan)]], FUN=function(x){x[[2]]}));
      
    } else {
      print("run here -0-2------------------------------")
      plot_function1_order <- which(new_command[[3]][[4]] == sapply(plan[[length(plan)]], FUN=function(x){x[[2]]}));
      
    }
    
    if ((plan[[length(plan)]][[data_folder_sample_order]] != plan[[length(plan)-1]][[data_folder_sample_order]]) | 
        plan$running.controller$operators[["operators_3"]]){
      # To make sure the data did change or not. If the data changed, re-run all.
      # c2 in others_1 is used to control the ImportRawMSData - Reading Part (c3 is used for plotting part)
      print("run here -0-3------------------------------");
      
      plan[["running.controller"]][["others_1"]][["c2"]] <- T;
      plan[["running.controller"]][["others_1"]][["c3"]] <- T;
      
    } else if (new_command == last_command) {
      # If the setting did not change, skip the import step - Reading Part (c3 is used for plotting part)
      print("run here -0-4------------------------------");
      
      plan[["running.controller"]][["others_1"]][["c2"]] <- F;
      
    } else {
      print("run here -0-5------------------------------")
      # If the setting did change, skip some steps in different cases
      for (i in 2:length(new_command[[3]])){
        for (j in 2:length(last_command[[3]])){
          if ((names(new_command[[3]])[i] == names(last_command[[3]])[j]) && 
              (new_command[[3]][[i]] != last_command[[3]][[j]])){
            
            if (names(new_command[[3]])[i] == "foldername" | (names(new_command[[3]])[i] == "" & i ==1)){
              
              plan[["running.controller"]][["others_1"]][["c2"]] <- T;
              
            };
            
            if (names(new_command[[3]])[i] == "mode" | (names(new_command[[3]])[i] == "" & i ==2)){
              
              plan[["running.controller"]][["others_1"]][["c2"]] <- T;
              
            };
            
            if (names(new_command[[3]])[i] == "ncores" | (names(new_command[[3]])[i] == "" & i ==3)){
              
              plan[["running.controller"]][["others_1"]][["c2"]] <- F;
              
            };
            
            if (names(new_command[[3]])[i] == "plotSettings" | (names(new_command[[3]])[i] == "" & i ==4)){
              
              plan[["running.controller"]][["others_1"]][["c2"]] <- F;
              plan[["running.controller"]][["others_1"]][["c3"]] <- T;
              
            };
          }
        }
      };
      
    }
    
    # others_1: c3 is the plotting part
    if (plan[[length(plan)]][[plot_function1_order]] != plan[[length(plan)-1]][[plot_function1_order]]) {
      
      print("run here -0-6------------------------------")
      plan[["running.controller"]][["others_1"]][["c3"]] <- T;
      
      
    } else if(plan$running.controller$operators[["operators_3"]]){
      print("run here -0-7------------------------------")
      plan[["running.controller"]][["others_1"]][["c3"]] <- T;
      
    } else {
      print("run here -0-8------------------------------")
      plan[["running.controller"]][["others_1"]][["c3"]] <- F;
      
    }
    
  }
  
  # Module 5 - Detect whether the params of peak profiling +  annotation has been changed
  if (new_command[[3]][[1]] == "PerformPeakProfiling"){
    print("run here PerformPeakProfiling + annotation parame changes detection !-");
    
    # TO get the position of functions defined the 'param'
    profiling_param_order <- which(sapply(1:length(plan[[length(plan)]]),FUN = function(x){plan[[length(plan)]][[x]][[2]]})==new_command[[3]][[3]] & 
                                     as.character(sapply(1:length(plan[[length(plan)]]),FUN = function(x){plan[[length(plan)]][[x]][[1]]})) == "`<-`")
    
    
    if (!identical(which(names(new_command[[3]])=="plotSettings"),integer(0))){
      print("run here -2------------------------------")
      plot_function2_order <- which(new_command[[3]][[which(names(new_command[[3]])=="plotSettings")]] == sapply(plan[[length(plan)]], FUN=function(x){x[[2]]}));
    } else {
      print("run here -3------------------------------")
      plot_function2_order <- which(new_command[[3]][[4]] == sapply(plan[[length(plan)]], FUN=function(x){x[[2]]}));
    };
    
    if (plan$running.controller$others_1[[2]]){ 
      # If data Import step (reading) was excecuted, the profiling step also has to be run/re-run;
      print("run here -4------------------------------")
      plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
      names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
      
    } else if (.on.public.web & new_command == last_command) {
      # If the param setting did change, run some profiling steps (web version)
      
      print("run here : web version - param change detection !-")
      
      envir.path <- paste0(getwd(),"/temp/envir");
      envir_tmp <- readRDS(paste0(envir.path,"/envir.rds"));
      last_param <- envir_tmp[["param"]];
      load("params.rda");
      new_param <- peakParams;
      
      if(is.null(last_param)){ # This is designed for the case that param have not been finsihed when killed
        
        print("run here : param phase killed last time !-")
        
        plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
        names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
        
        plan[["running.controller"]][["operators"]][[4]] <- T;
        names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
        
      } else {
        
        identifiers <- NULL;
        
        for (i in 1:length(new_param)){
          for (j in 1:length(last_param)){
            if (names(new_param[i])==names(last_param[j]) & new_param[[i]] != last_param[[j]]){
              identifiers <- c(identifiers, names(new_param[i]))
            }
          }
        }
        
        switch.path <- paste0(getwd(),"/temp/plan");      
        new_optimize_switch <- readRDS(paste0(switch.path,"/optimize_switch_",plan_count,".rds"));
        last_optimize_switch <- readRDS(paste0(switch.path,"/optimize_switch_",plan_count-1,".rds"));
        
        if(new_optimize_switch == T & 
           last_optimize_switch == T & 
           !plan$running.controller$operators[["operators_2"]] & 
           !plan$running.controller$operators[["operators_3"]]){
          
          identifiers <- NULL;
          
        } else if(plan$running.controller$operators[["operators_3"]] | plan$running.controller$operators[["operators_2"]]){
          
          # if data files included changed, re-run everything!
          print("run here : data files included changed !")
          
          plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
          names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
          
          plan[["running.controller"]][["operators"]][[4]] <- T;
          names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
          
        }
        
        if(is.null(identifiers)){
          #0. No parameters changed
          plan[["running.controller"]][["peak_profiling"]] <- c(F,F,F,F);
          names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
          
        } else if (any(identifiers %in% c("min_peakwidth","max_peakwidth","mzdiff","ppm","noise","prefilter","value_of_prefilter",
                                          "Peak_method","snthresh","fwhm","sigma","steps"))){
          
          print("run here : picking parameters change found !-")
          # 1. change picking parameters
          plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
          names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
          
          plan[["running.controller"]][["operators"]][[4]] <- T;
          names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
          
        } else if (any(identifiers %in% c("bw","RT_method","minFraction","minSamples","maxFeatures","family","smooth",
                                          "span","integrate","mzCenterFun","verbose.columns","fitgauss"))){
          
          print("run here : alignment parameters change found !-")
          # 2. change alignment parameters
          plan[["running.controller"]][["peak_profiling"]] <- c(F,T,T,T);
          names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
          
          plan[["running.controller"]][["operators"]][[4]] <- T;
          names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
          
        } else if(any(identifiers %in% c("polarity","perc_fwhm","mz_abs_iso","max_charge","max_iso","corr_eic_th","mz_abs_add"))){
          
          print("run here : Annotation parameters change found !-");
          
          # 3. change Annotation parameters
          plan[["running.controller"]][["peak_profiling"]] <- c(F,F,F,F);
          names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
          
          plan[["running.controller"]][["operators"]][[4]] <- T;
          names(plan[["running.controller"]][["operators"]][[4]]) <- "operators_4";
          
        };
        
      } 
      
    } else if (plan[[length(plan)]][[profiling_param_order]] != plan[[length(plan)-1]][[profiling_param_order]]) {
      # If the param setting did change, run some profiling steps (Package version)
      
      identifiers <- profiling_param_identifier(plan[[length(plan)]][[profiling_param_order]],
                                                plan[[length(plan)-1]][[profiling_param_order]]);
      print("run here -8------------------------------")
      # 1. change picking parameters
      
      if (any(identifiers %in% c("min_peakwidth","max_peakwidth","mzdiff","ppm","noise","prefilter","value_of_prefilter",
                                 "Peak_method","snthresh","fwhm","sigma","steps"))){
        plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
        names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
      };
      
      # 2. change alignment parameters
      
      if (any(identifiers %in% c("bw","RT_method","minFraction","minSamples","maxFeatures","family","smooth",
                                 "span","integrate","mzCenterFun","verbose.columns","fitgauss"))){
        plan[["running.controller"]][["peak_profiling"]] <- c(F,T,T,T);
        names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
      };
      
      
    } else if (new_command == last_command) {
      # If the setting did not change, skip the profiling step
      
      plan[["running.controller"]][["peak_profiling"]] <- c(F,F,F,F);
      names(plan[["running.controller"]][["peak_profiling"]]) <- c("c1","c2","c3","c4");
      
    } else {
      # Some unexcepted cases appears. Re-run everything to make sure correctness
      print("run here : unexcepted cases appears - Re-run everything !-")
      plan[["running.controller"]][["peak_profiling"]] <- c(T,T,T,T);
    };
    
    # peak_profiling: c4 is the plotting part
    if ((plan[[length(plan)]][[plot_function2_order]] != plan[[length(plan)-1]][[plot_function2_order]]) & 
        (!any(plan[["running.controller"]][["peak_profiling"]][c(1:3)]))) {
      
      print("run here -10------------------------------")
      plan[["running.controller"]][["peak_profiling"]][[4]] <- T;
      
    } else if(any(plan[["running.controller"]][["peak_profiling"]][c(1:3)])) {
      
      print("run here -11------------------------------")
      plan[["running.controller"]][["peak_profiling"]][[4]] <- T;
      
    } else {
      plan[["running.controller"]][["peak_profiling"]][[4]] <- F;
    };
    
  }
  
  # Module 6 - Detect whether annotation need to be re-run
  if (new_command[[3]][[1]] == "SetAnnotationParam"){
    
    if (new_command != last_command) {
      plan$running.controller$operators[[4]] <- T; # operators_4 for annotations
    };
    
    if (any(plan$running.controller$peak_profiling[1:3])){
      plan$running.controller$operators[[4]] <- T;
    };
    
  }
  
  # Return the prepared plan for excuting
  return(plan)
  
}

planType_identifier <- function(plan){
  
  switch.path <- paste0(getwd(),"/temp/plan");
  optimize_switch_last <- readRDS(paste0(switch.path,"/optimize_switch_",plan_count-1,".rds"));
  
  if  (optimize_switch_last == optimize_switch){
    # if the plan type is the same, not modification
    return(plan)
  } else {
    # otherwise, rerun everything
    plan[["running.controller"]][["data_trim"]][["c1"]] <- 
      plan[["running.controller"]][["data_trim"]][["c2"]] <- 
      plan[["running.controller"]][["data_trim"]][["c3"]] <- 
      plan[["running.controller"]][["data_trim"]][["c4"]] <- 
      plan[["running.controller"]][["peak_profiling"]][["c1"]] <-
      plan[["running.controller"]][["peak_profiling"]][["c2"]] <- 
      plan[["running.controller"]][["peak_profiling"]][["c3"]] <- 
      plan[["running.controller"]][["peak_profiling"]][["c4"]] <- 
      plan[["running.controller"]][["others_1"]][["c1"]] <-
      plan[["running.controller"]][["others_1"]][["c2"]] <- 
      plan[["running.controller"]][["others_1"]][["c3"]] <- 
      plan[["running.controller"]][["others_1"]][["c4"]] <- 
      plan[["running.controller"]][["operators"]][["operators_4"]] <-T; # special mark- TO DO: improve in future
    
    plan[["running.controller"]][["operators"]][["operators_1"]] <- 
      plan[["running.controller"]][["operators"]][["operators_2"]] <- 
      plan[["running.controller"]][["operators"]][["operators_3"]] <-
      plan[["running.controller"]][["operators"]][["operators_5"]] <-
      plan[["running.controller"]][["operators"]][["operators_6"]] <-
      plan[["running.controller"]][["operators"]][["operators_7"]] <-
      plan[["running.controller"]][["operators"]][["operators_8"]] <- F;
    
    return(plan);
    
  }
  
}

recording_identifier <- function(plan) {
  
  record.path <- paste0(getwd(),"/temp/records");
  
  if(file.exists(paste0(record.path,"/records_marker_",plan_count-1,".rds"))){
    record.marker_last <- readRDS(paste0(record.path,"/records_marker_",plan_count-1,".rds"));
  } else {
    recordMarker_resetter(plan_count-1);
    record.marker_last <- readRDS(paste0(record.path,"/records_marker_",plan_count-1,".rds"));
  }
  
  
  # If things were not run last time, forcedly start running the slot this time
  if(!as.logical(record.marker_last[1,2])){plan$running.controller$data_trim[1] <- T}
  if(!as.logical(record.marker_last[2,2])){plan$running.controller$data_trim[2] <- T}
  if(!as.logical(record.marker_last[3,2])){plan$running.controller$data_trim[3] <- T}
  if(!as.logical(record.marker_last[4,2])){plan$running.controller$data_trim[4] <- T}
  
  if(!as.logical(record.marker_last[5,2])){plan$running.controller$peak_profiling[1] <- T}
  if(!as.logical(record.marker_last[6,2])){plan$running.controller$peak_profiling[2] <- T}
  if(!as.logical(record.marker_last[7,2])){plan$running.controller$peak_profiling[3] <- T}
  if(!as.logical(record.marker_last[8,2])){plan$running.controller$peak_profiling[4] <- T}
  
  if(!as.logical(record.marker_last[9,2])){plan$running.controller$others_1[1] <- T}
  if(!as.logical(record.marker_last[10,2])){plan$running.controller$others_1[2] <- T}
  if(!as.logical(record.marker_last[11,2])){plan$running.controller$others_1[3] <- T}
  if(!as.logical(record.marker_last[12,2])){plan$running.controller$others_1[4] <- T}
  
  if(!as.logical(record.marker_last[13,2])){plan$running.controller$operators[1] <- T}
  if(!as.logical(record.marker_last[14,2])){plan$running.controller$operators[2] <- F} # Switch off this operator at the end of plan define
  if(!as.logical(record.marker_last[15,2])){plan$running.controller$operators[3] <- F} # Switch off this operator at the end of plan define
  if(!as.logical(record.marker_last[16,2])){plan$running.controller$operators[4] <- T}
  if(!as.logical(record.marker_last[17,2])){plan$running.controller$operators[5] <- T}
  if(!as.logical(record.marker_last[18,2])){plan$running.controller$operators[6] <- T}
  if(!as.logical(record.marker_last[19,2])){plan$running.controller$operators[7] <- T}
  if(!as.logical(record.marker_last[20,2])){plan$running.controller$operators[8] <- T}
  
  return(plan)
  
}

recordMarker_resetter <- function(plan_count){
  
  ## Recording initial markers about being ran or not
  record.marker <- matrix(nrow = 20,ncol = 2);
  record.marker[,1] <- c("trim_1","trim_2","trim_3","trim_4","profiling_1","profiling_2","profiling_3","profiling_4",
                         "others_1","others_2","others_3","others_4","operators_1","operators_2","operators_3","operators_4",
                         "operators_5","operators_6","operators_7","operators_8");
  record.marker[,2] <- rep(F, 20)
  record.path <- paste0(getwd(),"/temp/records");
  
  saveRDS(record.marker,file = paste0(record.path,"/records_marker_",plan_count,".rds"));
  
}

marker_record <- function(functionNM){
  record.path <- paste0(getwd(),"/temp/records");
  
  if (!file.exists(paste0(record.path,"/records_marker_",plan_count,".rds"))){
    record.marker <- matrix(nrow = 20,ncol = 2);
    record.marker[,1] <- c("trim_1","trim_2","trim_3","trim_4","profiling_1","profiling_2","profiling_3","profiling_4",
                           "others_1","others_2","others_3","others_4","operators_1","operators_2","operators_3","operators_4",
                           "operators_5","operators_6","operators_7","operators_8");
    record.marker[,2] <- rep(F, 20)
    saveRDS(record.marker,file = paste0(record.path,"/records_marker_",plan_count,".rds"))
  } else {
    record.marker <- readRDS(paste0(record.path,"/records_marker_",plan_count,".rds"));
  };
  
  
  # If this step has been run at "plan_count" time, is will be marked as T
  if(functionNM=="datatrim_c1"){record.marker[1,2] <- T};
  if(functionNM=="datatrim_c2"){record.marker[2,2] <- T};
  if(functionNM=="datatrim_c3"){record.marker[3,2] <- T};
  if(functionNM=="datatrim_c4"){record.marker[4,2] <- T};
  
  if(functionNM=="peak_profiling_c1"){record.marker[5,2] <- T};
  if(functionNM=="peak_profiling_c2"){record.marker[6,2] <- T};
  if(functionNM=="peak_profiling_c3"){record.marker[7,2] <- T};
  if(functionNM=="peak_profiling_c4"){record.marker[8,2] <- T};
  
  if(functionNM=="optimized_results_c1"){record.marker[9,2] <- T};
  if(functionNM=="raw_data_samples_c2"){record.marker[10,2] <- T};
  if(functionNM=="raw_data_samples_c3"){record.marker[11,2] <- T};
  if(functionNM=="others_c4"){record.marker[12,2] <- T};
  
  if(functionNM=="operators_operators_1"){record.marker[13,2] <- T};
  if(functionNM==""){record.marker[14,2] <- F};
  if(functionNM==""){record.marker[15,2] <- F};
  if(functionNM=="peak_annotation_0"){record.marker[16,2] <- T};
  if(functionNM==""){record.marker[17,2] <- F};
  if(functionNM==""){record.marker[18,2] <- F};
  if(functionNM==""){record.marker[19,2] <- F};
  if(functionNM==""){record.marker[20,2] <- F};
  
  saveRDS(record.marker,file = paste0(record.path,"/records_marker_",plan_count,".rds"));
}

controller.resetter <- function() {
  
  points <- list(rep(T,4));
  names(points[[1]]) <- c("c1","c2","c3","c4");
  running.controller <- rep(points,3);
  
  operators <- c(F,F,F,F,F,F,F,F);
  names(operators) <- c("operators_1","operators_2","operators_3","operators_4",
                        "operators_5","operators_6","operators_7","operators_8");
  running.controller[[4]] <- operators;
  
  names(running.controller) <- c("data_trim","peak_profiling","others_1","operators");
  
  return(running.controller);
  
}

perform.plan <- function(plan.set){
  
  for (i in plan.set){
    perform.command(i)
  }
  
}

perform.command <- function(command){
  
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

cache.save <- function(obj,funpartnm){
  
  tmp_path <- paste0(getwd(),"/temp/cache");
  if (!dir.exists(tmp_path)){dir.create(tmp_path,recursive = T)};
  temp <- tempfile(tmpdir=tmp_path,fileext = ".rds");
  saveRDS(obj,file = temp);
  
  tmp_path_r <- paste0(getwd(),"/temp/records");
  if (!dir.exists(tmp_path_r)){dir.create(tmp_path,recursive = T)};
  if (!file.exists(paste0(getwd(),"/temp/records/file_record.rds"))){
    
  };
  temp_file_name <- basename(temp)
  info.save(funpartnm,tmp_path_r,temp_file_name)
}

info.save <- function(funpartnm,tmp_path_r,temp_file_name){
  
  info.matrix <- readRDS(paste0(tmp_path_r,"/records.rds"));
  if (identical(which(info.matrix[,1]==funpartnm),integer(0))){
    record.row <- which(is.na(info.matrix[,1]))[1];
  } else {
    record.row <- which(info.matrix[,1]==funpartnm)[1];
  }
  
  
  info.matrix[record.row,1] <- funpartnm;
  info.matrix[record.row,2] <- temp_file_name;
  
  saveRDS(info.matrix,file = paste0(tmp_path_r,"/records.rds"))
}

cache.read <- function(function.name, point){
  
  info.matrix <- readRDS(paste0(getwd(),"/temp/records/records.rds"));
  temp_point <- paste0(function.name,"_",point);
  temp_file <- info.matrix[match(temp_point,info.matrix[,1]),2];
  
  obj <- readRDS(paste0(getwd(),"/temp/cache/",temp_file));
  
  return(obj)
  
}

profiling_param_identifier <- function(new_command,last_command){
  
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

CreateRawRscript <- function(guestName, planString, planString2, rawfilenms.vec){
  
  guestName <<- guestName;
  planString <<- planString;
  
  if(dir.exists("/home/glassfish/payara5/glassfish/domains/")){
    
    scr.path <- "/home/glassfish/payara5/glassfish/domains/domain1/applications/MetaboAnalyst/resources/rscripts/metaboanalystr/"
    users.path <-paste0("/data/glassfish/projects/metaboanalyst//", guestName)
    
  } else if(dir.exists("/media/zzggyy/disk/")){
    
    users.path <-getwd();
    scr.path <-"/media/zzggyy/disk/MetaboAnalyst/target/MetaboAnalyst-5-18/resources/rscripts/metaboanalystr/"
    
  }else {
    
    users.path <-getwd();
    scr.path <-paste0(strsplit(users.path, "users")[[1]][1], "rscripts/metaboanalystr/")
    
  }
  
  ## Prepare the script for running
  str <- paste0('scripts.path <- ','"', scr.path,'"');
  str <- paste0(str, '\n', 'general_files <- c("general_data_utils.R","general_misc_utils.R","general_lib_utils.R","generic_c_utils.R");')
  str <- paste0(str, '\n', 'spectra_files <- c("spectra_generic_utils.R","data_trimming.R","parameters_db.R","parameters_optimization.R","preproc_utils.R","spectra_resume.R","spectra_utils.R");')
  str <- paste0(str, '\n', 'rawfilenms <<-', rawfilenms.vec)
  str <- paste0(str, '\n', 'err.vec <<-', '""')
  str <- paste0(str, '\n', '.on.public.web <<- TRUE')
  str <- paste0(str, '\n', 'file.sources = c(general_files, spectra_files)')
  str <- paste0(str, '\n', 'for (f in file.sources) {source(paste0(scripts.path, f))}')
  
  ## Construct the opt pipeline
  if(planString2 == "opt"){
    str <- paste0(str, '\n', 'plan <- InitializaPlan("raw_opt","' , users.path,  '/")')
    str <- paste0(str, '\n', 'data_folder_QC <- "',users.path , '/upload/QC/"')
    str <- paste0(str, '\n',  planString)
    str <- paste0(str, '\n',  "ExecutePlan(plan)")
  }
  
  ## Construct the default pipeline
  if(planString2 == "default"){
    str <- paste0(str, '\n', 'plan <- InitializaPlan("raw_ms","' , users.path,  '/")')
    str <- paste0(str, '\n',  planString)
    str <- paste0(str, '\n',  "ExecutePlan(plan)")
  }
  
  # sink command for running
  sink("ExecuteRawSpec.R");
  cat(str);
  sink();
  return(1)
  
}

CreateRawRscript2 <- function(guestName, meth){
  
  guestName <<- guestName;  
  
  if(dir.exists("/home/glassfish/payara5/glassfish/domains/")){
    
    scr.path <- "/home/glassfish/payara5/glassfish/domains/domain1/applications/MetaboAnalyst/resources/rscripts/metaboanalystr/"
    users.path <-paste0("/data/glassfish/projects/metaboanalyst/", guestName)
    
  } else {
    
    users.path <-getwd();
    scr.path <-paste0(strsplit(users.path, "users")[[1]][1], "rscripts/metaboanalystr/")
    
  }
  
  str <- paste0('scripts.path <- ','"', scr.path,'"');
  str <- paste0(str, '\n', 'general_files <- c("general_data_utils.R","general_misc_utils.R","general_lib_utils.R","generic_c_utils.R");')
  str <- paste0(str, '\n', 'spectra_files <- c("spectra_generic_utils.R","data_trimming.R","parameters_db.R","parameters_optimization.R","preproc_utils.R","spectra_resume.R","spectra_utils.R");')
  str <- paste0(str, '\n', '.on.public.web <<- TRUE')
  str <- paste0(str, '\n', 'file.sources = c(general_files, spectra_files)')
  str <- paste0(str, '\n', 'for (f in file.sources) {source(paste0(scripts.path, f))}')
  
  
  ## Prepare the script for running
  if(meth == "auto"){
    str <- paste0(str, '\n', 'FastRunningShow_auto ("',users.path,'")');
  } else {
    str <- paste0(str, '\n', 'FastRunningShow_customized ("',users.path,'")');
  }
  
  # sink command for running
  sink("ExecuteRawSpec.R");
  cat(str);
  sink();
  return(1)
  
}



FastRunningShow_auto <- function(fullUserPath){
  
  setwd(fullUserPath);  
  
  time_interval1 <- 0.65;
  time_interval2 <- 2;
  # running to show the progress
  if (!file.exists("metaboanalyst_spec_proc.txt")){
    print_mes <- paste0("Running Status -- Job Submitted Successfully at: ", Sys.time(), "\nWaiting in queue to start...")
    write.table(print_mes, file = paste0(fullUserPath, "/metaboanalyst_spec_proc.txt"),row.names = F,col.names = F);
  }
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 0/6: Scanning ROIs for parameters optimization...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(1.0, file = "log_progress.txt",row.names = F,col.names = F); 
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Data Loading...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  write.table("Raw file import begin...",file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Importing ","QC_PREFA02.mzML",":");    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = " ");
  
  Sys.sleep(time_interval1);
  
  for(i in seq(0,100,20)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "");
    Sys.sleep(0.35);
    
  }
  
  print_mes <- paste0(" Done!");    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Importing ","QC_PREFB02.mzML",":");    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = " ");
  
  Sys.sleep(time_interval1);
  
  for(i in seq(0,100,20)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "");
    Sys.sleep(0.3);
    
  }
  
  print_mes <- paste0(" Done!");    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(2.0, file = "log_progress.txt",row.names = F,col.names = F);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Data Loaded !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Empty Scan scanning...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "No Empty scan found !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying regions of interest (ROI)...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying regions of potential contaminants ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "");
  
  Sys.sleep(time_interval1);
  
  for(i in 1:30){
    print_mes <- ".";    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "");
    Sys.sleep(0.75);
  }
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "300 potential contaminamts will not be used for parameters optimization !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "MS data Preparing...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "MS Data ready !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(3.0, file = "log_progress.txt",row.names = F,col.names = F);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying ROIs in m/z dimensions...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying ROIs in m/z dimensions Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying ROIs in RT dimensions...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identification on ROIs Finished! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(4.0, file = "log_progress.txt",row.names = F,col.names = F);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Optimization will be started soon...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 1/6: Start to optimize parameters! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "This step may take a long time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "DoE Optimization Starting Now...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Evaluating Noise level...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(5.00, file = "log_progress.txt",row.names = F,col.names = F);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Preparing Parameters for optimization finished !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Data Spliting Finished ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Peak Preparing Begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Peak Preparing Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Round:1";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "DoE Running Begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Finished:";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = " ");
  
  for(i in c(11,22,33,44,56,67,78,89,100)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "");
    write.table(5+i/50, file = "log_progress.txt",row.names = F,col.names = F);
    Sys.sleep(2);
    
  }
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  print_mes <- "Round 1 Finished !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Round:2";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "DoE Running Begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Finished:";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = " ");
  
  for(i in c(11,22,33,44,56,67,78,89,100)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "");
    write.table(8.0+1/50, file = "log_progress.txt",row.names = F,col.names = F);
    Sys.sleep(0.5);
    
  }
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  print_mes <- "Round 2 Finished !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  
  print_mes <- "Round:3";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "DoE Running Begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Finished:";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = " ");
  
  for(i in c(11,22,33,44,56,67,78,89,100)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "");
    write.table(8.0+1/50, file = "log_progress.txt",row.names = F,col.names = F);
    Sys.sleep(0.5);
    
  }
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  print_mes <- "Round 3 Finished !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 1/6: Parameters Optimization Finished ! (", Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(20.0, file = "log_progress.txt",row.names = F,col.names = F);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 2/6: Start to import the spectrum! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(21.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "This step will take a short time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Raw file import begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-6KUCT.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-77FXR.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-9OS5Y.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-9WOBP.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-9SNJ4.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-9X47O.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-AMR37.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-AUP8B.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "QC_PREFA02.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "QC_PREFB02.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(22.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(3);
  
  print_mes <- "Raw file initialized Successfully!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(24.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 2/6: Successfully imported raw MS data! (", Sys.time(),") \nGoing to the next step...");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Parameters for centWave have been successfully parsed!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "4 CPU Threads will be used for peak profiling !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(25, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 3/6: Started peak picking! This step will take some time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 3613 regions of interest ... OK: 551 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(30.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 3528 regions of interest ... OK: 591 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(30.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 5098 regions of interest ... OK: 598 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(35.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 3953 regions of interest ... OK: 805 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(35.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 4875 regions of interest ... OK: 676 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 5333 regions of interest ... OK: 717 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 4991 regions of interest ... OK: 805 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 4633 regions of interest ... OK: 847 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 4988 regions of interest ... OK: 791 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(47.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 5069 regions of interest ... OK: 755 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(47.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 3/6: Peak picking finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(50.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 4/6: Started peak alignment! This step is running...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(52.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Total of 2803 slices detected for processing... Done ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(55.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the next step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(56.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Retention time correction is running.....";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(60.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Performing retention time correction using  92  peak groups.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(65.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Applying retention time adjustment to the identified chromatographic peaks ... Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(66.0, file = "log_progress.txt",row.names = F,col.names = F);
  write.table(68.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Total of 2803 slices detected for processing... Done ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the next step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 4/6: Peak alignment finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(73.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 5/6: Started peak filling! This step may take some time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Starting peak filling!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(74.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Defining peak areas for filling-in...... OK";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Start integrating peak areas from original files...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(77.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 296 peaks from CD-6KUCT.mzML ...  got  150 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(78.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 420 peaks from CD-9WOBP.mzML ...  got  146 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(79.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 380 peaks from HC-9X47O.mzML ...  got  192 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(80.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 371 peaks from HC-AUP8B.mzML ...  got  197 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(81.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 35 peaks from QC_PREFA02.mzML ...  got  32 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 275 peaks from CD-77FXR.mzML ...  got  148 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(82.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 33 peaks from QC_PREFB02.mzML ...  got  30 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(83.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 358 peaks from HC-AMR37.mzML ...  got  173 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(84.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 356 peaks from HC-9SNJ4.mzML ...  got  189 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(85.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 283 peaks from CD-9OS5Y.mzML ...  got  160 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(86.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 5/6: Peak filing finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(88.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Peak Profiling finished successfully !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Begin to plotting figures...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(89.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 6/6: Starting Peak Annotation...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(91.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Start grouping after retention time.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Created 68 pseudospectra.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(92.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating peak matrix...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Run isotope peak annotation..";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(93.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Found isotopes:108";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Start grouping after correlation...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(94.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating EIC's....";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-6KUCT.mzML  ... 37 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-77FXR.mzML  ... 43 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-9OS5Y.mzML  ... 54 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-9WOBP.mzML  ... 124 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-9SNJ4.mzML  ... 28 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(95.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-9X47O.mzML  ... 0 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  print_mes <- "Detecting  HC-AMR37.mzML  ... 170 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-AUP8B.mzML  ... 20 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  QC_PREFA02.mzML  ... 64 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  QC_PREFB02.mzML  ... 87 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Calculating graph cross linking in 68 Groups...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(96.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "New number of ps-groups:  406";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "mSet has now 406 groups, instead of 68 !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating peak matrix for peak annotation!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Polarity is set in annotaParam: negative";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(97.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Ruleset could not read from object! Recalculating...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 6/6: Successfully performed peak annotation! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(98.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the final step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(99.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Everything has been finished Successfully ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(100.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
}

FastRunningShow_customized <- function(fullUserPath){
  print(fullUserPath);
  setwd(fullUserPath);    
  
  time_interval1 <- 0.80;
  time_interval2 <- 2;
  # running to show the progress
  if (!file.exists("metaboanalyst_spec_proc.txt")){
    print_mes <- paste0("Running Status -- Job Submitted Successfully at: ", Sys.time(), "\nWaiting in queue to start...")
    write.table(print_mes, file = paste0(fullUserPath, "/metaboanalyst_spec_proc.txt"),row.names = F,col.names = F);
  }
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 1/6: Internalize parameters!");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(1.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("This step will be finished soon...");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(2.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 1/6: Parameters Internalized Successfully! ");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(15.0, file = "log_progress.txt",row.names = F,col.names = F);
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Going to the next step...");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(20.0, file = "log_progress.txt",row.names = F,col.names = F);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 2/6: Start to import the spectrum! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(21.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "This step will take a short time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Raw file import begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-6KUCT.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-77FXR.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-9OS5Y.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-9WOBP.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-9SNJ4.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-9X47O.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-AMR37.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-AUP8B.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "QC_PREFA02.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "QC_PREFB02.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(22.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(3);
  
  print_mes <- "Raw file initialized Successfully!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(24.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 2/6: Successfully imported raw MS data! (", Sys.time(),") \nGoing to the next step...");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Parameters for centWave have been successfully parsed!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "4 CPU Threads will be used for peak profiling !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(25, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 3/6: Started peak picking! This step will take some time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5509 regions of interest ... OK: 1285 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(30.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5355 regions of interest ... OK: 1395 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(30.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6009 regions of interest ... OK: 1327 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(35.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5767 regions of interest ... OK: 1810 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(35.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5832 regions of interest ... OK: 1507 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6273 regions of interest ... OK: 1586 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6061 regions of interest ... OK: 1643 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5926 regions of interest ... OK: 1805 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6054 regions of interest ... OK: 1646 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(47.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6299 regions of interest ... OK: 1608 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(47.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 3/6: Peak picking finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(50.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 4/6: Started peak alignment! This step is running...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(52.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Total of 2805 slices detected for processing... Done ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(55.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the next step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(56.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Retention time correction is running.....";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(60.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Performing retention time correction using  54  peak groups.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(65.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Applying retention time adjustment to the identified chromatographic peaks ... Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(66.0, file = "log_progress.txt",row.names = F,col.names = F);
  write.table(68.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Total of 2805 slices detected for processing... Done ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the next step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 4/6: Peak alignment finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(73.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 5/6: Started peak filling! This step may take some time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Starting peak filling!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(74.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Defining peak areas for filling-in...... OK";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Start integrating peak areas from original files...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(77.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 139 peaks from CD-6KUCT.mzML ...  got  56 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(78.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 256 peaks from CD-9WOBP.mzML ...  got  55 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(79.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 261 peaks from HC-9X47O.mzML ...  got  101 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(80.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 234 peaks from HC-AUP8B.mzML ...  got  95 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(81.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 199 peaks from HC-AMR37.mzML ...  got  72 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 206 peaks from HC-9SNJ4.mzML ...  got  94 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(82.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 146 peaks from CD-77FXR.mzML ...  got  61 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(83.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 52 peaks from QC_PREFA02.mzML ...  got  30 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(84.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 188 peaks from CD-9OS5Y.mzML ...  got  74 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(85.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 54 peaks from QC_PREFB02.mzML ...  got  37 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(86.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 5/6: Peak filing finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(88.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Peak Profiling finished successfully !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Begin to plotting figures...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(89.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 6/6: Starting Peak Annotation...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(91.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Start grouping after retention time.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Created 62 pseudospectra.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(92.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating peak matrix...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Run isotope peak annotation..";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(93.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Found isotopes:28";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Start grouping after correlation...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(94.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating EIC's....";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-6KUCT.mzML  ... 66 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-77FXR.mzML  ... 35 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-9OS5Y.mzML  ... 219 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-9WOBP.mzML  ... 120 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-9SNJ4.mzML  ... 30 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(95.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-9X47O.mzML  ... 0 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(95.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-AMR37.mzML  ... 287 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-AUP8B.mzML  ... 21 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  QC_PREFA02.mzML  ... 89 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  QC_PREFB02.mzML  ... 39 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Calculating graph cross linking in 62 Groups...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(96.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "New number of ps-groups:  736";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "mSet has now 310 groups, instead of 62 !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating peak matrix for peak annotation!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Polarity is set in annotaParam: negative";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(97.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Ruleset could not read from object! Recalculating...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 6/6: Successfully performed peak annotation! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(98.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the final step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(99.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Everything has been finished Successfully ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  write.table(100.0, file = "log_progress.txt",row.names = F,col.names = F);
  Sys.sleep(time_interval1);
}


#' @param folderPath guest folder
#'
#' @title Cache Update
#' @author Zhiqiang Pang
#' @noRd
CachePathCorrection <- function(folderPath){
  
  cacheFiles <- list.files(paste0(folderPath,"/temp"),all.files = T, full.names = T, recursive = T);
  
  newfilepath <- list.files("/home/glassfish/projects/MetaboDemoRawData/upload",all.files = T, recursive = T,full.names = T)
  
  for (i in cacheFiles){
    tmp_rds <- readRDS(i);
    tmp_names <- names(tmp_rds);
    
    if(class(tmp_rds)[[1]] == "OnDiskMSnExp"){
      
      tmp_rds@processingData@files <- newfilepath;
      
    } else if (class(tmp_rds)[[1]] == "list"){
      
      for(j in tmp_names){
        tmp_obj <- tmp_rds[[j]];
        
        if(class(tmp_obj)[[1]] == "OnDiskMSnExp"){
          tmp_rds[[j]]@processingData@files <- newfilepath;
        }
        
        if(class(tmp_obj)[[1]] == "xcmsSet"){
          tmp_rds[[j]]@filepaths <- newfilepath;
        }
        
      } 
      
    } else if(class(tmp_rds)[[1]] == "environment"){
      
      tmp_rds$annotPeaks$onDiskData@processingData@files <- 
        tmp_rds$annotPeaks$xcmsSet@filepaths <- 
        tmp_rds$mSet$onDiskData@processingData@files <- 
        tmp_rds$mSet$xcmsSet@filepaths <- 
        tmp_rds$rawData@processingData@files <- newfilepath;
      
    }
    
    saveRDS(tmp_rds, file = i);
  }
  
  load(paste0(folderPath,"/","mSet.rda"));
  mSet$onDiskData@processingData@files <- mSet$xcmsSet@filepaths <- newfilepath;
  
  save(mSet, file = paste0(folderPath,"/","mSet.rda"));
  
}

#folderPath <- "/home/qiang/NetBeansProjects/MetaboAnalyst/target/MetaboAnalyst-5.18/resources/users/guest8386191689266673959tmp/"
#CachePathCorrection(folderPath)


## End of the resuming script ---
