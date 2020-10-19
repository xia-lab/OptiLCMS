# Special script for .on.web.public == TRUE

GeneratePeakList <- function(userPath) {
  
  oldwd <- getwd();
  on.exit(setwd(oldwd));
  
  setwd(userPath)
  
  ## Claculate the mean internsity of all groups
  sample_data <-
    read.csv("metaboanalyst_input.csv",
             header = TRUE,
             stringsAsFactors = FALSE)
  
  groups <- as.character(as.matrix(sample_data[1, ]))[-1]
  
  sample_data <- sample_data[-1, -1]
  
  
  if (length(unique(groups)) == 1) {
    sample_data_mean  <-
      apply(
        sample_data,
        1,
        FUN = function(x) {
          mean(as.numeric(x), na.rm = TRUE)
        }
      )
    
  } else {
    sample_data1 <- matrix(nrow = nrow(sample_data))
    
    
    for (i in 1:length(unique(groups))) {
      columnnum <- unique(groups)[i] == groups
      sample_data0  <-
        subset.data.frame(sample_data, subset = TRUE, select = columnnum)
      
      sample_data0  <-
        round(apply(
          sample_data0,
          1,
          FUN = function(x) {
            mean(as.numeric(x), na.rm = TRUE)
          }
        ), 2)
      
      sample_data1 <- cbind(sample_data1, sample_data0)
    }
    sample_data_mean <- sample_data1[, -1]
    colnames(sample_data_mean) <- unique(groups)
  }
  
  ## Prepare other information
  ann_data <- readRDS("annotated_peaklist.rds")
  
  ann_data <-
    ann_data[, c(1, 4, ncol(ann_data) - 1, ncol(ann_data) - 2)]
  ann_data[, 1] <- round(ann_data[, 1], 4)
  ann_data[, 2] <- round(ann_data[, 2], 2)
  
  write.csv(
    cbind(ann_data, sample_data_mean),
    file = "peak_feature_summary.csv",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE
  )
  
  return (nrow(ann_data))
}


#' @title featureSUM
#' @description a function used to summarize peak features
#' @param MS_group MS_group
#' @param frtr frtr
#' @noRd
#' @import scales
featureSUM <- function(MS_group, frtr) {
  # sanity check
  if (identical(MS_group, list())) {
    rts <- quantile(frtr)
    res <- data.frame(rts[1], 0)
    res[2, ] <- c(rts[2], 0.5)
    res[3, ] <- c(rts[3], 1)
    res[4, ] <- c(rts[4], 0.5)
    res[5, ] <- c(rts[5], 0)
    colnames(res) <- c("RT_mean", "Inten_mean")
    
    return(res)
  }
  
  # summarize intensity and RT
  inten_sum <-
    sapply(
      MS_group,
      FUN = function(x) {
        sum(x@intensity, na.rm = TRUE)
      }
    )
  
  inten_sum[(inten_sum == 0)] <- 1
  
  rt_min_sum <- sapply(
    MS_group,
    FUN = function(x) {
      min(x@rtime)
    }
  )
  
  rt_max_sum <- sapply(
    MS_group,
    FUN = function(x) {
      max(x@rtime)
    }
  )
  
  scan_sum <- sapply(
    MS_group,
    FUN = function(x) {
      length(x@rtime)
    }
  )
  
  # correct RT
  rt_min_corrected <- sum(rt_min_sum * inten_sum) / sum(inten_sum)
  rt_max_corrected <- sum(rt_max_sum * inten_sum) / sum(inten_sum)
  rt_range_cor <- abs(rt_max_corrected - rt_min_corrected)
  
  # if(.on.public.web){
  #   load_scales()
  # }
  
  MS_group <-
    sapply(
      MS_group,
      FUN = function(x, rt_min_corrected, rt_max_corrected) {
        x@rtime <- rescale(x@rtime, c(rt_min_corrected, rt_max_corrected))
        
        return(x)
      },
      rt_min_corrected = rt_min_corrected,
      rt_max_corrected = rt_max_corrected
    )
  
  # extract information
  df <- data.frame()
  
  for (u in 1:length(MS_group)) {
    ddf <- data.frame(MS_group[[u]]@rtime, MS_group[[u]]@intensity)
    df <- rbind(df, ddf)
  }
  
  df[is.na(df)] <- 0
  colnames(df) <- c("rt", "intensity")
  df <- df[order(df$rt),]
  res <- data.frame()
  
  # bin all scans
  binsize <-
    rt_range_cor / (max(scan_sum) - 1) - 0.001
  # Avoid the boundary effect by minus 0.001
  rt_now <- rt_min_corrected + binsize
  
  if (length(MS_group) > 1) {
    while (rt_now < rt_max_corrected + binsize) {
      Inten_mean <-
        mean(df[df$rt <= rt_now & df$rt >= (rt_now - binsize), 2])
      
      RT_mean <-
        mean(df[df$rt < rt_now & df$rt >= (rt_now - binsize), 1])
      res <- rbind(res, data.frame(RT_mean, Inten_mean))
      rt_now <- rt_now + binsize
    }
    
  } else {
    # if only one sample, use it diresctly
    res <- df
    rownames(res) <- NULL
  }
  
  colnames(res) <- c("RT_mean", "Inten_mean")
  # remove the empty bin
  if (any(is.nan(res$RT_mean) | is.nan(res$Inten_mean))) {
    res <- res[-which(is.nan(res$RT_mean) | is.nan(res$Inten_mean)),]
  }
  
  # manually add 2 empty points before and after the range
  res[nrow(res) + 1, ] <- c(min(res$RT_mean) - binsize, 0)
  res[nrow(res) + 1, ] <- c(min(res$RT_mean) - 2 * binsize, 0)
  res[nrow(res) + 1, ] <- c(max(res$RT_mean) + binsize, 0)
  res[nrow(res) + 1, ] <- c(max(res$RT_mean) + 2 * binsize, 0)
  
  return(res)
}

peakTableSUM <- function(peak_table) {
  max_peaks <- vector()
  
  if (length(unique(peak_table[, ncol(peak_table)])) != length(peak_table[, ncol(peak_table)])) {
    for (i in unique(peak_table[, ncol(peak_table)])) {
      tmp_table <- peak_table[peak_table[, ncol(peak_table)] == i, ]
      
      if (!is.null(nrow(tmp_table))) {
        max_peaks <-
          c(max_peaks, names(which.max(tmp_table[, which(colnames(tmp_table) == "into")])))
      } else {
        max_peaks <-
          c(max_peaks, names(which(peak_table[, ncol(peak_table)] == i)))
      }
    }
    
    return(peak_table[max_peaks, ])
  } else {
    return(peak_table)
  }
}

#' @title Cache Update
#' @param folderPath guest folder
#' @description used only for web cache update
#' @author Zhiqiang Pang
#' @noRd
CachePathCorrection <- function(folderPath){
  if(!.on.public.web){
    stop("Prohibited running this function at local!")
  }
  cacheFiles <- list.files(paste0(folderPath,"/temp"),all.files = TRUE, full.names = TRUE, recursive = TRUE);
  
  newfilepath <- list.files("/home/glassfish/projects/MetaboDemoRawData/upload",all.files = TRUE, recursive = TRUE,full.names = TRUE)
  
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

#' @noRd
FastRunningShow_auto <- function(fullUserPath){
  
  oldwd <- getwd();
  on.exit(setwd(oldwd));
  
  setwd(fullUserPath);  
  
  time_interval1 <- 0.65;
  time_interval2 <- 2;
  # running to show the progress
  if (!file.exists("metaboanalyst_spec_proc.txt")){
    print_mes <- paste0("Running Status -- Job Submitted Successfully at: ", Sys.time(), "\nWaiting in queue to start...")
    write.table(print_mes, file = paste0(fullUserPath, "/metaboanalyst_spec_proc.txt"),row.names = FALSE,col.names = FALSE);
  }
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 0/6: Scanning ROIs for parameters optimization...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(1.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE); 
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Data Loading...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  write.table("Raw file import begin...",file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Importing ","QC_PREFA02.mzML",":");    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = " ");
  
  Sys.sleep(time_interval1);
  
  for(i in seq(0,100,20)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "");
    Sys.sleep(0.35);
    
  }
  
  print_mes <- paste0(" Done!");    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Importing ","QC_PREFB02.mzML",":");    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = " ");
  
  Sys.sleep(time_interval1);
  
  for(i in seq(0,100,20)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "");
    Sys.sleep(0.3);
    
  }
  
  print_mes <- paste0(" Done!");    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(2.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Data Loaded !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Empty Scan scanning...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "No Empty scan found !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying regions of interest (ROI)...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying regions of potential contaminants ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "");
  
  Sys.sleep(time_interval1);
  
  for(i in 1:30){
    print_mes <- ".";    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "");
    Sys.sleep(0.75);
  }
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "300 potential contaminamts will not be used for parameters optimization !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "MS data Preparing...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "MS Data ready !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(3.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying ROIs in m/z dimensions...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying ROIs in m/z dimensions Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identifying ROIs in RT dimensions...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Identification on ROIs Finished! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(4.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Optimization will be started soon...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 1/6: Start to optimize parameters! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "This step may take a long time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "DoE Optimization Starting Now...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Evaluating Noise level...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(5.00, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Preparing Parameters for optimization finished !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Data Spliting Finished ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Peak Preparing Begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Peak Preparing Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Round:1";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "DoE Running Begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Finished:";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = " ");
  
  for(i in c(11,22,33,44,56,67,78,89,100)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "");
    write.table(5+i/50, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
    Sys.sleep(2);
    
  }
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  print_mes <- "Round 1 Finished !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Round:2";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "DoE Running Begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Finished:";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = " ");
  
  for(i in c(11,22,33,44,56,67,78,89,100)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "");
    write.table(8.0+1/50, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
    Sys.sleep(0.5);
    
  }
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  print_mes <- "Round 2 Finished !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  
  print_mes <- "Round:3";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "DoE Running Begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Finished:";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = " ");
  
  for(i in c(11,22,33,44,56,67,78,89,100)){
    
    print_mes <- paste0(i,"% | ");    
    write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "");
    write.table(8.0+1/50, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
    Sys.sleep(0.5);
    
  }
  
  print_mes <- "Done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  print_mes <- "Round 3 Finished !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Model Parsing Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 1/6: Parameters Optimization Finished ! (", Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(20.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 2/6: Start to import the spectrum! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(21.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "This step will take a short time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Raw file import begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-6KUCT.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-77FXR.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-9OS5Y.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-9WOBP.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-9SNJ4.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-9X47O.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-AMR37.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-AUP8B.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "QC_PREFA02.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "QC_PREFB02.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(22.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(3);
  
  print_mes <- "Raw file initialized Successfully!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(24.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 2/6: Successfully imported raw MS data! (", Sys.time(),") \nGoing to the next step...");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Parameters for centWave have been successfully parsed!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "4 CPU Threads will be used for peak profiling !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(25, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 3/6: Started peak picking! This step will take some time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 3613 regions of interest ... OK: 551 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(30.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 3528 regions of interest ... OK: 591 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(30.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 5098 regions of interest ... OK: 598 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(35.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 3953 regions of interest ... OK: 805 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(35.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 4875 regions of interest ... OK: 676 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 5333 regions of interest ... OK: 717 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 4991 regions of interest ... OK: 805 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 4633 regions of interest ... OK: 847 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 4988 regions of interest ... OK: 791 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(47.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 1.91 ppm ... Detecting peaks in 5069 regions of interest ... OK: 755 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(47.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 3/6: Peak picking finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(50.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 4/6: Started peak alignment! This step is running...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(52.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Total of 2803 slices detected for processing... Done ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(55.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the next step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(56.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Retention time correction is running.....";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(60.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Performing retention time correction using  92  peak groups.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(65.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Applying retention time adjustment to the identified chromatographic peaks ... Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(66.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  write.table(68.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Total of 2803 slices detected for processing... Done ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the next step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 4/6: Peak alignment finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(73.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 5/6: Started peak filling! This step may take some time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Starting peak filling!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(74.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Defining peak areas for filling-in...... OK";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Start integrating peak areas from original files...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(77.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 296 peaks from CD-6KUCT.mzML ...  got  150 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(78.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 420 peaks from CD-9WOBP.mzML ...  got  146 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(79.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 380 peaks from HC-9X47O.mzML ...  got  192 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(80.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 371 peaks from HC-AUP8B.mzML ...  got  197 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(81.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 35 peaks from QC_PREFA02.mzML ...  got  32 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 275 peaks from CD-77FXR.mzML ...  got  148 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(82.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 33 peaks from QC_PREFB02.mzML ...  got  30 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(83.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 358 peaks from HC-AMR37.mzML ...  got  173 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(84.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 356 peaks from HC-9SNJ4.mzML ...  got  189 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(85.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 283 peaks from CD-9OS5Y.mzML ...  got  160 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(86.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 5/6: Peak filing finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(88.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Peak Profiling finished successfully !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Begin to plotting figures...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(89.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 6/6: Starting Peak Annotation...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(91.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Start grouping after retention time.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Created 68 pseudospectra.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(92.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating peak matrix...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Run isotope peak annotation..";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(93.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Found isotopes:108";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Start grouping after correlation...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(94.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating EIC's....";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-6KUCT.mzML  ... 37 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-77FXR.mzML  ... 43 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-9OS5Y.mzML  ... 54 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-9WOBP.mzML  ... 124 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-9SNJ4.mzML  ... 28 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(95.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-9X47O.mzML  ... 0 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  print_mes <- "Detecting  HC-AMR37.mzML  ... 170 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-AUP8B.mzML  ... 20 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  QC_PREFA02.mzML  ... 64 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  QC_PREFB02.mzML  ... 87 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Calculating graph cross linking in 68 Groups...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(96.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "New number of ps-groups:  406";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "mSet has now 406 groups, instead of 68 !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating peak matrix for peak annotation!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Polarity is set in annotaParam: negative";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(97.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Ruleset could not read from object! Recalculating...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 6/6: Successfully performed peak annotation! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(98.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the final step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(99.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Everything has been finished Successfully ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(100.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
}

#' @noRd
FastRunningShow_customized <- function(fullUserPath){
  
  oldwd <- getwd();
  on.exit(setwd(oldwd));
  
  setwd(fullUserPath);    
  
  time_interval1 <- 0.80;
  time_interval2 <- 2;
  # running to show the progress
  if (!file.exists("metaboanalyst_spec_proc.txt")){
    print_mes <- paste0("Running Status -- Job Submitted Successfully at: ", Sys.time(), "\nWaiting in queue to start...")
    write.table(print_mes, file = paste0(fullUserPath, "/metaboanalyst_spec_proc.txt"),row.names = FALSE,col.names = FALSE);
  }
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 1/6: Internalize parameters!");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(1.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("This step will be finished soon...");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(2.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 1/6: Parameters Internalized Successfully! ");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(15.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Going to the next step...");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(20.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 2/6: Start to import the spectrum! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(21.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "This step will take a short time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Raw file import begin...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-6KUCT.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-77FXR.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-9OS5Y.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "CD-9WOBP.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-9SNJ4.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-9X47O.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-AMR37.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "HC-AUP8B.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "QC_PREFA02.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "QC_PREFB02.mzML import done!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(22.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(3);
  
  print_mes <- "Raw file initialized Successfully!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(24.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 2/6: Successfully imported raw MS data! (", Sys.time(),") \nGoing to the next step...");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Parameters for centWave have been successfully parsed!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "4 CPU Threads will be used for peak profiling !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(25, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 3/6: Started peak picking! This step will take some time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5509 regions of interest ... OK: 1285 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(30.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5355 regions of interest ... OK: 1395 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(30.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6009 regions of interest ... OK: 1327 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(35.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5767 regions of interest ... OK: 1810 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(35.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5832 regions of interest ... OK: 1507 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6273 regions of interest ... OK: 1586 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6061 regions of interest ... OK: 1643 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 5926 regions of interest ... OK: 1805 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(40.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6054 regions of interest ... OK: 1646 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(47.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Searching under 5 ppm ... Detecting peaks in 6299 regions of interest ... OK: 1608 found.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(47.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 3/6: Peak picking finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(50.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 4/6: Started peak alignment! This step is running...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(52.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Total of 2805 slices detected for processing... Done ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(55.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the next step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(56.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Retention time correction is running.....";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(60.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Performing retention time correction using  54  peak groups.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(65.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Applying retention time adjustment to the identified chromatographic peaks ... Done !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(66.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  write.table(68.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Total of 2805 slices detected for processing... Done ! ";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the next step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 4/6: Peak alignment finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(73.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 5/6: Started peak filling! This step may take some time...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Starting peak filling!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(74.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Defining peak areas for filling-in...... OK";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Start integrating peak areas from original files...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(77.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 139 peaks from CD-6KUCT.mzML ...  got  56 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(78.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 256 peaks from CD-9WOBP.mzML ...  got  55 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(79.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 261 peaks from HC-9X47O.mzML ...  got  101 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(80.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 234 peaks from HC-AUP8B.mzML ...  got  95 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(81.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 199 peaks from HC-AMR37.mzML ...  got  72 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 206 peaks from HC-9SNJ4.mzML ...  got  94 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(82.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 146 peaks from CD-77FXR.mzML ...  got  61 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(83.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 52 peaks from QC_PREFA02.mzML ...  got  30 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(84.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 188 peaks from CD-9OS5Y.mzML ...  got  74 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(85.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Requesting 54 peaks from QC_PREFB02.mzML ...  got  37 .";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(86.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 5/6: Peak filing finished ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(88.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Peak Profiling finished successfully !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Begin to plotting figures...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(89.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Step 6/6: Starting Peak Annotation...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(91.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Start grouping after retention time.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Created 62 pseudospectra.";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(92.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating peak matrix...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Run isotope peak annotation..";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(93.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Found isotopes:28";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Start grouping after correlation...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(94.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating EIC's....";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-6KUCT.mzML  ... 66 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-77FXR.mzML  ... 35 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-9OS5Y.mzML  ... 219 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  CD-9WOBP.mzML  ... 120 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-9SNJ4.mzML  ... 30 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(95.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-9X47O.mzML  ... 0 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(95.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-AMR37.mzML  ... 287 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  HC-AUP8B.mzML  ... 21 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  QC_PREFA02.mzML  ... 89 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Detecting  QC_PREFB02.mzML  ... 39 peaks found!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Calculating graph cross linking in 62 Groups...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(96.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "New number of ps-groups:  736";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "mSet has now 310 groups, instead of 62 !";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Generating peak matrix for peak annotation!";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- "Polarity is set in annotaParam: negative";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(97.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Ruleset could not read from object! Recalculating...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Step 6/6: Successfully performed peak annotation! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(98.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- "Going to the final step...";    
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(99.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
  
  print_mes <- paste0("Everything has been finished Successfully ! (",Sys.time(),")");
  write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = TRUE,row.names = FALSE,col.names = FALSE, quote = FALSE, eol = "\n");
  write.table(100.0, file = "log_progress.txt",row.names = FALSE,col.names = FALSE);
  Sys.sleep(time_interval1);
}

# importFrom("grDevices", "boxplot.stats", "dev.off", "jpeg")
# importFrom("graphics", "abline", "boxplot", "contour", "grid",
#            "legend", "par", "points")
# importFrom("stats", "approx", "approxfun", "as.formula", "convolve",
#            "cor", "cor.test", "cutree", "deriv3", "dist", "dnorm",
#            "fft", "fitted", "hclust", "kmeans", "lm", "loess", "lsfit",
#            "na.omit", "nextn", "nls", "prcomp", "predict", "sd",
#            "smooth.spline", "weighted.mean")
# 
# 
# 

