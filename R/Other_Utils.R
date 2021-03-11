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
  
  rt_min_sum <- rt_min_sum[is.finite(rt_min_sum)]
  rt_max_sum <- rt_max_sum[is.finite(rt_max_sum)]
  
  # correct RT
  rt_min_corrected <- sum(rt_min_sum * inten_sum) / sum(inten_sum)
  rt_max_corrected <- sum(rt_max_sum * inten_sum) / sum(inten_sum)
  rt_range_cor <- abs(rt_max_corrected - rt_min_corrected)

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

PeakGroupCV <- function(IntoLists, groupsInfo){
  allGroups <- names(groupsInfo)
  if(is.null(allGroups)){
    newGrIn <- list();
    allGroups <- colnames(groupsInfo);
    for(j in seq_along(allGroups)){
      newGrIn[[j]] <- as.vector(groupsInfo[,j]);
    }
    names(newGrIn) <- allGroups;
    groupsInfo <- newGrIn;
  }
  
  groups <- names(IntoLists);
  
  statsdf <- data.frame();
  count <- 0;
  
  for(i in allGroups){
    if(i %in% groups){
      count <- count +1;
      
      if(length(IntoLists[[i]]) == length(groupsInfo[[i]])){
        cv<- sd(as.numeric(IntoLists[[i]]))/mean(as.numeric(IntoLists[[i]]));
        statsdf[count, c(1,2)] <- c(i, round(cv, 2)) 
      } else {
        nmissing <- abs(length(IntoLists[[i]]) - length(groupsInfo[[i]]));
        ints <- c(as.numeric(IntoLists[[i]]), rep(0, nmissing));
        
        cv<- sd(ints)/mean(ints);
        statsdf[count, c(1,2)] <- c(i, round(cv, 2));
      }
    }
  }
  
  statsdf$V2 <- as.numeric(statsdf$V2)
  statsdf$V2[statsdf$V2 == 0] <- 0.05
  
  return(statsdf)
}

GetRGBColorGradient <- function(vals) {
  library(RColorBrewer)
  
  #seed.cols <- brewer.pal(3, "YlOrRd")
  #seed.cols <- brewer.pal(9, "Oranges")[c(2,5,7)]
  seed.cols <- c("#FCF5DF", "#FFEDA0", "#F03B20")
  cols <- colorRampPalette(seed.cols)(length(vals))
  # set alpha for
  my.alpha <-
    signif(seq(
      from = 0.3,
      to = 0.8,
      length.out = length(vals)
    ), 2)
  
  rgb.cols <- my.col2rgba(cols, alpha = my.alpha)
  
  # now need make sure values and colors are matched using names
  nms.orig <- names(vals)
  names(rgb.cols) <- names(sort(vals))
  ord.cols <- rgb.cols[nms.orig]
  return(as.vector(ord.cols))
  # note remove names
}

my.col2rgba <- function(cols, alpha) {
  rgbcols <- col2rgb(cols)
  rgbcols <- rbind(rgbcols, alpha)
  return(as.vector(apply(rgbcols, 2, function(x) {
    paste("rgba(", paste(x, collapse = ","), ")", sep = "")
  })))
}

GetDist3D <-function(mat, target=c(0,0,0)){
  dist.vec <- apply(mat, 2, function(x) dist(rbind(x, target)));
  return(dist.vec);
}

#' @title Cache Update
#' @param folderPath guest folder
#' @description used only for web cache update
#' @author Zhiqiang Pang
#' @noRd
CachePathCorrection <- function(folderPath){

  cacheFiles <- list.files(paste0(folderPath,"/specTemp"),all.files = TRUE, full.names = TRUE, recursive = TRUE);
  
  newfilepath <- list.files("/home/glassfish/projects/MetaboDemoRawData/upload",all.files = TRUE, recursive = TRUE,full.names = TRUE)
  
  for (i in cacheFiles){
    tmp_rds <- readRDS(i);
    tmp_names <- names(tmp_rds);
    
    if(class(tmp_rds)[1] == "mSet"){
      
      tmp_rds@rawOnDisk@processingData@files <- newfilepath;
      tmp_rds@rawInMemory@processingData@files <- newfilepath;
      tmp_rds@rawfiles <- newfilepath;
      
    } else if (class(tmp_rds)[1] == "OnDiskMSnExp"){
      
      tmp_rds@processingData@files <- newfilepath;
      
    } else if(class(tmp_rds)[1] == "environment"){
      tmp_rds[["mSet"]]@rawOnDisk@processingData@files <- 
        tmp_rds[["mSet"]]@rawInMemory@processingData@files <- 
        tmp_rds[["mSet"]]@rawfiles <- newfilepath;
    }
    
    saveRDS(tmp_rds, file = i);
  }
  
  load(paste0(folderPath,"/","mSet.rda"));
  mSet@rawOnDisk@processingData@files <- mSet@rawInMemory@processingData@files <- mSet@rawfiles <- newfilepath;
  save(mSet, file = paste0(folderPath,"/","mSet.rda"));
  
  load(paste0(folderPath,"/","raw_data_filt.rda"));
  raw_data_filt@processingData@files <- newfilepath;
  save(raw_data_filt, file = paste0(folderPath,"/","raw_data_filt.rda"));
  
}

 # folderPath <- "/home/qiang/NetBeansProjects/MetaboAnalyst/src/main/webapp/resources/data/spectra_example_cache/auto/"
 # CachePathCorrection(folderPath)

#' @noRd
FastRunningShow_auto <- function(fullUserPath){
  
  
  oldwd <- getwd();
  on.exit(setwd(oldwd));
  
  setwd(fullUserPath);    
  
  time_interval1 <- 0.80;
  time_interval2 <- 2;
  # running to show the progress
  MessageOutput(paste0("Running Status -- Plan Initialized Successfully at: ", 
                       Sys.time(), 
                       "\nCurrent OptiLCMS version is ",
                       packageVersion("OptiLCMS"),
                       "\nPlease define your running plan ..."), sleep = time_interval1);
  ##################### 0
  MessageOutput("Commands Origanization Finished!", progress = 1, sleep = time_interval1);
  MessageOutput("Step 0/6: Scanning ROIs for parameters optimization...", progress = 2, sleep = time_interval1);
  MessageOutput("QC_PREFA02.mzML", progress = 2, sleep = time_interval1);
  MessageOutput("QC_PREFB02.mzML", progress = 2, sleep = time_interval1);
  MessageOutput("Data Loading...", progress = 2, sleep = time_interval1);
  MessageOutput("Raw file import begin...", progress = 2, sleep = time_interval1);
  MessageOutput("Reading MS from QC_PREFA02.mzML ...", progress = 2, sleep = time_interval1);
  MessageOutput("Importing QC_PREFA02.mzML:", progress = 2, ecol = "", sleep = time_interval1);
  for(i in seq(0,100,20)){
    MessageOutput(paste0(i,"% | "), progress = 2, ecol = "", sleep = 0.35);
  }
  MessageOutput(" Done!", progress = 2, sleep = time_interval1);
  MessageOutput("Reading from QC_PREFA02.mzML finished successfully !", progress = 2, sleep = time_interval1);
  
  MessageOutput("Reading MS from QC_PREFB02.mzML ...", progress = 2, sleep = time_interval1);
  MessageOutput("Importing QC_PREFB02.mzML:", progress = 2, ecol = "", sleep = time_interval1);
  for(i in seq(0,100,20)){
    MessageOutput(paste0(i,"% | "), progress = 2, ecol = "", sleep = 0.35);
  }
  MessageOutput(" Done!", progress = 2, sleep = time_interval1);
  MessageOutput("Reading from QC_PREFB02.mzML finished successfully !", progress = 2, sleep = time_interval1);
  MessageOutput("Data Loaded !", progress = 2, sleep = time_interval1);
  MessageOutput("Empty Scan detecting...", progress = 2, sleep = time_interval1);
  MessageOutput("No Empty scan found !", progress = 2, sleep = time_interval1);
  MessageOutput("Identifying regions of interest (ROI)...", progress = 2, sleep = time_interval1);
  MessageOutput("Identifying regions of potential contaminants", ecol = "", progress = 2, sleep = time_interval1);
  for(i in 1:30){
    MessageOutput(".", ecol = "", progress = 3, sleep = 0.5);
  }
  MessageOutput(" Done!", progress = 3, sleep = time_interval1);
  
  MessageOutput("300 potential contaminamts will not be used for parameters optimization !", progress = 3, sleep = time_interval1);
  MessageOutput("Going to the next step...", progress = 3, sleep = time_interval1);
  MessageOutput("MS data Preparing...", progress = 3, sleep = time_interval1);
  MessageOutput("MS Data ready !", progress = 3, sleep = time_interval1);
  MessageOutput("Identifying ROIs in m/z dimensions...", progress = 4, sleep = time_interval1);
  MessageOutput("Identifying ROIs in m/z dimensions Done !", progress = 4, sleep = time_interval1);
  MessageOutput("Identifying ROIs in RT dimensions...", progress = 4, sleep = time_interval1);
  MessageOutput("Identification on ROIs Finished!", progress = 4, sleep = time_interval1);
  MessageOutput("Optimization will be started soon...", progress = 4, sleep = time_interval1);
  
  ##################### 1
  MessageOutput("\nStep 1/6: Start to optimize parameters! ", progress = 4, sleep = time_interval1);
  MessageOutput("This step may take a long time...", progress = 4, sleep = time_interval1);
  MessageOutput("DoE Optimization Starting Now...", progress = 4, sleep = time_interval1);
  MessageOutput("Evaluating Noise level...", progress = 4, sleep = time_interval1);
  MessageOutput("Done!", progress = 5, sleep = time_interval1);
  MessageOutput("Preparing Parameters for optimization finished !", progress = 5, sleep = time_interval1);
  MessageOutput("Data Spliting Finished !", progress = 5, sleep = time_interval1);
  MessageOutput("Peak Preparing Begin...", progress = 5, sleep = time_interval1);
  MessageOutput("Peak Preparing Done !", progress = 5, sleep = time_interval1);
  MessageOutput("Round:1", progress = 5, sleep = time_interval1);
  MessageOutput("DoE Running Begin...", progress = 5, sleep = time_interval1);
  MessageOutput("Finished: ", ecol = "", progress = 5, sleep = time_interval1);
  for(i in c(11,22,33,44,56,67,78,89,100)){
    MessageOutput(paste0(i,"% | "), progress = 5+i/50, ecol = "", sleep = 2);
  }
  MessageOutput("Done !", progress = 7, sleep = time_interval1);
  MessageOutput("Round 1 Finished !", progress = 7, sleep = time_interval1);
  MessageOutput("Model Parsing...", progress = 7, sleep = time_interval1);
  MessageOutput("Gaussian peak ratio (%): 59.7.", progress = 7, sleep = time_interval1);
  MessageOutput("Model Parsing Done !", progress = 8, sleep = time_interval1);
  
  
  MessageOutput("\nRound:2", progress = 8, sleep = time_interval1);
  MessageOutput("DoE Running Begin...", progress = 8, sleep = time_interval1);
  MessageOutput("Finished: ", ecol = "", progress = 8, sleep = time_interval1);
  for(i in c(11,22,33,44,56,67,78,89,100)){
    MessageOutput(paste0(i,"% | "), progress = 8+i/50, ecol = "", sleep = 2);
  }
  MessageOutput("Done !", progress = 10, sleep = time_interval1);
  MessageOutput("Round 2 Finished !", progress = 10, sleep = time_interval1);
  MessageOutput("Model Parsing...", progress = 10, sleep = time_interval1);
  MessageOutput("Gaussian peak ratio (%): 61.5.", progress = 10, sleep = time_interval1);
  MessageOutput("Model Parsing Done !", progress = 11, sleep = time_interval1);
  
  MessageOutput("\nNo Increase Stopping !", progress = 15, sleep = time_interval1);
  MessageOutput(paste0("Step 1/6: Parameters Optimization Finished ! (", Sys.time(),")"), progress = 16, sleep = time_interval1);
  
  ##################### 1
  MessageOutput("\nStep 2/6: Start to import the spectrum! ", progress = 16, sleep = time_interval1);
  MessageOutput("This step will take a short time...", progress = 16, sleep = time_interval1);
  MessageOutput("Raw file import begin...", progress = 16, sleep = time_interval1);
  MessageOutput("CD-6KUCT.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("CD-77FXR.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("CD-9OS5Y.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("CD-9WOBP.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("HC-9SNJ4.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("HC-9X47O.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("HC-AMR37.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("HC-AUP8B.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("QC_PREFA02.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("QC_PREFB02.mzML import done!", progress = 16, sleep = time_interval1);
  MessageOutput("Raw file initialized Successfully!", progress = 17, sleep = time_interval1);
  MessageOutput("Plotting BPIS and TICS.", progress = 17, sleep = time_interval1);
  MessageOutput(paste0("Step 1/6: Successfully imported raw MS data! (",
                       Sys.time(),
                       ") \nGoing to the next step..."), progress = 18, sleep = time_interval1);
  
  ##################### 3
  MessageOutput("\nStep 3/6: Started peak picking! This step will take some time...", progress = 22, sleep = time_interval1);
  MessageOutput("Detecting peaks in 3942 regions of interest ... OK: 877 found.", progress = 25, sleep = time_interval1);
  MessageOutput("Detecting peaks in 3968 regions of interest ... OK: 860 found.", progress = 28, sleep = time_interval1);
  MessageOutput("Detecting peaks in 4200 regions of interest ... OK: 835 found.", progress = 31, sleep = time_interval1);
  MessageOutput("Detecting peaks in 4239 regions of interest ... OK: 1121 found.", progress = 34, sleep = time_interval1);
  MessageOutput("Detecting peaks in 4235 regions of interest ... OK: 962 found.", progress = 37, sleep = time_interval1);
  MessageOutput("Detecting peaks in 4361 regions of interest ... OK: 1021 found.", progress = 40, sleep = time_interval1);
  MessageOutput("Detecting peaks in 4462 regions of interest ... OK: 1118 found.", progress = 43, sleep = time_interval1);
  MessageOutput("Detecting peaks in 4361 regions of interest ... OK: 1177 found.", progress = 46, sleep = time_interval1);
  MessageOutput("Detecting peaks in 4330 regions of interest ... OK: 1086 found.", progress = 47, sleep = time_interval1);
  MessageOutput("Detecting peaks in 4489 regions of interest ... OK: 1024 found.", progress = 50, sleep = time_interval1);
  MessageOutput(paste0("Step 3/6: Peak picking finished ! (",Sys.time(),")"), progress = 50, sleep = time_interval1);
  ##################### 4
  MessageOutput("\nStep 4/6: Started peak alignment! This step is running...", progress = 51, sleep = time_interval1);
  MessageOutput("Total of 2795 slices detected for processing... Done !", progress = 52, sleep = time_interval1);
  MessageOutput("Going to the next step...", progress = 55, sleep = time_interval2);
  MessageOutput("Retention time correction is running.PeakGroup -- loess is used for retention time correction.....", progress = 60, sleep = time_interval2);
  MessageOutput("Performing retention time correction using  174  peak groups.", progress = 61, sleep = time_interval2);
  MessageOutput("Applying retention time adjustment to the identified chromatographic peaks ... Done !", progress = 61, sleep = time_interval2);
  MessageOutput("Total of 2795 slices detected for processing... Done !", progress = 65, sleep = time_interval2);
  MessageOutput(paste0("Step 4/6: Peak alignment finished ! (",Sys.time(),")"), progress = 70, sleep = time_interval2);
  MessageOutput("Going to the next step...", progress = 73, sleep = time_interval2);
  ##################### 5
  MessageOutput("\nStep 5/6: Started peak filling! This step may take some time...", progress = 74, sleep = time_interval1);
  MessageOutput("Starting peak filling!", progress = 75, sleep = time_interval1);
  MessageOutput("Defining peak areas for filling-in....OK", progress = 76, sleep = time_interval1);
  MessageOutput("Start integrating peak areas from original files...", progress = 77, sleep = time_interval1);
  MessageOutput("Requesting 408 peaks from CD-6KUCT.mzML ... got 240.", progress = 78, sleep = time_interval1);
  MessageOutput("Requesting 580 peaks from CD-9WOBP.mzML ... got 267.", progress = 78.5, sleep = time_interval1);
  MessageOutput("Requesting 535 peaks from HC-9X47O.mzML ... got 299.", progress = 79, sleep = time_interval1);
  MessageOutput("Requesting 507 peaks from HC-AUP8B.mzML ... got 314.", progress = 79, sleep = time_interval1);
  MessageOutput("Requesting 367 peaks from CD-77FXR.mzML ... got 210.", progress = 80, sleep = time_interval1);
  MessageOutput("Requesting 53 peaks from QC_PREFA02.mzML ... got 46.", progress = 81, sleep = time_interval1);
  MessageOutput("Requesting 65 peaks from HC-AMR37.mzML ... got 58.", progress = 82, sleep = time_interval1);
  MessageOutput("Requesting 486 peaks from QC_PREFB02.mzML ... got 296.", progress = 83, sleep = time_interval1);
  MessageOutput("Requesting 479 peaks from HC-9SNJ4.mzML ... got 302.", progress = 84, sleep = time_interval1);
  MessageOutput("Requesting 402 peaks from CD-9OS5Y.mzML ... got 264.", progress = 85, sleep = time_interval1);
  MessageOutput(paste0("Step 5/6: Peak filing finished ! (",Sys.time(),")"), progress = 86, sleep = time_interval1);
  MessageOutput("Peak Profiling finished successfully !", progress = 87, sleep = time_interval1);
  MessageOutput("Begin to plotting figures...Done !", progress = 88, sleep = time_interval1);
  ##################### 6
  MessageOutput("\nStep 6/6: Starting Peak Annotation...", progress = 88, sleep = time_interval1);
  MessageOutput("Start grouping after retention time.", progress = 88, sleep = time_interval1);
  MessageOutput("Created 45 pseudospectra.", progress = 89, sleep = time_interval1);
  MessageOutput("Generating peak matrix...", progress = 89, sleep = time_interval1);
  MessageOutput("Run isotope peak annotation..", progress = 90, sleep = time_interval1);
  MessageOutput("Found isotopes:113", progress =90, sleep = time_interval1);
  MessageOutput("Generating EIC's....", progress = 91, sleep = time_interval1);
  MessageOutput("Detecting  CD-6KUCT.mzML  ... 67 peaks found!", progress = 91, sleep = time_interval1);
  MessageOutput("Detecting  CD-77FXR.mzML  ... 60 peaks found!", progress = 92, sleep = time_interval1);
  MessageOutput("Detecting  CD-9OS5Y.mzML  ... 152 peaks found!", progress = 92, sleep = time_interval1);
  MessageOutput("Detecting  CD-9WOBP.mzML  ... 166 peaks found!", progress = 93, sleep = time_interval1);
  MessageOutput("Detecting  HC-9SNJ4.mzML  ... 28 peaks found!", progress = 93, sleep = time_interval1);
  MessageOutput("Detecting  HC-9X47O.mzML  ... 1 peaks found!", progress = 93, sleep = time_interval1);
  MessageOutput("Detecting  HC-AMR37.mzML  ... 283 peaks found!", progress = 94, sleep = time_interval1);
  MessageOutput("Detecting  HC-AUP8B.mzML  ... 16 peaks found!", progress = 94, sleep = time_interval1);
  MessageOutput("Detecting  QC_PREFA02.mzML  ... 46 peaks found!", progress = 95, sleep = time_interval1);
  MessageOutput("Detecting  QC_PREFB02.mzML  ... 137 peaks found!", progress = 95, sleep = time_interval1);
  MessageOutput("Warning: Found NA peaks in selected sample.", progress = 95, sleep = time_interval1);
  MessageOutput("Calculating peak correlations in 54 Groups...", progress = 96, sleep = time_interval1);
  MessageOutput("Calculating graph cross linking in 54 Groups...", progress = 96, sleep = time_interval1);
  MessageOutput("New number of ps-groups:  614", progress = 96, sleep = time_interval1);
  MessageOutput("mSet has now 361 groups, instead of 54 !", progress = 97, sleep = time_interval1);
  MessageOutput("Generating peak matrix for peak annotation!", progress = 97, sleep = time_interval1);
  MessageOutput("Polarity is set in annotaParam: negative", progress = 98, sleep = time_interval1);
  MessageOutput("Ruleset could not read from object! Recalculating...", progress = 98, sleep = time_interval1);
  MessageOutput("Calculating possible adducts in 614 Groups... ", progress = 99, sleep = time_interval1);
  MessageOutput(paste0("Step 6/6: Successfully performed peak annotation! (",Sys.time(),")"), progress = 99, sleep = time_interval1);
  MessageOutput("Going to the final step...Done!", progress = 99, sleep = time_interval2);
  
  MessageOutput(paste0("\nEverything has been finished Successfully ! (",Sys.time(),")"), progress = 100, sleep = 0);
  
}

#' @noRd
FastRunningShow_customized <- function(fullUserPath){
  
  oldwd <- getwd();
  on.exit(setwd(oldwd));
  
  setwd(fullUserPath);    
  
  time_interval1 <- 0.80;
  time_interval2 <- 2;
  # running to show the progress
  MessageOutput(paste0("Running Status -- Plan Initialized Successfully at: ", 
                         Sys.time(), 
                         "\nCurrent OptiLCMS version is ",
                         packageVersion("OptiLCMS"),
                         "\nPlease define your running plan ..."), sleep = time_interval1);
  ##################### 1
  MessageOutput("Commands Origanization Finished!", progress = 1, sleep = time_interval1);
  MessageOutput("\nStep 1/6: Start to import the spectrum! ", progress = 2, sleep = time_interval1);
  MessageOutput("This step will take a short time...", progress = 3, sleep = time_interval1);
  MessageOutput("Raw file import begin...", progress = 4, sleep = time_interval1);
  MessageOutput("CD-6KUCT.mzML import done!", progress = 5, sleep = time_interval1);
  MessageOutput("CD-77FXR.mzML import done!", progress = 6, sleep = time_interval1);
  MessageOutput("CD-9OS5Y.mzML import done!", progress = 7, sleep = time_interval1);
  MessageOutput("CD-9WOBP.mzML import done!", progress = 8, sleep = time_interval1);
  MessageOutput("HC-9SNJ4.mzML import done!", progress = 9, sleep = time_interval1);
  MessageOutput("HC-9X47O.mzML import done!", progress = 10, sleep = time_interval1);
  MessageOutput("HC-AMR37.mzML import done!", progress = 11, sleep = time_interval1);
  MessageOutput("HC-AUP8B.mzML import done!", progress = 12, sleep = time_interval1);
  MessageOutput("QC_PREFA02.mzML import done!", progress = 13, sleep = time_interval1);
  MessageOutput("QC_PREFB02.mzML import done!", progress = 14, sleep = time_interval1);
  MessageOutput("Raw file initialized Successfully!", progress = 15, sleep = time_interval1);
  MessageOutput("Plotting BPIS and TICS.", progress = 16, sleep = time_interval1);
  MessageOutput(paste0("Step 1/6: Successfully imported raw MS data! (",
                       Sys.time(),
                       ") \nGoing to the next step..."), progress = 18, sleep = time_interval1);
  ##################### 2
  MessageOutput("\nStep 2/6: Internalize parameters!", progress = 18, sleep = time_interval1);
  MessageOutput("This step will be finished soon...", progress = 18, sleep = time_interval1);
  MessageOutput("Step 2/6: Parameters Internalized Successfully!", progress = 19, sleep = time_interval1);
  MessageOutput("Going to the next step...", progress = 19, sleep = time_interval1);
  MessageOutput("Parameters for centWave have been successfully parsed!", progress = 20, sleep = time_interval1);
  MessageOutput("4 CPU Threads will be used for peak profiling !", progress = 21, sleep = time_interval1);
  ##################### 3
  MessageOutput("\nStep 3/6: Started peak picking! This step will take some time...", progress = 22, sleep = time_interval1);
  MessageOutput("Detecting peaks in 5509 regions of interest ... OK: 1285 found.", progress = 25, sleep = time_interval1);
  MessageOutput("Detecting peaks in 6009 regions of interest ... OK: 1327 found.", progress = 28, sleep = time_interval1);
  MessageOutput("Detecting peaks in 5355 regions of interest ... OK: 1395 found.", progress = 31, sleep = time_interval1);
  MessageOutput("Detecting peaks in 5767 regions of interest ... OK: 1810 found.", progress = 34, sleep = time_interval1);
  MessageOutput("Detecting peaks in 5832 regions of interest ... OK: 1507 found.", progress = 37, sleep = time_interval1);
  MessageOutput("Detecting peaks in 6273 regions of interest ... OK: 1586 found.", progress = 40, sleep = time_interval1);
  MessageOutput("Detecting peaks in 6061 regions of interest ... OK: 1643 found.", progress = 43, sleep = time_interval1);
  MessageOutput("Detecting peaks in 5926 regions of interest ... OK: 1805 found.", progress = 46, sleep = time_interval1);
  MessageOutput("Detecting peaks in 6054 regions of interest ... OK: 1646 found.", progress = 47, sleep = time_interval1);
  MessageOutput("Detecting peaks in 6299 regions of interest ... OK: 1608 found.", progress = 50, sleep = time_interval1);
  MessageOutput(paste0("Step 3/6: Peak picking finished ! (",Sys.time(),")"), progress = 50, sleep = time_interval1);
  ##################### 4
  MessageOutput("\nStep 4/6: Started peak alignment! This step is running...", progress = 51, sleep = time_interval1);
  MessageOutput("Total of 2805 slices detected for processing... Done !", progress = 52, sleep = time_interval1);
  MessageOutput("Going to the next step...", progress = 55, sleep = time_interval2);
  MessageOutput("Retention time correction is running.PeakGroup -- loess is used for retention time correction.....", progress = 60, sleep = time_interval2);
  MessageOutput("Performing retention time correction using  54  peak groups.", progress = 61, sleep = time_interval2);
  MessageOutput("Applying retention time adjustment to the identified chromatographic peaks ... Done !", progress = 61, sleep = time_interval2);
  MessageOutput("Total of 2805 slices detected for processing... Done !", progress = 65, sleep = time_interval2);
  MessageOutput(paste0("Step 4/6: Peak alignment finished ! (",Sys.time(),")"), progress = 70, sleep = time_interval2);
  MessageOutput("Going to the next step...", progress = 73, sleep = time_interval2);
  ##################### 5
  MessageOutput("\nStep 5/6: Started peak filling! This step may take some time...", progress = 74, sleep = time_interval1);
  MessageOutput("Starting peak filling!", progress = 75, sleep = time_interval1);
  MessageOutput("Defining peak areas for filling-in....OK", progress = 76, sleep = time_interval1);
  MessageOutput("Start integrating peak areas from original files...", progress = 77, sleep = time_interval1);
  MessageOutput("Requesting 139 peaks from CD-6KUCT.mzML ... got 56.", progress = 78, sleep = time_interval1);
  MessageOutput("Requesting 256 peaks from CD-9WOBP.mzML ... got 55.", progress = 78.5, sleep = time_interval1);
  MessageOutput("Requesting 261 peaks from HC-9X47O.mzML ... got 101.", progress = 79, sleep = time_interval1);
  MessageOutput("Requesting 234 peaks from HC-AUP8B.mzML ... got 95.", progress = 79, sleep = time_interval1);
  MessageOutput("Requesting 146 peaks from CD-77FXR.mzML ... got 61.", progress = 80, sleep = time_interval1);
  MessageOutput("Requesting 52 peaks from QC_PREFA02.mzML ... got 30.", progress = 81, sleep = time_interval1);
  MessageOutput("Requesting 199 peaks from HC-AMR37.mzML ... got 72.", progress = 82, sleep = time_interval1);
  MessageOutput("Requesting 206 peaks from HC-9SNJ4.mzML ... got 94.", progress = 83, sleep = time_interval1);
  MessageOutput("Requesting 54 peaks from QC_PREFB02.mzML ... got 37.", progress = 84, sleep = time_interval1);
  MessageOutput("Requesting 188 peaks from CD-9OS5Y.mzML ... got 74.", progress = 85, sleep = time_interval1);
  MessageOutput(paste0("Step 5/6: Peak filing finished ! (",Sys.time(),")"), progress = 86, sleep = time_interval1);
  MessageOutput("Peak Profiling finished successfully !", progress = 87, sleep = time_interval1);
  MessageOutput("Begin to plotting figures...Done !", progress = 88, sleep = time_interval1);
  ##################### 6
  MessageOutput("\nStep 6/6: Starting Peak Annotation...", progress = 88, sleep = time_interval1);
  MessageOutput("Start grouping after retention time.", progress = 88, sleep = time_interval1);
  MessageOutput("Created 62 pseudospectra.", progress = 89, sleep = time_interval1);
  MessageOutput("Generating peak matrix...", progress = 89, sleep = time_interval1);
  MessageOutput("Run isotope peak annotation..", progress = 90, sleep = time_interval1);
  MessageOutput("Found isotopes:28", progress =90, sleep = time_interval1);
  MessageOutput("Generating EIC's....", progress = 91, sleep = time_interval1);
  MessageOutput("Detecting  CD-6KUCT.mzML  ... 66 peaks found!", progress = 91, sleep = time_interval1);
  MessageOutput("Detecting  CD-77FXR.mzML  ... 35 peaks found!", progress = 92, sleep = time_interval1);
  MessageOutput("Detecting  CD-9OS5Y.mzML  ... 219 peaks found!", progress = 92, sleep = time_interval1);
  MessageOutput("Detecting  CD-9WOBP.mzML  ... 120 peaks found!", progress = 93, sleep = time_interval1);
  MessageOutput("Detecting  HC-9SNJ4.mzML  ... 30 peaks found!", progress = 93, sleep = time_interval1);
  MessageOutput("Detecting  HC-9X47O.mzML  ... 0 peaks found!", progress = 93, sleep = time_interval1);
  MessageOutput("Detecting  HC-AMR37.mzML  ... 287 peaks found!", progress = 94, sleep = time_interval1);
  MessageOutput("Detecting  HC-AUP8B.mzML  ... 21 peaks found!", progress = 94, sleep = time_interval1);
  MessageOutput("Detecting  QC_PREFA02.mzML  ... 89 peaks found!", progress = 95, sleep = time_interval1);
  MessageOutput("Detecting  QC_PREFB02.mzML  ... 39 peaks found!", progress = 95, sleep = time_interval1);
  MessageOutput("Warning: Found NA peaks in selected sample.", progress = 95, sleep = time_interval1);
  MessageOutput("Calculating peak correlations in 62 Groups...", progress = 96, sleep = time_interval1);
  MessageOutput("Calculating graph cross linking in 62 Groups...", progress = 96, sleep = time_interval1);
  MessageOutput("New number of ps-groups:  736", progress = 96, sleep = time_interval1);
  MessageOutput("mSet has now 736 groups, instead of 62 !", progress = 97, sleep = time_interval1);
  MessageOutput("Generating peak matrix for peak annotation!", progress = 97, sleep = time_interval1);
  MessageOutput("Polarity is set in annotaParam: negative", progress = 98, sleep = time_interval1);
  MessageOutput("Ruleset could not read from object! Recalculating...", progress = 98, sleep = time_interval1);
  MessageOutput("Calculating possible adducts in 736 Groups... ", progress = 99, sleep = time_interval1);
  MessageOutput(paste0("Step 6/6: Successfully performed peak annotation! (",Sys.time(),")"), progress = 99, sleep = time_interval1);
  MessageOutput("Going to the final step...Done!", progress = 99, sleep = time_interval2);
  
  MessageOutput(paste0("\nEverything has been finished Successfully ! (",Sys.time(),")"), progress = 100, sleep = 0);
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

