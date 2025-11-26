#' @title Perform ROI Extraction from raw MS data (PerformDataTrimming)
#' @description This function performs the raw data trimming. This function will output 
#' an trimmed MSnExp file to memory or hardisk according to the choice of users must 
#' provide the data path for 'datapath', and optionally provide other corresponding parameters.
#' @param datapath Character, the path of the raw MS data files' or folder's path (.mzXML, .CDF and .mzML) 
#' for parameters training.
#' @param mode Character, mode for data trimming to select the chraracteristic peaks. 
#' Default is 'ssm'. Users could select random trimed according to mz value (mz_random) or 
#' RT value (rt_random). Besides, specific peaks at certain mz (mz_specific) or 
#' RT (rt_specific) could also be extracted. 'none' will not trim the data.
#' @param mz Numeric, mz value(s) for specific selection. Positive values means including (the values 
#' indicted) and negative value means excluding/removing.
#' @param mzdiff Numeric, the deviation (ppm) of mz value(s).
#' @param rt Numeric, rt value for specific selection. Positive values means including 
#' and negative value means excluding.
#' @param rtdiff Numeric, the deviation (seconds) of rt value(s).
#' @param rt.idx Numeric, the relative rt (retention time) range, from 0 to 1. 1 means all retention time
#' will be retained, while 0 means none. Default is 1/15. If default rt.idx produce too few peaks, 
#' please consider increasing this value.
#' @param write Logical, if true, will write the trimmed data to the directory 'trimmed' folder 
#' in the datapath. The data in memory will be kept.
#' @param plot Logical, if TRUE, will plot the chromatogram of the trimmed data.
#' @param rmConts LOgical, whether to exclude/remove the potential contamination for parameters optimization. Default is TRUE.
#' @param running.controller The resuming pipeline running controller. Optional. Don't need to define by hand.
#' @export
#' @import MSnbase
#' @import progress
#' @import Biobase
#' @import RColorBrewer
#' @import tools
#' @return will return an mSet objects with extracted ROI
#' @seealso \code{\link{PerformROIExtraction}} for the new version of this function.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
#' @examples
#' DataFiles <- dir(system.file("mzData", package = "mtbls2"), full.names = TRUE, recursive = TRUE)
#' # mSet <- PerformDataTrimming(datapath = DataFiles[1],rt.idx = 0.025, rmConts = FALSE);


PerformDataTrimming<- function(datapath, mode="ssm", write=FALSE, mz, mzdiff, rt, rtdiff, 
                                 rt.idx=1/15, rmConts = TRUE, plot=TRUE, running.controller=NULL){
  
  PerformROIExtraction(datapath, mode=mode, write, mz, mzdiff, rt, rtdiff, 
                       rt.idx, rmConts = rmConts, plot,running.controller);
  
}

#' @title Perform ROI Extraction from raw MS data
#' @description This function performs the raw data trimming. This function will output 
#' an trimmed MSnExp file to memory or hardisk according to the choice of users must 
#' provide the data path for 'datapath', and optionally provide other corresponding parameters.
#' @param datapath Character, the path of the raw MS data files' or folder's path (.mzXML, .CDF and .mzML) 
#' for parameters training.
#' @param mode Character, mode for data trimming to select the chraracteristic peaks. 
#' Default is 'ssm'. Users could select random trimed according to mz value (mz_random) or 
#' RT value (rt_random). Besides, specific peaks at certain mz (mz_specific) or 
#' RT (rt_specific) could also be extracted. 'none' will not trim the data.
#' @param mz Numeric, mz value(s) for specific selection. Positive values means including (the values 
#' indicted) and negative value means excluding/removing.
#' @param mzdiff Numeric, the deviation (ppm) of mz value(s).
#' @param rt Numeric, rt value for specific selection. Positive values means including 
#' and negative value means excluding.
#' @param rtdiff Numeric, the deviation (seconds) of rt value(s).
#' @param rt.idx Numeric, the relative rt (retention time) range, from 0 to 1. 1 means all retention time
#' will be retained, while 0 means none. Default is 1/15. If default rt.idx produce too few peaks, 
#' please consider increasing this value.
#' @param write Logical, if true, will write the trimmed data to the directory 'trimmed' folder 
#' in the datapath. The data in memory will be kept.
#' @param plot Logical, if TRUE, will plot the chromatogram of the trimmed data.
#' @param rmConts LOgical, whether to exclude/remove the potential contamination for parameters optimization. Default is TRUE.
#' @param running.controller The resuming pipeline running controller. Optional. Don't need to define by hand.
#' @export
#' @import MSnbase
#' @import progress
#' @import Biobase
#' @import RColorBrewer
#' @import tools
#' @return will return an mSet objects with extracted ROI
#' @seealso \code{\link{PerformDataTrimming}} for the old version of this function.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
#' @examples
#' DataFiles <- dir(system.file("mzData", package = "mtbls2"), full.names = TRUE, recursive = TRUE)
#' # mSet <- PerformROIExtraction(datapath = DataFiles[1],rt.idx = 0.025, rmConts = FALSE);


PerformROIExtraction <-
  function(datapath,
           mode = "ssm",
           write = FALSE,
           mz,
           mzdiff,
           rt,
           rtdiff,
           rt.idx = 1/15,
           rmConts = TRUE,
           plot = TRUE,
           running.controller = NULL) {
    
    if(!exists(".SwapEnv")){
      .SwapEnv <<- new.env(parent = .GlobalEnv);
      .SwapEnv$.optimize_switch <- TRUE;
      .SwapEnv$count_current_sample <- 0;
      .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
      .SwapEnv$envir <- new.env();
    }
    
    .SwapEnv$.optimize_switch <- TRUE;
    
    datapath <- normalizePath(datapath);
    
    if(!.on.public.web()){
      #Local R package running
      files <- Path2Files(path = datapath);
      files <- unlist(sapply(files, tools::file_path_as_absolute));
      
    } else {
      # code for web version
      if (!dir.exists(datapath)) {
        datapath <- "upload/";
      }
      if (!dir.exists(datapath)) {
        datapath <- "/home/glassfish/projects/MetaboDemoRawData/upload/QC/";
      }
      files <- Path2Files(path = datapath);
      files <- unlist(sapply(files, tools::file_path_as_absolute));
    }

    MessageOutput(
      "Step 0/6: Scanning ROIs for parameters optimization...",
      ecol = "\n",
      progress = 1.0
    )
    
    match.arg(mode,
              choices = c(
                "ssm",
                "mz_random",
                "rt_random",
                "mz_specific",
                "rt_specific"
              ))
    
    start.time <- Sys.time();
    function.name <- "ROI_extract";
    parallel_param <- BiocParallel::bpparam();
    
    combine_msnexp_list <- function(msnexp_list) {
      if (length(msnexp_list) == 1) {
        return(msnexp_list[[1]])
      }
      
      assayData_env <- new.env(parent = emptyenv())
      feature_data_list <- list()
      pheno_data_list <- list()
      assay_name_vec <- character()
      
      for (idx in seq_along(msnexp_list)) {
        obj <- msnexp_list[[idx]]
        ad_names <- ls(obj@assayData)
        new_ad_names <- sprintf("F%d.S%03d", idx, seq_along(ad_names))
        assay_name_vec <- c(assay_name_vec, new_ad_names)
        
        pheno_row <- obj@phenoData@data
        if (is.null(rownames(pheno_row))) {
          rownames(pheno_row) <- pheno_row$sample_name
        }
        pheno_data_list[[idx]] <- pheno_row
        
        fd <- obj@featureData@data
        if (nrow(fd) == length(new_ad_names)) {
          rownames(fd) <- new_ad_names
        }
        feature_data_list[[idx]] <- fd
        
        for (nm_idx in seq_along(ad_names)) {
          spectrum_obj <- obj@assayData[[ad_names[nm_idx]]]
          if (!is.null(slotNames(spectrum_obj)) && "fromFile" %in% slotNames(spectrum_obj)) {
            spectrum_obj@fromFile <- as.integer(idx)
          }
          assign(new_ad_names[nm_idx], spectrum_obj, envir = assayData_env)
        }
      }
      
      pd_data <- do.call(rbind, pheno_data_list)
      pd_data <- data.frame(pd_data,
                            row.names = as.character(seq_len(nrow(pd_data))),
                            stringsAsFactors = FALSE)
      combined_pd <- new("AnnotatedDataFrame", data = pd_data)
      attr(combined_pd@data, "row.names") <- seq_len(nrow(pd_data))
      
      combined_fd <- new("AnnotatedDataFrame",
                         data = do.call(rbind, feature_data_list),
                         varMetadata = msnexp_list[[1]]@featureData@varMetadata)
      combined_fd <- combined_fd[assay_name_vec, , drop = FALSE]
      
      processing_data <- msnexp_list[[1]]@processingData
      processing_data@files <- unname(unlist(lapply(msnexp_list, function(x) x@processingData@files)))
      
      combined_experiment <- msnexp_list[[1]]@experimentData
      exp_slots <- slotNames(combined_experiment)
      file_count <- length(processing_data@files)
      replicate_slots <- c("analyser", "detectorType",
                           "instrumentManufacturer", "instrumentModel", "ionSource")
      for (sn in exp_slots) {
        val <- slot(combined_experiment, sn)
        if (sn %in% replicate_slots && is.vector(val) && length(val) == 1 && file_count > 1) {
          slot(combined_experiment, sn) <- rep(val[1], file_count)
        }
      }
      
      protocol_df <- data.frame(row.names = character(0),
                                stringsAsFactors = FALSE)
      protocol_data <- new("AnnotatedDataFrame",
                           data = protocol_df,
                           dimLabels = c("sampleNames", "sampleColumns"))
      attr(protocol_data@data, "row.names") <- integer(0)
      
      new("MSnExp",
          assayData = assayData_env,
          phenoData = combined_pd,
          featureData = combined_fd,
          experimentData = combined_experiment,
          protocolData = protocol_data,
          processingData = processing_data)
    }
    
    if (is.null(running.controller)) {
      c1 <- c2 <- c3 <- c4 <- c5 <- TRUE;
      .running.as.plan <- FALSE;
    } else {
      c1 <- running.controller@ROI_extract[["c1"]]; # Data read part
      c2 <- running.controller@ROI_extract[["c2"]]; # Data trim part
      c3 <- running.controller@ROI_extract[["c3"]]; # Data write part
      c4 <- running.controller@ROI_extract[["c4"]]; # Data plot part
      c5 <- running.controller@ROI_extract[["c5"]]; # Data rmConts part
      .running.as.plan <- TRUE;
    }
    
    if (c1) {
      #Select first at least 2 QC data for optimization
      
      if (.on.public.web()) {
        #TODO: need to configure with the online pipeline
        mSet <- NULL;
        load("mSet.rda");
        
        dda_file <- files;
        rawfilenms <- basename(mSet@rawfiles);
        
        if (basename(datapath) == "QC") {
          QC_uploaded_list <- basename(dda_file);
          QC_index <- QC_uploaded_list %in% rawfilenms;
          QC_count <- length(which(QC_index));
          
          if (QC_count > 1) {
            # deal with the case: 2 or more QC provided/included
            if (QC_count > 2) {
              dda_file1 <- dda_file[QC_index][seq_len(3)]
            } else {
              dda_file1 <- dda_file[QC_index]
            }
            
          } else if (QC_count == 1) {
            # deal with the case: 1 QC provided/included
            # List all samples uploaded
            all_sample_list <-
              list.files(paste0(datapath, "/../"),
                         recursive = TRUE,
                         full.names = TRUE)
            
            # choose the included files
            included_samples <-
              all_sample_list[basename(all_sample_list) %in% rawfilenms]
            
            file_sizes <-
              sapply(
                included_samples,
                FUN = function(x) {
                  file.size(x)
                }
              )
            
            largest_file <-
              basename(names(sort(file_sizes, decreasing = TRUE)[seq_len(2)]))
            
            sample_index <-
              as.numeric(sapply(
                unique(c(
                  basename(dda_file[QC_index]), largest_file
                )),
                FUN = function(x) {
                  grep(x, included_samples)
                }
              ))
            
            dda_file1 <- included_samples[sample_index]
          } else {
            # deal with the case: 0 QC provided/included
            # List all samples uploaded
            all_sample_list <-
              list.files(paste0(datapath, "/../"),
                         recursive = TRUE,
                         full.names = TRUE)
            
            # choose the largest files
            included_samples <-
              all_sample_list[basename(all_sample_list) %in% rawfilenms]
            
            dda_file1 <-
              names(sort(
                sapply(
                  included_samples,
                  FUN = function(x) {
                    file.size(x)
                  }
                ),
                decreasing = TRUE
              )[seq_len(3)])
          }
          
        } else if (basename(datapath) == "upload") {
          included_files <- dda_file[basename(dda_file) %in% rawfilenms]
          QC_index <- grep("QC_*", included_files)
          
          if (!identical(QC_index, integer(0)) &
              length(QC_index) > 1) {
            dda_file1 <- included_files[QC_index][seq_len(3)]
            
          } else if (!identical(QC_index, integer(0)) &
                     length(QC_index) == 1) {
            file_sizes <-
              sapply(
                included_files,
                FUN = function(x) {
                  file.size(x)
                }
              )
            
            dda_file1 <-
              c(included_files[which(file_sizes == max(file_sizes))], included_files[QC_index])
            
          } else {
            
            file_sizes <-
              sapply(
                included_files,
                FUN = function(x) {
                  file.size(x)
                }
              )
            
            dda_file1 <-
              names(sort(file_sizes, decreasing = TRUE)[seq_len(3)])
          }
        }
        
      } else {
        mSet <- new("mSet");
        mSet@rawfiles <- dda_file1 <- files;
      }
      
      MessageOutput(paste0(basename(dda_file1), collapse = "\n"),"\n")
      
      pd <-
        data.frame(
          sample_name = tools::file_path_sans_ext(basename(dda_file1)),
          stringsAsFactors = FALSE
        )
      pdata_class <- if (isClass("NAnnotatedDataFrame")) "NAnnotatedDataFrame" else "AnnotatedDataFrame"
      
      MessageOutput("Data Loading in parallel...", "\n", NULL)
      
      raw_data_list <-
        BiocParallel::bplapply(seq_along(dda_file1),
                               FUN = function(idx) {
                                 tryCatch(
                                   read.MSdata(
                                     dda_file1[idx],
                                     pdata = new(pdata_class, pd[idx, , drop = FALSE]),
                                     msLevel. = 1L,
                                     mode = "inMemory"
                                   ),
                                   error = function(e) {
                                     e
                                   }
                                 )
                               },
                               BPPARAM = parallel_param)
      
      invalid_data <- which(!vapply(raw_data_list, function(x) {
        is(x, "MSnExp")
      }, logical(1)))
      
      if (length(invalid_data)) {
        first_error <- raw_data_list[[invalid_data[1]]]
        MessageOutput(
          mes = paste0(
            "<font color=\"red\">",
            "\nERROR:",
            first_error$message,
            "</font>"
          ),
          ecol = "\n",
          progress = 1
        )
        stop("EXCEPTION POINT CODE: RMS1")
      }
      
      MessageOutput("Data Loaded !", "\n", NULL)

      if (.running.as.plan) {
        cache.save(raw_data_list, paste0(function.name, "_c1"));
        marker_record(paste0(function.name, "_c1"));
      }
      
    } else {
      mSet <- new("mSet");
      mSet@rawfiles <- files;
      raw_data_list <- cache.read(function.name, "c1");
      pd <-
        data.frame(
          sample_name = tools::file_path_sans_ext(basename(mSet@rawfiles)),
          stringsAsFactors = FALSE
        )
      pdata_class <- if (isClass("NAnnotatedDataFrame")) "NAnnotatedDataFrame" else "AnnotatedDataFrame"
      marker_record(paste0(function.name, "_c1"));
    }
    
    if (c2) {
      if (!mode == "none") {
        if (missing(rt.idx)) {
          rt.idx <- 0.9
          # include the 90% RT range
        }
        
        if(.on.public.web()){
          rmConts <- FALSE;
          peakParams <- NULL;
          tmp_mes <- try(suppressWarnings(load("params.rda")),silent = TRUE);
        } else {
          tmp_mes <- 0;
          peakParams <- NULL;
        }
        
        MessageOutput("Empty Scan detecting...", "\n", NULL)

        prepare_single_file <- function(raw_data_single) {
          ms_list <-
            sapply(
              ls(raw_data_single@assayData),
              FUN = function(x)
                raw_data_single@assayData[[x]]@mz
            );
          
          emptyScan <- vector()
          
          for (i in seq_along(ms_list)) {
            if (identical(ms_list[[i]], numeric(0))) {
              emptyScan <- c(emptyScan, i)
            }
          }
          
          if (!identical(emptyScan, logical(0))) {
            raw_data_single <- .emptyscan.remove(raw_data_single, ms_list);
            ms_list <-
              sapply(
                ls(raw_data_single@assayData),
                FUN = function(x)
                  raw_data_single@assayData[[x]]@mz
              );
          }
          
          if(is(tmp_mes,"try-error") | rmConts){
            raw_data_single <- ContaminatsRemoval(raw_data_single, ms_list);
            raw_data_single <- .emptyscan.remove(raw_data_single, ms_list);
          } else if (.on.public.web()) {
            
            if(exists("peakParams")){
              if(peakParams[["rmConts"]]){
                raw_data_single <- ContaminatsRemoval(raw_data_single, ms_list);
                raw_data_single <- .emptyscan.remove(raw_data_single, ms_list);
              }
            } else {
              raw_data_single <- ContaminatsRemoval(raw_data_single, ms_list);
              raw_data_single <- .emptyscan.remove(raw_data_single, ms_list);
            }
          }
          
          return(raw_data_single)
        }
        
        if(c5){
          raw_data_list <- BiocParallel::bplapply(raw_data_list,
                                                  FUN = prepare_single_file,
                                                  BPPARAM = parallel_param)
          
          if (.running.as.plan) {
            cache.save(raw_data_list, paste0(function.name, "_c5"));
            marker_record(paste0(function.name, "_c5"));
          }
          
        } else {
          raw_data_list <- cache.read(function.name, "c5");
          marker_record(paste0(function.name, "_c5"));
        }

        compute_global_params <- function(raw_list, rt.idx) {
          spectra_mz_all <- unlist(lapply(raw_list, function(rd) unlist(lapply(MSnbase::spectra(rd), function(sp) MSnbase::mz(sp)))))
          spectra_int_all <- unlist(lapply(raw_list, function(rd) unlist(lapply(MSnbase::spectra(rd), function(sp) MSnbase::intensity(sp)))))
          
          bins.boundary <- seq(from = min(spectra_mz_all), to = max(spectra_mz_all), length.out = 5)
          
          if (length(spectra_mz_all) > 1000000) {
            rannum <- seq(from = 1, to = length(spectra_mz_all), by = 50)
          } else if (length(spectra_mz_all) > 100000) {
            rannum <- seq(from = 1, to = length(spectra_mz_all), by = 5)
          } else {
            rannum <- seq(length(spectra_mz_all))
          }
          
          spectra_mz_sample <- spectra_mz_all[rannum]
          spectra_abundance_sample <- spectra_int_all[rannum]
          
          spectra_mz_set <- lapply(seq_len(4), function(ii) {
            spectra_mz_sample[spectra_mz_sample > bins.boundary[ii] & spectra_mz_sample < bins.boundary[ii + 1]]
          })
          
          spectra_abundance_set <- lapply(seq_len(4), function(ii) {
            spectra_abundance_sample[spectra_mz_sample > bins.boundary[ii] & spectra_mz_sample < bins.boundary[ii + 1]]
          })
          
          good.bins.list <- lapply(
            seq_len(4),
            function(i) {
              binsb.low <- bins.boundary[i]
              binsb.high <- bins.boundary[i + 1]
              w <- 1
              bins.width <- (binsb.high - binsb.low) * 0.1
              wind.up <- wind.down <- 0
              inten.sum.set <- numeric()
              
              while (!wind.up > binsb.high) {
                wind.down <- binsb.low + (w - 1)
                wind.up <- binsb.low + bins.width + (w - 1)
                
                inten.sum.set <- c(
                  inten.sum.set,
                  sum(spectra_abundance_set[[i]][spectra_mz_set[[i]] > wind.down &
                                                  spectra_mz_set[[i]] < wind.up])
                )
                w <- w + 1
              }
              
              selected.bin.down <- binsb.low + which(max(inten.sum.set) == inten.sum.set) - 1
              selected.bin.up <- binsb.low + which(max(inten.sum.set) == inten.sum.set) - 1 + bins.width
              
              list(selected.bin.down, selected.bin.up)
            }
          )
          
          rt_vec <- unlist(lapply(raw_list, function(rd) {
            scan_names <- ls(rd@assayData)
            sapply(scan_names, function(nm) rd@assayData[[nm]]@rt)
          }))
          
          tic_vec <- unlist(lapply(raw_list, function(rd) {
            scan_names <- ls(rd@assayData)
            sapply(scan_names, function(nm) rd@assayData[[nm]]@tic)
          }))
          
          rt.bin.width <- (max(rt_vec) - min(rt_vec)) * rt.idx
          w <- 0
          rt.window.min <- min(rt_vec)
          tic.sum <- numeric()
          
          while (!rt.window.min > max(rt_vec)) {
            rt.window.min <- min(rt_vec) + w * 0.75
            rt.window.max <- min(rt_vec) + rt.bin.width + w * 0.75
            idx <- which(rt_vec > rt.window.min & rt_vec < rt.window.max)
            if (length(idx)) {
              tic.sum <- c(tic.sum, sum(tic_vec[idx]))
            }
            w <- w + 1
          }
          
          rt.boundary.lowlimit <- min(rt_vec) + which(max(tic.sum) == tic.sum)[ceiling(length(which(max(tic.sum) == tic.sum)) / 2)] * 0.75
          rt.boundary.uplimit <- min(rt_vec) + which(max(tic.sum) == tic.sum)[ceiling(length(which(max(tic.sum) == tic.sum)) / 2)] * 0.75 + rt.bin.width
          
          list(good.bins.list = good.bins.list,
               rt.boundary = c(rt.boundary.lowlimit, rt.boundary.uplimit))
        }

        global_params <- compute_global_params(raw_data_list, rt.idx)
        
        MessageOutput("Identifying regions of interest (ROI)...", "\n", NULL);
        
        trim_single_file <- function(raw_data_single) {
          ms_list <-
            sapply(
              ls(raw_data_single@assayData),
              FUN = function(x)
                raw_data_single@assayData[[x]]@mz
            );
          
          if (mode == "ssm") {
            trimed_MSnExp <- ssm_trim(raw_data_single,
                                      ms_list,
                                      rt.idx = rt.idx,
                                      good.bins.list = global_params$good.bins.list,
                                      rt.boundary = global_params$rt.boundary);
          }
          
          if (mode == "mz_random") {
            try(trimed_MSnExp <- mz.trim_random(raw_data_single, ms_list), silent = TRUE);
          }
          
          if (mode == "rt_random") {
            try(trimed_MSnExp <- rt.trim_random(raw_data_single, ms_list), silent = TRUE)
          }
          
          if (mode == "mz_specific") {
            trimed_MSnExp <- mz.trim_specific(raw_data_single, ms_list, mz, mzdiff = mzdiff)
          }
          
          if (mode == "rt_specific") {
            trimed_MSnExp <- rt.trim_specific(raw_data_single, ms_list, rt, rtdiff = rtdiff)
          }
          
          trimed_MSnExp <- .emptyscan.remove(trimed_MSnExp, ms_list);
          return(trimed_MSnExp)
        }
        
        trimed_list <- BiocParallel::bplapply(raw_data_list,
                                              FUN = trim_single_file,
                                              BPPARAM = parallel_param)
        
        invalid_trim <- which(!vapply(trimed_list, function(x) {
          is(x, "MSnExp")
        }, logical(1)))
        
        if (length(invalid_trim)) {
          stop("EXCEPTION POINT CODE: TRIM1")
        }
        
        trimed_MSnExp <- combine_msnexp_list(trimed_list);
        
      } else{
        trimed_MSnExp <- combine_msnexp_list(raw_data_list);
      }
      
      if (.running.as.plan) {
        cache.save(trimed_MSnExp, paste0(function.name, "_c2"));
        marker_record(paste0(function.name, "_c2"));
      }
      
    } else {
      trimed_MSnExp <- cache.read(function.name, "c2");
      marker_record(paste0(function.name, "_c2"));
    }
    
    MessageOutput(NULL, NULL, 4);

    if (c3) {
      if (write == TRUE) {
        MessageOutput("Data Writing...",ecol = "\n",NULL);
        writenames <-
          paste0(datapath,
                 "/trimmed/Trimmed_",
                 pd$sample_name,
                 ".mzML",
                 sep = "")
        dir.create(paste0(datapath, "/trimmed", collapse = ""))
        suppressMessages(writeMSData(trimed_MSnExp, writenames, outformat = "mzml"))
        MessageOutput("Data Writing Finished !",ecol = "\n",NULL);
      }
      
      if (.running.as.plan) {
        #cache.save(trimed_MSnExp,paste0(function.name,"_c2"));
        marker_record(paste0(function.name, "_c3"));
      }
    }
    
    if (c4) {
      if (plot == TRUE) {
        MessageOutput("Chromatogram Plotting Begin...",ecol = "\n",NULL);
        
        # if (.on.public.web()) {
        #   load_RColorBrewer();
        # } 
        
        ch.xdata <- chromatogram(trimed_MSnExp)
        group.col <-
          paste0(suppressWarnings(brewer.pal(length(
            trimed_MSnExp@processingData@files
          ), "Blues")));
        plot(ch.xdata, col = group.col[seq_along(trimed_MSnExp@processingData@files)])
      }
      
      if (.running.as.plan) {
        #cache.save(trimed_MSnExp,paste0(function.name,"_c2"));
        marker_record(paste0(function.name, "_c4"));
      }
      
    }

    MessageOutput(paste0("Identification on ROIs Finished!"),
                  ecol = "\n",
                  NULL);
    if(.running.as.plan){
      MessageOutput("Optimization will be started soon...", "\n", NULL);
    }

    pheno_tmp <- new("AnnotatedDataFrame",
                     data = pData(trimed_MSnExp),
                     varMetadata = varMetadata(phenoData(trimed_MSnExp)),
                     dimLabels = phenoData(trimed_MSnExp)@dimLabels)
    attr(pheno_tmp@data, "row.names") <- seq_len(nrow(pData(trimed_MSnExp)))
    trimed_MSnExp@phenoData <- pheno_tmp
    protocol_tmp <- new("AnnotatedDataFrame",
                        data = data.frame(row.names = character(0),
                                          stringsAsFactors = FALSE),
                        dimLabels = c("sampleNames", "sampleColumns"))
    trimed_MSnExp@protocolData <- protocol_tmp
    attr(trimed_MSnExp@protocolData@data, "row.names") <- integer(0)
    mSet@rawInMemory <- trimed_MSnExp;
    
    if(.on.public.web()){
      save(mSet, file = "mSet.rda");
    }
    
    return(mSet)
  }

#' @title Standards Simulation Method
#' @description Whole mass spectra will be divided as 4 bins according to the mz range. Trimming 
#' the raw with slide window method in every bins and retained the windows with highest scan intensity
#' and remove other scan signal in mz dimension. Then the data will be trimed again in the RT dimension
#' with slide window method. The window with highest intensity scans will be kept. After the timming
#' alongside mz and RT dimension, the peaks not only the high intensity peaks, but also the relatively 
#' low intensity peaks will also be retained as the 'simulated standards' data for parameters optimization.
#' @param raw_data MSnExp object, the raw data that has been read in memory.
#' @param ms_list List, the names list of all scans.
#' @param rt.idx Numeric, the retention time percentage, from 0 to 1. Default is 1/15.
#' @noRd
#' @import progress
#' @import BiocParallel
#' @import Biobase
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)

ssm_trim <- function(raw_data, ms_list, rt.idx, good.bins.list = NULL, rt.boundary = NULL){
  
  # mzdata_points<-unique(unlist(ms_list))  
  
  # Slide window to choose the high abundance bins in every districts
  
  MessageOutput("MS data Preparing...", "\n", NULL);
 
  #if(length(names(ms_list))>1000){
  #  ms_list_s<-sort(sample(names(ms_list),1000))
  #} else {
  #  ms_list_s<-sort(names(ms_list))
  #}

  # Update ms_lsit
  empty_toRemove <- which(sapply(names(ms_list),function(x){is.null(raw_data@assayData[[x]])}));
  if(length(empty_toRemove)){
    ms_list <- ms_list[-empty_toRemove];
  }
  
  if (is.null(good.bins.list)) {
    spectra_mz <- unlist(lapply(MSnbase::spectra(raw_data),mz))
    
    highestmz<-max(spectra_mz)
    lowestmz<-min(spectra_mz)
    
    # Split MS into 5 large districts
    bins.boundary<-seq(from=lowestmz, to= highestmz,length.out = 5)
    
    
    if (length(spectra_mz)>1000000){
      rannum <- seq(from = 1, 
                    to = length(spectra_mz), 
                    by = 50)
    } else if (length(spectra_mz)>100000){
      rannum <- seq(from = 1, 
                    to = length(spectra_mz), 
                    by = 5)
    } else {
      rannum<-seq(length(spectra_mz))
    }
    
    spectra_abundance <- unlist(lapply(MSnbase::spectra(raw_data),intensity))[rannum]
    spectra_mz <- spectra_mz[rannum]
    
    spectra_mz_set <- sapply(seq_len(4),FUN=function(ii){
      
      spectra_mz[spectra_mz > bins.boundary[ii] & spectra_mz <  bins.boundary[ii+1]]
      
    })
    
    spectra_abundance_set <- sapply(seq_len(4),FUN=function(ii){
      
      spectra_abundance[spectra_mz > bins.boundary[ii] &  spectra_mz  <  bins.boundary[ii+1]]
      
    })
    
    rm(spectra_abundance); rm(spectra_mz)
    
    good.bins.list<-bplapply(seq_len(4),
                             FUN = function(i, spectra_abundance_set, spectra_mz_set, bins.boundary){
                               
                               binsb.low<-bins.boundary[i];                             
                               binsb.high<-bins.boundary[i+1];                             
                               w<-1;                             
                               bins.width<-(binsb.high-binsb.low)*0.1;                             
                               wind.up<-wind.down<-0;
                               inten.sum.set<-numeric();
                               
                               while (!wind.up>binsb.high) {                               
                                 wind.down<-binsb.low+(w-1);                               
                                 wind.up<-binsb.low+bins.width+(w-1);                               
                                 
                                 inten.sum.set<-c(inten.sum.set,
                                                  sum (spectra_abundance_set[[i]]
                                                       [spectra_mz_set[[i]] > wind.down & spectra_mz_set[[i]] < wind.up]));
                                 w<-w+1;
                               }
                               
                               selected.bin.down<-binsb.low+which(max(inten.sum.set)==inten.sum.set)-1;
                               selected.bin.up<-binsb.low+which(max(inten.sum.set)==inten.sum.set)-1+bins.width
                               
                               return(list(selected.bin.down,selected.bin.up))
                               
                             },
                             spectra_abundance_set=spectra_abundance_set,
                             spectra_mz_set=spectra_mz_set,
                             bins.boundary=bins.boundary,
                             BPPARAM = MulticoreParam(4L))
  }
  
  MessageOutput("MS Data ready !", "\n", NULL);
  
  # Remove the peaks out of the bins
  MessageOutput("Identifying ROIs in m/z dimensions...", "\n", NULL);

  pb <- progress_bar$new(format = "ROIs Identification in MZ Dimension [:bar] :percent Time left: :eta", total = length(ms_list), clear = TRUE, width= 80)
  
  for (i in seq_along(ms_list)){
    pb$tick();
      ms.set<-raw_data@assayData[[names(ms_list)[i]]]@mz
      k<-which(sapply(ms.set,FUN = function(x){(x > good.bins.list[[1]][[1]] && x < good.bins.list[[1]][[2]]) | 
          (x > good.bins.list[[2]][[1]] && x < good.bins.list[[2]][[2]]) | 
          (x > good.bins.list[[3]][[1]] && x < good.bins.list[[3]][[2]]) | 
          (x > good.bins.list[[4]][[1]] && x < good.bins.list[[4]][[2]])}))
    
      raw_data@assayData[[names(ms_list)[i]]]@mz<-raw_data@assayData[[names(ms_list)[i]]]@mz[k];
      raw_data@assayData[[names(ms_list)[i]]]@intensity<-raw_data@assayData[[names(ms_list)[i]]]@intensity[k];
      raw_data@assayData[[names(ms_list)[i]]]@tic<-sum(raw_data@assayData[[names(ms_list)[i]]]@intensity);
      raw_data@assayData[[names(ms_list)[i]]]@peaksCount<-length(raw_data@assayData[[names(ms_list)[i]]]@mz)
  }
  MessageOutput("Identifying ROIs in m/z dimensions Done !", "\n", NULL);
  
  ## Trimed data outside the RT bins
  if (is.null(rt.boundary)) {
    rt_set<-sapply(names(ms_list),FUN=function(x) raw_data@assayData[[x]]@rt)
    rt.bin.width<-(max(rt_set)-min(rt_set))*rt.idx
    
    w<-0;rt.window.min<-min(rt_set);tic.sum<-numeric()
    
    while (!rt.window.min>max(rt_set)) {
      rt.window.min<-min(rt_set)+w*0.75
      rt.window.max<-min(rt_set)+rt.bin.width+w*0.75
      rt.name.set<-names(which(sapply(rt_set,FUN = function(x){x>rt.window.min && x<rt.window.max})))
      
      if(!identical(rt.name.set,character(0))){
        tic.sum<-c(tic.sum,sum(sapply(rt.name.set,FUN=function(x){raw_data@assayData[[x]]@tic})))
      }
      
      w<-w+1
    }
    
    rt.boundary.lowlimit<-min(rt_set)+which(max(tic.sum)==tic.sum)[ceiling(length(which(max(tic.sum)==tic.sum))/2)]*0.75
    rt.boundary.uplimit<-min(rt_set)+which(max(tic.sum)==tic.sum)[ceiling(length(which(max(tic.sum)==tic.sum))/2)]*0.75+rt.bin.width
  } else {
    rt.boundary.lowlimit <- rt.boundary[1]
    rt.boundary.uplimit <- rt.boundary[2]
  }
  
  MessageOutput("Identifying ROIs in RT dimensions...", "\n", NULL);
 
  pb <- progress_bar$new(format = "ROIs Identification in RT Dimension [:bar] :percent Time left: :eta", total = length(ms_list), clear = TRUE, width= 80)
  ncount<-numeric();
  
  for (w in seq_along(ms_list)){
    pb$tick();
    if (raw_data@assayData[[names(ms_list)[w]]]@rt>rt.boundary.lowlimit && raw_data@assayData[[names(ms_list)[w]]]@rt<rt.boundary.uplimit){
      ncount<-c(ncount,w)
    }
  }
  
  for (j in seq_along(ms_list)){
    if (!(j %in% ncount)){
      raw_data@assayData[[names(ms_list)[j]]]@mz<-raw_data@assayData[[names(ms_list)[j]]]@intensity<-as.double();
      raw_data@assayData[[names(ms_list)[j]]]@peaksCount<-as.integer(0);
      raw_data@assayData[[names(ms_list)[j]]]@tic<-as.double(0);
    }
  }
  
  return(raw_data)
}

#' @title Data trimming Method Based on Random MS
#' @description Trim raw data scan signal randomly in the mz dimension.
#' @param raw_data MSnExp object, the raw data that has been read in memory.
#' @param ms_list List, the names list of all scans.
#' @noRd
#' @import progress
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)

mz.trim_random <- function(raw_data, ms_list){
  
  mzdata_points <- unique(unlist(ms_list))
  highestmz <- max(mzdata_points)
  lowestmz <- min(mzdata_points)
  
  # Randomly select 100 mz center points
  #set.seed(1233);
  mzdata_points <- sample(unique(mzdata_points),size = 100)
  
  pb <- progress_bar$new(format = "Data Trimming [:bar] :percent Time left: :eta", total = length(raw_data@assayData), clear = TRUE, width= 60)
  
  for (i in seq_along(raw_data@assayData)){
    pb$tick();
    mz.data<-ms_list[[i]];
    remove.num.set<-numeric();k<-1
    
    for (j in seq_along(mz.data)){
      if (!(TRUE %in% (abs(mz.data[j]-mzdata_points)<=(highestmz-lowestmz)/10000))){
        remove.num.set[k]<-j;
        k<-k+1
      }
    }
    raw_data@assayData[[names(ms_list[i])]]@mz<-mz.data[-remove.num.set];
    raw_data@assayData[[names(ms_list[i])]]@intensity<-raw_data@assayData[[names(ms_list[i])]]@intensity[-remove.num.set];
    raw_data@assayData[[names(ms_list[i])]]@tic<-sum(raw_data@assayData[[names(ms_list[i])]]@intensity);
    raw_data@assayData[[names(ms_list[i])]]@peaksCount<-length(raw_data@assayData[[names(ms_list[i])]]@mz);
  }
  return(raw_data)
}

#' @title Data trimming Method Based on Random RT
#' @description Trim raw data scan signal randomly in the RT dimension.
#' @param raw_data MSnExp object, the raw data that has been read in memory.
#' @param ms_list List, the names list of all scans.
#' @noRd
#' @import progress
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)

rt.trim_random <- function(raw_data, ms_list){
  
  rt_set<-sapply(names(ms_list),FUN=function(x) raw_data@assayData[[x]]@rt)
  rt_begin<-min(rt_set);
  rt_end<-max(rt_set);
  
  xn2<-character();xn3<-xn<-xcount<-ncount<-rt_var_sum<-numeric();
  rt_width_aver<-(rt_end-rt_begin)/200
  xn_list<-list()
  
  pb <- progress_bar$new(format = "Data Trimming [:bar] :percent Time left: :eta", total = 200, clear = TRUE, width= 60)
  
  for (w in seq_len(200)){
    pb$tick();
    
    rt_var<-rt_begin+rt_width_aver*(w-0.5);
    rt_var_sum<-c(rt_var_sum,rt_var);
    
    xn_list[w][[1]]<-0
  }
  
  xn<-sapply(names(ms_list),FUN= function(x) which(abs(raw_data@assayData[[x]]@rt-rt_var_sum)<=(rt_width_aver/2+0.000001)))
  
  for (i in seq_along(xn)){xn_list[[xn[i]]]<-sum(xn_list[[xn[i]]],raw_data@assayData[[names(xn[i])]]@peaksCount)}
  
  rt_peakcounts<-unlist(xn_list)
  #rtslots<-rt_var_sum[sapply(tail(sort(unlist(xn_list)),10),FUN = function(x){which(x==rt_peakcounts)})]
  
  # randomly select 10 RT slots
  #set.seed(1244)
  rtslots<-rt_var_sum[sapply(sample(unlist(xn_list),10),
                             FUN = function(x){which(x==rt_peakcounts)})]
  
  for (j in seq_along(ms_list)){
    
    if (TRUE %in% (abs(raw_data@assayData[[names(ms_list)[j]]]@rt-rtslots)<=(rt_width_aver/2+0.000001))){
      xn2<-c(xn2,names(ms_list)[j]);xn3<-c(xn3,j)
    }
  }
  
  for (k in seq_along(ms_list))  {
    
    if (!(k %in% xn3)){
      raw_data@assayData[[names(ms_list)[k]]]@peaksCount<-as.integer(0);
      raw_data@assayData[[names(ms_list)[k]]]@mz<-as.double();
      raw_data@assayData[[names(ms_list)[k]]]@intensity<-as.double();
      raw_data@assayData[[names(ms_list)[k]]]@tic<-as.double(0);
    }
  }
  return(raw_data)
}

#' @title Data trimming Method Based on Specific MS
#' @description Trim data based on specific mz values. Positive values will be specially retained, 
#' while the negative values will be removed.
#' @param raw_data MSnExp object, the raw data that has been read in memory.
#' @param ms_list List, the names list of all scans.
#' @param mz Numeric, the specifric mz value that will be kept or removed.
#' @param mzdiff Numeric, the deviation (ppm) for the 'mz' values. Default is 100.
#' @noRd
#' @import progress
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)

mz.trim_specific<-function(raw_data, ms_list, mz, mzdiff=100){
  
  #mz<-c(72.323,100,120,240,360,480,720,780.524,1080);
  if (missing(mz)){
    stop("'mz' must be provided for mz_specific mode as a numeric vector!")
  }
  
  if (missing(mzdiff)){
    stop("'mzdiff' must be provided for mz_specific mode as a numeric value!")
  }
  
  if(length(mz) == 0){
    return(raw_data)
  }
  
  if (mz[1] == 0 | mzdiff==0){
    stop("'mz' or 'mzdiff' cannot be zero!")
  }
  
  mz.neg<-mz[which(mz<0)];
  mz.pos<-mz[which(mz>0)];
  
  if (!identical(mz.pos,numeric(0))){
    pb <- progress_bar$new(format = "Data Trimming_keeping [:bar] :percent Time left: :eta", total = length(ms_list), clear = TRUE, width= 60)
    
    for (w in seq_along(ms_list)){
      pb$tick();
      if(is.null(raw_data@assayData[[names(ms_list)[w]]])){next}
      xn<-unlist(sapply(mz.pos,FUN=function(x){which((abs(x-raw_data@assayData[[names(ms_list)[w]]]@mz)/(x)*1000000)<=mzdiff)}));
      
      raw_data@assayData[[names(ms_list)[w]]]@mz<-raw_data@assayData[[names(ms_list)[w]]]@mz[xn];
      raw_data@assayData[[names(ms_list)[w]]]@intensity<-raw_data@assayData[[names(ms_list)[w]]]@intensity[xn];
      raw_data@assayData[[names(ms_list)[w]]]@peaksCount<-length(raw_data@assayData[[names(ms_list)[w]]]@mz);
      raw_data@assayData[[names(ms_list)[w]]]@tic<-sum(raw_data@assayData[[names(ms_list)[w]]]@intensity);
    }
    
    for (w in seq_along(ms_list)){
      if (identical(raw_data@assayData[[names(ms_list)[w]]]@mz,numeric(0))){
        raw_data@assayData[[names(ms_list)[w]]]@mz<-as.double();
        raw_data@assayData[[names(ms_list)[w]]]@intensity<-as.double();
      }
    }
  } 
  
  if (!identical(mz.neg,numeric(0))){
    
    pb <- progress_bar$new(format = "Data Trimming_removing [:bar] :percent Time left: :eta", total = length(ms_list), clear = TRUE, width= 60)
    
    for (w in seq_along(ms_list)){
      pb$tick();
      if(is.null(raw_data@assayData[[names(ms_list)[w]]])){next}
      xn<-unlist(sapply(abs(mz.neg),FUN=function(x){which((abs(x-raw_data@assayData[[names(ms_list)[w]]]@mz)/(x)*1000000)<=mzdiff)}));
      
      if (!identical(xn, numeric(0))){
        raw_data@assayData[[names(ms_list)[w]]]@mz<-raw_data@assayData[[names(ms_list)[w]]]@mz[-xn];
        raw_data@assayData[[names(ms_list)[w]]]@intensity<-raw_data@assayData[[names(ms_list)[w]]]@intensity[-xn];
        raw_data@assayData[[names(ms_list)[w]]]@peaksCount<-length(raw_data@assayData[[names(ms_list)[w]]]@mz);
        raw_data@assayData[[names(ms_list)[w]]]@tic<-sum(raw_data@assayData[[names(ms_list)[w]]]@intensity);
      }
    }
    
    for (w in seq_along(ms_list)){
      if(is.null(raw_data@assayData[[names(ms_list)[w]]])){next}
      if (identical(raw_data@assayData[[names(ms_list)[w]]]@mz,numeric(0))){
        raw_data@assayData[[names(ms_list)[w]]]@mz<-as.double();
        raw_data@assayData[[names(ms_list)[w]]]@intensity<-as.double();
      }
    }
  }
  return(raw_data)
}

#' @title Data trimming Method Based on Specific RT
#' @description Trim data based on specific RT values. Positive values will be specially retained, 
#' while the negative values will be removed.
#' @param raw_data MSnExp object, the raw data that has been read in memory.
#' @param ms_list List, the names list of all scans.
#' @param rt Numeric, the specifric RT value that will be kept or removed.
#' @param rtdiff Numeric, the deviation (ppm) for the 'rt' values. Default is 100.
#' @noRd
#' @import progress
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)

rt.trim_specific<-function(raw_data,ms_list,rt,rtdiff=10){
  
  if (missing(rt)){
    stop("'rt' must be provided for rt_specific mode with a vector !")
  }
  
  if (missing(rtdiff)){
    stop("'rtdiff' must be provided for rt_specific mode with a numeric !")
  }
  
  if (rt[1] == 0 | rtdiff==0){
    stop("'rt' or 'rtdiff' can not be zero !")
  }
  
  if (rt[1] > 0){
    pb <- progress_bar$new(format = "Data Trimming [:bar] :percent Time left: :eta", total = length(ms_list), clear = TRUE, width= 60)
    
    ncount<-numeric();
    for (w in seq_along(ms_list)){
      pb$tick();
      if (TRUE %in% (abs(raw_data@assayData[[names(ms_list)[w]]]@rt-rt)<=rtdiff)){
        ncount<-c(ncount,w)
      }
    }
    
    for (j in seq_along(ms_list)){
      
      if (!(j %in% ncount)){
        raw_data@assayData[[names(ms_list)[j]]]@mz<-raw_data@assayData[[names(ms_list)[j]]]@intensity<-as.double();
        raw_data@assayData[[names(ms_list)[j]]]@peaksCount<-as.integer(0);
        raw_data@assayData[[names(ms_list)[j]]]@tic<-as.double(0);
      }
      
    }
  } else {
    
    pb <- progress_bar$new(format = "Data Trimming [:bar] :percent Time left: :eta", total = length(ms_list), clear = TRUE, width= 60)
    
    ncount<-numeric();
    for (w in seq_along(ms_list)){
      pb$tick();
      if (TRUE %in% (abs(raw_data@assayData[[names(ms_list)[w]]]@rt-abs(rt))<=rtdiff)){
        ncount<-c(ncount,w)
      }
      
    }
    for (j in seq_along(ms_list)){
      
      if (j %in% ncount){
        raw_data@assayData[[names(ms_list)[j]]]@mz<-raw_data@assayData[[names(ms_list)[j]]]@intensity<-as.double();
        raw_data@assayData[[names(ms_list)[j]]]@peaksCount<-as.integer(0);
        raw_data@assayData[[names(ms_list)[j]]]@tic<-as.double(0);
      }
    }
  }
  return(raw_data)
}

#' @title Function for 'Empty scan' removal
#' @description Function for 'Empty scan' removal (internal use only)
#' @noRd
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}
.emptyscan.remove<-function(raw_data, ms_list){
  
  name.set<-sort(names(raw_data@assayData))
  
  sample.set<-unique(sapply(seq_along(name.set),FUN=function(x){sub(".[S][0-9]+","",name.set[x])}))
  assayData<-list();assayData1<-new.env();ki<-k<-1
  
  df.data0<-raw_data@featureData@data;
  df.data<-data.frame();
  
  if(length(sample.set) > 9){
    nchar.num<-nchar(names(ms_list)[1])-5
  } else if (length(sample.set) > 99) {
    nchar.num<-nchar(names(ms_list)[1])-6
  } else {
    nchar.num<-nchar(names(ms_list)[1])-4
  }

  suppressWarnings(for (i in seq_along(raw_data@assayData)){
    
    if (!raw_data@assayData[[name.set[i]]]@tic==0){
      
      if (!k==1){
        if (!sample.set[which(sub(".[S][0-9]+","",name.set[i])==sample.set)]==sub(".[S][0-9]+","",names(assayData)[ki-1])){
          k<-1
        } 
      }
      
      name0<-paste0(sample.set[which(sub(".[S][0-9]+","",name.set[i])==sample.set)],
                    ".S",formatC(k, width=nchar.num, digits = nchar.num, flag="0"))
      
      assayData[name0]<-raw_data@assayData[[name.set[i]]];
      assayData[[name0]]@acquisitionNum<-as.integer(k);
      assayData[[name0]]@scanIndex<-as.integer(k);
      
      df.data[ki,1]<-k;
      row.names(df.data)[ki]<-name0;
      
      k<-k+1;ki<-ki+1;
    }
  })
  
  list2env(assayData,envir = assayData1)
  raw_data@assayData<-assayData1
  names(df.data)<-"spectrum"
  
  metaData<-raw_data@featureData@varMetadata;
  featureData<-AnnotatedDataFrame(data=df.data, varMetadata=metaData);
  raw_data@featureData<-featureData;
  
  return(raw_data)
}

#' @title Function MS Generation
#' @description Output the MS data. This function will generate .mzML MS data in the working dirctory.
#' @param raw_data MS data in R environment with "MSnExp" class.
#' @noRd
#' @import MSnbase
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
PerformMSDataOutput<-function(raw_data){
  message("Data Writing...")
  writenames<-paste0("new_",raw_data@phenoData@data[["sample_name"]],sep = "")
  #dir.create(paste0(datapath,"/trimed",collapse = ""))
  suppressMessages(writeMSData(raw_data, writenames, outformat = "mzML"))
  message("Output Data:","\n");
  message(writenames)
  message(c("Depositing Folder: ",getwd()))
  message("Data Writing Finished !")
}

#' @title Function for 3D ms plotting
#' @description Function for 3D ms plotting (internal use only)
#' @importFrom lattice cloud
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}
#' @noRd
plot.MS_3D<-function(object) {
  
  # if (.on.public.web()){load_lattice()};
  # 
  dd <- as(object, "data.frame")
  dd$rt <- dd$rt*60;
  ms <- NULL ## get rid of 'no visible global function definition' note
  par.set <- list(box.3d = list(lwd=.2))
  
  cloud(intensity ~ mz + rt , data = dd,
        type="h",
        scales= list(
          arrows=FALSE,
          cex=.65,
          draw=TRUE),
        aspect=c(.8, 1),
        group = ms,
        zoom = 1, 
        par.settings = par.set,
        axis.line = list(col = "transparent"),
        xlab="m/z", ylab="Retention time (sec)", zlab=NULL)
  
}

#' @noRd
ContaminatsRemoval <- function(raw_data, ms_list){
  
  scan_names <- sort(names(raw_data@assayData));
  fre_stats <- plyr::count(round(unname(unlist(ms_list)), 2))
  topScan_stats <- fre_stats[order(fre_stats$freq, decreasing = TRUE),][seq_len(300),] # get top 300 most frequent mz value
  scan_totalLength <- length(raw_data@assayData)
  top_mzs <- topScan_stats$x
  topScan_stats[,4] <- topScan_stats[,3] <- 0;
  colnames(topScan_stats)[c(3,4)] <- c("RT_ratio", "mz_mean");
  
  MessageOutput("Identifying regions of potential contaminants", "", NULL);
  
  k_count <- 0;
  
  for(i in seq_len(80)){
    
    mzrange <- c(top_mzs[i]-0.005, top_mzs[i]+0.005);
    RawData <- suppressWarnings(MSnbase::filterMz(raw_data, mzrange));
    scan_count <- 0;
    mz_vec <- vector();
    
    if ((i - k_count) > 0 ){
      MessageOutput(".", "", NULL);
      k_count <- k_count +10;
    }
    
    for(j in scan_names){
      if(!identical(RawData@assayData[[j]]@mz, numeric(0))){
        scan_count <- scan_count +1;
        mz_vec <- c(mz_vec, RawData@assayData[[j]]@mz)
      }
    }
    
    topScan_stats[i,3] <- scan_count/scan_totalLength;
    topScan_stats[i,4] <- mean(mz_vec);
    
    if(topScan_stats[i,3] < 0.5){
      break;
    }
    
  }
  
  MessageOutput("Done!", "\n", NULL)
  #save(topScan_stats, file = "topScan_stats.rda")
  mzs_toRemove <- topScan_stats[topScan_stats[,3] > 0.5, 4]
  
  raw_data_clean <- mz.trim_specific(raw_data, ms_list, -mzs_toRemove, mzdiff=10)
  MessageOutput(paste0(length(mzs_toRemove), " potential contaminamts will not be used for parameters optimization !\nGoing to the next step..."), 
                "\n",
                NULL)
  
  raw_data_clean <- .emptyscan.remove(raw_data_clean, ms_list);
  
  return(raw_data_clean);
  
}
