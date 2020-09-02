#' Set class information for MS data
#' @description This function sets the class information
#' for preprocessing MS data.
#' @author Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' Mai Yamamoto \email{yamamoto.mai@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
SetClass <- function(class) {
  groupInfo <<- class
}

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
        sum(x@intensity, na.rm = T)
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
  
  require(scales)
  
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

#' Perform peak profiling
#' This function performs feature extraction of user's raw MS data using
#' the rawData object created using the ImportRawMSData function.
#' @param rawData The object created using the ImportRawMSData function,
#' containing the raw MS data.
#' @param Params The object created using the SetPeakParam function,
#' containing user's specified or default parameters for downstream
#' raw MS data pre-processing.
#' @param plotSettings List, plotting parameters produced by SetPlotParam Function.
#' Defaut is set to true.
#' @param ncore Numeric, used to define the cores' number for Peak Profiling.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' Mai Yamamoto \email{yamamoto.mai@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import MSnbase
#' @import BiocParallel
#' @import ggplot2


PerformPeakProfiling <-
  function(mSet,
           Params,
           plotSettings,
           ncore,
           running.controller = NULL) {
    
    require(BiocParallel)
    
    #Build Running plan for data import - Indentify the controller
    function.name <- "peak_profiling";
    
    if (is.null(running.controller)) {
      c1 <- c2 <- c3 <- c4 <- T
    } else {
      c1 <-
        running.controller[["peak_profiling"]][["c1"]] # used to control peak picking
      c2 <-
        running.controller[["peak_profiling"]][["c2"]] # used to control peak alignment
      c3 <-
        running.controller[["peak_profiling"]][["c3"]] # used to control peak filing
      c4 <-
        running.controller[["peak_profiling"]][["c4"]] # used to control plotting
    }
    
    ### Update parameters' style
    write.table(
      unlist(Params),
      file = "param_optimized.txt",
      row.names = T,
      col.names = F,
      quote = F
    )
    
    mSet@params <- updateRawSpectraParam (Params)
    
    ### Setting the different parallel method for linux or windows
    MessageOutput(mes = NULL,
                  ecol = NULL,
                  progress = 25)
    
    if (.on.public.web) {
      load_biocparallel()
      total_threads <- 4
      
    } else if (missing(ncore)) {
      total_threads <- detectCores() * 2 / 3
    } else {
      total_threads <- ncore
    }
    
    if (.Platform$OS.type == "unix") {
      register(bpstart(MulticoreParam(ceiling(total_threads))))
    } else if (.Platform$OS.type == "windows") {
      register(bpstart(SnowParam(ceiling(total_threads))))
    }
    
    MessageOutput(
      mes = paste0(
        ceiling(total_threads),
        " CPU Threads will be used for peak profiling !"
      ),
      ecol = "\n",
      progress = NULL
    )
    
    #   ---------===========----- I. Peak picking -----===========------------
    
    if (c1) {
      MessageOutput(mes = "Step 3/6: Started peak picking! This step will take some time...",
                    ecol = "\n",
                    progress = NULL);
      
      mSet <-
        tryCatch(
          PerformPeakPicking(mSet),
          error = function(e) {
            e
          }
        )
      
      gc()
      
      if (running.as.plan & class(mSet)[1] != "simpleError") {
        cache.save(mSet, paste0(function.name, "_c1"))
        marker_record(paste0(function.name, "_c1"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c1")
      marker_record("raw_data_samples_c2")
    }
    
    if (.on.public.web) {
      if (class(mSet)[1] == "simpleError") {
        MessageOutput(
          mes =  paste0(
            "<font color=\"red\">",
            "\nERROR:",
            mSet$message,
            "</font>"
          ),
          ecol = "\n",
          progress = 50
        )
        stop("EXCEPTION POINT CODE: PU1")
      }
      
      MessageOutput(
        mes =  paste0("Step 3/6: Peak picking finished ! (", Sys.time(), ")"),
        ecol = "\n",
        progress = 50
      )
    }
    
    #   --------===========----- II. Peak alignment -----===========------------
    if (c2) {
      MessageOutput(
        mes = paste("Step 4/6: Started peak alignment! This step is running..."),
        ecol = "\n",
        progress = 50.1
      )
      
      mSet <-
        tryCatch(
          PerformPeakAlignment(mSet),
          error = function(e) {
            e
          }
        )
      
      gc()
      
      if (running.as.plan & class(mSet)[1] != "simpleError") {
        cache.save(mSet, paste0(function.name, "_c2"))
        marker_record(paste0(function.name, "_c2"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c2")
      marker_record("raw_data_samples_c2")
    }
    
    if (.on.public.web) {

      if (class(mSet)[1] == "simpleError") {
        MessageOutput(
          mes = paste0(
            "<font color=\"red\">",
            "\nERROR:",
            mSet$message,
            "</font>"
          ),
          ecol = "\n",
          progress = 73
        )
        stop("EXCEPTION POINT CODE: PU2")
      }
      
      MessageOutput(
        mes = paste0("Step 4/6: Peak alignment finished ! (", Sys.time(), ")"),
        ecol = "\n",
        progress = 73
      )
    }
    
    #   --------===========----- III. Peak filling -----===========------------
    if (c3) {
      MessageOutput(
        mes = paste(
          "Step 5/6: Started peak filling! This step may take some time..."
        ),
        ecol = "\n",
        progress = NULL
      )
      
      mSet <-
        tryCatch(
          PerformPeakFiling (mSet),
          error = function(e) {
            e
          }
        )
      
      gc()
      
      if (running.as.plan & class(mSet)[1] != "simpleError") {
        cache.save(mSet, paste0(function.name, "_c3"))
        marker_record(paste0(function.name, "_c3"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c3")
      marker_record("raw_data_samples_c2")
    }
    
    if (class(mSet)[1] == "simpleError") {
      MessageOutput(
        mes = paste0(
          "<font color=\"red\">",
          "\nERROR:",
          mSet$message,
          "</font>"
        ),
        ecol = "\n",
        progress = 88
      )
      
      stop("EXCEPTION POINT CODE: PU3")
    }
    
    MessageOutput(
      mes = paste0(
        "Step 5/6: Peak filing finished ! (",
        Sys.time(),
        ")",
        "\nPeak Profiling finished successfully !"
      ),
      ecol = "\n",
      progress = 88
    )
    
    save(mSet, file = "mSet.rda")
    
    MessageOutput(
      mes = paste0("Begin to plotting figures..."),
      ecol = "\n",
      progress = 89
    )
    
    #  ---------------====------IV. Plotting Results --------========-----------
    if (c4) {
      sample_idx <- mSet[["onDiskData"]]@phenoData@data[["sample_group"]]
      
      if (missing(plotSettings)) {
        plotSettings <- SetPlotParam(
          name_peak_in = "Peak_Intensity",
          name_PCA = "PCA",
          name_adj_RT = "Adjusted_RT",
          name_adj_BPI = "Adjusted_BPI"
        )
      } else {
        plotSettings$name_peak_in = "Peak_Intensity"
        plotSettings$name_PCA = "PCA"
        plotSettings$name_adj_RT = "Adjusted_RT"
        plotSettings$name_adj_BPI = "Adjusted_BPI"
      }
      
      if (plotSettings$Plot == T) {
        ### 1. Peak Intensity plotting -----
        PlotSpectraInsensityStistics(
          mSet,
          paste0(plotSettings$name_peak_in, ".", plotSettings$format),
          plotSettings$format,
          plotSettings$dpi,
          8
        )
        
        ### 2. PCA plotting -----
        if (.on.public.web) {
          load_ggplot()
        }
        
        PlotSpectraPCA(
          mSet,
          paste0(plotSettings$name_PCA, ".", plotSettings$format),
          plotSettings$format,
          plotSettings$dpi,
          8
        )
        
        ### 3. Adjusted RT plotting -----
        PlotSpectraRTadj(
          mSet,
          paste0(plotSettings$name_adj_RT, ".", plotSettings$format),
          plotSettings$format,
          plotSettings$dpi,
          8
        )
        
        ### 4. Chromatogram Generation -----
        PlotSpectraBPIadj(
          mSet,
          paste0(plotSettings$name_adj_BPI, ".", plotSettings$format),
          plotSettings$format,
          plotSettings$dpi,
          plotSettings$width
        )
      }
      
      marker_record(paste0(function.name, "_c4"))
      
    }
    
    MessageOutput(mes = NULL,
                  ecol = NULL,
                  progress = 90)
    
    return(mSet)
  }

#' Set annotation parameters
#' @description This function sets the parameters for peak annotation.
#' @param polarity Character, specify the polarity of the MS instrument. Either
#' "negative" or "positive".
#' @param perc_fwhm Numeric, set the percentage of the width of the FWHM for peak grouping.
#' Default is set to 0.6.
#' @param mz_abs_iso Numeric, set the allowed variance for the search (for isotope annotation).
#' The default is set to 0.005.
#' @param max_charge Numeric, set the maximum number of the isotope charge. For example,
#' the default is 2, therefore the max isotope charge is 2+/-.
#' @param max_iso Numeric, set the maximum number of isotope peaks. For example, the default
#' is 2, therefore the max number of isotopes per peaks is 2.
#' @param corr_eic_th Numeric, set the threshold for intensity correlations across samples.
#' Default is set to 0.85.
#' @param mz_abs_add Numeric, set the allowed variance for the search (for adduct annotation).
#' The default is set to 0.001.
#' @author Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' Mai Yamamoto \email{yamamoto.mai@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
SetAnnotationParam <-
  function(polarity = "positive",
           perc_fwhm = 0.6,
           mz_abs_iso = 0.005,
           max_charge = 2,
           max_iso = 2,
           corr_eic_th = 0.85,
           mz_abs_add = 0.001) {
    annParams <- list()
    
    if (.on.public.web) {
      load("params.rda")
      
      annParams$polarity <- peakParams$polarity
      annParams$perf.whm <- peakParams$perc_fwhm
      annParams$mz.abs.iso <- peakParams$mz_abs_iso
      annParams$max.charge <- peakParams$max_charge
      annParams$max.iso <- peakParams$max_iso
      annParams$corr.eic.th <- peakParams$corr_eic_th
      annParams$mz.abs.add <- peakParams$mz_abs_add
      
    } else {
      annParams$polarity <- polarity
      annParams$perf.whm <- perc_fwhm
      annParams$mz.abs.iso <- mz_abs_iso
      annParams$max.charge <- max_charge
      annParams$max.iso <- max_iso
      annParams$corr.eic.th <- corr_eic_th
      annParams$mz.abs.add <- mz_abs_add
      
    }
    
    return(annParams)
  }

#' Perform peak annotation
#' @description This function performs peak annotation on
#' the xset object created using the PerformPeakPicking function.
#' @param xset The object created using the PerformPeakPicking function,
#' containing the peak picked MS data.
#' @param annParams The object created using the SetAnnotationParam function,
#' containing user's specified or default parameters for downstream
#' raw MS data pre-processing.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' Mai Yamamoto \email{yamamoto.mai@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import MSnbase
#' @importFrom graph ftM2graphNEL
#' @importFrom RBGL highlyConnSG
#' @references Kuhl C, Tautenhahn R, Boettcher C, Larson TR, Neumann S (2012).
#' "CAMERA: an integrated strategy for compound spectra extraction and annotation of
#' liquid chromatography/mass spectrometry data sets." Analytical Chemistry, 84, 283-289.
#' http://pubs.acs.org/doi/abs/10.1021/ac202450g.
PerformPeakAnnotation <-
  function(mSet,
           annotaParam,
           ncore = 1,
           running.controller = NULL) {
    MessageOutput(
      mes = paste0("Step 6/6: Starting Peak Annotation..."),
      ecol = "\n",
      progress = 91
    )
    
    if (.on.public.web) {
      dyn.load(.getDynLoadPath())
    }
    
    if (ncore > 1) {
      print("Only single core mode is supported now. Parallel will be supported later !")
      ncore <- 1
    }
    
    load_progress()
    load_graph()
    load_RBGL()
    
    function.name <- "peak_annotation"
    
    if (is.null(running.controller) | plan_count == 1) {
      operators_4 <- T
    } else {
      operators_4 <- running.controller[["operators"]][["operators_4"]]
    }
    
    if (operators_4) {
      ## 1. Prepare the Annotation Object-------
      xs <- mSet$xcmsSet
      
      if (is.null(xs)) {
        stop("No xcmsSet object in 'mSet' was given !")
      } else if (!class(xs) == "xcmsSet") {
        stop("There is correct xcmsSet object in mSet !")
      }
      
      mSet$AnnotateObject <- list()
      
      if (length(xs@phenoData[["sample_name"]]) > 1 &&
          !nrow(xs@groups) > 0) {
        stop ('No group information found or contain only one sample.')
      }
      mSet$AnnotateObject$sample   <-  as.numeric(NA)
      mSet$AnnotateObject$groupInfo <- getPeaks_selection(xs)
      runParallel <- list()
      runParallel$enable   <-  0
      
      if (ncore > 1) {
        ## If MPI is available ...
        rmpi = "Rmpi"
        opt.warn <- options("warn")$warn
        options("warn" = -1)
        if ((Sys.info()["sysname"] != "Windows") &&
            require(rmpi, character.only = TRUE) &&
            !is.null(ncore)) {
          if (is.loaded('mpi_initialize')) {
            #test if not already slaves are running!
            if (mpi.comm.size() > 0) {
              warning(
                "There are already intialized mpi slaves on your machine.\nCamera will try to uses them!\n"
              )
              
              runParallel$enable <- 1
              runParallel$mode <- rmpi
              
            } else{
              mpi.spawn.Rslaves(ncore = ncore, needlog = FALSE)
              if (mpi.comm.size() > 1) {
                #Slaves have successfull spawned
                runParallel$enable <- 1
                runParallel$mode <- rmpi
                
              } else{
                warning(
                  "Spawning of mpi slaves have failed. CAMERA will run without parallelization.\n"
                )
              }
            }
          } else {
            #And now??
            warning("DLL mpi_initialize is not loaded. Run single core mode!\n")
            
          }
        } else {
          #try local sockets using snow package
          snow = "snow"
          if (try(require(snow, character.only = TRUE, quietly = TRUE))) {
            cat("Starting snow cluster with",
                ncore,
                "local sockets.\n")
            snowclust <- makeCluster(ncore, type = "SOCK")
            runParallel$enable <- 1
            runParallel$mode <- snow
            runParallel$cluster <- snowclust
          }
        }
        
        options("warn" = opt.warn)
        
        MessageOutput(mes = "Run cleanParallel after processing to remove the spawned slave processes!",
                      ecol = "\n",
                      progress = NULL)
        
      }
      
      if (!is.null(annotaParam[["polarity"]])) {
        if (is.na(match.arg(annotaParam[["polarity"]], c("positive", "negative")))) {
          stop("Parameter polarity has to be 'positive' or 'negative' !")
        } else{
          mSet$AnnotateObject$polarity <- annotaParam[["polarity"]]
        }
      }
      
      mSet$AnnotateObject$runParallel <- runParallel
      mSet$AnnotateObject$annoID <- matrix(ncol = 4, nrow = 0)
      mSet$AnnotateObject$annoGrp <- matrix(ncol = 4, nrow = 0)
      mSet$AnnotateObject$isoID <- matrix(ncol = 4, nrow = 0)
      
      colnames(mSet$AnnotateObject$annoID) <-
        c("id", "grpID", "ruleID", "parentID")
      
      colnames(mSet$AnnotateObject$annoGrp) <-
        c("id", "mass", "ips", "psgrp")
      
      colnames(mSet$AnnotateObject$isoID)  <-
        c("mpeak", "isopeak", "iso", "charge")
      
      MessageOutput(mes = NULL,
                    ecol = "\n",
                    progress = 92)
      
      ## 2. Group peaks according to their retention time into pseudospectra-groups-----
      
      intval <- "maxo"
      perfwhm <- annotaParam$perf.whm
      sigma <- 6
      sample    <- mSet$AnnotateObject$sample
      pspectra  <- list()
      psSamples <- NA
      
      MessageOutput(
        mes = paste0("Start grouping after retention time."),
        ecol = "\n",
        progress = NULL
      )
      
      if (mSet$AnnotateObject$groupInfo[1, "rt"] == -1) {
        # Like FTICR Data
        warning("Warning: no retention times avaiable. Do nothing\n")
        return(invisible(mSet$AnnotateObject))
        
      } else{
        if (is.na(sample[1]) || length(mSet$xcmsSet@filepaths) > 1) {
          # grouped peaktable within automatic selection or sub selection
          if (is.na(sample[1])) {
            index <- 1:length(mSet$xcmsSet@filepaths)
          } else{
            index <- sample
          }
          
          gvals    <- groupval(mSet$xcmsSet)[, index, drop = FALSE]
          peakmat  <- mSet$xcmsSet@peaks
          groupmat <- mSet$xcmsSet@groups
          
          #calculate highest peaks
          maxo      <-
            as.numeric(apply(gvals, 1, function(x, peakmat) {
              val <- na.omit(peakmat[x, intval])
              if (length(val) == 0) {
                return(NA)
              } else{
                return(max(val))
              }
            }, peakmat))
          
          maxo[which(is.na(maxo))] <- -1
          maxo      <- cbind(1:length(maxo), maxo)
          
          #highest peak index
          int.max   <-
            as.numeric(apply(gvals, 1, function(x, peakmat) {
              which.max(peakmat[x, intval])
            }, peakmat))
          
          peakrange <- matrix(apply(gvals, 1, function(x, peakmat) {
            val <- peakmat[x, intval]
            
            if (length(na.omit(val)) == 0) {
              return(c(0, 1))
              
            } else {
              return(peakmat[x[which.max(val)], c("rtmin", "rtmax")])
            }
          }, peakmat),
          ncol = 2,
          byrow = TRUE)
          
          colnames(peakrange) <- c("rtmin", "rtmax")
          
          while (length(maxo) > 0) {
            iint   <- which.max(maxo[, 2])
            rtmed  <-
              groupmat[iint, "rtmed"]
            #highest peak in whole spectra
            rt.min <- peakrange[iint, "rtmin"]
            
            rt.max <-
              peakrange[iint, "rtmax"]
            #begin and end of the highest peak
            hwhm   <-
              ((rt.max - rt.min) / sigma * 2.35 * perfwhm) / 2
            #fwhm of the highest peak
            #all other peaks whose retensiontimes are in the fwhm of the highest peak
            irt    <-
              which(groupmat[, 'rtmed'] > (rtmed - hwhm) &
                      groupmat[, 'rtmed'] < (rtmed + hwhm))
            if (length(irt) > 0) {
              #if peaks are found
              idx <- maxo[irt, 1]
              
              pspectra[[length(pspectra) + 1]] <- idx
              #create groups
              psSamples[length(pspectra)]  <-
                index[int.max[maxo[iint, 1]]] # saves the sample of the peak which is in charge for this pspectrum
              maxo <-
                maxo[-irt, , drop = FALSE]
              #set itensities of peaks to NA, due to not to be found in the next cycle
              groupmat   <- groupmat[-irt, , drop = FALSE]
              peakrange  <- peakrange[-irt, , drop = FALSE]
              
            } else{
              idx <- maxo[iint, 1]
              cat(
                "Warning: Feature ",
                idx,
                " looks odd for at least one peak. Please check afterwards.\n"
              )
              
              pspectra[[length(pspectra) + 1]] <- idx
              #create groups
              psSamples[length(pspectra)]  <-
                index[int.max[maxo[iint, 1]]] # saves the sample of the peak which is in charge for this pspectrum
              maxo <-
                maxo[-iint, , drop = FALSE]
              #set itensities of peaks to NA, due to not to be found in the next cycle
              groupmat   <- groupmat[-iint, , drop = FALSE]
              peakrange  <- peakrange[-iint, , drop = FALSE]
              
            }
          }
          
        } else{
          #One sample experiment
          peakmat <- mSet$xcmsSet@peaks
          maxo    <- peakmat[, intval]
          #max intensities of all peaks
          maxo    <- cbind(1:length(maxo), maxo)
          
          while (length(maxo) > 0) {
            iint   <- which.max(maxo[, 2])
            
            rtmed  <-
              peakmat[iint, "rt"]
            #highest peak in whole spectra
            rt.min <- peakmat[iint, "rtmin"]
            
            rt.max <-
              peakmat[iint, "rtmax"]
            #begin and end of the highest peak
            hwhm   <-
              ((rt.max - rt.min) / sigma * 2.35 * perfwhm) / 2
            #fwhm of the highest peak
            #all other peaks whose retensiontimes are in the fwhm of the highest peak
            irt    <-
              which(peakmat[, 'rt'] > (rtmed - hwhm) &
                      peakmat[, 'rt'] < (rtmed + hwhm))
            if (length(irt) > 0) {
              #if peaks are found
              idx <- maxo[irt, 1]
              
              pspectra[[length(pspectra) + 1]] <- idx
              #create groups
              maxo <-
                maxo[-irt, , drop = FALSE]
              #set itensities of peaks to NA, due to not to be found in the next cycle
              peakmat <- peakmat[-irt, , drop = FALSE]
              
            } else{
              idx <- maxo[iint, 1]
              cat(
                "Warning: Feature ",
                idx,
                " looks odd for at least one peak. Please check afterwards.\n"
              )
              
              pspectra[[length(pspectra) + 1]] <- idx
              #create groups
              maxo <-
                maxo[-iint, , drop = FALSE]
              #set itensities of peaks to NA, due to not to be found in the next cycle
              peakmat  <- peakmat[-iint, , drop = FALSE]
            }
          }
          psSamples <- rep(sample, length(pspectra))
        }
        
        mSet$AnnotateObject$pspectra  <- pspectra
        mSet$AnnotateObject$psSamples <- psSamples
        
        MessageOutput(
          mes = paste(
            "Created",
            length(mSet$AnnotateObject$pspectra),
            "pseudospectra."
          ),
          ecol = "\n",
          progress = NULL
        )
      }
      MessageOutput(mes = NULL,
                    ecol = "\n",
                    progress = 93)
      
      ## 3. Annotate isotope peaks -----
      
      maxcharge <-
        annotaParam$max.charge
      maxiso <- annotaParam$max.iso
      
      mzabs <- annotaParam$mz.abs.add
      intval = c("maxo")
      
      minfrac = 0.8
      isotopeMatrix = NULL
      filter = TRUE
      ppm <- 5
      
      if (!is.wholenumber(maxcharge) || maxcharge < 1) {
        stop("Invalid argument 'maxcharge'. Must be integer and > 0.\n")
      }
      if (!is.wholenumber(maxiso) || maxiso < 1) {
        stop("Invalid argument 'maxiso'. Must be integer and > 0.\n")
      }
      if (!is.numeric(mzabs) || mzabs < 0) {
        stop("Invalid argument 'mzabs'. Must be numeric and not negative.\n")
      }
      
      #intval <- match.arg(intval)
      
      if (!is.null(isotopeMatrix)) {
        if (!is.matrix(isotopeMatrix) ||
            ncol(isotopeMatrix) != 4 || nrow(isotopeMatrix) < 1
            || !is.numeric(isotopeMatrix)) {
          stop("Invalid argument 'isotopeMatrix'. Must be four column numeric matrix.\n")
        } else {
          colnames(isotopeMatrix) <- c("mzmin", "mzmax", "intmin", "intmax")
        }
      } else if (maxiso > 8) {
        stop("Invalid argument 'maxiso'. Must be lower 9 or provide your own isotopeMatrix.\n")
      } else{
        isotopeMatrix <- calcIsotopeMatrix(maxiso = maxiso)
      }
      ####End Test arguments
      
      npeaks.global <- 0
      #Counter for % bar
      npspectra <- length(mSet$AnnotateObject$pspectra)
      
      # scaling
      devppm <- ppm / 1000000
      filter <- filter
      
      #generate parameter list
      params <-
        list(
          maxiso = maxiso,
          maxcharge = maxcharge,
          devppm = devppm,
          mzabs = mzabs,
          IM = isotopeMatrix,
          minfrac = minfrac,
          filter = filter
        )
      
      #Check if object have been preprocessed with groupFWHM
      if (npspectra < 1) {
        cat("xsAnnotate contains no pseudospectra. Regroup all peaks into one!\n")
        npspectra <- 1
        
        mSet$AnnotateObject$pspectra[[1]] <-
          seq(1:nrow(mSet$AnnotateObject$groupInfo))
        mSet$AnnotateObject$psSamples  <- 1
      }
      
      #number of peaks in pseudospectra
      ncl <- sum(sapply(mSet$AnnotateObject$pspectra, length))
      
      # get mz,rt and intensity values from peaktable
      if (nrow(mSet$xcmsSet@groups) > 0) {
        ##multiple sample or grouped single sample
        if (is.na(mSet$AnnotateObject$sample[1])) {
          index <- 1:length(mSet$xcmsSet@filepaths)
        } else{
          index <- mSet$AnnotateObject$sample
        }
        
        MessageOutput(
          mes = paste0("Generating peak matrix..."),
          ecol = "\n",
          progress = NULL
        )
        
        mint <-
          groupval(mSet$xcmsSet, value = intval)[, index, drop = FALSE]
        imz <- mSet$AnnotateObject$groupInfo[, "mz", drop = FALSE]
        irt <- mSet$AnnotateObject$groupInfo[, "rt", drop = FALSE]
        
      } else{
        ##one sample case
        MessageOutput(
          mes = paste0("Generating peak matrix..."),
          ecol = "\n",
          progress = NULL
        )
        
        imz  <- mSet$AnnotateObject$groupInfo[, "mz", drop = FALSE]
        irt  <- mSet$AnnotateObject$groupInfo[, "rt", drop = FALSE]
        mint <-
          mSet$AnnotateObject$groupInfo[, intval, drop = FALSE]
        
      }
      
      isotope   <- vector("list", length(imz))
      isomatrix <- matrix(ncol = 5, nrow = 0)
      colnames(isomatrix) <-
        c("mpeak", "isopeak", "iso", "charge", "intrinsic")
      
      MessageOutput(
        mes = paste0("Run isotope peak annotation.."),
        ecol = "\n",
        progress = NULL
      )
      
      lp <- -1
      along = mSet$AnnotateObject$pspectra
      
      pb <-
        progress_bar$new(
          format = "Isotope [:bar] :percent Time left: :eta",
          total = length(along),
          clear = T,
          width = 75
        )
      
      #look for isotopes in every pseudospectra
      for (i in seq(along)) {
        #get peak indizes for i-th pseudospectrum
        ipeak <- mSet$AnnotateObject$pspectra[[i]]
        #Ouput counter
        pb$tick()
        
        #Pseudospectrum has more than one peak
        if (length(ipeak) > 1) {
          #peak mass and intensity for pseudospectrum
          mz  <- imz[ipeak]
          int <- mint[ipeak, , drop = FALSE]
          isomatrix <-
            findIsotopesPspec(isomatrix, mz, ipeak, int, params)
        }
      }
      
      #clean isotopes
      if (is.null(nrow(isomatrix))) {
        isomatrix = matrix(isomatrix,
                           byrow = F,
                           ncol = length(isomatrix))
      }
      
      #check if every isotope has only one annotation
      if (length(idx.duplicated <-
                 which(duplicated(isomatrix[, 2]))) > 0) {
        peak.idx <- unique(isomatrix[idx.duplicated, 2])
        
        for (i in 1:length(peak.idx)) {
          #peak.idx has two or more annotated charge
          #select the charge with the higher cardinality
          peak <- peak.idx[i]
          
          peak.mono.idx <- which(isomatrix[, 2] == peak)
          
          if (length(peak.mono.idx) < 2) {
            #peak has already been deleted
            next
          }
          
          peak.mono <- isomatrix[peak.mono.idx, 1]
          #which charges we have
          charges.list   <- isomatrix[peak.mono.idx, 4]
          tmp <- cbind(peak.mono, charges.list)
          charges.length <- apply(tmp, 1, function(x, isomatrix) {
            length(which(isomatrix[, 1] == x[1] & isomatrix[, 4] == x[2]))
          },
          isomatrix)
          
          idx <- which(charges.length == max(charges.length))
          
          if (length(idx) == 1) {
            #max is unique
            isomatrix <-
              isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] &
                                 isomatrix[, 4] %in% charges.list[-idx]), , drop = FALSE]
          } else{
            #select this one, which lower charge
            idx <- which.min(charges.list[idx])
            isomatrix <-
              isomatrix[-which(isomatrix[, 1] %in% peak.mono[-idx] &
                                 isomatrix[, 4] %in% charges.list[-idx]), , drop = FALSE]
          }
        }
      }
      
      #check if every isotope in one isotope grp, have the same charge
      if (length(idx.duplicated <-
                 which(duplicated(paste(
                   isomatrix[, 1], isomatrix[, 3]
                 )))) > 0) {
        #at least one pair of peakindex and number of isotopic peak is identical
        peak.idx <- unique(isomatrix[idx.duplicated, 1])
        
        for (i in 1:length(peak.idx)) {
          #peak.idx has two or more annotated charge
          #select the charge with the higher cardinality
          peak <- peak.idx[i]
          
          #which charges we have
          charges.list   <-
            unique(isomatrix[which(isomatrix[, 1] == peak), 4])
          
          #how many isotopes have been found, which this charges
          charges.length <-
            sapply(charges.list, function(x, isomatrix, peak) {
              length(which(isomatrix[, 1] == peak &
                             isomatrix[, 4] == x))
            }, isomatrix, peak)
          
          #select the charge which the highest cardinality
          idx <- which(charges.length == max(charges.length))
          
          if (length(idx) == 1) {
            #max is unique
            isomatrix <-
              isomatrix[-which(isomatrix[, 1] == peak &
                                 isomatrix[, 4] %in% charges.list[-idx]), , drop = FALSE]
          } else{
            #select this one, which lower charge
            idx <- which.min(charges.list[idx])
            isomatrix <-
              isomatrix[-which(isomatrix[, 1] == peak &
                                 isomatrix[, 4] %in% charges.list[-idx]), , drop = FALSE]
          }
        }
      }
      
      #Combine isotope cluster, if they overlap
      index2remove <- c()
      
      if (length(idx.duplicated <-
                 which(isomatrix[, 1] %in% isomatrix[, 2])) > 0) {
        for (i in 1:length(idx.duplicated)) {
          index <-  which(isomatrix[, 2] == isomatrix[idx.duplicated[i], 1])
          index2 <-
            sapply(index, function(x, isomatrix)
              which(isomatrix[, 1] == isomatrix[x, 1] &
                      isomatrix[, 3] == 1), isomatrix)
          if (length(index2) == 0) {
            index2remove <- c(index2remove, idx.duplicated[i])
          }
          max.index <- which.max(isomatrix[index, 4])
          
          isomatrix[idx.duplicated[i], 1] <-
            isomatrix[index[max.index], 1]
          
          isomatrix[idx.duplicated[i], 3] <-
            isomatrix[index[max.index], 3] + 1
        }
      }
      
      if (length(index <- which(isomatrix[, "iso"] > maxiso)) > 0) {
        index2remove <- c(index2remove, index)
      }
      
      if (length(index2remove) > 0) {
        isomatrix <- isomatrix[-index2remove, , drop = FALSE]
      }
      
      isomatrix <- isomatrix[order(isomatrix[, 1]), , drop = FALSE]
      #Create isotope matrix within object
      mSet$AnnotateObject$isoID <- matrix(nrow = 0, ncol = 4)
      
      colnames(mSet$AnnotateObject$isoID)  <-
        c("mpeak", "isopeak", "iso", "charge")
      
      #Add isomatrix to object
      mSet$AnnotateObject$isoID <-
        rbind(mSet$AnnotateObject$isoID, isomatrix[, 1:4])
      
      # counter for isotope groups
      globalcnt <- 0
      oldnum    <- 0
      
      if (nrow(isomatrix) > 0) {
        for (i in 1:nrow(isomatrix)) {
          if (!isomatrix[i, 1] == oldnum) {
            globalcnt <- globalcnt + 1
            
            isotope[[isomatrix[i, 1]]] <-
              list(
                y = globalcnt,
                iso = 0,
                charge = isomatrix[i, 4],
                val = isomatrix[i, 5]
              )
            
            oldnum <- isomatrix[i, 1]
          }
          
          isotope[[isomatrix[i, 2]]] <-
            list(
              y = globalcnt,
              iso = isomatrix[i, 3],
              charge = isomatrix[i, 4],
              val = isomatrix[i, 5]
            )
          
        }
      }
      cnt <- nrow(mSet$AnnotateObject$isoID)
      
      MessageOutput(
        mes = paste0("Found isotopes:", cnt),
        ecol = "\n",
        progress = 96
      )
      
      mSet$AnnotateObject$isotopes <- isotope
      
      ## 4. Peak grouping with information -----
      
      cor_eic_th <- annotaParam$corr.eic.th
      
      cor_exp_th <- 0.75
      pval = 0.05
      graphMethod = "hcs"
      calcIso = FALSE
      calcCaS = FALSE
      psg_list = NULL
      xraw = NULL
      intval = "into"
      
      if (!is.numeric(cor_eic_th) ||
          cor_eic_th < 0 || cor_eic_th > 1) {
        stop ("Parameter cor_eic_th must be numeric and between 0 and 1.\n")
      }
      
      npspectra <- length(mSet$AnnotateObject$pspectra)
      
      MessageOutput(
        mes = paste0("Start grouping after correlation..."),
        ecol = "\n",
        progress = NULL
      )
      
      #Data is not preprocessed with groupFWHM
      if (npspectra < 1) {
        MessageOutput(
          mes = paste0(
            "Data was not preprocessed with groupFWHM, creating one pseudospectrum with all peaks."
          ),
          ecol = "\n",
          progress = NULL
        )
        
        #Group all peaks into one group
        npspectra <- 1
        
        mSet$AnnotateObject$pspectra[[1]] <-
          seq(1:nrow(mSet$AnnotateObject$groupInfo))
        
        if (is.na(mSet$AnnotateObject$sample[1])) {
          mSet$AnnotateObject$psSamples <-
            rep(1, nrow(mSet$AnnotateObject$groupInfo))
          ##TODO: Change if sample=NA or sample=number
        } else{
          mSet$AnnotateObject$psSamples <-
            rep(mSet$AnnotateObject$sample,
                nrow(mSet$AnnotateObject$groupInfo))
          
        }
      }
      
      #save number of pspectra before groupCorr
      cnt <- length(mSet$AnnotateObject$pspectra)
      res <- list()
      
      # Check LC information and calcCorr was selected
      
      #Autoselect sample path for EIC correlation
      index <- rep(0, nrow(mSet$AnnotateObject$groupInfo))
      
      for (i in 1:npspectra) {
        index[mSet$AnnotateObject$pspectra[[i]]] <-
          mSet$AnnotateObject$psSamples[[i]]
      }
      
      #Generate EIC data
      tmp <- getAllPeakEICs(mSet, index = index)
      EIC <- tmp$EIC
      scantimes <- tmp$scantimes
      rm(tmp)
      gc()
      
      res[[1]] <- calcCiS(
        mSet,
        EIC = EIC,
        corval = cor_eic_th,
        pval = pval,
        psg_list = psg_list
      )
      
      #Check if we have at least 2 result matrixes
      if (length(res) > 2) {
        #combine the first two to create the result Table
        resMat <- combineCalc(res[[1]], res[[2]], method = "sum")
        for (i in 3:length(res)) {
          resMat <- combineCalc(resMat, res[[i]], method = "sum")
        }
        
      } else if (length(res) == 2) {
        #combine one time
        resMat <- combineCalc(res[[1]], res[[2]], method = "sum")
      } else {
        #Only one matrix
        resMat <- res[[1]]
      }
      
      #Perform graph seperation to seperate co-eluting pseudospectra
      mSet <-
        calcPC.hcs(mSet, ajc = resMat, psg_list = psg_list)
      
      #Create pc groups based on correlation results
      MessageOutput(
        mes = paste(
          "mSet has now",
          length(mSet$AnnotateObject$pspectra),
          "groups, instead of",
          cnt,
          "!"
        ),
        ecol = "\n",
        progress = 97
      )
      
      ## 5. Annotate adducts (and fragments) -----
      mSet$AnnotateObject$ruleset <- NULL
      mSet <-
        findAdducts (
          mSet,
          polarity = annotaParam$polarity,
          mzabs = annotaParam$mz.abs.add,
          maxcharge = annotaParam$max.charge
        )
      
      MessageOutput(mes = NULL,
                    ecol = "\n",
                    progress = 98)
      
      ## 6. Data Organization -----
      camera_output <- getPeaklist(mSet)
      
      sample_names <- mSet$xcmsSet@phenoData[[1]]
      sample_names_ed <-
        gsub(".mzML|.mzData|.mzXML|.cdf|.CDF", "", sample_names)
      
      # Account for multiple groups
      length <- ncol(camera_output)
      end <- length - 3
      camnames <- colnames(camera_output)
      groupNum <-
        nlevels(mSet[["xcmsSet"]]@phenoData[["sample_group"]])
      start <- groupNum + 8
      camnames[start:end] <- sample_names_ed
      colnames(camera_output) <- camnames
      
      endGroup <- 7 + groupNum
      camera_output <- camera_output[,-c(7:endGroup)]
      
      saveRDS(camera_output,
              paste0(fullUserPath, "annotated_peaklist.rds"))
      fast.write(camera_output,
                 paste0(fullUserPath, "annotated_peaklist.csv"))
      
      MessageOutput(
        mes = paste0(
          "Step 6/6: Successfully performed peak annotation! (",
          Sys.time(),
          ") \nGoing to the final step..."
        ),
        ecol = "\n",
        progress = 99
      )
      
      ## 0. Final Output -----
      
      if (running.as.plan) {
        cache.save(mSet, paste0(function.name, "_0"))
        marker_record(paste0(function.name, "_0"))
      }
      
    } else {
      mSet <- cache.read(function.name, "0")
      marker_record("raw_data_samples_c2")
    }
    
    return(mSet)
  }

#' Format Peak List
#' @description This function formats the CAMERA output to a usable format for MetaboAanlyst.
#' @param annotPeaks The object created using the PerformPeakAnnotation.
#' @param annParams The object created using the SetAnnotationParam function,
#' containing user's specified or default parameters for downstream
#' raw MS data pre-processing.
#' @param filtIso Logical, filter out all isotopes except for [M]+ for
#' positive ion mode and [M]- for negative ion mode. By default it is
#' set to true.
#' @param filtAdducts Logical, filter out all adducts except [M+H]+ for
#' positive ion more and [M-H]- for negative ion mode. By default it is set to false.
#' @param missPercent Numeric, specify the threshold to remove features
#' missing in X\% of samples. For instance, 0.5 specifies to remove features
#' that are missing from 50\% of all samples per group. Method is only valid
#' when there are two groups.
#' @author Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' Mai Yamamoto \email{yamamoto.mai@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
FormatPeakList <-
  function(annotPeaks,
           annParams,
           filtIso = TRUE,
           filtAdducts = FALSE,
           missPercent = 0.75) {
    camera_output <- readRDS("annotated_peaklist.rds")
    
    length <- ncol(camera_output)
    end <- length - 3
    
    # Format peaklist for MetaboAnalyst
    camera_ma <- camera_output[,-length]
    
    if (filtAdducts == TRUE) {
      if (annParams$polarity == "positive") {
        if (filtIso == TRUE) {
          camera_isotopes <-
            camera_ma[grepl("\\[M\\]\\+", camera_ma$isotopes), ]
        } else{
          camera_isotopes <- camera_ma[(camera_ma$isotopes != ""), ]
        }
        camera_adducts <-
          camera_ma[grepl("\\[M\\+H\\]\\+", camera_ma$adduct), ]
        camera_feats <-
          camera_ma[(camera_ma$isotopes == "" &
                       camera_ma$adduct == ""), ]
        unique_feats <-
          unique(rbind(camera_isotopes, camera_adducts, camera_feats))
        
      } else{
        # negative polarity
        if (filtIso == TRUE) {
          camera_isotopes <-
            camera_ma[grepl("\\[M\\]-", camera_ma$isotopes), ]
        } else{
          camera_isotopes <- camera_ma[(camera_ma$isotopes != ""), ]
        }
        camera_adducts <-
          camera_ma[grepl("\\[M-H\\]-", camera_ma$adduct), ]
        camera_feats <-
          camera_ma[(camera_ma$isotopes == "" &
                       camera_ma$adduct == ""), ]
        unique_feats <-
          unique(rbind(camera_isotopes, camera_adducts, camera_feats))
      }
    } else{
      if (annParams$polarity == "positive") {
        if (filtIso == TRUE) {
          camera_isotopes <-
            camera_ma[grepl("\\[M\\]\\+", camera_ma$isotopes), ]
        } else{
          camera_isotopes <- camera_ma[(camera_ma$isotopes != ""), ]
        }
        camera_adducts <- camera_ma[(camera_ma$adduct != ""), ]
        camera_feats <- camera_ma[(camera_ma$isotopes == ""), ]
        unique_feats <-
          unique(rbind(camera_isotopes, camera_adducts, camera_feats))
        
      } else{
        # negative polarity
        
        if (filtIso == TRUE) {
          camera_isotopes <-
            camera_ma[grepl("\\[M\\]-", camera_ma$isotopes), ]
        } else{
          camera_isotopes <- camera_ma[(camera_ma$isotopes != ""), ]
        }
        camera_adducts <- camera_ma[(camera_ma$adduct != ""), ]
        camera_feats <- camera_ma[(camera_ma$isotopes == ""), ]
        unique_feats <-
          unique(rbind(camera_isotopes, camera_adducts, camera_feats))
      }
    }
    
    unique_feats <- unique_feats[order(unique_feats[, 1]), ]
    
    # adjust decimal places, feats_info contains all samples
    feats_info <- unique_feats[, 7:end]
    feats_digits <- round(feats_info, 5)
    
    group_info <- annotPeaks$xcmsSet@phenoData[[2]]
    combo_info <- rbind(as.character(group_info), feats_digits)
    
    mzs_rd <-
      paste0(round(unique_feats[, 1], 4), "@", round(unique_feats[, 4], 2))
    mzs <- data.frame(c("Label", mzs_rd), stringsAsFactors = FALSE)
    
    # ensure features are unique
    mzs_unq <- mzs[duplicated(mzs), ]
    
    while (length(mzs_unq) > 0) {
      mzs[duplicated(mzs), ] <-
        sapply(mzs_unq, function(x)
          paste0(x, sample(1:999, 1, replace = FALSE)))
      
      mzs_unq <- mzs[duplicated(mzs), ]
    }
    
    colnames(mzs) <- "Sample"
    ma_feats <- cbind(mzs, combo_info)
    
    # remove features missing in over X% of samples per group
    # only valid for 2 group comparisons!!
    ma_feats_miss <-
      ma_feats[which(rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(unique(group_info[1])))]))
                     |
                       rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(unique(group_info[2])))])) < missPercent),]
    
    fast.write(ma_feats_miss,
               paste0(fullUserPath, "metaboanalyst_input.csv"),
               row.names = FALSE)
    
    # provide index for CAMERA output
    Pklist_inx <- row.names(ma_feats_miss)
    ma_feats_miss_inx <- cbind(ma_feats_miss, Pklist_inx)
    
    fast.write(ma_feats_miss_inx,
               paste0(fullUserPath, "filtered_peaklist.csv"),
               row.names = FALSE)
    
    # generate peak summary results
    peaksum <-
      camera_output[, -c(1:6, (ncol(camera_output) - 2):ncol(camera_output))]
    
    rt_info_min <-
      round(apply(
        peaksum,
        2,
        FUN = function(x) {
          min(camera_output[!is.na(x), 4])
        }
      ), digits = 2)
    
    rt_info_max <-
      round(apply(
        peaksum,
        2,
        FUN = function(x) {
          max(camera_output[!is.na(x), 4])
        }
      ), digits = 2)
    
    rt_range <- paste0(rt_info_min, "~", rt_info_max)
    
    mz_info_min <-
      round(apply(
        peaksum,
        2,
        FUN = function(x) {
          min(camera_output[!is.na(x), 1])
        }
      ), digits = 3)
    
    mz_info_max <-
      round(apply(
        peaksum,
        2,
        FUN = function(x) {
          max(camera_output[!is.na(x), 1])
        }
      ), digits = 3)
    
    mz_range <- paste0(mz_info_min, "~", mz_info_max)
    Group_detail <- ma_feats_miss[1, -1]
    sample_names <- names(Group_detail)
    
    peak_number <-
      apply(
        peaksum,
        2,
        FUN = function(x) {
          length(which(!is.na(x)))
        }
      )
    
    missing_perc <-
      round((nrow(peaksum) - peak_number) / nrow(peaksum) * 100, digits = 2)
    
    datam <- matrix(nrow = length(sample_names), ncol = 6)
    
    datam[, 1] <- sample_names
    datam[, 2] <- as.matrix(unname(Group_detail))
    datam[, 3] <- rt_range
    datam[, 4] <- mz_range
    datam[, 5] <- peak_number
    datam[, 6] <- missing_perc
    
    write.table(
      datam,
      file = "peak_result_summary.txt",
      row.names = F,
      col.names = F,
      quote = F
    )
    
    MessageOutput(
      mes = paste0("Everything has been finished Successfully ! (",
                   Sys.time(),
                   ")"),
      ecol = "\n",
      progress = 100
    )
    
    return(ma_feats_miss)
  }


GeneratePeakList <- function(userPath) {
  setwd(userPath)
  
  ## Claculate the mean internsity of all groups
  sample_data <-
    read.csv("metaboanalyst_input.csv",
             header = T,
             stringsAsFactors = F)
  
  groups <- as.character(as.matrix(sample_data[1, ]))[-1]
  
  sample_data <- sample_data[-1, -1]
  
  
  if (length(unique(groups)) == 1) {
    sample_data_mean  <-
      apply(
        sample_data,
        1,
        FUN = function(x) {
          mean(as.numeric(x), na.rm = T)
        }
      )
    
  } else {
    sample_data1 <- matrix(nrow = nrow(sample_data))
    
    
    for (i in 1:length(unique(groups))) {
      columnnum <- unique(groups)[i] == groups
      sample_data0  <-
        subset.data.frame(sample_data, subset = T, select = columnnum)
      
      sample_data0  <-
        round(apply(
          sample_data0,
          1,
          FUN = function(x) {
            mean(as.numeric(x), na.rm = T)
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
    row.names = F,
    quote = F
  )
  
  return (nrow(ann_data))
}

#' Verify the data is centroid or not
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
CentroidCheck <- function(filename) {
  fileh <- MSnbase:::.openMSfile(filename)
  
  allSpect <- mzR::peaks(fileh, c(1:10))
  
  nValues <- base::lengths(allSpect, use.names = FALSE) / 2
  allSpect <- do.call(rbind, allSpect)
  
  res <- MSnbase:::Spectra1_mz_sorted(
    peaksCount = nValues,
    rt = c(1:10),
    acquisitionNum = c(1:10),
    scanIndex = c(1:10),
    tic = c(1:10),
    mz = allSpect[, 1],
    intensity = allSpect[, 2],
    fromFile = c(1:10),
    centroided = rep(NA, 10),
    smoothed =  rep(NA, 10),
    polarity =  rep(-1, 10),
    nvalues = nValues
  )
  names(res) <-
    paste0("F1.s100",
           c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10"))
  
  
  mzR::close(fileh)
  rm(fileh)
  
  res <- lapply(
    res,
    FUN = function(z, APPLF, ...) {
      pk <- as.data.frame(list(z))
      
      k = 0.025
      qtl = 0.9
      .qtl <- quantile(pk[, 2], qtl)
      x <- pk[pk[, 2] > .qtl, 1]
      quantile(diff(x), 0.25) > k
      
    }
  )
  
  return(sum(unlist(res)) > 8)
}


#' Title
#'
#' @param param0_path 
#' @param users.path 
#'
#'
#' @examples
verifyParam <- function(param0_path, users.path) {
  load(paste0(users.path, "/params.rda"))
  load(param0_path)
  
  verify.vec <- NULL
  
  for (i in 1:length(peakParams)) {
    for (j in 1:length(peakParams0)) {
      if (names(peakParams[i]) == names(peakParams0[j])) {
        verify.vec_tmp <- peakParams[[i]] == peakParams0[[j]]
        verify.vec <- c(verify.vec, verify.vec_tmp)
      }
      
    }
  }
  
  if (all(verify.vec)) {
    return(1)
  } else{
    return(0)
  }
}


#' Title
#'
#' @param mes 
#' @param ecol 
#' @param progress 
#'
#' @return
#'
#' @examples
MessageOutput <- function(mes, ecol, progress) {
  if (!is.null(mes)) {
    if (.on.public.web) {
      # write down message
      write.table(
        mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = ecol
      )
    } else {
      # print message
      message(mes)
    }
  }
  
  # write down progress
  if (.on.public.web & !is.null(progress)) {
    progress <- as.numeric(progress)
    
    write.table(
      progress,
      file = paste0(fullUserPath, "log_progress.txt"),
      row.names = F,
      col.names = F
    )
  }
  
}
