#' @title Perform peak profiling
#' @description This function performs feature extraction of user's raw MS data using
#' the rawData object created using the ImportRawMSData function.
#' @param mSet The object created using the ImportRawMSData function,
#' containing the raw MS data.
#' @param Params The object created using the SetPeakParam function,
#' containing user's specified or default parameters for downstream
#' raw MS data pre-processing.
#' @param plotSettings List, plotting parameters produced by SetPlotParam Function.
#' Defaut is set to true.
#' @param ncore Numeric, used to define the cores' number for Peak Profiling.
#' @param running.controller The resuming pipeline running controller. Optional. Don't need to define by hand.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import MSnbase
#' @import BiocParallel
#' @import ggplot2
#' @examples 
#' ##' Load OptiLCMS package
#' library(OptiLCMS)
#' ##' Get raw spectra files
#' DataFiles <- dir(system.file("mzData", package = "mtbls2"), full.names = TRUE,
#'                  recursive = TRUE)[c(10:12, 14:16)]
#' ##' Create a phenodata data.frame
#' pd <- data.frame(sample_name = sub(basename(DataFiles), pattern = ".mzData",
#'                                    replacement = "", fixed = TRUE),
#'                  sample_group = c(rep("col0", 3), rep("cyp79", 3)),
#'                  stringsAsFactors = FALSE)
#' ##' Import raw spectra
#' mSet <- ImportRawMSData(path = DataFiles, metadata = pd);
#' 
#' ##' Perform spectra profiling
#' mSet <- PerformPeakProfiling(mSet, Params = SetPeakParam(ppm = 5, 
#'                                                          bw = 10, 
#'                                                          mzdiff = 0.001, 
#'                                                          max_peakwidth = 15, 
#'                                                          min_peakwidth = 10), 
#'                                                          ncore = 1, 
#'                              plotSettings = SetPlotParam(Plot = TRUE))
#' 
#' ##' Set peak annotation parameters
#' annParams <- SetAnnotationParam(polarity = 'positive',
#'                                 mz_abs_add = 0.035);
#' 
#' ##' Perform peak annotation
#' mSet <- PerformPeakAnnotation(mSet = mSet,
#'                               annotaParam = annParams,
#'                               ncore =1)
#' 
#' ##' Format the PeakList
#' mSet <- FormatPeakList(mSet = mSet,
#'                        annParams,
#'                        filtIso =FALSE,
#'                        filtAdducts = FALSE,
#'                        missPercent = 1)
#' 
#' ##' Export the annotation result
#' Export.Annotation(mSet, path = tempdir());
#' 
#' ##' Export the Peak Table
#' Export.PeakTable(mSet, path = tempdir());
#' 
#' ##' Export the Peak summary
#' Export.PeakSummary(mSet, path = tempdir())


PerformPeakProfiling <-
  function(mSet,
           Params = NULL,
           plotSettings,
           ncore,
           running.controller = NULL) {
    
    #Build Running plan for data import - Indentify the controller
    function.name <- "peak_profiling";
    
    if (is.null(running.controller)) {
      c1 <- c2 <- c3 <- c4 <- TRUE;
      .running.as.plan <- FALSE;
      
    } else {
      c1 <-
        running.controller@peak_profiling[["c1"]] # used to control peak picking
      c2 <-
        running.controller@peak_profiling[["c2"]] # used to control peak alignment
      c3 <-
        running.controller@peak_profiling[["c3"]] # used to control peak filing
      c4 <-
        running.controller@peak_profiling[["c4"]] # used to control plotting
      
      .running.as.plan <- TRUE;
    }
    
    if(!exists(".SwapEnv")){
      .SwapEnv <<- new.env(parent = .GlobalEnv);
      .SwapEnv$.optimize_switch <- FALSE;
      .SwapEnv$count_current_sample <- 0;
      .SwapEnv$count_total_sample <- 120; # maximum number for on.public.web
      .SwapEnv$envir <- new.env();
    }
    
    ## Clean the existing results
    mSet@peakpicking <- mSet@peakgrouping <- 
      mSet@peakRTcorrection <- mSet@peakfilling <- 
      mSet@peakAnnotation <- 
      list();
    
    .optimize_switch <- .SwapEnv$.optimize_switch <- FALSE;
    
    if(.on.public.web){
      ### Update parameters' style
      write.table(
        unlist(Params),
        file = "param_optimized.txt",
        row.names = TRUE,
        col.names = FALSE,
        quote = FALSE
      )
    }
    
    if(!is.null(Params)){
      mSet@params <- updateRawSpectraParam (Params);
    } else if (is.null(Params) & length(mSet@params) == 0){
      warning("Param is missing. Will use the default generic/general params, 
              but optimization with PerformParamsOptimization is strongly recommanded !");
      mSet@params <- updateRawSpectraParam (SetPeakParam());
    }

    ### Setting the different parallel method for linux or windows
    MessageOutput(mes = NULL,
                  ecol = NULL,
                  progress = 25)
    
    if (.on.public.web) {
      # load_biocparallel()
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
      MessageOutput(mes = "\nStep 3/6: Started peak picking! This step will take some time...",
                    ecol = "\n",
                    progress = NULL);
      
      if(length(mSet@rawOnDisk) == 0){
        stop("No MS found! Please ImportRawMSData first!")
      }
      
      mSet <-
        tryCatch(
          PerformPeakPicking(mSet),
          error = function(e) {
            e
          }
        )
      
      gc()
      
      if (.running.as.plan & class(mSet)[1] != "simpleError") {
        cache.save(mSet, paste0(function.name, "_c1"))
        marker_record(paste0(function.name, "_c1"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c1")
      marker_record("peak_profiling_c1")
    }
    
    .SwapEnv$envir$mSet <- mSet;
    
    if (.on.public.web) {
      if (class(mSet)[1] != "mSet") {
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

    }
    
    MessageOutput(
      mes =  paste0("Step 3/6: Peak picking finished ! (", Sys.time(), ")"),
      ecol = "\n",
      progress = 50
    )
    
    #   --------===========----- II. Peak alignment -----===========------------
    if (c2) {
      MessageOutput(
        mes = paste("\nStep 4/6: Started peak alignment! This step is running..."),
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
      
      if (.running.as.plan & class(mSet)[1] != "simpleError") {
        cache.save(mSet, paste0(function.name, "_c2"))
        marker_record(paste0(function.name, "_c2"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c2")
      marker_record("peak_profiling_c2")
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
    }
    
    .SwapEnv$envir$mSet <- mSet;
    
    MessageOutput(
      mes = paste0("Step 4/6: Peak alignment finished ! (", Sys.time(), ")"),
      ecol = "\n",
      progress = 73
    )
    
    #   --------===========----- III. Peak filling -----===========------------
    if (c3) {
      MessageOutput(
        mes = paste(
          "\nStep 5/6: Started peak filling! This step may take some time..."
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
      
      if (.running.as.plan & class(mSet)[1] != "simpleError") {
        cache.save(mSet, paste0(function.name, "_c3"))
        marker_record(paste0(function.name, "_c3"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c3")
      marker_record("peak_profiling_c3")
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

    if(.on.public.web){
      save(mSet, file = "mSet.rda");
    }
    
    .SwapEnv$envir$mSet <- mSet;

    MessageOutput(
      mes = paste0("Begin to plotting figures..."),
      ecol = "",
      progress = 89
    )
    register(bpstop());
    
    #  ---------------====------IV. Plotting Results --------========-----------
    if (c4) {
      
      sample_idx <- mSet@rawOnDisk@phenoData@data[["sample_group"]]
      
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
      
      if (plotSettings$Plot == TRUE) {
        ### 1. Peak Intensity plotting -----
        PlotSpectraInsensityStistics(
          mSet,
          paste0(plotSettings$name_peak_in, ".", plotSettings$format),
          plotSettings$format,
          plotSettings$dpi,
          8
        )
        
        ### 2. PCA plotting -----
        # if (.on.public.web) {
        #   load_ggplot()
        # }
        
        # PlotSpectraPCA(
        #   mSet,
        #   paste0(plotSettings$name_PCA, ".", plotSettings$format),
        #   plotSettings$format,
        #   plotSettings$dpi,
        #   9
        # )
        
        ### 3. Adjusted RT plotting -----
        PlotSpectraRTadj(
          mSet,
          paste0(plotSettings$name_adj_RT, ".", plotSettings$format),
          plotSettings$format,
          plotSettings$dpi,
          9
        )
        
        ### 4. Chromatogram Generation -----
        PlotSpectraBPIadj(
          mSet,
          paste0(plotSettings$name_adj_BPI, ".", plotSettings$format),
          plotSettings$format,
          plotSettings$dpi,
          9
        )
        
        MessageOutput(
          mes = paste0("Done !"),
          ecol = "\n",
          progress = 89
        )
      }
      
      if(.running.as.plan){
        marker_record(paste0(function.name, "_c4"))
      }
      
    }
    
    MessageOutput(mes = NULL,
                  ecol = NULL,
                  progress = 90);
    
    return(mSet)
  }

#' @title Set annotation parameters
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
#' @param adducts Character, specify the adducts based on your instrument settings.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @seealso \code{\link{ExecutePlan}} and \code{\link{PerformPeakProfiling}} for the whole pipeline.
#' @examples 
#' ##' Load OptiLCMS package
#' library(OptiLCMS)
#' 
#' ##' Set peak annotation parameters
#' annParams <- SetAnnotationParam(polarity = 'positive',
#'                                 mz_abs_add = 0.035);
#' 
#' ##' Please check the example of PerformPeakProfiling 
#' ##' and ExcutePlan for the whole running pipeline.

SetAnnotationParam <-
  function(polarity = "positive",
           perc_fwhm = 0.6,
           mz_abs_iso = 0.005,
           max_charge = 2,
           max_iso = 2,
           corr_eic_th = 0.85,
           mz_abs_add = 0.001,
           adducts = NULL) {
    
    annParams <- list()
    peakParams <- NULL;
    
    if (.on.public.web) {
      if(file.exists("params.rda")){
        load("params.rda")
      } else {
        peakParams <- SetPeakParam();
      }
      
      
      annParams$polarity <- peakParams$polarity
      annParams$perf.whm <- peakParams$perc_fwhm
      annParams$mz.abs.iso <- peakParams$mz_abs_iso
      annParams$max.charge <- peakParams$max_charge
      annParams$max.iso <- peakParams$max_iso
      annParams$corr.eic.th <- peakParams$corr_eic_th
      annParams$mz.abs.add <- peakParams$mz_abs_add;
      annParams$adducts <- unlist(strsplit(peakParams$adducts, "\\|"));
      
    } else {
      annParams$polarity <- polarity
      annParams$perf.whm <- perc_fwhm
      annParams$mz.abs.iso <- mz_abs_iso
      annParams$max.charge <- max_charge
      annParams$max.iso <- max_iso
      annParams$corr.eic.th <- corr_eic_th
      annParams$mz.abs.add <- mz_abs_add
      annParams$adducts <- NULL;
      
    }
    
    return(annParams)
  }

#' @title Perform peak annotation
#' @description This function performs peak annotation on
#' the xset object created using the PerformPeakPicking function.
#' @param mSet mSet object, usually generated by 'PerformPeakProfiling' here.
#' @param annotaParam The object created using the SetAnnotationParam function,
#' containing user's specified or default parameters for downstream
#' raw MS data pre-processing.
#' @param ncore annotation running core. Default is 1. Parallel running will be supported soon.
#' @param running.controller The resuming pipeline running controller. Optional. Don't need to define by hand.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' and Jeff Xia \email{jeff.xia@mcgill.ca}
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
#' @examples 
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' annParams <- SetAnnotationParam(polarity = 'positive',
#'                                 mz_abs_add = 0.035);
#' 
#' ## Perform peak annotation with newly deinfed annParams
#' # mSet <- PerformPeakAnnotation(mSet = mSet,
#' #                               annotaParam = annParams,
#' #                               ncore =1)
#' @seealso \code{\link{ExecutePlan}} and \code{\link{PerformPeakProfiling}} for the whole pipeline.
#' 
PerformPeakAnnotation <-
  function(mSet,
           annotaParam,
           ncore = 1,
           running.controller = NULL) {
    
    MessageOutput(
      mes = paste0("\nStep 6/6: Starting Peak Annotation..."),
      ecol = "\n",
      progress = 91
    )
    
    if(length(mSet@rawOnDisk) == 0){
      if(.on.public.web){
        MessageOutput("ERROR: No MS data imported, please import the MS data with 'ImportRawMSData' first !", NULL, NULL)
      } else {
        stop("No MS data Imported, please import the MS data with 'ImportRawMSData' first !")
      }
    }
    
    if(is.null(mSet@peakfilling$FeatureGroupTable)){
      if(.on.public.web){
        MessageOutput("ERROR: No features found for annotation, please do the \"ImportRawMSData\", \"PerformPeakProfiling\", first before annotation !");
        stop();
      } else {
        stop("No features found for annotation, please do the \"ImportRawMSData\", \"PerformPeakProfiling\", first before annotation !")
      }
    }
    
    # if (.on.public.web) {
    #   dyn.load(.getDynLoadPath());
    #   
    #   load_progress();
    #   load_graph();
    #   load_RBGL();
    # }
    
    if (ncore > 1) {
      MessageOutput("Only single core mode is supported now. Parallel will be supported later !\n")
      ncore <- 1
    }

    function.name <- "peak_annotation"
    
    if (is.null(running.controller)) {
      c1 <- TRUE;
      .running.as.plan <- FALSE;
    } else {
      c1 <- running.controller@peak_annotation[["c1"]];
      .running.as.plan <- TRUE;
    }
    
    if (c1) {
      ## 1. Prepare the Annotation Object-------
      #xs <- mSet$xcmsSet
      
      if (is.null(mSet)) {
        stop("No mSet object was given !")
      } else if (!class(mSet) == "mSet") {
        stop("There is correct mSet object !")
      }
      
      mSet@peakAnnotation$AnnotateObject <- list();
      fgs <- mSet@peakfilling$FeatureGroupTable;
      #xs@groups <- S4Vectors::as.matrix(fgs[, -ncol(fgs)])
      
      if (length(mSet@rawOnDisk@phenoData[["sample_name"]]) > 1 &&
          !nrow(fgs) > 0) {
        stop ('No group information found or contain only one sample.')
      }
      
      mSet@peakAnnotation$peaks <- mSet@peakfilling$msFeatureData$chromPeaks;
      #mSet@peakAnnotation$peaks <- mSet@peakRTcorrection$chromPeaks;
      mSet@peakAnnotation$groups <- S4Vectors::as.matrix(fgs[, -ncol(fgs)])
      rownames(mSet@peakAnnotation$groups) <- NULL;
      mSet@peakAnnotation$groupmat <- mSet@peakAnnotation$groups;
      
      mSet@peakAnnotation$AnnotateObject$sample   <-  as.numeric(NA);
      mSet@peakAnnotation$AnnotateObject$groupInfo <- getPeaks_selection(mSet);
      runParallel <- list();
      runParallel$enable <- 0;
      
      # if (ncore > 1) {
      #   ## If MPI is available ...
      #   rmpi = "Rmpi"
      #   opt.warn <- options("warn")$warn
      #   options("warn" = -1)
      #   if ((Sys.info()["sysname"] != "Windows") &&
      #       require(rmpi, character.only = TRUE) &&
      #       !is.null(ncore)) {
      #     if (is.loaded('mpi_initialize')) {
      #       #test if not already slaves are running!
      #       if (mpi.comm.size() > 0) {
      #         warning(
      #           "There are already intialized mpi slaves on your machine.\nCamera will try to uses them!\n"
      #         )
      #         
      #         runParallel$enable <- 1
      #         runParallel$mode <- rmpi
      #         
      #       } else{
      #         mpi.spawn.Rslaves(ncore = ncore, needlog = FALSE)
      #         if (mpi.comm.size() > 1) {
      #           #Slaves have successfull spawned
      #           runParallel$enable <- 1
      #           runParallel$mode <- rmpi
      #           
      #         } else{
      #           warning(
      #             "Spawning of mpi slaves have failed. CAMERA will run without parallelization.\n"
      #           )
      #         }
      #       }
      #     } else {
      #       #And now??
      #       warning("DLL mpi_initialize is not loaded. Run single core mode!\n")
      #     }
      #     
      #   } else {
      #     #try local sockets using snow package
      #     snow = "snow"
      #     if (try(require(snow, character.only = TRUE, quietly = TRUE))) {
      #       cat("Starting snow cluster with",
      #           ncore,
      #           "local sockets.\n")
      #       snowclust <- makeCluster(ncore, type = "SOCK")
      #       runParallel$enable <- 1
      #       runParallel$mode <- snow
      #       runParallel$cluster <- snowclust
      #     }
      #   }
      #   
      #   options("warn" = opt.warn);
      #   
      #   MessageOutput(mes = "Run cleanParallel after processing to remove the spawned slave processes!",
      #                 ecol = "\n",
      #                 progress = NULL);
      #   
      # }
      
      if (!is.null(annotaParam[["polarity"]])) {
        if (is.na(match.arg(annotaParam[["polarity"]], c("positive", "negative")))) {
          stop("Parameter polarity has to be 'positive' or 'negative' !")
        } else{
          mSet@peakAnnotation$AnnotateObject$polarity <- annotaParam[["polarity"]]
        }
      }
      
      mSet@peakAnnotation$AnnotateObject$runParallel <- runParallel;
      mSet@peakAnnotation$AnnotateObject$annoID <- matrix(ncol = 4, nrow = 0);
      mSet@peakAnnotation$AnnotateObject$annoGrp <- matrix(ncol = 4, nrow = 0);
      mSet@peakAnnotation$AnnotateObject$isoID <- matrix(ncol = 4, nrow = 0);
      
      colnames(mSet@peakAnnotation$AnnotateObject$annoID) <-
        c("id", "grpID", "ruleID", "parentID");
      
      colnames(mSet@peakAnnotation$AnnotateObject$annoGrp) <-
        c("id", "mass", "ips", "psgrp");
      
      colnames(mSet@peakAnnotation$AnnotateObject$isoID)  <-
        c("mpeak", "isopeak", "iso", "charge");
      
      MessageOutput(mes = NULL,
                    ecol = "\n",
                    progress = 92);
      
      ## 2. Group peaks according to their retention time into pseudospectra-groups-----
      
      intval <- "maxo";
      perfwhm <- annotaParam$perf.whm;
      sigma <- 6;
      sample    <- mSet@peakAnnotation$AnnotateObject$sample;
      pspectra  <- list();
      psSamples <- NA;
      
      MessageOutput(
        mes = paste0("Start grouping after retention time."),
        ecol = "\n",
        progress = NULL
      );
      
      if (mSet@peakAnnotation$AnnotateObject$groupInfo[1, "rt"] == -1) {
        # Like FTICR Data
        warning("Warning: no retention times avaiable. Do nothing\n")
        return(invisible(mSet@peakAnnotation$AnnotateObject))
        
      } else{
        if (is.na(sample[1]) || length(mSet@rawOnDisk@processingData@files) > 1) {
          # grouped peaktable within automatic selection or sub selection
          if (is.na(sample[1])) {
            index <- 1:length(mSet@rawOnDisk@processingData@files)
          } else{
            index <- sample
          }
          
          gvals    <- groupval(mSet)[, index, drop = FALSE];
          peakmat  <- mSet@peakAnnotation$peaks;
          groupmat <- mSet@peakAnnotation$groupmat;
          
          #calculate highest peaks
          maxo <-
            as.numeric(apply(gvals, 1, function(x, peakmat) {
              val <- na.omit(peakmat[x, intval])
              if (length(val) == 0) {
                return(NA)
              } else{
                return(max(val))
              }
            }, peakmat));
          
          maxo[which(is.na(maxo))] <- -1;
          maxo <- cbind(1:length(maxo), maxo);
          
          #highest peak index
          int.max <-
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
              MessageOutput(paste0(
                "Warning: Feature ",
                idx,
                " looks odd for at least one peak. Please check afterwards.\n"
              ))
              
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
          peakmat <- mSet@peakAnnotation$peaks
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
              MessageOutput(paste0(
                "Warning: Feature ",
                idx,
                " looks odd for at least one peak. Please check afterwards.\n"
              ))
              
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
        
        mSet@peakAnnotation$AnnotateObject$pspectra  <- pspectra
        mSet@peakAnnotation$AnnotateObject$psSamples <- psSamples
        
        MessageOutput(
          mes = paste(
            "Created",
            length(mSet@peakAnnotation$AnnotateObject$pspectra),
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
      npspectra <- length(mSet@peakAnnotation$AnnotateObject$pspectra)
      
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
        MessageOutput("xsAnnotate contains no pseudospectra. Regroup all peaks into one!\n")
        npspectra <- 1
        
        mSet@peakAnnotation$AnnotateObject$pspectra[[1]] <-
          seq(1:nrow(mSet@peakAnnotation$AnnotateObject$groupInfo))
        mSet@peakAnnotation$AnnotateObject$psSamples  <- 1
      }
      
      #number of peaks in pseudospectra
      ncl <- sum(sapply(mSet@peakAnnotation$AnnotateObject$pspectra, length))
      
      # get mz,rt and intensity values from peaktable
      if (nrow(mSet@peakAnnotation$groups) > 0) {
        ##multiple sample or grouped single sample
        if (is.na(mSet@peakAnnotation$AnnotateObject$sample[1])) {
          index <- 1:length(mSet@rawOnDisk@processingData@files)
        } else{
          index <- mSet@peakAnnotation$AnnotateObject$sample
        }
        
        MessageOutput(
          mes = paste0("Generating peak matrix..."),
          ecol = "\n",
          progress = NULL
        )
        
        mint <-
          groupval(mSet, value = intval)[, index, drop = FALSE]
        imz <- mSet@peakAnnotation$AnnotateObject$groupInfo[, "mz", drop = FALSE]
        irt <- mSet@peakAnnotation$AnnotateObject$groupInfo[, "rt", drop = FALSE]
        
      } else{
        ##one sample case
        MessageOutput(
          mes = paste0("Generating peak matrix..."),
          ecol = "\n",
          progress = NULL
        )
        
        imz  <- mSet@peakAnnotation$AnnotateObject$groupInfo[, "mz", drop = FALSE];
        irt  <- mSet@peakAnnotation$AnnotateObject$groupInfo[, "rt", drop = FALSE];
        mint <-
          mSet@peakAnnotation$AnnotateObject$groupInfo[, intval, drop = FALSE];
      }
      
      isotope   <- vector("list", length(imz));
      isomatrix <- matrix(ncol = 5, nrow = 0);
      colnames(isomatrix) <-
        c("mpeak", "isopeak", "iso", "charge", "intrinsic");
      
      MessageOutput(
        mes = paste0("Run isotope peak annotation.."),
        ecol = "\n",
        progress = NULL
      );
      
      lp <- -1;
      along = mSet@peakAnnotation$AnnotateObject$pspectra;
      
      pb <-
        progress_bar$new(
          format = "Isotope [:bar] :percent Time left: :eta",
          total = length(along),
          clear = TRUE,
          width = 75
        );
      
      #look for isotopes in every pseudospectra
      for (i in seq(along)) {
        #get peak indizes for i-th pseudospectrum
        ipeak <- mSet@peakAnnotation$AnnotateObject$pspectra[[i]]
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
                           byrow = FALSE,
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
      mSet@peakAnnotation$AnnotateObject$isoID <- matrix(nrow = 0, ncol = 4)
      
      colnames(mSet@peakAnnotation$AnnotateObject$isoID)  <-
        c("mpeak", "isopeak", "iso", "charge")
      
      #Add isomatrix to object
      mSet@peakAnnotation$AnnotateObject$isoID <-
        rbind(mSet@peakAnnotation$AnnotateObject$isoID, isomatrix[, 1:4])
      
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
      cnt <- nrow(mSet@peakAnnotation$AnnotateObject$isoID)
      
      MessageOutput(
        mes = paste0("Found isotopes:", cnt),
        ecol = "\n",
        progress = 96
      )
      
      mSet@peakAnnotation$AnnotateObject$isotopes <- isotope
      
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
      
      npspectra <- length(mSet@peakAnnotation$AnnotateObject$pspectra)
      
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
        
        mSet@peakAnnotation$AnnotateObject$pspectra[[1]] <-
          seq(1:nrow(mSet@peakAnnotation$AnnotateObject$groupInfo))
        
        if (is.na(mSet@peakAnnotation$AnnotateObject$sample[1])) {
          mSet@peakAnnotation$AnnotateObject$psSamples <-
            rep(1, nrow(mSet@peakAnnotation$AnnotateObject$groupInfo))
          ##TODO: Change if sample=NA or sample=number
        } else{
          mSet@peakAnnotation$AnnotateObject$psSamples <-
            rep(mSet@peakAnnotation$AnnotateObject$sample,
                nrow(mSet@peakAnnotation$AnnotateObject$groupInfo))
          
        }
      }
      
      #save number of pspectra before groupCorr
      cnt <- length(mSet@peakAnnotation$AnnotateObject$pspectra)
      res <- list()
      
      # Check LC information and calcCorr was selected
      
      #Autoselect sample path for EIC correlation
      index <- rep(0, nrow(mSet@peakAnnotation$AnnotateObject$groupInfo))
      
      for (i in 1:npspectra) {
        index[mSet@peakAnnotation$AnnotateObject$pspectra[[i]]] <-
          mSet@peakAnnotation$AnnotateObject$psSamples[[i]]
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
          length(mSet@peakAnnotation$AnnotateObject$pspectra),
          "groups, instead of",
          cnt,
          "!"
        ),
        ecol = "\n",
        progress = 97
      )
      
      ## 5. Annotate adducts (and fragments) -----
      mSet@peakAnnotation$AnnotateObject$ruleset <- NULL;
      rules <- annotaParam$adducts;
      
      if(rules == "NULL"){
        rules <- NULL;
      }
      
      mSet <-
        findAdducts (
          mSet,
          polarity = annotaParam$polarity,
          mzabs = annotaParam$mz.abs.add,
          maxcharge = annotaParam$max.charge,
          rules = rules
        )
      
      MessageOutput(mes = NULL,
                    ecol = "\n",
                    progress = 98)
      
      ## 6. Data Organization -----
      camera_output <- getPeaklist(mSet)

      sample_names <- pData(mSet@rawOnDisk)[[1]]
      sample_names_ed <-
        gsub(".mzML|.mzData|.mzXML|.cdf|.CDF", "", sample_names)

      # Account for multiple groups
      length <- ncol(camera_output)
      end <- length - 3
      camnames <- colnames(camera_output)
      groupNum <-
        length(unique(pData(mSet@rawOnDisk)[["sample_group"]]));
      start <- groupNum + 8
      camnames[start:end] <- sample_names_ed
      colnames(camera_output) <- camnames
      
      endGroup <- 7 + groupNum
      camera_output <- camera_output[,-c(7:endGroup)]
      fullUserPath <- mSet@WorkingDir;
      
      if(is.null(fullUserPath) | length(fullUserPath) == 0){
        fullUserPath <- getwd();
      }
      
      mSet@peakAnnotation$camera_output <- camera_output;
     
      MessageOutput(
        mes = paste0(
          "Step 6/6: Successfully performed peak annotation! (",
          Sys.time(),
          ") \nGoing to the final step..."
        ),
        ecol = "",
        progress = 99
      )
      
      ## 0. Final Output -----
      
      if (.running.as.plan) {
        cache.save(mSet, paste0(function.name, "_c1"));
        marker_record("peak_annotation_c1");
      }
      
    } else {
      mSet <- cache.read(function.name, "c1");
      marker_record("peak_annotation_c1");
    }
    
    if(.on.public.web){
      save(mSet, file = "mSet.rda");
    }
    MessageOutput("Done!")
    
    return(mSet)
  }

#' @title Format Peak List
#' @description This function formats the CAMERA output to a usable format for OptiLCMS.
#' @param mSet The mSet object generated by the PerformPeakAnnotation function.
#' @param annParams The object created using the SetAnnotationParam function,
#' containing user's specified or default parameters for downstream
#' raw MS data pre-processing.
#' @param filtIso Logical, filter out all isotopes except for `[M]`+ for
#' positive ion mode and `[M]`- for negative ion mode. By default it is
#' set to true.
#' @param filtAdducts Logical, filter out all adducts except `[M+H]`+ for
#' positive ion more and `[M-H]`- for negative ion mode. By default it is set to false.
#' @param missPercent Numeric, specify the threshold to remove features
#' missing in X\% of samples. For instance, 0.5 specifies to remove features
#' that are missing from 50\% of all samples per group. Method is only valid
#' when there are two groups.
#' @seealso \code{\link{ExecutePlan}} and \code{\link{PerformPeakProfiling}} for the whole pipeline.
#' @author Jasmine Chong \email{jasmine.chong@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @examples 
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' annParams <- SetAnnotationParam(polarity = 'positive',
#'                                 mz_abs_add = 0.035);
#' 
#' ## Perform peak annotation with newly deinfed annParams
#' # mSet <- PerformPeakAnnotation(mSet = mSet,
#' #                               annotaParam = annParams,
#' #                               ncore =1)
#' ## Format the PeakList
#' mSet <- FormatPeakList(mSet = mSet,
#'                        annParams,
#'                        filtIso =FALSE,
#'                        filtAdducts = FALSE,
#'                        missPercent = 1)

FormatPeakList <-
  function(mSet,
           annParams,
           filtIso = TRUE,
           filtAdducts = FALSE,
           missPercent = 0.75) {
    
    fullUserPath <- mSet@WorkingDir;
    if(is.null(fullUserPath) | length(fullUserPath) == 0){
      fullUserPath <- getwd();
    }
    
    camera_output <- mSet@peakAnnotation$camera_output;
    
    if(is.null(camera_output)){
      stop("No annotation results found! Please 'PerformPeakAnnotation' first !")
    }
    
    if(length(mSet@rawOnDisk) == 0){
      if(.on.public.web){
        MessageOutput("ERROR: No MS data imported, please import the MS data with 'ImportRawMSData' first !", NULL, NULL)
      } else {
        stop("No MS data Imported, please import the MS data with 'ImportRawMSData' first !")
      }
    }
    
    length <- ncol(camera_output);
    end <- length - 3;
    
    # Format peaklist for MetaboAnalyst
    camera_ma <- camera_output[,-length];
    
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
      
    };
    
    unique_feats <- unique_feats[order(unique_feats[, 1]), ];
    
    # adjust decimal places, feats_info contains all samples
    feats_info <- unique_feats[, 7:end];
    feats_digits <- round(feats_info, 5);
    
    group_info <- mSet@rawOnDisk@phenoData@data[["sample_group"]];
    combo_info <- rbind(as.character(group_info), feats_digits);
    
    mzs_rd <-
      paste0(round(unique_feats[, 1], 4), "__", round(unique_feats[, 4], 2));
    mzs <- data.frame(c("Label", mzs_rd), stringsAsFactors = FALSE);
    
    # ensure features are unique
    mzs_unq <- mzs[duplicated(mzs), ];
    
    while (length(mzs_unq) > 0) {
      mzs[duplicated(mzs), ] <-
        sapply(mzs_unq, function(x)
          paste0(x, sample(1:999, 1, replace = FALSE)))
      
      mzs_unq <- mzs[duplicated(mzs), ]
    };
    
    colnames(mzs) <- "Sample";
    ma_feats <- cbind(mzs, combo_info);
    
    # remove features missing in over X% of samples per group
    # only valid for 2 group comparisons!!
    ma_feats_miss <-
      ma_feats[which(rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(unique(group_info[1])))]))
                     | rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(unique(group_info[2])))])) < missPercent),];
    
    mSet@dataSet <- ma_feats_miss;
    
    # provide index for CAMERA output
    Pklist_inx <- row.names(ma_feats_miss);
    ma_feats_miss_inx <- cbind(ma_feats_miss, Pklist_inx);
    # 
    # fast.write.csv(ma_feats_miss_inx,
    #            paste0(fullUserPath, "/filtered_peaklist.csv"),
    #            row.names = FALSE);
    
    # generate peak summary results
    peaksum <-
      camera_output[, -c(1:6, (ncol(camera_output) - 2):ncol(camera_output))];
    
    rt_info_min <-
      round(apply(
        peaksum,
        2,
        FUN = function(x) {
          min(camera_output[!is.na(x), 4])
        }
      ), digits = 2);
    
    rt_info_max <-
      round(apply(
        peaksum,
        2,
        FUN = function(x) {
          max(camera_output[!is.na(x), 4])
        }
      ), digits = 2);
    
    rt_range <- paste0(rt_info_min, "~", rt_info_max);
    
    mz_info_min <-
      round(apply(
        peaksum,
        2,
        FUN = function(x) {
          min(camera_output[!is.na(x), 1])
        }
      ), digits = 3);
    
    mz_info_max <-
      round(apply(
        peaksum,
        2,
        FUN = function(x) {
          max(camera_output[!is.na(x), 1])
        }
      ), digits = 3);
    
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
      );
    
    missing_perc <-
      round((nrow(peaksum) - peak_number) / nrow(peaksum) * 100, digits = 2);
    datam <- matrix(nrow = length(sample_names), ncol = 6);
    
    datam[, 1] <- sample_names;
    datam[, 2] <- as.matrix(unname(Group_detail));
    datam[, 3] <- rt_range;
    datam[, 4] <- mz_range;
    datam[, 5] <- peak_number;
    datam[, 6] <- missing_perc;

    mSet@peakAnnotation$peak_result_summary <- datam;

    MessageOutput(
      mes = paste0("\nEverything has been finished Successfully ! (",
                   Sys.time(),
                   ")\n"),
      ecol = "",
      progress = 100
    );
    
    ma_feats_miss[is.na(ma_feats_miss)] <- 0;
    ma_feats_miss[ma_feats_miss == ""] <- 0;
    mSet@dataSet <- ma_feats_miss;
    
    PlotSpectraPCA(
      mSet,
      paste0("PCA.png"),
      "png",
      72,
      9
    )
    
    if(.on.public.web){
      save(mSet, file = "mSet.rda");
    }
    
    return(mSet)
  }


#' @title Export.Annotation
#' @description Export.Annotation is used to export the result of annotation
#' @param mSet mSet object, processed by FormatPeakList.
#' @param path character, used to specify the path for result rds and csv file. Default is the working directory.
#' @export
#' @seealso \code{\link{ExecutePlan}} and \code{\link{PerformPeakProfiling}} for the whole pipeline.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jeff Xia \email{jeff.xia@mcgill.ca}
#' @examples
#' data(mSet)
#' Export.Annotation(mSet, path = tempdir())
#' # delete the exported files from the tempdir with unlink
#' unlink(paste0(tempdir(),"/annotated_peaklist.csv"), recursive = TRUE, force = TRUE);
#' unlink(paste0(tempdir(),"/annotated_peaklist.rds"), recursive = TRUE, force = TRUE)

Export.Annotation <- function(mSet = NULL, path = getwd()){
  
  camera_output <- mSet@peakAnnotation$camera_output;
  
  if(is.null(camera_output)){
    stop("No annotation results found! Please 'PerformPeakAnnotation' first !")
  }
  
  saveRDS(camera_output,
          paste0(path, "/annotated_peaklist.rds"))
  fast.write.csv(camera_output,
                 paste0(path, "/annotated_peaklist.csv"))
}

#' @title Export.PeakTable
#' @description Export.PeakTable is used to export the table of peak
#' @param mSet mSet object, processed by FormatPeakList.
#' @param path character, used to specify the path for result rds and csv file. Default is the working directory.#'
#' @export
#' @seealso \code{\link{ExecutePlan}} and \code{\link{PerformPeakProfiling}} for the whole pipeline.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jeff Xia \email{jeff.xia@mcgill.ca}
#' @examples
#' data(mSet);
#' Export.PeakTable(mSet, path = tempdir());
#' # delete the exported files from the tempdir with unlink
#' unlink(paste0(tempdir(),"/metaboanalyst_input.csv"), recursive = TRUE, force = TRUE)

Export.PeakTable <- function(mSet = NULL, path = getwd()){
  
  mSet@dataSet -> ma_feats_miss;
  
  if(is.null(ma_feats_miss)){
    stop("No PeakTable found! Please 'FormatPeakList' first !")
  }
  
  fast.write.csv(ma_feats_miss,
                 paste0(path, "/metaboanalyst_input.csv"),
                 row.names = FALSE);
}

#' @title Export.PeakSummary
#' @description Export.PeakSummary is used to export the result of peak' summary
#' @param mSet mSet object, processed by FormatPeakList.
#' @param path character, used to specify the path for result rds and csv file. Default is the working directory.#'
#' @export
#' @seealso \code{\link{ExecutePlan}} and \code{\link{PerformPeakProfiling}} for the whole pipeline.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jeff Xia \email{jeff.xia@mcgill.ca}
#' @examples
#' data(mSet);
#' Export.PeakSummary(mSet, path = tempdir());
#' # delete the exported files from the tempdir with unlink
#' unlink(paste0(tempdir(),"/peak_result_summary.txt"), recursive = TRUE, force = TRUE)

Export.PeakSummary <- function(mSet = NULL, path = getwd()){
  
  mSet@peakAnnotation$peak_result_summary -> datam;
  
  if(is.null(datam)){
    stop("No PeakSummary found! Please 'FormatPeakList' first !")
  }
  
  write.table(
    datam,
    file = paste0(path, "/peak_result_summary.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  );
}



