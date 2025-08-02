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
#' @return will return a complete mSet object with the whole processes finished
#' @examples 
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
#' mSet <- PerformPeakProfiling(mSet, Params = SetPeakParam(ppm = 15, 
#'                                                          bw = 10, 
#'                                                          mzdiff = 0.001, 
#'                                                          max_peakwidth = 15, 
#'                                                          min_peakwidth = 10), 
#'                                                          ncore = 2, 
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
    
    if(.on.public.web()){
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
    
    if (.on.public.web()) {
      # load_biocparallel()
      total_threads <- 4
      
    } else if (missing(ncore)) {
      total_threads <- detectCores() * 2 / 3
    } else {
      total_threads <- ncore
    }
    
    #if (.Platform$OS.type == "unix") {
    #  register(bpstart(MulticoreParam(ceiling(total_threads))))
    #} else if (.Platform$OS.type == "windows") {
    #  register(bpstart(SnowParam(ceiling(total_threads))))
    #}
    
    if (.Platform$OS.type=="unix"){
      bp <- MulticoreParam(ceiling(total_threads));
    } else {
      bp <- SnowParam(ceiling(total_threads));
    };
    
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
          PerformPeakPicking(mSet, BPPARAM = bp),
          error = function(e) {
            e
          }
        )
      
      gc()
      
      if (.running.as.plan & !is(mSet,"simpleError")) {
        cache.save(mSet, paste0(function.name, "_c1"))
        marker_record(paste0(function.name, "_c1"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c1")
      marker_record("peak_profiling_c1")
    }
    
    .SwapEnv$envir$mSet <- mSet;
    
    if (.on.public.web()) {
      if (!is(mSet, "mSet")) {
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
        stop("EXCEPTION POINT CODE: Peak Picking. Please adjust your parameters and try again.")
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
      
      if (.running.as.plan & !is(mSet, "simpleError")) {
        cache.save(mSet, paste0(function.name, "_c2"))
        marker_record(paste0(function.name, "_c2"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c2")
      marker_record("peak_profiling_c2")
    }
    
    if (.on.public.web()) {

      if (is(mSet, "simpleError")) {
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
        stop("EXCEPTION POINT CODE: Peak Alignment. Please adjust your parameters and try again")
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
          PerformPeakFiling (mSet, BPPARAM = bp),
          error = function(e) {
            e
          }
        )
      
      gc();
      
      if (.running.as.plan & !is(mSet,"simpleError")) {
        cache.save(mSet, paste0(function.name, "_c3"))
        marker_record(paste0(function.name, "_c3"))
      }
      
    } else {
      mSet <- cache.read(function.name, "c3")
      marker_record("peak_profiling_c3")
    }
    
    if (is(mSet,"simpleError")) {
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

    if(.on.public.web()){
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
        # if (.on.public.web()) {
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
#' @return will return a annotation parameter set for following annotation steps
#' @seealso \code{\link{ExecutePlan}} and \code{\link{PerformPeakProfiling}} for the whole pipeline.
#' @examples 
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
    
    if (.on.public.web()) {
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
      annParams$adducts <- 'NULL';
      
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
#' @return will return an mSet object wirh annotation finished
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
      if(.on.public.web()){
        MessageOutput("ERROR: No MS data imported, please import the MS data with 'ImportRawMSData' first !", NULL, NULL)
      } else {
        stop("No MS data Imported, please import the MS data with 'ImportRawMSData' first !")
      }
    }
    
    if(is.null(mSet@peakfilling$FeatureGroupTable)){
      if(.on.public.web()){
        MessageOutput("ERROR: No features found for annotation, please do the \"ImportRawMSData\", \"PerformPeakProfiling\", first before annotation !");
        stop();
      } else {
        stop("No features found for annotation, please do the \"ImportRawMSData\", \"PerformPeakProfiling\", first before annotation !")
      }
    }
    
    # if (.on.public.web()) {
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
      } else if (!is(mSet,"mSet")) {
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
        warning("no retention times avaiable. Do nothing\n")
        return(invisible(mSet@peakAnnotation$AnnotateObject))
        
      } else{
        if (is.na(sample[1]) || length(mSet@rawOnDisk@processingData@files) > 1) {
          # grouped peaktable within automatic selection or sub selection
          if (is.na(sample[1])) {
            if(mSet@params[["BlankSub"]]){
              index <- seq_along(which(mSet@rawOnDisk@phenoData@data[["sample_group"]] != "BLANK"))
            } else {
              index <- seq_along(mSet@rawOnDisk@processingData@files)
            }
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
          maxo <- cbind(seq_along(maxo), maxo);
          
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
          maxo    <- cbind(seq_along(maxo), maxo)
          
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
          seq_len(nrow(mSet@peakAnnotation$AnnotateObject$groupInfo))
        mSet@peakAnnotation$AnnotateObject$psSamples  <- 1
      }
      
      #number of peaks in pseudospectra
      ncl <- sum(sapply(mSet@peakAnnotation$AnnotateObject$pspectra, length))
      
      # get mz,rt and intensity values from peaktable
      if (nrow(mSet@peakAnnotation$groups) > 0) {
        ##multiple sample or grouped single sample
        if (is.na(mSet@peakAnnotation$AnnotateObject$sample[1])) {
          if(mSet@params[["BlankSub"]]){
            index <- seq_along(which(mSet@rawOnDisk@phenoData@data[["sample_group"]] != "BLANK"))
          } else {
            index <- seq_along(mSet@rawOnDisk@processingData@files)
          }
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
        
        for (i in seq_along(peak.idx)) {
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
        
        for (i in seq_along(peak.idx)) {
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
        for (i in seq_along(idx.duplicated)) {
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
        rbind(mSet@peakAnnotation$AnnotateObject$isoID, isomatrix[, seq_len(4)])
      
      # counter for isotope groups
      globalcnt <- 0
      oldnum    <- 0
      
      if (nrow(isomatrix) > 0) {
        for (i in seq_len(nrow(isomatrix))) {
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
        progress = 95
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
          seq_len(nrow(mSet@peakAnnotation$AnnotateObject$groupInfo))
        
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
      
      for (i in seq_len(npspectra)) {
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
        progress = 96
      )
      
      ## 5. Annotate adducts (and fragments) -----
      mSet@peakAnnotation$AnnotateObject$ruleset <- NULL;
      rules <- annotaParam$adducts;
      
      if(rules[1] == "NULL"){
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
                    progress = 97)
      
      ## 6. Mass Matching
      if(.on.public.web()){
        mSet <- PerformMassMatching(mSet);
      }
      
      ## 7. Data Organization -----
      camera_output <- getPeaklist(mSet)

      sample_names <- pData(mSet@rawOnDisk)[[1]];
      if(mSet@params$BlankSub){
        sample_names <- sample_names[mSet@rawOnDisk@phenoData@data[["sample_group"]] != "BLANK"]
      }
      sample_names_ed <-
        gsub(".mzML|.mzData|.mzXML|.cdf|.CDF", "", sample_names)


      # Account for multiple groups
      length <- ncol(camera_output)
      end <- length - 3
      camnames <- colnames(camera_output)
      groupNum <-
        length(unique(pData(mSet@rawOnDisk)[["sample_group"]]));
      if(mSet@params$BlankSub & any(mSet@rawOnDisk@phenoData@data[["sample_group"]] == "BLANK")){
        groupNum <- groupNum -1;
      }
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
        progress = 98
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
    
    if(.on.public.web()){
      save(mSet, file = "mSet.rda");
    }
    
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
#' missing in X\\% of samples. For instance, 0.5 specifies to remove features
#' that are missing from 50\\% of all samples per group. Method is only valid
#' when there are two groups.
#' @return will return a mSet object with all result table formatted
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
      if(.on.public.web()){
        MessageOutput("ERROR: No MS data imported, please import the MS data with 'ImportRawMSData' first !", NULL, NULL)
      } else {
        stop("No MS data Imported, please import the MS data with 'ImportRawMSData' first !")
      }
    }
    
    length <- ncol(camera_output);
    end <- length - 3;
    
    # Format peaklist for MetaboAnalyst
    camera_ma <- camera_output[,-length];
    camera_ma <- cbind(camera_ma, Formula = mSet@peakAnnotation[["massMatching"]][["Formula"]])
    
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
        
      } else {
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
    if(mSet@params$BlankSub & any(group_info == "BLANK")){
      group_info <- group_info[group_info != "BLANK"]
    }
    combo_info <- rbind(as.character(group_info), feats_digits);
    
    mzs_rd <-
      paste0(round(unique_feats[, 1], 4), "__", round(unique_feats[, 4], 2));
    mzs <- data.frame(c("Label", mzs_rd), stringsAsFactors = FALSE);
    
    # ensure features are unique
    mzs_unq <- mzs[duplicated(mzs), ];
    
    while (length(mzs_unq) > 0) {
      mzs[duplicated(mzs), ] <-
        sapply(mzs_unq, function(x)
          paste0(x, sample(seq_len(999), 1, replace = FALSE)))
      
      mzs_unq <- mzs[duplicated(mzs), ]
    };
    
    colnames(mzs) <- "Sample";
    ma_feats <- cbind(mzs, combo_info);
    
    # remove features missing in over X% of samples per group
    # only valid for 2 group comparisons!!
    # ma_feats_miss <-
    #   ma_feats[which(rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(unique(group_info[1])))]))
    ##                  | rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(unique(group_info[2])))])) < missPercent),];
    # 
    
    grps_tb <- table(group_info)
    grps_info <- names(grps_tb[grps_tb > 2])
    if(length(grps_info) > 1){
      ma_feats_miss <-
        ma_feats[which(rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(grps_info[1]))]))
                       | rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(grps_info[2]))])) <= missPercent),];
    } else {
      ma_feats_miss <-
        ma_feats[which(rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(grps_info[1]))])) <= missPercent),];
    }
    
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
      camera_output[, -c(seq_len(6), (ncol(camera_output) - 2):ncol(camera_output))];
    
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
    if(.on.public.web()) {
      mSet@peakAnnotation[["massMatching"]][["Compound"]][is.na(mSet@peakAnnotation[["massMatching"]][["Compound"]])] <- 
        mSet@peakAnnotation[["massMatching"]][["Formula"]][is.na(mSet@peakAnnotation[["massMatching"]][["Formula"]])] <- 
        "";
      
      mSet@peakAnnotation$camera_output <- 
        cbind(camera_output, mSet@peakAnnotation[["massMatching"]]);
      
    } else {
      mSet@peakAnnotation[["massMatching"]][["Compound"]][is.na(mSet@peakAnnotation[["massMatching"]][["Compound"]])] <- 
        mSet@peakAnnotation[["massMatching"]][["Formula"]][is.na(mSet@peakAnnotation[["massMatching"]][["Formula"]])] <- 
        '';
    }

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
    
    if(.on.public.web()){
      # remove isotopes and adducts directly
      
      if (annParams$polarity == "positive") {
        isotopes_idx <- grepl("\\[M\\]\\+", camera_ma$isotopes)
        camera_isotopes <-
          camera_ma[isotopes_idx, ]
        adducts_idx <- grepl("\\[M\\+H\\]\\+", camera_ma$adduct)
        camera_adducts <-
          camera_ma[adducts_idx, ]
        
        idx1 <- unique(c(which(isotopes_idx), which(adducts_idx)))
        pc_grp_idx1 <- camera_output$pcgroup[idx1]
        
        parent_ion_idx <- which(camera_ma$isotopes == "" &
                             camera_ma$adduct == "")
        pc_grp_parent_ion_idx <- camera_output$pcgroup[parent_ion_idx]
        
        parent_ion_idx <- parent_ion_idx[!(pc_grp_parent_ion_idx %in% pc_grp_idx1)]
        
        camera_feats <-
          camera_ma[parent_ion_idx, ]
        camera_feats_out <-
          camera_output[parent_ion_idx, ]
        
        unique_pc_idx <- unique(camera_feats_out$pcgroup)
        unique_pc_inx <- vapply(unique_pc_idx, function(x){
          uid <- which(camera_feats_out$pcgroup == x)
          if(length(uid)==1){
            return(uid)
          } else {
            uuid <- findMostAbundantFeature(camera_feats_out[uid, ])
            return(uid[uuid])
          }
        }, FUN.VALUE = integer(1L))
        camera_feats_done <- camera_feats[unique_pc_inx, ]
        
        unique_feats <-
          unique(rbind(camera_isotopes, camera_adducts, camera_feats_done))
        
      } else{
        # negative polarity
        isotopes_idx <- grepl("\\[M\\]\\-", camera_ma$isotopes)
        camera_isotopes <-
          camera_ma[isotopes_idx, ]
        adducts_idx <- grepl("\\[M\\-H\\]\\-", camera_ma$adduct)
        camera_adducts <-
          camera_ma[adducts_idx, ]
        
        idx1 <- unique(c(which(isotopes_idx), which(adducts_idx)))
        pc_grp_idx1 <- camera_ma$pcgroup[idx1]
        
        parent_ion_idx <- which(camera_ma$isotopes == "" &
                                  camera_ma$adduct == "")
        pc_grp_parent_ion_idx <- camera_ma$pcgroup[parent_ion_idx]
        
        parent_ion_idx <- parent_ion_idx[!(pc_grp_parent_ion_idx %in% pc_grp_idx1)]
        
        camera_feats <-
          camera_ma[parent_ion_idx, ]
        camera_feats_out <-
          camera_output[parent_ion_idx, ]
        
        unique_pc_idx <- unique(camera_feats_out$pcgroup)
        unique_pc_inx <- vapply(unique_pc_idx, function(x){
          uid <- which(camera_feats_out$pcgroup == x)
          if(length(uid)==1){
            return(uid)
          } else {
            uuid <- findMostAbundantFeature(camera_feats_out[uid, ])
            return(uid[uuid])
          }
        }, FUN.VALUE = integer(1L))
        camera_feats_done <- camera_feats[unique_pc_inx, ]
        
        unique_feats <-
          unique(rbind(camera_isotopes, camera_adducts, camera_feats_done))
        
      }
      
      unique_feats <- unique_feats[order(unique_feats[, 1]), ];
      
      # adjust decimal places, feats_info contains all samples
      feats_info <- unique_feats[, 7:end];
      feats_digits <- round(feats_info, 5);
      
      group_info <- mSet@rawOnDisk@phenoData@data[["sample_group"]];
      if(mSet@params$BlankSub & any(group_info == "BLANK")){
        group_info <- group_info[group_info != "BLANK"]
      }
      combo_info <- rbind(as.character(group_info), feats_digits);
      
      mzs_rd <-
        paste0(round(unique_feats[, 1], 4), "__", round(unique_feats[, 4], 2));
      mzs <- data.frame(c("Label", mzs_rd), stringsAsFactors = FALSE);
      
      # ensure features are unique
      mzs_unq <- mzs[duplicated(mzs), ];
      
      while (length(mzs_unq) > 0) {
        mzs[duplicated(mzs), ] <-
          sapply(mzs_unq, function(x)
            paste0(x, sample(seq_len(999), 1, replace = FALSE)))
        
        mzs_unq <- mzs[duplicated(mzs), ]
      };
      
      colnames(mzs) <- "Sample";
      ma_feats <- cbind(mzs, combo_info);
      
      # remove features missing in over X% of samples per group
      # only valid for 2 group comparisons!!
      # ma_feats_miss <-
      #   ma_feats[which(rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(unique(group_info[1])))]))
      ##                  | rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(unique(group_info[2])))])) < missPercent),];
      # 

      grps_tb <- table(group_info)
      grps_info <- names(grps_tb[grps_tb > 2])
      if(length(grps_info) > 1){
        ma_feats_miss <-
          ma_feats[which(rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(grps_info[1]))]))
                         | rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(grps_info[2]))])) <= missPercent),];
      } else {
        ma_feats_miss <-
          ma_feats[which(rowMeans(is.na(ma_feats[, (ma_feats[1, ] == as.character(grps_info[1]))])) <= missPercent),];
      }
      
      
      # correction based on RSD of QC, if any
      if("QC" %in% group_info){
        qc_idx <-which(group_info == "QC")
        ma_feats_miss0 <- ma_feats_miss[-1,-1]
        res_idx <- apply(ma_feats_miss0, 1, function(x){
          qc_vec <- as.numeric(x[qc_idx])
          rsd <- sd(qc_vec, na.rm = T)/mean(qc_vec, na.rm = T)
          if(is.na(rsd)){
            return(FALSE)
          }
          return(rsd<=1)
        })
        ma_feats_miss <- ma_feats_miss[c(1, which(res_idx)+1),]
        unique_feats <- unique_feats[res_idx,]
      }
      
      # adding empirical compound information if any
      empirical_cmpds <- apply(unique_feats, 1, function(y){
        fres <- y[length(y)]
        strsplit(fres, "; ")[[1]][1]
      })
      
      new_labels <- vapply(1:length(empirical_cmpds), function(x){
        if(is.na(empirical_cmpds[x])){
          return(ma_feats_miss[x+1, 1])
        } else {
          return(paste0(empirical_cmpds[x], "_", ma_feats_miss[x+1, 1]))
        }
      }, character(1L))
      ma_feats_miss[2:nrow(ma_feats_miss), 1] <- new_labels
      
      ms1_res_dt <- list(ma_feats_miss, unique_feats)
      write.csv(ma_feats_miss, file = "metaboanalyst_input_clean.csv", row.names = F, quote = F)
      qs::qsave(ms1_res_dt, file = "metaboanalyst_input_clean_MS1.qs")
      # save the rda
      save(mSet, file = "mSet.rda");
    }
    
    return(mSet)
  }

findMostAbundantFeature <- function(ft){
  ft<- ft[-c(1:6, ncol(ft):(ncol(ft)-5))]
  res <- apply(ft, 1, function(x){mean(x, na.rm = T)})
  return(which.max(res))
}

#' @title Export.Annotation
#' @description Export.Annotation is used to export the result of annotation
#' @param mSet mSet object, processed by FormatPeakList.
#' @param path character, used to specify the path for result rds and csv file. Default is the working directory.
#' @return will save annotated_peaklist.rds and annotated_peaklist.csv into working path
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
#' @return will save metaboanalyst_input.csv into working path
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
  MessageOutput(
    mes = paste0("Done!"),
    ecol = "\n",
    progress = 99
  );
  fast.write.csv(ma_feats_miss,
                 paste0(path, "/metaboanalyst_input.csv"),
                 row.names = FALSE);
  MessageOutput(
    mes = paste0("\nEverything for MS1 spectra has been finished successfully ! (",
                 Sys.time(),
                 ")\n"),
    ecol = "",
    progress = 100
  );
}

#' @title Export.PeakSummary
#' @description Export.PeakSummary is used to export the result of peak' summary
#' @param mSet mSet object, processed by FormatPeakList.
#' @param path character, used to specify the path for result rds and csv file. Default is the working directory.#'
#' @export
#' @return will save peak_result_summary.txt into working path
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


PerformMassMatching <- function(mSet = NA) {
  if(!.on.public.web()) {
    message("PerformMassMatching is not supported for local version now!")
    return(0)
  }
  
  MessageOutput(mes = "Performing mass matching to HMDB database...", 
                ecol = "")
  libPath <- system.file('hmdb/hmdb_all.rds', package = "OptiLCMS",
                         lib.loc=.libPaths());
  hmdb <- readRDS(file = libPath);
  hmdb <- hmdb[!is.na(hmdb$mono_weight),]
  ### running matching
  res <- sapply(seq(nrow(mSet@peakAnnotation[["AnnotateObject"]][["groupInfo"]])),
         FUN = function(x) {
           thismz0 <- mSet@peakAnnotation[["AnnotateObject"]][["groupInfo"]][x,1];
           thismz1 <- thismz2 <- thismz3 <- numeric();
           isoAnnotated <- !is.null(mSet@peakAnnotation[["AnnotateObject"]][["isotopes"]][[x]])
           addAnnotated <- !is.null(mSet@peakAnnotation[["AnnotateObject"]][["derivativeIons"]][[x]])
           
           if(isoAnnotated) {
             chargeVal <- mSet@peakAnnotation[["AnnotateObject"]][["isotopes"]][[x]]$charge;
             isoVal <- mSet@peakAnnotation[["AnnotateObject"]][["isotopes"]][[x]]$iso;
             thismz1 <- thismz0*chargeVal - isoVal*1.007825;
           }
           
           if(addAnnotated) {
             addlists <- mSet@peakAnnotation[["AnnotateObject"]][["derivativeIons"]][[x]];
             for(i in seq(addlists)){
               thismz2 <- c(thismz2, addlists[[i]]$mass)
             }
           }
           
           if(!isoAnnotated & !addAnnotated){
             if(mSet@peakAnnotation[["AnnotateObject"]][["polarity"]] == "negative"){
               thismz3 <- thismz0 + 1.007825;
             } else {
               thismz3 <- thismz0 - 1.007825;
             }
           }
           
           thismz <- c(thismz1, thismz2, thismz3)
           return(thismz)
         })
  ppmValue <- mSet@params$ppm
  matchingRes <- lapply(res, FUN = function(x){
    mRes <- numeric()
    for(i in x){
      iup <- i + i*ppmValue*1e-6;
      idw <- i - i*ppmValue*1e-6;
      mRes <- c(mRes, head(which(hmdb$mono_weight > idw & hmdb$mono_weight < iup),5))
    }
    
    if(length(mRes) == 0){
      iup <- i + i*ppmValue*2e-6;
      idw <- i - i*ppmValue*2e-6;
      mRes <- c(mRes, head(which(hmdb$mono_weight > idw & hmdb$mono_weight < iup), 5))
    }
    formulaList <- hmdb$formula[mRes];
    compoundList <- hmdb$compound_name[mRes];
    HMDBIDList <- hmdb$HMDB_ID[mRes];
    corRes <- list()
    for(f in unique(formulaList)){
      corRes[[f]]$cmpd <- compoundList[formulaList == f]
      corRes[[f]]$hmdb <- HMDBIDList[formulaList == f]
    }
    
    MixRes <- list(compoundList = compoundList, formulaList = unique(formulaList), HMDBIDList = HMDBIDList)
    return(list(corRes = corRes, MixRes = MixRes))
  })
  
  MassMatchingTable <- data.frame()
  for(i in seq(length(matchingRes))){
    if(length(matchingRes[[i]]$MixRes$compoundList) == 0){
      tmpTable <- data.frame(Compound = NA, Formula = NA, HMDBID = NA)
    } else {
      tmpTable <- data.frame(Compound = gsub("\t|\\s+", " ",paste0(matchingRes[[i]]$MixRes$compoundList, sep = "; ", collapse = "")), 
                             Formula = gsub("\t|\\s+", " ",paste0(matchingRes[[i]]$MixRes$formulaList, sep = "; ", collapse = "")),
                             HMDBID = gsub("\t|\\s+", " ",paste0(matchingRes[[i]]$MixRes$HMDBIDList, sep = "; ", collapse = "")))
    }
    
    MassMatchingTable <- rbind(MassMatchingTable, tmpTable)
  }
  
  matchCoreRes <- lapply(matchingRes, FUN = function(x){
    x$corRes
  })
  mSet@peakAnnotation$massMatching <- MassMatchingTable;
  mSet@peakAnnotation$Formula2Cmpd <- matchCoreRes;
  ###
  
  MessageOutput(mes = "Done!", 
                ecol = "\n")
  return(mSet)
}

PerformFormulaGeneration <- function(mSet = NULL, ms_instr = NULL, 
                                     halogen = FALSE, ppm = TRUE, 
                                     ms1_tol = 5, ms2_tol = 10, charge = 1,
                                     parallel = TRUE, ncore = 4,  
                                     timeout_secs = 600, batch_size = 1000, 
                                     c_range = 100, h_range = 150, ms2_denoise = TRUE){
  
  # Check if ms_instr is "qtof", "orbitrap" or "fticr"
  if (!ms_instr %in% c("qtof", "orbitrap", "fticr")) {
    stop("Error: ms_instr must be one of 'qtof', 'orbitrap', or 'fticr'")
  }
  
  # Check if there is Concensus spectra
  if (length(mSet@MSnResults[["Concensus_spec"]][[1]]) == 0) {
    stop("Error: There is no Concensus spectra in the mSet object.")
  }
  
  # Check the charge
  if (!is.numeric(charge)) {
    stop("'charge' must be a numeric value or numeric vector.")
  }
  
  len_charge <- length(charge)
  
  if (mSet@MSnData[["acquisitionMode"]] == "DDA") {
    len_targeted_peaks <- nrow(mSet@MSnData[["peak_mtx"]])
  }
  
  if (mSet@MSnData[["acquisitionMode"]] == "DIA") {
    len_targeted_peaks <- length(mSet@MSnData[["peak_mtx"]])
  }
  
  if (len_charge != 1 && len_charge != len_targeted_peaks) {
    stop(sprintf("'charge' should either be a single integer or a numeric vector of length %d.", len_targeted_peaks))
  }
  
  ##importing data from mSet (with columns ‘mz’, ‘intensity’) containing an MS/MS spectrum 
  ##Import the required Python modules
  msbuddy <- reticulate::import("msbuddy")
  base <- msbuddy$base
  np <- reticulate::import("numpy")
  
  # Instantiate a Msbuddy object
  config <- msbuddy$MsbuddyConfig(ms_instr= ms_instr, 
                                  halogen=halogen, ppm = ppm, 
                                  ms1_tol = ms1_tol, ms2_tol = ms2_tol, 
                                  parallel = parallel, n_cpu = ncore, 
                                  timeout_secs = timeout_secs, batch_size = batch_size, 
                                  c_range = c(0, c_range), h_range = c(0, h_range))
  engine <- msbuddy$Msbuddy(config)
  
  # generate the features_list
  features_list <- list()
  
  if (mSet@MSnData[["acquisitionMode"]] == "DDA") {
    # import all the concensus spectra
    if (len_charge == 1) {
      for(i in mSet@MSnResults[["Concensus_spec"]][[1]]) {
        
        # get MS2 spectrum from the current mSet object
        position <- which(mSet@MSnResults[["Concensus_spec"]][[1]] == i)
        ms2_df <- mSet@MSnResults[["Concensus_spec"]][[2]][[position]][[1]]
        colnames(ms2_df) <- c("mz", "intensity")
        
        # Create a Spectrum object for MS2
        ms2_spec <- base$Spectrum(mz_array = np$array(ms2_df['mz']),
                                  int_array = np$array(ms2_df['intensity']))
        
        # Create a MetaFeature object (assuming the other parameters remain constant for all files)
        metafeature <- base$MetaFeature(identifier = i+1,
                                        mz = mean(mSet@MSnData[["peak_mtx"]][i+1,1], mSet@MSnData[["peak_mtx"]][i+1,2]),
                                        rt = mean(mSet@MSnData[["peak_mtx"]][i+1,3], mSet@MSnData[["peak_mtx"]][i+1,4]),
                                        charge = charge,
                                        ms2 = ms2_spec)
        
        # Append metafeature to the list
        features_list <- c(features_list, list(metafeature))
      }
    } else{
      for(i in mSet@MSnResults[["Concensus_spec"]][[1]]) {
        
        # get MS2 spectrum from the current mSet object
        position <- which(mSet@MSnResults[["Concensus_spec"]][[1]] == i)
        ms2_df <- mSet@MSnResults[["Concensus_spec"]][[2]][[position]][[1]]
        colnames(ms2_df) <- c("mz", "intensity")
        
        # Create a Spectrum object for MS2
        ms2_spec <- base$Spectrum(mz_array = np$array(ms2_df['mz']),
                                  int_array = np$array(ms2_df['intensity']))
        
        # Create a MetaFeature object (assuming the other parameters remain constant for all files)
        metafeature <- base$MetaFeature(identifier = i+1,
                                        mz = mean(mSet@MSnData[["peak_mtx"]][i+1,1], mSet@MSnData[["peak_mtx"]][i+1,2]),
                                        rt = mean(mSet@MSnData[["peak_mtx"]][i+1,3], mSet@MSnData[["peak_mtx"]][i+1,4]),
                                        charge = charge[i+1],
                                        ms2 = ms2_spec)
        
        # Append metafeature to the list
        features_list <- c(features_list, list(metafeature))
      }
    }
  }
  
  
  if (mSet@MSnData[["acquisitionMode"]] == "DIA") {
    # import all the concensus spectra
    if (len_charge == 1) {
      for(i in mSet@MSnResults[["Concensus_spec"]][[1]]) {
        
        # get MS2 spectrum from the current mSet object
        position <- which(mSet@MSnResults[["Concensus_spec"]][[1]] == i)
        ms2_df <- mSet@MSnResults[["Concensus_spec"]][[2]][[position]][[1]]
        colnames(ms2_df) <- c("mz", "intensity")
        
        # Create a Spectrum object for MS2
        ms2_spec <- base$Spectrum(mz_array = np$array(ms2_df['mz']),
                                  int_array = np$array(ms2_df['intensity']))
        
        # Create a MetaFeature object (assuming the other parameters remain constant for all files)
        metafeature <- base$MetaFeature(identifier = i+1,
                                        mz = mSet@MSnData[["peak_mtx"]][[i+1]][[1]],
                                        rt = mSet@MSnData[["peak_mtx"]][[i+1]][[4]],
                                        charge = charge,
                                        ms2 = ms2_spec)
        
        # Append metafeature to the list
        features_list <- c(features_list, list(metafeature))
      }
    } else{
      for(i in mSet@MSnResults[["Concensus_spec"]][[1]]) {
        
        # get MS2 spectrum from the current mSet object
        position <- which(mSet@MSnResults[["Concensus_spec"]][[1]] == i)
        ms2_df <- mSet@MSnResults[["Concensus_spec"]][[2]][[position]][[1]]
        colnames(ms2_df) <- c("mz", "intensity")
        
        # Create a Spectrum object for MS2
        ms2_spec <- base$Spectrum(mz_array = np$array(ms2_df['mz']),
                                  int_array = np$array(ms2_df['intensity']))
        
        # Create a MetaFeature object (assuming the other parameters remain constant for all files)
        metafeature <- base$MetaFeature(identifier = i+1,
                                        mz = mSet@MSnData[["peak_mtx"]][[i+1]][[1]],
                                        rt = mSet@MSnData[["peak_mtx"]][[i+1]][[4]],
                                        charge = charge[i+1],
                                        ms2 = ms2_spec)
        
        # Append metafeature to the list
        features_list <- c(features_list, list(metafeature))
      }
    }
  }
  
  
  # Add features_list the Msbuddy object
  engine$add_data(features_list)
  
  # Annotate molecular formula
  engine$annotate_formula()
  
  # Retrieve the annotation result summary
  formula_result <-engine$get_summary()
  #formula_result <-as.data.frame(do.call(rbind, engine$get_summary()))
  
  # Function to replace NULL with NA and unlist
  prepare_row <- function(x) {
    x[sapply(x, is.null)] <- NA  # Replace NULL with NA
    unlist(x)  # Unlist the sublist into a vector
  }
  
  # Apply the function to each sublist and combine them into a data frame
  formula_result <- do.call(rbind, lapply(formula_result, prepare_row))
  
  # Convert the result into a data frame
  formula_result <- as.data.frame(formula_result, stringsAsFactors = FALSE)
  
  # Convert factors to characters if there are any
  formula_result[] <- lapply(formula_result, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  # Convert columns that should be numeric
  formula_result$identifier <- as.numeric(formula_result$identifier)
  formula_result$mz <- as.double(formula_result$mz)
  formula_result$rt <- as.double(formula_result$rt)
  formula_result$estimated_fdr <- as.double(formula_result$estimated_fdr)
  mSet@MSnResults[["FormulaPre"]] <- formula_result
  return(mSet)
}



##run asari process in R
##need to install asari lirary first
##then, need to replace two files 
run_process <- function(parameters, args, sample_class = NULL) {
  # main process function
  ##list_input_files <- asari$workflow$read_project_dir(args$input)
  if (class(parameters$mz_tolerance_ppm) != "integer") {
    parameters$mz_tolerance_ppm <- as.integer(round(parameters$mz_tolerance_ppm))
  }
  if (class(parameters$min_intensity_threshold) != "integer") {
    parameters$min_intensity_threshold <- as.integer(round(parameters$min_intensity_threshold))
  }
  if (class(parameters$min_peak_height) != "integer") {
    parameters$min_peak_height <- as.integer(round(parameters$min_peak_height))
  }
  if (class(parameters$min_timepoints) != "integer") {
    parameters$min_timepoints <- as.integer(round(parameters$min_timepoints))
  }
  if (class(parameters$signal_noise_ratio) != "integer") {
    parameters$signal_noise_ratio <- as.integer(round(parameters$signal_noise_ratio))
  }
  if (class(parameters$wlen) != "integer") {
    parameters$wlen <- as.integer(round(parameters$wlen))
  }
  
  list_input_files <- list.files(args$input, full.names = TRUE, recursive = TRUE, pattern = ".mzML")
  if(length(list_input_files) == 0) {
    cat("No valid mzML files are found in the input directory :(\n")
  } else {
    if(args$autoheight){
      tryCatch({
        parameters[['min_peak_height']] <- asari$analyze$estimate_min_peak_height(list_input_files)
      }, error = function(err) {
        cat(paste("Problems with input files:", err$message, "Back to default min_peak_height.\n"))
      })
    }
    parameters[['min_prominence_threshold']] <- as.integer(0.33 * parameters[['min_peak_height']])
    ##process project
    result <- process_project(list_input_files, parameters)
  }
  
  if (!is.null(sample_class)) {
    sample_class <- read.csv(sample_class_path)
    rownames(sample_class) <- sample_class[,1]
    peak_table <- result[["feature_tables"]][["preferred_table"]][,c(1,12:ncol(result[["feature_tables"]][["preferred_table"]]))]
    nrow <- nrow(peak_table)
    for (i in (2:ncol(peak_table))) {
      peak_table[nrow+1,i] <- sample_class[colnames(peak_table)[i],2]
    }
    peak_table <- rbind(peak_table[nrow(peak_table),], peak_table[-nrow(peak_table),])
    peak_table[1,1] <- "class"
    result$peak_table <- peak_table 
  }
  return(result)
}

## main function of run_process
process_project <- function(list_input_files, parameters) {
  
  sample_registry <- asari$workflow$register_samples(list_input_files)
  
  if (parameters['database_mode'] == 'auto') {
    if (length(list_input_files) <= parameters['project_sample_number_small']) {
      parameters['database_mode'] <- 'memory'
    } else {
      parameters['database_mode'] <- 'ondisk'
    }
  }
  
  time_stamp <- format(Sys.time(), "%m%d%H%M%S")
  parameters['time_stamp'] <- gsub("","",time_stamp)
  #asari$workflow$create_export_folders(parameters, time_stamp)
  
  shared_dict <- asari$workflow$batch_EIC_from_samples_(sample_registry, parameters)
  
  for (sid in names(sample_registry)) {
    sam <- sample_registry[[sid]]
    temp <- shared_dict[[sid]]
    names(temp) <- c('status:mzml_parsing', 'status:eic', 'data_location', 'max_scan_number', 'list_scan_numbers', 'list_retention_time', 'track_mzs', 'number_anchor_mz_pairs', 'anchor_mz_pairs', 'sample_data')
    sam <- c(sam, temp)
    sam['name'] <- gsub(".mzML", "", basename(sam[["input_file"]]))
    sample_registry[[sid]] <- sam
  }
  ##choose and change the reference sample based on the number_anchor_mz_pairs
  ##parameters$reference_sample_id <- get_reference_sample_id(parameters, sample_registry)
  parameters$reference <- get_reference_sample_id(parameters, sample_registry)
  

  
  
  EE <- asari$experiment$ext_Experiment(sample_registry, parameters)
  ##asari$experiment$ext_Experiment$process_all()
  result <- list()
  result <- EE$process_all()
  result <- EE$export_all(anno=parameters["anno"])
  
  if (!parameters[["pickle"]] & parameters[["database_mode"]] != 'memory') {
    remove_intermediate_pickles(parameters)
  }
  result$feature_tables$full_table <- py_to_r(result$feature_tables$full_table)
  result$feature_tables$preferred_table <- py_to_r(result$feature_tables$preferred_table)
  if (!is.null(result$feature_tables$targeted_table)) {
    result$feature_tables$targeted_table <- py_to_r(result$feature_tables$targeted_table)
    }
  result$feature_tables$unique_compound_table <- py_to_r(result$feature_tables$unique_compound_table)
  
  return(result)
}

## get the reference sample id from all the samples based on the number_anchor_mz_pairs
get_reference_sample_id <- function(parameters, sample_registry) {
  if(!is.null(parameters$reference)){
    for (k in names(sample_registry)){
      v <- sample_registry[[k]]
      if(basename(parameters$reference) == basename(v$input_file)){
        return(parameters$reference)
      }
    }
  } else if(length(sample_registry) > 0) {
    L <- lapply(names(sample_registry), function(k) 
      c(sample_registry[[k]]$number_anchor_mz_pairs, k)
    )
    L <- do.call(rbind, L)
    L <- L[order(L[,1], decreasing = TRUE), ]
    ref <- sample_registry[[L[1, 2]]]
    #parameters$reference <- ref$name
    ##parameters$reference <- ref$sample_id
    #parameters$reference <- ref$input_file
    cat("\n    The reference sample is:\n    ||* ", ref$name, " *||\n")
    cat("Max reference retention time is ", max(ref$list_retention_time), 
        " at scan number ", ref$max_scan_number, ".\n")
  }
  return(paste0(ref$name, ".mzML"))
}

PerformAsariResultsFormating <- function(minFrac = 0.7){
  
  alfs <- list.files(".", pattern = "results_metabo_asari_res")
  alfs_idx <- as.numeric(gsub("results_metabo_asari_res_","",alfs))
  result_folder <- alfs[which.max(alfs_idx)];

  # generate metaboanalyst table
  load("mSet.rda")
  load("params.rda")
  mSet@params <- updateRawSpectraParam (peakParams);
  ftable_ori <- ftable <- read.csv(paste0(result_folder, "/preferred_Feature_table.tsv"), sep = "\t", check.names=FALSE)
  features <- paste0(ftable$mz, "__", ftable$rtime)
  ftable1 <- ftable[,c(12:ncol(ftable))]
  allSamples <- colnames(ftable1)
  allGroups <- 
    vapply(allSamples, FUN = function(x){
      idx <- which(mSet@rawOnDisk@phenoData@data[["sample_name"]] == x)
      mSet@rawOnDisk@phenoData@data[["sample_group"]][idx]
      }, character(1L))
  ftable2 <- t(data.frame(Groups = allGroups))
  ftable3 <- data.frame(Samples = c("Groups", features))
  ftable0 <- rbind(ftable2, ftable1)
  ftable0 <- cbind(ftable3, ftable0)
  rownames(ftable0) <- NULL;
  mSet@dataSet <- ftable0;
  write.csv(ftable0, file = "metaboanalyst_input.csv", row.names = F, quote = F)
  
  # generate peak_feature_summary
  ftab_annotation <- read.csv(paste0(result_folder, "/Feature_annotation.tsv"), sep = "\t", check.names=FALSE)
  idx_num <- ftable$id_number
  idx_row <- vapply(idx_num, FUN = function(x){
    which(ftab_annotation[,1] == x)[1]
  }, FUN.VALUE = integer(1L))
  ftab_annotation <- ftab_annotation[idx_row, ]
  
  annots <- strsplit(ftab_annotation[,6], ",")
  adds <- vapply(annots, function(x){
    if(length(x) == 0){
      return("")
    } else {
      return(x[2])
    }
  }, FUN.VALUE = character(1L))
  isos <- vapply(annots, function(x){
    if(length(x) == 0){
      return("")
    } else {
      return(x[1])
    }
  }, FUN.VALUE = character(1L))
  #all_annots <- rjson::fromJSON(file = paste0(result_folder, "/Annotated_empiricalCompounds.json"))
  
  all_recrds <- ftab_annotation$matched_DB_records
  all_forumus <- sapply(all_recrds, function(x){
    if(is.na(x)){return("")}
    res <- strsplit(x, "\\), \\(")
    if(length(res[[1]]) == 0){
      return("")
    } else {
      res <- gsub("\\(|\\)", "", res[[1]])
      res2 <- vapply(res, FUN = function(y){
        strsplit(strsplit(y, "'")[[1]][2], "_")[[1]][1]
      }, character(1L), USE.NAMES = F)
      if(length(res2) == 1){
        res2 <- paste0(res2, ";")
      } else {
        res2 <- paste0(res2, collapse = "; ")
      }
      return(res2)
    }
  })
  
  all_cmpds <- ftab_annotation$matched_DB_shorts
  all_cmpd <- sapply(all_cmpds, function(x){
    if(is.na(x)){return("")}
    res <- strsplit(x, "\\), \\(")
    if(length(res[[1]]) == 0){
      return("")
    } else {
      res <- gsub("\\(|\\)", "", res[[1]])
      res2 <- lapply(res, FUN = function(y){
        w <- strsplit(y, ";")[[1]]
        strsplit(w, "\\$")
      })
      
      res2_done <- vapply(res2, function(nn){
        if(length(nn)==1){
          nn[[1]][2]
        } else {
          res2x <- vector(mode = "character", length = length(nn))
          for(kk in 1:length(nn)){
            res2x[kk] <- nn[[kk]][2]
          }
          return(paste0(res2x, collapse = "; "))
        }
      }, FUN.VALUE = character(1L))
      
      if(length(res2_done) == 1){
        res2_done <- paste0(res2_done, ";")
      } else {
        res2_done <- paste0(res2_done, collapse = "; ")
      }
      return(res2_done)
    }
  })
  all_hmdb <- sapply(all_cmpds, function(x){
    if(is.na(x)){return("")}
    res <- strsplit(x, "\\), \\(")
    if(length(res[[1]]) == 0){
      return("")
    } else {
      res <- gsub("\\(|\\)", "", res[[1]])
      res2 <- lapply(res, FUN = function(y){
        w <- strsplit(y, ";")[[1]]
        strsplit(w, "\\$")
      })
      
      res2_done <- vapply(res2, function(nn){
        if(length(nn)==1){
          nn[[1]][1]
        } else {
          res2x <- vector(mode = "character", length = length(nn))
          for(kk in 1:length(nn)){
            res2x[kk] <- nn[[kk]][1]
          }
          return(paste0(res2x, collapse = "; "))
        }
      }, FUN.VALUE = character(1L))
      
      if(length(res2_done) == 1){
        res2_done <- paste0(res2_done, ";")
      } else {
        res2_done <- paste0(res2_done, collapse = "; ")
      }
      return(res2_done)
    }
  })
  
  ftab_annotation$matched_DB_records[is.na(ftab_annotation$matched_DB_records)] <- ""
  
  Formula2Cmpd_list <- lapply(1:length(all_recrds), function(x){
    res <- list()
    if(ftab_annotation$matched_DB_records[x] == ""){
      return(res)
    }
    
    fms <- ftab_annotation$matched_DB_records[x]
    res <- strsplit(fms, "\\), \\(")
    res <- gsub("\\(|\\)", "", res[[1]])
    res2 <- vapply(res, FUN = function(y){
      strsplit(strsplit(y, "'")[[1]][2], "_")[[1]][1]
    }, character(1L), USE.NAMES = F)
    res2_ori <- res2
    res2 <- unique(res2)
    res <- rep(list(list(cmpd = vector(mode = "character"),
                    hmdb = vector(mode = "character"))), length(res2))
    
    names(res) <- res2
    
    this_cmpd <- ftab_annotation$matched_DB_shorts[x]
    #cmpd
    resx <- strsplit(this_cmpd, "\\), \\(")
    if(length(resx[[1]]) == 0){
      
    } else {
      resx <- gsub("\\(|\\)", "", resx[[1]])
      res2cmpd <- lapply(resx, FUN = function(y){
        w <- strsplit(y, ";")[[1]]
        strsplit(w, "\\$")
      })
      
      res2_done <- lapply(res2cmpd, function(nn){
        if(length(nn)==1){
          nn[[1]][2]
        } else {
          res2x <- vector(mode = "character", length = length(nn))
          for(kk in 1:length(nn)){
            res2x[kk] <- nn[[kk]][2]
          }
          return(res2x)
        }
      })
      
      for(nn in 1:length(res2_done)){
        res[[res2_ori[nn]]][["cmpd"]] <- unique(c(res[[res2_ori[nn]]][["cmpd"]], res2_done[[nn]]))
      }
    }
    #hmdb
    resy <- strsplit(this_cmpd, "\\), \\(")
    if(length(resy[[1]]) == 0){
      
    } else {
      resy <- gsub("\\(|\\)", "", resy[[1]])
      res2hmdb <- lapply(resy, FUN = function(y){
        w <- strsplit(y, ";")[[1]]
        strsplit(w, "\\$")
      })
      
      res2_done <- lapply(res2hmdb, function(nn){
        if(length(nn)==1){
          nn[[1]][1]
        } else {
          res2x <- vector(mode = "character", length = length(nn))
          for(kk in 1:length(nn)){
            res2x[kk] <- nn[[kk]][1]
          }
          return(res2x)
        }
      })
      
      for(nn in 1:length(res2_done)){
        res[[res2_ori[nn]]][["hmdb"]] <- unique(c(res[[res2_ori[nn]]][["hmdb"]], res2_done[[nn]]))
      }
    }
    
    return(res)
  })
  
  mSet@peakAnnotation[["Formula2Cmpd"]] <- Formula2Cmpd_list

  # peak_feature_summary
  LogNorm<-function(x, min.val){
    log10((x + sqrt(x^2 + min.val^2))/2)
  }
  CalCV<- function(x){
    x <- as.numeric(x)
    sd(x)/mean(x)
  }
  ftable1 -> data
  allGroups -> groups
  
  min.val <- min(abs(data[data!=0]))/10;
  data<-apply(data, 2, LogNorm, min.val);
  
  sample_data_log <- data;
  cvs <- round(apply(data, 1,FUN = CalCV),4)*100
  lvls <- groups[groups != "QC"];
  sample_data_log <- sample_data_log[,groups != "QC"];
  groups <- as.factor(lvls);
  
  ttest_res <- PerformFastUnivTests(t(sample_data_log), as.factor(groups))
  pvals <- ttest_res[,2]
  pfdr <-p.adjust(pvals, method = "fdr")
  pvals <- signif(pvals, 8)
  pfdr <- round(signif(pfdr, 8), 8)
  
  pvals[is.nan(pvals)] = 1
  pfdr[is.nan(pfdr)] = 1
  
  my.dat <- data.frame(mz = ftab_annotation$mz,
                            rt = ftab_annotation$rtime,
                            adduct = adds,
                            isotopes = isos,
                            Formula = all_forumus,
                            Compound = all_cmpd,
                            HMDBID = all_hmdb,
                            intenVale = ftable$peak_area,
                            pvals = pvals,
                            pfdr = pfdr,
                            cvs = cvs)
  
  HMDBIDs <- my.dat$HMDBID
  HMDBIDs[HMDBIDs==""] <- NA
  mSet@peakAnnotation[["massMatching"]] <- data.frame(Compound = my.dat$Compound, 
                                                      Formula = my.dat$Formula, 
                                                      HMDBID = HMDBIDs)
  FeatureOrder <- order(pvals)
  my.dat0 <- my.dat[FeatureOrder,]
  
  write.table(my.dat0, sep = "\t",
              file = "peak_feature_summary.tsv",
              row.names = FALSE,
              quote = FALSE);
  write.csv(my.dat0, 
            file = "peak_feature_summary.csv", 
            quote = TRUE, 
            row.names = FALSE);
  
  qs::qsave(mSet@peakAnnotation[["Formula2Cmpd"]][FeatureOrder], 
            file = "formula2cmpd.qs")
  
  # camera_output
  camera_output_df <- data.frame(mz = ftable$mz, 
                                 mzmin = ftable$mz-10*ftable$mz*1e-6,
                                 mzmax = ftable$mz+10*ftable$mz*1e-6,
                                 rt = ftable$rtime,
                                 rtmin = ftable$rtime_left_base,
                                 rtmax = ftable$rtime_right_base,
                                 ftable1,
                                 isotopes = my.dat$isotopes,
                                 adduct = my.dat$adduct,
                                 pcgroup = 1,
                                 Compound = my.dat$Compound, 
                                 Formula = my.dat$Formula, 
                                 HMDBID = HMDBIDs, check.names = FALSE)
  
  
  mSet@peakAnnotation[["camera_output"]] <- camera_output_df
  
  # generate peak_result_summary
  mzrange <- vapply(allSamples, FUN = function(x){
    ftab_annotation$mz -> mzs;
    idx <- which(colnames(ftable1) == x)
    min_mz <- round(min(mzs[ftable1[,idx] != 0]), 4)
    max_mz <- round(max(mzs[ftable1[,idx] != 0]), 4)
    paste0(min_mz, "~", max_mz)
  }, FUN.VALUE = character(1L));
  
  rtrange <- vapply(allSamples, FUN = function(x){
    ftab_annotation$rtime -> rts;
    idx <- which(colnames(ftable1) == x)
    min_rt <- round(min(rts[ftable1[,idx] != 0]), 2)
    max_rt <- round(max(rts[ftable1[,idx] != 0]), 2)
    paste0(min_rt, "~", max_rt)
  }, FUN.VALUE = character(1L))
  
  peak_num <- vapply(allSamples, FUN = function(x){
    idx <- which(colnames(ftable1) == x)
    length(which(ftable1[,idx] != 0))
  }, FUN.VALUE = integer(1L))
  
  peak_ratio_missing <- vapply(allSamples, FUN = function(x){
    idx <- which(colnames(ftable1) == x)
    round(1-length(which(ftable1[,idx] != 0))/nrow(ftable1),4)*100
  }, FUN.VALUE = double(1L))
  
  datam <- data.frame(samples = allSamples,
                      groups = allGroups,
                      rtrange = rtrange,
                      mzrange = mzrange,
                      peak_num = peak_num,
                      missing = peak_ratio_missing)
  
  mSet@peakAnnotation[["peak_result_summary"]] <- as.matrix(datam)
  
  write.table(
    datam,
    file = paste0("peak_result_summary.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  );
  save(mSet, file = "mSet.rda")
  
  # Generate PCA figure and intensity stats
  {
    imgName <- "PCA.png";
    dpi = 72;
    width = 8;
    format = "png"
    if (.on.public.web()) {
      Cairo::Cairo(
        file = imgName,
        unit = "in",
        dpi = dpi,
        width = width,
        height = width * 0.80,
        type = format,
        bg = "white"
      )
    }
    
    sample_idx <- allGroups;
    
    if(any(sample_idx == "BLANK")){
      sample_idx_nonblk <- sample_idx != "BLANK"
      sample_idx <- sample_idx[sample_idx_nonblk]
      ftable1 <- ftable1[,sample_idx_nonblk]
    }
    
    feature_value0 <- ftable1;
    rownames(feature_value0) <- features
    feature_value <- feature_value0
    feature_value[is.na(feature_value)] <- 0;
    
    int.mat <- as.matrix(feature_value)
    rowNms <- rownames(int.mat);
    colNms <- colnames(int.mat);
    int.mat <- t(apply(int.mat, 1, function(x) .replace.by.lod(as.numeric(x))));
    int.mat[int.mat == Inf] <- 0
    rownames(int.mat) <- rowNms;
    colnames(int.mat) <- colNms; 
    feature_value <- int.mat;
    feature_value[feature_value==0] <- 1;
    
    pca_feats <- log10(feature_value);
    
    if (nrow(feature_value) < 2) {
      MessageOutput(
        mes =  paste0(
          "<font color=\"red\">",
          "\nERROR: No enough peaks detected, please adjust your parameters or use other Peak/Alignment method",
          "</font>"
        ),
        ecol = "\n",
        progress = 65
      );
      
      if (.on.public.web()) {
        dev.off()
      }
      
      return(NULL)
    }
    
    pca_feats[is.na(pca_feats)] <- 0;
    df0 <- na.omit(pca_feats);
    df1 <- df0[is.finite(rowSums(df0)),];
    df <- t(df1);
    
    mSet_pca <- prcomp(df, center = TRUE, scale = FALSE);
    sum.pca <- summary(mSet_pca);
    var.pca <-
      sum.pca$importance[2,]; # variance explained by each PCA
    
    xlabel <- paste("PC1", "(", round(100 * var.pca[1], 1), "%)");
    ylabel <- paste("PC2", "(", round(100 * var.pca[2], 1), "%)");
    zlabel <- paste("PC3", "(", round(100 * var.pca[3], 1), "%)");
    # using ggplot2
    df <- as.data.frame(mSet_pca$x);
    df$group <- sample_idx;
    
    ## Handle to generate json file for PCA3D online
    if(.on.public.web()){
      ## For score plot
      pca3d <- list();
      pca3d$score$axis <- c(xlabel, ylabel, zlabel);
      xyz0 <- df[,seq_len(3)];
      colnames(xyz0) <- rownames(xyz0) <- NULL;
      pca3d$score$xyz <- data.frame(t(xyz0));
      colnames(pca3d$score$xyz) <- NULL;
      pca3d$score$name <- rownames(df);
      pca3d$score$facA <- df$group;
      
      if(length(unique(df$group)) < 9){
        col.fun <-
          grDevices::colorRampPalette(RColorBrewer::brewer.pal(length(unique(df$group)), "Set1"));
      } else {
        col.fun <-
          grDevices::colorRampPalette(RColorBrewer::brewer.pal(length(unique(df$group)), "Set3"));
      }
      
      pca3d$score$colors <- col.fun(length(unique(df$group)));
      
      ## For loading plot
      
      pca3d$loading$axis <- paste("Loading ", seq_len(3), sep="");
      coords0 <- coords <- data.frame(t(signif(mSet_pca$rotation[,seq_len(3)], 5)));
      colnames(coords) <- NULL; 
      pca3d$loading$xyz <- coords;
      pca3d$loading$name <- rownames(mSet_pca$rotation);
      pca3d$loading$entrez <- paste0(gsub("__", "@", features));
      
      dists <- GetDist3D(coords0);
      pca3d$loading$cols <- GetRGBColorGradient(dists);
      
      pca3d$cls =  df$group;
      
      # json.obj <- RJSONIO::toJSON(pca3d, .na='null');
      # sink("spectra_3d_loading.json");
      # cat(json.obj);
      # sink();
      qs::qsave(pca3d$score, "score3d.qs");
      qs::qsave(pca3d$loading, "loading3d.qs");
      fileNm <- paste0("spectra_3d_loading.json");
      
      my.json.scatter(fileNm, T);
      
    }
    
    if (nrow(df) < 30) {
      if (length(unique(sample_idx)) > 9) {
        col.fun <-
          grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"));
        
        p <-
          ggplot2::ggplot(df, aes_string(
            x = "PC1",
            y = "PC2",
            color = "group",
            label = "row.names(df)"
          )) +
          geom_text_repel(force = 1.5) + 
          geom_point(size = 5,  fill = col.fun(length(unique(sample_idx)))) + 
          theme(axis.text = element_text(size = 12))
        
      } else{
        p <-
          ggplot2::ggplot(df, aes_string(
            x = "PC1",
            y = "PC2",
            color = "group",
            label = "row.names(df)"
          )) +
          geom_text_repel(force = 1.5) + 
          geom_point(size = 5) + 
          scale_color_brewer(palette = "Set1") + 
          theme(axis.text = element_text(size = 12))
      }
      
    } else {
      if (length(unique(sample_idx)) > 9) {
        p <-
          ggplot2::ggplot(df, aes(x = "PC1",
                                  y = "PC2",
                                  color = "group")) + geom_point(size = 5)
        
      } else{
        p <-
          ggplot2::ggplot(df, aes(x = "PC1",
                                  y = "PC2",
                                  color = "group")) + geom_point(size = 5) + scale_color_brewer(palette = "Set1");
      }
    }
    
    p <-
      p + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(axis.title = element_text(size = 12));
    
    print(p)
    
    if (.on.public.web()) {
      dev.off()
    }
    
    # process peak intensity
    sample_idx <- mSet@rawOnDisk@phenoData@data[["sample_group"]];
    
    sample_num <-
      mSet@rawOnDisk@phenoData@data[["sample_name"]];
    
    sample_num <- sample_num[sample_idx!="MS2"]
    sample_idx <- sample_idx[sample_idx!="MS2"]
    
    if (length(unique(sample_idx)) > 9) {
      col.fun <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      group_colors <- col.fun(length(unique(sample_idx)))
      
    } else{
      group_colors <-
        paste0(RColorBrewer::brewer.pal(9, "Set1")[seq_along(unique(sample_idx))], "60")
    }
    if(any(sample_idx == "BLANK")){
      sample_num <- sample_num[sample_idx != "BLANK"]
      sample_idx <- sample_idx[sample_idx!="BLANK"]
    }
    ftable1[ftable1 == 0] <- 1
    ints <- lapply(sample_num, function(x){
      log2(ftable1[,x])
    })
    #ints <- ints[names(ints) != 0];
    names(ints) <- as.character(sample_num)
    
    sample_idx <- as.factor(sample_idx)
    group_colors <-
      sapply(
        seq(length(levels(sample_idx))),
        FUN = function(x) {
          rep(group_colors[x], length(sample_idx[sample_idx == levels(sample_idx)[x]]))
        }
      )
    
    if(!.on.public.web()){
      oldpar <- par(no.readonly = TRUE);
      on.exit(par(oldpar));
    }
    imgName <- "Peak_Intensity.png"
    if (.on.public.web()) {
      Cairo::Cairo(
        file = imgName,
        unit = "in",
        dpi = dpi,
        width = width,
        height = length(sample_num) * 0.65 * (width/8),
        type = format,
        bg = "white"
      )
    }
    
    #op <- 
    par(mar = c(3.5, 10, 4, 1.5), xaxt = "s")
    
    sampleNMs <- names(ints);
    len_nms <- nchar(sampleNMs);
    if(any(len_nms > 15)){
      names(ints) <- 
        unname(unlist(sapply(sampleNMs, function(x){
          LEN_x <- nchar(x);
          if(LEN_x > 15){
            substring(x, LEN_x-14,LEN_x)
          } else {
            x
          }
        })))
    }
    
    boxplot(
      ints,
      varwidth = TRUE,
      col = as.character(unlist(group_colors)),
      ylab = "",
      horizontal = TRUE,
      las = 2,
      main = expression(log[2] ~ intensity),
      cex.lab = 0.8,
      cex.main = 1.25
    )
    
    #title(ylab=expression(log[2]~intensity), line=7.5, cex.lab=1.2)
    grid(nx = NA, ny = NULL)
    
    if (.on.public.web()) {
      dev.off()
    }
    
  }

  # generate clean table and object (optionally) for ms2
  {
    annotation_list <- rjson::fromJSON(file = paste0(alfs, "/Annotated_empiricalCompounds.json"))
    unique_ft_anns <- unique(ftab_annotation$`[EmpCpd]interim_id`)
    
    # filter adducts/isotopes
    ft_res <- vapply(unique_ft_anns, function(x){
      idxs <- which(ftab_annotation$`[EmpCpd]interim_id` == x)
      rltions <- ftab_annotation$`[EmpCpd]ion_relation`[idxs]
      
      parental_ion_idx <- which(grepl("M0,M\\+H\\+|M0,M\\-H\\-", rltions))
      ft_final_idx <- integer(1L)
      if(length(parental_ion_idx)== 0){
        return(ftab_annotation$`[peak]id_number`[idxs[1]])
      } else if(length(parental_ion_idx)>1){
        parental_ion_idx <- parental_ion_idx[1]
      }
      ft_final_idx <- idxs[parental_ion_idx]
      return(ftab_annotation$`[peak]id_number`[ft_final_idx])
    }, FUN.VALUE = character(1L))
    idx_keep1 <- which(ftable_ori$id_number %in% ft_res)
    ftable_f1 <- ftable_ori[idx_keep1, ]
    ftab_annotation_f1 <- ftab_annotation[idx_keep1, ]
    
    # filter QC
    group_info <- mSet@rawOnDisk@phenoData@data[["sample_group"]];
    
    if("QC" %in% group_info){
      qc_names <- names(which(allGroups == "QC"))
      qc_idx <- vapply(qc_names, function(y) {which(y == colnames(ftable_f1))}, integer(1L))
      ma_feats_miss0 <- ftable_f1[,qc_idx]
      idx_keep2 <- apply(ma_feats_miss0, 1, function(x){
        qc_vec <- as.numeric(x)
        rsd <- sd(qc_vec, na.rm = T)/mean(qc_vec, na.rm = T)
        if(is.na(rsd)){
          return(FALSE)
        }
        return(rsd<=1)
      })
      ftable_f2 <- ftable_f1[idx_keep2, ]
      ftab_annotation_f2 <- ftab_annotation_f1[idx_keep2, ]
    }

    # add empirical label
    features <- paste0(ftable_f2$mz, "__", ftable_f2$rtime)
    ftable_f2_1 <- ftable_f2[,c(12:ncol(ftable_f2))]
    
    empir_cmpds <- apply(ftab_annotation_f2, 1, function(x){
      emp1 <- as.character(x[length(x)])
      if(emp1==""){
        return("")
      } else {
        emp_id <- as.character(x[5])
        
        emp_res <- annotation_list[[emp_id]][["list_matches"]][[1]][[1]]
        if(is.null(emp_res)){
          return("")
        } else {
          return(strsplit(emp_res,"_")[[1]][1])
        }
      }
    })
    features[empir_cmpds != ""] <- paste0(empir_cmpds[empir_cmpds !=""],"_",features[empir_cmpds != ""])
    
    allSamples <- colnames(ftable_f2_1)
    allGroups <- 
      vapply(allSamples, FUN = function(x){
        idx <- which(mSet@rawOnDisk@phenoData@data[["sample_name"]] == x)
        mSet@rawOnDisk@phenoData@data[["sample_group"]][idx]
      }, character(1L))
    ftable_f2_2 <- t(data.frame(Groups = allGroups))
    ftable_f2_3 <- data.frame(Samples = c("Groups", features))
    ftable_f2_0 <- rbind(ftable_f2_2, ftable_f2_1)
    ftable_f2_0 <- cbind(ftable_f2_3, ftable_f2_0)
    rownames(ftable_f2_0) <- NULL;
    
    ms1_res_dt <- list(ftable_f2_0, ftab_annotation_f2, ftable_f2)
    write.csv(ftable_f2_0, file = "metaboanalyst_input_clean.csv", row.names = F, quote = F)
    qs::qsave(ms1_res_dt, file = "metaboanalyst_input_clean_MS1_asari.qs")
    
  }
  
  
  # restore the mSet Obj for other tasks
  save(mSet, file = "mSet.rda")
  MessageOutput(
    mes =  paste0(
      "<font color=\"blue\">",
      "\nRaw spectral data processing has been finished completely!",
      "</font>"
    ),
    ecol = "\n",
    progress = 100
  );
  
  
}

findMostAbundantAsariFeature <- function(){
  
  
}

PerformFastUnivTests <- function(data, cls, var.equal=TRUE){
  if(!exists("mem.univ")){
    require("memoise");
    mem.univ <<- memoise(.perform.fast.univ.tests);
  }
  return(mem.univ(data, cls, var.equal));
}

.perform.fast.univ.tests <- function(data, cls, var.equal=TRUE){
  
  print("Performing fast univariate tests ....");
  # note, feature in rows for gene expression
  data <- t(as.matrix(data));
  if(length(levels(cls)) > 2){
    res <- try(rowcolFt(data, cls, var.equal = var.equal));
  }else{
    res <- try(rowcoltt(data, cls, FALSE, 1L, FALSE));
  }  
  
  if(class(res) == "try-error") {
    res <- cbind(NA, NA);
  }else{
    # res <- cbind(res$statistic, res$p.value);
    # make sure row names are kept
    res <- res[, c("statistic", "p.value")];
  }
  
  return(res);
}

rowcolFt =  function(x, fac, var.equal, which = 1L) {
  
  if(!(which %in% c(1L, 2L)))
    stop(sQuote("which"), " must be 1L or 2L.")
  
  if(which==2L)
    x = t(x)
  
  if (typeof(x) == "integer")
    x[] <- as.numeric(x)
  
  sqr = function(x) x*x
  
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]
  
  ## Number of levels (groups)
  k <- nlevels(fac)
  
  ## xm: a nrow(x) x nlevels(fac) matrix with the means of each factor
  ## level
  xm <- matrix(
    sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE])),
    nrow = nrow(x),
    ncol = nlevels(fac))
  
  ## x1: a matrix of group means, with as many rows as x, columns correspond to groups 
  x1 <- xm[,fac, drop=FALSE]
  
  ## degree of freedom 1
  dff    <- k - 1
  
  if(var.equal){
    ## x0: a matrix of same size as x with overall means
    x0 <- matrix(rowMeans(x), ncol=ncol(x), nrow=nrow(x))
    
    ## degree of freedom 2
    dfr    <- ncol(x) - dff - 1
    
    ## mean sum of squares
    mssf   <- rowSums(sqr(x1 - x0)) / dff
    mssr   <- rowSums(sqr( x - x1)) / dfr
    
    ## F statistic
    fstat  <- mssf/mssr
    
  } else{
    
    ## a nrow(x) x nlevels(fac) matrix with the group size  of each factor
    ## level
    ni <- t(matrix(tapply(fac,fac,length),ncol=nrow(x),nrow=k))
    
    ## wi: a nrow(x) x nlevels(fac) matrix with the variance * group size of each factor
    ## level
    sss <- sqr(x-x1)
    x5 <- matrix(
      sapply(levels(fac), function(fl) rowSums(sss[,which(fac==fl), drop=FALSE])),
      nrow = nrow(sss),
      ncol = nlevels(fac))          
    wi <- ni*(ni-1) /x5
    
    ## u : Sum of wi
    u  <- rowSums(wi)
    
    ## F statistic
    MR <- rowSums(sqr((1 - wi/u)) * 1/(ni-1))*1/(sqr(k)-1)
    fsno <- 1/dff * rowSums(sqr(xm - rowSums(wi*xm)/u) * wi)
    fsdeno <- 1+ 2* (k-2)*MR
    fstat <- fsno/fsdeno
    
    ## degree of freedom 2: Vector with length nrow(x)
    dfr <- 1/(3 * MR)
    
  }
  
  res = data.frame(statistic = fstat,
                   p.value   = pf(fstat, dff, dfr, lower.tail=FALSE),
                   row.names = rownames(x))
  
  attr(res, "df") = c(dff=dff, dfr=dfr)
  return(res)
}

rowcoltt =  function(x, fac, tstatOnly, which, na.rm) {
  
  #if(.on.public.web){
  #  dyn.load(.getDynLoadPath());
  #}
  
  if (!missing(tstatOnly) && (!is.logical(tstatOnly) || is.na(tstatOnly)))
    stop(sQuote("tstatOnly"), " must be TRUE or FALSE.")
  
  f = checkfac(fac)
  if ((f$nrgrp > 2) || (f$nrgrp <= 0))
    stop("Number of groups is ", f$nrgrp, ", but must be >0 and <=2 for 'rowttests'.")
  
  if (typeof(x) == "integer")
    x[] <- as.numeric(x)
  
  #cc = .Call("rowcolttests", x, f$fac, f$nrgrp, which-1L, na.rm)
  cc = XiaLabCppLib::rowcolttestsR(x, f$fac, f$nrgrp, which-1L, na.rm)
  
  res = data.frame(statistic = cc$statistic,
                   dm        = cc$dm,
                   row.names = dimnames(x)[[which]])
  
  if (!tstatOnly)
    res = cbind(res, p.value = 2*pt(abs(res$statistic), cc$df, lower.tail=FALSE))
  
  attr(res, "df") = cc$df    
  return(res)
}

checkfac = function(fac) {
  
  if(is.numeric(fac)) {
    nrgrp = as.integer(max(fac, na.rm=TRUE)+1)
    fac   = as.integer(fac)
  }
  ## this must precede the factor test
  if(is.character(fac))
    fac = factor(fac)
  
  if (is.factor(fac)) {
    nrgrp = nlevels(fac)
    fac   = as.integer(as.integer(fac)-1)
  } 
  if(!is.integer(fac))
    stop("'fac' must be factor, character, numeric, or integer.")
  
  if(any(fac<0, na.rm=TRUE))
    stop("'fac' must not be negative.")
  
  return(list(fac=fac, nrgrp=nrgrp))
}

jobFinished <- function(mSet = NA){
  try(exportCompoundTable(mSet), silent = T)
  MessageOutput(
    mes = paste0('<b>Everything of this LC-MS/MS dataset has been completed successfully! </b>\n\n'),
      ecol = '',
      progress = 200
    )
}

exportCompoundTable <- function(mSet){
  dt <- mSet@dataSet
  ms1_idxs <- mSet@MSnData[["peak_mtx_idx"]]
  dt_cmpd <- mSet@MSnResults[["ResTable"]]
  cmpd_idx <- vapply(1:nrow(dt_cmpd), function(x){
    mz_min <- dt_cmpd$mzmin[x]
    mz_max <- dt_cmpd$mzmax[x]
    rt_min <- dt_cmpd$rtmin[x]
    rt_max <- dt_cmpd$rtmax[x]
    which(mSet@peakAnnotation[["camera_output"]]$mzmin >= mz_min-0.0001 & 
            mSet@peakAnnotation[["camera_output"]]$mzmax <= mz_max+0.0001 & 
            mSet@peakAnnotation[["camera_output"]]$rtmin >= rt_min-0.1 & 
            mSet@peakAnnotation[["camera_output"]]$rtmax <= rt_max+0.1)[1]
  }, integer(1L))
  dt_cpmd_ms1 <- dt[c(1, cmpd_idx+1), ]
  cmpd_labels <- dt_cmpd$Compound_1
  dt_cpmd_ms1[-1, 1] <- cmpd_labels
  write.csv(dt_cpmd_ms1, file = "metaboanalyst_input_compound.csv", row.names = F, quote = T)
  dt[c(cmpd_idx+1), 1] <- cmpd_labels
  write.csv(dt, file = "metaboanalyst_input_integrated.csv", row.names = F, quote = T)
}

PerformExpsomeClassify <- function(mSet, path_repo = ""){
  if(path_repo == ""){
    stop("No classification database provided!")
  }
  
  exposome_repo <- qs::qread("/home/glassfish/projects/exposome_lib/complte_exposome_categories_lib.qs")
  anno_res <- mSet@MSnResults[["DBAnnoteRes"]]
  exposome_repo <- as.data.frame(exposome_repo)
  
  exposome_res <- lapply(anno_res, function(x){
    res1 <- x[[1]]
    if(length(res1$IDs)==0){
      return(NA)
    }
    all_inchikeys <- res1[["InchiKeys"]]
    nms <- colnames(exposome_repo)[-1]
    all_cls <- lapply(all_inchikeys, function(u){
      idx <- which(exposome_repo$InChiKeys == u)
      if(length(idx)!=0){
        nms[as.logical(exposome_repo[idx[1], -1])]
      } else {
        return(NA)
      }
    })
    names(all_cls) <- all_inchikeys
    return(all_cls)
  })
  
  mSet@MSnResults[["ExposomeRes"]] <- exposome_res
  
  exposome_res_clean <- lapply(exposome_res, function(x){
    if(length(x)==1){
      if(is.na(x)){
        return(NA)
      } else {
        return(x[[1]])
      }
    } else {
      return(x[[1]])
    }
  })
  
  # plot a summary histogram
  
  # prepare df_all
  if(nrow(mSet@dataSet) == 0){
    # reload mSet if MS1 information is missing
    mSet <- reloadMS1mSet(mSet);
  }
  #save(mSet, file = "mSet_PerformExpsomeClassify.rda")
  
  ft_idx <- mSet@MSnResults[["Concensus_spec"]][[1]]+1
  meta_info0 <- meta_info <- as.character(mSet@dataSet[1, -1])
  meta_info <- unique(meta_info)
  meta_info <- meta_info[meta_info!="QC" & meta_info!="MS2"]
  
  idx2check <- vapply(exposome_res_clean, function(m){
    is.na(m[1])
  }, FUN.VALUE = logical(1L))
  ft_idx <- ft_idx[!idx2check]
  ft_idx_seq <- seq(1:length(ft_idx))
  
  exposome_res_clean <- exposome_res_clean[!idx2check]
  
  dt <- mSet@dataSet[-1, -1]
  res_exp_class_by_group <- lapply(meta_info, FUN = function(n){
    idx_grp_col <- which(meta_info0 == n)
    bool_idxs <- vapply(ft_idx, function(x){
      this_ft_grp <- dt[x, idx_grp_col]
      length(which(this_ft_grp==0))/length(this_ft_grp)<= 0.75
    }, FUN.VALUE = logical(1L))
    this_good_ft_idx <- ft_idx_seq[bool_idxs]
    res_all_this_pho <- exposome_res_clean[this_good_ft_idx]
    res_all_this_pho <- unlist(res_all_this_pho)
    return(res_all_this_pho)
  })
  
  names(res_exp_class_by_group) <- meta_info
  
  all_cls <- c("Biocides", "Drugs", "Environment_Contaminantes", "Foods", 
               "Indoor_Environment", "Industrial_Toxins", "Microbes", "Natural_Toxins",
               "Neuro_Toxins", "Other_Hazards", "Other_Healths_Related", "Others",
               "Personal_Cares", "PFAS", "Phenols", "Plants",
               "Plastics", "PMTs", "Smokes", "Surfactants",
               "Waters")
  
  
  all_cls_grps <- lapply(res_exp_class_by_group, function(x){
    rs1 <- vapply(all_cls, function(y){
      length(which(y==x))
    }, FUN.VALUE = integer(1L))
    return(rs1)
  })
  
  df_all1 <- lapply(1:length(all_cls_grps), function(z) {data.frame(Categories = all_cls, Number = as.numeric(all_cls_grps[[z]]), Group = names(all_cls_grps)[z])})
  df_all <- do.call(rbind, df_all1)
  
  qs::qsave(df_all, file = "exposome_classification_summary.qs")
  
  require("viridis") 
  require("ggplot2")
  if(length(all_cls_grps)==2){
    fill <- c("#0000FF", "#FF0000")
  } else {
    fill <- RColorBrewer::brewer.pal(length(all_cls_grps), "Set2")
  }
  
  p4 <- ggplot(data = df_all) + 
    geom_bar(aes(y = Number, x = Categories, fill = Group), stat="identity", width = 0.5, position = "dodge") + theme_light() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 8.5))
  p4 <- p4 + scale_fill_manual(values=fill) 
  #p4 <- p4 + scale_color_viridis(discrete = TRUE, option = "D")+scale_fill_viridis(discrete = TRUE)
  p4
  
  Cairo::Cairo(820, 500,file = paste0("exposome_cat_comparison.png"),dpi = 90,bg = "white")
  print(p4)
  dev.off()
  
  Cairo::CairoPDF(width = 8, height = 5,file = paste0("exposome_cat_comparison.pdf"))
  print(p4)
  dev.off()
  
  return(mSet)
}

PerformMetabolomeClassify <- function(mSet, path_repo = ""){
  #save(mSet, file = "PerformMetabolomeClassify_mSet.rda")
  if(path_repo == ""){
    stop("No classification database provided!")
  }
  
  metabolome_repo <- qs::qread("/home/glassfish/projects/metabolome_lib/complete_metabolome_taxonomies_lib.qs")
  anno_res <- mSet@MSnResults[["DBAnnoteRes"]]
  metabolome_repo <- as.data.frame(metabolome_repo)
  
  metabolome_res <- lapply(anno_res, function(x){
    res1 <- x[[1]]
    if(length(res1$IDs)==0){
      return(NA)
    }
    all_inchikeys <- res1[["InchiKeys"]]
    nms <- colnames(metabolome_repo)[-1]
    all_cls <- lapply(all_inchikeys, function(u){
      idx <- which(metabolome_repo$InChiKeys == u)
      if(length(idx)!=0){
       metabolome_repo[idx[1], -1]
      } else {
        return(NA)
      }
    })
    names(all_cls) <- all_inchikeys
    return(all_cls)
  })

  mSet@MSnResults[["MetabolomeRes"]] <- metabolome_res

  metabolome_res_clean <- lapply(metabolome_res, function(x){
    if(length(x)==1){
      if(is.na(x)){
        return(NA)
      } else {
        return(x[[1]])
      }
    } else {
      return(x[[1]])
    }
  })
  
  # prepare df_all
  if(nrow(mSet@dataSet) == 0){
    # reload mSet if MS1 information is missing
    mSet <- reloadMS1mSet(mSet);
  }
  
  #
  ft_idx <- mSet@MSnResults[["Concensus_spec"]][[1]]+1
  meta_info0 <- meta_info <- as.character(mSet@dataSet[1, -1])
  meta_info <- unique(meta_info)
  meta_info <- meta_info[meta_info!="QC" & meta_info!="MS2"]
  
  idx2check <- vapply(metabolome_res_clean, function(m){
    is.na(m[1])
  }, FUN.VALUE = logical(1L))
  ft_idx <- ft_idx[!idx2check]
  ft_idx_seq <- seq(1:length(ft_idx))
  
  metabolome_res_clean <- metabolome_res_clean[!idx2check]
  
  dt <- mSet@dataSet[-1, -1]
  res_exp_class_by_group <- lapply(meta_info, FUN = function(n){
    
    idx_grp_col <- which(meta_info0 == n)
    bool_idxs <- vapply(ft_idx, function(x){
      this_ft_grp <- dt[x, idx_grp_col]
      length(which(this_ft_grp==0))/length(this_ft_grp)<= 0.75
    }, FUN.VALUE = logical(1L))
    
    this_good_ft_idx <- ft_idx_seq[bool_idxs]
    res_all_this_pho <- metabolome_res_clean[this_good_ft_idx]
    res_all_this_pho_done <- do.call(rbind, res_all_this_pho)
    return(res_all_this_pho_done)
  })
  
  names(res_exp_class_by_group) <- meta_info
  qs::qsave(res_exp_class_by_group, file = "metabolome_classification_summary.qs")
  
  return(mSet)
  
}

SummarizeAllResults4Reference <- function(mSet){
  save(mSet, file = "mSet_SummarizeAllResults4Reference.rda");
  
  ms1_dt <- mSet@peakAnnotation[["camera_output"]];
  ms2_dt <- qs::qread("compound_msn_results_index.qs")
  ms2_dtx <- as.data.frame(matrix(NA, nrow = nrow(ms1_dt), ncol = 25))
  for(i in 1:nrow(ms2_dt)){
    ms2_dtx[ms2_dt$peak_idx_vals[i], ] <- ms2_dt[i, -c(1:5)]
  }
  colnames(ms2_dtx) <- colnames(ms2_dt)[ -c(1:5)]
  
  ms1_2_dt <- cbind(ms1_dt, ms2_dtx)
  
  ft_idx <- mSet@MSnResults[["Concensus_spec"]][[1]]+1
  
  # add metabolome info
  mSet@MSnResults[["MetabolomeRes"]] -> metabolome_res
  metabolome_dtx <- as.data.frame(matrix(NA, nrow = nrow(ms1_dt), ncol = 4))
  for(i in 1:nrow(ms2_dt)){
    idx <- which(ms2_dt$peak_idx_vals[i] == ft_idx)
    metabolome_dtx[ms2_dt$peak_idx_vals[i], ] <- metabolome_res[[idx]][[1]][2:5]
  }
  colnames(metabolome_dtx) <- c("Kingdoms", "Super_classes", "Classes", "Sub_classes")
  # add metabolome info
  mSet@MSnResults[["ExposomeRes"]] -> exposome_res
  exposome_dtx <- as.data.frame(matrix(NA, nrow = nrow(ms1_dt), ncol = 1))
  for(i in 1:nrow(ms2_dt)){
    idx <- which(ms2_dt$peak_idx_vals[i] == ft_idx)
    exposome_dtx[ms2_dt$peak_idx_vals[i], 1] <- paste(exposome_res[[idx]][[1]], collapse = "; ")
  }
  colnames(exposome_dtx) <- c("Exposome_classes")
  
  metabo_feat_ref <- cbind(ms1_2_dt, metabolome_dtx, exposome_dtx)
  write.csv(metabo_feat_ref, file = "metaboanalyst_feature_reference.csv", row.names = F, quote = T)
  
  return(mSet)
}


reloadMS1mSet <- function(mSet1){
  
  load("mSet.rda")
  mSet1@peakAnnotation <- mSet@peakAnnotation
  mSet1@dataSet <- mSet@dataSet
  return(mSet1)
}