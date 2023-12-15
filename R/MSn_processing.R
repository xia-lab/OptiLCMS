#' PerformMSnImport
#'
#' @param mSet mSet
#' @param filesPath filesPath
#' @param targetFeatures targetFeatures
#' @param acquisitionMode acquisitionMode
#' @param SWATH_file SWATH_file
#'
#' @return mSet Object
#' @export
#'
#' @examples to add
PerformMSnImport <- function(mSet = NULL, 
                             filesPath = NULL, 
                             targetFeatures = NULL,
                             acquisitionMode = "DDA",
                             SWATH_file = NULL){
  require(mzR)
  
  if(is.null(mSet) & is.null(targetFeatures)){
    warning("Neither mSet object nor targetFeatures matrix provided! We will initialize one!")
  }
  if(!is.null(mSet) & !is.null(targetFeatures)){
    warning("Both mSet object and targetFeatures matrix provided! We will use targetFeatures matrix here!")
  }
  if(acquisitionMode != "DDA" & acquisitionMode != "DIA"){
    stop("acquisitionMode must be one of the \"DDA\" and \"DIA\".")
  }
  if(acquisitionMode == "DIA" & is.null(mSet) & is.null(targetFeatures)){
    stop("For DIA data, you must provide a targeted feature list (targetFeatures) or MS1 data processing results (mSet).")
  }
  
  ## Read raw data by using mzR package
  message("Importing raw MS data .. ", appendLF = F);
  # rawData1 <- readMSData(files = filesPath, msLevel. = 1L, mode = "inMemory");
  # rawData2 <- readMSData(files = filesPath, msLevel. = 2L, mode = "inMemory");
  # 
  msFiles_list <- lapply(filesPath, function(x){mzR::openMSfile(x, backend = "pwiz")})
  pks_list <- lapply(msFiles_list, mzR::spectra)
  hdr_list <- lapply(msFiles_list, mzR::header)
  message("done.")
  # Process ion mode (1 is positive, 0 is negative)
  all_polarities <- sapply(1:length(filesPath), function(x){
    unique(hdr_list[[x]][["polarity"]])
  })
  ion_mode <- unique(unlist(all_polarities))
  if(length(ion_mode)>1){
    stop("Scans with mutiple polarities are not allowed! Please keep one polarity at a time!")
  }
  
  # Split scanrts_ms1, scanrts_ms2, scan_ms1, scan_ms2 based on files' information
  scanrts_ms1_List = 
    scanrts_ms2_List = 
    scan_ms1_List = 
    scan_ms2_List = 
    prec_mzs_List = list();
  
  # Process the data of MS1 level
  scanrts_ms1_List <- lapply(1:length(filesPath), function(x){
    hdr_list[[x]][["retentionTime"]][hdr_list[[x]][["msLevel"]] == 1]
  })
  scan_ms1_idx <- lapply(1:length(filesPath), function(x){
    hdr_list[[x]][["msLevel"]] == 1
  })
  scan_ms1_List <- lapply(1:length(filesPath), function(x){
    pks_list[[x]][scan_ms1_idx[[x]]]
  })
  
  
  # Process the data of MS2 level
  scan_ms2_nonEmpty_idx <- lapply(1:length(filesPath), function(x){
    vapply(pks_list[[x]], nrow, integer(1L)) != 0
  })
  scanrts_ms2_List <- lapply(1:length(filesPath), function(x){
    hdr_list[[x]][["retentionTime"]][(hdr_list[[x]][["msLevel"]] == 2) & (scan_ms2_nonEmpty_idx[[x]])]
  })
  scan_ms2_idx <- lapply(1:length(filesPath), function(x){
    hdr_list[[x]][["msLevel"]] == 2 & (scan_ms2_nonEmpty_idx[[x]])
  })
  scan_ms2_List <- lapply(1:length(filesPath), function(x){
    pks_list[[x]][scan_ms2_idx[[x]]]
  })
  
  
  # Process the precursors OR swath information
  swath <- matrix()
  if(acquisitionMode == "DDA"){
    prec_mzs_List <- lapply(1:length(filesPath), function(x){
      mtx <- cbind(hdr_list[[x]][["precursorMZ"]][scan_ms2_idx[[x]]],
                   hdr_list[[x]][["precursorIntensity"]][scan_ms2_idx[[x]]])
    })
    checkPrecInt <- vapply(1:length(filesPath), function(x){
      all(prec_mzs_List[[x]][,2] == 0)
    }, FUN.VALUE = logical(1L))
    if(any(checkPrecInt)){
      noPreIntData_idx <- which(checkPrecInt)
      pseduInts <- lapply(noPreIntData_idx, function(x){
        vapply(scan_ms2_List[[x]], function(y) {max(y[,2])}, double(1L))
      })
      for(i in noPreIntData_idx){prec_mzs_List[[i]][,2] == pseduInts[[i]]}
    }
    
  } else if(acquisitionMode == "DIA"){
    swath <- read.table(SWATH_file, header = F)
  }

  # Extract the feature matrix from MS1 level
  if(acquisitionMode == "DDA"){
    if(!is.null(mSet)){
      feature_data <- mSet@peakAnnotation[["camera_output"]];
      mz_min <- feature_data[["mzmin"]]
      mz_max <- feature_data[["mzmax"]]
      rt_min <- feature_data[["rtmin"]]
      rt_max <- feature_data[["rtmax"]]
      peak_mtx <- cbind(mz_min, mz_max, rt_min, rt_max)
    } else if(!is.null(targetFeatures)) {
      peak_mtx <- targetFeatures
    } else {
      peak_mtx <- matrix();
    }
  } else {
    if(!is.null(mSet)){
      peak_mtx <- dia_feature_preparation(mSet@peakAnnotation[["groupmat"]],
                                          mSet@peakfilling[["msFeatureData"]][["chromPeaks"]],
                                          mSet@peakfilling[["FeatureGroupTable"]]@listData[["peakidx"]])
    } else if(!is.null(targetFeatures)) {
      peak_mtx <- lapply(1:nrow(targetFeatures), function(x){
        mz_min <- targetFeatures[x,1];
        mz_max <- targetFeatures[x,2];
        mz_med <- mean(c(mz_min, mz_max));
        rt_min <- targetFeatures[x,3];
        rt_max <- targetFeatures[x,4];
        rt_med <- mean(c(rt_min, rt_max));
        baseline_val <- 0.5;
        c(mz_med, mz_min, mz_max, rt_med, rt_min, rt_max, baseline_val)
      })
    } else {
      peak_mtx <- NA;
    }
  }
  

  # idxVec <- seq_along(filesPath);
  # if(acquisitionMode == "DDA"){
  #   for(i in idxVec){ # scanrts_ms1 scan_ms1
  #     idxBoolVec1 <- (scanidx_ms1 == i);
  #     idxBoolVec2 <- (scanidx_ms2 == i);
  #     scanrts_ms1_List[[i]] = scanrts_ms1[idxBoolVec1];
  #     scanrts_ms2_List[[i]] = scanrts_ms2[idxBoolVec2];
  #     scan_ms1_List[[i]] = scan_ms1[idxBoolVec1];
  #     scan_ms2_List[[i]] = scan_ms2[idxBoolVec2];
  #     prec_mzs_List[[i]] = prec_mzs[idxBoolVec2, ];
  #   }
  # } else if(acquisitionMode == "DIA"){
  #   for(i in idxVec){ # scanrts_ms1 scan_ms1
  #     idxBoolVec1 <- (scanidx_ms1 == i);
  #     idxBoolVec2 <- (scanidx_ms2 == i);
  #     scanrts_ms1_List[[i]] = scanrts_ms1[idxBoolVec1];
  #     scanrts_ms2_List[[i]] = scanrts_ms2[idxBoolVec2];
  #     scan_ms1_List[[i]] = scan_ms1[idxBoolVec1];
  #     scan_ms2_List[[i]] = scan_ms2[idxBoolVec2];
  #   }
  # }
  
  MSn_data <- list(peak_mtx = peak_mtx, 
                   precursors = prec_mzs_List,
                   scanrts_ms1 = scanrts_ms1_List,
                   scanrts_ms2 = scanrts_ms2_List,
                   scan_ms1 = scan_ms1_List,
                   scan_ms2 = scan_ms2_List,
                   swath = swath,
                   ion_mode = ion_mode,
                   acquisitionMode = acquisitionMode,
                   fileNames = basename(filesPath))
  
  #if(is.null(mSet)){
    mSet <- new("mSet")
    mSet@MSnData <- MSn_data
  #}
  
  return(mSet);
}



#' PerformDDADeconvolution
#'
#' @param mSet mSet
#' @param ppm1 ppm1
#' @param ppm2 ppm2
#' @param sn sn
#' @param filtering filtering
#' @param window_size window_size
#' @param intensity_thresh intensity_thresh
#' @param database_path database_path
#' @param ncores ncores
#' @param decoOn decoOn
#' @param useEntropy useEntropy
#'
#' @return mSet Object
#' @export
#'
#' @examples to add
PerformDDADeconvolution <- function(mSet= NULL, 
                                    ppm1 = 5, 
                                    ppm2 = 15,
                                    sn = 12,
                                    filtering = 2000,
                                    window_size = 1,
                                    intensity_thresh = 1e3,
                                    database_path = "",
                                    ncores = 1L,
                                    decoOn = TRUE,
                                    useEntropy = FALSE) {
  
  require(parallel);
  if(is.null(mSet)){
    stop("mSet is missing! Please import your MSn data first with function \"PerformMSnImport\" !")
  }
  if(ncores <= 0){
    warning("At least 1 cpu core required for processing!")
  }
  if(!is.integer(ncores)){
    stop("\"ncores\" must be an integer!")
  }
  if(database_path == ""){
    stop("Database is required, please define a valid database via \"database_path\".")
  }
  if(!file.exists(database_path)){
    stop("Your database doesnot exist in this path: ", database_path, ", please check!")
  }
  MSn_data <- mSet@MSnData;
  idxVec <- seq_along(MSn_data[["scanrts_ms2"]]);
  
  # find ion mode - ion_mode
  ion_mode <- MSn_data$ion_mode;
  
  #register(bpstart(SnowParam(2))); # for windows
  #register(bpstart(MulticoreParam(ncores)));
  
  # Orchestrating bpapply
  nPeaks <- nrow(MSn_data[["peak_mtx"]]);
  if(ncores<= nPeaks){
    Npbatch <- ceiling(nPeaks/ncores);
    batches <- split(c(1:nPeaks), ceiling(c(1:nPeaks)/Npbatch));
    if(length(batches) < ncores){
      batches <- split(c(1:nPeaks), rep(c(1:ncores),Npbatch)[1:nPeaks]);
    }
    col1 <- as.integer(vapply(idxVec, function(x) {rep(x, ncores)}, FUN.VALUE = integer(length = ncores)));
    col2 <- rep(1:ncores, length(idxVec));
  } else {
    batches <- list(seq(nPeaks))
    col1 <- idxVec
    col2 <- rep(1, length(idxVec))
  }
  
  ls_orches <- list(col1, col2, batches);
  
  
  if(ncores>1){
    cl <- makeCluster(getOption("cl.cores", ncores))
    clusterExport(cl, c("ls_orches", "MSn_data", 
                        "window_size", "ppm1", "ppm2", 
                        "sn", "filtering", "intensity_thresh", 
                        "ion_mode", "database_path", "decoOn", "useEntropy",
                        "PerformDDADeco"), envir = environment())
    DecRes <- list()
    DecRes <- parLapply(cl, 
                        seq(col1), 
                        function(x, ls_orches, winSize, ppm1, ppm2, sn, filt, 
                                 intensity_thresh, ionmode, db_path, decoOn, useEntropy){
                          x1 <- ls_orches[[1]][x];
                          x2 <- ls_orches[[3]][[ls_orches[[2]][x]]];
                          if(any(x2 == 1)){showOutput = TRUE} else{showOutput = FALSE}
                          if(length(x2) == 1){
                            peak_matrix <- matrix(MSn_data[["peak_mtx"]][x2,], 
                                                  ncol = length(MSn_data[["peak_mtx"]][x2,]));
                          } else {
                            peak_matrix <- MSn_data[["peak_mtx"]][x2,];
                          }
                          prec_mzs <- MSn_data[["precursors"]][[x1]];
                          scant1 <- MSn_data[["scanrts_ms1"]][[x1]];
                          scant2 <- MSn_data[["scanrts_ms2"]][[x1]];
                          scanms1 <- MSn_data[["scan_ms1"]][[x1]];
                          scanms2 <- MSn_data[["scan_ms2"]][[x1]];
                          fileName <- MSn_data[["fileNames"]][x1];
                          thread_int <- x;
                          res <- PerformDDADeco(peak_matrix,
                                                scant1, scant2,
                                                scanms1, scanms2,
                                                prec_mzs, winSize,
                                                ppm1, ppm2,
                                                sn, filt, intensity_thresh,
                                                ionmode, db_path,
                                                decoOn, useEntropy,
                                                showOutput, thread_int, fileName);
                          res},
                        ls_orches = ls_orches,
                        winSize = window_size,
                        ppm1 = ppm1,
                        ppm2 = ppm2,
                        sn = sn,
                        filt = filtering,
                        intensity_thresh = intensity_thresh,
                        ionmode = ion_mode,
                        db_path = database_path,
                        decoOn = decoOn,
                        useEntropy = useEntropy)
    
    stopCluster(cl)
    
  } else {
    
    DecRes <- list()
    DecRes <- lapply(seq(col1), 
                     function(x, ls_orches, winSize, ppm1, ppm2, sn, filt, 
                              intensity_thresh, ionmode, db_path, decoOn){
                       x1 <- ls_orches[[1]][x];
                       x2 <- ls_orches[[3]][[ls_orches[[2]][x]]];
                       if(any(x2 == 1)){showOutput = TRUE} else{showOutput = FALSE}
                       if(length(x2) == 1){
                         peak_matrix <- matrix(MSn_data[["peak_mtx"]][x2,], 
                                               ncol = length(MSn_data[["peak_mtx"]][x2,]));
                       } else {
                         peak_matrix <- MSn_data[["peak_mtx"]][x2,];
                       }
                       prec_mzs <- MSn_data[["precursors"]][[x1]];
                       scant1 <- MSn_data[["scanrts_ms1"]][[x1]];
                       scant2 <- MSn_data[["scanrts_ms2"]][[x1]];
                       scanms1 <- MSn_data[["scan_ms1"]][[x1]];
                       scanms2 <- MSn_data[["scan_ms2"]][[x1]];
                       fileName <- MSn_data[["fileNames"]][x1];
                       thread_int <- x;
                       res <- PerformDDADeco(peak_matrix,
                                             scant1, scant2,
                                             scanms1, scanms2,
                                             prec_mzs, winSize,
                                             ppm1, ppm2,
                                             sn, filt, intensity_thresh,
                                             ionmode, db_path,
                                             decoOn,
                                             showOutput, thread_int, fileName);
                       res},
                     ls_orches = ls_orches,
                     winSize = window_size,
                     ppm1 = ppm1,
                     ppm2 = ppm2,
                     sn = sn,
                     filt = filtering,
                     intensity_thresh = intensity_thresh,
                     ionmode = ion_mode,
                     db_path = database_path,
                     decoOn = decoOn)
  }

  
  # DecRes <- bplapply(seq(col1), 
  #                    FUN = function(x, ls_orches, winSize, ppm1, ppm2, sn, filt, intensity_thresh, ionmode, db_path, decoOn){
  #                      x1 <- ls_orches[[1]][x];
  #                      x2 <- ls_orches[[3]][[ls_orches[[2]][x]]];
  #                      if(any(x2 == 1)){showOutput = TRUE} else{showOutput = FALSE}
  #                      peak_matrix <- MSn_data[["peak_mtx"]][x2,];
  #                      prec_mzs <- MSn_data[["precursors"]][[x1]];
  #                      scant1 <- MSn_data[["scanrts_ms1"]][[x1]];
  #                      scant2 <- MSn_data[["scanrts_ms2"]][[x1]];
  #                      scanms1 <- MSn_data[["scan_ms1"]][[x1]];
  #                      scanms2 <- MSn_data[["scan_ms2"]][[x1]];
  #                      fileName <- MSn_data[["fileNames"]][x1];
  #                      thread_int <- x;
  #                      
  #                      res <- PerformDDADeco(peak_matrix, 
  #                                            scant1, scant2, 
  #                                            scanms1, scanms2, 
  #                                            prec_mzs, winSize, 
  #                                            ppm1, ppm2, 
  #                                            sn, filt, intensity_thresh, 
  #                                            ionmode, db_path, 
  #                                            decoOn, 
  #                                            showOutput, thread_int, fileName);
  #                      res},
  #                    ls_orches = ls_orches,
  #                    winSize = window_size,
  #                    ppm1 = ppm1,
  #                    ppm2 = ppm2,
  #                    sn = sn,
  #                    filt = filtering,
  #                    intensity_thresh = intensity_thresh,
  #                    ionmode = ion_mode,
  #                    db_path = database_path,
  #                    decoOn = decoOn,
  #                    BPPARAM = MulticoreParam(ncores));
  cat("\nDeconvolution executed successfully!")
  # Organize DecRes
  DecResx <- lapply(idxVec, function(x){
    idx0 <- c(ls_orches[[1]] == x);
    DecRes0 <- DecRes[which(idx0)];
    Spectra <- Indicator <- list(); 
    FeatureIdx <- vector();
    # spectra
    Spectra <- lapply(DecRes0, function(u){u[[1]]})
    Spectra <- do.call(c, Spectra)
    # Indicator
    Indicator <- lapply(DecRes0, function(u){u[[2]]})
    Indicator <- do.call(c, Indicator)
    # FeatureIdx
    FeatureIdx <- lapply(1:length(batches), function(u){
      d0 <- batches[[u]]
      d0[DecRes0[[u]][["FeatureIdx"]]+1]-1
    })
    FeatureIdx <- do.call(c, FeatureIdx)
    # return
    list(Spectra = Spectra, 
         Indicator = Indicator,
         FeatureIdx = FeatureIdx)
  })
  
  mSet@MSnResults$DecRes <- DecResx;
  return(mSet)
}



#' PerformDIADeconvolution
#'
#' @param mSet mSet Object contains raw spectral data from *PerformMSnImport*
#' @param min_width minimum peak width value, in seconds
#' @param ppm2 
#' @param sn 
#' @param span 
#' @param filtering 
#' @param ncores 
#'
#' @return mSet Object
#' @export
#'
#' @examples to add
PerformDIADeconvolution <- function(mSet= NULL, 
                                    min_width = 5,
                                    ppm2,
                                    sn = 12,
                                    span = 0.3,
                                    filtering = 2000,
                                    ncores = 1L) {
  
  # load("/ultraData/new_wb/HILICpos/mzML/MS1/mSet.rda")
  # load("data/MSn_data_dia.rda")
  require(BiocParallel);
  if(is.null(mSet)){
    stop("mSet is missing! Please import your MSn data first with function \"PerformMSnImport\" !")
  }
  if(missing(ppm2)){
    warning("ppm_ms2 is not defined, will use the default value 15 here!")
    ppm2 <- 15;
  }
  if(ncores <= 0){
    warning("At least 1 cpu core required for processing!")
  }
  if(!is.integer(ncores)){
    stop("\"ncores\" must be an integer!")
  }
  
  MSn_data <- mSet@MSnData;
  idxVec <- seq_along(MSn_data[["scanrts_ms2"]]);
  #register(bpstart(SnowParam(2))); # for windows
  #register(bpstart(MulticoreParam(ncores)));
  
  #Rcpp::sourceCpp("src/export_interfece.cpp")
  #x= 1; min_width = 2; ppm2 = 25; sn = 4;span = 0.3; filt = 100
  DecRes <- bplapply(idxVec, 
                     FUN = function(x, min_width, ppm2, sn, span, filt){
                       peak_list <- MSn_data[["peak_mtx"]];
                       swath <- as.matrix(MSn_data[["swath"]]);
                       scant1 <- MSn_data[["scanrts_ms1"]][[x]];
                       scant2 <- MSn_data[["scanrts_ms2"]][[x]];
                       scanms1 <- MSn_data[["scan_ms1"]][[x]];
                       scanms2 <- MSn_data[["scan_ms2"]][[x]];
                       
                       res <- PerformDIADeco(peak_list, 
                                             swath,
                                             scant1, scant2, 
                                             scanms1, scanms2, 
                                             min_width, 
                                             ppm2, 
                                             sn, 
                                             span,
                                             filt);
                       
                       
                       res},
                     min_width = min_width,
                     ppm2 = ppm2,
                     sn = sn,
                     span = span,
                     filt = filtering,
                     BPPARAM = MulticoreParam(ncores))
  register(bpstop());
  mSet@MSnResults$DecRes <- DecRes;
  return(mSet)
}


#' PerformSpectrumConsenus
#'
#' @param mSet 
#' @param ppm2 
#' @param concensus_fraction 
#' @param database_path 
#' @param use_rt 
#' @param user_dbCorrection 
#'
#' @return mSet Object
#' @export
#'
#' @examples to add
PerformSpectrumConsenus <- function(mSet = NULL, 
                                    ppm2,
                                    concensus_fraction = 0.5,
                                    database_path = "",
                                    use_rt = FALSE,
                                    user_dbCorrection = TRUE,
                                    useEntropy = FALSE) {
  
  
  if(is.null(mSet)){
    stop("mSet is missing! Please import your MSn data first with function \"PerformMSnImport\" and then run the data deconvolution, and then run this function!")
  }
  
  DecoedRes <- mSet@MSnResults$DecRes;
  MSn_data <- mSet@MSnData;
  ionmode <- MSn_data$ion_mode;
  
  if(is.null(DecoedRes)){
    stop("DecoedRes is missing! Please perform data deconvolution with \"PerformDIADeconvolution\" (for DIA data) or \"PerformDDADeconvolution\" (for DDA data)!")
  }

  if(missing(ppm2)){
    warning("ppm2 is not defined, will use the default value 15 here!")
    ppm2 <- 15;
  }

  # need to process/format peak_mtx here for DIA case
  if(MSn_data$acquisitionMode == "DIA"){
    peak_list <- MSn_data$peak_mtx;
    peak_mtx_res <- vapply(peak_list, function(x){
      c(x[2], x[3], x[5], x[6])
    }, FUN.VALUE = vector(mode = "double", length = 4))
    peak_mtx <- t(peak_mtx_res)
  } else {
    peak_mtx <- MSn_data$peak_mtx;
  }

  Concensus_spec = list();
  Concensus_spec = SpectrumConsensus(DecoedRes, 
                                     peak_mtx, 
                                     ppm2, 
                                     concensus_fraction, 
                                     user_dbCorrection, 
                                     database_path = database_path, 
                                     use_rt,
                                     ionmode, 
                                     useEntropy);
  
  ## filter out all empty spectrum
  idx_empty <- vapply(Concensus_spec[[2]], function(x){
    nrow(x[[1]]) == 0
  }, logical(length = 1L))
  Concensus_spec[[1]] <- Concensus_spec[[1]][!idx_empty]
  Concensus_spec[[2]] <- Concensus_spec[[2]][!idx_empty]
  mSet@MSnResults$Concensus_spec <- Concensus_spec;
  return(mSet)
}

#' PerformDBSearchingBatch
#'
#' @param mSet mSet Object contains raw spectral data after results consensus from *PerformSpectrumConsenus*
#' @param ppm1 numeric, ppm value of m/z for precursours;
#' @param ppm2 numeric, ppm value of m/z for ms/ms fragments matching;
#' @param rt_tol numeric, retention time tolerance, in seconds. Only effective when use_rt is TRUE;
#' @param database_path character, specify the path of database (.sqlite format);
#' @param use_rt logical, to use retention time if TRUE;
#' @param enableNL logical, to enable use Neutral Loss matching for unmatched features if TRUE;
#' @param NLdatabase_path path of neutral loss database. Must be specified to a valid neutral loss database when enableNL is TRUE.
#' @param ncores 
#'
#' @return mSet Object
#' @export
#'
PerformDBSearchingBatch <- function(mSet = NULL, 
                                    ppm1 = 5,
                                    ppm2 = 15,
                                    rt_tol = 0,
                                    database_path = "",
                                    use_rt = FALSE,
                                    enableNL = FALSE,
                                    NLdatabase_path = NULL,
                                    ncores = 1,
                                    useEntropy = FALSE
                                    ){
  
  # This function is a parallel function to enable the database searching in a parallel way
  require(parallel);
  if(is.null(mSet)){
    stop("mSet is required! Please run \"PerformSpectrumConsenus\" before this step!")
  }
  if(database_path == ""){
    stop("Database is required, please define a valid database via \"database_path\".")
  }
  if(!file.exists(database_path)){
    stop("Your database doesnot exist in this path: ", database_path, ", please check!")
  } else if(sub(pattern = "^(.*\\.|[^.]+)(?=[^.]*)", replacement = "", database_path, perl = TRUE) != "sqlite") {
    stop("Your database is not a valid sqlite file!")
  } else {
    mSet@MSnData[["database_path"]] <- database_path;
  }
  if(!is.numeric(ppm1) | !is.numeric(ppm2)){
    stop("ppm1 or ppm2 is not a numeric value!")
  }
  if(!is.numeric(rt_tol)){
    stop("rt_tol is not a numeric value!")
  }
  ConsensusRes <- mSet@MSnResults[["Concensus_spec"]];
  if(is.null(ConsensusRes)){
    stop("Missing ConsensusRes is not allowed. Please try to run \"PerformSpectrumConsenus\" first!")
  }
  ionmode <- mSet@MSnData[["ion_mode"]];
  if(ionmode != 0 & ionmode != 1){
    stop("Ion mode must be 0 (ESI-) or 1 (ESI+). Please don't modify the value in mSet@MSnData[[\"ion_mode\"]],");
  }
  if(ncores <= 0){
    warning("At least 1 cpu core required for processing!")
  }
  if(!is.integer(ncores)){
    stop("\"ncores\" must be an integer!")
  }
  if(is.null(NLdatabase_path) & enableNL){
    stop("'NLdatabase_path' must be specified because you are trying to use neutral loss database");
  } else if(is.null(NLdatabase_path)) {
    NLdatabase_path <- "";
  }
  if(enableNL){
    if(!file.exists(NLdatabase_path)) stop("Your neutral loss database does not exist!");
    if(sub(pattern = "^(.*\\.|[^.]+)(?=[^.]*)", replacement = "", NLdatabase_path, perl = TRUE) != "sqlite"){
      stop("Your neutral loss database is not a valid sqlite file!");
    } 
  }
  num_fts <- length(ConsensusRes[[1]]);
  seq_idx <- seq(0, num_fts-1);
  rt_ms1 <- mSet@MSnData[["scanrts_ms1"]];
  scan_ms1 <- mSet@MSnData[["scan_ms1"]];
  
  # need to process/format peak_mtx here for DIA case
  if(mSet@MSnData$acquisitionMode == "DIA"){
    peak_list <- mSet@MSnData[["peak_mtx"]];
    peak_mtx_res <- vapply(peak_list, function(x){
      c(x[2], x[3], x[5], x[6])
    }, FUN.VALUE = vector(mode = "double", length = 4))
    peak_matrix <- t(peak_mtx_res)
  } else {
    peak_matrix <- mSet@MSnData[["peak_mtx"]];
  }
  cat("==== Database searching against", basename(database_path), "started ====\n");
  if(ncores == 1){
    res <- SpectraSearching(ConsensusRes, 
                            seq_idx, 
                            peak_matrix,
                            ppm1, 
                            ppm2, 
                            rt_tol, 
                            rt_ms1, 
                            scan_ms1, 
                            ion_mode = ionmode, 
                            database_path = database_path, 
                            enableNL = enableNL,
                            NLdatabase_path = NLdatabase_path,
                            useEntropy = useEntropy)
  } else if (ncores > 1) {
    # perform parallel searching
    mem_num <- ceiling(num_fts/ncores); # number of members in each parallel
    ft_idx_grps <- lapply(1:ncores, function(x){
      rs <- seq_idx[(1 + (x-1)*mem_num):(x*mem_num)]
      rs[!is.na(rs)]
    })
    
    cl <- makeCluster(getOption("cl.cores", ncores))
    clusterExport(cl, c("ft_idx_grps", 
                        "ConsensusRes", "peak_matrix", 
                        "ppm1", "ppm2", "rt_tol", "rt_ms1", 
                        "scan_ms1", "ionmode", "database_path", 
                        "enableNL", "useEntropy", "NLdatabase_path",
                        "SpectraSearching"), envir = environment())
    res1 <- list()
    res1 <- parLapply(cl, 
                      1:ncores, 
                      function(x, ft_idx_grps, ConsensusRes, peak_matrix, 
                               ppm1, ppm2, rt_tol, rt_ms1, 
                               scan_ms1, ionmode, database_path, enableNL, 
                               NLdatabase_path, useEntropy){
                        res0 <- SpectraSearching(ConsensusRes, 
                                                 ft_idx_grps[[x]], 
                                                 peak_matrix,
                                                 ppm1, 
                                                 ppm2, 
                                                 rt_tol, 
                                                 rt_ms1, 
                                                 scan_ms1, 
                                                 ion_mode = ionmode, 
                                                 database_path = database_path, 
                                                 enableNL = enableNL,
                                                 NLdatabase_path = NLdatabase_path,
                                                 useEntropy = useEntropy)
                        res0},
                      ft_idx_grps = ft_idx_grps,
                      ConsensusRes = ConsensusRes,
                      peak_matrix = peak_matrix,
                      ppm1 = ppm1,
                      ppm2 = ppm2,
                      rt_tol = rt_tol,
                      rt_ms1 = rt_ms1,
                      scan_ms1 = scan_ms1,
                      ionmode = ionmode,
                      database_path = database_path,
                      enableNL = enableNL,
                      NLdatabase_path = NLdatabase_path,
                      useEntropy = useEntropy)
    stopCluster(cl)
    res <- list()
    for(i in 1:ncores){
      res <- c(res, res1[[i]])
    }
    
  } else {
    # perform parallel searching
    mem_num <- ceiling(num_fts/ncores); # number of members in each parallel
    ft_idx_grps <- lapply(1:ncores, function(x){
      rs <- seq_idx[(1 + (x-1)*mem_num):(x*mem_num)]
      rs[!is.na(rs)]
    })
    
    #register(bpstart(MulticoreParam(ncores)));
    #register(bpstart(SnowParam(6)));
    res1 <- bplapply(1:ncores, function(x, 
                                        ConsensusRes, 
                                        peak_matrix, 
                                        ppm1, 
                                        ppm2, 
                                        rt_tol, 
                                        rt_ms1, 
                                        scan_ms1,
                                        ionmode,
                                        database_path,
                                        enableNL){
      res0 <- SpectraSearching(ConsensusRes, 
                              ft_idx_grps[[x]], 
                              peak_matrix,
                              ppm1, 
                              ppm2, 
                              rt_tol, 
                              rt_ms1, 
                              scan_ms1, 
                              ion_mode = ionmode, 
                              database_path = database_path, 
                              enableNL = enableNL)
      res0},
      ConsensusRes = ConsensusRes,
      peak_matrix = peak_matrix,
      ppm1 = ppm1,
      ppm2 = ppm2,
      rt_tol = rt_tol,
      rt_ms1 = rt_ms1,
      scan_ms1 = scan_ms1,
      ionmode = ionmode,
      database_path = database_path,
      enableNL = enableNL,
      BPPARAM = MulticoreParam(ncores));
    #register(bpstop());
    res <- list()
    for(i in 1:ncores){
      res <- c(res, res1[[i]])
    }
  }
  cat("\n============ Searching finished successfully! ============");
  mSet@MSnResults[["DBMatchRes"]] <- res;
  return(mSet)
}

#' PerformResultsExport
#'
#' @param mSet 
#' @param type 
#' @param topN 
#' @param ncores 
#' @param lipids 
#'
#' @return mSet Object
#' @export
#'
#' @examples to add
PerformResultsExport <- function(mSet = NULL,
                                 type = 0L,
                                 topN = 10L,
                                 ncores = 1L,
                                 lipids = F){
  
  require(BiocParallel);
  if(is.null(mSet)){
    stop("mSet is required! Please run \"PerformSpectrumConsenus\" before this step!")
  }
  if(!(type %in% c(0L,1L,2L,3L))){
    stop("\"type\" must be one of 0L, 1L, 2L, 3L.")
  }
  if(!is.integer(topN)){
    stop("\"topN\" must be an integer!")
  }
  if(ncores <= 0){
    warning("At least 1 cpu core required for processing!")
  }
  if(!is.integer(ncores)){
    stop("\"ncores\" must be an integer!")
  }
  database_path <- mSet@MSnData[["database_path"]];
  if(database_path == ""){
    stop("Database is required, please define a valid database via \"database_path\".")
  }
  if(!file.exists(database_path)){
    stop("Your database is missing in this path: ", database_path, ", please check!")
  }
  
  Match_res <- mSet@MSnResults[["DBMatchRes"]];
  if(is.null(Match_res)){
    stop("Missing Matching result is not allowed. Please try to run \"PerformDBSearchingBatch\" first!")
  }
  ionmode <- mSet@MSnData[["ion_mode"]];
  if(ionmode != 0 & ionmode != 1){
    stop("Ion mode must be 0 (ESI-) or 1 (ESI+). Please don't modify the value in mSet@MSnData[[\"ion_mode\"]],");
  }
  cat("======= Converting started. See progress below ======= \n");
  if(ncores == 1){
    # single core
    annote_res <- annotation_export(Match_res, 
                                    type, 
                                    topN,
                                    ion_mode = ionmode,
                                    database_path = database_path,
                                    lipidsClass = lipids)
  } else {
    # multiple cores, need to split Match_res first
    num_res <- ceiling(length(Match_res)/ncores);
    seq_idx <- seq(1, length(Match_res));
    ft_idx_grps <- lapply(1:ncores, function(x){
      rs <- seq_idx[(1 + (x-1)*num_res):(x*num_res)]
      rs[!is.na(rs)]
    })
    #register(bpstart(MulticoreParam(ncores)));
    #register(bpstart(SnowParam(6)));
    if(.Platform$OS.type == "windows"){
      bppm <- SnowParam(ncores, progressbar = TRUE);
    } else {
      bppm <- MulticoreParam(ncores);
    }
    res1 <- bplapply(1:ncores, function(x, Match_res, type, topN, ionmode, database_path, lipidsClass){
      res0 <- annotation_export(Match_res[ft_idx_grps[[x]]], 
                                type, topN,
                                ion_mode = ionmode,
                                database_path = database_path, lipidsClass = lipidsClass);
    }, 
    Match_res = Match_res,
    type = type,
    topN = topN, 
    ionmode = ionmode,
    database_path = database_path, 
    lipidsClass = lipids,
    BPPARAM = bppm)
    
    #register(bpstop());
    annote_res <- list()
    for(i in 1:ncores){
      annote_res <- c(annote_res, res1[[i]])
    }
  }
  cat("\n============ Exporting finished successfully! ============");
  mSet@MSnResults[["DBAnnoteRes"]] <- annote_res;
  return(mSet)
}


#' Title
#'
#' @param mSet 
#' @param topN 
#' @param isLipidomics 
#'
#' @return mSet Object
#' @export
#'
#' @examples to add
FormatMSnAnnotation <- function(mSet = NULL,
                                topN = 5L,
                                isLipidomics = FALSE){
  
  if(is.null(mSet)){
    stop("mSet is required! Please run \"PerformResultsExport\" before this step!")
  }
  if(!is.integer(topN)){
    stop("\"topN\" must be an integer!")
  }
  if(is.null(mSet@MSnResults[["DBAnnoteRes"]])){
    stop("Annotation need to be converted at first. Please run \"PerformResultsExport\" before this step!")
  }
  if(length(mSet@MSnResults[["DBAnnoteRes"]]) == 0){
    warning("No annotated peaks exported. Will return nothing.")
    return(NULL)
  }
  ## the mSet is valid, start processing
  ## formatting MS1 peaks
  if(mSet@MSnData[["acquisitionMode"]] == "DIA"){
    peak_list <- mSet@MSnData[["peak_mtx"]]
    peak_mtx <- do.call(rbind, peak_list)[,c(2,3,5,6)]
    colnames(peak_mtx) <- c("mzmin", "mzmax", "rtmin", "rtmax")
  } else {
    peak_mtx <- mSet@MSnData[["peak_mtx"]]
    colnames(peak_mtx) <- c("mzmin", "mzmax", "rtmin", "rtmax")
  }
  peak_idx <- mSet@MSnResults[["Concensus_spec"]][[1]]
  peak_mtx_identified <- peak_mtx[peak_idx+1,]
  ## formarring MS2 results [topN, compounds, inchikey and scores]
  maxN <- max(vapply(1:length(mSet@MSnResults$DBAnnoteRes), function(x){
    length(mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][[1]])
  }, integer(1L)))
  if(topN > maxN){
    warning("Current MSn annotation results contain at most ", maxN, " annotations.")
    topN <- maxN;
  }
  nmss <- names(mSet@MSnResults[["DBAnnoteRes"]][[1]][[1]])
  if(all(c("InchiKeys", "Compounds", "Formula") %in% nmss)){
    topNResults <- lapply(1:length(mSet@MSnResults$DBAnnoteRes), function(x){
      
      cmpds <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Compounds"]];
      inchikeys <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["InchiKeys"]];
      formulas <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Formula"]];
      scores <-  mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Scores"]];
      databases <-  mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Database"]]
      
      uincks <- unique(inchikeys)
      uincks <- uincks[uincks != ""]
      if(length(uincks) > topN){
        uincks <- uincks[1:topN]; #original order kept
      }
      idx <- vapply(uincks, function(x){which(x==inchikeys)[1]}, FUN.VALUE = integer(1L),USE.NAMES = F)
      if(isLipidomics){
        if(!is.null(mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Super_Class"]])){
          supcls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Super_Class"]];
          mancls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Main_Class"]];
          subcls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Sub_Class"]];
          data.frame(cmpds = cmpds[idx],
                     inchikeys = inchikeys[idx],
                     formulas = formulas[idx],
                     scores = scores[idx],
                     databases = databases[idx],
                     supclass = supcls[idx],
                     mainclass = mancls[idx],
                     subclass = subcls[idx])
        } else {
          data.frame(cmpds = cmpds[idx],
                     inchikeys = inchikeys[idx],
                     formulas = formulas[idx],
                     scores = scores[idx],
                     databases = databases[idx],
                     supclass = NA,
                     mainclass = NA,
                     subclass = NA)
        }
      } else {
        data.frame(cmpds = cmpds[idx],
                   inchikeys = inchikeys[idx],
                   formulas = formulas[idx],
                   scores = scores[idx],
                   databases = databases[idx])
      }
    })
    
  } else if("InchiKeys" %in% nmss){
    topNResults <- lapply(1:length(mSet@MSnResults$DBAnnoteRes), function(x){
      
      inchikeys <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["InchiKeys"]];
      scores <-  mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Scores"]];
      
      uincks <- unique(inchikeys)
      uincks <- uincks[uincks != ""]
      if(length(uincks) > topN){
        uincks <- uincks[1:topN]; #original order kept
      }
      idx <- vapply(uincks, function(x){which(x==inchikeys)[1]}, FUN.VALUE = integer(1L),USE.NAMES = F)
      if(isLipidomics){
        if(!is.null(mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Super_Class"]])){
          supcls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Super_Class"]];
          mancls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Main_Class"]];
          subcls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Sub_Class"]];
          data.frame(inchikeys = inchikeys[idx],
                     scores = scores[idx],
                     supclass = supcls[idx],
                     mainclass = mancls[idx],
                     subclass = subcls[idx])
        } else {
          data.frame(inchikeys = inchikeys[idx],
                     scores = scores[idx],
                     supclass = NA,
                     mainclass = NA,
                     subclass = NA)
        }
      } else {
        data.frame(inchikeys = inchikeys[idx],
                   scores = scores[idx])
      }
    })
    
  } else if("Compounds" %in% nmss){
    topNResults <- lapply(1:length(mSet@MSnResults$DBAnnoteRes), function(x){
      
      cmpds <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Compounds"]];
      scores <-  mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Scores"]];
      
      uincks <- unique(cmpds)
      uincks <- uincks[uincks != ""]
      if(length(uincks) > topN){
        uincks <- uincks[1:topN]; #original order kept
      }
      idx <- vapply(uincks, function(x){which(x==cmpds)[1]}, FUN.VALUE = integer(1L),USE.NAMES = F)
      if(isLipidomics){
        if(!is.null(mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Super_Class"]])){
          supcls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Super_Class"]];
          mancls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Main_Class"]];
          subcls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Sub_Class"]];
          data.frame(cmpds = cmpds[idx],
                     scores = scores[idx],
                     supclass = supcls[idx],
                     mainclass = mancls[idx],
                     subclass = subcls[idx])
        } else {
          data.frame(cmpds = cmpds[idx],
                     scores = scores[idx],
                     supclass = NA,
                     mainclass = NA,
                     subclass = NA)
        }
      } else {
        data.frame(cmpds = cmpds[idx],
                   scores = scores[idx])
      }
    })
    
  }  else if("Formula" %in% nmss){
    topNResults <- lapply(1:length(mSet@MSnResults$DBAnnoteRes), function(x){
      
      formulas <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Formula"]];
      scores <-  mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Scores"]];
      
      uincks <- unique(formulas)
      uincks <- uincks[uincks != ""]
      if(length(uincks) > topN){
        uincks <- uincks[1:topN]; #original order kept
      }
      idx <- vapply(uincks, function(x){which(x==formulas)[1]}, FUN.VALUE = integer(1L),USE.NAMES = F)
      if(isLipidomics){
        if(!is.null(mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Super_Class"]])){
          supcls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Super_Class"]];
          mancls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Main_Class"]];
          subcls <- mSet@MSnResults[["DBAnnoteRes"]][[x]][[1]][["Sub_Class"]];
          data.frame(formulas = formulas[idx],
                     scores = scores[idx],
                     supclass = supcls[idx],
                     mainclass = mancls[idx],
                     subclass = subcls[idx])
        } else {
          data.frame(formulas = formulas[idx],
                     scores = scores[idx],
                     supclass = NA,
                     mainclass = NA,
                     subclass = NA)
        }
      } else {
        data.frame(formulas = formulas[idx],
                   scores = scores[idx])
      }
    })
    
  }

  keep_peak_idx <- vapply(topNResults, function(x){nrow(x) != 0}, logical(length = 1L))
  peak_mtx_identified <- peak_mtx_identified[keep_peak_idx, ]
  ## formatting lipidomics results
  if(isLipidomics){
    # generate colnames
    colnms <- vapply(1:topN, function(x){
      paste0(c("Compound_", "InchiKey_", "Formula_", "Score_", "Database_", "SuperClass_", "MainClass_", "SubClass_"), x)
    }, FUN.VALUE = character(8L))
    dim(colnms) <- c(1,topN*8);
    
    res <- lapply(topNResults, function(x){
      if(nrow(x) == 0){
        df_res <- data.frame(matrix(ncol = 8*topN, nrow = 0))
      } else {
        df_res <- data.frame(matrix(ncol = 8*topN, nrow = 1))
        if(!is.null(x$cmpds)){
          df_res[1,seq(1,(nrow(x)-1)*8+1, by = 8)] <- x$cmpds
        }
        if(!is.null(x$inchikeys)){
          df_res[1,seq(2,(nrow(x)-1)*8+2, by = 8)] <- x$inchikeys
        }
        if(!is.null(x$formulas)){
          df_res[1,seq(3,(nrow(x)-1)*8+3, by = 8)] <- x$formulas
        }
        df_res[1,seq(4,(nrow(x)-1)*8+4, by = 8)] <- ifelse(x$scores<0.1, yes = 0, no = round(x$scores,2))
        df_res[1,seq(5,(nrow(x)-1)*8+5, by = 8)] <- x$databases
        
        df_res[1,seq(6,(nrow(x)-1)*8+6, by = 8)] <- x$supclass
        df_res[1,seq(7,(nrow(x)-1)*8+7, by = 8)] <- x$mainclass
        df_res[1,seq(8,(nrow(x)-1)*8+8, by = 8)] <- x$subclass
      }
      colnames(df_res) <- colnms
      df_res
    })
    
  } else {
    # no need to extract class information
    # generate colnames
    colnms <- vapply(1:topN, function(x){
      paste0(c("Compound_", "InchiKey_", "Formula_", "Score_", "Database_"), x)
    }, FUN.VALUE = character(5L))
    dim(colnms) <- c(1,topN*5);
    
    res <- lapply(topNResults, function(x){
      if(nrow(x) == 0){
        df_res <- data.frame(matrix(ncol = 5*topN, nrow = 0))
      } else {
        df_res <- data.frame(matrix(ncol = 5*topN, nrow = 1))
        if(!is.null(x$cmpds)){
          df_res[1,seq(1,(nrow(x)-1)*5+1, by = 5)] <- x$cmpds
        }
        if(!is.null(x$inchikeys)){
          df_res[1,seq(2,(nrow(x)-1)*5+2, by = 5)] <- x$inchikeys
        }
        if(!is.null(x$formulas)){
          df_res[1,seq(3,(nrow(x)-1)*5+3, by = 5)] <- x$formulas
        }
        df_res[1,seq(4,(nrow(x)-1)*5+4, by = 5)] <- ifelse(x$scores<0.1, yes = 0, no = round(x$scores,2))
        df_res[1,seq(5,(nrow(x)-1)*5+5, by = 5)] <- x$databases
      }
      colnames(df_res) <- colnms
      df_res
    })
  }
  ## generating results table
  res_dt <- do.call(rbind, res)
  if(all(!is.na(res_dt$InchiKey_1))){
    uniInks <- unique(res_dt$InchiKey_1);
    row_idx <- vapply(uniInks, function(x){
      score_nums <- res_dt$Score_1[res_dt$InchiKey_1 == x]
      which(res_dt$InchiKey_1 == x)[which.max(score_nums)]
    }, integer(1L))
    peak_mtx_res <- peak_mtx_identified[row_idx,]
    res_dtx <- res_dt[row_idx,]
  } else if(all(!is.na(res_dt$Compound_1))){
    uniInks <- unique(tolower(res_dt$Compound_1));
    row_idx <- vapply(uniInks, function(x){
      score_nums <- res_dt$Score_1[tolower(res_dt$Compound_1) == x]
      which(res_dt$Compound_1 == x)[which.max(score_nums)]
    }, integer(1L))
    peak_mtx_res <- peak_mtx_identified[row_idx,]
    res_dtx <- res_dt[row_idx,]
  } else {
    peak_mtx_res <- peak_mtx_identified;
    res_dtx <- res_dt;
  }
  
  peak_mtx_res[,1] <- round(peak_mtx_res[,1], 4)
  peak_mtx_res[,2] <- round(peak_mtx_res[,2], 4)
  peak_mtx_res[,3] <- round(peak_mtx_res[,3], 2)
  peak_mtx_res[,4] <- round(peak_mtx_res[,4], 2)
  
  res_final <- cbind(peak_mtx_res, res_dtx)
  
  return(res_final)
}


# Function to parse the MS2Peaks string into a two-column matrix
parse_ms2peaks <- function(ms2peaks_str) {
  # Split the string into lines
  lines <- strsplit(ms2peaks_str, "\n")[[1]]
  
  # Split each line into mz and intensity, and convert to numeric
  # need to determine mz and intensity values are separated by " " or "\t", because there are two types in the sqlite database even in one entry
  # Process each line individually to handle different delimiters
  peaks <- lapply(lines, function(line) {
    # Determine delimiter by checking for tabs or spaces
    delimiter <- ifelse(grepl("\t", line), "\t", " ")
    # Split the line using the determined delimiter and convert to numeric
    numeric_line <- as.numeric(strsplit(trimws(line), delimiter)[[1]])
    return(numeric_line)
  })
  
  # Combine all the peaks into a matrix, handling potential uneven lengths by padding with NA
  max_length <- max(sapply(peaks, length))
  peaks_matrix <- do.call(rbind, lapply(peaks, function(x) {
    c(x, rep(NA, max_length - length(x)))
  }))
  
  # Remove rows that contain NAs (which may have been due to empty or incorrect format lines)
  peaks_matrix <- peaks_matrix[!apply(is.na(peaks_matrix), 1, any), ]
  
  # Ensure the result is a matrix with two columns
  if (length(peaks_matrix) == 2) {
    peaks_matrix <- matrix(peaks_matrix, nrow = 1)  # Convert single row vector to a matrix
  }
  return(peaks_matrix)
}




## normalize the spectrum
normalize_spectrum <- function(spec, cutoff_relative){
  ##cutoff_relative: relative noise cutoff
  tmp <- data.frame(mz = spec[, 1], intensity = spec[, 2])
  tmp$normalized <- round((tmp$intensity/max(tmp$intensity)) * 100)
  subset(tmp, tmp$normalized >= cutoff_relative)
}

# Calculate the ppm difference
ppm_diff <- function(mz1, mz2, ppm) {
  abs(mz1 - mz2) / mz2 * 1e6 <= ppm
}

# Find matching m/z values within the ppm tolerance
find_matches <- function(top_mz, bottom_mz, ppm) {
  matched <- vector("list", length(top_mz))
  for (i in seq_along(top_mz)) {
    matched[[i]] <- which(ppm_diff(top_mz[i], bottom_mz, ppm))
  }
  matched
}


MirrorPlotting <- function(spec.top, 
                           spec.bottom, 
                           ppm = 25, 
                           title= NULL, 
                           subtitle = NULL,
                           cutoff_relative = 5) {
  
  if (is.null(spec.top) || is.null(spec.bottom)) {
    stop("spec.top and spec.bottom cannot be NULL.")
  }
  
  top <- normalize_spectrum(spec.top, cutoff_relative)
  bottom <- normalize_spectrum(spec.bottom, cutoff_relative)
  
  ### Calculate the ppm difference
  ##ppm_diff <- function(mz1, mz2) {
  ##  abs(mz1 - mz2) / mz2 * 1e6
  ##}
  ##
  ### Find matching m/z values within the ppm tolerance
  ##matches <- which(sapply(top$mz, function(mz_top) {
  ##  sapply(bottom$mz, function(mz_bottom) {
  ##    ppm_diff(mz_top, mz_bottom) <= ppm
  ##  }) %>% any
  ##}))
  
  matches <- find_matches(top$mz, bottom$mz, ppm)
  
  # Add an additional column for plotting stars
  library(dplyr)
  top$star_intensity <- ifelse(matches %in% 1:nrow(bottom), top$normalized + 5, NA)
  #top$star_intensity <- ifelse(1:nrow(top) %in% matches, top$normalized + 5, NA)

  # Plotting with ggplot2
  library(ggplot2)
  
  ## determine the x axis range
  range_top <- range(top$mz, na.rm = TRUE)
  range_bottom <- range(bottom$mz, na.rm = TRUE)
  # The actual limits will be the minimum of the combined ranges and the maximum of the combined ranges
  x_axis_limits <- range(c(range_top, range_bottom), na.rm = TRUE)
  ## modify the range if it's too narrow
  if (x_axis_limits[2] - x_axis_limits[1] <= 50) {
    mean_x_axis_limits <- mean(x_axis_limits)
    x_axis_limits <- c(mean_x_axis_limits-25, mean_x_axis_limits+25)
  }
  
  # Check if there are any non-NA values in star_intensity
  #could add a layer for structure 
  #inchikey <- "IKGXIBQEEMLURG-NVPNHPEKSA-N"
  #compound <- ChemmineR::pubchemInchikey2sdf(inchikey)
  #plot(compound$sdf_set)
  if (any(!is.na(top$star_intensity))) {
    p <- ggplot() +
      geom_segment(data=top, aes(x = mz, xend = mz, y = 0, yend = normalized), color = "blue") +
      geom_segment(data=bottom, aes(x = mz, xend = mz, y = 0, yend = -normalized), color = "red") +
      geom_point(data=top, aes(x = mz, y = star_intensity), shape=18, color="red", size=3, na.rm = TRUE) +
      labs(title = title, subtitle = subtitle, y = "Relative Intensity (%)", x = "m/z") +
      theme_minimal() +
      xlim(x_axis_limits) +
      ylim(-105, 105)
  }else{
    p <- ggplot() +
      geom_segment(data=top, aes(x = mz, xend = mz, y = 0, yend = normalized), color = "blue") +
      geom_segment(data=bottom, aes(x = mz, xend = mz, y = 0, yend = -normalized), color = "red") +
      labs(title = title, subtitle = subtitle, y = "Relative Intensity (%)", x = "m/z") +
      theme_minimal() +
      xlim(x_axis_limits) +
      ylim(-105, 105)
  }
  
  return(p)
}




## to be done
##mark some mz values
##cutoff_relative or cutoff_absolute
#' PerformMirrorPlotting
#'
#' @param mSet mSet
#' @param cutoff_relative cutoff value for relative intensity
#' @param ppm ppm of mz of msms fragment
#' @param display_plot display mirror plot or not, TRUE or FALSE, default is FALSE
#' @param width width of the image, default is 8
#' @param height height of the image, default is 6
#' @param dpi dpi value, default is 300
#' @param interactive to make figure interactive, default is FALSE. TRUE or FALSE.
#' @export
#'
PerformMirrorPlotting <- function(mSet = NULL, 
                                  cutoff_relative = 5, 
                                  ppm = 25,
                                  display_plot = FALSE,
                                  width = 8, 
                                  height = 6, 
                                  dpi = 300,
                                  format = "png", 
                                  interactive = FALSE) {
  # Retrieve the DBAnnoteRes list
  dbAnnoteResList <- mSet@MSnResults[["DBAnnoteRes"]]
  
  if (is.null(dbAnnoteResList)) {
    stop("Process mSet with PerformResultsExport function (if you want to plot mirror plot, set type as 0L).")
  }
  
  # determine if all dbAnnoteResList is null
  # Initialize a flag to track if any non-NULL is found
  allNull <- TRUE
  
  # Loop through the list to check for non-NULL content
  for (i in seq_along(dbAnnoteResList)) {
    # Access the first element and then the MS2Pekas of each sublist
    ms2Pekas <- dbAnnoteResList[[i]][[1]][["MS2Pekas"]]
    
    # Check if MS2Pekas is not NULL
    if (!is.null(ms2Pekas)) {
      allNull <- FALSE
      #cat(sprintf("Non-NULL content found in sublist %d, stopping the check.\n", i))
      break # Stop the loop if non-NULL is found
    }
  }
  
  # If allNull is TRUE, it means all sublists had NULL in MS2Pekas
  if (allNull) {
    stop("Process mSet with PerformResultsExport function with the type set as 0L, then run PerformMirrorPlotting).")
  }
  
  # Iterate over each 'spec.top'
  spec.top.list <- mSet@MSnResults[["Concensus_spec"]][[2]]
  
  for (i in seq_along(spec.top.list)) {
    spec.top.m <- spec.top.list[[i]][[1]]
    ##spec.top.m <- parse_ms2peaks(spec.top)
    
    # Check if there are corresponding 'spec.bottom' spectra
    if (i <= length(dbAnnoteResList) && !is.null(dbAnnoteResList[[i]]) && length(dbAnnoteResList[[i]][[1]][["MS2Pekas"]]) != 0) {
      
      # Extract spec.bottom spectra
      spec.bottom.list <- dbAnnoteResList[[i]][[1]][["MS2Pekas"]]
      
      ## for title and file name
      if (mSet@MSnData[["acquisitionMode"]] == "DIA") {
        mz <- mSet@MSnData[["peak_mtx"]][[mSet@MSnResults[["Concensus_spec"]][[1]][[i]]+1]][[1]]
        rt <- mSet@MSnData[["peak_mtx"]][[mSet@MSnResults[["Concensus_spec"]][[1]][[i]]+1]][[4]]
      }
      if (mSet@MSnData[["acquisitionMode"]] == "DDA") {
        mz = mean(mSet@MSnData[["peak_mtx"]][mSet@MSnResults[["Concensus_spec"]][[1]][[i]]+1,1], mSet@MSnData[["peak_mtx"]][mSet@MSnResults[["Concensus_spec"]][[1]][[i]]+1,2])
        rt = mean(mSet@MSnData[["peak_mtx"]][mSet@MSnResults[["Concensus_spec"]][[1]][[i]]+1,3], mSet@MSnData[["peak_mtx"]][mSet@MSnResults[["Concensus_spec"]][[1]][[i]]+1,4])
      }
      

      # Loop through each spec.bottom
      for (j in seq_along(spec.bottom.list)) {
        spec.bottom <- spec.bottom.list[[j]]
        
        # Parse strings to numeric matrix
        spec.bottom.m <- parse_ms2peaks(spec.bottom)
        
        # Perform mirror plotting
        compound_name <- dbAnnoteResList[[i]][[1]][["Compounds"]][[j]]
        score <- round(dbAnnoteResList[[i]][[1]][["Scores"]][[j]], 2)
        database <- dbAnnoteResList[[i]][[1]][["Database"]][[j]]
        title <- paste0(mz, "__", rt)
        subtitle <- paste0(compound_name, "\n", score, "\n", database)
        p1 <- MirrorPlotting(spec.top.m, 
                             spec.bottom.m, 
                             ppm = ppm,
                             title= title, 
                             subtitle = subtitle,
                             cutoff_relative = cutoff_relative)
        
        # Construct the plot file name
        plot_filename <- paste0(sub_dir, "/", as.character(paste0(mz, "__", rt)), "_", j, ".png")
        
        # Save the plot with Cairo
        # ggsave(plot_filename, plot = last_plot(), width = width, height = height, dpi = dpi)
        if(display_plot){
          print(p1)
        } else {
          Cairo::Cairo(
            file = paste0("mirrorplot_", mz, "_", rt, "_", dpi, "_", j, ".", format),
            unit = "in", dpi = dpi, width = width, height = height, type = format, bg = "white")
          print(p1)
          dev.off()
        }
      }
    }
  }
}
