#' OptiLCMS: A package for computating the notorious bar statistic.
#'
#' The OptiLCMS package provides a pipeline for metabolomics processing.
#' 
#' @section OptiLCMS functions:
#' The OptiLCMS functions ...
#'
#' @docType package
#' @name OptiLCMS
#' @useDynLib OptiLCMS, .registration=TRUE
NULL
#> NULL



######## =-------- Call C by .C -------= ########

#' continuousPtsAboveThresholdR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
continuousPtsAboveThresholdR <-
  function(y, threshold, num, istart = 1) {
    if (!is.double(y)){
      y <- as.double(y)
    }

    n <- continuousPtsAboveThreshold(y, as.integer(istart - 1), 
                                     length(y), 
                                     threshold = as.double(threshold), 
                                     num = as.integer(num));
    if (n > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
}

#' continuousPtsAboveThresholdIdxR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
continuousPtsAboveThresholdIdxR <-
  function(y, threshold, num, istart = 1) {
    if (!is.double(y)){
      y <- as.double(y)
    }
    
    as.logical(continuousPtsAboveThresholdIdx(y, as.integer(istart - 1), length(y), threshold = as.double(threshold), num = as.integer(num)))
}

#' descendMinR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
descendMinR <- function(y, istart = which.max(y)) {
  if (!is.double(y))
    y <- as.double(y)
  DescendMin(y, length(y), as.integer(istart - 1))+1
}

#' which.colMax
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
which.colMax <- function (x, na.rm = FALSE, dims = 1) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.array(x) || length(dn <- dim(x)) < 2)
    stop("`x' must be an array of at least two dimensions")
  if (dims < 1 || dims > length(dn) - 1)
    stop("invalid `dims'")
  n <- prod(dn[seq_len(dims)])
  dn <- dn[-(seq_len(dims))]
  if (!is.double(x))
    x <- as.double(x)

  z <- WhichColMax(x, as.integer(n), as.integer(prod(dn)))
  if (length(dn) > 1) {
    dim(z) <- dn
    dimnames(z) <- dimnames(x)[-(seq_len(dims))]
  }
  else
    names(z) <- dimnames(x)[[dims + 1]]
  z
}


#' DescendZeroR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
DescendZeroR <- function(y, istart = which.max(y)){
  if (!is.double(y)) y <- as.double(y)
  DescendZero( y, length(y), as.integer(istart - 1))+1
}

#' colMaxR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
colMaxR <- function (x, na.rm = FALSE, dims = 1) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.array(x) || length(dn <- dim(x)) < 2)
    stop("`x' must be an array of at least two dimensions")
  if (dims < 1 || dims > length(dn) - 1)
    stop("invalid `dims'")
  n <- prod(dn[seq_len(dims)])
  dn <- dn[-(seq_len(dims))]
  if (!is.double(x))
    x <- as.double(x)

  z <- ColMax(x, as.integer(n), as.integer(prod(dn)))
  
  if (length(dn) > 1) {
    dim(z) <- dn
    dimnames(z) <- dimnames(x)[-(seq_len(dims))]
  }
  else
    names(z) <- dimnames(x)[[dims + 1]]
  z
}

#' rectUniqueR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
rectUniqueR <-
  function(m,
           order = seq(length = nrow(m)),
           xdiff = 0,
           ydiff = 0) {
    nr <- nrow(m)
    nc <- ncol(m)
    if (!is.double(m))
      m <- as.double(m)

    RectUnique(m, as.integer(order - 1), nr, nc, as.double(xdiff), as.double(ydiff)) == 1
  }

#' findEqualGreaterMR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
findEqualGreaterMR <- function(x, values) {
  if (!is.double(x))
    x <- as.double(x)
  if (!is.double(values))
    values <- as.double(values)
  FindEqualGreaterM(x, length(x), values, length(values)) + 1;
}


######## = ------- Call C by .Call -----= #########
#' R_set_obiwarpR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
R_set_obiwarpR <-
  function(valscantime1,
           scantime1,
           mzvals,
           mzs,
           cntrPr,
           valscantime2,
           scantime2,
           curP,
           parms) {

      rtadj <-
        .Call(
          R_set_obiwarp,
          valscantime1,
          scantime1,
          mzvals,
          mzs,
          cntrPr$profMat,
          valscantime2,
          scantime2,
          mzvals,
          mzs,
          curP$profMat,
          parms$response,
          parms$distFun,
          parms$gapInit,
          parms$gapExtend,
          parms$factorDiag,
          parms$factorGap,
          as.numeric(parms$localAlignment),
          parms$initPenalty,
          PACKAGE = 'OptiLCMS'
        )
 
    return (rtadj)
  }

#' getEIC4Peaks
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
getEIC4Peaks <-
  function(mset, peaks, maxscans) {
    mset$env$mz <-
      lapply(MSnbase::spectra(mset$onDiskData, BPPARAM = SerialParam()),
             mz)
    mset$env$intensity <-
      MSnbase::intensity(mset$onDiskData, BPPARAM = SerialParam())
    
    mset$scantime <- MSnbase::rtime(mset$onDiskData)
    
    valCount <- cumsum(lengths(mset$env$mz, FALSE))
    ## Get index vector for C calls
    mset$scanindex <-
      as.integer(c(0, valCount[-length(valCount)])) 
    
    npeaks <- dim(peaks)[1]
    
    scans  <- length(mset$scantime)
    eics <- matrix(NA, npeaks, maxscans)
    
    mset$env$mz <- as.double(unlist(mset$env$mz))
    mset$env$intensity <- as.double(unlist(mset$env$intensity))
    scan_time_leng <- as.integer(length(mset$scantime))
    
    MessageOutput(paste0(npeaks, " peaks found!"), "\n", NULL)
    
    for (p in seq_len(npeaks)) {
      timerange <- c(peaks[p, "rtmin"], peaks[p, "rtmax"])
      tidx <-
        which((mset$scantime >= timerange[1]) &
                (mset$scantime <= timerange[2]))
      
      if (length(tidx) > 0) {
        scanrange <- range(tidx)
      } else{
        scanrange <- seq_len(scans)
      }
      massrange <- c(peaks[p, "mzmin"], peaks[p, "mzmax"])

        eic <- .Call(
          getEIC,
          mset$env$mz,
          mset$env$intensity,
          as.integer(mset$scanindex),
          as.double(massrange),
          as.integer(scanrange),
          scan_time_leng,
          PACKAGE = 'OptiLCMS'
        )$intensity
        
      eic[eic == 0] <- NA;
      eics[p, scanrange[1]:scanrange[2]] <- eic;
    }
    
    eics
  }

#' breaks_on_nBinsR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
breaks_on_nBinsR <-
  function(fromX, toX, nBins, shiftByHalfBinSize = FALSE) {
    if (missing(fromX) | missing(toX) | missing(nBins))
      stop("'fromX', 'toX' and 'nBins' are required!")
    if (!is.double(fromX))
      fromX <- as.double(fromX)
    if (!is.double(toX))
      toX <- as.double(toX)
    if (!is.integer(nBins))
      nBins <- as.integer(nBins)
    shiftIt <- 0L
    
    if (shiftByHalfBinSize) {
      shiftIt <- 1L
    }
    
      res <- .Call(breaks_on_nBins, fromX, toX, nBins, shiftIt,
                   PACKAGE = "OptiLCMS")
    
    return(res)
  }

#' imputeLinInterpol
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
imputeLinInterpol <-
  function(x,
           baseValue,
           method = "lin",
           distance = 1L,
           noInterpolAtEnds = FALSE) {
    method <-
      match.arg(method, c("none", "lin", "linbase")) ## interDist? distance = 2
    if (method == "none") {
      return(x)
    }
    if (!is.double(x))
      x <- as.double(x)
    
    if (method == "lin") {
      noInter <- 0L
      if (noInterpolAtEnds)
        noInter <- 1L
        res <- .Call(impute_with_linear_interpolation,
                     x,
                     noInter,
                     PACKAGE = "OptiLCMS")
      
      return(res)
    }
    
    if (method == "linbase") {
      if (missing(baseValue))
        baseValue <- min(x, na.rm = TRUE) / 2
      if (!is.double(baseValue))
        baseValue <- as.double(baseValue)
      if (!is.integer(distance))
        distance <- as.integer(distance)

        res <- .Call(impute_with_linear_interpolation_base,
                     x,
                     baseValue,
                     distance,
                     PACKAGE = "OptiLCMS")
      
      return(res)
    }
  }

#' binYonXR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
binYonXR <- function(x,
                    y,
                    breaks,
                    nBins,
                    binSize,
                    binFromX,
                    binToX,
                    fromIdx = 1L,
                    toIdx = length(x),
                    method = "max",
                    baseValue,
                    sortedX = !is.unsorted(x),
                    shiftByHalfBinSize = FALSE,
                    returnIndex = FALSE,
                    returnX = TRUE) {
  if (!missing(x) & missing(y))
    y <- x
  if (missing(x) | missing(y))
    stop("Arguments 'x' and 'y' are mandatory!")
  if (missing(nBins) & missing(binSize) & missing(breaks))
    stop("One of 'breaks', 'nBins' or 'binSize' has to be defined!")
  if (!sortedX) {
    message("'x' is not sorted, will sort 'x' and 'y'.")
    ## Sort method; see issue #180 for MSnbase
    ## Note: order method = "radix" is considerably faster - but there is no
    ## method argument for older R versions.
    o <- order(x)
    x <- x[o]
    y <- y[o]
  }
  ## Check fromIdx and toIdx
  if (any(fromIdx < 1) | any(toIdx > length(x)))
    stop("'fromIdx' and 'toIdx' have to be within 1 and lenght(x)!")
  if (length(toIdx) != length(fromIdx))
    stop("'fromIdx' and 'toIdx' have to have the same length!")
  if (missing(binFromX))
    binFromX <- as.double(NA)
  if (missing(binToX))
    binToX <- as.double(NA)
  ## For now we don't allow NAs in x
  if (anyNA(x))
    stop("No 'NA' values are allowed in 'x'!")
  ## Check that input values have the correct data types.
  if (!is.double(x))
    x <- as.double(x)
  if (!is.double(y))
    y <- as.double(y)
  if (!is.double(binFromX))
    binFromX <- as.double(binFromX)
  if (!is.double(binToX))
    binToX <- as.double(binToX)
  if (!is.integer(fromIdx))
    fromIdx <- as.integer(fromIdx)
  if (!is.integer(toIdx))
    toIdx <- as.integer(toIdx)
  ## breaks has precedence over nBins over binSize.
  shiftIt <- 0L
  if (!missing(breaks)) {
    if (shiftByHalfBinSize)
      warning("Argument 'shiftByHalfBinSize' is ignored if 'breaks'",
              " are provided.")
    if (!is.double(breaks))
      breaks <- as.double(nBins)
    nBins <- NA_integer_
    binSize <- as.double(NA)
  } else {
    if (!missing(nBins)) {
      breaks <- as.double(NA)
      if (!is.integer(nBins))
        nBins <- as.integer(nBins)
      binSize <- as.double(NA)
    } else{
      breaks <- as.double(NA)
      nBins <- NA_integer_
      if (!is.double(binSize))
        binSize <- as.double(binSize)
    }
  }
  if (shiftByHalfBinSize)
    shiftIt <- 1L
  ## Define default value for baseValue
  if (missing(baseValue)) {
    baseValue = as.double(NA)
  } else {
    if (!is.double(baseValue))
      baseValue <- as.double(baseValue)
  }
  
  getIndex <- 0L
  if (returnIndex)
    getIndex <- 1L
  getX <- 0L
  if (returnX)
    getX <- 1L
  if (length(toIdx) > 1) {
      .Call(
        binYonX_multi,
        x,
        y,
        breaks,
        nBins,
        binSize,
        binFromX,
        binToX,
        force(fromIdx - 1L),
        force(toIdx - 1L),
        shiftIt,
        as.integer(.aggregateMethod2int(method)),
        baseValue,
        getIndex,
        getX,
        PACKAGE = "OptiLCMS"
      )

  } else {

      .Call(
        binYonX,
        x,
        y,
        breaks,
        nBins,
        binSize,
        binFromX,
        binToX,
        force(fromIdx - 1L),
        force(toIdx - 1L),
        shiftIt,
        as.integer(.aggregateMethod2int(method)),
        baseValue,
        getIndex,
        getX,
        PACKAGE = "OptiLCMS"
      )
  }
}


#' massifquantROIs
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
massifquantROIs <- function(mz,
                            int,
                            scanindex,
                            scantime,
                            mzrange,
                            scanrange,
                            minIntensity,
                            minCentroids,
                            consecMissedLim,
                            ppm,
                            criticalVal,
                            segs,
                            scanBack) {

  ROIs <-
      .Call(
        massifquant,
        mz,
        int,
        scanindex,
        scantime,
        as.double(mzrange),
        as.integer(scanrange),
        as.integer(length(scantime)),
        as.double(minIntensity),
        as.integer(minCentroids),
        as.double(consecMissedLim),
        as.double(ppm),
        as.double(criticalVal),
        as.integer(segs),
        as.integer(scanBack),
        PACKAGE = 'OptiLCMS'
      )

  return(ROIs)
  
}

#' getEICR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
getEICR <-
  function(mz,
           int,
           scanindex,
           mzrange,
           scanrange) {

      noised <- .Call(
        getEIC,
        mz,
        int,
        scanindex,
        as.double(mzrange),
        as.integer(scanrange),
        as.integer(length(scanindex)),
        PACKAGE = "OptiLCMS"
      )

    return(noised)
  }

#' getMZR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
getMZR <- function(mz,
                  int,
                  scanindex,
                  mzrange,
                  scrange,
                  scantime) {

    omz <- .Call(
      getMZ,
      mz,
      int,
      scanindex,
      as.double(mzrange),
      as.integer(scrange),
      as.integer(length(scantime))
    )

  return(omz)
}


#' findmzROIR
#' @noRd
#' @references Smith, C.A. et al. 2006. {Analytical Chemistry}, 78, 779-787
findmzROIR <- function(mz,
                      int,
                      scanindex,
                      scanrange,
                      scantime,
                      param,
                      minCentroids) {

    roiList <- .Call(
      findmzROI,
      mz,
      int,
      scanindex,
      as.double(c(0.0, 0.0)),
      as.integer(scanrange),
      as.integer(length(scantime)),
      as.double(param$ppm * 1e-6),
      as.integer(minCentroids),
      as.integer(param$prefilter),
      as.integer(param$noise),
      PACKAGE = "OptiLCMS"
    )
    
  return(roiList)
  
}
