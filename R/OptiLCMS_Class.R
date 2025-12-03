# .onAttatch
#
.onAttach <- function (libname, pkgname){
  k1 <- paste("OptiLCMS",
              utils::packageVersion( "OptiLCMS"),
              "initialized Successfully !")
  k0 <- "\n"
  #k2 <- paste("https://github.com/xia-lab/OptiLCMS")
  packageStartupMessage(c(k1,k0))
}

# Internal mSet Class for OptiLCMS
# Internal mSet Class for OptiLCMS
setClass("mSet",
         representation = representation(rawfiles = "character",
                                         rawOnDisk = "OnDiskMSnExp",
                                         rawInMemory = "MSnExp",
                                         params = "list",
                                         peakpicking ="list",
                                         peakgrouping = "list",
                                         peakRTcorrection = "list",
                                         peakfilling = "list",
                                         MSnData = "list",
                                         MSnResults = "list",
                                         peakAnnotation = "list",
                                         dataSet = "data.frame",
                                         runningplan = "list",
                                         msgSet = "list",
                                         userpath = "character",
                                         WorkingDir = "character"),
         prototype = prototype(rawfiles = vector("character"),
                               rawOnDisk = new("OnDiskMSnExp"),
                               rawInMemory = new("MSnExp"),
                               params = vector("list"),
                               peakpicking = vector("list"),
                               peakgrouping = vector("list"),
                               peakRTcorrection = vector("list"),
                               peakfilling = vector("list"),
                               MSnData = vector("list"),
                               MSnResults = vector("list"),
                               peakAnnotation = vector("list"),
                               dataSet = data.frame(),
                               runningplan = vector("list"),
                               msgSet = vector("list"),
                               userpath = character(0),
                               WorkingDir = vector("character")),
         validity = function(object) {
           return(TRUE)
         })


### Generic DataClass Defination----------

# 1. mSetRule were set for CAMERA processing.
##' This class is highly referenced from CAMERA but will improve later
##' @references Kuhl C, Tautenhahn R, Boettcher C, Larson TR, Neumann S (2012). 
##' "CAMERA: an integrated strategy for compound spectra extraction 
##' and annotation of liquid chromatography/mass spectrometry data sets.
##' " Analytical Chemistry, 84, 283-289. 
##' http://pubs.acs.org/doi/abs/10.1021/ac202450g.
#' @noRd
#' 
setClass("mSetRule",
         representation(ionlistfile="character",
                        ionlist="data.frame", 
                        neutrallossfile="character",
                        neutralloss="data.frame", 
                        neutraladditionfile="character",
                        neutraladdition="data.frame",
                        maxcharge="numeric",
                        mol="numeric",
                        nion="numeric",
                        nnloss="numeric",
                        nnadd="numeric",
                        nh="numeric",
                        polarity="character",
                        rules="data.frame",
                        lib.loc="character"),        
         contains=c("Versioned"),
         prototype=prototype(
           ionlistfile="",
           ionlist=data.frame(),
           neutrallossfile="",
           neutralloss=data.frame(),           
           neutraladditionfile="",
           neutraladdition=data.frame(),
           maxcharge=numeric(),
           mol=numeric(),
           nion=numeric(),
           nnloss=numeric(),
           nnadd=numeric(),
           nh=numeric(),
           polarity=NULL,
           rules=data.frame(),
           lib.loc=NULL,
           new("Versioned", versions=c(ruleSet="0.1.0"))),
         validity=function(object) {
           TRUE
         })

setClass(
  "OptiCommandSet",
  representation(
    ROIExtraction = "language",
    ParamsOptimization = "language",
    ImportRawMSData = "language",
    PeakProfiling = "language",
    PeakAnnotation = "language",
    FormatPeakList = "language"
  ),
  prototype = prototype(
    ROIExtraction = new("call"),
    ParamsOptimization = new("call"),
    ImportRawMSData = new("call"),
    PeakProfiling = new("call"),
    PeakAnnotation = new("call"),
    FormatPeakList = new("call")
  ),
  validity = function(object) {
    TRUE
  }
)

setClass(
  "CustCommandSet",
  representation(
    ImportRawMSData = "language",
    PeakProfiling = "language",
    PeakAnnotation = "language",
    FormatPeakList = "language"
  ),
  prototype = prototype(
    ImportRawMSData = new("call"),
    PeakProfiling = new("call"),
    PeakAnnotation = new("call"),
    FormatPeakList = new("call")
  ),
  validity = function(object) {
    TRUE
  }
)

setClass(
  "controller",
  representation(
    ROI_extract = "logical",
    data_import = "logical",
    peak_profiling = "logical",
    peak_annotation = "logical",
    others_1 = "logical",
    operators = "logical"
  ),
  prototype = prototype(
    ROI_extract = logical(),
    data_import = logical(),
    peak_profiling = logical(),
    peak_annotation = logical(),
    others_1 = logical(),
    operators = logical()
  ),
  validity = function(object) {
    TRUE
  }
)

setClass(
  "ResumeHistory",
  representation(setNO = "numeric",
                 FunishedPosition = "numeric"),
  prototype = prototype(
    setNO = vector("numeric"),
    FunishedPosition = vector("numeric")
  ),
  validity = function(object) {
    TRUE
  }
)

#' Running Plan Class is designed for fast resuming
#' @author Zhiqiang Pang
#' @noRd
#' 
setClass(
  "ResumingePlan",
  representation(
    running.controller = "controller",
    CommandSet = "list",
    RunningHistory = "ResumeHistory",
    PlanNumber = "numeric",
    WorkingDir = "character"
  ),
  prototype = prototype(
    running.controller = new("controller"),
    CommandSet = vector("list"),
    RunningHistory = new("ResumeHistory"),
    PlanNumber = numeric(0),
    WorkingDir = NULL,
    new("Versioned", versions = c(OptiLCMS = "0.1.0"))
  ),
  validity = function(object) {
    TRUE
  }
)

#' Running Plan Class is designed for fast resuming
#' @author Zhiqiang Pang
#' @noRd
#' 
setClass(
  "cacheEnvi",
  representation(
    cache = "list",
    envir = "list",
    plan = "list",
    records = "list"
  ),
  prototype = prototype(
    cache = vector("list"),
    envir = vector("list"),
    plan = vector("list"),
    records = vector("list"),
    new("Versioned", versions = c(OptiLCMS = "0.1.0"))
  ),
  validity = function(object) {
    TRUE
  }
)

#' Running Group Class
#' @author Zhiqiang Pang
#' @noRd
#' 
setClass("NAnnotatedDataFrame",
         representation(multiplex = "numeric",
                        multiLabels = "character"),
         contains = c("AnnotatedDataFrame"),
         prototype = prototype(
           new("Versioned", versions = list(NAnnotatedDataFrame="0.0.3")),
           multiplex = 1,
           multiLabels = "Single run"),
         validity = function(object) {
           msg <- validMsg(NULL, NULL)
           if (length(object@multiLabels) != object@multiplex)
             msg <- validMsg(msg, "Number of multiplex does not match it's labels.")
           if (is.null(msg)) TRUE
           else msg
         })

