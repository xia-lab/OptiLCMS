#'Constructs a dataSet object for storing data 
#'@description This functions handles the construction of a mSetObj object for storing data for further processing and analysis.
#'It is necessary to utilize this function to specify to MetaboAnalystR the type of data and the type of analysis you will perform. 
#'@usage InitDataObjects(data.type, anal.type, paired=FALSE)
#'@param data.type The type of data, either list (Compound lists), conc (Compound concentration data), 
#'specbin (Binned spectra data), pktable (Peak intensity table), nmrpeak (NMR peak lists), mspeak (MS peak lists), 
#'or msspec (MS spectra data)
#'@param anal.type Indicate the analysis module to be performed: stat, pathora, pathqea, msetora, msetssp, msetqea, ts, 
#'cmpdmap, smpmap, or pathinteg
#'@param paired Indicate if the data is paired or not. Logical, default set to FALSE
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import methods
#'@importFrom  Cairo CairoFonts

InitDataObjects <- function(data.type, anal.type, paired=FALSE){
  
  if(anal.type == "raw" & data.type == "spec") {
    cat("OptiLCMS R objects initialized ...\n");
    return(new("mSet"))
  }
  
  if(!.on.public.web){
    if(exists("mSet")){
      mSetObj <- .get.mSet(mSet);
      mSetObj$dataSet$type <- data.type;
      mSetObj$analSet$type <- anal.type;
      return(.set.mSet(mSetObj));
    }
  }
  
  dataSet <- list();
  dataSet$type <- data.type;
  dataSet$design.type <- "regular"; # one factor to two factor
  dataSet$cls.type <- "disc"; # default until specified otherwise
  dataSet$format <- "rowu";
  dataSet$paired <- paired;
  analSet <- list();
  analSet$type <- anal.type;
  
  mSetObj <- list();
  mSetObj$dataSet <- dataSet;
  mSetObj$analSet <- analSet;
  mSetObj$imgSet <- list();
  mSetObj$msgSet <- list(); # store various message during data processing
  mSetObj$msgSet$msg.vec <- vector(mode="character");     # store error messages
  mSetObj$cmdSet <- vector(mode="character"); # store R command
  
  # other global variables
  msg.vec <<- "";
  err.vec <<- "";
  
  # for network analysis
  module.count <<- 0;
  # counter for naming different json file (pathway viewer)
  smpdbpw.count <<- 0; 
  # for mummichog
  peakFormat <<- "mpt"  
  
  # for meta-analysis
  mdata.all <<- list(); 
  mdata.siggenes <<- vector("list");
  meta.selected <<- TRUE;
  anal.type <<- anal.type;
  
  if(.on.public.web){
    # disable parallel prcessing for public server
    load_BiocParallel();
    register(SerialParam());
  }else{
    if("stat" %in% anal.type | "msetqea" %in% anal.type | "pathqea" %in% anal.type | "roc" %in% anal.type)
      # start Rserve engine for Rpackage
      load_Rserve();
  }
  
  # plotting required by all
  Cairo::CairoFonts(regular="Arial:style=Regular",bold="Arial:style=Bold",italic="Arial:style=Italic",bolditalic = "Arial:style=Bold Italic",symbol = "Symbol")
  
  # sqlite db path for gene annotation
  if(file.exists("/home/glassfish/sqlite/")){ #.on.public.web
    url.pre <<- "/home/glassfish/sqlite/";
  }else if(file.exists("/home/jasmine/Downloads/sqlite/")){ #jasmine's local
    url.pre <<- "/home/jasmine/Downloads/sqlite/";
  }else if(file.exists("/Users/soufanom/Documents/Projects/gene-id-mapping/")){ # soufan laptop
    url.pre <<- "/Users/soufanom/Documents/Projects/gene-id-mapping/";
  }else if(file.exists("~/Documents/Projects/gene-id-mapping/")){
    url.pre <<- "~/Documents/Projects/gene-id-mapping/"
  }else if(file.exists("/Users/xia/Dropbox/sqlite/")){ # xia local
    url.pre <<- "/Users/xia/Dropbox/sqlite/";
  }else if(file.exists("/media/zzggyy/disk/sqlite/")){
    url.pre <<-"/media/zzggyy/disk/sqlite/"; #zgy local)
  }else if(file.exists("/home/zgy/sqlite/")){
    url.pre <<-"/home/zgy/sqlite/"; #zgy local)
  } else if(file.exists("/home/le/sqlite/GeneID_25Species_JE/")){# le local
    url.pre <<-"/home/le/sqlite/GeneID_25Species_JE/";
  }else{
    url.pre <<- paste0(dirname(system.file("database", "sqlite/GeneID_25Species_JE/ath_genes.sqlite", package="MetaboAnalystR")), "/")
  }
  
  print("MetaboAnalyst R objects initialized ...");
  return(.set.mSet(mSetObj));
}

#' Import raw MS data
#' @description This function handles the reading in of
#' raw MS data (.mzML, .CDF and .mzXML). Users must provide
#' a matrix with meta information about file such that each file has the name,
#' file path, group class and extension type.
#' The function will output two chromatograms into the user's working directory, a
#' base peak intensity chromatogram (BPIC) and a total ion
#' chromatogram (TIC). Further, this function sets the number of cores
#' to be used for parallel processing. It first determines the number of cores
#' within a user's computer and then sets it that number/2.
#' @param dataset.meta Matrix, input the meta data for files containing
#' the raw MS spectra to be processed.
#' @param format Character, input the format of the image to create.
#' @param dpi Numeric, input the dpi of the image to create.
#' @param width Numeric, input the width of the image to create.
#' @param par.cores Logical, if true, the function will automatically
#' set the number of parallel cores. If false, it will not.
#' @param plot Logical, if true the function will create BPIS and TICS plots.
#' @param bpis_name Character, input the name of the BPIS image to create.
#' @param tics_name Character, input the name of the TICS image to create.
#' @author Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' Mai Yamamoto \email{yamamoto.mai@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import MSnbase
#' @import BiocParallel
#' @import parallel
#' @importFrom Cairo Cairo

ImportRawMSDataList <-
  function(dataset.meta,
           format = "png",
           dpi = 72,
           width = 9,
           par.cores = TRUE,
           plot = TRUE,
           bpis_name = "BPIS_",
           tics_name = "TICS_") {
    
    msg.vec <<- vector(mode = "character")
    
    if (bpis_name == "BPIS_") {
      bpis_name = paste("BPIS_", dpi, ".", format, sep = "")
    }
    if (tics_name == "TICS_") {
      tics_name <- paste("TICS_", dpi, ".", format, sep = "")
    }
    
    msg <- c("The uploaded files are raw MS spectra.")
    
    # The dataset.meta should be a matrix such that each row has the following:
    #   1- file name, 2- file path, 3- group and 4- file extension type
    # provided. The accepted file extensions are (.mzML/.CDF/.mzXML files)
    
    if (nrow(dataset.meta) == 0 || is.na(dataset.meta)) {
      AddErrMsg("No spectra were found!")
      return(0)
    }
    
    compfile.types <-
      sum(sapply(dataset.meta[, 3], function(x) {
        x %in% c("mzml", "cdf", "mzxml")
      }))
    if (compfile.types < nrow(dataset.meta)) {
      AddErrMsg("Only mzML, cdf and mzXML input types can be handled!")
      return(0)
    }
    
    snames <- dataset.meta[, 1]
    files <-
      dataset.meta[, 2]
    files <-
      as.character(files)
    # Otherwise, a factor form of files will cause an error
    sclass <- dataset.meta[, 4]
    
    # some sanity check before proceeds
    sclass <- as.factor(sclass)
    
    if (length(levels(sclass)) < 2) {
      AddErrMsg("You must provide classes labels (at least two classes)!")
      return(0)
    }
    
    SetClass(sclass)
    
    # check for unique sample names
    if (length(unique(snames)) != length(snames)) {
      AddErrMsg("Duplicate sample names are not allowed!")
      dup.nm <- paste(snames[duplicated(snames)], collapse = " ")
      AddErrMsg("Duplicate sample names are not allowed!")
      AddErrMsg(dup.nm)
      return(0)
    }
    
    pd <- data.frame(
      sample_name = snames,
      sample_group = sclass,
      stringsAsFactors = FALSE
    )
    
    if (!.on.public.web & par.cores == TRUE) {
      cores <- parallel::detectCores()
      num_cores <- ceiling(cores / 2)
      cat(paste0("The number of CPU cores to be used is set to ", num_cores, ".","\n"))
      
      if (.Platform$OS.type == "unix") {
        BiocParallel::register(BiocParallel::bpstart(BiocParallel::MulticoreParam(num_cores)))
      } else {
        # for windows
        BiocParallel::register(BiocParallel::bpstart(BiocParallel::SnowParam(num_cores)))
      }
    }
    
    raw_data <-
      suppressMessages(read.MSdata(
        files = files,
        pdata = new("NAnnotatedDataFrame", pd),
        mode = "onDisk"
      ))
    
    if (plot == TRUE) {
      # Plotting functions to see entire chromatogram
      bpis <- chromatogram(raw_data, aggregationFun = "max")
      tics <- chromatogram(raw_data, aggregationFun = "sum")
      
      groupNum <- nlevels(groupInfo)
      
      if (groupNum > 9) {
        col.fun <-
          grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
        group_colors <- col.fun(groupNum)
      } else{
        group_colors <-
          paste0(RColorBrewer::brewer.pal(9, "Set1")[1:groupNum], "60")
      }
      
      names(group_colors) <- levels(groupInfo)
      
      Cairo::Cairo(
        file = bpis_name,
        unit = "in",
        dpi = dpi,
        width = width,
        height = width * 5 / 9,
        type = format,
        bg = "white"
      )
      
      plot(bpis, col = group_colors[raw_data$sample_group])
      legend(
        "topright",
        legend = levels(groupInfo),
        pch = 15,
        col = group_colors
      )
      
      dev.off()
      
      Cairo::Cairo(
        file = tics_name,
        unit = "in",
        dpi = dpi,
        width = width,
        height = width * 5 / 9,
        type = format,
        bg = "white"
      )
      
      plot(tics, col = group_colors[raw_data$sample_group])
      legend(
        "topright",
        legend = levels(groupInfo),
        pch = 15,
        col = group_colors
      )
      
      dev.off()
    }
    
    cat("Successfully imported raw MS data!\n")
    return(raw_data)
  }

#' Import raw MS data
#' @description This function handles the reading in of
#' raw MS data (.mzML, .CDF and .mzXML). Users must set
#' their working directory to the folder containing their raw
#' data, divided into two subfolders named their desired group labels. The
#' function will output two chromatograms into the user's working directory, a
#' base peak intensity chromatogram (BPIC) and a total ion
#' chromatogram (TIC). Further, this function sets the number of cores
#' to be used for parallel processing. It first determines the number of cores
#' within a user's computer and then sets it that number/2.
#' @param mSet mSet Object, can be optional. Usually generated by InitDataObjects("spec", "raw", FALSE) before the data import.
#' @param foldername Character, input the file path to the folder containing
#' the raw MS spectra to be processed.
#' @param mode Character, the data input mode. Default is "onDisk" to avoid memory crash. "inMemory" will
#' absorb data into the memory.
#' @param ncores Numeric, a value used to defined the parallel cores.
#' @param plotSettings List, plotting parameters produced by SetPlotParam Function. "plot.opts" can be added through this
#' function for samples numbers for plotting. Defalut is "default", "all" will apply all samples for plotting and may cause
#' memory crash, especially for large sample dataset.
#' @param running.controller The resuming pipeline running controller. Optional. Don't need to define by hand.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' Mai Yamamoto \email{yamamoto.mai@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import MSnbase
#' @import BiocParallel
#' @import parallel
#' @importFrom Cairo Cairo

ImportRawMSData <-
  function(mSet = NULL,
           foldername,
           mode = "onDisk",
           ncores = 4,
           plotSettings,
           running.controller = NULL) {
    
    foldername <- tools::file_path_as_absolute(foldername);
    
    if(.on.public.web){
      load_msnbase();
    } 
    
    if(missing(mSet)){
      mSet <- new("mSet")
    } else if(is.null(mSet)){
      mSet <- new("mSet")
    }
    
    #Build Running plan for data import - Indentify the controller
    if (is.null(running.controller)) {
      c1 <- TRUE;
      c2 <- TRUE;
      plan_switch <- FALSE;
    } else {
      plan_switch <- TRUE;
      c1 <-
        running.controller@data_import[["c1"]] # used to control data import
      c2 <-
        running.controller@data_import[["c2"]] # used to control plotting option
    }
    
    .optimize_switch <<- FALSE;
    
    if (!dir.exists(foldername) & .on.public.web) {
      foldername <- "/home/glassfish/projects/MetaboDemoRawData/upload"
    }
    
    start.time <- Sys.time()
    msg.vec <<- vector(mode = "character")
    msg <- c("The uploaded files are raw MS spectra.")
    
    # the "upload" folder should contain two subfolders (groups, i.e. Healthy vs. Disease)
    # each subfolder must contain samples (.mzML/.CDF/.mzXML files)
    
    MessageOutput(mes = "Step 2/6: Start to import the spectrum! \nThis step will take a short time...",
                  ecol = "\n",
                  progress = 21.0)
    
    files <-
      dir(
        foldername,
        pattern = ".mzML|.mzml|.cdf|.mzXML|.mzxml|.mzData|.CDF",
        recursive = T,
        full.names = TRUE
      )
    
    if (length(files) == 0) {
      MessageOutput(
        mes = paste0(
          "<font color=\"red\">",
          "\nERROR: No standard MS file found ! Please check the extension of your data.",
          "</font>"
        ),
        ecol = "\n",
        progress = NULL
      )
      stop()
      
    } else if (length(files) < 3) {
      MessageOutput(
        mes = paste0(
          "<font color=\"red\">",
          "\nERROR: At least 3 samples should be provided.",
          "</font>"
        ),
        ecol = "\n",
        progress = NULL
      )
    }
    
    count_total_sample <<- length(files)
    count_current_sample <<- 0
    toRemove = vector();
    
    # Update first
    if(length(mSet@rawfiles) == 0){
      rawfilenms <- basename(files);
    } else {
      rawfilenms <- basename(mSet@rawfiles);
    }
    
    for (i in 1:length(files)) {
      file = basename(files[i])
      if (!(file %in% rawfilenms)) {
        toRemove = c(toRemove, files[i])
      }
    }
    
    toKeepInx = !(files %in% toRemove)
    files = files[toKeepInx]
    
    if (length(files) == 0) {
      AddErrMsg("No spectra were found!")
      return(0)
    }
    
    snames <- gsub("\\.[^.]*$", "", basename(files))
    msg <- c(msg, paste("A total of ", length(files), "samples were found."))
    sclass <- gsub("^\\.$", "sample", dirname(files))
    
    scomp <- strsplit(substr(sclass, 1, min(nchar(sclass))), "", fixed = TRUE)
    scomp <- matrix(c(scomp, recursive = TRUE), ncol = length(scomp))
    
    i <- 1
    while (all(scomp[i, 1] == scomp[i, -1]) && i < nrow(scomp)) {
      i <- i + 1
    }
    
    i <-
      min(i, tail(c(0, which(
        scomp[1:i, 1] == .Platform$file.sep
      )), n = 1) + 1)
    
    if (i > 1 && i <= nrow(scomp)) {
      sclass <- substr(sclass, i, max(nchar(sclass)))
    }
    
    if (.on.public.web &
        unique(sclass)[1] == "upload" &
        length(unique(sclass)) == 1) {
      sclass <- rep("Unknown", length(sclass))
    }
    # some sanity check before proceeds
    sclass <- as.factor(sclass);
    SetClass(sclass);
    
    # check for unique sample names
    if (length(unique(snames)) != length(snames)) {
      AddErrMsg("Duplicate sample names are not allowed!")
      dup.nm <- paste(snames[duplicated(snames)], collapse = " ")
      AddErrMsg("Duplicate sample names are not allowed!")
      AddErrMsg(dup.nm)
      return(0)
    }
    
    pd <- data.frame(
      sample_name = snames,
      sample_group = sclass,
      stringsAsFactors = FALSE
    )
    
    if (!.on.public.web) {
      cores <- parallel::detectCores()
      if (missing(ncores)) {
        num_cores <- ceiling(cores * 2 / 3)
      } else{
        ncores -> num_cores
      }
      
      cat(paste0("The number of CPU cores to be used is set to ", num_cores, ".","\n"))
      
      if (.Platform$OS.type == "unix") {
        BiocParallel::register(BiocParallel::bpstart(BiocParallel::MulticoreParam(num_cores)))
      } else {
        # for windows
        BiocParallel::register(BiocParallel::bpstart(BiocParallel::SnowParam(num_cores)))
      }
      
    } else {
      num_cores <- 2
      cat(paste0("The number of CPU cores to be used is set to ", num_cores, ".","\n"))
      BiocParallel::register(BiocParallel::bpstart(BiocParallel::MulticoreParam(num_cores)))
    }
    
    if (c1) {
      raw_data <-
        suppressMessages(read.MSdata(
          files = files,
          pdata = new("NAnnotatedDataFrame", pd),
          mode = mode,
          msLevel. = 1
        ))
      
      if(plan_switch){
        cache.save(raw_data, funpartnm = "data_import_c1");
        marker_record("data_import_c1");
      }
 
    } else {
      raw_data <- cache.read ("data_import", "c1")
      marker_record("data_import_1")
    }
    
    MessageOutput(NULL, NULL, 22)
    
    if (c2) {
      if (plotSettings$Plot == TRUE) {
        if (is.null(plotSettings$plot.opts)) {
          plot.opts <- "default"
        } else {
          plot.opts <- plotSettings$plot.opts
        }
        
        if (plot.opts == "default") {
          #subset raw_data to first 50 samples
          cat("To reduce memory usage BPIS and TICS plots will be created using only 10 samples per group.\n")
          
          grp_nms <- names(table(pd$sample_group))
          files <- NA
          
          for (i in 1:length(grp_nms)) {
            numb2ext <- min(table(pd$sample_group)[i], 10)
            filt_df <- pd[pd$sample_group == grp_nms[i], ]
            files.inx <- sample(nrow(filt_df), numb2ext)
            sel.samples <- filt_df$sample_name[files.inx]
            files <-
              c(files, which(pd$sample_name %in% sel.samples))
          }
          
          raw_data_filt <-
            filterFile(raw_data, file = na.omit(files))
          
        } else{
          raw_data_filt <- raw_data
          # just for plotting
        }
        
        save(raw_data_filt, file = "raw_data_filt.rda")
        
        if (plot.opts == "all") {
          h <-
            readline(prompt = "Using all samples to create BPIS and TICS plots may cause severe memory issues! Press [0] to continue, or [1] to cancel: ")
          h <- as.integer(h)
          
          if (h == 1) {
            cat("ImportRawMSData function aborted!\n")
            return(0)
          }
        }
        
        cat("Plotting BPIS and TICS.\n")
        # Plotting functions to see entire chromatogram
        bpis <- chromatogram(raw_data_filt, aggregationFun = "max")
        tics <- chromatogram(raw_data_filt, aggregationFun = "sum")
        
        groupNum <- nlevels(groupInfo)
        
        if (groupNum > 9) {
          col.fun <-
            grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
          group_colors <- col.fun(groupNum)
          
        } else{
          group_colors <-
            paste0(RColorBrewer::brewer.pal(9, "Set1")[1:groupNum], "60")
        }
        
        names(group_colors) <- levels(groupInfo)
        
        bpis_name <-
          paste("BPIS_",
                plotSettings$dpi,
                ".",
                plotSettings$format,
                sep = "")
        
        tics_name <-
          paste("TICS_",
                plotSettings$dpi,
                ".",
                plotSettings$format,
                sep = "")
        
        #save(bpis, file = "bpis.rda"); # Don't need bpis for now.
        save(tics, file = "tics.rda")
        
        Cairo::Cairo(
          file = bpis_name,
          unit = "in",
          dpi = plotSettings$dpi,
          width = 8,
          height = 6,
          type = plotSettings$format,
          bg = "white"
        )
        
        plot(bpis, col = group_colors[raw_data_filt$sample_group])
        legend(
          "topright",
          legend = levels(groupInfo),
          pch = 15,
          col = group_colors
        )
        
        dev.off()
        
        Cairo::Cairo(
          file = tics_name,
          unit = "in",
          dpi = plotSettings$dpi,
          width = 8,
          height = 6,
          type = plotSettings$format,
          bg = "white"
        )
        
        plot(tics, col = group_colors[raw_data_filt$sample_group])
        legend(
          "topright",
          legend = levels(groupInfo),
          pch = 15,
          col = group_colors
        )
        
        dev.off()
      }
      if (plan_switch) {
        marker_record("data_import_c2")
      }
    }
    
    MessageOutput(
      mes = paste0(
        "Step 2/6: Successfully imported raw MS data! (",
        Sys.time(),
        ") \nGoing to the next step..."
      ),
      ecol = "\n",
      progress = 24
    )
    
    if(mode == "onDisk"){
      mSet@rawOnDisk <- raw_data;
    } else if( mode == "inMemory"){
      mSet@rawInMemory <- raw_data;
    }
      
    return(mSet)
  }


updateSpectraFiles <- function(mSet, workingDir, filesNames){
  
  if(missing(mSet) | .on.public.web){
    load("mSet.rda");
  }
  
  if(missing(workingDir)){
    mSet@WorkingDir <- getwd();
  }
  
  if(!missing(filesNames)){
    mSet@rawfiles <- filesNames;  
  }
  
  if(.on.public.web){
    save(mSet, file = "mSet.rda")
  } else {
    return(mSet)
  }
  
}

read.MSdata <- function(files, 
                        pdata = NULL, 
                        msLevel. = NULL, 
                        centroided. = NA,
                        smoothed. = NA, 
                        cache. = 1L,
                        mode = c("inMemory", "onDisk")) {
  mode <- match.arg(mode)
  ## o normalize the file path, i.e. replace relative path with absolute
  ##   path. That fixes possible problems on Windows with SNOW parallel
  ##   processing and also proteowizard problems on unis system with ~ paths.
  files <- normalizePath(files)
  suppressWarnings(.hasChroms <- MSnbase::hasChromatograms(files))
  
  MessageOutput ("Raw file import begin...", "\n", NULL)
  
  if (!length(files)) {
    process <- new("MSnProcess",
                   processing = paste("No data loaded:", date()))
    if (mode == "inMemory")
      res <- new("MSnExp",
                 processingData = process)
    else res <- new("OnDiskMSnExp",
                    processingData = process)
  } else {
    if (mode == "inMemory") {
      if (is.null(msLevel.)) msLevel. <- 2L
      res <- read.InMemMSd.data(files, pdata = pdata, msLevel. = msLevel.,
                                centroided. = centroided., smoothed. = smoothed., cache. = cache.)
    } else { ## onDisk
      res <- read.OnDiskMS.data(files = files, pdata = pdata,
                                msLevel. = msLevel., centroided. = centroided., smoothed. = smoothed.)
    }
  }
  res
}

#' @import utils
read.InMemMSd.data <- function(files, 
                               pdata, 
                               msLevel., 
                               centroided., 
                               smoothed., 
                               cache. = 1) {
  MSnbase:::.testReadMSDataInput(environment())
  if (MSnbase:::isCdfFile(files)) {
    #message("Polarity can not be extracted from netCDF files, please set ",
    #        "manually the polarity with the 'polarity' method.")
    msLevel. <- 1;
  }
  
  if (msLevel. == 1) cache. <- 0
  msLevel. <- as.integer(msLevel.)
  ## Creating environment with Spectra objects
  assaydata <- new.env(parent = emptyenv())
  ioncount <- c()
  ioncounter <- 1
  filenams <- filenums <- c()
  fullhd2 <- fullhdorder <- c()
  fullhdordercounter <- 1
  .instrumentInfo <- list()
  ## List eventual limitations
  
  ## ## Idea:
  ## ## o initialize a featureData-data.frame,
  ## ## o for each file, extract header info and put that into
  ##      featureData;
  
  count.idx <- 0;
  
  for (f in files) {
    cat(paste("Reading MS from",basename(f),"begin !\n"))
    
    filen <- match(f, files)
    filenums <- c(filenums, filen)
    filenams <- c(filenams, f)
    ## issue #214: define backend based on file format.
    msdata <- mzR::openMSfile(f,backend = NULL)
    .instrumentInfo <- c(.instrumentInfo, list(mzR::instrumentInfo(msdata)))
    fullhd <- mzR::header(msdata)
    ## Issue #325: get centroided information from file, but overwrite if
    ## specified with centroided. parameter.
    if (!is.na(centroided.))
      fullhd$centroided <- as.logical(centroided.)
    spidx <- which(fullhd$msLevel == msLevel.)
    ## increase vectors as needed
    ioncount <- c(ioncount, numeric(length(spidx)))
    fullhdorder <- c(fullhdorder, numeric(length(spidx)))
    if (msLevel. == 1) {
      if (length(spidx) == 0)
        stop("No MS1 spectra in file",f)
      
      
      if (.on.public.web){   
        print_mes <- paste0("Importing ",basename(f),":");    
        write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = " ");
      }
      else {
        pb <- progress_bar$new(format = "Reading [:bar] :percent Time left: :eta", 
                               total = length(spidx), clear = T, width= 75)
      }
      
      k_count <- 0;
      for (i in 1:length(spidx)) {
        
        if (!.on.public.web){   
          pb$tick();
        }
        
        if (.on.public.web){  
          if (round(i/length(spidx),digits = 4)*100 - k_count > -0.2){
            print_mes <- paste0(k_count,"% | ");    
            write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = " ");
            k_count <- k_count +20;
          }
          
        }
        
        j <- spidx[i]
        hd <- fullhd[j, ]
        ## Fix missing polarity from netCDF
        pol <- hd$polarity
        if (length(pol) == 0)
          pol <- NA
        .p <- mzR::peaks(msdata, j)
        sp <- new("Spectrum1",
                  rt = hd$retentionTime,
                  acquisitionNum = as.integer(hd$acquisitionNum),
                  scanIndex = as.integer(hd$seqNum),
                  tic = hd$totIonCurrent,
                  mz = .p[, 1],
                  intensity = .p[, 2],
                  fromFile = as.integer(filen),
                  centroided = as.logical(hd$centroided),
                  smoothed = as.logical(smoothed.),
                  polarity = as.integer(pol))
        ## peaksCount
        ioncount[ioncounter] <- sum(.p[, 2])
        ioncounter <- ioncounter + 1
        .fname <-MSnbase:::formatFileSpectrumNames(fileIds=filen,
                                                   spectrumIds=i,
                                                   nSpectra=length(spidx),
                                                   nFiles=length(files))
        assign(.fname, sp, assaydata)
        fullhdorder[fullhdordercounter] <- .fname
        fullhdordercounter <- fullhdordercounter + 1
      }
    } else { ## .msLevel != 1
      if (length(spidx) == 0)
        stop("No MS(n>1) spectra in file", f)
      cat(paste("Reading ", length(spidx), " MS", msLevel.,
                  " spectra from file ", basename(f),"\n"))
      
      scanNums <- fullhd[fullhd$msLevel == msLevel., "precursorScanNum"]
      if (length(scanNums) != length(spidx))
        stop("Number of spectra and precursor scan number do not match!")
      
      pb <- progress_bar$new(format = "Reading [:bar] :percent Time left: :eta", 
                             total = length(spidx), clear = T, width= 75)
      
      for (i in 1:length(spidx)) {
        
        pb$tick();
        
        j <- spidx[i]
        hd <- fullhd[j, ]
        .p <- mzR::peaks(msdata, j)
        sp <- new("Spectrum2",
                  msLevel = as.integer(hd$msLevel),
                  merged = as.numeric(hd$mergedScan),
                  precScanNum = as.integer(scanNums[i]),
                  precursorMz = hd$precursorMZ,
                  precursorIntensity = hd$precursorIntensity,
                  precursorCharge = as.integer(hd$precursorCharge),
                  collisionEnergy = hd$collisionEnergy,
                  rt = hd$retentionTime,
                  acquisitionNum = as.integer(hd$acquisitionNum),
                  scanIndex = as.integer(hd$seqNum),
                  tic = hd$totIonCurrent,
                  mz = .p[, 1],
                  intensity = .p[, 2],
                  fromFile = as.integer(filen),
                  centroided = as.logical(hd$centroided),
                  smoothed = as.logical(smoothed.),
                  polarity = as.integer(hd$polarity))
        ## peaksCount
        ioncount[ioncounter] <- sum(.p[, 2])
        ioncounter <- ioncounter + 1
        .fname <- MSnbase:::formatFileSpectrumNames(fileIds=filen,
                                                    spectrumIds=i,
                                                    nSpectra=length(spidx),
                                                    nFiles=length(files))
        assign(.fname, sp, assaydata)
        fullhdorder[fullhdordercounter] <- .fname
        fullhdordercounter <- fullhdordercounter + 1
      }
    }
    if (cache. >= 1)
      fullhd2 <- rbind(fullhd2, fullhd[spidx, ])
    
    gc()
    mzR::close(msdata)
    rm(msdata);
    
    if (.on.public.web){ 
      print_mes <- paste0("Done!");    
      write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
      
      count.idx <- count.idx + 1;  
      write.table(1.0 + count.idx/length(files)*3, file = paste0(fullUserPath, "log_progress.txt"),row.names = F,col.names = F);
    }
    
    cat(paste("This reading finished !\n"))
    
  }
  
  ## cache level 2 yet implemented
  cache. <- MSnbase:::testCacheArg(cache., maxCache = 2)
  if (cache. >= 1) {
    fl <- sapply(assaydata, function(x) x@fromFile)
    featnms <- ls(assaydata) ## feature names in final MSnExp
    fl <- fl[featnms] ## reorder file numbers
    stopifnot(all(base::sort(featnms) == base::sort(fullhdorder)))
    fullhdorder <- match(featnms, fullhdorder)
    tmphd <- fullhd2[fullhdorder, ] ## reorder
    ioncount <- ioncount[fullhdorder]
    newhd <- data.frame(fileIdx = fl,
                        retention.time = tmphd$retentionTime,
                        precursor.mz = tmphd$precursorMZ,
                        precursor.intensity = tmphd$precursorIntensity,
                        charge = tmphd$precursorCharge,
                        peaks.count = tmphd$peaksCount,
                        tic = tmphd$totIonCurrent,
                        ionCount = ioncount,
                        ms.level = tmphd$msLevel,
                        acquisition.number = tmphd$acquisitionNum,
                        collision.energy = tmphd$collisionEnergy)
  } else {
    newhd <- NULL ## not used anyway
  }
  .cacheEnv <- MSnbase:::setCacheEnv(list("assaydata" = assaydata,
                                          "hd" = newhd),
                                     cache., lock = TRUE)
  ## CACHING AS BEEN SUPERSEDED BY THE OnDiskMSnExp IMPLEMENTATION
  ## if cache==2, do not lock assign msdata in .cacheEnv then lock
  ## it and do not close(msdata) above; rm(msdata) is OK
  
  ## Create 'MSnProcess' object
  process <- new("MSnProcess",
                 processing = paste("Data loaded:", date()),
                 files = files,
                 smoothed = smoothed.)
  ## Create 'fdata' and 'pdata' objects
  nms <- ls(assaydata)
  if (is.null(pdata)) {
    .pd <- data.frame(sampleNames = basename(files))
    rownames(.pd) <- .pd$sampleNames
    pdata <- new("AnnotatedDataFrame",
                 data = .pd)
  }
  fdata <- new("AnnotatedDataFrame",
               data = data.frame(
                 spectrum = 1:length(nms),
                 row.names = nms))
  fdata <- fdata[ls(assaydata)] ## reorder features
  ## expriment data slot
  if (length(.instrumentInfo) > 1) {
    cmp <- length(unique(sapply(.instrumentInfo, "[[", 1)))
    if (cmp > 1)
      message("According to the instrument information in the files,\n",
              "the data has been acquired on different instruments!")
    for (nm in names(.instrumentInfo[[1]]))
      .instrumentInfo[[1]][[nm]] <- sapply(.instrumentInfo, "[[", nm)
  }
  expdata <- new("MIAPE",
                 instrumentManufacturer = .instrumentInfo[[1]]$manufacturer,
                 instrumentModel = .instrumentInfo[[1]]$model,
                 ionSource = .instrumentInfo[[1]]$ionisation,
                 analyser = as.character(.instrumentInfo[[1]]$analyzer),
                 detectorType = .instrumentInfo[[1]]$detector)
  ## Create and return 'MSnExp' object
  
  toReturn <- new("MSnExp",
                  assayData = assaydata,
                  phenoData = pdata,
                  featureData = fdata,
                  processingData = process,
                  experimentData = expdata,
                  .cache = .cacheEnv)
  return(toReturn)
}

read.OnDiskMS.data <- function(files, 
                               pdata, 
                               msLevel., 
                               centroided., 
                               smoothed.) {
  
  MSnbase:::.testReadMSDataInput(environment())
  stopifnot(is.logical(centroided.))
  
  ## Creating environment with Spectra objects
  assaydata <- new.env(parent = emptyenv())
  filenams <- filenums <- c()
  fullhd2 <- fullhdorder <- c()
  fullhdordercounter <- 1
  .instrumentInfo <- list()
  ## List eventual limitations
  if (MSnbase:::isCdfFile(files)) {
    message("Polarity can not be extracted from netCDF files, please set ",
            "manually the polarity with the 'polarity' method.")
  }
  ## Idea:
  ## o initialize a featureData-data.frame,
  featureDataList <- list()
  ## o for each file, extract header info and put that into featureData
  ##pb <- progress_bar$new(format = "Reading [:bar] :percent Time left: :eta", 
  #                       total = length(spidx), clear = T, width= 75);
  
  count_mark <- 0;
  
  for (f in files) {
    
    if (.on.public.web){
      
      print_mes <- paste(basename(f),"import done!");    
      write.table(print_mes,file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
      
    } else {
      #pb$tick();
    }    
    
    filen <- match(f, files)
    filenums <- c(filenums, filen)
    filenams <- c(filenams, f)
    ## issue #214: define backend based on file format.
    msdata <- mzR::openMSfile(f,backend = NULL)
    .instrumentInfo <- c(.instrumentInfo, list(mzR::instrumentInfo(msdata)))
    fullhd <- mzR::header(msdata)
    spidx <- seq_len(nrow(fullhd))
    
    ## Don't read the individual spectra, just define the names of
    ## the spectra.
    fullhdorder <- c(fullhdorder,
                     MSnbase:::formatFileSpectrumNames(fileIds=filen,
                                                       spectrumIds=seq_along(spidx),
                                                       nSpectra=length(spidx),
                                                       nFiles=length(files)))
    ## Extract all Spectrum info from the header and put it into the featureData
    fdData <- fullhd[spidx, , drop = FALSE]
    ## rename totIonCurrent and peaksCount, as detailed in
    ## https://github.com/lgatto/MSnbase/issues/105#issuecomment-229503816
    names(fdData) <- sub("peaksCount", "originalPeaksCount", names(fdData))
    ## Add also:
    ## o fileIdx -> links to fileNames property
    ## o spIdx -> the index of the spectrum in the file.
    fdData <- cbind(fileIdx = rep(filen, nrow(fdData)),
                    spIdx = spidx,
                    smoothed = rep(as.logical(smoothed.), nrow(fdData)),
                    fdData, stringsAsFactors = FALSE)
    if (MSnbase:::isCdfFile(f)) {
      ## Add the polarity columns if missing in netCDF
      if (!any(colnames(fdData) == "polarity"))
        fdData <- cbind(fdData, polarity = rep(as.integer(NA),
                                               nrow(fdData)))
    }
    ## Order the fdData by acquisitionNum to force use of acquisitionNum
    
    ## as unique ID for the spectrum (issue #103). That way we can use
    ## the spIdx (is the index of the spectrum within the file) for
    ## subsetting and extracting.
    if (!all(sort(fdData$acquisitionNum) == fdData$acquisitionNum))
      warning(paste("Unexpected acquisition number order detected.",
                    "Please contact the maintainers or open an issue",
                    "on https://github.com/lgatto/MSnbase.",
                    sep = "\n")) ## see issue #160
    fdData <- fdData[order(fdData$acquisitionNum), ]
    featureDataList <- c(featureDataList, list(fdData))
    ## Fix for #151; would be nice if we could remove that at some point.
    gc()
    mzR::close(msdata)
    rm(msdata)
  }
  ## new in version 1.9.8
  lockEnvironment(assaydata, bindings = TRUE)
  .cacheEnv <- MSnbase:::setCacheEnv(list("assaydata" = assaydata,
                                          "hd" = NULL),
                                     level = 0,
                                     lock = TRUE)
  
  ## Create 'MSnProcess' object
  process <- new("MSnProcess",
                 processing = paste0("Data loaded [", date(), "]"),
                 files = files,
                 smoothed = NA)
  ## Create 'fdata' and 'pdata' objects
  if (is.null(pdata)) {
    .pd <- data.frame(sampleNames = basename(files))
    rownames(.pd) <- .pd$sampleNames
    pdata <- new("AnnotatedDataFrame",
                 data = .pd)
  }
  ## If we've got a featureDataList, use it
  if (length(featureDataList) > 0) {
    fdata <- do.call(rbind, featureDataList)
    fdata <- cbind(fdata, spectrum = 1:nrow(fdata),
                   stringsAsFactors = FALSE)
    ## Setting rownames on the data.frame not on the AnnotatedDataFrame;
    ## did get strange errors otherwise.
    rownames(fdata) <- fullhdorder
    ## Re-order them
    fdata <- fdata[base::sort(fullhdorder), ]
    fdata <- new("AnnotatedDataFrame", data = fdata)
    ## Re-order the features.
    ## fdata <- fdata[ls(assaydata), ]
  } else fdata <- new("AnnotatedDataFrame")
  
  ## expriment data slot
  if (length(.instrumentInfo) > 1) {
    cmp <- length(unique(sapply(.instrumentInfo, "[[", 1)))
    if (cmp > 1)
      message("According to the instrument information in the files,\n",
              "the data has been acquired on different instruments!")
    for (nm in names(.instrumentInfo[[1]]))
      .instrumentInfo[[1]][[nm]] <- sapply(.instrumentInfo, "[[", nm)
  }
  expdata <- new("MIAPE",
                 instrumentManufacturer = .instrumentInfo[[1]]$manufacturer,
                 instrumentModel = .instrumentInfo[[1]]$model,
                 ionSource = .instrumentInfo[[1]]$ionisation,
                 analyser = as.character(.instrumentInfo[[1]]$analyzer),
                 detectorType = .instrumentInfo[[1]]$detector)
  ## Create ProcessingStep if needed.
  ## Create the OnDiskMSnExp object.
  res <- new("OnDiskMSnExp",
             assayData = assaydata,
             phenoData = pdata,
             featureData = fdata,
             processingData = process,
             experimentData = expdata,
             .cache  =  .cacheEnv)
  if (!is.null(msLevel.)) {
    msLevel. <- as.integer(msLevel.)
    res <- filterMsLevel(res, msLevel.)
  }
  if (any(!is.na(centroided.))) {
    if (length(centroided.) == 1) {
      centroided(res) <- centroided.
    } else {
      for (i in seq_along(centroided.))
        centroided(res, msLevel. = i) <- centroided.[i]
    }
  }
  
  if (.on.public.web){
    write.table("Raw file initialized Successfully!",file="metaboanalyst_spec_proc.txt",append = T,row.names = F,col.names = F, quote = F, eol = "\n");
  }
  
  return(res)
}


