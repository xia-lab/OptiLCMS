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
      print(paste0("The number of CPU cores to be used is set to ", num_cores, "."))
      
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
    
    print("Successfully imported raw MS data!")
    
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
#' @param foldername Character, input the file path to the folder containing
#' the raw MS spectra to be processed.
#' @param mode Character, the data input mode. Default is "onDisk" to avoid memory crash. "inMemory" will
#' absorb data into the memory.
#' @param plotSettings List, plotting parameters produced by SetPlotParam Function. "plot.opts" can be added through this
#' function for samples numbers for plotting. Defalut is "default", "all" will apply all samples for plotting and may cause
#' memory crash, especially for large sample dataset.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, Jasmine Chong \email{jasmine.chong@mail.mcgill.ca},
#' Mai Yamamoto \email{yamamoto.mai@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import MSnbase
#' @import BiocParallel
#' @import parallel

ImportRawMSData <-
  function(foldername,
           mode = "onDisk",
           ncores = 4,
           plotSettings,
           running.controller = NULL) {
    foldername <- paste0(fullUserPath, foldername)
    
    
    if (.on.public.web) {
      load_msnbase()
      write.table(
        21.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
      
      print_mes <-
        "Step 2/6: Start to import the spectrum! \nThis step will take a short time..."
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
      
    }
    
    
    #Build Running plan for data import - Indentify the controller
    if (is.null(running.controller)) {
      c2 <- T
      c3 <- T
    } else {
      c2 <-
        running.controller[["others_1"]][["c2"]] # used to control data import
      c3 <-
        running.controller[["others_1"]][["c3"]] # used to control plotting option
    }
    
    if (!dir.exists(foldername)) {
      foldername <- "/home/glassfish/projects/MetaboDemoRawData/upload"
    }
    
    start.time <- Sys.time()
    
    msg.vec <<- vector(mode = "character")
    
    msg <- c("The uploaded files are raw MS spectra.")
    
    # the "upload" folder should contain two subfolders (groups, i.e. Healthy vs. Disease)
    # each subfolder must contain samples (.mzML/.CDF/.mzXML files)
    
    files <-
      dir(
        foldername,
        pattern = ".mzML|.mzml|.cdf|.mzXML|.mzxml|.mzData|.CDF",
        recursive = T,
        full.name = TRUE
      )
    
    if(length(files) == 0){
      print_mes <-  
        paste0("<font color=\"red\">",
               "\nERROR: No standard MS file found ! Please check the extension of your data.",
               "</font>")
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      stop();
    } else if(length(files) < 3) {
      print_mes <-
        paste0("<font color=\"red\">",
               "\nERROR: At least 3 samples should be provided.",
               "</font>")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
    }
    
    
    count_total_sample <<- length(files)
    
    count_current_sample <<- 0
    
    toRemove = vector()
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
    msg <-
      c(msg, paste("A total of ", length(files), "samples were found."))
    
    sclass <- gsub("^\\.$", "sample", dirname(files))
    
    scomp <-
      strsplit(substr(sclass, 1, min(nchar(sclass))), "", fixed = TRUE)
    
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
        unique(sclass)[1] == "upload" & length(unique(sclass)) == 1) {
      sclass <- rep("Unknown", length(sclass))
    }
    # some sanity check before proceeds
    sclass <- as.factor(sclass)
    
    #if(length(levels(sclass))<2){
    #  AddErrMsg("You must provide classes labels (at least two classes)!");
    #  return(0);
    #}
    
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
    
    if (!.on.public.web) {
      cores <- parallel::detectCores()
      
      if (missing(ncores)) {
        num_cores <- ceiling(cores * 2 / 3)
      } else{
        ncores -> num_cores
      }
      
      print(paste0("The number of CPU cores to be used is set to ", num_cores, "."))
      
      if (.Platform$OS.type == "unix") {
        BiocParallel::register(BiocParallel::bpstart(BiocParallel::MulticoreParam(num_cores)))
      } else {
        # for windows
        BiocParallel::register(BiocParallel::bpstart(BiocParallel::SnowParam(num_cores)))
      }
      
    } else {
      num_cores <- 2
      
      print(paste0("The number of CPU cores to be used is set to ", num_cores, "."))
      BiocParallel::register(BiocParallel::bpstart(BiocParallel::MulticoreParam(num_cores)))
      
    }
    
    
    if (c2) {
      raw_data <-
        suppressMessages(read.MSdata(
          files = files,
          pdata = new("NAnnotatedDataFrame", pd),
          mode = mode,
          msLevel = 1
        ))
      
      cache.save(raw_data, funpartnm = "raw_data_samples_c2")
      
      marker_record("raw_data_samples_c2")
      
      
    } else {
      raw_data <- cache.read ("raw_data_samples", "c2")
      
      marker_record("raw_data_samples_c2")
      
    }
    
    
    if (.on.public.web) {
      write.table(
        22.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
    }
    
    if (c3) {
      if (plotSettings$Plot == TRUE) {
        if (is.null(plotSettings$plot.opts)) {
          plot.opts <- "default"
          
        } else {
          plot.opts <- plotSettings$plot.opts
          
        }
        
        if (plot.opts == "default") {
          #subset raw_data to first 50 samples
          print(
            "To reduce memory usage BPIS and TICS plots will be created using only 10 samples per group."
          )
          
          grp_nms <- names(table(pd$sample_group))
          files <- NA
          
          for (i in 1:length(grp_nms)) {
            numb2ext <- min(table(pd$sample_group)[i], 10)
            filt_df <- pd[pd$sample_group == grp_nms[i], ]
            files.inx <- sample(nrow(filt_df), numb2ext)
            sel.samples <- filt_df$sample_name[files.inx]
            files <- c(files, which(pd$sample_name %in% sel.samples))
          }
          
          raw_data_filt <- filterFile(raw_data, file = na.omit(files))
          
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
            print("ImportRawMSData function aborted!")
            return(0)
          }
        }
        
        print("Plotting BPIS and TICS.")
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
      
      marker_record("raw_data_samples_c3")
      
      
    }
    
    
    print("Successfully imported raw MS data!")
    
    if (.on.public.web) {
      write.table(
        24.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
      
      print_mes <-
        paste0(
          "Step 2/6: Successfully imported raw MS data! (",
          Sys.time(),
          ") \nGoing to the next step..."
        )
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
    }
    
    
    return(raw_data)
  }


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

#' Plot EIC
#' @description This functionn creates an extracted ion chromatogram (EIC) for a specific
#' m/z and retention time. This is used for quality-control of raw m/s data.
#' @param raw_data The object created using the ImportRawMSData function,
#' containing the raw MS data.
#' @param rt_mn Numeric, specify the minimum bound of the retention time range.
#' @param rt_mx Numeric, specify the maximum bound of the retention time range.
#' @param mz_mn Numeric, specify the minimum bound of the m/z range.
#' @param mz_mx Numeric, specify the maximum bound of the m/z range.
#' @param aggreg Character, if "sum", creates a total ion chromatogram.
#' If "max", creates a base peak chromatogram. By default it is set
#' to "sum".
#' @param format Character, input the format of the image to create.
#' @param dpi Numeric, input the dpi of the image to create.
#' @param width Numeric, input the width of the image to create.
#' @export

PlotEIC <-
  function(raw_data,
           rt_mn,
           rt_mx,
           mz_mn,
           mz_mx,
           aggreg = "sum",
           format = "png",
           dpi = 72,
           width = 9) {
    filt_data <- filterRt(raw_data, rt = c(rt_mn, rt_mx))
    filt_mz <- filterMz(filt_data, mz = c(mz_mn, mz_mx))
    mz_slice <- chromatogram(filt_mz, aggregationFun = aggreg)
    
    groupNum <- nlevels(groupInfo)
    
    if (groupNum > 9) {
      col.fun <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      group_colors <- col.fun(groupNum)
    } else{
      group_colors <-
        paste0(RColorBrewer::brewer.pal(9, "Set1")[1:groupNum], "60")
    }
    
    names(group_colors) <- unique(raw_data$sample_group)
    
    eic_name <- paste("EIC_", dpi, ".", format, sep = "")
    
    
    Cairo::Cairo(
      file = eic_name,
      unit = "in",
      dpi = dpi,
      width = width,
      height = width * 5 / 9,
      type = format,
      bg = "white"
    )
    
    plot(mz_slice, col = group_colors[raw_data$sample_group])
    legend(
      "topright",
      legend = unique(raw_data$sample_group),
      pch = 15,
      col = group_colors
    )
    
    dev.off()
    
    
    print("EIC created!")
    return(1)
  }


PlotXIC <-
  function(featureNum,
           sample_labeled,
           Group_labeled,
           format,
           dpi,
           width,
           height) {
    # Parameters initialization
    if (missing(sample_labeled)) {
      sample_labeled <- FALSE
    }
    if (missing(Group_labeled)) {
      Group_labeled <- TRUE
    }
    if (missing(format)) {
      format <- "png"
    }
    if (missing(dpi)) {
      dpi <- 72
    }
    if (missing(width)) {
      width <- 6
    }
    if (missing(height)) {
      height <- width * 1.05
    }
    
    # Load data results
    load("mSet.rda");
    require("MSnbase");
    
    raw_data <- mSet[["onDiskData"]]
    
    raw_data@featureData$retentionTime <-
      unlist(mSet[["msFeatureData"]][["adjustedRT"]])
    
    # Get groups information
    groupsInfo <- raw_data@phenoData@data[["sample_group"]]
    
    group_levels <- levels(groupsInfo)
    
    groupsInfo <-
      sapply(
        group_levels,
        FUN = function(x) {
          which(x == groupsInfo)
        }
      )
    
    samples_names <- raw_data@phenoData@data[["sample_name"]]
    
    
    # Get current feature information
    peak_idx_current <-
      mSet[["FeatureGroupTable"]]@listData[["peakidx"]][[featureNum]]
    
    #peak_table <- mSet[["msFeatureData"]][["chromPeaks"]][peak_idx_current,];
    peak_table <- mSet[["xcmsSet"]]@peaks[peak_idx_current, ]
    
    peak_table <- peakTableSUM(peak_table)
    
    
    rtrange <-
      c(mSet[["FeatureGroupTable"]]@listData[["rtmin"]][featureNum] - 10,
        mSet[["FeatureGroupTable"]]@listData[["rtmax"]][featureNum] + 10)
    
    mzrange <-
      c(mSet[["FeatureGroupTable"]]@listData[["mzmin"]][featureNum] - 0.2,
        mSet[["FeatureGroupTable"]]@listData[["mzmax"]][featureNum] + 0.2)
    
    RawData <- filterMz(filterRt(raw_data, rtrange), mzrange)
    
    title <-
      paste0(round(mSet[["FeatureGroupTable"]]@listData[["mzmed"]][[featureNum]], 4),
             "mz@",
             round(mSet[["FeatureGroupTable"]]@listData[["rtmed"]][[featureNum]], 2),
             "s")
    
    
    # Get feature details as a list in MSdb
    MSdb0 <- MSdb <- IntoLists <- IntoListg <- list()    
    
    for (i in group_levels) {
      
      IntoListg[[i]] <- 0
      
      #Different groups
      if (class(groupsInfo)[1] == "list") {
        array_sample <- groupsInfo[[i]]
        
      } else {
        array_sample <- groupsInfo[, i]
        
      }
      nc <- 0
      
      for (x in array_sample) {
        #Different samples per group
        sample_idx <- which(x == peak_table[, "sample"])
        
        #print(paste0("sample_idx: ",x, sample_idx))
        for (j in sample_idx) {
          #Multiple features in corresponding sample  - extract a little bit more
          mzr <-
            c(peak_table[j, 2] - 0.001, peak_table[j, 3] + 0.001)
          
          rtr <- c(peak_table[j, 5] - 0.1, peak_table[j, 6] + 0.1)
          
          IntoListg[[i]] <- IntoListg[[i]] +  peak_table[j, 7]
          
          IntoLists[[i]][[samples_names[x]]] <- peak_table[j, 7]
          
          
          MSdb0[[samples_names[x]]] <-
            chromatogram(
              RawData,
              mz = mzr,
              rt = rtr,
              aggregationFun = "max"
            )[[x]]
          
          nc <- nc + 1
          
        }
      }
      # Do the average on the same group
      if (nc != 0) {
        IntoListg[[i]] <- IntoListg[[i]] / nc
        
      } else {
        IntoListg[[i]] <- 0.0
        
      }
      if (!identical(MSdb0, list())) {
        MSdb[[i]] <- MSdb0
        
        MSdb0 <- list()
        
      }
      
    }
    
    ## PLotting sample EIC begin -
    
    res0 <- lapply(
      MSdb,
      FUN = function(x) {
        res_sample <- data.frame()
        
        if (length(x) != 0) {
          for (j in 1:length(x)) {
            res_sample <- rbind(res_sample,
                                data.frame(x[[j]]@rtime, x[[j]]@intensity, rep(names(x)[j], length(x[[j]]@rtime))))
          }
        }
        
        res_sample[is.na(res_sample)] <- 0
        
        return(res_sample)
      }
    )
    
    res <- data.frame()
    
    for (w in 1:length(res0)) {
      ncouts <- nrow(res0[[w]])
      
      resds <-
        cbind(res0[[w]], rep(names(res0)[w], nrow(res0[[w]])), rep(NA, ncouts))
      
      
      for (ww in 1:length(unique(resds[, 3]))) {
        cn <- unique(resds[, 3])[ww]
        
        cint <- IntoLists[[w]][cn == names(IntoLists[[w]])]
        
        resds[resds[, 3] == cn, ][which(resds[resds[, 3] == cn, 2] == max(resds[resds[, 3] == cn, 2])), 5] <-
          formatC(cint[[1]], format = "e", digits = 2)
        
      }
      
      res <- rbind(res, resds)
      
    }
    
    colnames(res) <-
      c("RT", "Intensity", "Samples", "Groups", "Labels")
    
    peak_width <- max(res$RT) - min(res$RT)
    
    require(ggplot2)
    require(ggrepel)
    
    
    Cairo::Cairo(
      file = paste0("EIC_", title, "_sample_", dpi, ".", format),
      unit = "in",
      dpi = dpi,
      width = width,
      height = height,
      type = format,
      bg = "white"
    )
    
    s_image <-
      ggplot(res, aes(x = RT, y = Intensity, color = Samples)) + #geom_line() #+
      stat_smooth(
        geom = 'area',
        method = "loess",
        se = F,
        span = 0.4,
        size = 0.35,
        formula = "y ~ x",
        alpha = 1 / 4,
        aes(fill = Samples)
      ) +
      theme_bw() +
      ylim(0, NA) +
      xlim(min(res$RT) - peak_width * 0.05 , max(res$RT) + peak_width * 0.35) +
      theme(
        legend.position = c(0.87, 0.85),
        plot.title = element_text(
          size = 12,
          face = "bold",
          hjust = 0.5
        ),
        legend.title = element_text(size = 8.5),
        legend.text = element_text(size = 8),
        legend.key.size = unit(3, "mm")
      ) +
      ggtitle(title)
    
    if (sample_labeled) {
      s_image <-
        s_image + geom_text_repel(aes(y = Intensity * 0.2, label = Labels),
                                  force = 1.5,
                                  show.legend = FALSE)
    }
    
    print(s_image)
    
    dev.off()
    ## PLotting sample EIC finished -
    
    
    ## PLotting group EIC begin --
    res <- NULL
    
    res <- lapply(MSdb, featureSUM, frtr = rtrange)
    
    if (any(sapply(res, is.null))) {
      res <- res[-which(sapply(res, is.null))]
    }
    
    res_data <- data.frame()
    
    for (k in 1:length(res)) {
      ncout <- nrow(res[[k]])
      
      resd <-
        cbind(res[[k]], rep(names(res[k]), ncout), rep(NA, ncout))
      
      resd[which(resd$Inten_mean == max(resd$Inten_mean)), 4] <-
        formatC(IntoListg[[k]], format = "e", digits = 2)
      
      res_data <- rbind(res_data, resd)
      
    }
    
    colnames(res_data) <- c("RT", "Intensity", "Groups", "Labels")
    
    rownames(res_data) <- NULL
    
    peak_width <- max(res_data$RT) - min(res_data$RT)
    
    Cairo::Cairo(
      file = paste0("EIC_", title, "_group_", dpi, ".", format),
      unit = "in",
      dpi = dpi,
      width = width,
      height = height,
      type = format,
      bg = "white"
    )
    
    g_image <-
      ggplot(res_data, aes(x = RT, y = Intensity, color = Groups)) + #geom_line() +
      stat_smooth(
        geom = 'area',
        method = "loess",
        se = F,
        span = 0.4,
        size = 0.35,
        formula = "y ~ x",
        alpha = 1 / 4,
        aes(fill = Groups)
      ) +
      theme_bw() +
      ylim(0, NA) +
      xlim(min(res_data$RT) - peak_width * 0.75 ,
           max(res_data$RT) + peak_width * 0.75) +
      theme(
        legend.position = c(0.85, 0.88),
        plot.title = element_text(
          size = 12,
          face = "bold",
          hjust = 0.5
        ),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.key.size = unit(4, "mm")
      ) +
      ggtitle(title)
    
    if (Group_labeled) {
      g_image <-
        g_image + geom_text_repel(aes(y = Intensity * 0.2, label = Labels),
                                  force = 1.5,
                                  show.legend = FALSE)
      
    }
    
    print(g_image)
    
    dev.off()
    ## PLotting group EIC finished --
    
    return(title)
  }


featureSUM <- function(MS_group, frtr) {
  # sanity check
  if (identical(MS_group, list())) {
    rts <- quantile(frtr)
    
    res <- data.frame(rts[1], 0)
    
    res[2,] <- c(rts[2], 0.5)
    
    res[3,] <- c(rts[3], 1)
    
    res[4,] <- c(rts[4], 0.5)
    
    res[5,] <- c(rts[5], 0)
    
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
  df <- df[order(df$rt), ]
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
    res <- res[-which(is.nan(res$RT_mean) | is.nan(res$Inten_mean)), ]
  }
  
  # manually add 2 empty points before and after the range
  res[nrow(res) + 1,] <- c(min(res$RT_mean) - binsize, 0)
  
  res[nrow(res) + 1,] <- c(min(res$RT_mean) - 2 * binsize, 0)
  
  res[nrow(res) + 1,] <- c(max(res$RT_mean) + binsize, 0)
  
  res[nrow(res) + 1,] <- c(max(res$RT_mean) + 2 * binsize, 0)
  
  
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
#' @import stats
#' @import MSnbase
#' @import BiocParallel
#' @import ggplot2

PerformPeakProfiling <-
  function(rawData,
           Params,
           plotSettings,
           ncore,
           running.controller = NULL) {
    #Build Running plan for data import - Indentify the controller
    function.name <- "peak_profiling"
    
    
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
    
    param <- updateRawSpectraParam (Params)
    
    ### Setting the different parallel method for linux or windows
    if (.on.public.web) {
      write.table(
        25,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
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
    
    if (.on.public.web) {
      print_mes <-
        paste0(ceiling(total_threads),
               " CPU Threads will be used for peak profiling !")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
      
    } else {
      print(paste0(
        ceiling(total_threads),
        " CPU Threads will be used for peak profiling !"
      ))
    }
    
    #   ---------===========----- I. Peak picking -----===========------------
    
    if (c1) {
      if (.on.public.web) {
        print_mes <-
          "Step 3/6: Started peak picking! This step will take some time..."
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
      } else {
        print("Step 1/3: Started peak picking! This step will take some time...")
      }
      
      mSet <-
        tryCatch(
          PerformPeakPicking(rawData, param = param),
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
      write.table(
        50.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
      if (class(mSet)[1] == "simpleError") {
        #write.table(0.0, file = paste0(fullUserPath, "log_progress.txt"),row.names = F,col.names = F);
        print_mes <-
          paste0("<font color=\"red\">",
                 "\nERROR:",
                 mSet$message,
                 "</font>")
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        stop("EXCEPTION POINT CODE: PU1")
      }
      
      
      print_mes <-
        paste0("Step 3/6: Peak picking finished ! (", Sys.time(), ")")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
    }
    
    #   --------===========----- II. Peak alignment -----===========------------
    if (c2) {
      if (.on.public.web) {
        print_mes <-
          paste("Step 4/6: Started peak alignment! This step is running...")
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        
      } else {
        print("Step 2/3: Started peak alignment! This step is running and will take some time...")
        
      }
      
      mSet <-
        tryCatch(
          PerformPeakAlignment(mSet, param),
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
      write.table(
        73.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
      if (class(mSet)[1] == "simpleError") {
        #write.table(0.0, file = paste0(fullUserPath, "log_progress.txt"),row.names = F,col.names = F);
        print_mes <-
          paste0("<font color=\"red\">",
                 "\nERROR:",
                 mSet$message,
                 "</font>")
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        stop("EXCEPTION POINT CODE: PU2")
        
      }
      
      print_mes <-
        paste0("Step 4/6: Peak alignment finished ! (", Sys.time(), ")")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
    }
    
    
    #   --------===========----- III. Peak filling -----===========------------
    if (c3) {
      if (.on.public.web) {
        print_mes <-
          paste("Step 5/6: Started peak filling! This step may take some time...")
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        
      } else {
        print("Step 3/3: Started peak filling! This step may take some time...")
      }
      
      mSet <-
        tryCatch(
          PerformPeakFiling (mSet, param),
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
    
    
    if (.on.public.web) {
      write.table(
        88.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
      
      if (class(mSet)[1] == "simpleError") {
        #write.table(0.0, file = paste0(fullUserPath, "log_progress.txt"),row.names = F,col.names = F);
        print_mes <-
          paste0("<font color=\"red\">",
                 "\nERROR:",
                 mSet$message,
                 "</font>")
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        stop("EXCEPTION POINT CODE: PU3")
        
      }
      
      print_mes <-
        paste0("Step 5/6: Peak filing finished ! (", Sys.time(), ")")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
      
      print_mes <-
        paste0("Peak Profiling finished successfully !")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
    } else {
      print("Peak picking finished successfully !")
      
    }
    
    save(mSet, file = "mSet.rda")
    
    
    if (.on.public.web) {
      print_mes <- paste0("Begin to plotting figures...")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
      write.table(
        89.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
    }
    
    
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
    
    if (.on.public.web) {
      write.table(
        90.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
    }
    
    return(mSet)
    
  }

#' Set generic Plotting Parameters
#' @description This function sets the generic Plotting Parameters for different functions
#' @param Plot Logical, if true, the function will plot internal figures for different functions.
#' @param labels Logical, if true, the labels in the plot will be added.
#' @param format Numeric, input the format of the image to create.
#' @param dpi Numeric, input the dpi of the image to create.
#' @param width Numeric, input the width of the image to create.
#' @param ... Other specific parameters for specific function. Please set them according to the corresponding function.
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca}, and Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export

SetPlotParam <-
  function(Plot = F,
           labels = TRUE,
           format = "png",
           dpi = 72,
           width = 9,
           ...) {
    return(list(
      Plot = Plot,
      labels = labels,
      format = format,
      dpi = dpi,
      width = width,
      ...
    ))
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
    if (.on.public.web) {
      dyn.load(.getDynLoadPath())
      
      write.table(
        91.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
      print_mes <-
        paste0("Step 6/6: Starting Peak Annotation...")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
    }
    
    
    if (ncore > 1) {
      print("Only single core mode is supported now. Parallel will be supported later !")
      ncore <- 1
    }
    
    if (.on.public.web) {
      load_progress()
      
      load_graph()
      
      load_RBGL()
      
    }
    
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
            require(rmpi, character.only = TRUE) && !is.null(ncore)) {
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
          if (try(require(snow, character.only = TRUE, quietly = TRUE))
          ) {
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
        cat("Run cleanParallel after processing to remove the spawned slave processes!\n")
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
      
      if (.on.public.web) {
        write.table(
          92.0,
          file = paste0(fullUserPath, "log_progress.txt"),
          row.names = F,
          col.names = F
        )
        
      }
      
      ## 2. Group peaks according to their retention time into pseudospectra-groups-----
      
      intval <- "maxo"
      perfwhm <- annotaParam$perf.whm
      sigma <- 6
      
      sample    <- mSet$AnnotateObject$sample
      
      pspectra  <- list()
      
      psSamples <- NA
      
      
      if (.on.public.web) {
        print_mes <- paste0("Start grouping after retention time.")
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        
      } else {
        print("Start grouping after retention time.")
        
        
      }
      
      
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
              maxo       <-
                maxo[-iint, , drop = FALSE]
              #set itensities of peaks to NA, due to not to be found in the next cycle
              peakmat  <- peakmat[-iint, , drop = FALSE]
              
            }
          }
          psSamples <- rep(sample, length(pspectra))
        }
        
        mSet$AnnotateObject$pspectra  <- pspectra
        
        mSet$AnnotateObject$psSamples <- psSamples
        
        
        if (.on.public.web) {
          print_mes <-
            paste("Created",
                  length(mSet$AnnotateObject$pspectra),
                  "pseudospectra.")
          
          write.table(
            print_mes,
            file = "metaboanalyst_spec_proc.txt",
            append = T,
            row.names = F,
            col.names = F,
            quote = F,
            eol = "\n"
          )
          
          
        } else {
          message (paste(
            "Created",
            length(mSet$AnnotateObject$pspectra),
            "pseudospectra."
          ))
          
        }
        
        
      }
      
      if (.on.public.web) {
        write.table(
          93.0,
          file = paste0(fullUserPath, "log_progress.txt"),
          row.names = F,
          col.names = F
        )
        
      }
      
      
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
        
        if (.on.public.web) {
          print_mes <- paste0("Generating peak matrix...")
          
          write.table(
            print_mes,
            file = "metaboanalyst_spec_proc.txt",
            append = T,
            row.names = F,
            col.names = F,
            quote = F,
            eol = "\n"
          )
          
          
        } else {
          message("Generating peak matrix...")
          
          
        }
        
        
        mint     <-
          groupval(mSet$xcmsSet, value = intval)[, index, drop = FALSE]
        
        imz <- mSet$AnnotateObject$groupInfo[, "mz", drop = FALSE]
        
        irt <- mSet$AnnotateObject$groupInfo[, "rt", drop = FALSE]
        
      } else{
        ##one sample case
        message("Generating peak matrix...")
        
        imz  <- mSet$AnnotateObject$groupInfo[, "mz", drop = FALSE]
        
        irt  <- mSet$AnnotateObject$groupInfo[, "rt", drop = FALSE]
        
        mint <-
          mSet$AnnotateObject$groupInfo[, intval, drop = FALSE]
        
      }
      
      isotope   <- vector("list", length(imz))
      
      
      isomatrix <- matrix(ncol = 5, nrow = 0)
      
      colnames(isomatrix) <-
        c("mpeak", "isopeak", "iso", "charge", "intrinsic")
      
      if (.on.public.web) {
        print_mes <- paste0("Run isotope peak annotation..")
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        
      } else {
        message("Run isotope peak annotation")
        
        
      }
      
      
      
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
      
      
      if (.on.public.web) {
        print_mes <- paste0("Found isotopes:", cnt)
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        
      } else {
        cat("Found isotopes:", cnt, "\n")
        
        
      }
      
      mSet$AnnotateObject$isotopes <- isotope
      
      
      if (.on.public.web) {
        write.table(
          96.0,
          file = paste0(fullUserPath, "log_progress.txt"),
          row.names = F,
          col.names = F
        )
        
      }
      
      
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
      
      
      if (.on.public.web) {
        print_mes <- paste0("Start grouping after correlation...")
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        
      } else {
        print("Start grouping after correlation.")
        
      }
      
      
      #Data is not preprocessed with groupFWHM
      if (npspectra < 1) {
        cat(
          "Data was not preprocessed with groupFWHM, creating one pseudospectrum with all peaks.\n"
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
      if (.on.public.web) {
        print_mes <-
          paste(
            "mSet has now",
            length(mSet$AnnotateObject$pspectra),
            "groups, instead of",
            cnt,
            "!"
          )
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
      } else {
        message (paste(
          "mSet has now",
          length(mSet$AnnotateObject$pspectra),
          "groups, instead of",
          cnt,
          "!"
        ))
        
      }
      
      
      if (.on.public.web) {
        write.table(
          97.0,
          file = paste0(fullUserPath, "log_progress.txt"),
          row.names = F,
          col.names = F
        )
        
      }
      
      
      
      ## 5. Annotate adducts (and fragments) -----
      mSet$AnnotateObject$ruleset <- NULL
      mSet <-
        findAdducts (
          mSet,
          polarity = annotaParam$polarity,
          mzabs = annotaParam$mz.abs.add,
          maxcharge = annotaParam$max.charge
        )
      
      if (.on.public.web) {
        write.table(
          98.0,
          file = paste0(fullUserPath, "log_progress.txt"),
          row.names = F,
          col.names = F
        )
        
      }
      
      
      
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
      camera_output <- camera_output[, -c(7:endGroup)]
      
      saveRDS(camera_output,
              paste0(fullUserPath, "annotated_peaklist.rds"))
      fast.write(camera_output,
                 paste0(fullUserPath, "annotated_peaklist.csv"))
      
      if (.on.public.web) {
        print_mes <-
          paste0(
            "Step 6/6: Successfully performed peak annotation! (",
            Sys.time(),
            ") \nGoing to the final step..."
          )
        
        write.table(
          print_mes,
          file = "metaboanalyst_spec_proc.txt",
          append = T,
          row.names = F,
          col.names = F,
          quote = F,
          eol = "\n"
        )
        
        
      } else {
        message("Successfully performed peak annotation!")
        
      }
      
      
      if (.on.public.web) {
        write.table(
          99.0,
          file = paste0(fullUserPath, "log_progress.txt"),
          row.names = F,
          col.names = F
        )
        
      }
      
      
      
      ## 0. Final Output -----
      
      if (running.as.plan) {
        cache.save(mSet, paste0(function.name, "_0"))
        
        marker_record(paste0(function.name, "_0"))
        
      }
      
      
    } else {
      mSet <- cache.read(function.name, "0")
      
      marker_record("raw_data_samples_c2")
      
      
    }
    
    
    
    if (.on.public.web) {
      write.table(
        99.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
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
          camera_ma[(camera_ma$isotopes == "" & camera_ma$adduct == ""), ]
        unique_feats <-
          unique(rbind(camera_isotopes, camera_adducts, camera_feats))
        
      } else{
        # negative polarity
        if (filtIso == TRUE) {
          camera_isotopes <- camera_ma[grepl("\\[M\\]-", camera_ma$isotopes), ]
        } else{
          camera_isotopes <- camera_ma[(camera_ma$isotopes != ""), ]
        }
        camera_adducts <-
          camera_ma[grepl("\\[M-H\\]-", camera_ma$adduct), ]
        camera_feats <-
          camera_ma[(camera_ma$isotopes == "" & camera_ma$adduct == ""), ]
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
          camera_isotopes <- camera_ma[grepl("\\[M\\]-", camera_ma$isotopes), ]
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
    
    mzs_rd <- paste0(round(unique_feats[, 1], 4),"@",round(unique_feats[, 4], 2))
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
    
    
    if (.on.public.web) {
      write.table(
        100.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
      print_mes <-
        paste0("Everything has been finished Successfully ! (",
               Sys.time(),
               ")")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
    }
    
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
  
  ann_data <- ann_data[, c(1, 4, ncol(ann_data) - 1, ncol(ann_data) - 2)]
  
  ann_data[, 1] <- round(ann_data[, 1], 4)
  
  ann_data[, 2] <- round(ann_data[, 2], 2)
  
  print("finished the feature summary")
  
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

SanityCheckRawSpectra <- function() {
  mesObj <- rep(list(NA), 8)
  
  names(mesObj) <-
    c(
      "groupInfo",
      "sampleCount",
      "dataformats",
      "QCinclude",
      "centroid",
      "junkfile",
      "decision",
      "messageSum"
    )
  #---------1. Group Information--------/
  mesObj$groupInfo <- dir("upload")
  
  mes1 <- paste0(mesObj$groupInfo, " Uploaded successfully !")
  
  #---------2. Samples Count------------/
  mesObj$sampleCount <-
    sapply(
      mesObj$groupInfo,
      FUN = function(i) {
        length(list.files(paste0("upload/", i)))
      }
    )
  
  mes2 <-
    paste0(unname(mesObj$sampleCount),
           " samples in ",
           names(mesObj$sampleCount),
           " founded !")
  
  #---------3. Data Format--------------/
  require(tools)
  
  mesObj$dataformats <-
    sapply(
      mesObj$groupInfo,
      FUN = function(i) {
        file_ext(list.files(paste0("upload/", i)))
      }
    )
  
  mes3 <-
    paste0("Formats: ", unique(unname(unlist(
      mesObj$dataformats
    ))), " founded ! ")
  
  #---------4. QCs included-------------/
  mesObj$QCinclude <- any(names(mesObj$sampleCount) == "QC")
  
  if (mesObj$QCinclude) {
    mes4 <- paste0("QC samples founded !")
  } else {
    mes4 <-
      paste0("QC samples not founded, will randomly selected some samples for optimization !")
  }
  
  #---------5. Data CCentroid-----------/
  require(mzR)
  
  mesObj$centroid <- verifyCentroid(paste0("upload"))
  
  if (is.null(mesObj$centroid)) {
    mes5 <- paste0("All data have been centroided !")
  } else {
    mes5 <-
      paste0(
        "Files: ",
        mesObj$centroid,
        " have not been centroided ! They will not be included for processing  !"
      )
  }
  
  
  #---------6. Junk File----------------/
  mesObj$junkfile <-
    sapply(
      unique(unname(unlist(
        mesObj$dataformats
      ))),
      FUN = function(x) {
        !any (x == "mzML",
              x == "mzml",
              x == "mzXML",
              x == "mzxml",
              x == "mzData",
              x == "mzdata")
      }
    )
  if (!mesObj$junkfile) {
    mes6 <- paste0("No junk files founded in side all zips !")
  } else {
    mes6 <-
      paste0("Some files except MS data have been found, please remove them first !")
  }
  
  #---------7. Decision ----------------/
  if (any(mesObj$junkfile) |
      !is.null(mesObj$centroid) | any(mesObj$sampleCount < 2)) {
    mesObj$decision <- F
  } else {
    mesObj$decision <- T
  }
  
  
  #---------8. messagesum ----------------/
  mesObj$messageSum <- c(mes1, mes2, mes3, mes4, mes5, mes6)
  
  return(mesObj)
}

## plot single TIC - useful for web only
plotSingleTIC <- function(filename, imagename) {
  load_msnbase()
  
  
  load("raw_data_filt.rda")
  
  load("tics.rda")
  
  
  file_order <-
    which(filename == raw_data_filt@phenoData@data[["sample_name"]])
  
  
  Cairo::Cairo(
    file = paste0(imagename, ".png"),
    unit = "in",
    dpi = 72,
    width = 7,
    height = 5,
    type = "png",
    bg = "white"
  )
  
  
  plot(tics[[file_order]], col = "#0080FF", main = filename)
  
  
  dev.off()
  
  
}

## Verify the parameters changed or not
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

## Functions for other plotting

PlotSpectraInsensityStistics <-
  function(mSet,
           imgName,
           format = "png",
           dpi = 72,
           width = NA) {
    sample_idx <- mSet[["onDiskData"]]@phenoData@data[["sample_group"]]
    
    sample_num <-
      mSet[["onDiskData"]]@phenoData@data[["sample_name"]]
    
    
    Cairo::Cairo(
      file = imgName,
      unit = "in",
      dpi = dpi,
      width = width,
      height = length(sample_num) * 0.65,
      type = format,
      bg = "white"
    )
    
    if (length(unique(sample_idx)) > 9) {
      col.fun <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      group_colors <- col.fun(length(unique(sample_idx)))
      
    } else{
      group_colors <-
        paste0(RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(sample_idx))], "60")
    }
    
    ints <-
      split(log2(mSet[["msFeatureData"]][["chromPeaks"]][, "into"]),
            f = mSet[["msFeatureData"]][["chromPeaks"]][, "sample"])
    
    names(ints) <- as.character(sample_num)
    
    group_colors <-
      sapply(
        seq(length(levels(sample_idx))),
        FUN = function(x) {
          rep(group_colors[x], length(sample_idx[sample_idx == levels(sample_idx)[x]]))
        }
      )
    
    op <- par(mar = c(3.5, 10, 4, 1.5), xaxt = "s")
    
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
    
    dev.off()
    
  }


PlotSpectraPCA <-
  function(mSet,
           imgName,
           format = "png",
           dpi = 72,
           width = NA) {
    Cairo::Cairo(
      file = imgName,
      unit = "in",
      dpi = dpi,
      width = width,
      height = width * 0.80,
      type = format,
      bg = "white"
    )
    
    sample_idx <-
      mSet[["onDiskData"]]@phenoData@data[["sample_group"]]
    
    
    feature_value <-
      .feature_values(
        pks = mSet[["msFeatureData"]][["chromPeaks"]],
        fts = mSet[["FeatureGroupTable"]],
        method = "medret",
        value = "into",
        intensity = "into",
        colnames = mSet[["onDiskData"]]@phenoData@data[["sample_name"]],
        missing = NA
      )
    
    
    pca_feats <- log2(feature_value)
    
    
    if (nrow(feature_value) < 2) {
      print_mes_tmp <-
        paste(
          "\nERROR: No enough peaks detected, please adjust your parameters or use other Peak/Alignment method"
        )
      
      print_mes <-
        paste0("<font color=\"red\">",
               "\nERROR:",
               print_mes_tmp,
               "</font>")
      
      write.table(
        print_mes,
        file = "metaboanalyst_spec_proc.txt",
        append = T,
        row.names = F,
        col.names = F,
        quote = F,
        eol = "\n"
      )
      
      write.table(
        65.0,
        file = paste0(fullUserPath, "log_progress.txt"),
        row.names = F,
        col.names = F
      )
      
      
      dev.off()
      
      return(NULL)
      
    }
    
    df0 <- na.omit(pca_feats)
    
    df1 <- df0[is.finite(rowSums(df0)), ]
    df <- t(df1)
    
    
    mSet_pca <- prcomp(df, center = TRUE, scale = T)
    sum.pca <- summary(mSet_pca)
    var.pca <- sum.pca$importance[2, ] # variance explained by each PCA
    
    xlabel <- paste("PC1", "(", round(100 * var.pca[1], 1), "%)")
    
    ylabel <- paste("PC2", "(", round(100 * var.pca[2], 1), "%)")
    
    
    # using ggplot2
    df <- as.data.frame(mSet_pca$x)
    df$group <- sample_idx
    
    if (.on.public.web) {
      require("ggrepel")
      
      
      if (nrow(df) < 30) {
        if (length(unique(sample_idx)) > 9) {
          col.fun <-
            grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
          p <-
            ggplot2::ggplot(df, aes(
              x = PC1,
              y = PC2,
              color = group,
              label = row.names(df)
            )) +
            geom_text_repel(force = 1.5) + geom_point(size = 5) + fill = col.fun(length(unique(sample_idx))) + theme(axis.text =
                                                                                                                       element_text(size = 12))
          
        } else{
          p <-
            ggplot2::ggplot(df, aes(
              x = PC1,
              y = PC2,
              color = group,
              label = row.names(df)
            )) +
            geom_text_repel(force = 1.5) + geom_point(size = 5) + scale_color_brewer(palette =
                                                                                       "Set1") + theme(axis.text = element_text(size = 12))
        }
        
      } else {
        if (length(unique(sample_idx)) > 9) {
          
          p <-
            ggplot2::ggplot(df, aes(
              x = PC1,
              y = PC2,
              color = group
            )) + geom_point(size = 5)
          
        } else{
          p <-
            ggplot2::ggplot(df, aes(
              x = PC1,
              y = PC2,
              color = group
            )) + geom_point(size = 5) + scale_color_brewer(palette = "Set1")
          
        }
        
      }
      
    } else{
      if (length(unique(sample_idx)) > 9) {
        
        p <-
          ggplot2::ggplot(df, aes(x = PC1, y = PC2, color = group)) +
          geom_point(size = 3) + scale_color_brewer(palette = "Set1")
        
      } else{
        p <- ggplot2::ggplot(df, aes(x = PC1, y = PC2, color = group)) +
          geom_point(size = 3) + scale_color_brewer(palette = "Set1")
      }
      
    }
    
    p <-
      p + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(axis.title = element_text(size =
                                                                                       12))
    print(p)
    
    dev.off()
    
  }



PlotSpectraRTadj <-
  function(mSet,
           imgName,
           format = "png",
           dpi = 72,
           width = NA) {
    sample_idx <- mSet[["onDiskData"]]@phenoData@data[["sample_group"]]
    
    
    Cairo::Cairo(
      file = imgName,
      unit = "in",
      dpi = dpi,
      width = width,
      height = width * 0.75,
      type = format,
      bg = "white"
    )
    
    if (length(unique(sample_idx)) > 9) {
      col.fun <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      group_colors <- col.fun(length(unique(sample_idx)))
    } else{
      group_colors <-
        paste0(RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(sample_idx))], "60")
    }
    
    
    names(group_colors) <- unique(sample_idx)
    
    
    ## Extract RT information
    
    rt.set <-
      list(mSet[["xcmsSet"]]@rt$raw, unlist(mSet[["xcmsSet"]]@rt$corrected))
    
    diffRt <- rt.set[[2]] - rt.set[[1]]
    
    diffRt <- split(diffRt, fromFile(mSet$onDiskData))
    
    
    xRt <- mSet[["msFeatureData"]][["adjustedRT"]]
    
    col = group_colors[sample_idx]
    
    lty = 1
    
    lwd = 1
    
    col <- rep(col, length(diffRt))
    
    lty <- rep(lty, length(diffRt))
    
    lwd <- rep(lwd, length(diffRt))
    
    
    ylim <- range(diffRt, na.rm = TRUE)
    
    plot(
      3,
      3,
      pch = NA,
      xlim = range(xRt, na.rm = TRUE),
      ylim = ylim,
      xlab = "RT_adj",
      ylab = "RT_diff"
    )
    
    
    for (i in 1:length(diffRt)) {
      points(
        x = xRt[[i]],
        y = diffRt[[i]],
        col = col[i],
        lty = lty[i],
        type = "l",
        lwd = lwd[i]
      )
    }
    
    rawRt <-
      split(mSet[["xcmsSet"]]@rt$raw, fromFile(mSet$onDiskData))
    
    adjRt <- xRt
    
    
    ####
    peaks_0 <- mSet[["msFeatureData"]][["chromPeaks"]]
    subs <-
      seq_along(mSet[["onDiskData"]]@phenoData@data[["sample_name"]])
    ####
    pkGroup <- mSet[["msFeatureData"]][["pkGrpMat_Raw"]]
    
    ####
    
    rawRt <- rawRt[subs]
    adjRt <- adjRt[subs]
    ## Have to "adjust" these:
    pkGroupAdj <- pkGroup
    if (!is.null(pkGroup)) {
      for (i in 1:ncol(pkGroup)) {
        pkGroupAdj[, i] <-
          .applyRtAdjustment(pkGroup[, i], rawRt[[i]], adjRt[[i]])
      }
      
      diffRt <- pkGroupAdj - pkGroup
      
      xRt <- pkGroupAdj
      
      ## Loop through the rows and plot points - ordered by diffRt!
      for (i in 1:nrow(xRt)) {
        idx <- order(diffRt[i,])
        points(
          x = xRt[i,][idx],
          diffRt[i,][idx],
          col = "#00000080",
          type = "b",
          pch = 16,
          lty = 3
        )
      }
      
    }
    
    legend(
      "topright",
      legend = unique(sample_idx),
      pch = 15,
      col = group_colors
    )
    
    
    dev.off()
    
  }

PlotSpectraBPIadj <-
  function(mSet,
           imgName,
           format = "png",
           dpi = 72,
           width = NA) {
    Cairo::Cairo(
      file = imgName,
      unit = "in",
      dpi = dpi,
      width = width,
      height = width * 0.618,
      type = format,
      bg = "white"
    )
    
    sample_idx <-
      mSet[["onDiskData"]]@phenoData@data[["sample_group"]]
    
    
    if (length(unique(sample_idx)) > 9) {
      col.fun <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      group_colors <- col.fun(length(unique(sample_idx)))
      
    } else{
      group_colors <-
        paste0(RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(sample_idx))], "60")
    }
    group_colors2 <- group_colors
    
    names(group_colors2) <- unique(sample_idx)
    
    
    group_colors <-
      sapply(
        seq(length(levels(sample_idx))),
        FUN = function(x) {
          rep(group_colors[x], length(sample_idx[sample_idx == levels(sample_idx)[x]]))
        }
      )
    
    object_od <- mSet$onDiskData
    adj_rt <- unlist(mSet$msFeatureData$adjustedRT)
    
    object_od <- selectFeatureData(
      object_od,
      fcol = c(
        "fileIdx",
        "spIdx",
        "seqNum",
        "acquisitionNum",
        "msLevel",
        "polarity",
        "retentionTime",
        "precursorScanNum"
      )
    )
    object_od@featureData$retentionTime <- adj_rt
    
    res <- MSnbase::chromatogram(
      object_od,
      aggregationFun = "max",
      missing = NA_real_,
      msLevel = 1,
      BPPARAM = bpparam()
    )
    
    plot(res, col = as.character(unlist(group_colors)))
    
    legend(
      "topright",
      legend = unique(sample_idx),
      pch = 15,
      col = group_colors2
    )
    
    dev.off()
  }

plotMSfeature <- function(FeatureNM,
                          dpi = 72,
                          format = "png") {
  require(ggplot2)
  
  require(RColorBrewer)
  
  
  peakdata <- readRDS("annotated_peaklist.rds")
  
  peakdata1 <-
    peakdata[, c(-1:-6,
                 -ncol(peakdata),
                 -ncol(peakdata) + 1,
                 -ncol(peakdata) + 2)]
  
  
  peakdata1[is.na(peakdata1)] <- 0
  
  
  group_info <- read.csv("metaboanalyst_input.csv")[c(1), -1]
  
  data_table <-
    as.data.frame(t(rbind(
      round(as.numeric(peakdata1[FeatureNM, ]), 1), as.character(as.matrix(group_info))
    )))
  
  data_table[, 1] <- as.numeric(data_table[, 1])
  
  colnames(data_table) <- c("value", "Group")
  
  rownames(data_table) <- NULL
  
  title = paste0(round(peakdata[FeatureNM, 1], 4), "mz@", round(peakdata[FeatureNM, 4], 2), "s")
  
  
  Cairo::Cairo(
    file = paste0(title, ".png"),
    unit = "in",
    dpi = dpi,
    width = 6,
    height = 6.18,
    type = format,
    bg = "white"
  )
  
  p1 <-
    ggplot(data_table, aes(
      x = Group,
      y = log2(value + 1),
      fill = Group
    )) + # geom_violin(trim = T,draw_quantiles = T) +
    stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
    geom_boxplot(
      size = 0.35,
      width = 0.5,
      fill = "white",
      outlier.fill = "white",
      outlier.color = "white"
    ) +
    geom_jitter(aes(fill = Group),
                width = 0.2,
                shape = 21,
                size = 2.5) +
    scale_color_manual(
      values = c(
        "black",
        "black",
        "black",
        "black",
        "black",
        "black",
        "black",
        "black",
        "black",
        "black",
        "black",
        "black",
        "black"
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        family = "Arial",
        size = 12,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        family = "Arial",
        size = 12,
        face = "plain"
      ),
      axis.title.y = element_text(
        family = "Arial",
        size = 16,
        face = "plain"
      ),
      axis.title.x = element_text(
        family = "Arial",
        size = 16,
        face = "plain"
      ),
      plot.title = element_text(
        family = "Arial",
        size = 16,
        face = "bold",
        hjust = 0.5
      )
    ) +
    ylab(expression(log[2] ~ intensity)) + xlab("Groups") + ggtitle(title)
  
  print(p1);
  
  dev.off();
  
  return(paste0(title, ".png"));
  
}
