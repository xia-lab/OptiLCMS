#' @title PerformDataInspect
#' @description This functions provide a path for users to visually inspect their raw data before the data
#' trimming so as to remove the dirty or significantly uneluted peaks.
#' @param datapath Character, the path of the raw MS data files (.mzXML, .CDF and .mzML)
#' for the visual and intuitive data inspectation or the file folder (if only a folder path provided, the first file will
#' be inspected).
#' @param rt.range Numerics, a congregation of two values to define the lower and upper RT range (seconds) for
#' users to inspect. This is an optional parameter, if absent, will display the MS of the whole RT range.
#' @param mz.range Numerics, a congregation of two values to define the lower and upper mz range for
#' users to inspect. This is an optional parameter, if absent, will display the MS of the whole mz range.
#' @param dimension Character, the dimension for sample to display, including '2D' or '3D'. The default is '3D'.
#' @param res Numeric, the resolution for data inspectation. The larger the value, the higher the resolution.
#' The default value is 100. This value is usually clearly enough and also give consideration to the speed.
#' @export
#' @import MSnbase
#' @importFrom Cairo Cairo
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
#' @examples 
#' library(OptiLCMS)
#' ## Get raw spectra files
#' DataFiles <- dir(system.file("mzData", package = "mtbls2"), full.names = TRUE,
#'                  recursive = TRUE)[c(10:12, 14:16)]
#' PerformDataInspect(DataFiles[1])

PerformDataInspect <-
  function(datapath = NULL,
           rt.range = c(0,0),
           mz.range = c(0,0),
           dimension = "3D",
           res = 100) {
    
    if(.on.public.web){
      
      fullUserPath <- getwd();
      
      if (datapath == "null" | is.null(datapath)) {
        if (.on.public.web & dir.exists("/upload/QC")) {
          datapath <- "/upload/QC"
          datapath <- paste0(fullUserPath, datapath)
        } else if (.on.public.web) {
          datapath <- paste0("/upload/", list.files("upload", recursive = TRUE)[1])
          datapath <- paste0(fullUserPath, datapath)
        } else {
          MessageOutput("Local Inspectation !\n")
        }
      } else {
        files <-
          list.files(paste0(getwd(), "/upload/"),
                     full.names = TRUE,
                     recursive = TRUE)
        
        
        if (isEmpty(files)) {
          # Handle the example issue - show other files
          files <-
            list.files(
              "/home/glassfish/projects/MetaboDemoRawData/upload/",
              full.names = TRUE,
              recursive = TRUE
            )
          datapath =  files[datapath == basename(files)]
          
        } else {
          # Handle the regular data insepctation
          datapath =  files[datapath == basename(files)]
        }
        
      }
      
      if (basename(datapath) == "NA") {
        # Handle the example issue - default showing
        datapath <-
          "/home/glassfish/projects/MetaboDemoRawData/upload/QC"
      }
    }
    
    if (!grepl(pattern = c("*.mzXML"), basename(datapath)) &
        !grepl(pattern = c("*.mzxml"), basename(datapath)) &
        !grepl(pattern = c("*.mzData"), basename(datapath)) &
        !grepl(pattern = c("*.mzdata"), basename(datapath)) &
        !grepl(pattern = c("*.mzML"), basename(datapath)) &
        !grepl(pattern = c("*.mzml"), basename(datapath)) &
        !grepl(pattern = c("*.CDF"), basename(datapath)) &
        !grepl(pattern = c("*.cdf"), basename(datapath))) {
      MessageOutput(paste("First file in ", datapath, " will be inspected !\n"),SuppressWeb = TRUE)
      mzf <- list.files(datapath, recursive = FALSE, full.names = TRUE)[1]
    } else {
      mzf <- datapath
    }
    
    if (grepl(pattern = c("*.CDF"), basename(mzf)) |
        grepl(pattern = c("*.cdf"), basename(mzf))) {
      # centroid_cdf();
      return()
    }
    
    ms <- mzR::openMSfile(mzf)
    hd <- header(ms)
    
    ms1 <- which(hd$msLevel == 1)
    
    if (missing(rt.range) | (rt.range[1] == 0 & rt.range[2] == 0)) {
      rtsel <-
        hd$retentionTime[ms1] > min(hd$retentionTime) &
        hd$retentionTime[ms1] < max(hd$retentionTime)
      
      rt.extension <- F
      MessageOutput(paste(
        "RT range is:",
        min(hd$retentionTime),
        "and",
        max(hd$retentionTime),
        "seconds !\n"
      ), SuppressWeb = TRUE)
      
    } else{
      if (rt.range[2] < rt.range[1]) {
        a1 <- rt.range[1]
        rt.range[1] <- rt.range[2]
        rt.range[2] <- a1
      } else if (rt.range[2] == rt.range[1]) {
        rt.range[2] <- rt.range[2] + 1
        rt.range[1] <- rt.range[1] - 1
      }
      MessageOutput(paste("RT range is:", rt.range[1], "and", rt.range[2], "seconds !\n"),SuppressWeb = TRUE)
      
      rtsel <-
        hd$retentionTime[ms1] > rt.range[1] &
        hd$retentionTime[ms1] < rt.range[2]
      
      rt.min <- min(hd$retentionTime[ms1])
      rt.max <- max(hd$retentionTime[ms1])
      
      if (rt.range[1] < rt.min | rt.range[2] > rt.max) {
        rt.extension <- T
      } else {
        rt.extension <- F
      }
    }
    if (missing(mz.range) | (mz.range[1] == 0 & mz.range[2] == 0)) {
      min.mz <- min(hd$lowMZ)
      max.mz <- max(hd$highMZ)
      
      if (min.mz == 0 & max.mz == 0 | min.mz == max.mz) {
        MessageOutput(
          "mz.range information is missing in your data file. mz between 100 and 1200 will be shown here !\n"
        )
        min.mz <- 100
        max.mz <- 1200
        
      } else if (is.infinite(min.mz) |
                 is.infinite(max.mz) | min.mz == -1 | max.mz == -1) {
        MessageOutput(
          "mz.range information is missing in your data file. mz between 50 and 2000 will be shown here !\n"
        )
        min.mz <- 50
        max.mz <- 2000
      }
      
      MessageOutput(paste("MZ range is:", min.mz, "and", max.mz, "Thomson !\n"),SuppressWeb = TRUE)
      
      res.mz <- (max.mz - min.mz) / res
      M <- MSmap(ms,
                 ms1[rtsel],
                 min.mz,
                 max.mz,
                 res.mz,
                 hd,
                 zeroIsNA = TRUE)
      
    } else{
      if (mz.range[2] < mz.range[1]) {
        a1 <- mz.range[1]
        mz.range[1] <- mz.range[2]
        mz.range[2] <- a1
      } else if (mz.range[2] == mz.range[1]) {
        mz.range[2] <- mz.range[2] + 0.01
        mz.range[1] <- mz.range[1] - 0.01
      }
      
      MessageOutput(paste("MZ range is:", mz.range[1], "and", mz.range[2], "Thomson !\n"),SuppressWeb = TRUE)
      
      res.mz <- (mz.range[2] - mz.range[1]) / res
      M <- MSmap(ms,
                 ms1[rtsel],
                 mz.range[1],
                 mz.range[2],
                 res.mz,
                 hd,
                 zeroIsNA = TRUE)
    }
    
    if (rt.extension) {
      if (min(M@rt) > rt.range[1]) {
        M@rt <- c(rt.range[1], M@rt)
        M@map <- rbind(rep(NA, dim(M@map)[2]), M@map)
        M@ms <- c(M@ms, 1)
      }
      
      if (max(M@rt) < rt.range[2]) {
        M@rt <- c(M@rt, rt.range[2])
        M@map <- rbind(M@map, rep(NA, dim(M@map)[2]))
        M@ms <- c(M@ms, 1)
      }
    }
    
    if (missing(dimension)) {
      dimension = "3D"
    }
    
    FN <- tools::file_path_sans_ext(basename(mzf))
    
    filename <-
      paste0(
        FN,
        "_mz_",
        mz.range[1],
        "_",
        mz.range[2],
        "_RT_",
        rt.range[1],
        "_",
        rt.range[2],
        dimension,
        ".png"
      )

    if(.on.public.web){
      Cairo::Cairo(
        file = filename,
        unit = "in",
        dpi = 72,
        width = 7,
        height = 7,
        type = "png",
        bg = "white"
      )
    }

    if (dimension == "3D") {
      print(plot.MS_3D(M))
    } else {
      print(plot(M, aspect = 1, allTicks = FALSE))
    }
    
    if(.on.public.web){
      dev.off()
    }
    
    return(filename)
  }

#' @title SetPlotParam
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
#' @examples
#' SetPlotParam(Plot = TRUE, dpi = 144, width = 12)

SetPlotParam <-
  function(Plot = FALSE,
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

#' @title PlotXIC/EIC
#' @description This functionn creates an extracted ion chromatogram (XIC/EIC) for a specific
#' m/z and retention time. This is used for quality-control of raw m/s data.
#' @param mSet mSet Object. Should contain the spectra processing result.
#' @param featureNum Numeric, Feature number in the feature table.
#' @param sample_labeled Logical, whether to lable the sample name.
#' @param Group_labeled Logical, whether to lable the group name.
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png".
#' @param dpi Numeric, to define the dpi of the figures. Default is 72.
#' @param width Numeric, to define the width of the figure.
#' @param height Numeric, to define the height of the figure.
#' @importFrom Cairo Cairo
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid viewport
#' @export
#' @examples 
#' library(OptiLCMS)
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' PlotXIC(mSet, 1, TRUE, TRUE);

PlotXIC <-
  function(mSet = NULL,
           featureNum,
           sample_labeled,
           Group_labeled,
           format,
           dpi,
           width,
           height,
           sample_filled,
           group_filled) {
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
    
    if(missing(sample_filled)){
      sample_filled <- FALSE;
    }
    
    if(missing(group_filled)){
      group_filled <- TRUE;
    }
    
    # Load data results
    if(.on.public.web){
      load("mSet.rda")
    } else if(is.null(mSet)) {
      stop("mSet is missing!")
    }

    raw_data <- mSet@rawOnDisk;
    raw_data@featureData$retentionTime <-
      unlist(mSet@peakRTcorrection$adjustedRT);
    
    # Get groups information
    groupsInfo <- raw_data@phenoData@data[["sample_group"]]
    group_levels <- levels(as.factor(groupsInfo))
    
    groupsInfo <-
      sapply(
        group_levels,
        FUN = function(x) {
          which(x == groupsInfo)
        }
      )
    
    samples_names <- raw_data@phenoData@data[["sample_name"]]
    
    # Get current feature information
    peak_idx_current <-mSet@peakfilling$FeatureGroupTable@listData$peakidx[[featureNum]]
    
    #peak_table <- mSet[["msFeatureData"]][["chromPeaks"]][peak_idx_current,];
    peak_table <- mSet@peakfilling[["msFeatureData"]][["chromPeaks"]][peak_idx_current,]
    peak_table <- peakTableSUM(peak_table)
    
    rtrange <-
      c(mSet@peakfilling[["FeatureGroupTable"]]@listData[["rtmin"]][featureNum] - 10,
        mSet@peakfilling[["FeatureGroupTable"]]@listData[["rtmax"]][featureNum] + 10)
    
    mzrange <-
      c(mSet@peakfilling[["FeatureGroupTable"]]@listData[["mzmin"]][featureNum] - 0.2,
        mSet@peakfilling[["FeatureGroupTable"]]@listData[["mzmax"]][featureNum] + 0.2)
    
    RawData <- filterMz(filterRt(raw_data, rtrange), mzrange)
    
    title <-
      paste0(round(mSet@peakfilling[["FeatureGroupTable"]]@listData[["mzmed"]][[featureNum]], 4),
             "mz@",
             round(mSet@peakfilling[["FeatureGroupTable"]]@listData[["rtmed"]][[featureNum]], 2),
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
          rtr <- c(peak_table[j, which(colnames(peak_table) == "rtmin")] - 0.1, 
                   peak_table[j, which(colnames(peak_table) == "rtmax")] + 0.1)
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
    
    ## Get Group CV Table
    CVTable <- PeakGroupCV(IntoLists, groupsInfo);
    CVTable <- data.frame(x = CVTable$V1,y=CVTable$V2);
    pCV <- ggplot(data=CVTable, aes(x=x, y=y, fill =x, alpha = 0.5))  + 
      geom_col(width = 1.32/ncol(CVTable)) + 
      theme_bw() + 
      theme(text = element_text(size=10), 
            axis.title.x = element_blank(),
            axis.text.x=element_blank(), 
            legend.position = "none",
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
      labs(y = expression(italic("Coefficient of variation"))) + 
      geom_hline(yintercept=1, linetype="dashed", color = "red");

    ## PLotting sample EIC begin -
    res0 <- lapply(
      MSdb,
      FUN = function(x) {
        res_sample <- data.frame()
        
        if (length(x) != 0) {
          for (j in 1:length(x)) {
            sampleRes <- data.frame(x[[j]]@rtime, x[[j]]@intensity, rep(names(x)[j], length(x[[j]]@rtime)));
            rtdiff <- mean(diff(sampleRes$x..j...rtime));
            insertEmptyScan <- data.frame(c(unname(x[[j]]@rtime[1]) - 4*rtdiff, unname(x[[j]]@rtime[1]) - 3*rtdiff, 
                                            unname(x[[j]]@rtime[1]) - 2*rtdiff, unname(x[[j]]@rtime[1]) - rtdiff,
                                            unname(tail(x[[j]]@rtime,1)) + rtdiff, unname(tail(x[[j]]@rtime,1)) + 2*rtdiff,
                                            unname(tail(x[[j]]@rtime,1)) + 3*rtdiff, unname(tail(x[[j]]@rtime,1)) + 4*rtdiff),
                                            rep(0, 8), 
                                            rep(names(x)[j], 8));
            colnames(sampleRes) <- colnames(insertEmptyScan) <- colnames(sampleRes);
            sampleRes <- rbind(insertEmptyScan, sampleRes)
            res_sample <- rbind(res_sample,sampleRes);
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
        resds[resds[, 3] == cn,][which(resds[resds[, 3] == cn, 2] == max(resds[resds[, 3] == cn, 2])), 5] <-
          formatC(cint[[1]], format = "e", digits = 2)
      }
      
      res <- rbind(res, resds)
    }
    
    colnames(res) <-
      c("RT", "Intensity", "Samples", "Groups", "Labels")
    
    peak_width <- max(res$RT) - min(res$RT);
    
    if(!is.null(mSet@params[["Peak_method"]])){
      if(mSet@params[["Peak_method"]] == "Massifquant"){
        spanValue <- 0.8
      }else {
        spanValue <- 0.6
      }
    } else {
      spanValue <- 0.55
    }
    print(paste0("EIC_", title, "_sample_", dpi, ".", format));
    if (.on.public.web) {
      Cairo::Cairo(
        file = paste0("EIC_", title, "_sample_", dpi, ".", format),
        unit = "in",
        dpi = dpi,
        width = width,
        height = height,
        type = format,
        bg = "white"
      )
    }
    
    if(sample_filled){
      s_image0 <- ggplot(res, aes_string(x = "RT", y = "Intensity", color = "Samples")) +
        stat_smooth(
          geom = 'area',
          method = "loess",
          se = FALSE,
          span = spanValue,
          size = 0.35,
          formula = "y ~ x",
          alpha = 1 / 4,
          aes_string(fill = "Samples")
        )
    } else {
      s_image0 <- ggplot(res, aes_string(x = "RT", y = "Intensity", color = "Samples")) +
        stat_smooth(
          method = "loess",
          se = FALSE,
          span = spanValue,
          size = 0.35,
          formula = "y ~ x",
          alpha = 1 / 4,
          aes_string(color = "Samples"),
          fill = NA
        )
    }
    
    s_image <-s_image0 +
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
        legend.key.size = unit(3, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
      ) +
      ggtitle(title)
    
    if (sample_labeled) {
      s_image <-
        s_image + geom_text_repel(aes_string(y = "Intensity * 0.2", label = "Labels"),
                                  force = 1.5,
                                  show.legend = FALSE)
    }
    
    print(s_image);
    
    if (.on.public.web) {
      dev.off()
    }
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
    
    if (.on.public.web) {
      Cairo::Cairo(
        file = paste0("EIC_", title, "_group_", dpi, ".", format),
        unit = "in",
        dpi = dpi,
        width = width,
        height = height,
        type = format,
        bg = "white"
      )
    }
    
    if(group_filled){
      g_image0 <- ggplot(res_data, aes_string(x = "RT", y = "Intensity", color = "Groups")) + #geom_line() +
        stat_smooth(
          geom = 'area',
          method = "loess",
          se = FALSE,
          span = spanValue,
          size = 0.35,
          formula = "y ~ x",
          alpha = 1 / 4,
          aes_string(fill = "Groups")
        ) 
    } else {
      g_image0 <- ggplot(res_data, aes_string(x = "RT", y = "Intensity", color = "Groups")) + #geom_line() +
        stat_smooth(
          geom = 'area',
          method = "loess",
          se = FALSE,
          span = spanValue,
          size = 0.35,
          formula = "y ~ x",
          alpha = 1 / 4,
          aes_string(color = "Groups"),
          fill = NA
        ) 
    }
    
    g_image <- g_image0 +
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
        legend.key.size = unit(4, "mm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
      ) +
      ggtitle(title)
    
    if (Group_labeled) {
      g_image <-
        g_image + geom_text_repel(aes_string(y = "Intensity * 0.2", label = "Labels"),
                                  force = 1.5,
                                  show.legend = FALSE)
    }
    require(grid);
    print(g_image); print(pCV,vp=viewport(.255, .82, .24, 0.24));
    
    if (.on.public.web) {
      dev.off()
    }
    ## PLotting group EIC finished --
    return(title)
  }


#' @title PlotSpectraInsensityStistics
#' @description This function is used to do the statistics on the spectra intensity
#' @param mSet mSet object, usually generated after the peakannotaion finished here.
#' @param imgName Character, to give the name of BPI figures ploted.
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png".
#' @param dpi Numeric, to define the dpi of the figures. Default is 72.
#' @param width Numeric, to define the width of the figure. Height = width * 0.618. 
#' @export
#' @importFrom Cairo Cairo
#' @examples 
#' library(OptiLCMS)
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' PlotSpectraInsensityStistics(mSet);

PlotSpectraInsensityStistics <-
  function(mSet = NULL,
           imgName,
           format = "png",
           dpi = 72,
           width = NA) {
    
    if(is.null(mSet) & .on.public.web){
      load("mSet.rda");
    } else if(is.null(mSet)) {
      stop("mSet is missing!")
    }
    
    sample_idx <- mSet@rawOnDisk@phenoData@data[["sample_group"]];
    
    sample_num <-
      mSet@rawOnDisk@phenoData@data[["sample_name"]];
    
    if (length(unique(sample_idx)) > 9) {
      col.fun <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      group_colors <- col.fun(length(unique(sample_idx)))
      
    } else{
      group_colors <-
        paste0(RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(sample_idx))], "60")
    }
    
    ints <-
      split(log2(mSet@peakfilling$msFeatureData$chromPeaks[, "into"]),
            f = mSet@peakfilling$msFeatureData$chromPeaks[, "sample"])
    
    names(ints) <- as.character(sample_num)
    
    sample_idx <- as.factor(sample_idx)
    group_colors <-
      sapply(
        seq(length(levels(sample_idx))),
        FUN = function(x) {
          rep(group_colors[x], length(sample_idx[sample_idx == levels(sample_idx)[x]]))
        }
      )
    
    if(!.on.public.web){
      oldpar <- par(no.readonly = TRUE);
      on.exit(par(oldpar));
    }
    
    if (.on.public.web) {
      Cairo::Cairo(
        file = imgName,
        unit = "in",
        dpi = dpi,
        width = width,
        height = length(sample_num) * 0.65,
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
    
    if (.on.public.web) {
      dev.off()
    }
  }


#' @title PlotSpectraPCA
#' @description This function is used to plot the PCA of all spectra
#' @param mSet mSet object, usually generated after the peakannotaion finished here.
#' @param imgName Character, to give the name of BPI figures ploted.
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png".
#' @param dpi Numeric, to define the dpi of the figures. Default is 72.
#' @param width Numeric, to define the width of the figure. Height = width * 0.618. 
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples 
#' library(OptiLCMS)
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' PlotSpectraPCA(mSet);

PlotSpectraPCA <-
  function(mSet = NULL,
           imgName,
           format = "png",
           dpi = 72,
           width = NA) {
    
    if(missing(mSet) || is.null(mSet)){
      load("mSet.rda")
    }
    
    if (.on.public.web) {
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
    
    sample_idx <-
      mSet@rawOnDisk@phenoData@data[["sample_group"]];
    
    # feature_value <-
    #   .feature_values(
    #     pks = mSet@peakfilling$msFeatureData$chromPeaks,
    #     fts = mSet@peakfilling$FeatureGroupTable,
    #     method = "medret",
    #     value = "into",
    #     intensity = "into",
    #     colnames = mSet@rawOnDisk@phenoData@data[["sample_name"]],
    #     missing = NA
    #   );
    
    feature_value0 <- mSet@dataSet[-1,];
    rownames(feature_value0) <- feature_value0[,1];
    feature_value <- feature_value0[,-1];
    feature_value[is.na(feature_value)] <- 0;
    
    int.mat <- as.matrix(feature_value)
    rowNms <- rownames(int.mat);
    colNms <- colnames(int.mat);
    int.mat <- t(apply(int.mat, 1, function(x) .replace.by.lod(as.numeric(x))));
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
      
      if (.on.public.web) {
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
    zlabel <- paste("PC2", "(", round(100 * var.pca[3], 1), "%)");
    
    # using ggplot2
    df <- as.data.frame(mSet_pca$x);
    df$group <- sample_idx;
    
    ## Handle to generate json file for PCA3D online
    if(.on.public.web){
      ## For score plot
      pca3d <- list();
      pca3d$score$axis <- c(xlabel, ylabel, zlabel);
      xyz0 <- df[,c(1:3)];
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
      
      json.obj <- RJSONIO::toJSON(pca3d, .na='null');
      sink("spectra_3d_score.json");
      cat(json.obj);
      sink();
      ## For loading plot
      pca3d <- list();
      
      pca3d$loading$axis <- paste("Loading ", c(1:3), sep="");
      coords0 <- coords <- data.frame(t(signif(mSet_pca$rotation[,1:3], 5)));
      colnames(coords) <- NULL; 
      pca3d$loading$xyz <- coords;
      pca3d$loading$name <- rownames(mSet_pca$rotation);
      pca3d$loading$entrez <- paste0(round(mSet@peakfilling[["FeatureGroupTable"]]@listData$mzmed, 4), 
                                     "@", 
                                     round(mSet@peakfilling[["FeatureGroupTable"]]@listData$rtmed, 2));
      save(coords0, file = "coords0_ori.rda")
      dists <- GetDist3D(coords0);
      pca3d$loading$cols <- GetRGBColorGradient(dists);
      
      pca3d$cls =  df$group;
      
      json.obj <- RJSONIO::toJSON(pca3d, .na='null');
      sink("spectra_3d_loading.json");
      cat(json.obj);
      sink();
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
          ggplot2::ggplot(df, aes_string(x = "PC1",
                                  y = "PC2",
                                  color = "group")) + geom_point(size = 5)

      } else{
        p <-
          ggplot2::ggplot(df, aes_string(x = "PC1",
                                  y = "PC2",
                                  color = "group")) + geom_point(size = 5) + scale_color_brewer(palette = "Set1");
      }
    }

    p <-
      p + xlab(xlabel) + ylab(ylabel) + theme_bw() + theme(axis.title = element_text(size = 12));

    print(p)
    
    if (.on.public.web) {
      dev.off()
    }
  }


#' @title PlotSpectraRTadj
#' @description This function is used to plot the adjustment of retention time of all spectra
#' @param mSet mSet object, usually generated after the peakannotaion finished here.
#' @param imgName Character, to give the name of BPI figures ploted.
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png".
#' @param dpi Numeric, to define the dpi of the figures. Default is 72.
#' @param width Numeric, to define the width of the figure. Height = width * 0.618. 
#' @export
#' @importFrom Cairo Cairo
#' @examples 
#' library(OptiLCMS)
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' PlotSpectraRTadj(mSet);

PlotSpectraRTadj <-
  function(mSet = NULL,
           imgName,
           format = "png",
           dpi = 72,
           width = NA) {
    
    if(is.null(mSet) & .on.public.web){
      load("mSet.rda")
    } else if(is.null(mSet)) {
      stop("mSet is missing!")
    }
    
    sample_idx <- mSet@rawOnDisk@phenoData@data[["sample_group"]];
    
    if (.on.public.web) {
      Cairo::Cairo(
        file = imgName,
        unit = "in",
        dpi = dpi,
        width = width,
        height = width * 0.75,
        type = format,
        bg = "white"
      )
    }
    
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
      list(rtime(mSet@rawOnDisk), unlist(mSet@peakRTcorrection$adjustedRT))
    
    diffRt <- rt.set[[2]] - rt.set[[1]]
    diffRt <- split(diffRt, fromFile(mSet@rawOnDisk))
    xRt <- mSet@peakRTcorrection$adjustedRT
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
      xlab = "RT_Adjustment",
      ylab = "RT_Difference"
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
      split(rtime(mSet@rawOnDisk), fromFile(mSet@rawOnDisk))
    
    adjRt <- xRt
    
    ####
    peaks_0 <- mSet@peakfilling$msFeatureData$chromPeaks;
    subs <-
      seq_along(mSet@rawOnDisk@phenoData@data[["sample_name"]]);
    ####
    pkGroup <- mSet@peakRTcorrection[["pkGrpMat_Raw"]];
    ####
    
    rawRt <- rawRt[subs]
    adjRt <- adjRt[subs]
    ## Have to "adjust" these for peakgroup only:
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
        idx <- order(diffRt[i, ])
        points(
          x = xRt[i, ][idx],
          diffRt[i, ][idx],
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
    
    if (.on.public.web) {
      dev.off()
    }
  }


#' @title PlotSpectraBPIadj
#' @description This function is used to plot the adjust BPI (Base Peak Ion)
#' @param mSet mSet object, usually generated after the peakannotaion finished here.
#' @param imgName Character, to give the name of BPI figures ploted.
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png".
#' @param dpi Numeric, to define the dpi of the figures. Default is 72.
#' @param width Numeric, to define the width of the figure. Height = width * 0.618. 
#' @export
#' @importFrom Cairo Cairo
#' @examples 
#' library(OptiLCMS)
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' PlotSpectraBPIadj(mSet);

PlotSpectraBPIadj <-
  function(mSet = NULL,
           imgName,
           format = "png",
           dpi = 72,
           width = NA) {
    
    if(is.null(mSet) & .on.public.web){
      load("mSet.rda")
    } else if(is.null(mSet)) {
      stop("mSet is missing!")
    }
    
    if (.on.public.web) {
      Cairo::Cairo(
        file = imgName,
        unit = "in",
        dpi = dpi,
        width = width,
        height = width * 0.75,
        type = format,
        bg = "white"
      )
    }
    
    sample_idx <-
      mSet@rawOnDisk@phenoData@data[["sample_group"]]
    sample_idx <- as.factor(sample_idx);
    
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
    
    object_od <- mSet@rawOnDisk
    adj_rt <- unlist(mSet@peakRTcorrection$adjustedRT)
    
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
      msLevel = 1L,
      BPPARAM = bpparam()
    )
    
    plot(res, col = as.character(unlist(group_colors)))
    
    legend(
      "topright",
      legend = unique(sample_idx),
      pch = 15,
      col = group_colors2
    )
    
    if (.on.public.web) {
      dev.off()
    }
  }

#' @title plotMSfeature
#' @description plotMSfeature is used to plot the feature intensity of different groups
#' @param mSet mSet Object, should be processed aby 'PerformPeakProfiling'.
#' @param FeatureNM Numeric, feature number in the feature table.
#' @param dpi Numeric, to define the dpi of the figures. Default is 72. (only works for web version)
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png". (only works for web version)
#' @export
#' @importFrom Cairo Cairo
#' @examples 
#' library(OptiLCMS)
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' plotMSfeature (mSet, 1); # Here is only one group

plotMSfeature <- function(mSet = NULL, 
                          FeatureNM = 1,
                          dpi = 72,
                          format = "png",
                          width = NA) {
  
  if(is.null(mSet)){
    if(.on.public.web){
      load("mSet.rda")
    } else {
      stop("No mSet Object found !")
    }
  }

  peakdata <- mSet@peakAnnotation$camera_output;
  peakdata1 <-
    peakdata[, c(-1:-6,-ncol(peakdata),-ncol(peakdata) + 1,-ncol(peakdata) + 2)]
  
  peakdata1[is.na(peakdata1)] <- 0
  
  group_info <- mSet@dataSet[c(1),-1]
  
  data_table <-
    as.data.frame(t(rbind(
      round(as.numeric(peakdata1[FeatureNM,]), 1), as.character(as.matrix(group_info))
    )))
  
  data_table[, 1] <- as.numeric(data_table[, 1])
  colnames(data_table) <- c("value", "Group")
  rownames(data_table) <- NULL
  
  if(is.na(width)){
    width = 6;
    title = paste0(round(peakdata[FeatureNM, 1], 4), "mz@", round(peakdata[FeatureNM, 4], 2), "s")
  } else {
    title = paste0(round(peakdata[FeatureNM, 1], 4), "mz@", round(peakdata[FeatureNM, 4], 2), "s_", dpi, "_", width)
  }
  
  if (.on.public.web) {
    Cairo::Cairo(
      file = paste0(title, ".", format),
      unit = "in",
      dpi = dpi,
      width = width,
      height = width/6*6.18,
      type = format,
      bg = "white"
    )
  }
    
  p1 <-
    ggplot(data_table, aes_string(
      x = "Group",
      y = "log2(value + 1)",
      fill = "Group"
    )) + # geom_violin(trim = TRUE,draw_quantiles = TRUE) +
    stat_boxplot(geom = "errorbar", width = 0.15, aes(color = "black")) +
    geom_boxplot(
      size = 0.35,
      width = 0.5,
      fill = "white",
      outlier.fill = "white",
      outlier.color = "white"
    ) +
    geom_jitter(aes_string(fill = "Group"),
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
        #family = "Arial",
        size = 12,
        angle = 45,
        hjust = 1
      ),
      axis.text.y = element_text(
        #family = "Arial",
        size = 12,
        face = "plain"
      ),
      axis.title.y = element_text(
        #family = "Arial",
        size = 16,
        face = "plain"
      ),
      axis.title.x = element_text(
        #family = "Arial",
        size = 16,
        face = "plain"
      ),
      plot.title = element_text(
        #family = "Arial",
        size = 16,
        face = "bold",
        hjust = 0.5
      )
    ) +
    ylab(expression(log[2] ~ intensity)) + xlab("Groups") + ggtitle(title)
  
  print(p1)
  
  if (.on.public.web) {
    dev.off()
    return(paste0(title, ".", format))
  }
}

#' @title plotSingleTIC
#' @description plotSingleTIC is used to plot the TIC of a certain spectra
#' @param mSet mSet Object, should be processed by ImportMSData.
#' @param filename Character, to give the filename for the TIC plotting.
#' @param imagename Character, to give the filename of the TIC plotted. (only works for web version)
#' @export
#' @importFrom Cairo Cairo
#' @examples
#' library(OptiLCMS)
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' plotSingleTIC(mSet, "MSpos-Ex2-Col0-48h-Ag-2_1-A,3_01_9829.mzData", 
#'               "MSpos-Ex2-Col0-48h-Ag-2_1-A,3_01_9829.png")

plotSingleTIC <- function(mSet = NULL, filename, imagename, dpi = 72, width = 7, format = "png") {
  
  if(is.null(mSet)){
    if(.on.public.web){
      load("mSet.rda");
    } else {
      stop("No mSet found !")
    }
  }
  
  AllFileNMs <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(mSet@rawOnDisk@processingData@files));
  file_order <-
    which(filename == AllFileNMs)
  
  raw_data_filt <-
    filterFile(mSet@rawOnDisk, file = file_order);
  
  tics <- chromatogram(raw_data_filt, aggregationFun = "sum");

  if (.on.public.web) {
    Cairo::Cairo(
      file = imagename,
      unit = "in",
      dpi = dpi,
      width = width,
      height = width/7*5,
      type = format,
      bg = "white"
    )
  }
  
  plot(tics, col = "#0080FF", main = filename);
  
  if (.on.public.web) {
    dev.off()
  }
}

#' @title plotTICs
#' @description plotTICs is used to plot the TIC of all files
#' @param mSet mSet Object, should be processed by ImportMSData.
#' @param imgName Character, to name the imgName for the TIC plotting.
#' @param dpi Numeric, to define the dpi of the figures. Default is 72. (only works for web version)
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png". (only works for web version)
#' @export
#' @importFrom Cairo Cairo
#' @examples

plotTICs <-function(mSet = NULL,
                   imgName,
                   format = "png",
                   dpi = 72,
                   width = NA){
  #need to extract the plotting part from import function
  if(is.null(mSet) & !.on.public.web){
    load("mSet.rda")
    raw_data_filt <- mSet@rawOnDisk;
  } else if(!.on.public.web & (class(mSet) == "mSet")) {
    raw_data_filt <- mSet@rawOnDisk;
  } else if(.on.public.web) {
    load("raw_data_filt.rda");
    load("tics.rda");
  } else if(is.null(mSet)) {
    stop("mSet is missing!")
  }
  
  samplegroup <- raw_data_filt@phenoData@data$sample_group;
  samplename <- raw_data_filt@phenoData@data$sample_name;
  
  groupInfo <-as.factor(samplegroup);
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
  if (.on.public.web) {
    Cairo::Cairo(
      file = imgName,
      unit = "in",
      dpi = dpi,
      width = width,
      height = width* 0.75,
      type = format,
      bg = "white"
    )
  }
  
  plot(tics, col = group_colors[raw_data_filt$sample_group])
  legend(
    "topright",
    legend = levels(groupInfo),
    pch = 15,
    col = group_colors
  )
  
  if (.on.public.web) {
    dev.off()
  }

}

#' @title plotBPIs
#' @description plotBPIs is used to plot the BPI of all files
#' @param mSet mSet Object, should be processed by ImportMSData.
#' @param filename Character, to give the filename for the TIC plotting.
#' @param dpi Numeric, to define the dpi of the figures. Default is 72. (only works for web version)
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png". (only works for web version)
#' @export
#' @importFrom Cairo Cairo
#' @examples
plotBPIs <-function(mSet = NULL,
                    imgName,
                    format = "png",
                    dpi = 72,
                    width = NA){
  #need to extract the plotting part from import function
  if(is.null(mSet) & !.on.public.web){
    load("mSet.rda")
    raw_data_filt <- mSet@rawOnDisk;
  } else if(!.on.public.web & (class(mSet) == "mSet")) {
    raw_data_filt <- mSet@rawOnDisk;
  } else if(.on.public.web) {
    load("raw_data_filt.rda");
    load("bpis.rda");
  } else if(is.null(mSet)) {
    stop("mSet is missing!")
  }
  
  samplegroup <- raw_data_filt@phenoData@data$sample_group;
  samplename <- raw_data_filt@phenoData@data$sample_name;
  
  groupInfo <-as.factor(samplegroup);
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
  if (.on.public.web) {
    Cairo::Cairo(
      file = imgName,
      unit = "in",
      dpi = dpi,
      width = width,
      height = width* 0.75,
      type = format,
      bg = "white"
    )
  }
  
  plot(bpis, col = group_colors[raw_data_filt$sample_group])
  legend(
    "topright",
    legend = levels(groupInfo),
    pch = 15,
    col = group_colors
  )
  
  if (.on.public.web) {
    dev.off()
  }
  
}


