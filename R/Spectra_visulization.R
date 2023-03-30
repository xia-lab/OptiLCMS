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
#' @return will output a figure for viewing the data structure
#' @import MSnbase
#' @importFrom Cairo Cairo
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)
#' @examples 
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
        if (.on.public.web & dir.exists("upload/QC")) {
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
        hd$retentionTime[ms1] >= min(hd$retentionTime) &
        hd$retentionTime[ms1] <= max(hd$retentionTime)
      
      rt.extension <- FALSE
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
        rt.extension <- TRUE
      } else {
        rt.extension <- FALSE
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

    seed.cols <- c("#fafafa", "#e0e0e0");
    cols <- colorRampPalette(seed.cols)(8)
    seed.cols22 <- c("#ffee00","#ffc400", "#ff0000","#a600ff", "#0006b0","#000363")
    cols2 <- colorRampPalette(seed.cols22)(8)
    
    if (dimension == "3D") {
      print(plot.MS_3D(M))
    } else {
      M@mz <- round(M@mz,3)
      print(plot(M, aspect = 1, allTicks = FALSE, col.regions = c(cols, cols2)))
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
#' @return will return a plotting parameters set
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
#' @param sample_filled Logical, to determine the EIC/XIC is filled or not for sample EIC
#' @param group_filled Logical, to determine the EIC/XIC is filled or not for group EIC
#' @importFrom Cairo Cairo
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid viewport
#' @importFrom MSnbase filterMz
#' @export
#' @return will return a figure of EIC/XIC
#' @examples 
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' #PlotXIC(mSet, 1, TRUE, TRUE);


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
    groupsInfo0 <- groupsInfo <- raw_data@phenoData@data[["sample_group"]]
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
    
    RawData <- MSnbase::filterMz(MSnbase::filterRt(raw_data, rtrange), mzrange)
    
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
      if (is(groupsInfo,"list")) {
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
    y <- NULL;
    CVTable <- PeakGroupCV(IntoLists, groupsInfo);
    
    CVTable <- data.frame(x = CVTable$V1,y=CVTable$V2);
    pCV <- ggplot(data=CVTable, aes(x=x, y=y, fill =x, alpha = 0.5))  + 
      geom_col(width = 1.32/ncol(CVTable)) + 
      theme_bw() + 
      theme(text = element_text(size=8.5), 
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
          for (j in seq_along(x)) {
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
    
    for (w in seq_along(res0)) {
      ncouts <- nrow(res0[[w]])
      
      resds <-
        cbind(res0[[w]], rep(names(res0)[w], nrow(res0[[w]])), rep(NA, ncouts))
      
      for (ww in seq_along(unique(resds[, 3]))) {
        cn <- unique(resds[, 3])[ww]
        cint <- IntoLists[[w]][cn == names(IntoLists[[w]])]
        resds[resds[, 3] == cn,][which(resds[resds[, 3] == cn, 2] == max(resds[resds[, 3] == cn, 2])), 5] <-
          formatC(cint[[1]], format = "e", digits = 2)
      }
      
      res <- rbind(res, resds)
    }
    
    colnames(res) <-
      c("RT", "Intensity", "Samples", "Groups", "Labels");
    peak_width <- max(res$RT) - min(res$RT);
    
    ## Resume for the missing peaks
    minRT <- min(res$RT);
    maxRT <- max(res$RT);
    RTvec <- c(minRT, minRT + (maxRT - minRT)/4, minRT + (maxRT - minRT)/2, minRT + (maxRT - minRT)*3/4, maxRT)
    
    missingIdx <- which(!(samples_names %in% unique(res$Samples)));

    for(m in missingIdx){
      res <- rbind(res, data.frame(RT=RTvec, 
                                   Intensity = rep(0, 5), 
                                   Samples = rep(samples_names[m], 5), 
                                   Groups = rep(groupsInfo0[m], 5), 
                                   Labels = c(NA, NA, 0, NA, NA)))
    }

    
    
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
    
    res_data <- data.frame();
    
    for (k in seq_along(res)) {
      ncout <- nrow(res[[k]]);
      
      resd <-
        cbind(res[[k]], rep(names(res[k]), ncout), rep(NA, ncout))
      
      resd[which(resd$Inten_mean == max(resd$Inten_mean)), 4] <-
        formatC(IntoListg[[names(res[k])]], format = "e", digits = 2)
      
      res_data <- rbind(res_data, resd)
    }
    
    colnames(res_data) <- c("RT", "Intensity", "Groups", "Labels");
    rownames(res_data) <- NULL;
    peak_width <- max(res_data$RT) - min(res_data$RT);
    
    ## ADD the missing groups
    minRT <- maxRT <- RTvec <- missingIdx <- m <- NULL;
    minRT <- min(res_data$RT);
    maxRT <- max(res_data$RT);
    RTvec <- c(minRT, minRT + (maxRT - minRT)/4, minRT + (maxRT - minRT)/2, minRT + (maxRT - minRT)*3/4, maxRT)
    
    missingIdx <- which(!(unique(groupsInfo0) %in% unique(res_data$Groups)));
    
    for(m in missingIdx){
      res_data <- rbind(res_data, data.frame(RT=RTvec,
                                             Intensity = rep(0, 5),
                                             Groups = rep(unique(groupsInfo0)[m], 5), 
                                             Labels = c(NA, NA, 0, NA, NA)))
    }
    
    
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
    #require(grid);
    print(g_image); print(pCV,vp=viewport(.3, .76, .24, 0.24));
    
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
#' @return will return a figure of spectral peak intensity
#' @importFrom Cairo Cairo
#' @importFrom graphics par
#' @importFrom graphics boxplot
#' @importFrom graphics grid
#' @import RColorBrewer
#' @examples 
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
        paste0(RColorBrewer::brewer.pal(9, "Set1")[seq_along(unique(sample_idx))], "60")
    }
    
    ints <-
      split(log2(mSet@peakfilling$msFeatureData$chromPeaks[, "into"]),
            f = mSet@peakfilling$msFeatureData$chromPeaks[, "sample"])
    
    if(mSet@params$BlankSub & any(sample_idx == "BLANK")){
      ints <- ints[names(ints) != 0];
      sample_num <- sample_num[sample_idx != "BLANK"]
    }
    
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
#' @importFrom RJSONIO toJSON
#' @importFrom stats prcomp
#' @import RColorBrewer
#' @return will return a figure of PCA after log tranformation (log2)
#' @export
#' @examples 
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
    
    if(mSet@params$BlankSub & any(sample_idx == "BLANK")){
      sample_idx <- sample_idx[sample_idx != "BLANK"]
    }
    
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
    zlabel <- paste("PC3", "(", round(100 * var.pca[3], 1), "%)");
    
    # using ggplot2
    df <- as.data.frame(mSet_pca$x);
    df$group <- sample_idx;
    
    ## Handle to generate json file for PCA3D online
    if(.on.public.web){
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
      
      # json.obj <- RJSONIO::toJSON(pca3d, .na='null');
      # sink("spectra_3d_score.json");
      # cat(json.obj);
      # sink();
      
      ## For loading plot
      #pca3d <- list();
      
      pca3d$loading$axis <- paste("Loading ", seq_len(3), sep="");
      coords0 <- coords <- data.frame(t(signif(mSet_pca$rotation[,seq_len(3)], 5)));
      colnames(coords) <- NULL; 
      pca3d$loading$xyz <- coords;
      pca3d$loading$name <- rownames(mSet_pca$rotation);
      pca3d$loading$entrez <- paste0(round(mSet@peakfilling[["FeatureGroupTable"]]@listData$mzmed, 4), 
                                     "@", 
                                     round(mSet@peakfilling[["FeatureGroupTable"]]@listData$rtmed, 2));
      
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
#' @return will return a figure of spectral adjustment of retention time
#' @export
#' @import RColorBrewer
#' @importFrom Cairo Cairo
#' @importFrom graphics points
#' @examples 
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
    
    if(mSet@params$BlankSub & any(sample_idx == "BLANK")){
      sample_idx <- sample_idx[sample_idx != "BLANK"];
      
      specdata0 <- mSet@rawOnDisk;
      blk2rms <- which(specdata0@phenoData@data[["sample_group"]] == "BLANK");
      fdnew <- specdata0@featureData
      fdnew@data <- fdnew@data[!(fdnew@data[["fileIdx"]] %in% blk2rms),];
      #fdnew@data[["fileIdx"]] <- fdnew@data[["fileIdx"]] - length(blk2rms)
      ii <- 1;
      for(fii in unique(fdnew@data[["fileIdx"]])){
        fdnew@data[["fileIdx"]][fdnew@data[["fileIdx"]] == fii] <- ii;
        ii <- ii + 1
      }
      
      pdnew <- specdata0@phenoData;
      pdnew@data <- pdnew@data[-blk2rms,];
      
      numfiles <- length(specdata0@experimentData@instrumentModel) - length(blk2rms)

      specdata <- new(
        "OnDiskMSnExp",
        processingData = new("MSnProcess",
                             files = specdata0@processingData@files[-blk2rms]),
        featureData = fdnew,
        phenoData = pdnew,
        experimentData = new("MIAPE",
                             instrumentManufacturer = rep("a",numfiles),
                             instrumentModel = rep("a",numfiles),
                             ionSource = rep("a",numfiles),
                             analyser = rep("a",numfiles),
                             detectorType = rep("a",numfiles)))
      
    } else {
      specdata <-  mSet@rawOnDisk;
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
    
    if (length(unique(sample_idx)) > 9) {
      col.fun <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      group_colors <- col.fun(length(unique(sample_idx)))
    } else{
      group_colors <-
        paste0(RColorBrewer::brewer.pal(9, "Set1")[seq_along(unique(sample_idx))], "60")
    }
    
    names(group_colors) <- unique(sample_idx)
    
    ## Extract RT information
    rt.set <-
      list(rtime(specdata), unlist(mSet@peakRTcorrection$adjustedRT))
    
    diffRt <- rt.set[[2]] - rt.set[[1]]
    diffRt <- split(diffRt, fromFile(specdata))
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
    
    for (i in seq_along(diffRt)) {
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
      split(rtime(specdata), fromFile(specdata))
    
    adjRt <- xRt
    
    ####
    peaks_0 <- mSet@peakfilling$msFeatureData$chromPeaks;
    subs <-
      seq_along(specdata@phenoData@data[["sample_name"]]);
    ####
    pkGroup <- mSet@peakRTcorrection[["pkGrpMat_Raw"]];
    ####
    
    rawRt <- rawRt[subs]
    adjRt <- adjRt[subs]
    ## Have to "adjust" these for peakgroup only:
    pkGroupAdj <- pkGroup
    if (!is.null(pkGroup)) {
      for (i in seq_len(ncol(pkGroup))){
        pkGroupAdj[, i] <-
          .applyRtAdjustment(pkGroup[, i], rawRt[[i]], adjRt[[i]])
      }
      
      diffRt <- pkGroupAdj - pkGroup
      xRt <- pkGroupAdj
      
      ## Loop through the rows and plot points - ordered by diffRt!
      for (i in seq_len(nrow(xRt))){
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
#' @return will return a figure of adjusted BPIs
#' @export
#' @import RColorBrewer
#' @importFrom Cairo Cairo
#' @examples 
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
      mSet@rawOnDisk@phenoData@data[["sample_group"]];
    
    if(mSet@params$BlankSub & any(sample_idx == "BLANK")){
      sample_idx <- sample_idx[sample_idx != "BLANK"];
      
      specdata0 <- mSet@rawOnDisk;
      blk2rms <- which(specdata0@phenoData@data[["sample_group"]] == "BLANK");
      fdnew <- specdata0@featureData
      fdnew@data <- fdnew@data[!(fdnew@data[["fileIdx"]] %in% blk2rms),];
      #fdnew@data[["fileIdx"]] <- fdnew@data[["fileIdx"]] - length(blk2rms)
      ii <- 1;
      for(fii in unique(fdnew@data[["fileIdx"]])){
        fdnew@data[["fileIdx"]][fdnew@data[["fileIdx"]] == fii] <- ii;
        ii <- ii + 1
      }
      
      pdnew <- specdata0@phenoData;
      pdnew@data <- pdnew@data[-blk2rms,];
      
      numfiles <- length(specdata0@experimentData@instrumentModel) - length(blk2rms)
      
      object_od <- new(
        "OnDiskMSnExp",
        processingData = new("MSnProcess",
                             files = specdata0@processingData@files[-blk2rms]),
        featureData = fdnew,
        phenoData = pdnew,
        experimentData = new("MIAPE",
                             instrumentManufacturer = rep("a",numfiles),
                             instrumentModel = rep("a",numfiles),
                             ionSource = rep("a",numfiles),
                             analyser = rep("a",numfiles),
                             detectorType = rep("a",numfiles)))
      
      
    } else {
      object_od <- mSet@rawOnDisk;
    }
    sample_idx <- as.factor(sample_idx);
    
    if (length(unique(sample_idx)) > 9) {
      col.fun <-
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      group_colors <- col.fun(length(unique(sample_idx)))
      
    } else{
      group_colors <-
        paste0(RColorBrewer::brewer.pal(9, "Set1")[seq_along(unique(sample_idx))], "60")
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
    
    #object_od <- mSet@rawOnDisk
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
      BPPARAM = MulticoreParam(4)
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
#' @param width Numeric, width of the figure (default is NA, usually set it as 6~12)
#' @export
#' @return will return a figure of ms stats
#' @import RColorBrewer
#' @importFrom Cairo Cairo
#' @examples 
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
    peakdata[, c((-1):-6,-ncol(peakdata),-ncol(peakdata) + 1,-ncol(peakdata) + 2,
                 -ncol(peakdata) + 3, -ncol(peakdata) + 4, -ncol(peakdata) + 5)]
  
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
#' @param dpi Numeric, dpi of the figure (default is 72, usually set it as 72, 144, 360)
#' @param width Numeric, width of the figure (default is 7, usually set it as 6~12)
#' @param format Character, format of the figure (default is 'png', usually can be 'png', 'pdf','tiff','svg','eps','jpg')
#' @export
#' @return will return a figure of a single TIC
#' @importFrom Cairo Cairo
#' @examples
#' data(mSet);
#' newPath <- dir(system.file("mzData", package = "mtbls2"),
#'                full.names = TRUE, recursive = TRUE)[c(10, 11, 12)]
#' mSet <- updateRawSpectraPath(mSet, newPath);
#' plotSingleTIC(mSet, "MSpos-Ex2-Col0-48h-Ag-2_1-A,3_01_9829.mzData", 
#'               "MSpos-Ex2-Col0-48h-Ag-2_1-A,3_01_9829.png")

plotSingleTIC <- function(mSet = NULL, 
                          filename, 
                          imagename, 
                          dpi = 72, 
                          width = 7, 
                          format = "png") {
  
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
  
  tics <- chromatogram(raw_data_filt, aggregationFun = "sum", MulticoreParam(4));

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
#' @param width Numeric, width of the figure (default is NA, usually set it as 6~12)
#' @export
#' @return will return a figure of TICs
#' @importFrom Cairo Cairo
#' @import RColorBrewer
#' @examples
#' newPath <- dir(system.file("mzData", package = "mtbls2"), 
#' full.names = TRUE, recursive = TRUE)[c(10:12)]
#' data(mSet)
#' mSet <- updateRawSpectraPath(mSet, newPath)
#' plotTICs(mSet)
plotTICs <-function(mSet = NULL,
                   imgName,
                   format = "png",
                   dpi = 72,
                   width = NA){
  tics <- NULL;
  #need to extract the plotting part from import function
  if(is.null(mSet) & !.on.public.web){
    load("mSet.rda")
    raw_data_filt <- mSet@rawOnDisk;
  } else if(!.on.public.web & (is(mSet, "mSet"))) {
    raw_data_filt <- mSet@rawOnDisk;
    tics <- chromatogram(raw_data_filt, aggregationFun = "sum", MulticoreParam(4))
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
      paste0(RColorBrewer::brewer.pal(9, "Set1")[seq_len(groupNum)], "60")
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
#' @param imgName Character, to give the filename for the TIC plotting.
#' @param dpi Numeric, to define the dpi of the figures. Default is 72. (only works for web version)
#' @param format Character, to give the format of BPI figures ploted. Can be "jpeg", "png", "pdf", "svg",
#'  "tiff" or "ps". Default is "png". (only works for web version)
#' @param width Numeric, width of the figure (default is NA, usually set it as 6~12)
#' @export
#' @import RColorBrewer
#' @importFrom Cairo Cairo
#' @return will return a figure of BPIs
#' @examples
#' newPath <- dir(system.file("mzData", package = "mtbls2"), 
#' full.names = TRUE, recursive = TRUE)[c(10:12)]
#' data(mSet)
#' mSet <- updateRawSpectraPath(mSet, newPath)
#' plotBPIs(mSet)
plotBPIs <-function(mSet = NULL,
                    imgName,
                    format = "png",
                    dpi = 72,
                    width = NA){
  bpis <- NULL;
  #need to extract the plotting part from import function
  if(is.null(mSet) & !.on.public.web){
    load("mSet.rda")
    raw_data_filt <- mSet@rawOnDisk;
  } else if(!.on.public.web & (is(mSet, "mSet"))) {
    raw_data_filt <- mSet@rawOnDisk;
    bpis <- chromatogram(raw_data_filt, aggregationFun = "max", MulticoreParam(4))
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
      paste0(RColorBrewer::brewer.pal(9, "Set1")[seq_len(groupNum)], "60")
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




my.json.scatter <- function(filenm, containsLoading=F){
  library(igraph);
  res <- qs::qread("score3d.qs")
  nodes <- vector(mode="list");
  names <- res$name;
  if(ncol(res$xyz) > nrow(res$xyz)){
    orig.pos.xyz <- t(res$xyz);
  }else{
    orig.pos.xyz <- res$xyz;
  }
  ticksX <- pretty(range(orig.pos.xyz[,1]*1.2),10, bound=F);
  ticksY <- pretty(range(orig.pos.xyz[,2]*1.2),10, bound=F);
  ticksZ <- pretty(range(orig.pos.xyz[,3]*1.2),10, bound=F);
  #add two nodes, 1 with all min values and another with all max values. For scaling purposes
  minNode <-  c(min(ticksX), min(ticksY), min(ticksZ));
  maxNode <-  c(max(ticksX), max(ticksY), max(ticksZ));
  # Add the new rows to the data frame
  orig.pos.xyz.mod <- rbind(orig.pos.xyz, minNode)
  orig.pos.xyz.mod <- rbind(orig.pos.xyz.mod, maxNode)
  #scaling
  pos.xyz <- orig.pos.xyz.mod;
  pos.xyz[,1] <- scale_range(orig.pos.xyz.mod[,1], -1,1);
  pos.xyz[,2] <- scale_range(orig.pos.xyz.mod[,2], -1,1);
  pos.xyz[,3] <- scale_range(orig.pos.xyz.mod[,3], -1,1);
  #remove last two rows
  pos.xyz <- pos.xyz[1:(nrow(pos.xyz) - 2), ]
  metadf <- res$facA
  col = vector();
  meta.vec = as.vector(metadf)
  meta.vec.num = as.integer(as.factor(metadf))
  col.s <- res$colors
  for(i in 1:length(meta.vec.num)){
    col[i] = col.s[meta.vec.num[i]];
  }
  legendData <- list(label=unique(meta.vec),color=col.s)
  #for IPCA in multifactor
  if("facB" %in% names(res)){
    meta.vec2 <- res$facB
    metadf <- res$metadata_list;
    shape <- vector();
    meta.vec.num <- as.integer(as.factor(meta.vec2))
    shape.s <- c("circle", "triangle", "square", "diamond")[1:length(unique(meta.vec2))];
    for(i in 1:length(meta.vec.num)){
      shape[i] = shape.s[meta.vec.num[i]];
    }
    legendData2 <- list(label=unique(meta.vec2),shape=shape.s);
  }
  nodeSize = 18;
  for(i in 1:length(names)){
    nodes[[i]] <- list(
      id=names[i],
      label=names[i],
      size=nodeSize,
      meta=meta.vec[i],
      fx = unname(pos.xyz[i,1])*1000,
      fy = unname(pos.xyz[i,2])*1000,
      fz = unname(pos.xyz[i,3])*1000,
      origX = unname(orig.pos.xyz[i,1]),
      origY = unname(orig.pos.xyz[i,2]),
      origZ = unname(orig.pos.xyz[i,3]),
      colorb=col[i]
    );
    if("facB" %in% names(res)){
      nodes[[i]][["meta2"]] <- meta.vec2[i]
      nodes[[i]][["shape"]] <- shape[i]
    }
  }
  ticks <- list(x=ticksX, y=ticksY, z=ticksZ);
  library(RJSONIO)
  if(!containsLoading){
    netData <- list(nodes=nodes,
                    edges="NA",
                    meta=metadf,
                    loading="NA",
                    axis=res$axis,
                    ticks=ticks,
                    metaCol = legendData);
  }else{
    res2 <- qs::qread("loading3d.qs");
    if(ncol(res2$xyz) > nrow(res2$xyz)){
      orig.load.xyz <- t(res2$xyz);
    }else{
      orig.load.xyz <- res2$xyz;
    }
    ticksX <- pretty(range(orig.load.xyz[,1]),10);
    ticksY <- pretty(range(orig.load.xyz[,2]),10);
    ticksZ <- pretty(range(orig.load.xyz[,3]),10);
    #add two nodes, 1 with all min values and another with all max values. For scaling purposes
    minNode <-  c(min(ticksX), min(ticksY), min(ticksZ));
    maxNode <-  c(max(ticksX), max(ticksY), max(ticksZ));
    # Add the new rows to the data frame
    orig.load.xyz.mod <- rbind(orig.load.xyz, minNode)
    orig.load.xyz.mod <- rbind(orig.load.xyz.mod, maxNode)
    load.xyz <- orig.load.xyz.mod;
    load.xyz[,1] <- scale_range(orig.load.xyz.mod[,1], -1,1);
    load.xyz[,2] <- scale_range(orig.load.xyz.mod[,2], -1,1);
    load.xyz[,3] <- scale_range(orig.load.xyz.mod[,3], -1,1);
    #remove last two rows
    load.xyz <- load.xyz[1:(nrow(load.xyz) - 2), ];
    names <- res2$name;
    if("entrez" %in% names(res2)){
      ids <- res2$entrez;
    }else{
      ids <- res2$name;
    }
    colres <- rgba_to_hex_opacity(res2$cols);
    colorb <- colres[[1]];
    opacity_array <- colres[[2]];
    nodes2 <- vector(mode="list");
    for(i in 1:length(names)){
      nodes2[[i]] <- list(
        id=ids[i],
        label=names[i],
        size=24,
        opacity=opacity_array[i],
        fx = unname(load.xyz[i,1])*1000,
        fy = unname(load.xyz[i,2])*1000,
        fz = unname(load.xyz[i,3])*1000,
        origX = unname(orig.load.xyz[i,1]),
        origY = unname(orig.load.xyz[i,2]),
        origZ = unname(orig.load.xyz[i,3]),
        seedArr = "notSeed",
        colorb=colorb[i],
        colorw=colorb[i]
      );
    }
    ticksLoading <- list(x=ticksX, y=ticksY, z=ticksZ);
    netData <- list(omicstype="",
                    nodes=nodes,
                    meta=metadf,
                    loading=nodes2,
                    misc="", ticks=ticks,
                    ticksLoading=ticksLoading,
                    axis=res$axis,
                    axisLoading=res2$axis,
                    metaCol = legendData);
  }
  if("facB" %in% names(res)){
    netData$metaShape <- legendData2;
  }
  rownames(pos.xyz) <- res$name;
  res$pos.xyz <- pos.xyz;
  .set.rdt.set(res);
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(1);
}



# Define a function to convert RGBA to Hex and opacity values
rgba_to_hex_opacity <- function(rgba_array) {
  # Define an empty vector to store the hex color values
  hex_values <- vector("character", length = length(rgba_array))
  # Define an empty vector to store the opacity values
  opacity_values <- vector("numeric", length = length(rgba_array))
  # Loop through each element in the input array
  for (i in seq_along(rgba_array)) {
    rgba <- rgba_array[i]
    rgba <- gsub("rgba\\(", "",rgba);
    rgba <- gsub("\\)", "",rgba);
    # Convert the RGBA value to a list of numeric values
    rgba <- strsplit(rgba, ",")[[1]]
    rgba <- as.numeric(rgba)
    # Convert the RGB values to hexadecimal notation
    hex <- rgb(rgba[1], rgba[2], rgba[3], maxColorValue = 255)
    hex_values[i] <- hex
    # Extract the opacity value from the RGBA string
    opacity_values[i] <- rgba[4]
  }
  # Return a list containing the hex color values and opacity values
  return(list(hex_values = hex_values, opacity_values = opacity_values))
}

scale_range <- function(x, new_min = 0, new_max = 1) {
  range <- pretty(x,10);
  old_min <- min(range);
  old_max <- max(range);
  (x - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
}

ComputeEncasing <- function(filenm, type, names.vec, level=0.95, omics="NA"){
  level <- as.numeric(level)
  names = strsplit(names.vec, "; ")[[1]]
  res <- .get.rdt.set();
  res <- res$pos.xyz
  pos.xyz <- res
  inx = rownames(pos.xyz) %in% names;
  coords = as.matrix(pos.xyz[inx,c(1:3)])
  mesh = list()
  if(type == "alpha"){
    library(alphashape3d)
    library(rgl)
    sh=ashape3d(coords, 1.0, pert = FALSE, eps = 1e-09);
    mesh[[1]] = as.mesh3d(sh, triangles=T);
  }else if(type == "ellipse"){
    library(rgl);
    pos=cov(coords, y = NULL, use = "everything");
    mesh[[1]] = ellipse3d(x=as.matrix(pos), level=level);
  }else{
    library(ks);
    res=kde(coords);
    r = plot(res, cont=level*100, display="rgl");
    sc = scene3d();
    mesh = sc$objects;
  }
  library(RJSONIO);
  sink(filenm);
  cat(toJSON(mesh));
  sink();
  return(filenm);
}

.get.rdt.set <- function(){
  return(qs::qread("rdt.set.qs"));
}

.set.rdt.set <- function(my.set){
  qs::qsave(my.set, file="rdt.set.qs");
}




