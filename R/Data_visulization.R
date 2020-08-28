#' Data inspectation
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
#' @author Zhiqiang Pang \email{zhiqiang.pang@mail.mcgill.ca} Jeff Xia \email{jeff.xia@mcgill.ca}
#' Mcgill University
#' License: GNU GPL (>= 2)

PerformDataInspect<-function(datapath=NULL, rt.range, mz.range, dimension="3D", res=100){
  
  
  if(datapath == "null" | is.null(datapath)){
    if (.on.public.web & dir.exists("upload/QC")){
      datapath <- "upload/QC";
      datapath <- paste0(fullUserPath, datapath);
    } else if (.on.public.web){
      datapath <- paste0("upload/",list.files("upload", recursive=T)[1]);
      datapath <- paste0(fullUserPath, datapath);
    } else {
      print("Local Inspectation !");
    }
  } else {
    
    files <- list.files(paste0(getwd(),"/upload/"), full.names = T, recursive = T);
    
    if(isEmpty(files)){ # Handle the example issue - show other files
      files <- list.files("/home/glassfish/projects/MetaboDemoRawData/upload/", full.names = T, recursive = T);
      datapath =  files[datapath == basename(files)];
    } else { # Handle the regular data insepctation
      datapath =  files[datapath == basename(files)];
    }
    
  }
  
  suppressMessages(require("MSnbase"))
  
  
  if(basename(datapath) == "NA"){ # Handle the example issue - default showing
    datapath <- "/home/glassfish/projects/MetaboDemoRawData/upload/QC";
  }
  
  if (!grepl(pattern = c("*.mzXML"),basename(datapath)) & 
      !grepl(pattern = c("*.mzML"),basename(datapath)) & 
      !grepl(pattern = c("*.CDF"),basename(datapath))& 
      !grepl(pattern = c("*.cdf"),basename(datapath))){
    
    print(paste("First file in ", datapath, " will be inspected !"))
    
    mzf<-list.files(datapath,recursive = F,full.names = T)[1]
    
  } else {
    mzf <- datapath;
  }
  
  if (grepl(pattern = c("*.CDF"),basename(mzf)) | grepl(pattern = c("*.cdf"),basename(mzf))){
    # centroid_cdf();
    return();
  }
  
  cat(mzf,"\n");
  cat(basename(mzf),"\n");
  ms <- openMSfile(mzf)
  hd <- header(ms)
  
  ms1 <- which(hd$msLevel == 1)
  
  if (missing(rt.range) | (rt.range[1] == 0 & rt.range[2] == 0)){
    
    rtsel <- hd$retentionTime[ms1] > min(hd$retentionTime) & hd$retentionTime[ms1] < max(hd$retentionTime);
    rt.extension <- F
    print(paste("RT range is:",min(hd$retentionTime), "and",max(hd$retentionTime),"seconds !"))
    
    
  } else{
    
    if (rt.range[2]<rt.range[1]){
      a1<-rt.range[1];
      rt.range[1]<-rt.range[2];
      rt.range[2]<-a1;
    } else if (rt.range[2] == rt.range[1]){
      rt.range[2] <- rt.range[2] + 1;
      rt.range[1] <- rt.range[1] - 1;
    }
    print(paste("RT range is:",rt.range[1], "and",rt.range[2],"seconds !"))
    
    rtsel <- hd$retentionTime[ms1] > rt.range[1] & hd$retentionTime[ms1] < rt.range[2]
    
    rt.min <- min(hd$retentionTime[ms1]);
    rt.max <- max(hd$retentionTime[ms1]);
    
    if (rt.range[1] < rt.min | rt.range[2] > rt.max){
      rt.extension <- T
    } else {
      rt.extension <- F
    }
    
  }
  
  if (missing(mz.range) | (mz.range[1] == 0 & mz.range[2] == 0)){
    
    min.mz<-min(hd$lowMZ);
    max.mz<-max(hd$highMZ);
    
    if (min.mz==0 & max.mz==0 | min.mz==max.mz){
      print("mz.range information is missing in your data file. mz between 100 and 1200 will be shown here !")
      min.mz <- 100;
      max.mz <- 1200;
    } else if (is.infinite(min.mz) | is.infinite(max.mz) | min.mz == -1 | max.mz == -1){
      print("mz.range information is missing in your data file. mz between 50 and 2000 will be shown here !")
      min.mz <- 50;
      max.mz <- 2000;
    }
    
    print(paste("MZ range is:",min.mz, "and",max.mz,"Thomson !"))
    
    res.mz<-(max.mz-min.mz)/res
    M <- MSmap(ms,
               ms1[rtsel],
               min.mz,
               max.mz,
               res.mz,
               hd,
               zeroIsNA = TRUE)
    
  } else{
    
    if (mz.range[2]<mz.range[1]){
      a1<-mz.range[1];
      mz.range[1]<-mz.range[2];
      mz.range[2]<-a1;
    } else if (mz.range[2] == mz.range[1]){
      mz.range[2] <- mz.range[2] + 0.01;
      mz.range[1] <- mz.range[1] - 0.01;
    }
    
    print(paste("MZ range is:",mz.range[1], "and",mz.range[2],"Thomson !"))
    
    res.mz<-(mz.range[2]-mz.range[1])/res
    M <- MSmap(ms,
               ms1[rtsel],
               mz.range[1],
               mz.range[2],
               res.mz,
               hd,
               zeroIsNA = TRUE)
  }
  
  
  if (rt.extension){
    
    if (min(M@rt) > rt.range[1]){
      M@rt <- c(rt.range[1],M@rt);
      M@map <- rbind(rep(NA,dim(M@map)[2]),M@map);
      M@ms <- c(M@ms,1)
    }
    
    if (max(M@rt) < rt.range[2]){
      M@rt <- c(M@rt,rt.range[2]);
      M@map <- rbind(M@map,rep(NA,dim(M@map)[2]));
      M@ms <- c(M@ms,1)
    }
  }
  
  if(missing(dimension)){
    dimension="3D";
  }
  
  FN <- tools::file_path_sans_ext(basename(mzf));
  
  filename <- paste0(FN,"_mz_",mz.range[1],"_",mz.range[2],"_RT_",rt.range[1],"_",rt.range[2],dimension,".png");
  
  Cairo::Cairo(file = filename, unit="in", dpi=72, width=7, height= 7, 
               type="png", bg="white");
  par(mfrow=c(1,2));
  if (dimension=="3D"){
    print(plot.MS_3D(M))
  } else {
    print(plot(M, aspect = 1, allTicks = FALSE))
  }
  
  dev.off();
  
  return(filename);
  
}

