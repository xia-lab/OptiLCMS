#' @title testOptiLCMS
#' @description testOptiLCMS
#' @author zhiqiang Pang
#' @export
testOptiLCMS <- function(){
  print("Running Successfully !")
}


MessageOutput(
  mes = paste0("Run isotope peak annotation.."),
  ecol = "\n",
  progress = NULL
)

.optimize_switch <<- .on.public.web <<- F;
Params <- SetPeakParam();

mSet<-InitDataObjects("spec", "raw", FALSE)
mSet <- ImportRawMSData(mSet, "/home/glassfish/projects/MetaboDemoRawData/upload/", mode = "onDisk",plotSettings = SetPlotParam(Plot = F))
